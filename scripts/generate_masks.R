# # # # # # # #
# NDSI FUNCTION
compute_ndsi <- function(folder_path, region=NULL, scale_factor = 0.0000275, offset = -0.2, return.binary=TRUE) {
  #
  # scale and offset set for Landsat OLI (8)
  # binary threshold of 0.4 by default
  #
  # file paths (SR product)
  if(grepl("4|5", strsplit(basename(folder_path), "_")[[1]][1])){
    pattern_ext="SR_B(2|5).TIF$"; b1="SR_B2";b2="SR_B5"
  }else if(grepl("8|9", strsplit(basename(folder_path), "_")[[1]][1])){
    pattern_ext="SR_B(3|6).TIF$"; b1="SR_B3";b2="SR_B6"}else{stop("unrecognized sensor in folder_path")}
  files <- list.files(folder_path, pattern = pattern_ext, full.names = TRUE)
  # band assignment
  green <- raster(files[grep(b1, files)])
  swir <- raster(files[grep(b2, files)])
  if(!is.null(region)){green <- crop(green, region);swir <- crop(swir, region)}
  # scale factor 
  green <- green * scale_factor + offset
  swir <- swir * scale_factor + offset
  # compute and filter for realistic values
  ndsi <- (green - swir) / (green + swir)
  # ndsi[ndsi<(-1) | ndsi>1] = NA
  if(return.binary){ndsi <- ndsi >= 0.4}
  return(ndsi)
}

# # # # # # # #
# NDWI FUNCTION
compute_ndwi <- function(folder_path, region=NULL, scale_factor = 0.0000275, offset = -0.2, return.binary=TRUE){
  # same as above but for water
  if(grepl("4|5", strsplit(basename(folder_path), "_")[[1]][1])){
    pattern_ext="SR_B(2|4).TIF$"; b1="SR_B2";b2="SR_B4"
  }else if(grepl("8|9", strsplit(basename(folder_path), "_")[[1]][1])){
    pattern_ext="SR_B(3|5).TIF$"; b1="SR_B3";b2="SR_B5"}else{stop("unrecognized sensor in folder_path")}

  files <- list.files(folder_path, pattern = pattern_ext, full.names = TRUE)
  # band assignment
  green <- raster(files[grep(b1, files)])
  nir <- raster(files[grep(b2, files)])
  if(!is.null(region)){green <- crop(green, region);nir <- crop(nir, region)}
  # scale factor 
  green <- green * scale_factor + offset
  nir <- nir * scale_factor + offset
  # compute and filter for realistic values
  ndwi <- (green - nir) / (green + nir)
  # ndwi[ndwi<(-1) | ndwi>1] = NA
  # > 0.3 high probability of open water | 0-0.3 shallow water or wetlands | < 0 non-water features
  if(return.binary){ndwi <- ndwi >= 0.3} 
  return(ndwi)
}

# # # # # # # # # # #
# CFMASK BIT FUNCTION
classify_from_bits <- function(folder_path, region=NULL) {
  # 
  # https://www.usgs.gov/landsat-missions/landsat-collection-2-quality-assessment-bands
  # compute high confidence (clouds, etc.)
  #
  qa_file <- list.files(folder_path, pattern = "QA_PIXEL.TIF$", full.names = TRUE)
  if(!is.null(region)){raster_layer <- crop(raster(qa_file), region)}else{raster_layer <- raster(qa_file)}
  
  # Define the reclassification function
  reclassify_fun <- function(x) {
    # Extract individual bits (1-8)
    fill  <- bitwAnd(bitwShiftR(x, 0), 1)   # Bit 1
    dilate <- bitwAnd(bitwShiftR(x, 1), 1)  # Bit 2
    cir   <- bitwAnd(bitwShiftR(x, 2), 1)   # Bit 3
    cld   <- bitwAnd(bitwShiftR(x, 3), 1)   # Bit 4
    cshd  <- bitwAnd(bitwShiftR(x, 4), 1)   # Bit 5
    sno   <- bitwAnd(bitwShiftR(x, 5), 1)   # Bit 6
    clear <- bitwAnd(bitwShiftR(x, 6), 1)   # Bit 7
    wat   <- bitwAnd(bitwShiftR(x, 7), 1)   # Bit 8
    
    # Extract the 2-bit confidence values (bits 9-16)
    cld_conf    <- bitwAnd(bitwShiftR(x, 8), 3)  # Bits 9-10
    cldshd_conf <- bitwAnd(bitwShiftR(x, 10), 3) # Bits 11-12
    sno_conf    <- bitwAnd(bitwShiftR(x, 12), 3) # Bits 13-14
    cir_conf    <- bitwAnd(bitwShiftR(x, 14), 3) # Bits 15-16
    
    result <- ifelse(clear == 1 & cld_conf == 1 & cldshd_conf == 1 & cir_conf == 1 & sno_conf == 1, 1,
                     ifelse(wat == 1 & cld_conf == 1 & cldshd_conf == 1 & cir_conf == 1 & sno_conf == 1, 2,
                            ifelse( cld == 1 & (cld_conf == 3), 3,
                                    ifelse(cshd == 1 & (cldshd_conf == 3), 4,
                                           ifelse(cir == 1 & (cir_conf == 3), 5,
                                                  ifelse( cld == 1 & (cld_conf == 2), 6,
                                                          ifelse(cshd == 1 & (cldshd_conf == 2), 7,
                                                                 ifelse(cir == 1 & (cir_conf == 2), 8,
                                                                        ifelse(sno == 1 & sno_conf == 3, 9, NA)))))))))
    
    # 1:clear, 2:water, 3:cld_high, 4:cldshd_high, 5:cir_high, 6:cld_mid, 7:cldshd_mid, 8:cir_mid, 9:snow_high
    return(result)
  }
  # Apply to raster using terra::app (vectorized)
  out <- calc(raster_layer, fun = reclassify_fun)
  return(out)
}

# # # # # # # #
# LST FUNCTION
process_LST <- function(folder_path, region=NULL, save_temp = NULL){
  #
  # https://www.usgs.gov/landsat-missions/landsat-collection-2-level-2-science-products
  # scaling factors for OLI 
  #
  if(grepl("4|5", strsplit(basename(folder_path), "_")[[1]][1])){
    t_band="B6.TIF$"
  }else if(grepl("8|9", strsplit(basename(folder_path), "_")[[1]][1])){
    t_band="B10.TIF$"}else{stop("unrecognized sensor in folder_path")}
  band_10 <- raster(list.files(folder_path, pattern = t_band, full.names = TRUE))
  band_error <- raster(list.files(folder_path, pattern = "ST_QA", full.names = TRUE))
  
  if(!is.null(region)){
    band_10 <- crop(band_10, region) 
    band_error <- crop(band_error, region) 
  }
  band_error_scaled <- (band_error*0.01) 
  band_10_scaled <- (band_10*0.00341802+149)
  
  # scale offset - temperatures at or below freezing
  lst0C <- band_10_scaled < (273.15+band_error_scaled)
  if(!is.null(save_temp)){
    writeRaster(lst0C, filename = save_temp)
    return(save_temp)
  }else{
    return(lst0C)
  }
}
