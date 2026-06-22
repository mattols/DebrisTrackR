# update does not return path



# # # # # # # # # # # # #
# CLASSIFICATION FUNCTION
process_snowdebris <- function(folder_path, region = NULL, out_dir = tempdir()) {
  # NDSI | NDWI
  binary_snow <- compute_ndsi(folder_path, region)
  binary_water <- compute_ndwi(folder_path, region)
  # CFMask | 1:clear, 2:water, 3:cld_high, 4:cldshd_high, 5:cir_high, 6:cld_mid, 7:cldshd_mid, 8:cir_mid, 9:snow_high
  fmask <- classify_from_bits(folder_path, region)
  # LST <=0
  lst_binary <- process_LST(folder_path, region)
  # Shading | ray tracing
  shadow_mask <- make.shade(folder_path, region)
  
  # Classification logic
  final_class_fun <- function(snow, water, fmask_val, lst, shd) {
    # test 02
    class01 <- (snow == 1) & (water != 1) & (fmask_val != 3) & (lst == 1) & shd
    class02 <- (lst != 1) & (snow != 1) & (water != 1) & shd
    # test 04
    # class01 <- (snow == 1) & (water != 1) & (fmask_val != 3) & (fmask_val != 4) & (lst == 1)
    # class02 <- (lst != 1) & (snow != 1) & (water != 1) & (fmask_val != 4)
    
    out <- rep(NA_integer_, length(snow))
    out[class01] <- 1
    out[class02] <- 2
    return(out)
  }
  
  # File path to save
  scene_id <- basename(folder_path)
  out_file <- file.path(out_dir, paste0("final_class02_", scene_id, ".tif"))
  
  # Use overlay with on-disk write
  final_class02 <- overlay(binary_snow, binary_water, fmask, lst_binary, shadow_mask,
                           fun = final_class_fun,
                           filename = out_file,
                           datatype = "INT1U", overwrite = TRUE)
  
  
  # # # # # # # # # #
  # # # # # SAVE
  # Extract folder name and date
  folder_name <- basename(folder_path)
  date_str <- sub(".*_(\\d{8})_\\d{8}_.*", "\\1", folder_name)
  acquisition_date <- as.Date(date_str, format="%Y%m%d")
  
  percent_in_mask <- function(r, mask_idx, total_masked, target_vals) {
    vals <- getValues(r)
    masked_vals <- vals[mask_idx]
    sum(masked_vals == target_vals, na.rm = TRUE) / total_masked * 100
  }
  
  # Try-catch to avoid crashes from bad files
  try({
    
    # Load gmask and get masked indices
    # gmask <- raster(file.path(folder, "gmask.tif"))
    gmask_vals <- getValues(gmask)
    gmask_idx <- which(gmask_vals == 1)
    total_masked <- length(gmask_idx)
    
    if (total_masked == 0) next  # skip if no masked area

    # Load other rasters and compute
    snow_pct  <- percent_in_mask(binary_snow,  gmask_idx, total_masked, 1)
    water_pct <- percent_in_mask(binary_water, gmask_idx, total_masked, 1)
    shd_pct <- percent_in_mask(shadow_mask, gmask_idx, total_masked, 0)
    
    fmask_3_pct <- percent_in_mask(fmask, gmask_idx, total_masked, 3)
    fmask_4_pct <- percent_in_mask(fmask, gmask_idx, total_masked, 4)
    fmask_5_pct <- percent_in_mask(fmask, gmask_idx, total_masked, 5)
    fmask_9_pct <- percent_in_mask(fmask, gmask_idx, total_masked, 9)
    
    lst_pct <- percent_in_mask(lst_binary, gmask_idx, total_masked, 1)

    # Compute final_class02 percentages
    final_vals <- getValues(final_class02)
    masked_final_vals <- final_vals[gmask_idx]
    
    class02_1_pct <- sum(masked_final_vals == 1, na.rm = TRUE) / total_masked * 100
    class02_2_pct <- sum(masked_final_vals == 2, na.rm = TRUE) / total_masked * 100
    class02_unclass <- sum(is.na(masked_final_vals)) / total_masked * 100
    
    safe_extract_value <- function(pattern, lines, convert_fn = identity) {
      line <- grep(pattern, lines, value = TRUE)
      if (length(line) == 0) return(NA)
      value <- sub(".*=\\s*", "", line[1])
      value <- gsub('"', '', value)
      return(convert_fn(value))
    }
    
    lines <- readLines(list.files(folder_path, pattern = "MTL.txt$", full.names = TRUE), warn = FALSE)
    scene_cloud <- safe_extract_value("CLOUD_COVER\\s*=\\s*", lines, as.numeric)
    
    df_folder_i <- data.frame(
      folder = folder_name,
      scene_cloud,
      snow_pct,
      water_pct,
      shd_pct,
      fmask_3_pct,
      fmask_4_pct,
      fmask_5_pct,
      fmask_9_pct,
      lst_pct,
      class02_1_pct,
      class02_2_pct,
      class02_unclass,
      stringsAsFactors = FALSE
    )
    
    didx <- match(date_str, df_tile_info$date)
    df_tile_info[didx,] <- cbind(df_tile_info[didx,1:2], df_folder_i)
    df_tile_info <<- df_tile_info
    
  }, silent = TRUE)
  
  
  # # # # # # # # # # #
  
  return(out_file)  # return path, not raster object
}
