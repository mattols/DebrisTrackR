###############################################
make.shade <- function(folder_path, region = NULL, shade_size_filter = 10, #40, 5?
                       dilate_shadows=TRUE){
  ## SHADE
  # check dem
  if(!exists('dem')){stop('dem not loaded for tile\n', basename(dirname(folder_path)))}
  # crop
  if(!is.null(region)){
    dem <- crop(dem, region)
    dem <- resample(dem, crop(raster(list.files(folder_path, pattern = "SR_B2.TIF$",full.names = TRUE)), region) )
  }
  # s_elv, s_az, datetime
  scene_vars <- extract.solar.vars(folder_path)
  # ZEN = 90 - s_elv ; AZ = s_az
  ZEN=90-scene_vars[[1]];AZ=scene_vars[[2]]
  sol_vect <- insol::normalvector(ZEN,AZ)
  shd <- insol::doshade(dem,sol_vect)
  shd2 = shd;shd2[is.na(dem)] = NA
  # filter pixels based on size
  shclump <-  clump(!shd2, directions=8)
  f<-as.data.frame(freq(shclump))
  exludeShade <- f$value[which(f$count <= shade_size_filter)]
  shfilter <- shclump
  shfilter[shclump %in% exludeShade] <- NA 
  # dilate?
  if(dilate_shadows){shfilter <- shade.dilation(folder_path, shfilter, window_size = 9, region = region)}
  shf <- (!is.na(shfilter))*!is.na(dem)
  shdd <- !shf
  # # ## # # #
  if(FALSE){
    # plot shade dilation
    plotRGB(crop(stack(rev(list.files(folder_path, pattern = "SR_B[2-4].TIF$",full.names = TRUE))), p2), stretch='lin')
    plot(crop(shd, p2), add=T, col=adjustcolor(c("green",NA), 0.6))
    plot(crop(shdd, p2), add=T, col=adjustcolor(c("red",NA), 0.6))
  }
  return(shdd)
}

shade.dilation <- function(folder_path, shfilter, window_size = 9, save_file=NULL, region=tile_roi){
  # # # # #
  # shade dilation - uses 9x9 buffer and other shaded pixel values to find additional shade
  shdilate <- focal((!is.na(shfilter)), w=matrix(1,nrow=window_size,ncol=window_size))>0
  tband <- resample(crop(raster(list.files(folder_path,pattern="B2",full.names=T)), region), shdilate) # blv
  bluemasks <- mask(tband,shdilate, maskvalue=0)
  #
  # using a gaussian mixture model
  # mean of shade distribution +2 standard deviations
  bmix <- mixtools::normalmixEM(getValues(bluemasks)%>%na.omit(), lambda=0.5) 
  bmax_new <- (bmix$mu[1] + (bmix$sigma[1]*2))
  
  shdilate_r <-(!is.na(shfilter)) | (bluemasks<bmax_new); shdilate_r[shdilate_r==0]<-NA
  
  if (!is.null(save_file)) {
    writeRaster(shdilate_r, filename = save_filename, format = "raster", overwrite = TRUE)
  }
  if(FALSE){
    # plot shade dilation
    plotRGB(crop(stack(rev(list.files(folder_path, pattern = "SR_B[2-4].TIF$",
                                      full.names = TRUE))), p2), stretch='lin')
    # too slow
    # feature_line1 = rasterToPolygons(crop(shdilate, p2)==1,dissolve=T)
    # lines(feature_line1, col=adjustcolor('yellow', 0.8), lwd=1.2)
    # plot(crop(shdilate, p2), add=T, col=adjustcolor(c(NA, "yellow"), 0.6))
    plot(crop(shdilate_r, p2), add=T, col=adjustcolor(c(NA, "red"), 0.6))
    # plot(crop(shd, p2), add=T, col=adjustcolor(c("green",NA), 0.6))
    plot(crop(shdd, p2), add=T, col=adjustcolor(c("green",NA), 0.8))
  }
  return(shdilate_r)
}


###############################################
extract.solar.vars <- function(folder_path){
  # GET SOLAR VARS
  ## METADATA
  mtl_path <- list.files(folder_path, pattern = "MTL.txt", full.names = TRUE)
  if (length(mtl_path) > 0) {
    con <- file(mtl_path[1], open = "r")  # use mtl_path[1] to get the actual file path
    line <- readLines(con)
    close(con)
  } else {warning("No matching 'MTL.txt' in folder_path")}
  ## SOLAR POSITION AND DATETIME
  s_elv <- as.numeric(strsplit(trimws(grep("SUN_ELEVATION", line, value = T)), "= ")[[1]][2])
  s_az <- as.numeric(strsplit(trimws(grep("SUN_AZ", line, value = T)), "= ")[[1]][2])
  date_time <- paste(strsplit(trimws(grep("DATE_ACQUI", line, value = T)), "= ")[[1]][2],
                     substr(strsplit(strsplit(trimws(grep("TIME", line, value = T)), "= ")[[1]][2], "\"")[[1]][2], 1, 8))
  datetime <- as.POSIXct(date_time, tz = "UTC")
  return(list(s_elv, s_az, datetime))
}
