land_temp_273_faster <- function(years, months, parent_dir, region=tile_roi, T1=TRUE){
  #
  # save files to tempdir()
  # 
  message("  ...processing LST min\n     years: ", years, " | months: ", months, "\n")
  
  scene_folders <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
  if(T1){scene_folders <- scene_folders[grepl("_T1", scene_folders)]}
  
  date_str <- stringr::str_extract(basename(scene_folders), "\\d{8}")
  year <- as.integer(substr(date_str, 1, 4))
  month <- as.integer(substr(date_str, 5, 6))
  idx_subset <- (year >= min(years) & year <= max(years)) & (month >= min(months) & month <= max(months))
  
  scene_folders <- scene_folders[idx_subset]
  # date_sub <- date_str[idx_subset]
  output_files <- character()
  # message("   ntiles: ", length(scene_folders))
  
  for (folder in scene_folders) {
    idx <- match(folder, scene_folders)
    # message(idx, " of ", length(scene_folders))
    out_file <- try(
      process_LST(folder, region, 
                  save_temp = file.path(tempdir(),
                                        paste0(basename(parent_dir),'_lnd_', 
                                               paste0(years, collapse=""),'_i',idx,'_',paste0(sample(letters, 5), collapse=''),'.grd'))), 
      silent = TRUE)
    if (!inherits(out_file, "try-error")) {
      output_files <- c(output_files, out_file)
    }
  }
  
  out_file_stk <- try(
    lndstk <- calc(stack(output_files), min,na.rm=T), 
    silent = TRUE)
  if (inherits(out_file_stk, "try-error")) {
    tile_ex <- do.call(raster::merge, lapply(lapply(output_files, raster::raster), raster::extent))
    lndstk <- calc(stack(lapply(lapply(output_files, raster::raster), function(x)  extend(x, tile_ex))), min, na.rm=T)
  }
  unlink(output_files)
  # remove .gri files as well 
  unlink(gsub('.grd', '.gri', output_files))
  return(lndstk)
}

land_temp_diff <- function(years1, years2, months, parent_dir, region=tile_roi, save_file=NULL){
  #
  # Stack thermal data
  # 
  message("Converting time-series land surface temperature minimum... \n")
  # lnd1 <- land_temp_273(years1,months,parent_dir,region)
  # lndmin1 <- calc(brick(lnd1),min,na.rm=T)
  # years1=time1;years2=time2
  # YEAR 1
  if(max(years1)<2013){pdir1=gsub('Landsat08', 'LTM', parent_dir)}else{pdir1=parent_dir}
  lndmin1 <- land_temp_273_faster(years1,months,pdir1,region=region)
  # YEAR 2
  if(max(years2)<2013){pdir2=gsub('Landsat08', 'LTM', parent_dir)}else{pdir2=parent_dir}
  lndmin2 <- land_temp_273_faster(years2,months,pdir2,region=region)
  # lnd_diff <- lndmin1 - lndmin2
  out_file_stk <- try(
    lnd_stk <- stack(lndmin1, lndmin2), 
    silent = TRUE)
  if (inherits(out_file_stk, "try-error")) {
    max_ext <- raster::union(extent(lndmin1), extent(lndmin2))
    r1_ext <- extend(lndmin1, max_ext)
    r2_ext <- extend(lndmin2, max_ext)
    lnd_stk <- stack(r1_ext, r2_ext)
  }
  
  if(!is.null(save_file)){
    writeRaster(lnd_stk, filename = save_file, format = "raster", overwrite = TRUE)
  }
  # return(lnd_diff)
  return(save_file)
}
