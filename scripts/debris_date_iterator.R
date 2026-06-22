# # # # # # # # # # # # #
# DATE ITERATOR
process_debris_date_range <- function(parent_dir, T1=TRUE, region = NULL,
                                      start_year = 2013,
                                      end_year = 2014,
                                      start_month = 7,
                                      end_month = 10,
                                      save_filename = NULL,
                                      dest_path) {
  strt <- Sys.time()
  cat("...processing dates\n", start_year, "to", end_year, "\n")
  
  scene_folders <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
  if(T1){scene_folders <- scene_folders[grepl("_T1", scene_folders)]}
  
  date_str <- stringr::str_extract(basename(scene_folders), "\\d{8}")
  year <- as.integer(substr(date_str, 1, 4))
  month <- as.integer(substr(date_str, 5, 6))
  idx_subset <- (year >= start_year & year <= end_year) & (month >= start_month & month <= end_month)
  
  scene_folders <- scene_folders[idx_subset]
  date_sub <- date_str[idx_subset]
  output_files <- character()
  message("   ntiles: ", length(scene_folders))
  
  # # # # # # SAVE DF # # # # #
  
  # if(exists('df_tile_info')){
  #   df_tile_info1 <<- df_tile_info
  # }else{}
  df_tile_info <- data.frame(tile = rep(basename(parent_dir), length(scene_folders)),
                             date = date_sub)
  df_tile_info[,3:(3+12)] <- NA
  names(df_tile_info) <- c( "tile", 'date', "folder", "scene_cloud",   
                            "snow_pct","water_pct","shd_pct","fmask_cld",    
                            "fmask_cldshd", "fmask_cir", "fmask_snow", "lst_pct",        
                            "class02_1_pct", "class02_2_pct", "class02_unclass")
  df_tile_info <<- df_tile_info
  
  # # # # # # # #
  
  for (folder in scene_folders) {
    idx <- match(folder, scene_folders)
    message("Processing: ", folder, " | ", idx, " of ", length(scene_folders))
    
    out_file <- try(process_snowdebris(folder, region, out_dir = tempdir()), silent = TRUE)
    if (!inherits(out_file, "try-error")) {
      output_files <- c(output_files, out_file)
    }
  }
  
  # STACK DFs
  # df_tile_info_all <<- rbind(df_tile_info1, df_tile_info)
  # save scene info
  save_csv_file = file.path(dest_path, basename(parent_dir), 
                            paste0('df_tile_info_', basename(parent_dir), 
                                   '_T1_max_', 
                                   start_year, end_year, "_m",
                                   start_month, end_month, ".csv"))
  write.csv(df_tile_info, 
            save_csv_file, row.names = FALSE)
                                                                       
  
  # Stack or brick the rasters from disk
  out_file_stk <- try(
    raster_stack <- stack(output_files), 
    silent = TRUE)
  if (inherits(out_file_stk, "try-error")) {
    tile_ex <- do.call(raster::merge, lapply(lapply(output_files, raster::raster), raster::extent))
    raster_stack <- stack(lapply(lapply(output_files, raster::raster), function(x)  extend(x, tile_ex)))
  }
  
  
  # CURRENT
  sno_prob = calc(raster_stack, fun=max, na.rm=T)
  # sno_prob = calc(raster_stack, fun=modal, na.rm=T)
  
  # Compute snow presence probability (max of 1/2 class codes)
  # sno_prob <- calc(raster_stack, fun = function(x) max(x == 1, na.rm = TRUE))
  
  # my_brick = brick(new_tiles)
  # my_brick[my_brick == 0] <- NA 
  # sno_prob2 = calc(my_brick, fun=modal, na.rm=T)
  
  # Save to file
  if (!is.null(save_filename)) {
    # writeRaster(sno_prob, filename = save_filename, format = "GTiff", overwrite = TRUE)
    writeRaster(sno_prob, filename = save_filename, format = "raster", overwrite = TRUE)
  }
  
  # remove temp files
  unlink(output_files)
  
  cat("Finished\n")
  print(Sys.time() - strt)
  return(sno_prob)
}
