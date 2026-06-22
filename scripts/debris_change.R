# # # # # # # # # # # #
# TIME WRAPPER FUNCTION 
debris_change_wrapper <- function(tile_region, T1=TRUE,
                                  years1, years2, months, region = tile_roi,
                                  dest_path = '/uufs/chpc.utah.edu/common/home/u1037042/olson/glacier/hma-debris/output/class-new-072025'
                                  ) {
  strt <- Sys.time()
  cat("Debris change for tile", tile_region, "\n", min(years1), max(years1), "to", min(years2), max(years2), "\n")
  
  # run tile manager
  # if(!exists('gmask')){}
  path_tile_manager(tile_path=tile_region, gpoly_type="SAMRGI")
  parent_dir <<- file.path("~/olson/glacier/hma-debris/Landsat/Landsat08", tile_region)
  # parent_dirTM <<- file.path('~/olson/glacier/hma-debris/Landsat/LTM', tile_path)
  
  if(!dir.exists(dest_path)){dir.create(dest_path)}
  dest_tile_path <- file.path(dest_path, tile_region)
  if(!dir.exists(dest_tile_path)){dir.create(dest_tile_path)}
  
  # Build filenames for sno_prob
  # SAVE NAMES
  ctype <<- "max" # 'mode'
  prob1_file <- file.path(dest_tile_path, paste0(tile_region, "_T1_class_", ctype, "_", min(years1), "", max(years1), "_m", min(months), max(months), ".grd"))
  prob2_file <- file.path(dest_tile_path, paste0(tile_region, "_T1_class_", ctype, "_", min(years2), "", max(years2), "_m", min(months), max(months), ".grd"))
  change_filename <- file.path(dest_tile_path, paste0(tile_region, "_T1_change_", ctype, "_", min(years1), 
                                                      max(years1), "_", min(years2), max(years2), "_m", min(months), max(months),".grd"))
  
  if(!file.exists(prob1_file)){
    # Time 1
    if(max(years1)<2013){pdir1=gsub('Landsat08', 'LTM', parent_dir)}else{pdir1=parent_dir}
    layer1 <- process_debris_date_range(
      parent_dir = pdir1,
      T1 = T1,
      region = region,
      start_year = min(years1),
      end_year = max(years1),
      start_month = min(months),
      end_month = max(months),
      save_filename = prob1_file,
      dest_path = dest_path
    )
  }else{layer1 <- raster(prob1_file)}
  if(!file.exists(prob2_file)){
    # Time 2
    if(max(years2)<2013){pdir2=gsub('Landsat08', 'LTM', parent_dir)}else{pdir2=parent_dir}
    layer2 <- process_debris_date_range(
      parent_dir = pdir2,
      T1 = T1,
      region = region,
      start_year = min(years2),
      end_year = max(years2),
      start_month = min(months),
      end_month = max(months),
      save_filename = prob2_file,
      dest_path = dest_path
    )
  }else{layer2 <- raster(prob2_file)}
  # CHANGE CLASS
  if(file.exists(change_filename)){
    return(raster(change_filename))
  }else{
    class_change <- reclassify(layer1 - (layer2*4), 
                               rcl=matrix(c(-3, 1,
                                            -7, 3,
                                            -2, 4,
                                            -6, 2), ncol=2,byrow=T),
                               filename = change_filename)
    
    cat("Change classification complete.\n")
    print(Sys.time() - strt)
    return(class_change)
  }
  
}
