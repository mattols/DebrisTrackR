##########################################
# FUSION FUNCTION!
glacier_debris_fusion <- function(tile_region = '140041', T1=TRUE, 
                                  time1=2013:2015, time2=2022:2024, months=7:9, region = tile_roi,
                                  dest_path){
  strt <- Sys.time()
  message('Initialize glacier debris fusion classification code...')
  
  # TEST IF EXISTS
  fusion_save_name <- file.path(dest_path, tile_region, paste0(tile_region, '_T1_max_', 
                                           min(time1), max(time1), '_',min(time2), max(time2), "_m",
                                           min(months), max(months), "_newFusion1.grd"))
  if(file.exists(fusion_save_name)){message("Classification fusion file exists for fime period\n");print(fusion_save_name);return(NULL)}
  
  # calculate NDSI mask class
  debris_change_wrapper(tile_region=tile_region, T1=T1,
                        years1=time1, years2=time2, months=months,
                        region = region,
                        dest_path = dest_path)
  # determine 
  
  # run land temp minimum
  temp_filename = file.path(file.path(dest_path, tile_region), 
                            paste0(tile_region, "_minLand_stk_", min(time1), max(time1), '_',min(time2), max(time2), "_m",
                                   min(months), max(months), ".grd"))
  if(!file.exists(temp_filename)){
    if(!exists('parent_dir')){
      parent_dir <- file.path("~/olson/glacier/hma-debris/Landsat/Landsat08", tile_region)
    }
    lst_min = land_temp_diff(time1, time2, months, parent_dir, region,
                             save_file = temp_filename)
  }
  # returns 2 stack
  lst_min = stack(temp_filename)
  lnd_diff <- lst_min[[1]] - lst_min[[2]]
  lnd_mult <- lst_min[[1]] * lst_min[[2]]
  
  # fusion logic
  # class_t1 <- list.files(file.path(dest_path, tile_region), full.names = TRUE,
  #                        pattern = paste0("class_mode.*", min(time1), ".*grd"))
  class_t1 <- raster(list.files(file.path(dest_path, tile_region), full.names = TRUE,
                         pattern = paste0("class_max_", min(time1),max(time1), '_m',min(months),max(months),".grd")))
  class_t2 <- raster(list.files(file.path(dest_path, tile_region), full.names = TRUE,
                         pattern = paste0("class_max_", min(time2),max(time2), '_m',min(months),max(months),".grd")))
  change_class <- raster(list.files(file.path(dest_path, tile_region), full.names = TRUE,
                                pattern = paste0("change_max_",min(time1),max(time1),
                                                 '_',min(time2),max(time2), '_m',min(months),max(months),".grd")))
  
  if(FALSE){
    # SUBSET GMASK BOUNDARY - Temporary fix
    # FOR NGOZUMPA & EVEREST
    r <- raster(list.files(list.files(parent_dir, full.names = T, pattern ='.*T[1-2]$')[1], pattern = 'SR_B4', full.names = T))
    rroi <- crop(r, region)
    gmask <- sf::st_as_sf(gpoly[c(55,78),]) %>% 
      fasterize::fasterize(.,rroi) %>% raster::crop(.,rroi) %>% raster::mask(.,rroi)
  }
  
  # boundaries
  bd<-gmask%>%boundaries(.,type="inner")
  # thick border (60-meters from edge removed)
  bd2=bd;bd2[bd2==1]<-NA
  bd2<-bd2%>%boundaries(.,type="inner")
  # bd3=bd|bd2
  bd3=bd2;bd3[bd3==1]<-NA
  bd4=!is.na(bd3)
  
  print("Starting fusion classification...")
  
  # FINAL CLASSIFICATION
  result <- try({
    cstk <- stack(change_class, lnd_diff, lnd_mult, bd4, gmask)
  }, silent = TRUE)
  
  if (inherits(result, "try-error")) {
    # do something else if it failed
    message("Realigning glacier mask extent to image classes...")
    result2 <- try({
      cstk <- stack(change_class, lnd_diff, lnd_mult, extend(bd4, change_class), extend(gmask, change_class))
    }, silent = TRUE)
    if (inherits(result2, "try-error")) {
      # do something else if it failed
      message("Realigning all to main change stack...")
      # COULD INCLUDE LOGIC TO DETERMINE WHETHER TO CROP OR EXPAND
      result3 <- try({
        cstk <- stack(change_class, crop(lnd_diff, change_class), crop(lnd_mult, change_class), extend(bd4, change_class), extend(gmask, change_class))
      }, silent = TRUE)
      if (inherits(result3, "try-error")) {
        message("Switch extend to crop for glacier mask...")
        result4 <- try({
          cstk <- stack(change_class, crop(lnd_diff, change_class), crop(lnd_mult, change_class), crop(bd4, change_class), crop(gmask, change_class))
        }, silent = TRUE)
        if (inherits(result4, "try-error")) {
          message("Final attempt to find common extent...")
          e1 <- extent(change_class);e2 <- extent(lnd_diff);e3 <- extent(gmask)
          xmin_union <- min(e1@xmin, e2@xmin, e3@xmin);xmax_union <- max(e1@xmax, e2@xmax, e3@xmax)
          ymin_union <- min(e1@ymin, e2@ymin, e3@ymin);ymax_union <- max(e1@ymax, e2@ymax, e3@ymax)
          # Create a new extent object
          common_extent <- extent(xmin_union, xmax_union, ymin_union, ymax_union)
          cstk <- stack(extend(change_class, common_extent), extend(lnd_diff, common_extent), extend(lnd_mult, common_extent), 
                        extend(bd4, common_extent), extend(gmask, common_extent))
        }
      }
    }
  }
  overlay(cstk, 
          fun=function(a,b,c,d,e) ifelse((a==1 | c==1) & e==1, 1, 
                                         ifelse(a==3 & b==1 & d==1, 3, 
                                                ifelse(a==3 & b==1 & e==1, 3.5,
                                                       ifelse((a==2 | c==0) & d==1, 2,
                                                              ifelse(a==4 & b<0 & d==1, 4, NA))))),
  filename=fusion_save_name)
  
  message("Completed fusion")
  print(Sys.time() - strt)
  cat("")
}
