##########################################
# FUSION FUNCTION - MODIFIED
glacier_debris_fusion <- function(tile_region = '140041', T1=TRUE, 
                                  time1=2013:2015, time2=2023:2024, months=7:10, 
                                  dest_path){
  strt <- Sys.time()
  message('Initialize glacier debris fusion classification code...')
  
  # TEST IF EXISTS
  fusion_save_name <- file.path(dest_path, tile_region, paste0(tile_region, '_T1_max_', 
                                                               min(time1), max(time1), '_',min(time2), max(time2), "_m",
                                                               min(months), max(months), "_newFusionXb0.grd"))
  # if(file.exists(fusion_save_name)){message("Classification fusion file exists for fime period\n");print(fusion_save_name);return(NULL)}
  
  # calculate NDSI mask class
  debris_change_wrapper(tile_region=tile_region, T1=T1,
                        years1=time1, years2=time2, months=months,
                        dest_path = dest_path)
  # determine 
  
  # run land temp minimum
  temp_filename = file.path(file.path(dest_path, tile_region), 
                            paste0(tile_region, "_minLand_stk_", min(time1), max(time1), '_',min(time2), max(time2), "_m",
                                   min(months), max(months), ".grd"))
  if(!file.exists(temp_filename)){
    if(!exists('parent_dir')){parent_dir <- file.path("~/olson/glacier/hma-debris/Landsat/Landsat08", tile_region)}
    lst_min = land_temp_diff(time1, time2, months, parent_dir, tile_roi,
                             save_file = temp_filename)
  }
  # returns 2 stack
  lst_min = stack(temp_filename)
  lnd_diff <- lst_min[[1]] - lst_min[[2]]
  lnd_mult <- lst_min[[1]] + lst_min[[2]]*2 # 0:both debris; 3:both ice; 2:only b is 1; 1:only a is 1
  
  # fusion logic
  # class_t1 <- list.files(file.path(dest_path, tile_region), full.names = TRUE,
  #                        pattern = paste0("class_mode.*", min(time1), ".*grd"))
  #
  # class_t1 <- raster(list.files(file.path(dest_path, tile_region), full.names = TRUE,
  #                               pattern = paste0("class_max_", min(time1),max(time1), '_m',min(months),max(months),".grd")))
  # class_t2 <- raster(list.files(file.path(dest_path, tile_region), full.names = TRUE,
  #                               pattern = paste0("class_max_", min(time2),max(time2), '_m',min(months),max(months),".grd")))
  change_class <- raster(list.files(file.path(dest_path, tile_region), full.names = TRUE,
                                    pattern = paste0("change_max_",min(time1),max(time1),
                                                     '_',min(time2),max(time2), '_m',min(months),max(months),".grd")))
  
  # DEM
  if(!exists('dem')){
    path_tile_manager(tile_path=tile_region, gpoly_type="SAMRGI")
  }
  message('resampling DEM for slope correction...')
  dem2 = resample(dem, change_class)
  sl_thresh <- raster::terrain(dem2, opt='slope', unit='degrees') #< 25
  # sl_thresh <- raster::terrain(dem2, opt='slope', unit='degrees') < 35
  # Fusion2 and Fusion3 (ice pixel speckle) use slope angles less than 25 degrees
  # Fusion 4 and 5 use 30 degrees
  # Fusion 6 and 7 use 40 degree threshold (angle of repose)
  # Fusion 8 and 9 use 25 degrees - but also apply this to debris cover classification
    
  # boundaries
  # bd<-gmask%>%boundaries(.,type="inner")
  # border (30-meters from edge removed)
  # bd2=bd;bd2[bd2==1]<-NA
  # bd2<-bd2%>%boundaries(.,type="inner")
  
  # OLD (to 60m buffer only)
  # bd3=bd2;bd3[bd3==1]<-NA # turn off for 30
  # bd4=!is.na(bd3) # turn off for 30
  
  # bd4=!is.na(bd2) # turn ON for 30
  
  bd4 = gmask # turn ON for 0 buffer
  
  # NEW
  # bd2<-bd3%>%boundaries(.,type="inner")
  # bd3=bd2;bd3[bd3==1]<-NA # 90 meters
  # bd2<-bd3%>%boundaries(.,type="inner")
  # bd3=bd2;bd3[bd3==1]<-NA # 120 meters
  # bd2<-bd3%>%boundaries(.,type="inner")
  # bd3=bd2;bd3[bd3==1]<-NA # 150 meters
  # bd4=!is.na(bd3)
  
  
  
  print("Starting fusion classification...")
  
  # FINAL CLASSIFICATION
  result <- try({
    cstk <- stack(change_class, lnd_diff, lnd_mult, bd4, gmask, sl_thresh)
  }, silent = TRUE)
  
  if (inherits(result, "try-error")) {
    # do something else if it failed
    message("Realigning glacier mask extent to image classes...")
    
    result2 <- try({
      cstk <- stack(change_class, lnd_diff, lnd_mult, extend(bd4, change_class), extend(gmask, change_class), sl_thresh)
    }, silent = TRUE)
    
    if (inherits(result2, "try-error")) {
      # do something else if it failed
      message("Realigning all to main change stack...")
      # COULD INCLUDE LOGIC TO DETERMINE WHETHER TO CROP OR EXPAND
      cstk <- stack(change_class, crop(lnd_diff, change_class), crop(lnd_mult, change_class), extend(bd4, change_class), extend(gmask, change_class), sl_thresh )
    }
    
  }
  overlay(cstk, 
          fun=function(a,b,c,d,e,f) ifelse((a==1 | c==3) & e==1, 1, # ice
                                         ifelse(a==3 & b==1 & d==1 & (f<25), 3, # debris gain (f==1)
                                                ifelse(a==3 & b==1 & e==1 & (f<25), 3.5, # border correction
                                                       ifelse((a==2 | c==0) & e==1 & (f<35) & !(a==3 & b==1), 2, # debris
                                                              ifelse(a==4 & b<0 & d==1 & (f<25), 4, NA))))), # ice gain
          filename=fusion_save_name, overwrite=T)
  
  
  # filter ICE GAIN pixels based on size
  #
  cs2_mod <- raster(fusion_save_name)
  cs2_ig <- cs2_mod==4
  cs2_ig[cs2_ig==0]=NA
  # ICE gain limit
  csclump <-  clump(cs2_ig, directions=8)
  f<-as.data.frame(freq(csclump))
  exludePIX <- f$value[which(f$count < 5)] # filter areas with less than x pixels (9)
  cs2_mod[csclump %in% exludePIX] <- NA 
  # Boundary limit
  cs2_b <- cs2_mod==3.5
  cs2_b[cs2_b==0]=NA
  csclump <-  clump(cs2_b, directions=8)
  f<-as.data.frame(freq(csclump))
  exludePIX <- f$value[which(f$count < 5)] # filter areas with less than x pixels
  cs2_mod[csclump %in% exludePIX] <- NA 
  # Debris
  cs2_d <- cs2_mod==2
  cs2_d[cs2_d==0]=NA
  csclump <-  clump(cs2_d, directions=8)
  f<-as.data.frame(freq(csclump))
  exludePIX <- f$value[which(f$count < 5)] # filter areas with less than x pixels
  cs2_mod[csclump %in% exludePIX] <- NA 
  
  # SAVE
  # writeRaster(cs2_mod, gsub('Fusion2','Fusion3',fusion_save_name), overwrite=T)
  writeRaster(cs2_mod, gsub('Fusion8','Fusion9',fusion_save_name), overwrite=T)
  
  message("Completed NEW fusion")
  print(Sys.time() - strt)
  cat("")
}
