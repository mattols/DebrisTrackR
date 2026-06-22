#
# Supraglacial debris cover classification 
#  Thermal-optical multi-image fusion with shadow (DOSTS) and modified CFMask cloud handling
#  Olson M. - update 06/2026
# [Paper in review]
#  Utilizes time-series fusion of cloud-, temperature- (LST), and shadow-corrected NDSI 
#  to generate glacial debris gain over specified period
# Uses:
#    - Landsat OLI, OLI2, TM
#

# DEPENDENCIES
library(raster);library(sp);library(tools);library(stringr);library(dplyr);library(sf)     
# library(insol) # R 4.0.3 is required to use the doshade function from the insol package (called below)
#
# "R version 4.0.3 (2020-10-10)"
# platform       x86_64-pc-linux-gnu         
# arch           x86_64                      
# os             linux-gnu                   
# system         x86_64, linux-gnu           
# major          4                           
# minor          0.3                         
# year           2020                        
# month          10                          
# day            10                          
# svn rev        79318                       
# language       R                           
# version.string R version 4.0.3 (2020-10-10)
# nickname       Bunny-Wunnies Freak Out 
#
# GEOS 3.10.2, GDAL 3.4.1, PROJ 8.2.1

# other packages `fasterize`; 


# what this does:
# - looks for downloaded tile information
# - existing shapefiles -> creates masks & creates ROI 
# - reads in DEM

# replace shapefiles w/ RGI 7.0 for second round

# # # # # # # # # # # # # #
# DATAPATH LOADER FUNCTION
path_tile_manager <- function(tile_path='140041', gpoly_type="SAMRGI"){
  #
  # regional information & path locations
  # gpoly_type: "SAMRGI" or "glacierPolygons"
  #
  # !UPDATE >>>>> NEW REGIONS 
  region_names <- c("Bhutan", "Everest", "Annapurna", "Gangotri", "Himachal Pradesh", "Karakoram")
  tile_paths <-  c("138041", "140041", "142040", "145039", "148037", "148035")

  # TEST
  if(length(tile_path)>1 | !tile_path%in%tile_paths){stop("incorrect tile path specified")}
  print(paste("Loading file paths and ancillary data for:", region_names[match(tile_path, tile_paths)],
              "| tile no.", tile_path), quote=FALSE)
  #
  # PARENT DIR
  project_home_path <- '~/olson/glacier/hma-debris'
  parent_dir <<- file.path(project_home_path, 'Landsat/Landsat08', tile_path)
  parent_dirTM <<- file.path(project_home_path, 'Landsat/LTM', tile_path)
  if(!dir.exists(parent_dir)){stop('parent_dir does not exist')}else{print(parent_dir, quote=FALSE)}
  #

  # !UPDATE >>>>> NEW RGI 7.0 
  # Glacier polygons - default SAMRGI
  g_poly_path <- file.path(project_home_path, 'shp/gpoly', tile_path)
  g_path <- list.files(g_poly_path, pattern = 'SAMRGI', full.names = T)
  if(length(g_path)<1){
    print('Creating glacier polygons for region', quote=F)
    # clip polygons for region
    shp_paths_old <- '~/molson/DebrisCover/data/shp'
    # rgi_path = file.path(shp_paths_old, "RGI60/rgi60_HMA.shp")
    rgi_path = file.path(shp_paths_old, 'HerreidPelliccotti/Snew/SamRGI2km2.shp')
    rgi0 <- raster::shapefile(rgi_path)
    r <- raster(list.files(list.files(parent_dir, full.names = T, pattern ='.*T[1-2]$')[1], pattern = 'SR_B4', full.names = T))
    p <- as(raster::extent( r ), 'SpatialPolygons'); raster::crs(p) <- raster::crs( r )
    bf <- 2e4
    shproj <- spTransform(raster::buffer(p ,width=bf),crs(rgi0))
    gpoly <<-  rgi0%>%rgeos::gBuffer(byid=TRUE, width=0)%>%crop(shproj)%>%spTransform(crs(r)) %>% sf::st_as_sf()
    # rgi1 <- rgi1[rgi1$Area > area_limit,] # apply area limit
    # saveRDS(gpoly, file.path(dirname(g_poly_path) ,paste0(tile_path, "_glacierPolygons.rds")) )
    saveRDS(gpoly, file.path(g_poly_path ,paste0(tile_path, "__Herreid2km_SAMRGI.rds")) )
  }else{gpoly <<- readRDS(list.files(g_poly_path, pattern = gpoly_type, full.names = T))}
  #
  # Region of interest
  region_path <- list.files(g_poly_path, pattern = '_roi.rds', full.names = T)
  if(length(region_path)<1){
    print('Creating ROI', quote=F)
    # too slow
    # toi <- raster::rasterToPolygons(!is.na(crop(r, gpoly)), dissolve=T)
    tile_roi <<- as(raster::extent( gpoly ), 'SpatialPolygons'); raster::crs(tile_roi) <- raster::crs( gpoly )
    saveRDS(tile_roi, file.path(g_poly_path ,paste0(tile_path, '_roi.rds')) )
  }else{tile_roi <<- readRDS(list.files(g_poly_path, pattern = '_roi.rds', full.names = T))}
  #
  # Glacier mask
  gmask_path <- list.files(g_poly_path, pattern = 'gMask', full.names = T)
  if(length(gmask_path)<1){
    print('Creating glacier mask for region', quote=F)
    bf <- 2e4 # buffer radius
    r <- raster(list.files(list.files(parent_dir, full.names = T, pattern ='.*T[1-2]$')[1], pattern = 'SR_B4', full.names = T))
    rroi <- crop(r, tile_roi)
    sf::st_as_sf(gpoly) %>% 
      fasterize::fasterize(.,rroi) %>% raster::crop(.,rroi) %>% raster::mask(.,rroi) %>% 
      saveRDS(., file.path(g_poly_path ,paste0(tile_path, "_gMask.rds")))
    gmask <<- readRDS(list.files(g_poly_path, pattern = "_gMask.rds", full.names = T))
  }else{gmask <<- readRDS(list.files(g_poly_path, pattern = "_gMask.rds", full.names = T))}
  
  # DEM - terrain needed
  # dem_type="NASADEM"
  dem_path <- list.files(file.path(project_home_path, 'dems', tile_path), pattern = '.grd', full.names = T)
  dem <<- raster::raster(dem_path)
  
  # other?
  tile_path <<- tile_path
  print('Loaded: \'tile_path\', \'parent_dir\', \'gpoly\', \'tile_roi\', \'gmask\', \'dem\' into memory', quote=F)
}
