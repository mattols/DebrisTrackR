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

# INFORMATION about required libraries and dependencies see `path_manager.R`

##################
# MAIN CALL SCRIPT

# # # # # # # # # # #
# LOAD CUSTOM SCRIPTS
files <- list.files(
  path = "scripts",
  pattern = "\\.R$",
  full.names = TRUE
)
files <- files[!basename(files) %in% "run_main.R"]
lapply(files, source)

# SAVE DIR
## >>> UPDATE!
dp0 <- '~/olson/glacier/hma-debris/output/new-debris-class-X/'
if(!dir.exists(dp0)){dir.create(dp0)}

# VALID TILES
tile_paths00 <- c("138041", "140041", "142040", "145039", "148037", "148035")

# SPECIFY TILE!
tile_01 = tile_paths00[4]


### NEEDED?
# MINMAL
path_tile_manager(tile_01)

# CHANGE WORKFLOW
# >>>>>>>

# ITERATIVE CHANGE
strt_ranges = list(1993:1995, 2003:2005, 2013:2015, 2023:2024)
var_cache <- c(ls(), lsf.str(), "var_cache")
for(ti in 1:(length(strt_ranges)-1)){
  glacier_debris_fusion(tile_region = tile_01, T1=TRUE,
                        time1=strt_ranges[[ti]], 
                        time2=strt_ranges[[ti+1]], months=7:10, 
                        region = tile_roi, 
                        dest_path = dp0)
  rm(list = setdiff(ls(), var_cache));gc()
}

# FULL TIME REGION
glacier_debris_fusion(tile_region = tile_01, T1=TRUE,
                      time1=strt_ranges[[1]], 
                      time2=strt_ranges[[length(strt_ranges)]], months=7:10,
                      region = tile_roi, 
                      dest_path = dp0)
