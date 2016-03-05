require(msm)
require(raster)
require(rgdal)

##########################
#I. MODEL DEFINITIONS#####
##########################

#1.) Global parameters

#1.a) Paths and file lists
path_climate <- "/Users/Simon/Studium/MSC/Masterarbeit/data/climate/BOM/"
path_shp <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/fire_data_reprojected_1970-2015"
path_fire_rasters <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/fire rasterized 1970-2015"
path_tsf_rasters <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/tsf rasterized 1970-2015"
files_shp <- list.files(path_shp, pattern = "*.shp$")
files_tsf <- list.files(path_tsf_rasters, pattern = "*.asc$")
files_fs <- list.files(path_fire_rasters, pattern = "*.asc$")
output <- "/Users/Simon/Studium/MSC/Masterarbeit/data/new_model_out2" #output path for fire record and rasters

#1.b) Simulation parameters
scenario <- c("low", "high", "neut") #scenario
simle <- 75 #simulation length (cannot be changed!)
simno <- 20 #number of simulations
ini_year <- 2016

#1.c) Spatial constants
cell_size <- raster("/Users/Simon/Studium/MSC/Masterarbeit/data/cell_size.asc")
elev <- raster(paste0("/Users/Simon/Studium/MSC/Masterarbeit/data//Elevation/elevation9secNH.asc")) #elev raster
cell_size <- cell_size[] #cell sizes in km^2
elev <- elev[] #elev in m
maskrast <- raster("/Users/Simon/Studium/MSC/Masterarbeit/data/maskraster.asc")
study_a <- maskrast[] #study area
inds <- 1:length(maskrast[]) #cell indices
ncol <- ncol(maskrast) #raster dim
nrow <- nrow(maskrast) #raster dim
