require(raster)
path.fire <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire"
path.fire.rasters <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/fire rasterized 1970-2015"

files <- list.files(path.fire.rasters, pattern = "*.asc$")
mask <- raster(paste0(path.fire, "/mask.asc"), format = "ascii")

#load fire scar layers
fire.scar.ts <- stack()
for (i in 1:length(files)){
  r <- raster(paste0(path.fire.rasters, "/",files[i]))
  r[which(is.na(mask[]))] <- NA
  fire.scar.ts <- stack(fire.scar.ts,r)
}

#make data frame from fs data
fs.df <- data.frame(fire.scar.ts[])
# #write mfn per land type into raster
# mfn.lt <- lt.rast
# for (i in 1:length(mfn.per.lt[,1])){
#   mfn.lt[which(mfn.lt[] == mfn.per.lt[i,1])] <- mfn.per.lt[i,2]
# }
