require(rgeos)
require(rgdal)
#fire shp paths
setwd("/Users/Simon/Studium/MSC/Masterarbeit/data/fire/fire raw data 1970-2015/")
pshp12 <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/2012 fire scar shapefiles/"
pshp13 <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/2013 fire scar shapefiles/"
pshp14 <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/2014 fire scar shapefiles/"
pshp15 <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/2015 fire scar shapefiles/"



#read new fire shapefiles
y12_15 <- list()
y12_15[[1]] <- readOGR(pshp12, "fs12_mths_gda")
y12_15[[2]] <- readOGR(pshp13, "fs13_mths_gda")
y12_15[[3]]<- readOGR(pshp14, "fs2014_mths_gda")
y12_15[[4]]<- readOGR(pshp15, "fs2015shp")
fs.shp
#crop, calculate shape areas and lengths, assign cal year
y12_15cropped <- list()
y0 <- 2011
for (i in 1:4){
  x <- crop(y12_15[[i]], mask)
  if (!is.null(x)){
    names(x)[3] <- "BurnMonth"
    names(x) <- strtrim(names(x),10)
    x@data$Shape_Area <- gArea(x, byid = T)
    x@data$Shape_Length <- gLength(x, byid = T)
    x@data$CalanderYe <- y0 + i
    writeOGR(x, dsn = path.shp, layer = y0 + i, driver = "ESRI Shapefile", overwrite_layer = T) #function needs to tweaked
  }
}
require(rgeos)
setwd("/Users/Simon/Studium/MSC/Masterarbeit/data/fire/raw data/cropped/")
raw.data <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/fire raw data 1970-2015/"


#read new fire shapefiles
files <- list.files(raw.data, pattern = "*.shp$")
substr(files[1], 1, 4)

#crop, calculate shape areas and lengths, assign cal year
for (i in 1:16){
  shp <- readOGR(raw.data, substr(files[i], 1, 4))
  x <- crop(shp, mask)
  if (!is.null(x)){
    names(x) <- strtrim(names(x),10)
    x@data$Shape_Area <- gArea(x, byid = T)
    x@data$Shape_Length <- gLength(x, byid = T)
    x@data$CalanderYe <- 1999 + i
    writeOGR(x, dsn = getwd(), layer = paste0(stri_sub(files[i], 1, -5),"_cropped"), driver = "ESRI Shapefile", overwrite_layer = T) #function needs to tweaked
  }
}


#load shp data
fs.shp<- list()
for (i in 1:length(files)){
  temp <- readOGR(raw.data, substr(files[i], 1, 4))
  names(temp)[1] <- c("OBJECTID")
  fs.shp[[i]] <- temp
}

#defuckisize 1970
target <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/fire_data_reprojected_1970-2015"
crs(fs.shp[[1]]) <- crs(fs.shp[[2]])

fs.shp.repr <- list()
for (i in 1:29){
   temp <- spTransform(fs.shp[[i]], crs("+proj=utm +zone=52 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
   #temp <- crop(temp, mask)
   if (!is.null(temp)){
   temp@data$Shape_Area <- gArea(temp, byid = T) * 1E-6
   fs.shp.repr[[i]] <- temp
   writeOGR(temp, dsn = target, layer = paste0(unique(temp@data$CalanderYe), "_repr" ), driver = "ESRI Shapefile", overwrite_layer = T)
   }
}
