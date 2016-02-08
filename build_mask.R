require(rgeos)
#build a mask for area with available veg data
path.shp <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/"
files <- list.files(path.shp, pattern = "*.asc$")

template <- raster("/Users/Simon/Studium/MSC/Masterarbeit/data/Elevation/elevation9secNH.asc")
maskrast <- raster(vals = values(template), extent(template), res = res(template), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
maskrast[which(!is.na(maskrast[]))] <- 1

writeRaster(maskrast, paste0(path.fire,"maskraster"), format = "ascii", overwrite = T)

mask <- rasterToPolygons(maskrast, dissolve = T)
mask <- spTransform(mask, crs("+proj=utm +zone=52 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

writeOGR(mask, path.fire,  "mask", driver = "ESRI Shapefile", overwrite_layer = T)
cell_size <- raster::area(template, na.rm=TRUE, weights=FALSE)
writeRaster(cell_size, paste0(path.fire, "cell_size.asc"), format = "ascii", overwrite = T)
sum(na.omit(cell_size[]))
gArea(mask)/1000000

#determine approximate extended area
plot(extent(fs.shp[[26]]))
plot(fs.shp[[26]], add = T)
topright <- drawPoly()
crs(topright) <- crs(fs.shp[[26]])
leftbottom <- drawPoly()
crs(leftbottom) <- crs(fs.shp[[26]])

topright <- spTransform(topright, crs(mask))
leftbottom <- spTransform(leftbottom, crs(mask))
cutout <- gUnion(topright, leftbottom)
as.data.frame(cutout)
complete <- drawPoly()
crs(complete) <- crs(fs.shp[[26]])
complete <- spTransform(complete, crs(mask))
writeOGR(cutout, path.fire, "cutout_extended_a", driver = "ESRI Shapefile", overwrite_layer = T)
wide_area <- (gArea(complete) - gArea(cutout)) * 1E-6
#wide_area approx 5503.63 km^2