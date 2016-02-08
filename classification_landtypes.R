require("raster")
require("rgdal")
path.fire <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire"
veg.path <- "/Users/Simon/Studium/MSC/Masterarbeit/data/land types/Flora/"
geo.path <- "/Users/Simon/Studium/MSC/Masterarbeit/data/land types/Landform/"
veg.map <- readOGR(veg.path, "vegmap_v3_all_veg_types")
geo.map <- readOGR(geo.path, "Geology")
mask <- raster(paste0(path.fire, "/mask.asc"), format = "ascii")
#convert to longlat
#veg.map <- spTransform(veg.map, CRS = "+proj=longlat +ellps=GRS80 +no_defs")
#geo.map <- spTransform(geo.map, CRS = "+proj=longlat +ellps=GRS80 +no_defs")

#rasterize
veg.rast <- rasterize(veg.map, field = veg.map@data$CLASS, mask)
geo.rast <- rasterize(geo.map, field = as.integer(geo.map@data$UNITNAME), mask)
plot(veg.rast)
veg.leg
#set geo raster NA where veg is not available
geo.rast[which(is.na(mask[]))] <- NA

#correct map ( 10 BW and 7 MUW)
veg.map@data$CLASS[which(veg.map@data$CLASS == 10 & veg.map@data$COMMENTS == "BW")] <- 5
veg.map@data$CLASS[which(veg.map@data$CLASS == 7 & veg.map@data$COMMENTS == "MUW")] <- 17
  
#extract legs
#vegetation legend
veg.temp <- data.frame("ID" = veg.map@data$CLASS, "NAME" = veg.map@data$COMMENTS)
veg.temp <- veg.temp[order(veg.temp$ID),]
veg.leg <- unique(veg.temp)

#geoclasses legend
geo.leg <- unique(geo.map@data$UNITNAME)

geo.leg <- data.frame("ID" = as.integer(geo.leg), "NAME" = as.character(geo.leg))
geo.leg <- geo.leg[order(geo.leg$ID),]

geo.df <- as.data.frame(geo.rast)
veg.df <- as.data.frame(veg.rast)

#combine rasters into new land type raster
lt <- geo.df[,1]*100 + veg.df[,1]

#create land type legend
lt.leg <- unique(lt)
lt.leg <- lt.leg[order(lt.leg)]
lt.leg <- lt.leg[-which(is.na(lt.leg))]
geo.classes <- round(lt.leg/100)
veg.classes <- lt.leg-round(lt.leg, digits = -2) 
ID <- c(1:length(lt.leg))
lt.leg <- data.frame(lt.leg, geo.classes, veg.classes, ID)


for(i in 1:length(lt.leg[,1])){
  lt[which(lt == lt.leg$lt.leg[i])] <- lt.leg$ID[i]
}

lt.rast <- mask
lt.rast[] <- lt

#output new land type raster and legend files
path.land.types <- "/Users/Simon/Studium/MSC/Masterarbeit/data/land types"
writeRaster(lt.rast, paste0(path.land.types, "/landtypes.asc"), format = "ascii", overwrite = T)


#spinifex types and area they cover
veg.cells <- veg.rast[which(!is.na(veg.rast[]))]
perc.ecosyst <- numeric()
areas <- aggregate(veg.map@data$NEWAREA, by = list(veg.map@data$CLASS),  FUN = "sum")
for (i in 1:23){
veg.leg$PERC[i] <- (areas$x[i]/colSums(areas)[2]) * 100
}
#matching long names and abbirvations from file into legetnd
veg.names <- read.csv("/Users/Simon/Studium/MSC/Masterarbeit/data/Veg classifications (easy to read).csv", header = F)
pos <- match(veg.leg$NAME, veg.names$V1)
for (i in 1:23){
  veg.leg$LONG[i] <- as.character(veg.names$V2[pos[i]])
}

ssp.cover <- sum(veg.leg$PERC[which(grep("spinifex", veg.leg$LONG) > 0)])

veg.map@data$COL <- NA
cols <- terrain.colors(23)
for (i in 1:23){
veg.map@data$COL[which(veg.map@data$CLASS == i)] <- cols[i]
}

spinifex <- veg.map[veg.map@data$CLASS%in%c(3,8,9,16),]
plot(spinifex, col = spinifex@data$COL)
spinifex@data

write.csv(lt.leg, paste0(path.land.types, "/legend_land_types.csv"))
write.csv(geo.leg, paste0(path.land.types, "/legend_geo.csv"))
write.csv(veg.leg, paste0(path.land.types, "/legend_veg.csv"))