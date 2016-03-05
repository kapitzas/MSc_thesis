require(raster)
require(rgdal)
require(fitdistrplus)

#calibration of firesize distribution function
#Extraction of total area burned
path.shp <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/fire raw data 1970-2015/"
files <- list.files(path.shp, pattern = "*.shp$")

fs.shp <- list()

# read fire data
for (i in 1:length(files)){
  fs.shp[[i]] <- readOGR(path.shp, substr(files[i], 1, 4)) #-5 if 2000-2015 data are used
}

#Get fire sizes, calender year and number of fires from fire scar shapefiles
all <- data.frame("area" = fs.shp.repr[[1]]@data$Shape_Area, "year" = fs.shp.repr[[1]]@data$CalanderYe, "no_of_fires" = length(fs.shp.repr[[1]]@data$OBJECTID))
for (i in 2:length(files)){
  allnext <- data.frame("area" = fs.shp.repr[[i]]@data$Shape_Area, "year" = fs.shp.repr[[i]]@data$CalanderYe, "no_of_fires" = length(fs.shp.repr[[i]]@data$OBJECTID))
  all <- rbind(all, allnext)
  print(i)
}

#determine annual fire stats
total.burned <- aggregate(all$Shape_Area, by = list(all$CalanderYe), FUN = "sum")
total.count <- aggregate(all$Shape_Area, by = list(all$CalanderYe), FUN = "length")
average.size <- aggregate(all$Shape_Area, by = list(all$CalanderYe), FUN = "mean")
min.size <- aggregate(all$Shape_Area, by = list(all$CalanderYe), FUN = "min")
max.size <- aggregate(all$Shape_Area, by = list(all$CalanderYe), FUN = "max")


#annual burnt distr
ab_ln <- fitdist(total.burned[,2], distr = "lnorm")

#Filter cumulative precipitation for the modelled years
pos <- match(unique(all$year), prec.stats$year) 
precs <- prec.stats$yyycumu[pos]

#colour ramp
indices <- data.frame(precs, 1:29)
indices <- indices[order(precs),]
pal <- colorRampPalette(c("lightblue", "black"))
cl <- pal(29)
names(cl) <- indices$precs

hist <- ggplot(data = all) + geom_histogram(data = all, aes(Shape_Area), bins = 200)  + xlim(0,0.01) + ylim(0,1000)
# 
# 
# for(i in indices$X1.19){
# x <- seq(fs_pl[[i]]$xmin, 0.1, length.out = 1000)
# hist <- hist + geom_line(aes(x,y, colour = cols), data = data.frame("x" = x, "y" = dplcon(x, xmin = fs_pl[[i]]$xmin, alpha = fs_pl[[i]]$pars), "cols" = names(cl)[i]))
# j = j +1
# }
# hist + scale_colour_manual(name = "Cumu. Prec. [mm]", values = cl) +   
#   theme(legend.position = c(0.9, 0.51), legend.key.size = unit(1, "cm")) +
#   xlab("Fire size [%]") + ylab("Frequency")

##analysis of relationships
#burned vs count
dev.off()
df <- data.frame("count" = total.count$x, "burned" = total.burned$x)
df <- df[order(df$count),]
lm.burntvscount <- lm(burned ~ count, data = df)
anova(lm.burntvscount)
preds.burntvscount <- predict(lm.burntvscount, se.fit = T)
plot(df$count, df$burned)
lines(df$count, preds.burntvscount$fit)
lines(df$count, preds.burntvscount$fit + preds.burntvscount$se.fit, col = 'red', lty = 2)
lines(df$count, preds.burntvscount$fit - preds.burntvscount$se.fit, col = 'red', lty = 2)

prec.stats <- aggregate(prec$Rainfall.amount..millimetres., by = list(prec$Year), FUN = sum, na.rm = T)
shapiro.test(log(prec.stats$x))
shapiro.test(log(prec.stats$cumu))


####TRENDS
dev.off()
#number of fires per year
plot(total.count$Group.1, total.count$x)
lm.count <- lm(total.count$x ~ total.count$Group.1)
preds.count <- predict(lm.count, se.fit = T)
lines(total.count$Group.1, preds.count$fit, col = 'red')
lines(total.count$Group.1, preds.count$fit + preds.count$se.fit, col = 'red', lty = 2)
lines(total.count$Group.1, preds.count$fit - preds.count$se.fit, col = 'red', lty = 2)

anova(lm.count)
hist(total.count$x, breaks = 20)

#total area burned per year
plot(total.burned$Group.1, total.burned$x)
lm.burned <- lm(total.burned$x ~ total.burned$Group.1)
preds.burned <- predict(lm.burned, se.fit = T)
lines(total.burned$Group.1, preds.burned$fit, col = 'red')
lines(total.burned$Group.1, preds.burned$fit + preds.burned$se.fit, col = 'red', lty = 2)
lines(total.burned$Group.1, preds.burned$fit - preds.burned$se.fit, col = 'red', lty = 2)
anova(lm.burned)


# fire size per year
plot(all$CalanderYe, all$Shape_Area)
lm.size <- lm(all$Shape_Area ~ all$CalanderYe)
preds.size <- predict(lm.size, se.fit = T)
lines(all$CalanderYe, preds.size$fit, col = 'red')
lines(all$CalanderYe, preds.size$fit + preds.size$se.fit, col = 'red', lty = 2)
lines(all$CalanderYe, preds.size$fit - preds.size$se.fit, col = 'red', lty = 2)
anova(lm.size)


# size <- length(na.omit(fire.probs[]))
# 
# ty <- round(rlnorm(100000, ab_ln$estimate[1], ab_ln$estimate[2])*size)
# tf <- round(rplcon(100000, xmin = fs_pl$getXmin(), alpha = fs_pl$getPars())*size)
# 
# poly.mask <- rasterToPolygons(mask, dissolve = T)
# plot(poly.mask)
# crs(poly.mask) <- "+init=epsg:3112"
# getwd()
#writeOGR(poly.mask, driver = "ESRI Shapefile", layer = "mask_newhaven", dsn = getwd())