require(raster)
require(rgdal)

#shp, convert to fire scar rasters and output
path.shp <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/fire raw data 1970-2015"
path.fire <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire"
files <- list.files(path.shp, pattern = "*.shp$")
path.fire.rasters <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/fire rasterized 1970-2015"

n <- paste0("y", substr(files,1,4))
fire.scar.ts <- stack()

for (i in 1:length(files)){
  shp <- readOGR(path.shp, substr(files[i],1,4))
  crs(shp) <- "+proj=longlat +ellps=GRS80 +no_defs"
  r.f <- rasterize(shp, maskrast, background = NA)
  r.f[which(r.f[]>0)] <- 1
  r.f[is.na(r.f[])] <- 0
  r.f[is.na(maskrast[])] <- NA
  r.f@data@names <- paste0("fs_", n[i])
#  writeRaster(r.f, paste0(path.fire.rasters, "/fs_", n[i]), format = "ascii", overwrite = T)
  fire.scar.ts <- stack(fire.scar.ts, r.f)
}
#calc tsf for each year, starting 1970, output
tsf.ts <- stack()
r.tsf <- fire.scar.ts[[1]]
r.tsf[which(r.tsf[] == 1)] <- 0 #first time step set everything to 0
path.tsf <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/tsf rasterized 1970-2015"
for (i in 1:(length(files)-1)){
  y1 <- as.numeric(substr(files[i],1,4))
  y2 <- as.numeric(substr(files[i+1],1,4))
  r.tsf[] <- r.tsf[] + (y2-y1)
  r.tsf@data@names <- paste0("tsf_", n[i])
  writeRaster(r.tsf, paste0(path.tsf, "/tsf_", n[i+1]), format = "ascii", overwrite = T)
  tsf.ts <- stack(tsf.ts, r.tsf)
  r.tsf[which(fire.scar.ts[[i]][] == 1)] <- 0
}

##################################################################################
##################################################################################

#find cells that have never burnt and exclude from analysis
fire.scars <- numeric()
for (i in 15:nlayers(fire.scar.ts)){
  fire.scars <- cbind(fire.scars, fire.scar.ts[[i]][])
}

lc <- numeric()
for (i in 15:nlayers(fire.scar.ts)){
  lc <- c(lc, veg.rast[which(burned >0)])
}
length(lc)
burned <- rowSums(fire.scars)
fs <- numeric()
for (i in 15:nlayers(fire.scar.ts)){
fs <- c(fs, fire.scar.ts[[i]][which(burned > 0)])
}
length(tsf.v)
#extract tsf.v data frmo layers
tsf.v <- numeric()
for (i in 14:nlayers(tsf.ts)){
  tsf.v <- c(tsf.v, tsf.ts[[i]][which(burned > 0)])
}
tsf.ts[[14:28]]
prec.cumu <- numeric()
for(i in c(38:51,53)){
  cumu <- prec.stats$cumu[i]
  prec.cumu <- c(prec.cumu, rep(cumu, length(tsf.ts[[1]][]))[which(burned > 0)])
}


#make model
fire.data <- data.frame("lc" = lc, "fs" = fs, "tsf.v" = tsf.v, "cumu" = prec.cumu)
fit <- glm(fs ~ tsf.v + cumu + lc, data = na.omit(fire.data), family = "binomial")

summary(prec.cumu[which(fs == 0)])

hist(tsf.)
plot(tsf.v, fs)
# Logistic Regression
# where F is a binary factor and 
# x1-x3 are continuous predictors 
summary(fit) # display results
confint(fit) # 95% CI for the coefficients
exp(coef(fit)) # exponentiated coefficients
exp(confint(fit)) # 95% CI for exponentiated coefficients
preds <- plogispredict(fit, newdata = data.frame("tsf.v" = c(1:200), "cumu" = seq(1, 1000, by = 5)), type="response") # predicted values
resids <- residuals(fit, type="deviance") # residuals
hist(resids)

#effect plots
dev.off()
#tsf.v
#plot(fire.data$fs ~ fire.data$tsf.v)
plot(seq(0,1, len = 10) ~ seq(min(fire.data$tsf.v), max(fire.data$tsf.v), len = 10), type = "n")
new.data.tsf.v <- data.frame("cumu" = mean(fire.data$cumu), "tsf.v" = seq(min(fire.data$tsf.v), max(fire.data$tsf.v), len = 100), "lc" = median(lc))
preds.tsf.v <- predict(fit, newdata = new.data.tsf.v,type = "response", se.fit=T)
lines(new.data.tsf.v$tsf.v, preds.tsf.v$fit)

#cumu
#plot(fire.data$fs ~ fire.data$cumu)
plot(seq(0,1, len = 10) ~ seq(min(fire.data$cumu), max(fire.data$cumu), len = 10), type = "n")
new.data.cumu <- data.frame("tsf.v" = mean(fire.data$tsf.v), "cumu" = seq(min(fire.data$cumu), max(fire.data$cumu), len = 100), "lc" = median(lc))
preds.cumu <- predict(fit, newdata = new.data.cumu, type="response", se.fit=T)
lines(new.data.cumu$cumu,preds.cumu$fit)


hist(preds.tsf.v$fit)
length(new.data$tsf.v)
length(preds.tsf.v$fit)
?plogis
hist(resids)
dev.off()
summary(fire.data$tsf.v)

# #landcover data subsets
# path.dlcd.input <- "/Users/Simon/Studium/MSC/Masterarbeit/data/land cover/DLCDV2"
# path.dlcd.output <- "/Users/Simon/Studium/MSC/Masterarbeit/data/land cover/subsets"
# dlcd.files <- list.files(path.dlcd.input, pattern = "*.tif$")
# temp <- gregexpr("[0-9]+",dlcd.files)
# n <- paste0("lc_", substr(sapply(regmatches(dlcd.files, temp), function (x) x[3]), 1, 4))
# dlcd.ts <- stack()
# for (i in 1:length(dlcd.files)){
#   r <- raster(paste0(path.dlcd.input, "/", dlcd.files[i]))
#   r <- projectRaster(r, mask, method = "ngb")
#   r <- crop(r,mask)
#   writeRaster(r, paste0(path.dlcd.output, "/", n[i]), format = "ascii", overwrite = T)
#   dlcd.ts <- stack(dlcd.ts,r)
# }
# 
# names(dlcd.ts) <- n
# 
# dlcd <- as.data.frame(rasterToPoints(dlcd.ts))
# #change data structure, extract 2002-2009
# dlcd <- dlcd[sort(names(dlcd))]
# dlcd <- dlcd[,c(4:5)]
# head(dlcd)
# dlcd.stack <- stack(dlcd)
# names(dlcd.stack) <- c("class", "year")
# classes <- unique(dlcd.stack$class)
# for (i in 1:length(classes)){
#   cl.n <- paste0("cls", classes[i])
#   cl <- (rep(0, length(dlcd.stack[,1])))
#   cl[which(dlcd.stack$class == classes[i])] <- 1
#   dlcd.stack <- cbind(dlcd.stack, cl)
#   names(dlcd.stack)[i+2] <- cl.n
# }
# head(dlcd.stack)
# tsf <- as.data.frame(rasterToPoints(tsf.ts[[20:21]]))[,-c(1:2)]
# tsf <- stack(tsf)
# names(tsf) <- c("tsf", "year")
# fs <- as.data.frame(rasterToPoints(fire.scar.ts[[20:21]]))[,-c(1:2)]
# fs <- stack(fs)
# names(fs) <- c("fs","year")
# 
# raw.data <- data.frame(dlcd.stack, tsf, fs)
# hist(log(raw.data$tsf))
# 
# dev.off()
# summary(glm)
# AIC(glm)
# AIC(glm2)
# summary(glm3)
# gam <- gam(fs ~ (cls19 + cls16 + cls24 + cls25 + cls4 + cls33 + cls34 + cls14 + cls18 + cls11 + tsf), data = raw.data.scaled)
# 
# #create new data to predict to
# dlcd.new <- as.data.frame(rasterToPoints(dlcd.ts))
# 
# #change data structure, extract year to simulate
# dlcd.new <- dlcd.new[sort(names(dlcd.new))]
# dlcd.new <-dlcd.new[,10]
# class <- dlcd.new
# classes.new <- unique(class)
# for (i in 1:length(classes)){
#   cl.n <- paste0("cls", classes[i])
#   cl <- (rep(0, length(class)))
#   cl[which(dlcd.new == classes.new[i])] <- 1
#   dlcd.new <- data.frame(dlcd.new, cl)
#   names(dlcd.new)[i+1] <- cl.n
# }
# 
# tsf.new <- as.data.frame(rasterToPoints(tsf.ts[[26]]))[3]
# tsf.new <- stack(tsf.new)
# names(tsf.new) <- c("tsf", "year")
# fs.new <- as.data.frame(rasterToPoints(fire.scar.ts[[26]]))[3]
# fs.new <- stack(fs.new)
# names(fs.new) <- c("fs","year")
# 
# raw.data.new <- data.frame(tsf.new, fs.new, dlcd.new)
