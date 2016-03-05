require(msm)
require(raster)
require(rgdal)
require(mixtools)
require(MASS)
#########################
#II. Model estimation####
#########################

#1.) Import data

#1.a) Fire scar
fire.scar.ts <- stack()
old <- 1970
for (i in 1:length(files_fs)){
  r <- raster(paste(path_fire_rasters, files_fs[i], sep ="/"))
  new <- as.numeric(substr(names(r@data), 5, 8))
  while(new-old > 1){
    s <- maskrast - 1
    names(s) <- paste0("fs_y", old + 1)
  crs(s) <- "+proj=longlat +ellps=GRS80 +no_defs"
  fire.scar.ts <- stack(fire.scar.ts, s)
  old <- as.numeric(substr(names(s@data), 5, 8))
  }
  crs(r) <- "+proj=longlat +ellps=GRS80 +no_defs"
  fire.scar.ts <- stack(fire.scar.ts, r)
  old <- as.numeric(substr(names(r@data), 5, 8))
}

#1.b) Time since fire (tsf)

tsf.ts <- stack()
for (i in 1:length(files_tsf)){
  r <- raster(paste(path_tsf_rasters, files_tsf[i], sep ="/"))
  crs(r) <- "+proj=longlat +ellps=GRS80 +no_defs"
  tsf.ts <- stack(tsf.ts, r)
}


#1.c) Combine fire and precipitation statistics

#i) Read fire sizes from shape files and write into data frame
wide_area <- 5503.63 #the total area of available data
all <- data.frame()
for (i in 1:length(files_shp)){
  
  #read shapefiles
  shp <- readOGR(path_shp, substr(files_shp[i], 1, 9))
  names(shp)[1] <- c("OBJECTID")
  
    #extract data
  f <- data.frame("area" = shp@data$Shape_Area, 
                  "year" = shp@data$CalanderYe, 
                  "no_of_fires" = length(shp@data$OBJECTID))
  all <- rbind(all, f)
}

#ii) read and process precipitation and enso data and attach to data frame
enso <- read.csv("/Users/Simon/Studium/MSC/Masterarbeit/data/climate/BOM/southern_oscillation_index_1960_2015.csv")
enso_means <- data.frame("enso" = rowSums(enso[,8:12])/4) #aggregate to July - November means
enso_means$year <- enso$Year

pre <- read.csv(paste0(path_climate, "IDCJAC0009_015611_1800_Data.csv"))
prec.stats <- aggregate(pre$Rainfall.amount..millimetres., by = list(pre$Year), FUN = sum, na.rm = T)

prec.stats$yycumu <- NA
prec.stats$yyycumu <- NA
for (i in 3:nrow(prec.stats)){
  prec.stats$yycumu[i] <-  prec.stats$x[i-1] + prec.stats$x[i-2]   
  prec.stats$yyycumu[i] <-  prec.stats$x[i] + prec.stats$x[i-1] + prec.stats$x[i-2] 
}
names(prec.stats) <- c("year", "prec", "yycumu", "yyycumu")
pos <- match(all$year, prec.stats$year)
pos.ens <- match(all$year, enso_means$year)
all$prec <- prec.stats$prec[pos]
all$cumuyy <- prec.stats$yycumu[pos]
all$cumuyyy <- prec.stats$yyycumu[pos]
all$enso <- enso_means$enso[pos.ens-1] #esno of previous year
#2.) Model Total Area Burned vs Precipitation
rm(enso)
attach(all)
#2.a) Aggregate Data

total <- round(tapply(area/wide_area * 100, year, sum), digits = 3) #Total area burned per year [%]
cumu2 <- tapply(cumuyy, year, function(x) unique(x)) #Cumulative Precipitation per year
cumu3 <- tapply(cumuyyy, year, function(x) unique(x)) #Cumulative Precipitation per year
precan <- tapply(prec, year, function(x) unique(x)) #Cumulative Precipitation per year
ensoan <- tapply(enso, year, function(x) unique(x)) #Cumulative Precipitation per year
#2.b) Model comparison
tab_mod1 <- lm(log(total) ~ cumu3)
tab_mod2 <- lm(sqrt(total) ~ sqrt(cumu3))
tab_mod3 <- lm(sqrt(total) ~ sqrt(cumu3) + ensoan)
tab_mod4 <- lm(sqrt(total) ~ sqrt(cumu2))
stepmod <- stepAIC(mod4)
summary(tab_mod1)
summary(tab_mod2)
summary(tab_mod3)
summary(tab_mod4)
summary(stepmod)

tab_mod <- tab_mod2

#3.) Model Fire Distribution Parameters vs Precipitation

#3.a) Aggregate and subset data (only years with sufficeitn data)
s <- tapply(area, year, FUN = length)
subyrs <- c(1976, 1982,1983, 1984, 2000:2005, 2008:2012)
area.subs <- area[which(year%in%subyrs)]
years.subs <- year[which(year%in%subyrs)]
precan.subs <- unique(prec[which(year%in%subyrs)])
cumu2.subs <- unique(cumuyy[which(year%in%subyrs)])
cumu3.subs <- unique(cumuyyy[which(year%in%subyrs)])
#3.b) Determine unfitted lnorm pars
scales <- tapply(area.subs, years.subs, function(x) lnscale(x))
locs <- tapply(area.subs, years.subs, function(x) lnmean(x))

#3.c) Model
locs_mod1 <- lm(locs ~ cumu3.subs)
locs_mod2 <- lm(locs ~ sqrt(cumu3.subs))
scales_mod1 <- lm(scales ~ cumu3.subs)
scales_mod2 <- lm(scales ~ sqrt(cumu3.subs))
summary(locs_mod1)
summary(locs_mod2)
summary(scales_mod1)
summary(scales_mod2)

plot(locs ~ cumu3.subs)
points(cumu3.subs, predict(locs_mod1), col = "red")
points(cumu3.subs, predict(locs_mod2), col = "green")
plot(scales ~ cumu3.subs)
points(cumu3.subs, predict(scales_mod1), col = "red")
points(cumu3.subs, predict(scales_mod2), col = "green")

locs_mod <- locs_mod2
scales_mod <- scales_mod1
summary(mod_probs[[5]])
#4.) Model fire occurance vs precipitation
summary(mod_occ)
#4.a) Aggregation of data
yrs <- 1970:2015
oc <- 1:length(yrs)
oc[match(unique(year), yrs)] <- 1
oc[-match(unique(year), yrs)] <- 0
pr <- prec.stats$yyycumu[which(prec.stats$year%in%yrs)]

#4.b) Model
mod_occ <- glm(oc ~ pr, family = "binomial")

#5.) Model Fire probabilities vs Time Since Fire

#5.a) Aggregation of Data

#i) Identification of pyric classes
haz_time <-as.numeric(dimnames(total)[[1]]) ##for estimation of hazard function
haz_fs <- data.frame(fire.scar.ts[])
pyric_class <- rowSums(haz_fs) #determine number of fires on each cell
pyric_class[which(pyric_class%in%c(6,7))] <- 5
total_class <- tapply(pyric_class, pyric_class, "length")[-1]

#ii) Estimation of fire probability vs tsf mean curves by pyric class
Ftl <- data.frame(matrix(0, nrow = 20, ncol = 6))
for (i in 1:10){
  Ft <- data.frame(matrix(nrow = length(haz_fs[1,i:(i + 19)]), ncol = length(total_class)))
  
  Ft$yr <- rev(2006+i - haz_time[i:(i + 19)])
  for (j in 1: length(total_class)){
    
    sub <- haz_fs[which(pyric_class == j),i:(i + 19)] #starting yr 6 because big fire in yr 5
    Ft[1,j] <- length(sub[,1][which(sub[,1] == 1)])
    
    for (k in 2:(length(sub)-1)){
      
      ch <- sub[,k+1] - rowSums(sub[,1:k])
      Ft[k,j] <- length(ch[which(ch == 1)])
    }
    Ft[,j] <- cumsum(Ft[,j]/length(sub[,1])) #build cumulative sum of p
  }
  Ftl <- Ftl + Ft
}
Ftl <- Ftl/10

#5.b.) Models

mod_probs <- list()
for (i in 1: length(total_class)){
  mod_probs[[i]] <- NA
  while(length(mod_probs[[i]]) == 1){
    py <- data.frame(cumu = Ftl[,i], yr = Ftl[,6])
    mod_probs[[i]] <- tryCatch(nls(cumu ~ c * exp(-(yr/b)^a), data = py, start = list(a = sample(-5:5, 1), b = sample(0:30,1), c = sample(0:1,1))), error = function(e) NA)
  }
}
summary(mod_probs[[5]])
#6.) Estimation of precipitation denisty distribution in reference time psan
attach(prec.stats)

mixmdl_prec <- normalmixEM(prec[which(year <= 2006 & year >= 1987)])
components <- sample(1:2,prob=c(mixmdl_prec$lambda[1], mixmdl_prec$lambda[2]),size=10000,replace=TRUE)
mus <- c(mixmdl_prec$mu[1], mixmdl_prec$mu[2])
sds <- c(mixmdl_prec$sigma[1], mixmdl_prec$sigma[2])

plot(mixmdl_prec,2)
hist(log(prec[which(year <= 2006 & year >= 1987)]), breaks = 10)

#7.) Estimation of climate scenario variables for length of time series (projected annual precipiation change)

#7.a) Estimating 87 - 2005 reference data
reference_prec <- prec[year >= 1986 & year <= 2005]
mean_prec <- mean(reference_prec)
high_prec_30 <- mean_prec *108/100
low_prec_30 <- mean_prec * 89/100
high_prec_90 <- mean_prec *119/100
low_prec_90 <- mean_prec * 69/100
neut_prec_30 <- mean_prec

#7.b) Calculation of annual increments/decrements for high and low scenario
years_30 <- 2030 - 2015
years_90 <- 2090 -2030
incr_15_30 <- (high_prec_30 - mean_prec) / years_30
decr_15_30 <- (low_prec_30 - mean_prec) / years_30
incr_30_90 <- (high_prec_90 - high_prec_30) / years_90
decr_30_90 <- (low_prec_90 - low_prec_30) / years_90

incr <- c(rep(incr_15_30, times = years_30), rep(incr_30_90, times = years_90))
decr <- c(rep(decr_15_30, times = years_30), rep(decr_30_90, times = years_90))
neut <- rep(0, 75)
summary(mod_occ)
detach(prec.stats)
detach(all)