require(raster)
require(grid)
require(rgdal)

#####Climate Data#####
climate.path <- "/Users/Simon/Studium/MSC/Masterarbeit/data/climate/BOM/"

#prec newhaven, temp Yuendumu
prec <- read.csv(paste0(climate.path, "IDCJAC0009_015611_1800_Data.csv"))
temp <- read.csv(paste0(climate.path, "IDCJAC0010_015528_1800_Data.csv"))
enso <- read.csv(paste0(climate.path, "southern_oscillation_index_1960_2015.csv"))

#aggregate prec data by year (annual rainfall sum)
prec.stats <- aggregate(prec$Rainfall.amount..millimetres., by = list(prec$Year), FUN = sum, na.rm = T)
temp.stats <- aggregate(temp$Maximum.temperature..Degree.C., by = list(temp$Year), FUN =  function(x){length(na.omit(x[x > 41]))})
enso.stats <- data.frame("year" = enso$Year, "mean_enso_index" = rowMeans(enso[,2:13]))
clim.stats <- merge(prec.stats, enso.stats, by.x = "Group.1", by.y = "year")
clim.stats <- merge(clim.stats, temp.stats, by.x = "Group.1", by.y = "Group.1")

#Calculate cumulative prec
clim.stats$yycumu <- NA
for (i in 3:51){
  clim.stats$yycumu[i] <-  clim.stats$x.x[i] + clim.stats$x.x[i-1] + clim.stats$x.x[i-2]
}
names(clim.stats) <- c("year", "prec", "enso_ind", "temp", "yycumu")

#####Fire Data#####
path.shp <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/fire_data_reprojected_1970-2015"
files <- list.files(path.shp, pattern = "*.shp$")

# read fire data
fs.shp.repr <- list()
for (i in 1:length(files)){
  temp <- readOGR(path.shp, substr(files[i], 1, 9))
  names(temp)[1] <- c("OBJECTID")
  fs.shp.repr[[i]] <- temp
  }

#Get fire sizes, calender year and number of fires from fire scar shapefiles
all <- data.frame("area" = fs.shp.repr[[1]]@data$Shape_Area, "year" = fs.shp.repr[[1]]@data$CalanderYe, "no_of_fires" = length(fs.shp.repr[[1]]@data$OBJECTID))
for (i in 2:length(files)){
  allnext <- data.frame("area" = fs.shp.repr[[i]]@data$Shape_Area, "year" = fs.shp.repr[[i]]@data$CalanderYe, "no_of_fires" = length(fs.shp.repr[[i]]@data$OBJECTID))
  all <- rbind(all, allnext)
  print(i)
}

#attach prec and temp data
pos <- match(all$year, clim.stats$year)
all$cumu_prec <- clim.stats$yycumu[pos]
all$temp <- clim.stats$temp[pos]

######################################################
###PREC VS TOTAL AREA BURNED
######################################################
rm(year, temp)
detach(clim.stats)
detach(all)
attach(all)
wide_area <- 5503.63
total <- round((tapply(area, year, function(x) sum(x)))/wide_area * 100, digits = 3)

cumu <- tapply(cumu_prec, year, function(x) unique(x))

#model
mod_total <- lm(log(total) ~ cumu)
summary(mod_total)
df <- data.frame(total, cumu)

n.cumu <- seq(0,2300, length.out = 50)
preds <- predict(mod_total, newdata = data.frame("cumu" = n.cumu), se = T)

n.df <- data.frame("cumu" = n.cumu, "fit" = exp(preds$fit), "se1" = exp(preds$fit + preds$se.fit), "se2" = exp(preds$fit - preds$se.fit))

#plot of data points, fit line and standard errors
  ggplot(df, aes(cumu, total)) +
    geom_point() +
  labs(y = "Total Area Burned [km^2]", x = "Cumulative Rain") +
  geom_line(aes(cumu, fit), n.df, linetype = 1) +
  geom_line(aes(cumu, se1), n.df, linetype = 2) +
  geom_line(aes(cumu, se2), n.df, linetype = 2)

plot(residuals(mod_total))
  
######################################################
####PREC VS FIRE SIZES####
######################################################
detach(all)
detach(clim.stats)
rm(year, temp)
attach(all)

#use years with a lot of data
#data points per year
fn <- tapply(area, year, function (x) length(x))

area.subs <- area[which(year%in%c(1976, 1982,1983, 1984, 2000:2012))]
years.subs <- year[which(year%in%c(1976, 1982,1983, 1984, 2000:2012))]

#extract person moments
library(PearsonDS)
MM <-tapply(area.subs, years.subs, function(x) empMoments(x))
dpearson(x, moments = empMoments(x))

#Analyse pearson moments for correlation
mean <- vector(); variance <- vector() ;skewness <- vector(); kurtosis <- vector()
for (i in 1:length(MM)){
  mean[i] <- MM[[i]][[1]]
  variance[i] <- MM[[i]][[2]]
  skewness[i] <- MM[[i]][[3]]
  kurtosis[i] <- MM[[i]][[4]]
}

#isolate years with data
cumu.subs <- unique(cumu_prec[which(year%in%c(1976, 1982,1983, 1984, 2000:2012))])

#loglik and pars for Pearson functions
loglikpears <- function(x){
  temp <- pearsonMSC(x)
  max(temp$logLik)
}

pearsFUN <- function(x){
  temp <- pearsonMSC(x)
  temp$Best$AIC
}

pearsfuncs<- tapply(area.subs, years.subs, function(x) pearsFUN(x))

lnscales <- tapply(area.subs, years.subs, function(x) lnscale(x))
lnmeans <- tapply(area.subs, years.subs, function(x) lnmean(x))
total <- tapply(area.subs, years.subs, function(x) sum(x))
#unfitted lnorm parameters
lnscale <- function(x){
  variance <- (sd(x))^2
  mean <- mean(x)
  scale <- sqrt(log(1+ variance/mean^2))
  return(scale)
}

lnmean <- function(x){
  variance <- (sd(x))^2
  mean <- mean(x)
  scale <- sqrt(log(1+ variance/mean^2))
  loc <- log(mean) - (0.5 * scale^2)
  return(loc)
}
summary(lm(lnmeans ~ total))
plot(lnmeans ~ total)
  
  fitdist(x, "lnorm" , start = list(loc, scale))
}
x <- area
curve(dlnorm(x, loc[1], scale[1]))
?fitdist

library(ggplot2)
library(fitdistrplus)
library(MASS)

#calculate loglikelihood for each dist
Lnorm <- tapply(area.subs, years.subs, function(x) logLik(fitdist(x, "lnorm")))
Lgamma <- tapply(area.subs, years.subs, function(x) logLik(fitdist(x, "gamma", start=list(shape = 1, rate = 0.1),lower=0.1)))
Lpears<- tapply(area.subs, years.subs, function(x) loglikpears(x))


Lnorm2 <- tapply(area.subs, years.subs, function(x) logLik(fitdist(x, "lnorm")))

#Pearson moments
df <- data.frame(cumu.subs, mean, variance, skewness, kurtosis, lmean = log(mean), lvariance = log(variance))

a <- ggplot(df, aes(x = cumu.subs, y = mean)) + 
  geom_point(size = 1) +
  labs(x="Cumulative precipitation", y="Mean")

b <- ggplot(df, aes(x = cumu.subs, y =variance)) + 
  geom_point(size = 1) +
  labs(x="Cumulative precipitation", y="Variance")

c <- ggplot(df, aes(x = cumu.subs, y =skewness)) + 
  geom_point(size = 1) +
  labs(x="Cumulative precipitation", y="Skewness")

d <- ggplot(df, aes(x = cumu.subs, y =kurtosis)) + 
  geom_point(size = 1) +
  labs(x="Cumulative precipitation", y="Kurtosis")

grid.newpage()
pushViewport(viewport(layout = grid.layout(2,2)))
loc <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(a, vp = loc(1,1))
print(b, vp = loc(1,2))
print(c, vp = loc(2,1))
print(d, vp = loc(2,2))

df <- data.frame(unique(years.subs), Lnorm,Lpears,  Lgamma, mean, variance, skewness, kurtosis, cumu)
ggplot(data = df) +
  geom_point(aes(x = df[,1], y = Lnorm, col = "lognorm")) +
  geom_point(aes(df[,1], y = Lgamma, col = "gamma")) +
  geom_point(aes(df[,1], y = Lpears, col = "pearson")) +
  theme(legend.title=element_blank()) +
  labs(x = "year", y = "LogLik")

#it seems like fitted preason or fitted log normal are best. Since fitted lognormal means that the function parameters are likelihood maximised, we can't calculate parameters from empiric mean and sd and have to correlate MLE with prec

lnormfuncs <- tapply(area.subs, years.subs, function(x) fitdistr(x, "lognormal"))
meanlogs <- vector()
sdlogs <- vector()
for (i in 1:length(lnormfuncs)){
  meanlogs[i] <- lnormfuncs[[i]]$estimate[1]
  sdlogs[i] <- lnormfuncs[[i]]$estimate[2]
}

plot(meanlogs,sdlogs)
mean(meanlogs)
summary(lm(meanlogs~cumu))
plot(cumu, meanlogs)
plot(cumu, sdlogs)

require(lattice)
densityplot(~log(area),groups = year,
            plot.points = FALSE, ref = TRUE, 
            auto.key = list(space = "right"))


#plot qq plots of fitted gamma functions
require(grid)
dat <- tapply(area.subs, years.subs, function(x) list(x))
params <- tapply(area.subs, years.subs, function(x) fitdist(x, "gamma", method = "mme"))
xloc <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)
yloc <- c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6)

grid.newpage()
pushViewport(viewport(layout = grid.layout(6,3)))

for (i in 1:length(dat)){
  df <- data.frame("y" = dat[[i]])
  x <- quantile(params[[i]], c(0.25, 0.75))
  x <- unlist(x$quantiles)
  y <- quantile(df$y, c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y - slope * x
  a <- ggplot(df, aes(sample = y)) +
    stat_qq(distribution = qgamma, dparams = as.list(params[[i]]$estimate)) +
    geom_abline(slope = slope, intercept = int, lty = 2, col = "red") +
    labs(x = names(params[i]), y = "")
  print(a, vp = viewport(layout.pos.row = yloc[i], layout.pos.col = xloc[i]))
}

#plot qq plots of fitted lnorm functions
require(grid)
dat <- tapply(area.subs, years.subs, function(x) list(x))
params <- tapply(area.subs, years.subs, function(x) fitdist(x, "lnorm"))
grid.newpage()
pushViewport(viewport(layout = grid.layout(6,3)))
for (i in 1:length(dat)){
  df <- data.frame("y" = dat[[i]])
  x <- quantile(params[[i]], c(0.25, 0.75))
  x <- unlist(x$quantiles)
  y <- quantile(df$y, c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y - slope * x
  a <- ggplot(df, aes(sample = y)) +
    stat_qq(distribution = qlnorm, dparams = as.list(params[[i]]$estimate)) +
    geom_abline(slope = slope, intercept = int, lty = 2, col = "red") +
    labs(x = names(params[i]), y = "")
  print(a, vp = viewport(layout.pos.row = yloc[i], layout.pos.col = xloc[i]))
}

dev.off()
require(grid)
dev.off()


