#####CLIMATE DATA PROCESSING#####

#####DATA PREP#####
climate.path <- "/Users/Simon/Studium/MSC/Masterarbeit/data/climate/BOM/"

#prec newhaven, temp Yuendumu
prec <- read.csv(paste0(climate.path, "IDCJAC0009_015611_1800_Data.csv"))
temp <- read.csv(paste0(climate.path, "IDCJAC0010_015528_1800_Data.csv"))
enso <- read.csv(paste0(climate.path, "southern_oscillation_index_1960_2015.csv"))

#aggregate prec data by year (annual rainfall sum)
prec.stats <- aggregate(prec$Rainfall.amount..millimetres., by = list(prec$Year), FUN = sum, na.rm = T)
temp.stats <- aggregate(temp$Maximum.temperature..Degree.C., by = list(temp$Year), FUN =  function(x){length(na.omit(x[x > 41]))})
prec.stats.month <- aggregate(prec$Rainfall.amount..millimetres., by = list(prec$Month), FUN = sum, na.rm = T)
enso.stats <- data.frame("year" = enso$Year, "mean_enso_index" = rowMeans(enso[,2:13]))
clim.stats <- merge(prec.stats, enso.stats, by.x = "Group.1", by.y = "year")
clim.stats <- merge(clim.stats, temp.stats, by.x = "Group.1", by.y = "Group.1")

#Calculate cumulative prec
clim.stats$yycumu <- numeric(length = 51)
clim.stats$yyycumu <- numeric(length = 51)
for (i in 3:51){
  clim.stats$yycumu[i] <-  clim.stats$x.x[i] + clim.stats$x.x[i-1] + clim.stats$x.x[i-2]
}

clim.stats$yycumuoct <- numeric(length = 51)
for(i in 4:51){
  clim.stats$yyycumu[i] <-  clim.stats$x.x[i-1] + clim.stats$x.x[i-2] + clim.stats$x.x[i-3]
}

#precipitation 24m before april
year <- vector()

for (i in 1:51){
  y <- sum(prec$Rainfall.amount..millimetres.[which(prec$Month%in%c(5:12) & prec$Year == 1962+i)], na.rm = T)
  yy <- sum(prec$Rainfall.amount..millimetres.[which(prec$Month%in%c(1:12) & prec$Year == 1962 + i + 1)], na.rm = T)
  yyy <- sum(prec$Rainfall.amount..millimetres.[which(prec$Month%in%c(1:4) & prec$Year == 1962 + i + 2)], na.rm = T)
  clim.stats$yycumuoct[i] <- y +yy + yyy
  year[i] <- 1962 + i + 2
}
clim.stats
names(clim.stats) <- c("year", "prec", "enso_ind", "temp", "yycumu", "yyycumu", "yycumuoct")

#enso vs prec
dev.off()
plot(clim.stats$year, clim.stats$enso_ind, type = "l")
par(new = T)
plot(clim.stats$year, clim.stats$prec, type = "l", col = "red")
par(new = T)
plot(clim.stats$year, clim.stats$yyycumu, type = "l", col = "green")
axis(4)
plot(clim.stats$enso_ind ~clim.stats$prec)
mod <- glm(clim.stats$prec ~clim.stats$enso_ind)
summary(mod)

fire.stats <- cbind(average.size, total.burned$x, total.count$x, max.size$x, min.size$x)
names(fire.stats) <- c("year", "average.size", "total.burned", "total.count", "max.size", "min.size")

fire_clim <- merge(fire.stats, clim.stats, by = "year", all.y = T)

#add data point labels
fire_clim$lab <- substr(fire_clim$year, 3,4)

write.csv(fire_clim, file = paste0(climate.path, "fire_clim_newhaven_1965-2015.csv"))

#prec vs total burned [significant at 5%]
fire_clim_subs <- na.omit(fire_clim[-which(fire_clim$year == 1975),])
mod <- glm(log(total.burned) ~ yycumu, data = fire_clim)
mod2 <- glm(log(total.burned) ~ yycumu, data = fire_clim_subs)

require(ggrepel)
require(grid)
require(ggplot2)

#calculate predictors and make new data frame with preds and SEs
preds <- predict(mod, newdata = data.frame("yycumu" = seq(0,2300, length.out = 50)), se = T)
preds2 <- predict(mod2, newdata = data.frame("yycumu" = seq(0,2300, length.out = 50)), se = T)
newdata <- data.frame("yycumu" = seq(0, 2300, length.out = 50), "total.burned" = exp(preds$fit), "se1" = exp(preds$fit + preds$se.fit), "se2" = exp(preds$fit - preds$se.fit))
newdata2 <- data.frame("yycumu" = seq(0, 2300, length.out = 50), "total.burned" = exp(preds2$fit), "se1" = exp(preds2$fit + preds2$se.fit), "se2" = exp(preds2$fit - preds2$se.fit))

#plot of data points, fit line and standard errors
ggplot(fire_clim) + 
  geom_point(aes(yycumu, total.burned), size = 1) +
  geom_point(data = fire_clim[11,], aes(yycumu, total.burned,  color = "Influential Point")) +
  ylab("Total Area Burned [%]") +
  xlab("Cumulative rainfall 24 months before April") +
  geom_text_repel(aes(yycumu, total.burned, label = lab)) +
  geom_line(data = newdata, aes(yycumu, total.burned, color = "excl. 1975")) +
  geom_line(data = newdata, aes(yycumu, se1, color = "excl. 1975 SE"), lty = 2) +
  geom_line(data = newdata, aes(yycumu, se2, color = "excl. 1975 SE"), lty = 2) +
  geom_line(data = newdata2, aes(yycumu, total.burned, color = "incl. 1975")) +
  geom_line(data = newdata2, aes(yycumu, se1,  color = "incl. 1975 SE"),lty = 2) +
  geom_line(data = newdata2, aes(yycumu, se2,  color = "incl. 1975 SE"), lty = 2) +
  scale_colour_manual(name = "", values = c("excl. 1975" = "grey", "incl. 1975" = "black", "excl. 1975 SE" = "grey", "incl. 1975 SE" = "black", "Influential Point" = "red")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(1,2, 1, 2,0), shape=c( NA, NA, NA, NA, 16)))) +
  theme(legend.position = c(0.2, 0.8), legend.key.size = unit(1.3, "cm"))


preds2 <- preds1$fit + rtnorm(length(preds1$fit), mean(residuals(mod)), sd(residuals(mod)), lower = min(residuals(mod)), upper = max(residuals(mod)))

#leverage vs cook's distance plot
p6<-ggplot(mod2, aes(.hat, .cooksd))+geom_point(na.rm=TRUE)+stat_smooth(method="loess", na.rm=TRUE)

p6<-p6+xlab("Leverage hii")+ylab("Cook's Distance")
p6<-p6+ggtitle("Cook's dist vs Leverage hii/(1-hii)")
p6<-p6+geom_abline(slope=seq(0,3,0.5), color="gray", linetype="dashed")
p6<-p6+theme_bw()
p6
plot(mod)

plot(fire_clim$temp, fire_clim$total.count)
mod <- glm(fire_clim$total.count ~ fire_clim$temp)
plot(mod2)


#prec vs average fire size [significant at 1%]
mod <- glm(log(average.size) ~ yycumu, data = fire_clim)
mod2 <- glm(log(average.size) ~ yycumu, data = fire_clim_subs)

preds <- predict(mod, newdata = data.frame("yycumu" = seq(0,2300, length.out = 50)), se = T)
preds2 <- predict(mod2, newdata = data.frame("yycumu" = seq(0,2300, length.out = 50)), se = T)
newdata <- data.frame("yycumu" = seq(0, 2300, length.out = 50), "average.size" = exp(preds$fit), "se1" = exp(preds$fit + preds$se.fit), "se2" = exp(preds$fit - preds$se.fit))
newdata2 <- data.frame("yycumu" = seq(0, 2300, length.out = 50), "average.size" = exp(preds2$fit), "se1" = exp(preds2$fit + preds2$se.fit), "se2" = exp(preds2$fit - preds2$se.fit))
ggplot(fire_clim) + 

  geom_point(aes(yycumu, average.size), size = 1) +
  geom_point(data = fire_clim[11,], aes(yycumu, average.size,  color = "Influential Point")) +
  ylab("Average fire size [%]") +
  xlab("Cumulative rainfall 24 months before April") +
  geom_text_repel(aes(yycumu, average.size, label = lab)) +
  geom_line(data = newdata, aes(yycumu, average.size, color = "excl. 1975")) +
  geom_line(data = newdata, aes(yycumu, se1, color = "excl. 1975 SE"), lty = 2) +
  geom_line(data = newdata, aes(yycumu, se2, color = "excl. 1975 SE"), lty = 2) +
  geom_line(data = newdata2, aes(yycumu, average.size, color = "incl. 1975")) +
  geom_line(data = newdata2, aes(yycumu, se1,  color = "incl. 1975 SE"),lty = 2) +
  geom_line(data = newdata2, aes(yycumu, se2,  color = "incl. 1975 SE"), lty = 2) +
  scale_colour_manual(name = "", values = c("excl. 1975" = "grey", "incl. 1975" = "black", "excl. 1975 SE" = "grey", "incl. 1975 SE" = "black", "Influential Point" = "red")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(1,2, 1, 2,0), shape=c( NA, NA, NA, NA, 16)))) +
  theme(legend.position = c(0.2, 0.8), legend.key.size = unit(1.3, "cm"))

#prec vs fire number [no significant relationship]
mod <- glm(log(total.count) ~ yycumu, data = fire_clim[fire_clim$year > 1975,])
summary(mod)
plot(fire_clim$yycumu[fire_clim$year > 1975], fire_clim$total.count[fire_clim$year > 1975])
lines(0:2500, exp(predict(mod, newdata = data.frame("yycumu" = 0:2500))), col = "red")

#prec vs max fire [significant at 0.1%]
mod <- glm(log(max.size) ~ yycumu, data = fire_clim)
mod2 <- glm(log(max.size) ~ yycumu, data = fire_clim_subs)
preds <- predict(mod, newdata = data.frame("yycumu" = seq(0,1600, length.out = 50)), se = T)
preds2 <- predict(mod2, newdata = data.frame("yycumu" = seq(0,1600, length.out = 50)), se = T)
newdata <- data.frame("yycumu" = seq(0, 1600, length.out = 50), "max.size" = exp(preds$fit), "se1" = exp(preds$fit + preds$se.fit), "se2" = exp(preds$fit - preds$se.fit))
newdata2 <- data.frame("yycumu" = seq(0, 1600, length.out = 50), "max.size" = exp(preds2$fit), "se1" = exp(preds2$fit + preds2$se.fit), "se2" = exp(preds2$fit - preds2$se.fit))

ggplot(fire_clim) + 
  geom_point(aes(yycumu, max.size), size = 1) +
  geom_point(data = fire_clim[11,], aes(yycumu, max.size,  color = "Influential Point")) +
  ylab("Total Area Burned [%]") +
  xlab("Cumulative rainfall 24 months before April") +
  geom_text_repel(aes(yycumu, max.size, label = lab)) +
  geom_line(data = newdata, aes(yycumu, max.size, color = "excl. 1975")) +
  geom_line(data = newdata, aes(yycumu, se1, color = "excl. 1975 SE"), lty = 2) +
  geom_line(data = newdata, aes(yycumu, se2, color = "excl. 1975 SE"), lty = 2) +
  geom_line(data = newdata2, aes(yycumu, max.size, color = "incl. 1975")) +
  geom_line(data = newdata2, aes(yycumu, se1,  color = "incl. 1975 SE"),lty = 2) +
  geom_line(data = newdata2, aes(yycumu, se2,  color = "incl. 1975 SE"), lty = 2) +
  scale_colour_manual(name = "", values = c("excl. 1975" = "grey", "incl. 1975" = "black", "excl. 1975 SE" = "grey", "incl. 1975 SE" = "black", "Influential Point" = "red")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(1,2, 1, 2,0), shape=c( NA, NA, NA, NA, 16)))) +
  theme(legend.position = c(0.2, 0.8), legend.key.size = unit(1.3, "cm"))

sum((fs.shp[[j]]$Shape_Area / sum(fs.shp[[j]]$Shape_Area)) * fs.shp[[j]]$BurnMonth)

fire.months <- vector()
for (i in 1:length(fs.shp)){
  fire.months <- c(fire.months, fs.shp[[i]]$BurnMonth)
}