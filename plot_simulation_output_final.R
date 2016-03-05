require(ggplot2)
sf <- read.csv("/Users/Simon/Studium/MSC/Masterarbeit/data/new_model_out2/fires_high.csv")
stb <- read.csv("/Users/Simon/Studium/MSC/Masterarbeit/data/new_model_out2/tab_high.csv")

#tab simulations
ci_stb <- aggregate(stb$tab_target, by = list(c(stb$yr)), FUN =  function(x) ci(x))
quants <- aggregate(stb$tab_target, by = list(c(stb$yr)), FUN =  quantile, c(0.025, 0.5,0.975))

#Plot tab
df = data.frame(year = quants$Group.1, median = quants$x[,2], lower = quants$x[,1], upper = quants$x[,3], mean = ci_stb$x[,2])

ggplot(data = df,aes(year, median)) +
  geom_polygon(data = data.frame(y = c(0,0, 2000 ,2000), x = c(2090,2016, 2016, 2090)), aes(x = x, y = y), fill = "grey", alpha = .5) +
  geom_ribbon(aes(x=year, ymax=upper, ymin=lower), fill="blue", alpha=.3) +
  geom_line(aes(y = mean), col = "black", alpha = .3) +
  geom_line(aes(y = median), col = "black", alpha = .3)

#fire size simulations
sf_agr <- aggregate(sf$fire_target ~ sf$sim_no + sf$yr, FUN = "mean")
ci_sf <- aggregate(sf_agr[,3], by = list(c(sf_agr[,2])), FUN =  function(x) ci(x))
quants <- aggregate(sf_agr[,3], by = list(c(sf_agr[,2])), FUN =  quantile, c(0.1, 0.5, 0.9))

#Plot fire sizes
df = data.frame(year = quants$Group.1, median = quants$x[,2], lower = quants$x[,1], upper = quants$x[,3], mean = ci_sf$x[,2])
ggplot(data = df,aes(year, median)) +
  geom_polygon(data = data.frame(y = c(0,0, 120 ,120), x = c(2090,2016, 2016, 2090)), aes(x = x, y = y), fill = "grey", alpha = .5) +
  geom_ribbon(aes(x=year, ymax=upper, ymin=lower), fill="red", alpha=.2) +
  geom_line(aes(y = median), col = "black", alpha = .4) +
  geom_line(aes(y = mean), col = "red", alpha = .4) +
  stat_smooth(span = 20, se = F, alpha = .4, col = "black", size = .4)

#fire occurance simulations
fire_occ <- tapply(stb$sim_no, stb$yr, length)/100

#plot fire occuance
df <- data.frame(fire_occ = as.numeric(fire_occ), year = as.numeric(dimnames(fire_occ)[[1]]))
ggplot(df, aes(x = year, y = fire_occ)) +
  geom_point() +
  stat_smooth(span = 4)

#prec
a <- mean_prec + rnorm(1000, 0, sd(reference_prec))
b <- data
c <- reference_prec
lines(density(a))
plot(density(b))
lines(density(c), col = "red")

ci_prec <- aggregate(stb$precip, by = list(c(stb$yr)), FUN =  function(x) ci(x))
quants <- aggregate(stb$precip, by = list(c(stb$yr)), FUN =  quantile, c(0.025, 0.5,0.975))

#Plot prec
df = data.frame(year = quants$Group.1, median = quants$x[,2], lower = quants$x[,1], upper = quants$x[,3], mean = ci_prec$x[,2])
ggplot(data = df,aes(year, median)) +
  geom_polygon(data = data.frame(y = c(0,0, 2000 ,2000), x = c(2090,2016, 2016, 2090)), aes(x = x, y = y), fill = "grey", alpha = .5) +
  geom_ribbon(aes(x=year, ymax=upper, ymin=lower), fill="blue", alpha=.3) +
  geom_line(aes(y = mean), col = "black", alpha = .3) +
  geom_line(aes(y = median), col = "blue", alpha = .3) +
  geom_smooth(span = 1)


# density plots of real and predicted fire sizes (logs)
# s <- sample(1:length(sf$precip[which(sf$yr < 2015)]), 10000)

ggplot() +
  geom_density(aes(log(area[which(year < 2015)])), alpha = .2, fill = "red") +
  geom_density(aes(log(sf$fire_target[which(sf$yr < 2015)])), fill = "blue", alpha = 0.2)

#density plots of real and predicted tab values
s <- 1:length(stb$precip[which(stb$yr < 2015)])
ggplot() +
  geom_density(aes(log(total)), alpha = .2, fill = "red") +
  geom_density(aes(log(stb$tab_target[which(stb$yr < 2016)]/sum(na.omit(cell_size[]))*100)), fill = "blue", alpha = 0.2)

#density plots of real and predicted prec values
s <- 1:length(stb$precip[which(stb$yr < 2015)])
ggplot() +
  geom_density(aes(log(prec.stats$yycumu[which(prec.stats$year < 2006 & prec.stats$year > 1985)])) , alpha = .2, fill = "red") +
  geom_density(aes(log(test$precip[which(test$yr < 2006 & test$yr > 1985)])), fill = "blue", alpha = 0.2)

install.packages("mixtools")


#Mixed model comparison
test  <- read.csv("/Users/Simon/Studium/MSC/Masterarbeit/data/non_spatial_new_model_out/ns_total_area_low_test.csv")

#Plot cumu prec
ggplot() +
  geom_density(aes(stb$precip[which(stb$yr < 2005 & stb$yr >= 1986)]), fill = "blue", alpha = 0.2) +
  geom_density(aes(prec.stats$yycumu[which(prec.stats$year < 2005 & prec.stats$year >= 1986)]) , alpha = .2, fill = "red") +
  geom_density(aes(test$precip[which(test$yr < 2005 & test$yr >= 1986)]), fill = "green", alpha = 0.2)

#Plot tab
attach(all)
ggplot() +
  geom_density(aes(log(total[which(unique(year) <= 2005 & unique(year) >= 1986)]/100 * sum(cell_size, na.rm = T))) , alpha = .2, fill = "red") +
  geom_density(aes(log(stb$tab_target[which(stb$yr < 2016)])), fill = "blue", alpha = 0.2)

# rm(list = ls(all = T))

df <- data.frame(total, cumu3)
preds1 <- predict(mod1, newdata = data.frame("cumu3" = 600:2500), se.fit = T)
preds3 <- predict(mod3, newdata = data.frame("cumu3" = 600:2500), se.fit = T)
cumu3 <- 600:2500

dfmod <- data.frame(preds_mod1 =exp(preds1$fit), preds_mod3 = (preds3$fit)^2, cumu3 = cumu3)
dfmod$upper_mod1 <- exp(preds1$fit + preds1$se.fit)
dfmod$lower_mod1 <- exp(preds1$fit - preds1$se.fit)
dfmod$upper_mod3 <- (preds3$fit + preds3$se.fit)^2
dfmod$lower_mod3 <- (preds3$fit - preds3$se.fit)^2
ggplot() +
  #mod1
  geom_ribbon(data = dfmod, aes(x=cumu3, ymax=upper_mod1, ymin=lower_mod1), fill="blue", alpha=.3) +
  geom_ribbon(data = dfmod, aes(x=cumu3, ymax=upper_mod3, ymin=lower_mod3), fill="red", alpha=.3) +
  geom_line(data = dfmod, aes(x = cumu3, y = preds_mod1), col = "black", alpha = .4) + 
  geom_line(data = dfmod, aes(x = cumu3, y = preds_mod3), col = "black", alpha = .4) +
  geom_point(data = df, aes(cumu3, total))

