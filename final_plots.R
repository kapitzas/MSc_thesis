#PLOT TOTAL AREA BURNED VS PREC
df <- data.frame(total, cumu)
plot(cumu, total)
n.cumu <- seq(0,2300, length.out = 50)
preds <- predict(mod_total, newdata = data.frame("cumu" = n.cumu), se = T)

n.df <- data.frame("cumu" = n.cumu, "fit" = exp(preds$fit), "se1" = exp(preds$fit + preds$se.fit), "se2" = exp(preds$fit - preds$se.fit))

ggplot(df, aes(cumu, total)) +
  geom_point() +
  labs(y = "Total Area Burned [km^2]", x = "Cumulative Rain") +
  geom_line(aes(cumu, fit), n.df, linetype = 1) +
  geom_line(aes(cumu, se1), n.df, linetype = 2) +
  geom_line(aes(cumu, se2), n.df, linetype = 2)


#PLOT PEARSON MOMENTS VS PREC
df <- data.frame(cumu.subs, mean = MMs[,1], variance = MMs[,2], skewness = MMs[,3], kurtosis = MMs[,3])

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

#PLOT UNFITTED LN MOMENTS VS PREC
plot(lnmeans ~ cumu.subs)
lnmeans_mod <- lm(lnmeans ~ cumu.subs)
lines(cumu.subs, predict(lnmeans_mod))


plot(lnscales ~ cumu.subs)
lnscales_mod <- lm(lnscales ~ cumu.subs)
lines(cumu.subs, predict(lnscales_mod))

plot(lnscales ~lnmeans)
plot(residuals(lnscales_mod) ~ residuals(lnmeans_mod))

###QQPLOTS

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


#plot qq plots of unfitted lnorm functions
require(grid)
dat <- tapply(area.subs, years.subs, function(x) list(x))
grid.newpage()
pushViewport(viewport(layout = grid.layout(6,3)))
for (i in 1:length(dat)){
  df <- data.frame("y" = dat[[i]])
  x <- qlnorm(c(0.25, 0.75), sdlog = lnscales[i], meanlog = lnmeans[[i]])
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
