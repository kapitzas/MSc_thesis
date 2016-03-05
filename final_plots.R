require(fitdistrplus)
require(MASS)
image_path <- "/Users/Simon/Studium/MSC/Masterarbeit/write-up/Final/images/"
# Plot loglikelihoods of fire distribution fits
df <- data.frame(unique(years.subs), "lognorm fitted" = ll_lnorm_fit, "lognorm unfitted" = ll_lnorm_unfit, "fitted pearson" = ll_pears, "fitted gamma" = ll_gamma)
df <- melt(df, id.vars = "unique.years.subs.")
pdf(paste0(image_path,"loglik_fd.pdf"), 7,7)
ggplot(data = df, aes(x = unique.years.subs., y = value, group = variable)) +
  geom_point(aes(shape=variable, col = variable)) +
  theme(legend.position = c(0.1, 0.9), legend.title=element_blank()) +
  labs(x = "year", y = "LogLik")
#PLOT TOTAL AREA BURNED VS PREC
dev.off()
df <- data.frame(total, cumu3)
plot(cumu3, total)
n.cumu <- seq(0,2300, length.out = 50)
preds <- predict(tab_mod2, newdata = data.frame("cumu3" = n.cumu), se = T)

n.df <- data.frame("cumu3" = n.cumu, "fit" = exp(preds$fit), "se1" = exp(preds$fit + preds$se.fit), "se2" = exp(preds$fit - preds$se.fit))

ggplot(df, aes(cumu3, total)) +
  geom_point() +
  labs(y = "Total Area Burned [km^2]", x = "Cumulative Rain") +
  geom_line(aes(cumu3, fit), n.df, linetype = 1) +
  geom_line(aes(cumu3, se1), n.df, linetype = 2) +
  geom_line(aes(cumu3, se2), n.df, linetype = 2)


#PLOT PEARSON MOMENTS VS PREC
df <- data.frame(cumu3.subs, mean = MMs[,1], variance = MMs[,2], skewness = MMs[,3], kurtosis = MMs[,3])

a <- ggplot(df, aes(x = cumu3.subs, y = mean)) + 
  geom_point(size = 1) +
  labs(x="Cumulative precipitation", y="Mean")

b <- ggplot(df, aes(x = cumu3.subs, y =variance)) + 
  geom_point(size = 1) +
  labs(x="Cumulative precipitation", y="Variance")

c <- ggplot(df, aes(x = cumu3.subs, y =skewness)) + 
  geom_point(size = 1) +
  labs(x="Cumulative precipitation", y="Skewness")

d <- ggplot(df, aes(x = cumu3.subs, y =kurtosis)) + 
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
plot(lnmeans ~ cumu3.subs)
lnmeans_mod <- lm(lnmeans ~ cumu3.subs)
lines(cumu3.subs, predict(lnmeans_mod))


plot(lnscales ~ cumu3.subs)
lnscales_mod <- lm(lnscales ~ cumu3.subs)
lines(cumu3.subs, predict(lnscales_mod))

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
pdf(paste0(image_path,"qq_lnorm_fit.pdf"), width =7, height = 12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(5,3)))
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
#plot qq plots of unfitted lnorm functions
require(grid)
dat <- tapply(area.subs, years.subs, function(x) list(x))
pdf(paste0(image_path,"qq_lnorm_unfit.pdf"), width =7, height = 12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(5,3)))
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
###Plot fire probability curves
pdf(paste0(image_path,"prob_weib.pdf"), width =7, height = 12)
grid.newpage()
yloc <- c(1,1,2,2,3)
xloc <- c(1,2,1,2,1)
df <- data.frame()
pushViewport(viewport(layout = grid.layout(3,2)))
for (i in 1:length(mod_probs)){
  df <- data.frame(cumu = Ftl[,i], yr = Ftl[,6])
  nd <- seq(0,35, by = 0.1)
  df2 <- data.frame(nd = nd, preds = predict(mod_probs[[i]], newdata = data.frame(yr = nd)))
  a <- ggplot()+
    geom_line(data = df2, aes(nd, preds)) +
    geom_point(data = df, aes(yr, cumu)) + 
    labs(x = "t", y = "p") +
    ggtitle(paste("class", i))
  print(a, vp = viewport(layout.pos.row = yloc[i], layout.pos.col = xloc[i]))
}
dev.off()
mod_probs[[i]]
