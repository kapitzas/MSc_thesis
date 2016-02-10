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
