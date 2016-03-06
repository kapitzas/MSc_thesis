attach(all)
head(all)
no_of_fires <- tapply(area, year, length)
yrs <- 1970:2015
oc <- 1:length(yrs)
oc[match(unique(year), yrs)] <- 1
oc[-match(unique(year), yrs)] <- 0


pr <- prec.stats$yycumu[which(prec.stats$year%in%yrs)]

mod_occ <- glm(oc ~ pr, family = "binomial")
lines(0:2000, predict(mod, newdata = data.frame("pr" = 0:2000), type = "response"))
plot(oc ~ pr, xlim = c(0,2000))
summary(mod)
predict(mod, newdata = data.frame("pr" = 700), type = "response")

-log(1- At[,2])
At[,2]

plot(Ft[,1] ~ unique(all$cumu_prec)[5:28])
