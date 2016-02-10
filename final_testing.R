attach(all)
install.packages("truncdist")
library(truncdist)

minmax <- tapply(area.subs, years.subs, function(x) c(min(x), max(x)))
trunc.lnorm <- function (x) fitdist(x, "truncated_log_normal", start = c(meanlog=0, sdlog=1), method = "mle")

j <- tapply(area.subs, years.subs)
trunc.lnorm.pars <- list()
for (i in 1:length(mini)){
  dtruncated_log_normal <- function(x, meanlog, sdlog) 
    dtrunc(x, "lnorm", a=minmax[[i]][1], b=minmax[[i]][2], coef = list( meanlog=meanlog, sdlog= sdlog))
  ptruncated_log_normal <- function(x, meanlog, sdlog) 
    ptrunc(x, "lnorm", a=minmax[[i]][1], b=minmax[[i]][2], coef = list( meanlog=meanlog, sdlog=sdlog))
  trunc.lnorm.pars[[i]] <-  fitdist(area[which(j == i)] * 100, "truncated_log_normal", start = c(-1,0.1), method = "mle", control=list(trace=1, REPORT=1))
}

?ptrunc
t <- sample(1:length(lnmeans), 1) # choose distribution
ml <- lnmeans[t] #mean log
sdl <- lnscales[t] # sd log
target_distr <- rlnorm(10000,ml, sdl) # sample target distr
target_distr_sub <- target_distr[target_distr < max(area.subs[which(yr_inds == t)]) & target_distr > min(area.subs[which(yr_inds == t)])] #limit distribution to min max (will move the mean by a little bit)