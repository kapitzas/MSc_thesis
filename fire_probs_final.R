#TABLE
subs <- total
time <-as.numeric(dimnames(subs)[[1]])

subs <- c(subs, 200)
time <- c(time, 2018)
At <- data.frame(subs, time)
At$perc <- At$subs/sum(At$subs)
At$yr <- rev(2015 + 1 - time)
At$cumu <- cumsum(At$perc)


#PLOT
dev.off()

a <- seq(-10,10, length.out = 2000)
b <- seq(-10,10, length.out = 2000)
c <- seq(-10,10, length.out = 2000)
mod_probs
while(is.null(mod_probs)){
mod_probs <- tryCatch(nls(cumu ~  c * exp(-(yr/b)^a), data = At, start = list(a = sample(a, 1), b = sample(b,1), c = sample(c,1))), error = function(e) NULL)
}


summary(mod_probs)
dev.off()
plot(At$yr,At[,5], xlim = c(min(At$yr), max(At$yr)))
lines(1:120, predict(mod_probs[[5]], newdata = data.frame(yr = 1:120), type = "l"), col = "green")

p <- predict(newdata = data.frame("yr" = tsf_cur[]), mod_probs) #probabilities
test <- maskrast
values(test) <- p
plot(test)
-log(1-At$cumu)
sum(!is.na(p))
length(unique(mfi))

dt <- function(y, t){
  f <- numeric()
  for(i in 1:(length(y)-1)){
    f[i] <- (y[i + 1] - y[i]) / (t[i+1] - t[i])
  }
  return(f)
}
plot(cumu[1:28], dt(At[,4], At$yr))


plot(At$yr, At[,6])
lines(1:60, predict(mod_probs[[5]], newdata = data.frame(yr = 1:60)), type = "l")

pyric_class <- list()
for (i in 1:length(unique(mfi))){
  pyric_class[[i]] <- which(mfi == unique(mfi)[i])
}
values(test) <- mfi
plot(test)
