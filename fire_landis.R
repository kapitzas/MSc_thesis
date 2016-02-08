#LANDIS Fire probability using lf (time since fire), MI (mean fire return interval of a given cell) and B (fire prbability coefficient, B = MI by default)

#function for fire probability equation
P <- function(B, lf, MI){
  e <- 2.71828182
  P <- B * lf * MI^-(e+2)
 return(P)
}

lf <- seq(0.1, 1, 0.05)
mean(lf)
MI = 150
B <- MI

fireP <- P(B, lf, MI)

plot(fireP ~ lf, pch = 0)
lines(fireP ~ lf)
summary(fireP)

#function for fire size equation



S <- function(A, MS, r){
  S <- A * 10^r * MS
  return(S)
}

A <- 0.34
MS = 1000
a1 <- sample(runif(100, 0,1),1)
a2 <- sample(runif(100, 0,1),1)
r <- sqrt(-0.75 * log(a1) * sin(pi^2 * a2))
fireS <- S(A,MS,r)
fireS

