#FUNCTION TO CALCULATE PLANAR, UNPROJECTED DISTANCE (PLANAR UNITS) BETWEEN CELLS
dist <- function(from, to, ncol, nrow){
  row_from <- ceiling(from/ncol)
  col_from <- from - ((row_from - 1) * ncol)
  row_to <- ceiling(to/ncol)
  col_to <- to - ((row_to - 1) * ncol)
  dist_y <- row_from - row_to
  dist_x <- col_from - col_to
  distance <- sqrt(dist_x^2 + dist_y^2)
  return(distance)
}

#FUNCTION TO CLAMP VALUES
clamp <- function(x, min, max) {
  x[which(x<min)] <- min
  x[which(x>max)] <- max
  return(x)
}

#LOGLIK OF PEARSON FITS
loglikpears <- function(x){
  temp <- pearsonMSC(x)
  max(temp$logLik)
}

#PEARSON FUNCTION WITH HIGHEST AIC
pearsFUN <- function(x){
  temp <- pearsonMSC(x)
  temp$Best$AIC
}

#FUNCTION TO ESTIMATE EMPIRIC LOGSD
lnscale <- function(x){
  variance <- (sd(x))^2
  mean <- mean(x)
  scale <- sqrt(log(1+ variance/mean^2))
  return(scale)
}

#FUNCTION TO ESTIMATE EMPIRIC LOGMEAN
lnmean <- function(x){
  variance <- (sd(x))^2
  mean <- mean(x)
  scale <- sqrt(log(1+ variance/mean^2))
  loc <- log(mean) - (0.5 * scale^2)
  return(loc)
}
