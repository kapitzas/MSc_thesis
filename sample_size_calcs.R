#create 20 m res grid
#determine average shape size
# grow shapes into grid
# - for each shape, until target shape size is reached
# - until total area occupied by shapes is reached
nls <- list()
ta.burned.real <- vector()
deviation <- vector()
ta.burned <- seq(1,50,by = 0.5)/100
for (k in 1:length(ta.burned)){
  r <- raster()
  extent(r) <- c(0, 250, 0, 250)
  res(r) <- 1
  crs(r) <- NULL
  r[] <- 1
  #input params
  
  ta <- ta.burned[k]
  unburned <- 0
  unburned.clus <- 0
  
  #target area occupied by polygons
  target.area <- ncell(r) * ta
  
  b <- numeric()
  while(length(b) < target.area){
    seed <- sampleRandom(r, size = 1, cells=T)[,1]
    cells <- seed
    cluster.size <- round(runif(1,1,1000))
    bc <- numeric()
    
    while(length(bc) < cluster.size){
      cells <- adjacent(r, cells, directions = 16, pairs = F)
      bc <- unique(c(cells, bc))
    }
    b <- unique(c(bc, b))
  }
  
  r[b] <- 0
  
  real <- length(r[which(r[] == 0)])/length(r[])
  
  
  #create random transects
  samples.2.5.v <- vector()
  #samples.2.5.h <- vector()
  
  transects <- round(runif(100, 1, length(r[,1])))
  
for (i in 1:length(transects)){
    trans <- transects[i]
    
    dist <- 20
    #sample observation points on transect
    points <- round(seq(dist,250, by = dist) + runif(250/dist,1,dist-1))
    points[which(points > 250)] <- 250
    
    #simulate sampling every 20 m (sampling points are stored in points)
    samples.2.5.v[i] <- length(r[points,trans][which(r[points,trans] == 0)]) / length(points)
    #samples.2.5.h[i] <- length(r[trans,points][which(r[trans,points] == 0)]) / length(points)
  }
  
  
  #calculate RMSE of sampled data vs the real value
  #different transect lengths to predict to
  trans.length <- c(seq(0.5,3.5, by = 0.25), 4.0, 4.5, 5, 7.5, 10)
  obs <- vector()
  obs.list <- list()
  residuals <- list()
  rmse <- vector()
  
  #number of transects
  s <- c(2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,30,40)
  for (j in 1:18) {
    for (i in 1:100) {
      obs[i] <- mean(sample(samples.2.5.v, size = s[j]))
    }
    obs.list[[j]] <- obs
    residuals[[j]] <- obs.list[[j]] - real
    rmse[j] <- sqrt(mean((obs.list[[j]]-real)^2))
  }
  
  #make nls from rmse data
  nls[[k]] <- nls(rmse ~ i*trans.length^-z, start= list(i=-3,z=-2))
  print(paste(k, "%"))
  ta.burned.real[k] <- real
  deviation[k] <- 10/100 * real
}


#inverse function
predict.length <- function(fun, value){
  i <- coefficients(fun)[[1]]
  z <- coefficients(fun)[[2]]
  return((value/i)^(-1/z))
}

#Plot
dev.off()
layout(t(1:2), widths=c(5,0.2))

par(mar=rep(0.5, 4), oma=rep(3, 4), las=1)

plot(0:100, seq(0, 0.2, by = 0.01/5), type = "n", ylab = "RMSE")

colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))
cl<- colfunc(99)

trans.length.new <- c(seq(0,100, by = .01))
for (i in 1:length(nls)){
lines(trans.length.new, predict(nls[[i]], newdata = data.frame("trans.length" = trans.length.new)), col = cl[i])
  }

length <- vector()
for (i in 1:length(nls)){
length[i] <- predict.length(nls[[i]], deviation[i])
}


image(1, ta.burned, t(seq_along(ta.burned)), col=cl, axes=FALSE)
axis(4)

plot(ta.burned.real, length, ylim = c(0, max(length) + 50))
warnings()
quantile(length, c(.1, .5, .90))

plot(r)
