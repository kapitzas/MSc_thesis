#create grid
#determine average shape size
# grow shapes into grid
# - for each shape, until target shape size is reached
# - until total area occupied by shapes is reached

###################
#1. Definitions####
###################

#1.a) Storage variables
nls_a <- list() #For non linear models
nls_b <- list()
ta_unburned_real <- vector() #realised total area burned
deviation <- vector() #deviation 

#1.b)
total_clusters <- seq(1,50,by = 0.5)/100

##########################################
#2. Grow unburned patches into raster#####
##########################################
for (k in 1:length(total_clusters)){
  #2.a) Create empty raster at 1m resolution
  r <- raster()
  extent(r) <- c(0, 250, 0, 250)
  res(r) <- 1
  crs(r) <- NULL
  r[] <- 0 #0 means "burned", since we're looking at a presumably burned cell
  
  #2.b) Calculate number of unburned cells
  target_area <- ncell(r) * total_clusters[k]
  
  #Storage of unburned cells
  b <- numeric()
  
  #2.c) Initialise unburned clusters on raster
  while(length(b) < target_area){
    
    #Random seed from which to grow unburned cluster
    seed <- sampleRandom(r, size = 1, cells=T)[,1]
    cells <- seed
    
    #Sample the size of  unburned cluster
    cluster_size <- round(runif(1,1,200))
    
    #Storage of unburned cluster cells
    bc <- numeric()
    
    #2.d) Grow clusters on raster
    while(length(bc) < cluster_size && length(b) < target_area){
      
      #Sample adjacent cells
      cells <- adjacent(r, cells, directions = 16, pairs = F)
      
      #Store unburned cluster cells
      bc <- unique(c(cells, bc))
      b <- unique(c(bc, b))
    }
  }
  
  #2.d) Identify realised ratio of landscape occupied by unburned patches
  r[b] <- 1
  real <- length(r[which(r[] == 1)])/length(r[])
  
  ##################################################
  #3.)Sampling of landcover on created landscape####
  ##################################################
  
  #3.a) Create 100 random transects
  transects <- round(runif(100, 3, length(r[1,])-2))
  
  #3.b) Sampling of landcover on random points on transect
  
  samples_a <- vector()
  samples_b <- vector()
  for (i in 1:length(transects)){

    #Determine observation points on transect (random point per 25m section)
    points <- round(seq(0, 225, by = 25) + runif(250/25,1,25-1))
    
    #Sample land cover on each transect point
    #1m transect
    samples_a[i] <- sum(r[points,transects[i]][which(r[points,transects[i]] == 1)]) / length(points)
    #5m transect
    samples_b[i] <- sum(r[points, transects[i]] +
    r[points, transects[i]-1] +
    r[points, transects[i]-2] +
    r[points, transects[i]+1] +
    r[points, transects[i]+2])/(5*length(points))
  }
  
  ##############################################
  #4.) Calculation of RMSE of sampled values####
  ##############################################
  
  #4.a) Determine different transect lengths
  #Transect lengths and the number of transects necessary to achieve the length
  trans_length <- c(seq(0.5,3.5, by = 0.25), 4.0, 4.5, 5, 7.5, 10)
  s <- c(2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,30,40)
  
  #4.b) Calculate simulated observed ratio of unburned patches on different transect lengths
  obs_a <- vector()
  obs_b <- vector()
  obs_list_a <- list()
  obs_list_b <- list()
  residuals_a <- list()
  residuals_b <- list()
  rmse_a <- vector()
  rmse_b <- vector()
  
  for (j in 1:18) {
    for (i in 1:100) {
      #Drawing random sampled transects and building the mean of observations
      obs_a[i] <- mean(sample(samples_a, size = s[j]))
      obs_b[i] <- mean(sample(samples_b, size = s[j]))
    }
    #Calculating the residuals and rmse for each transect length
    obs_list_a[[j]] <- obs_a
    obs_list_b[[j]] <- obs_b
    residuals_a[[j]] <- obs_list_a[[j]] - real
    residuals_b[[j]] <- obs_list_b[[j]] - real
    rmse_a[j] <- sqrt(mean((obs_list_a[[j]]-real)^2))
    rmse_b[j] <- sqrt(mean((obs_list_b[[j]]-real)^2))
  }
  
  #4.c) Make nls of RMSE vs trans_length for plot and estimate 10% deviation interval from the realised value
  nls_a[[k]] <- nls(rmse_a ~ i*trans_length^-z, start= list(i=-3,z=-2))
  nls_b[[k]] <- nls(rmse_b ~ i*trans_length^-z, start= list(i=-3,z=-2))
  print(paste(k, "%"))
  ta_unburned_real[k] <- real
  deviation[k] <- 10/100 * real
}


#inverse function
predictLength <- function(fun, value){
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

trans_length_new <- c(seq(0,100, by = .01))

#1m transect
for (i in 1:length(nls_a)){
  lines(trans_length_new, predict(nls_a[[i]], newdata = data.frame("trans_length" = trans_length_new)), col = cl[i])
}

#5m transect
for (i in 1:length(nls_b)){
  lines(trans_length_new, predict(nls_b[[i]], newdata = data.frame("trans_length" = trans_length_new)), col = cl[i])
}

length_a <- vector()
length_b <- vector()

for (i in 1:length(nls_a)){
  length_a[i] <- predictLength(nls_a[[i]], deviation[i])
  length_b[i] <- predictLength(nls_b[[i]], deviation[i])
}


image(1, ta_unburned, t(seq_along(ta_unburned)), col=cl, axes=FALSE)
axis(4)
dev.off()
plot(ta_unburned_real, length_a, ylim = c(0, max(length_a) + 50), type = "l")
lines(ta_unburned_real, length_b, ylim = c(0, max(length_b) + 50), col = "red")

quantile(length_a, c(.075, .5, .95))
quantile(length_b, c(.075, .5, .95))
hist(sqrt(length_a))
plot(r)
