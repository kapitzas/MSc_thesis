#Simulation of precepitation patterns

#87 - 2005 reference
reference.prec <- fire_clim$prec[fire_clim$year > 1987 & fire_clim$year < 2005]
reference.average <- mean(reference.prec)
high_prec_30 <- reference.average *108/100
low_prec_30 <- reference.average * 89/100
high_prec_90 <- reference.average *119/100
low_prec_90 <- reference.average * 69/100

#mean of prec of reference time span
mean_prec <- mean(reference.average)

#calc of increments
years_30 <- 2030 - 2015
years_90 <- 2090 -2030
incr_15_30 <- (high_prec_30 - mean_prec) / years_30
decr_15_30 <- (low_prec_30 - mean_prec) / years_30
incr_30_90 <- (high_prec_90 - high_prec_30) / years_90
decr_30_90 <- (low_prec_90 - low_prec_30) / years_90

incr <- c(rep(incr_15_30, times = years_30), rep(incr_30_90, times = years_90))
decr <- c(rep(decr_15_30, times = years_30), rep(decr_30_90, times = years_90))

#means 2015 - 2090 (first value for 2016)
sim_prec_high <- numeric()
sim_prec_low <- numeric()
for (i in 1:length(incr)){
sim_prec_high[i] <- exp(log(reference.average + sum(incr[1:i])) + rnorm(1, 0, sd(log(reference.prec))))
sim_prec_low[i] <- exp(log(reference.average + sum(decr[1:i])) + rnorm(1, 0, sd(log(reference.prec))))
}

#simulated cumulative prec
for (i in 3:length(sim_prec_high)){
  sim_cumu_high[i] <-  sim_prec_high[i] + sim_prec_high[i-1] + sim_prec_high[i-2]
  sim_cumu_low[i] <-  sim_prec_low[i] + sim_prec_low[i-1] + sim_prec_low[i-2]
}

sim_high_total_burned <- exp(predict(mod, newdata = (data.frame("yycumu" = sim_cumu_high))))
sim_low_total_burned <- exp(predict(mod, newdata = (data.frame("yycumu" = sim_cumu_low))))
plot(exp(predict(mod, newdata = (data.frame("yycumu" = sim_cumu_high)))))

all$prec <- precs
names(all) <- c("area", "year", "fire_number", "yy_cumu_prec")
fires <- all[,-3]

save(fires, file = "fires_newhaven.RData")

?save
sim.length <- 75
#sim.reps <- 20

#specify output path
output <- "/Users/Simon/Studium/MSC/Masterarbeit/data/new_model_out"

elev.df <- as.data.frame(elev)
tsf.df <- data.frame("tsf.ini" = tsf.ts[[28]][])
tsf <- tsf.df[,1]
fs.df <- data.frame(fire.scar.ts[])
year <- 2015

for (i in 1:sim.length){

  #update fire probability map
  #add new fire scar to fire scar df (only chagnes something for i > 1)
  names(fs.df)[which(names(fs.df) == "mask")] <- paste0("fs_s", year)
  #determine number of fires on each cell
  
  firenumber <- rowSums(fs.df[c(3:length(fs.df))])
  
  #find mean fire return interval per land type
  #number of years between first record and current time step
  no_years <- as.numeric(tail(substr(names(fs.df), 5, 8),1)) - as.numeric(head(substr(names(fs.df), 5, 8),1))
  
  # average fire number per landtype
  mfn.per.lt <- aggregate(firenumber, list(lt), mean)
  
  #average fire interval per landtype
  mfi.per.lt <- cbind(mfn.per.lt[,1], no_years/mfn.per.lt[,2])
  
  #write mfi per land type into df
  mfi.lt <- data.frame(lt = lt, mfi = lt)
  for (i in 1:length(mfi.per.lt[,1])){
    mfi.lt[which(mfi.lt[,2] == mfi.per.lt[i,1]),2] <- mfi.per.lt[i,2]
  }
  
  #probability map
  #P * B * tsf * MI^(e+2)
  B <- mfi.lt[,2]
  MI <- mfi.lt[,2]
  fire.probs <- B * tsf *  MI^-(exp(1)+2)
  
  target.year <- sample(ty, 1)
  if (target.year > size){
    while(target.year > size){
      target.year <- sample(ty, 1)
    }
  }
  
  wind.dir <- 190
  windiness <- 5
  total.fire <- 1
  
  #normalising fire.probs
  fp.df <- data.frame((fire.probs - min(na.omit(fire.probs)))/(max(na.omit(fire.probs)) - min(na.omit(fire.probs))))
  
  fire.scar <- as.data.frame(mask-1) *-1
  
  while(total.fire < target.year){
    cells.burned <- 1
    #sample fire size
    target.fire <- sample(tf, 1)
    if (target.fire > target.year){
      while(target.fire > target.year){
        target.fire <- sample(tf, 1)
      }
    }
    
    #sampling ignition cell
    ignition.cell <- sample(which(!is.na(fp.df)), size = 1, prob = fp.df[which(!is.na(fp.df)),])
    fire.scar[ignition.cell,] <- 2
    ac <- ignition.cell
    fp.df[ac,] <- 0
    fp <- data.frame(ac = numeric(), sc = numeric(), iP = numeric(), tsf = numeric(), dist = numeric(), xdist = numeric(), disty = numeric(), dir = numeric(), slope = numeric(), pburn = numeric(), burn = numeric())
    
    while(cells.burned < target.fire & length(ac) > 0 & total.fire < target.year & length(ac) > 0){
      sc <- as.data.frame(adjacent(mask, cells = ac[which(!is.na(elev.df[ac,]))], pairs = T, directions = 8, include = F, sorted =T))
      
      #filter "only burnable" cells as spreading cells
      sc <- sc[which(fp.df[sc[,2],] > 0 & !is.na(elev.df[sc[,2],])),]
      
#       #distance
#       dist <- distance(sc[,1], sc[,2], mask)
#       
#       #wind
#       spr.dir <- direction(dist[[1]], dist[[2]])
#       dow <- ifelse(wind.dir < spr.dir, spr.dir - wind.dir, wind.dir - spr.dir)
#       dow[which(dow > 180)] <- 360 - dow[which(dow > 180)]
#       wind <- (180-dow)*(pi/180)
#       
#       #rescale wind
#       NewMax <- 1 + 0.2 * windiness
#       NewMin <- 1 - 0.2 * windiness
#       wind <- ((wind * (NewMax - NewMin)) / 3.141593) + NewMin
#       
#       #slope
#       diff.elev <- elev.df[sc[,2],] - elev.df[sc[,1],]
#       slope <- diff.elev/(dist[[1]] *100)
#       slope[which(slope > 0.5)] <- 0.5
#       slope[which(slope < -0.5)] <- - 0.5
#       
#       #calculate modulation variable
#       m <- ((tsf[sc[,1]]/15 * wind)  / (dist[[1]]/0.0025)) * 1.5^(slope + 0.5)
#       
      #create iteration dataframe
      fp.new <-data.frame(ac = sc[,1], sc = sc[,2], iP = fp.df[sc[,2],], tsf = tsf[sc[,2]], dist = dist[[1]], distx = dist[[2]][,1], disty = dist[[2]][,2], dir = spr.dir, m = m)
      fp.new$pburn <- (fp.new$iP * fp.new$m)
      fp.new$burn <- ifelse(runif(length(fp.new$pburn),0,1) < fp.new$pburn, 1, 0)
      
      ac <- fp.new$sc[which(fp.new$burn == 1)]
      fp <- rbind(fp, fp.new)
      cells.burned <- cells.burned + length(ac)
      total.fire <- total.fire + length(ac)
      #restart fire on fire front, if size not reached
      if (length(fp[,1]) > 0 & cells.burned < fire.size & length(ac) == 0){
        ff <- which(fp$burn == 0)
        ac <- fp$sc[sample(ff, size = round(length(ff)/20), prob = fp$iP[ff])]
      }
      #update fire prob and fire scar maps
      fp.df[ac,] <- 0
      fire.scar[ac,] <- 1
      print(cells.burned)
    }
    total.fire <- total.fire
    print(paste("total:", total.fire))
  }



#increment year
year <- year + 1

#update fs
fs.df <- cbind(fs.df, fire.scar)

#update tsf
tsf <- tsf +1
tsf[which(fire.scar == 1)] <- 1
tsf.df <- cbind(tsf.df, tsf)
}

fire.scar.map <- raster(mask)
fire.scar.map[] <- fire.scar[,1]
#writeRaster(paste0(output,"/", tsf, sim, year), format = "ascii")
fire.scar.map[which(fire.scar.map[] == 0)] <- NA
plot(fire.scar.ts[[25]])
plot(fire.scar.map, col = c("red", "black"), add = T, legend = F)

