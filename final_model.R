require(msm)
rm(prec,year)
detach(clim.stats)
detach(all)
attach(clim.stats)

####################################################
#1) Preprocessing of data and simulation definitions
####################################################

#1.a) Preprocessing precipitation reference data

#87 - 2005 reference
reference_prec <- prec[year > 1987 & year < 2005]
mean_prec <- mean(reference_prec)
high_prec_30 <- mean_prec *108/100
low_prec_30 <- mean_prec * 89/100
high_prec_90 <- mean_prec *119/100
low_prec_90 <- mean_prec * 69/100

#calculation of annual increments/decrements
years_30 <- 2030 - 2015
years_90 <- 2090 -2030
incr_15_30 <- (high_prec_30 - mean_prec) / years_30
decr_15_30 <- (low_prec_30 - mean_prec) / years_30
incr_30_90 <- (high_prec_90 - high_prec_30) / years_90
decr_30_90 <- (low_prec_90 - low_prec_30) / years_90

incr <- c(rep(incr_15_30, times = years_30), rep(incr_30_90, times = years_90))
decr <- c(rep(decr_15_30, times = years_30), rep(decr_30_90, times = years_90))



#1.b) Initial parameter values and definitions for entire simulation

simle <- 75 #simulation length
simno <- 10 #number of simulations
sim_data <- list() #data frame for all sim data
output <- "/Users/Simon/Studium/MSC/Masterarbeit/data/new_model_out"
cell_size <- raster(paste0(path.fire, "cell_size.asc"))
cell_size <- cell_size[] #cell sizes in km^2
maskrast <- raster(paste0(path.fire, "maskraster.asc"))
study_a <- maskrast[] #study area
inds <- 1:length(maskrast[]) #cell indices

#########################
#2) SIMULATION LEVEL#####
#########################
for (j in 1:simno){
  
  #2.a) Simulation of precipitation regime for time span
  
  #Estimation of simulated precipitation means 2015 - 2090 (first value for 2016), with random terms
  sim_prec_high <- numeric(); sim_prec_low <- numeric()
  for (i in 1:length(incr)){
    sim_prec_low[i] <- exp(log(mean_prec + sum(decr[1:i])) + rnorm(1, 0, sd(log(reference_prec))))
    sim_prec_high[i] <- exp(log(mean_prec + sum(incr[1:i])) + rnorm(1, 0, sd(log(reference_prec))))
  }
  
  #Attach last two measured precipitation values (2014 and 2015) for calculation of first time step
  sim_prec_high <- c(tail(prec, 2), sim_prec_high)
  sim_prec_low <- c(tail(prec, 2), sim_prec_low)
  
  sim_cumu_high <- numeric(); sim_cumu_low <- numeric()
  
  #Calculate simulated cumulative precipitation
  for (i in 1:(length(sim_prec_high)-2)){
    sim_cumu_high[i] <-  sim_prec_high[i+2] + sim_prec_high[i+1] + sim_prec_high[i]
    sim_cumu_low[i] <-  sim_prec_low[i+2] + sim_prec_low[i+1] + sim_prec_low[i]
  }
  
  #2.b) Simulation of total area burned based on estimated precipitation regime
  
  #Simulate total area burned
  sim_high_total_burned <- exp(predict(mod_total, newdata = (data.frame("cumu" = sim_cumu_high))) + rnorm(length(sim_cumu_high), 0, sd(residuals(mod_total))))
  sim_low_total_burned <- exp(predict(mod_total, newdata = (data.frame("cumu" = sim_cumu_low))) + rnorm(length(sim_cumu_low), 0, sd(residuals(mod_total))))
  
  # 2.c) Initial parameter settings and definitions for current simulation run
  
  #store fire maps and summary stats annually
  tsf.df <- data.frame("tsf.ini" = tsf.ts[[28]][]) # store annual tsf
  fs.df <- data.frame(fire.scar.ts[]) # store annual fire scar
  fires <- data.frame()   #fire sizes
  
  #initial tsf map
  tsf_cur <- tsf.df[,1]
  
  # counter for fires that couldnt fully burn
  not_complete <- 0
  
  #is a fire restarted? Default F (important for rescheduling of a fire that couldnt reach its size)
  restart = F
  
  #start year
  yr <- 2016
  
  ######################
  #3.) TIME STEP LEVEL##
  ######################
  
  for (i in 1:simle){
    
    #3.a) Sample if fire will occur during this time step
    occurance <- runif(1, 0, 1) < length(unique(all$year))/(2015-1970)
    
    if(occurance == F){
      target_year <- 0
    }
    
    #3.b) Set initial time step parameters
    
    # i) General time step definitions
    precip <- sim_cumu_high[i] #draw precipitation value (for storage later on)
    fire_scar <- values(maskrast-1) #fire scar
    total_fire <- 0 # initial realised total burned area
    
    #ii) If fires occur in this year
    if(occurance == T){
      
      #Calculation of fire probability map
      firenumber <- rowSums(fs.df[]) #determine number of fires on each cell
      no_years <- as.numeric(tail(substr(names(fs.df), 5, 8),1)) - as.numeric(head(substr(names(fs.df), 5, 8),1))
      mfi <- no_years/firenumber #cell level fire interval, produces inf, but 1/inf is 0
      p <- tsf_cur/mfi #probabilities
      p[p > 1] <- 1
      
      #Sample annual area burned
      target_year <- sim_high_total_burned[i]/100 * sum(na.omit(cell_size[]))
      if (target_year > sum(na.omit(cell_size[])) * 60 /100){
        target_year <- sum(na.omit(cell_size[])) * 60 /100 # limit to empiric 60% of study area
      }
      
      #Sample time step fire size distribitution from empiric distriubtions
      t <- sample(1:length(lnmeans), 1) # choose distribution
      ml <- lnmeans[t] #mean log
      sdl <- lnscales[t] # sd log
      target_distr <- rlnorm(10000,ml, sdl) # sample target distr
      target_distr_sub <- target_distr[target_distr < max(all$area) & target_distr > min(all$area)] #limit distribution to min max (will move the mean by a little bit)
    }
    
    #####################
    #4.) FIRE LEVEL######
    #####################
    while(total_fire < target_year){
      
      #4.a) Storage of probability map, fire scar map and total achieved fire size before a new fire is started
      p_old <- p
      fire_scar_old <- fire_scar
      total_fire_old <- total_fire
      
      #4.b) Sample initial fire level parameters
      
      #Sample fire size
      if (restart == F){
        target_fire <- sample(target_distr_sub, size = 1)
      }
      
      #Sample ignition cell
      ini_pburn <- 0 #sum of fire probabilities in adjacent cells, ignition cell only selected when it is > 0 (at least one burnable cell)
      while(ini_pburn == 0){
        
        if (restart == F){ #Ignition probabilty based on probability map
          ic <- sample(inds[which(!is.na(p))], size = 1 , prob = p[which(!is.na(p))])
        }
        
        if (restart == T){ #when a fire is resampled, all unburned cells
          ic <- sample(inds[which(!is.na(p))], size = 1)
          restart = F
        }
        
        sc <- adjacent(maskrast, cells = ic, pairs = F, directions = 8)
        ini_pburn <- sum(na.omit(p[sc]))
      }
      
      ac <- ic
      
      fire_scar[ac] <- 2
      p[ac] <- 0
      
      
      fire_achieved <- cell_size[ac]
      total_fire <- total_fire + cell_size[ac]
      rc <- numeric()
      
      #iteration counter
      it <- 0
      
      ##fire loop
      while(fire_achieved < target_fire & total_fire < target_year){
        if (it > 0){
          sc <- adjacent(maskrast, cells = ac, pairs = F, directions = 8)
        }
        
        #filter cells inside study area as spreading cells
        sc <- sc[!is.na(study_a[sc])]
        
        #burn cells based on mfi/tsf ratio,e.g. cells for which tsf approaches mfi are most likelz to burn
        burn <- ifelse(runif(length(sc),0,1) < p[sc], 1, 0)
        ac <- sc[which(burn == 1)]
        
        #store burning front for potential restart of the fire
        rc <- c(rc, sc[which(burn == 0 & p[sc] > 0)])
        
        #restart fire if prematurely off, based on p
        #will be the case when: fire did not spread to spreading cells although it can
        if (length(ac) == 0 && length(sc) > 1 && sum(p[sc]) > 0){
          ac <- sample(sc, size = 1, p = p[sc])
        }
        
        # fire did not spread to spreading cells because none of them can burn: restart from rc when rc > 1
        if(length(ac) == 0 && length(sc) > 1 && sum(p[sc]) == 0 && sum(p[rc]) > 0 && length(rc) > 1){
          ac <- sample(rc, size = 1, p = p[rc])
        }
        # fire did not spread to spreading cells, because none of them can burn: restart from rc if rc = 1
        if(length(ac) == 0 && length(sc) > 1 && sum(p[sc]) == 0 && sum(p[rc]) > 0 && length(rc) == 1){
          ac <- rc
        }
        
        #fire did not spread to spreading cells, because none of them can burn and none of the rc can burn: fire stops and new fire is started
        if(length(ac) == 0 && length(sc) > 1 && sum(p[sc]) == 0 && sum(p[rc]) == 0){
          not_complete <- not_complete + 1
          total_fire <- total_fire_old
          fire_achieved <- 0
          p <- p_old
          fire_scar <- fire_scar_old
          restart = T
          break
        }
        
        #number of cells to burn till target fire
        diff <- target_fire - fire_achieved
        no <- ceiling(diff/mean(cell_size[sc]))
        
        #randomly choose current active cells, to which fire will not spread
        if(no < length(ac)){
          ac <- sample(ac, size = no, prob = 1-p[ac] + 0.00001, replace = T)
        }
        
        fire_achieved <- fire_achieved + sum(cell_size[ac])
        total_fire <- total_fire + sum(cell_size[ac])
        
        #update fire prob and fire scar maps
        p[ac] <- 0
        fire_scar[ac] <- 1
        
        #increment iteration counter
        it <- it + 1
      }
      #store fires
      if(restart == F){
        fires <- rbind(fires, data.frame(yr, fire_achieved, target_fire, precip))
      }
      cat(yr, round(total_fire/target_year * 100, 2),"%","target:", target_year, "mean size:", mean(target_distr_sub), not_complete, "\r")
    }
    
    #update fs
    fs.df <- cbind(fs.df, fire_scar)
    names(fs.df)[which(names(fs.df) == "fire_scar")] <- paste0("fs_s", yr)
    
    #update tsf
    tsf_cur <- tsf_cur +1
    tsf_cur[which(fire_scar == 1)] <- 1
    tsf.df <- cbind(tsf.df, tsf_cur)
    
    #increment year
    yr <- yr + 1
  }
  
  sim_data[[j]] <- list(fires, tsf.df, fs.df)
}




plot(sim_data[[1]][[1]]$fire_achieved ~ sim_data[[1]][[1]]$target_fire)
head(fs.df)
str(sim_data)
plot(fires$fire_achieved ~ fires$target_fire)
test <- maskrast
values(test) <- tsf.df[,30]
plot(test)
length(fire_scar[which(fire_scar%in%c(1,2))])/length(fire_scar[!is.na(fire_scar)])
sum(na.omit(cell_size[which(fire_scar%in%c(1))]))/sum(na.omit(cell_size))

total_fire/target_year
length(fire_scar)
values(test) <- fs.df[,30]
head(fs.df)
plot(test)
head(tsf.df)
mean(fires$fire_achieved[which(fires$year == 2059)])
#writeRaster(paste0(output,"/", tsf, sim, year), format = "ascii")
fire_scar.map[which(fire_scar.map[] == 0)] <- NA
plot(fire_scar.ts[[25]])
plot(fire_scar.map, col = c("red", "black"), add = T, legend = F)
length(tsf.df)
