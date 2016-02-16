require(msm)
require(raster)
require(rgdal)

###################################
#1.) LOADING ALL NECESSARY DATA####
###################################

#1.a) Paths and file lists
path_climate <- "/Users/Simon/Studium/MSC/Masterarbeit/data/climate/BOM/"
path_shp <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/fire_data_reprojected_1970-2015"
path_fire_rasters <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/fire rasterized 1970-2015"
path_tsf_rasters <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/tsf rasterized 1970-2015"
files_shp <- list.files(path_shp, pattern = "*.shp$")
files_tsf <- list.files(path_tsf_rasters, pattern = "*.asc$")
files_fs <- list.files(path_fire_rasters, pattern = "*.asc$")

#1.b) LOAD RASTERS

#i) Fire scar data

fire.scar.ts <- stack()
for (i in 1:length(files_fs)){
  r <- raster(paste(path_fire_rasters, files_fs[i], sep ="/"))
  crs(r) <- "+proj=longlat +ellps=GRS80 +no_defs"
  fire.scar.ts <- stack(fire.scar.ts, r)
}

#ii) Time since fire (tsf) data
tsf.ts <- stack()
for (i in 1:length(files_tsf)){
  r <- raster(paste(path_tsf_rasters, files_tsf[i], sep ="/"))
  crs(r) <- "+proj=longlat +ellps=GRS80 +no_defs"
  tsf.ts <- stack(tsf.ts, r)
}

#1.b) Combine fire and precipitation statistics

#i) Read fire sizes from shape files and write into data frame

all <- data.frame()
for (i in 1:length(files_shp)){
  
  #read shapefiles
  shp <- readOGR(path_shp, substr(files_shp[i], 1, 9))
  names(shp)[1] <- c("OBJECTID")
  
  #extract data
  f <- data.frame("area" = shp@data$Shape_Area, 
                  "year" = shp@data$CalanderYe, 
                  "no_of_fires" = length(shp@data$OBJECTID))
  all <- rbind(all, f)
}

#ii) read and process precipitation data and attach to data frame

pre <- read.csv(paste0(path_climate, "IDCJAC0009_015611_1800_Data.csv"))
prec.stats <- aggregate(pre$Rainfall.amount..millimetres., by = list(pre$Year), FUN = sum, na.rm = T)

prec.stats$yycumu <- NA
for (i in 3:nrow(prec.stats)){
  prec.stats$yycumu[i] <-  prec.stats$x[i] + prec.stats$x[i-1] + prec.stats$x[i-2] 
}
names(prec.stats) <- c("year", "prec", "yycumu")
pos <- match(all$year, prec.stats$year)
all$cumu_prec <- prec.stats$yycumu[pos]

#1.c) Find relationship between precipitation and total area burned
attach(all)
wide_area <- 5503.63 #the total area of available data

#i) Total area burned per year [%]
total <- round((tapply(area, year, function(x) sum(x)))/wide_area * 100, digits = 3)

#ii) Cumulative Precipitation per year
cumu <- tapply(cumu_prec, year, function(x) unique(x))

#iii) Model
mod_total <- lm(log(total) ~ cumu)

#1.d)

#i) Subset data
area.subs <- area[which(year%in%c(1976, 1982,1983, 1984, 2000:2012))]
years.subs <- year[which(year%in%c(1976, 1982,1983, 1984, 2000:2012))]
cumu.subs <- unique(cumu_prec[which(year%in%c(1976, 1982,1983, 1984, 2000:2012))])

#ii) Determine unfitted lnorm pars
lnscales <- tapply(area.subs, years.subs, function(x) lnscale(x))
lnmeans <- tapply(area.subs, years.subs, function(x) lnmean(x))


#iii) Models
lnmeans_mod <- lm(lnmeans ~ cumu.subs)
lnscales_mod <- lm(lnscales ~ cumu.subs)

########################################################
#1.) Preprocessing of data and simulation definitions####
########################################################
detach(all)
attach(prec.stats)

#1.a) Preprocessing precipitation reference data

#87 - 2005 reference
reference_prec <- prec[year >= 1987 & year <= 2005]
mean_prec <- mean(reference_prec)
high_prec_30 <- mean_prec *108/100
low_prec_30 <- mean_prec * 89/100
high_prec_90 <- mean_prec *119/100
low_prec_90 <- mean_prec * 69/100
neut_prec_30 <- mean_prec

#calculation of annual increments/decrements for high and low scenario
years_30 <- 2030 - 2015
years_90 <- 2090 -2030
incr_15_30 <- (high_prec_30 - mean_prec) / years_30
decr_15_30 <- (low_prec_30 - mean_prec) / years_30
incr_30_90 <- (high_prec_90 - high_prec_30) / years_90
decr_30_90 <- (low_prec_90 - low_prec_30) / years_90

incr <- c(rep(incr_15_30, times = years_30), rep(incr_30_90, times = years_90))
decr <- c(rep(decr_15_30, times = years_30), rep(decr_30_90, times = years_90))
neut <- rep(0, 75)

#1.b) Initial parameter values and definitions for entire simulation
scen <- "high"
simle <- 75 #simulation length
simno <- 1 #number of simulations
sim_data <- list() #data frame for all sim data
output <- "/Users/Simon/Studium/MSC/Masterarbeit/data/new_model_out"
cell_size <- raster("/Users/Simon/Studium/MSC/Masterarbeit/data/cell_size.asc")
cell_size <- cell_size[] #cell sizes in km^2
maskrast <- raster("/Users/Simon/Studium/MSC/Masterarbeit/data/maskraster.asc")
study_a <- maskrast[] #study area
inds <- 1:length(maskrast[]) #cell indices
elev <- raster(paste0("/Users/Simon/Studium/MSC/Masterarbeit/data//Elevation/elevation9secNH.asc"))
elev <- elev[]
ncol <- ncol(maskrast)
nrow <- nrow(maskrast)

#########################
#2.) SIMULATION LEVEL#####
#########################
for (j in 1:simno){
  
  #2.a) Simulation of precipitation regime for time span
  
  #Estimation of simulated precipitation means 2015 - 2090 (first value for 2016), with random terms
  if (scen == "high"){
    change <- incr
  }else if (scen == "low"){
    change <- decr
  }else if (scen == "neut"){
    change <- neut
  }
  
  sim_prec <- numeric()
  for (i in 1:length(change)){
    sim_prec[i] <- exp(log(mean_prec + sum(change[1:i])) + rnorm(1, 0, sd(log(reference_prec))))
  }
  
  #Attach last two measured precipitation values (2014 and 2015) for calculation of first time step
  sim_prec <- c(tail(prec, 2), sim_prec)
  
  #Calculate simulated cumulative precipitation
  sim_cumu <- numeric()
  for (i in 1:(length(sim_prec)-2)){
    sim_cumu[i] <-  sim_prec[i+2] + sim_prec[i+1] + sim_prec[i]
  }
  
  #2.b) Simulation of total area burned based on estimated precipitation regime
  
  #Simulate total area burned
  sim_total_burned <- exp(predict(mod_total, newdata = (data.frame("cumu" = sim_cumu))) + rnorm(length(sim_cumu), 0, sd(residuals(mod_total))))
  
  # 2.c) Initial parameter settings and definitions for current simulation run
  
  #store fire maps and summary stats annually
  tsf.df <- data.frame("tsf.ini" = tsf.ts[[28]][]) #store annual tsf
  fs.df <- data.frame(fire.scar.ts[]) #store annual fire scar
  fires <- data.frame() #fire sizes
  
  #initial tsf map
  tsf_cur <- tsf.df[,1]
  
  # counter for fires that couldn't fully burn
  not_complete <- 0
  
  #is a fire restarted? Default F (important for rescheduling of a fire that couldnt reach its size)
  restart = F
  
  #Start year
  yr <- 2016
  
  ######################
  #3.) TIME STEP LEVEL##
  ######################
  
  for (i in 1:simle){
    
    #3.a) Sample if fire will occur during this time step (based on historic fire frequency)
    occurance <- runif(1, 0, 1) < length(unique(all$year))/(2015-1970)
    
    if(occurance == F){
      target_year <- 0
    }
    
    #3.b) Set initial time step parameters
    
    # i) General time step definitions
    precip <- sim_cumu[i] #draw precipitation value (for storage later on)
    fire_scar <- maskrast[]-1 #fire scar
    total_fire <- 0 # initial realised total burned area
    
    #ii) If fires occur in this year
    if(occurance == T){
      
      #Calculation of fire probability map
      no_fires <- rowSums(fs.df[]) #determine number of fires on each cell
      no_years <- as.numeric(tail(substr(names(fs.df), 5, 8),1)) - as.numeric(head(substr(names(fs.df), 5, 8),1))
      mfi <- no_years/no_fires #cell level fire interval, produces inf, but 1/inf is 0
      p <- tsf_cur/mfi #probabilities
      p[p > 1] <- 1
      
      #Sample annual area burned
      target_year <- sim_total_burned[i]/100 * sum(na.omit(cell_size[]))
      if (target_year > sum(na.omit(cell_size[])) * 60 /100){
        target_year <- sum(na.omit(cell_size[])) * 60 /100 # limit to empiric 60% of study area
      }
      
      #Sample time step fire size distribitution from empiric distriubtions
      loc <- predict(lnmeans_mod, newdata = data.frame(cumu.subs = precip)) + rnorm(1, 0, sd(residuals(lnmeans_mod)))
      scale <- predict(lnscales_mod, newdata = data.frame(cumu.subs = precip)) + rnorm(1, 0, sd(residuals(lnscales_mod)))
      target_distr <- rlnorm(10000,loc, scale)
      target_distr_sub <- target_distr[target_distr < max(area.subs) & target_distr > min(area.subs)] #limit distribution to min max (will move the mean by a little bit)
    }
    
    ####################
    #4.) FIRE LEVEL#####
    ####################
    while(total_fire < target_year){
      #4.a) Storage of probability map, fire scar map and total achieved fire size before a new fire is started
      p_old <- p
      fire_scar_old <- fire_scar
      total_fire_old <- total_fire
      
      #4.b) Sample fire size
      if (restart == F){
        target_fire <- sample(target_distr_sub, size = 1)
      }
      
      if (restart == T){
        restart <- F
      }
      
      #4.b) Sample ignition cell
      ini_pburn <- 0 #sum of fire probabilities in adjacent cells, ignition cell only selected when it is > 0 (at least one burnable cell)
      while(ini_pburn == 0){
        
        ic <- sample(inds[which(!is.na(p))], size = 1 , prob = p[which(!is.na(p))]) #Ignition probabilty based on probability map
        pairs <- adjacent(maskrast, cells = ic, pairs = T, directions = 8) #test if adjacent cells are burnable
        ini_pburn <- sum(na.omit(p[pairs[,2]]))
      }
      
      #Update fire scar, probability map achieved fire size and achieved total fire size
      fire_scar[ic] <- 1
      p[ic] <- 0
      fire_achieved <- cell_size[ic]
      total_fire <- total_fire + cell_size[ic]
      
      #4.c) Other initial fire level parameters and definitions
      ac <- ic #Ignition cell becomes first active cell
      rc <- numeric() # vector for storage of fire front
      
      #####################
      #5.) CELL LEVEL######
      #####################
      while(fire_achieved < target_fire && total_fire < target_year){
        
        #5.a) Determine spreading cells
        
        #Identify all cells that are adjacent to current active cells
        if (fire_achieved > cell_size[ic]){
          pairs <- adjacent(maskrast, cells = ac, pairs = T, directions = 8)
        }
        
        #Filter cells in study area
        pairs <- pairs[!is.na(study_a[pairs[,2]]),]
        sc <- pairs[,2]
        
        #5.b) Burn cells
        
        #adjust for distance (diagonal cells are further away, thus less likely to burn)
        dists_a <- dist(from = pairs[,1], to = pairs[,2], ncol = ncol, nrow = nrow)
        d <- min(dists_a)/dists_a
        
        #adjust for slope
        elevDiff <- (elev[pairs[,1]] - elev[pairs[,2]])/1000
        slopes <- elevDiff/dists_a
        slopes <- clamp(slopes, -1, 1)
        s <- 1 + slopes
        
        #Based on mfi/tsf ratio, e.g. cells for which tsf approaches mfi are most likely to burn
        burn <- ifelse(runif(length(sc),0,1) < (p[sc] * d * s), 1, 0)
        ac <- sc[which(burn == 1)] # burned cells become new active cells
        
        #Store burning front for potential restart of the fire
        rc <- c(rc, sc[which(burn == 0 & p[sc] > 0)])
        
        #5.c) Restart a fire if no spreading cells burned:
        
        #Fire did not spread to spreading cells although they can potentially burn: Restart from current spreading cells
        if (length(ac) == 0 && length(sc) > 1 && sum(p[sc]) > 0){
          ac <- sample(sc, size = 1, p = p[sc])
        }
        
        #Fire did not spread to spreading cells, because none of them can burn: restart from fire front (rc) when length(rc) > 1
        if(length(ac) == 0 && length(sc) > 1 && sum(p[sc]) == 0 && sum(p[rc]) > 0){
          if(length(rc) > 1){
            ac <- sample(rc, size = 1, p = p[rc])
          }else if(length(rc) == 1){
            ac <- rc
          }
        }
        
        #fire did not spread to spreading cells, because none of them can burn and fire front can't burn: fire stops and new fire is started
        if(length(ac) == 0 && length(sc) > 1 && sum(p[sc]) == 0 && sum(p[rc]) == 0){
          not_complete <- not_complete + 1
          total_fire <- total_fire_old
          fire_achieved <- 0
          p <- p_old
          fire_scar <- fire_scar_old
          restart = T
          break
        }
        
        #5.d) Decrement currently active cells to stay within target fire size
        
        #area to burn till target fire size is reached
        diff <- ifelse(target_year - total_fire < target_fire - fire_achieved, target_year - total_fire,target_fire - fire_achieved) # how much area is left to burn
        no <- ceiling(diff/mean(cell_size[ac])) # approximitely how many cells will that make
        
        #Choose current active cells, to which fire will not spread, based on p (+ 0.00001 to assure sufficient positive probabilities)
        if(no < length(ac)){
          ac <- sample(ac, size = no, prob = 1-p[ac] + 0.00001, replace = T)
        }
        
        #5.e) Update fire parameters after each itereation 
        fire_scar[ac] <- 1
        p[ac] <- 0
        fire_achieved <- sum(cell_size[which(fire_scar - fire_scar_old == 1)])
        total_fire <- sum(cell_size[which(fire_scar == 1)])
        cat(yr, round(total_fire/target_year * 100, 2),"%","target:", target_year, not_complete, "\r")
      }
      
      ##########################
      #6.)STORAGE OF SIM DATA###
      ##########################
      
      #6.a) fire level data
      if(fire_achieved != 0){
        fires <- rbind(fires, data.frame(yr, fire_achieved, target_fire, precip, total_fire, target_year))
      }
    }
    
    #6.b) Time-step level data
    
    #fire scar
    fs.df <- cbind(fs.df, fire_scar)
    names(fs.df)[which(names(fs.df) == "fire_scar")] <- paste0("fs_s", yr)
    
    #time since fire
    tsf_cur <- tsf_cur +1
    tsf_cur[which(fire_scar == 1)] <- 1
    tsf.df <- cbind(tsf.df, tsf_cur)
    names(tsf.df)[which(names(tsf.df) == "tsf_cur")] <- paste0("tsf_s", yr)
    
    #increment year
    yr <- yr + 1
  }
  
  #6.c) Simulation run level data
  sim_data[[j]] <- list("fires" = fires, "tsf" = tsf.df, "fs" = fs.df)
}

test <- maskrast
i = 1
values(test) <- tsf.df[,i]
plot(test)
i <- i + 1
