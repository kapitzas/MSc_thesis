require(msm)
require(raster)
require(rgdal)

#######################
#1.) SCENARIO LEVEL####
#######################
for(k in 1:length(scenario)){
  
  #1.a) Global scenario variables
  scen <- scenario[k]
  sim_no <- 1 #Simulation counter
  fires_all <- matrix(ncol = 4, nrow = 0) #dataframe for fire level fire record, will be written for each scenario
  tab_all <- matrix(ncol = 4, nrow = 0) #dataframe for year level fire record
  
  ##########################
  #2.) TIME SERIES LEVEL####
  ##########################
  for (j in 1:simno){
    
    #2.a) Simulation of precipitation regime for time span
    
    if (scen == "high"){
      change <- incr
    }else if (scen == "low"){
      change <- decr
    }else if (scen == "neut"){
      change <- neut
    }
    
    #i) Estimation of bimodal means based on projected annual change
    sim_prec <- numeric()
    for (i in 1:length(change)){
      data <- rnorm(n=10000,mean=mus[components] + sum(change[1:i]),sd=sds[components])
      sim_prec[i] <- sample(data, 1)
      while(sim_prec[i] < 0){
        sim_prec[i] <- sample(data, 1)
      }
    }
    # 
    # for (i in 1:length(change)){
    #   data <- rnorm(n=10000, mean= mean_prec + sum(change[1:i]),sd=sd(reference_prec))
    #   sim_prec[i] <- sample(data, 1)
    #   while(sim_prec[i] < 0){
    #     sim_prec[i] <- sample(data, 1)
    #   }
    # }
    
    #ii) Attach last two measured precipitation values (2014 and 2015) for calculation of first time step
    sim_prec <- c(tail(prec.stats$prec, 2), sim_prec)
    
    #iii) Cumulative precipitation
    sim_cumu <- numeric()
    for (i in 1:(length(sim_prec)-2)){
      sim_cumu[i] <-  sim_prec[i+2] + sim_prec[i+1] + sim_prec[i]
    }
    
    #2.b) Simulation of total area burned based on estimated precipitation regime
    sim_tab <- numeric()
    loc <- numeric()
    scale <- numeric()
    for (i in 1:length(sim_cumu)){
      sim_tab[i] <- exp(rtnorm(1, predict(mod_total, newdata = (data.frame("cumu" = sim_cumu[i]))), sd(residuals(mod_total)),lower = log(min(total)), upper = log(max(total))))
      loc[i] <- rtnorm(1, predict(lnmeans_mod, newdata = (data.frame("cumu" = sim_cumu[i]))), sd(residuals(lnmeans_mod)), lower = min(lnmeans), upper = max(lnmeans))
      scale[i] <- rtnorm(1, predict(lnscales_mod, newdata = (data.frame("cumu" = sim_cumu[i]))), sd(residuals(lnscales_mod)), lower = min(lnscales), upper = max(lnscales))
    }
    
    #2.d) Simulation of fire occurrence
    occu <- round(predict(mod_occ, newdata = data.frame("pr" = sim_cumu), type = "response"),3)
    
    # 2.e) Initialisation of dynamic time series variables
    
    #i) Spatial variables
    tsf.df <- data.frame("tsf.ini" = tsf.ts[[28]][]) #Storage annual tsf
    tsf_cur <- tsf.df[,1] #Initial TSF map
    fs.df <- data.frame(fire.scar.ts[]) #Store annual fire scar
    nf_cur <- rowSums(fs.df[,3:46])
    nf.df <- data.frame("nf.ini" = nf_cur)
    
    #ii) Other time series variables
    
    yr <- ini_year #Start year
    fires <- matrix(ncol = 4, nrow = 0) ##Storage of fire record
    tab <- matrix(ncol = 4, nrow = 0)
    
    ########################
    #3.) TIME STEP LEVEL####
    ########################
    
    for (i in 1:simle){
      
      #3.a) Initialisation of time step variables
      restart = F #no fire is initially restarted
      p <- study_a-1 #reset p map
      precip <- sim_cumu[i] #draw precipitation value (for storage later on)
      fire_scar <- maskrast[]-1 #fire scar
      tab_achieved <- 0 # initial realised total burned area  
      occurrence <- runif(1, 0, 1) < occu[i]  #will fire occur?
      if(occurrence == F){
        tab_target <- 0 #total area burned
      }
      
      #3.b) Initialisation of time step variables if fire occurs
      
      if(occurrence == T){
        
        #i) Predict burn probabilities
        for (l in 1: length(total_class)){
          p[which(pyric_class == l)] <- round(predict(mod_probs[[l]],newdata = data.frame("yr" = tsf_cur[which(pyric_class == l)])),3)
        }
        p[which(p > 1)] <- 1
        
        #ii) Sample annual area burned
        tab_target <- round(sim_tab[i]/100 * sum(na.omit(cell_size[])),5)
        
        #iii) Sample time step fire size distribitution
        fire_distr <- rlnorm(10000,loc[i], scale[i])
        fire_distr_sub <- round(fire_distr[fire_distr < max(area.subs) & fire_distr > min(area.subs)],5) #limit distribution to min max (will move the mean by a little bit)
      }
      
      ####################
      #4.) FIRE LEVEL#####
      ####################
      while(tab_achieved < tab_target){
        #4.a) Storage of probability map, fire scar map and total achieved fire size before a new fire is started
        p_old <- p
        fire_scar_old <- fire_scar
        tab_achieved_old <- tab_achieved
        
        #4.b) Sample fire size
        if (restart == F){
          fire_target <- sample(fire_distr_sub, size = 1)
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
        
        fire_scar[ic] <- 1 #Update fire scar
        p[ic] <- 0 #Update probability map
        fire_achieved <- cell_size[ic] #Update fire achieved
        tab_achieved <- tab_achieved + cell_size[ic] #update total fire
        
        #4.c) Initialisation of fire level variables
        ac <- ic #Ignition cell becomes first active cell
        rc <- numeric() # vector for storage of fire front
        
        #####################
        #5.) CELL LEVEL######
        #####################
        while(fire_achieved < fire_target && tab_achieved < tab_target){
          
          #5.a) Determine spreading cells
          
          #Identify all cells that are adjacent to current active cells
          if (fire_achieved > cell_size[ic]){
            pairs <- adjacent(maskrast, cells = ac, pairs = T, directions = 8)
          }
          
          #Filter cells in study area
          pairs <- pairs[!is.na(study_a[pairs[,2]]),]
          sc <- pairs[,2]
          
          #5.b) Burn cells
          
          #i) adjust for distance (diagonal cells are further away, thus less likely to burn)
          dists_a <- round(dist(from = pairs[,1], to = pairs[,2], ncol = ncol, nrow = nrow),3)
          d <- min(dists_a)/dists_a
          
          #ii) adjust for slope
          elevDiff <- round((elev[pairs[,1]] - elev[pairs[,2]])/1000,3)
          slopes <- round(elevDiff/dists_a,4)
          slopes <- clamp(slopes, -1, 1)
          s <- 1 + slopes
          
          #iii) Based on mfi/tsf ratio, e.g. cells for which tsf approaches mfi are most likely to burn
          burn <- ifelse(runif(length(sc),0,1) < (p[sc] * d * s), 1, 0)
          ac <- sc[which(burn == 1)] # burned cells become new active cells
          
          #ix) Store burning front for potential restart of the fire
          rc <- c(rc, sc[which(burn == 0 & p[sc] > 0)])
          
          #5.c) Restart a fire if no spreading cells burned:
          
          #i) Fire did not spread to spreading cells although they can potentially burn: Restart from current spreading cells
          if (length(ac) == 0 && length(sc) > 1 && sum(p[sc]) > 0){
            ac <- sample(sc, size = 1, p = p[sc])
          }
          
          #ii) Fire did not spread to spreading cells, because none of them can burn: restart from fire front (rc) when length(rc) > 1
          if(length(ac) == 0 && length(sc) > 1 && sum(p[sc]) == 0 && sum(p[rc]) > 0){
            if(length(rc) > 1){
              ac <- sample(rc, size = 1, p = p[rc])
            }else if(length(rc) == 1){
              ac <- rc
            }
          }
          
          #iii) fire did not spread to spreading cells, because none of them can burn and fire front can't burn: fire stops and new fire is started
          if(length(ac) == 0 && length(sc) > 1 && sum(p[sc]) == 0 && sum(p[rc]) == 0){
            tab_achieved <- tab_achieved_old
            fire_achieved <- 0
            p <- p_old
            fire_scar <- fire_scar_old
            restart = T
            break
          }
          
          #5.d) Decrement currently active cells to stay within target fire size
          
          diff <- ifelse(tab_target - tab_achieved < fire_target - fire_achieved, tab_target - tab_achieved,fire_target - fire_achieved) # how much area is left to burn
          no <- ceiling(diff/mean(cell_size[ac])) # approximitely how many cells will that make
          
          if(no < length(ac)){ #randomly extinguish current active cells
            ac <- sample(ac, size = no, prob = p[ac], replace = T)
          }
          
          #5.e) Update time step variables after each itereation 
          fire_scar[ac] <- 1
          p[ac] <- 0
          fire_achieved <- sum(cell_size[which(fire_scar - fire_scar_old == 1)])
          tab_achieved <- sum(cell_size[which(fire_scar == 1)])
          cat(yr, round(tab_achieved/tab_target * 100, 2),"%","target:", tab_target, "Scenario:",scen, ";", "Simulation run:", sim_no, "\r")
        }
        
        ##################################
        #6.)STORE AND UPDATE VARIABLES####
        ##################################
        
        #6.a) fire level data
        if(fire_achieved != 0){
          fires <- rbind(fires, c(sim_no, yr, precip, fire_achieved, fire_target))
        }
      }
      
      #6.b) Time series variables
      if(tab_target != 0){
        tab <- rbind(tab, c(sim_no, yr, precip, tab_achieved ,tab_target))
      }
      
      #i) Fire scar
      fs.df <- cbind(fs.df, fire_scar)
      names(fs.df)[which(names(fs.df) == "fire_scar")] <- paste0("fs_s", yr)
      
      #ii) Time since fire
      tsf_cur <- tsf_cur +1
      tsf_cur[which(fire_scar == 1)] <- 1
      tsf.df <- cbind(tsf.df, tsf_cur)
      names(tsf.df)[which(names(tsf.df) == "tsf_cur")] <- paste0("tsf_s", yr)
      
      #iii) Fire Number
      nf_cur <- rowSums(fs.df[,(i+3):(i+46)])
      nf.df <- cbind(nf.df, nf_cur)
      names(nf.df)[which(names(nf.df) == "nf_cur")] <- paste0("nf_s", yr)
      
      #iii) Year
      yr <- yr + 1
      
      #6.c) Output tsf raster
      if(length(tsf.df[]) == 2){
        r <- maskrast
        r[] <- tsf.df[,1]
        name_ras <- as.character(paste0("tsf_", scen, "_", sim_no,"_", 0))
        writeRaster(r, paste(output, name_ras, sep = "/"), format = "ascii")  
      }
      
      r <- maskrast
      r[] <- tsf.df[,i + 1]
      name_ras <- as.character(paste0("tsf_",scen, "_", sim_no,"_", i))
      writeRaster(r, paste(output, name_ras, sep = "/"), format = "ascii")
      
      #6.d)
      if(length(nf.df[]) == 2){
      r <- maskrast
      r[] <- nf.df[,1]
      name_ras <- as.character(paste0("nf_",scen, "_", sim_no,"_", 0))
      writeRaster(r, paste(output, name_ras, sep = "/"), format = "ascii")
      }
      
      r <- maskrast
      r[] <- nf.df[,i]
      name_ras <- as.character(paste0("nf_", scen, "_", sim_no,"_", i))
      writeRaster(r, paste(output, name_ras, sep = "/"), format = "ascii")
    }
    #6.d) Update Global variables
    fires_all <- rbind(fires_all, fires)
    sim_no <- sim_no + 1
    tab_all <- rbind(tab_all, tab)
  }
  #6.e) Output Fire data
  fires_all <- as.data.frame(fires_all)
  tab_all <- as.data.frame(tab_all)
  colnames(fires_all) <- c("sim_no", "yr", "precip", "fire_target")
  colnames(tab_all) <- c("sim_no", "yr", "precip", "tab_target")
  write.csv(fires_all, paste0(output,"/", "fires_",scen,".csv"), row.names = F)
  write.csv(tab_all, paste0(output,"/", "tab_",scen,".csv"), row.names = F)
}