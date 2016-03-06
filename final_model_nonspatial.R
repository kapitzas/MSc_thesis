#Output path
output <- "/Users/Simon/Studium/MSC/Masterarbeit/data/non_spatial_new_model_out" #output path for fire record and rasters
test <- F
#Simulation parameters
scenario <- c("low", "high", "neut") #scenario
simle <- 120 #simulation length (cannot be changed!)
simno <- 100 #number of simulations
ini_year <- 1971

#######################
#1.) SCENARIO LEVEL####
#######################
for(k in 1:length(scenario)){
  #1.a) Global scenario variables
  scen <- scenario[k]
  sim_no <- 1 #Simulation counter
  fires_sim_all <- matrix(ncol = 4, nrow = 0) #matrix for fire level fire record, will be written for each scenario
  tab_sim_all <- matrix(ncol = 4, nrow = 0) #matrix for year level fire record
  
  ##############################
  #2.) FOOTPRINT SIMULATIONS####
  ##############################
  for (j in 1:simno){
    
    #2.a) Simulation of precipitation time series
    
    if (scen == "high"){
      change <- c(rep(0, 45), incr) #add observation time span for model validation
    }else if (scen == "low"){
      change <- c(rep(0, 45),decr)
    }else if (scen == "neut"){
      change <- c(rep(0, 45),neut)
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
    for (i in 1:(length(sim_prec))-2){
      sim_cumu[i] <-  sim_prec[i+2] + sim_prec[i+1] + sim_prec[i]
    }
    
    #2.b) Simulation of total area burned based on estimated precipitation regime
    sim_tab <- numeric()
    sim_loc <- numeric()
    sim_scale <- numeric()
    for (i in 1:length(sim_cumu)){
      sim_tab[i] <- (rtnorm(1, predict(tab_mod, newdata = (data.frame("cumu3" = sim_cumu[i]))), sd(residuals(tab_mod)), lower = sqrt(min(tab)), upper = sqrt(max(tab))))^2
      sim_loc[i] <- rtnorm(1, predict(locs_mod, newdata = (data.frame("cumu3.subs" = sim_cumu[i]))), sd(residuals(locs_mod)), lower = min(locs), upper = max(locs))
      sim_scale[i] <- rtnorm(1, predict(scales_mod, newdata = (data.frame("cumu3.subs" = sim_cumu[i]))), sd(residuals(scales_mod)), lower = min(scales), upper = max(scales))
    }
    
    occu <- round(predict(fo_mod, newdata = data.frame("pr" = sim_cumu), type = "response"),3)
    
    yr <- ini_year #Start year
    fires <- matrix(ncol = 4, nrow = 0) ##Storage of fire record
    tab <- matrix(ncol = 4, nrow = 0)
    
    #TIME STEP LEVEL
    for (i in 1:simle){
      precip <- sim_cumu[i] #draw precipitation value (for storage later on)
      occurrence <- runif(1, 0, 1) < occu[i] #will fire occur?
      tab_achieved <- 0 #initial realised total burned area  
      if(occurrence == F){
        tab_target <- 0 #total area burned
      }
      
      if(occurrence == T){
        
        #ii) Sample annual area burned
        tab_target <- sim_tab[i]/100 * sum(na.omit(cell_size[])) 
        
        #iii) Sample time step fire size distribitution
        fire_distr <- rlnorm(10000,sim_loc[i], sim_scale[i])
        fire_distr_sub <- fire_distr[fire_distr < max(area.subs) & fire_distr > min(area.subs)]
        }
      
      ####################
      #4.) FIRE LEVEL#####
      ####################      
      while(tab_achieved < tab_target){
        fire_target <- sample(fire_distr_sub, size = 1)
        
        if(tab_achieved + fire_target > tab_target){
          fire_target <- tab_target - tab_achieved
        }
        tab_achieved <- tab_achieved + fire_target
        fires <- rbind(fires, c(sim_no, yr, precip, fire_target))
        cat(yr, tab_target, tab_achieved, "Scenario:",scen, ";", "Simulation run:", sim_no, "\r")
      }
      if(tab_target != 0){
        tab <- rbind(tab, c(sim_no, yr, precip, tab_target))
      }
        yr <- yr + 1
    }
    
    #STORAGE
    fires_all <- rbind(fires_all, fires)
    tab_all <- rbind(tab_all, tab)
    sim_no <- sim_no + 1
  }
  fires_all <- as.data.frame(fires_all)
  tab_all <- as.data.frame(tab_all)
  colnames(fires_all) <- c("sim_no", "yr", "precip", "fire_target")
  colnames(tab_all) <- c("sim_no", "yr", "precip", "tab_target")
  if (test == F){
    write.csv(fires_all, paste0(output,"/", "ns_fires_",scen,".csv"), row.names = F)
    write.csv(tab_all, paste0(output,"/", "ns_tab_",scen,".csv"), row.names = F)  
  }
  if (test == T){
    write.csv(fires_all, paste0(output,"/", "ns_fires_",scen,"_test",".csv"), row.names = F)
    write.csv(tab_all, paste0(output,"/", "ns_tab_",scen,"_test",".csv"), row.names = F)  
  }
}
