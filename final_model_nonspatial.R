#Output path
output <- "/Users/Simon/Studium/MSC/Masterarbeit/data/non_spatial_new_model_out" #output path for fire record and rasters
test <- F
#Simulation parameters
scenario <- c("low", "high", "neut") #scenario
simle <- 120 #simulation length (cannot be changed!)
simno <- 100 #number of simulations
ini_year <- 1971

#SCENARIO LEVEL
for(k in 1:length(scenario)){
  #1.a) Global scenario variables
  scen <- scenario[k]
  sim_no <- 1 #Simulation counter
  fires_all <- matrix(ncol = 4, nrow = 0) #matrix for fire level fire record, will be written for each scenario
  tab_all <- matrix(ncol = 4, nrow = 0) #matrix for year level fire record
  
  #TIME SERIES LEVEL
  for (j in 1:simno){
    
    #Simulation of precipitation regime for time span
    
    if (scen == "high"){
      change <- c(rep(0, 45), incr)
    }else if (scen == "low"){
      change <- c(rep(0, 45),decr)
    }else if (scen == "neut"){
      change <- c(rep(0, 45),neut)
    }
    
    # Estimation of means based on projected annual change
    # i) Estimation of means based on projected annual change
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
    
    sim_prec <- c(tail(prec.stats$prec, 2), sim_prec)
    
    #Cumulative precipitation
    sim_cumu <- numeric()
    for (i in 1:(length(sim_prec))-2){
      sim_cumu[i] <-  sim_prec[i+2] + sim_prec[i+1] + sim_prec[i]
    }
    
    #Simulation of total area burned based on estimated precipitation regime
    sim_tab <- numeric()
    loc <- numeric()
    scale <- numeric()
    for (i in 1:length(sim_cumu)){
      sim_tab[i] <- (rtnorm(1, predict(mod_total, newdata = (data.frame("cumu3" = sim_cumu[i]))), sd(residuals(mod_total)), lower = sqrt(min(total)), upper = sqrt(max(total))))^2
      loc[i] <- rtnorm(1, predict(lnmeans_mod, newdata = (data.frame("cumu3.subs" = sim_cumu[i]))), sd(residuals(lnmeans_mod)), lower = min(lnmeans), upper = max(lnmeans))
      scale[i] <- rtnorm(1, predict(lnscales_mod, newdata = (data.frame("cumu3.subs" = sim_cumu[i]))), sd(residuals(lnscales_mod)), lower = min(lnscales), upper = max(lnscales))
    }
    
    occu <- predict(mod_occ, newdata = data.frame("pr" = sim_cumu), type = "response")
    
    yr <- ini_year #Start year
    fires <- matrix(ncol = 4, nrow = 0) ##Storage of fire record
    tab <- matrix(ncol = 4, nrow = 0)
    
    #TIME STEP LEVEL
    for (i in 1:simle){
      precip <- sim_cumu[i] #draw precipitation value (for storage later on)
      occurrence <- runif(1, 0, 1) < occu[i] #will fire occur?
      if(occurrence == F){
        tab_target <- 0 #total area burned
      }
      tab_achieved <- 0 #initial realised total burned area  
      
      if(occurrence == T){
        
        #Sample annual area burned
        tab_target <- sim_tab[i]/100 * sum(na.omit(cell_size[])) 
        
        #iii) Sample time step fire size distribitution
        fire_distr <- rlnorm(10000,loc[i], scale[i])
        fire_distr_sub <- fire_distr[fire_distr < max(area.subs) & fire_distr > min(area.subs)]
        }
      
      #FIRE LEVEL
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
