

traj_path <- "/Users/Simon/Studium/MSC/Masterarbeit/data/Pop_traj_results/"

sc <- c("low", "neut", "high")
trajectories <- list()
traj_files <- list.files(traj_path)
for (j in 1:3){
  traj_all <- matrix(0, nrow = 76, ncol = 6)
    w <- grep(sc[j], traj_files)
  for (i in w){
    traj <- as.matrix(read.table(paste0(traj_path,traj_files[i]), skip = 15, nrows = 76, header = F))/length(w)
    traj_all <- traj_all + traj
  }
    traj_all[,1] <- traj_all[,1] + 2015
  trajectories[[j]] <- traj_all
}
