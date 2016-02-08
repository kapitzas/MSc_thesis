require(raster)
require(rgdal)

path.fire <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire"
path.tsf <- "/Users/Simon/Studium/MSC/Masterarbeit/data/fire/tsf rasterized 1970-2015"
files <- list.files(path.tsf, pattern = "*.asc$")
tsf.ts <- stack()
for (i in 1:length(files)){
  r <- raster(paste0(path.tsf, "/",files[i]))
  tsf.ts <- stack(tsf.ts,r)
}

fire.scar.ts[[1]]
tsf <- tsf.ts[[28]][]
#probability maps for each year
#P * B * tsf * MI^(e+2)
B <- mfi.per.lt
MI <- mfi.per.lt
fire.probs <- B * tsf *  MI^-(exp(1)+2)

# #calibration of B
# while(length(which(fire.probs > x)) < av.fire.size){
#   B = B + 10
#   fire.probs <- B * tsf *  MI^-(exp(1)+2)
#   x <- runif(length(fire.probs),0,1)
#   print(length(which(fire.probs > x)))
# }
# 
# Bcal <- B/mfi.lt[]
# #average number of burned cells per year
# av.fire.size <- length(which(fire.scar.ts[] == 1))/length(files)

