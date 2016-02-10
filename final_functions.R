distance <- function(from, to, layer){
  dist <- numeric()
  for (i in 1:length(from)){
    distxy <- xyFromCell(layer, to[i]) - xyFromCell(layer, from[i])
    dist[i] <- sqrt((distxy[1])^2 + (distxy[2])^2)  
  }
  return(dist)
}