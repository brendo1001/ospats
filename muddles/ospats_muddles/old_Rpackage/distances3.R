sx<- matrix(NA, nrow=length(s2), ncol=length(s2))

system.time(for (i in 1:length(s2)){
  sx[i,]<- s2[i]+ s2})
  
sourceCpp("fastVsum.cpp") #C++ function for calculating distances
sx2<- system.time(fastVsum(Ar = as.matrix(s2))
