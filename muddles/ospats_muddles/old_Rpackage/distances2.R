library(Rcpp);library(RcppArmadillo)
#inputs
setwd("/home/malone/Dropbox/2016/rcode/")
data<- read.table("NowleyCStock2014_75m.txt", sep=",",header=TRUE) #soil data
xy_dat<- as.matrix(data[,c("x","y")])


sourceCpp("fastPdist.cpp")
system.time(fastPdist2(Ar = xy_dat, Br = xy_dat))
distCalc<- fastPdist2(Ar = xy_dat, Br = xy_dat)
