#Run Ospats
library(raster);library(rasterVis);library(gstat)
#inputs
setwd("/home/malone/Dropbox/2016/rcode/")
data<- read.table("NowleyCStock2014_75m.txt", sep=",",header=TRUE) #soil data
xy_dat<- as.matrix(data[,c("x","y")])

###
naive_pdist <- function(A,B) {
  # A: matrix with obersvation vectors 
  #         (nrow = number of observations)
  #
  # B: matrix with another set of vectors
  #          (e.g. cluster centers)
  result = matrix(ncol=nrow(B), nrow=nrow(A))
  for (i in 1:nrow(A))
    for (j in 1:nrow(B))
      result[i,j] = sqrt(sum( (A[i,] - B[j,])^2 ))
    
    result
}

system.time(naive_pdist(A= xy_dat, B=xy_dat))

####
install.packages("pdist")
library(pdist)
system.time(as.matrix(pdist(X = xy_dat, Y = xy_dat)))



############
vectorized_pdist <- function(A,B) {
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))

m = nrow(A)
n = nrow(B)

tmp = matrix(rep(an, n), nrow=m) 
tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
sqrt( tmp - 2 * tcrossprod(A,B) )
}


system.time(vectorized_pdist(A= xy_dat, B=xy_dat))



library(Rcpp)
library(RcppArmadillo)
sourceCpp("fastPdist.cpp")
system.time(fastPdist2(Ar = xy_dat, Br = xy_dat))
distCalc<- fastPdist2(Ar = xy_dat, Br = xy_dat)
