
setwd("/home//brendo/myWork/ospats//sandbox/rPackage/ospats//muddles")
#Run Ospats
library(raster);library(gstat);library(Rcpp);library(RcppArmadillo)



#inputs

data<- read.table("NowleyCStock2014_75m.txt", sep=",",header=TRUE) #soil data
data<- data[,c(1,2,4,5)]
names(data)<- c("x", "y", "pred", "s2")


##rasters
p.map<- rasterFromXYZ(data[,c(1:3)])
plot(p.map)

u.map<- rasterFromXYZ(data[,c(1,2,4)])
plot(u.map)



#variogram
vgm1 <- variogram(s2~1, ~x+y, data, width= 100)
mod1<-vgm(psill= var(data$s2), "Exp", range= 1000, nugget = 0) 
model_globe<-fit.variogram(vgm1, mod1)
model_globe
plot(vgm1,model=model_globe)




ra<- 575 #spatial 'range' for error (from variogram)
ncy<- 400 
start<- 1
maxrun<- 1
R2<-  0.4 #R2 of regression
H<- 5 # no strata



#Run Function
tester<- ospatsF(data= data, 
                dRange= ra, 
                nCycles= ncy, 
                dStart= 0, 
                ClusterStart = c() , # external for dStart == 3
                dMaxrun=1, 
                dRSquare=0.4, 
                dStrata=4, 
                initialTemperature = 1,
                coolingRate = 0.95, 
                debug=FALSE, 
                verbose=TRUE)



##Plot outputs
#plot ospats stratification
r1<- rasterFromXYZ(tester[[2]][,c(1,2,5)])
plot(r1, main="ospats stratification")


#Cum-sqrt-f
nclass<-  100
strat0<- cumsqfDel(z= data[,3], nclass= nclass, ns= 5)[[1]]  
r2<- rasterFromXYZ(cbind(tester[[2]][,c(1,2)],strat0))
plot(r2, main="Cum-sqrt-f stratification")


#use kmeans of coordinates
k1<- kmeans(x= data[,1:2], centers= 5, nstart=10, iter.max = 500)
strat0<- as.matrix(k1$cluster)
r3<- rasterFromXYZ(cbind(tester[[2]][,c(1,2)],strat0))
plot(r3, main="kmeans coordinate stratification")



