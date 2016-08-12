##Stratifcation Allocation
#Initial muddles
library(raster)
library("rasterVis", lib.loc="~/R/win-library/3.2")

setwd("I:/rcodes/opspats/allocation")
#test data set
dat<- read.table("myFile.txt", header = F,sep=",")
names(dat)<- c("x", "y", "pred", "s2")

##ospats parameters
ra=4000;
ncy=400;
start=1;
maxrun=1;
R2= 0.42;
H=5;

library(Ospats)
#Run Function
tester_1<- ospatsF(data = dat, # input data (data frame). 4 columns [ X, Y, Pred, Var]
                   dRange = ra, # spatial structure of prediction variance
                   nCycles = ncy, # Number of allocation iterations 
                   dStart = start, # choose between kMeans (0) or CumrootSquare(1) or external (3)
                   ClusterStart = c() , # external for dStart == 3 (Saby Input)
                   dMaxrun = maxrun, # Number of runs the algorithm will go through to find optimal allocation
                   dRSquare = R2, # Used for compensation of leveling
                   dStrata = H, # Number of strata
                   initialTemperature = 1, #simulated annealing parameter
                   coolingRate = 0.95,  # simulated annealing parameter 
                   debug=F, # Useful during development for troubleshooting issues
                   verbose=T) # Prints messages during the running of function
str(tester_1)


##Plot outputs
#plot ospats stratification
r1<- rasterFromXYZ(tester_1[[2]][,c(1,2,5)])
plot(r1, main="ospats stratification")

map.c <- as.factor(r1)
rat <- levels(map.c)[[1]]
rat[["strata"]] <- c("HVT_001", "HVT_002", "HVT_003", "HVT_004", "HVT_005")
levels(map.c) <- rat

area_colors <- c("#FF0000", "#38A800", "#73DFFF", "#FFEBAF", "#800000")
levelplot(map.c, col.regions = area_colors, xlab = "", ylab = "")






#####ALLOCATION
str(tester_1)

#file to do allocation on
datA<- read.table("myFile_all.txt", header = F,sep=",")
names(datA)<- c("x", "y", "pred", "s2")
Sd2<- tester_1$sumSq
ospatIN<- tester_1$stratFrame


#some output matirces
strat_t<- matrix(NA, nrow=nrow(datA), ncol=2)


#######
library(doParallel)
# Begin a parallel cluster and register it with foreach:
cpus = 4  # The number of nodes/cores to use in the cluster
cl <- makeCluster(spec = cpus)
# Register the cluster with foreach:
registerDoParallel(cl)
#######


#parra processing
#oper1 <- foreach(i=1:500) %dopar% { 
for (i in 1:nrow(datA)){
  Sd2_1<- Sd2
  StratBest<- as.matrix(ospatIN$ospats_strata)
  n<- nrow(ospatIN)
  subsA<- datA[i,]
  #subsA
  #Distance  of new point to existing data
  Zt2<- (subsA$pred-ospatIN$pred)^2 #prediciton
  St2<-  subsA$s2 + ospatIN$s2 #variance
  Lagt<- sqrt(((subsA$x - ospatIN$x)^2) + ((subsA$y - ospatIN$y)^2))
  Zt2<- Zt2/R2 
  
  # Calculate matrix of maximum of covariance
  Cov_max<- 0.5*St2
  Cov<- Cov_max*exp((-3*Lagt)/ra)
  dt2<- Zt2 + St2 -2*Cov
  
  #cleanup
  rm(Zt2, St2, Lagt, Cov)
  
  # Insert the new point to d2
  # Add a new row and column
  t1<- cbind(tester_1$covMat, 0) #
  t2<- rbind(t1, 0)
  rm(t1)
  # Insert the dt2 data
  t2[n+1, 1:n]<- dt2
  t2[1:n, 1+n]<- dt2
  
  A<- 1 # try with strata 1
  n<- n + 1
  StratBest<- rbind(StratBest,1)
  
  # Calculate sums of d2's within strata (Sd2) and contributions of strata to Obj (cbObj).
  ij<- which(StratBest==A)
  dA<- sum(t2[n,ij])
  Sd2_1[1,1]<- Sd2_1[1,1] + dA
  cbObj<-  sqrt(Sd2_1)
  
  for (strat in 2:H){
    # remove t from A
    A<- strat-1
    t<- n
    ij<- which(StratBest==A)
    dA<- sum(t2[t,ij]) - t2[t,t]
    sumd2tinA <- dA
    Sd2Amint = Sd2_1[1,A]- sumd2tinA
    cbObjA = sqrt(Sd2Amint)
    
    #Add to B
    B<-  strat
    ij<- which(StratBest==B)
    dB<- sum(t2[t,ij])
    sumd2plus<-  dB
    
    cbObjB = sqrt(Sd2_1[1,B] + sumd2plus)
    Delta = cbObjA + cbObjB - cbObj[1,A] - cbObj[1,B]
    #print(Delta)
    
    if(Delta < 0){ # Onbj function decrease? then  transfer
      StratBest[t,1] = B #update stratification
      Sd2_1[1,A] = Sd2_1[1,A] -sumd2tinA
      Sd2_1[1,B] = Sd2_1[1,B] + sumd2plus
      cbObj = sqrt(Sd2_1)
      Obj = sum(cbObj)
      Delta = 0}}
  strat_t[i,1]<- StratBest[t,1]
  Obj<- sum(cbObj)
  Obar<- Obj/n
  strat_t[i,2]<- Obar
  print(c(i,StratBest[t,1]))}
  
  






