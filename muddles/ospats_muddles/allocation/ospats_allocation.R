##Stratifcation Allocation
#Initial muddles
library(raster);library(Ospats);library(parallel);library(doParallel)
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
#plot(r1, main="ospats stratification")

map.c <- as.factor(r1)
rat <- levels(map.c)[[1]]
rat[["strata"]] <- c("HVT_001", "HVT_002", "HVT_003", "HVT_004", "HVT_005")
levels(map.c) <- rat

area_colors <- c("#FF0000", "#38A800", "#73DFFF", "#FFEBAF", "#800000")
levelplot(map.c, col.regions = area_colors, xlab = "", ylab = "")



#####ALLOCATION
#file to do allocation on
##small subset of the exisitng data but a much finer resolution
datA<- read.table("Nowley2014_sub.txt", header = T,sep=",")
names(datA)<- c("x", "y", "pred", "s2")
datA$s2<- datA$s2^2 *2

#rasterise
pred.r<- rasterFromXYZ(datA[,c(1,2,3)])
s2.r<- rasterFromXYZ(datA[,c(1,2,4)])
stacks<- stack(pred.r, s2.r)
stacks
#####################################################

temp1<- ospats_allocate(grids = stacks , ospat.obj= tester_1, cores = 4)
plot(temp1[[1]])






##OSPATS ALLOCATION FUNCTION
ospats_allocate<- function(grids = stacks , ospat.obj= tester_1, cores = 4){

#inputs from ospats
tester_1<- ospat.obj
stack.ins<- grids

Sd2<- tester_1$sumSq #strata sums of squares
ospatIN<- tester_1$stratFrame #statification data frame

#some output rasters
all.r<- raster(stack.ins[[1]])
all.r<- writeStart(all.r,filename='allocated_strata.tif', format='GTiff', overwrite=TRUE) #allocated strata
all.o<- raster(stack.ins[[1]])
all.o<- writeStart(all.o,filename='objectiveF.tif', format='GTiff', overwrite=TRUE) #objective function

##############################
# Begin a parallel cluster and register it with foreach:
cpus = cores # The number of nodes/cores to use in the cluster
cl <- makeCluster(spec = cpus)
# Register the cluster with foreach:
registerDoParallel(cl)

##For each row of the raster
pb <- txtProgressBar(min=0, max=dim(all.r)[1], style=3)
for(i in 1:dim(all.r)[1]){ #for each line of each input raster
#for(i in 1:25){ #for each line of each input raster
#oper1 <- foreach(i=1:25,.packages=c("raster", "rgdal")) %dopar% {
  cov.Frame<- getValues(stack.ins,i) #get raster values
  cov.cell<- cellFromRow(stack.ins, i)
  sub.frame<- cov.Frame[which(complete.cases(cov.Frame)),] #get the complete cases
  sub.index<- which(complete.cases(cov.Frame))
  sub.cell<- cov.cell[sub.index]
  sub.xy<- xyFromCell(stack.ins, sub.cell, spatial=FALSE)
  sub.frame<- as.data.frame(cbind(sub.xy,sub.frame))
  a.matrix <- matrix(NA, nrow = nrow(cov.Frame), ncol = 2)
  
  ##For each element of the row
  oper1<- foreach(qq=1:nrow(sub.frame),.combine = rbind) %dopar% {  #element by element in paralel
  #oper1<- for (qq in 1:nrow(sub.frame)){
    Sd2_1<- Sd2
    StratBest<- as.matrix(ospatIN$ospats_strata)
    n<- nrow(ospatIN)
    subsA<- sub.frame[qq,]
    subsB<-sub.index[qq]
    
    #subsA
    #Distance  of new point to existing data
    Zt2<- (subsA$pred-ospatIN$pred)^2 #prediciton
    St2<-  subsA$s2 + ospatIN$s2 #variance
    Lagt<- sqrt(((subsA$x - ospatIN$x)^2) + ((subsA$y - ospatIN$y)^2))
    Zt2<- Zt2/tester_1$ospatsParams[5] 
    
    # Calculate matrix of maximum of covariance
    Cov_max<- 0.5*St2
    Cov<- Cov_max*exp((-3*Lagt)/tester_1$ospatsParams[1])
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
    
    for (strat in 2:tester_1$ospatsParams[6]){
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
    
  Obj<- sum(cbObj)
  Obar<- Obj/n
  c(StratBest[t,1], Obar)}
  
  a.matrix[sub.index,]<- oper1
  all.r <- writeValues(all.r, a.matrix[, 1], i)
  all.o <- writeValues(all.o, a.matrix[, 2], i)
  setTxtProgressBar(pb, i)}
all.r<- writeStop(all.r)
all.o<- writeStop(all.o)
stopCluster(cl)     
message(paste("OSPATS outputs can be located at:",getwd(), sep=" "))
out.stack<- stack(all.r, all.o)
return(out.stack)
}
















