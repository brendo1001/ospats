#' opstats algorithm WITH SIMULATED ANNEALING option
#'
#'  run opstats algorithm
#'  @import Rcpp
#'  @import RcppArmadillo
#'  @author Brendan Malone
#'  @useDynLib ospats
#'  @importFrom Rcpp sourceCpp RcppArmadillo
NULL

#opstats algorithm WITH SIMULATED ANNEALING option

# 9 Juillet add external starting point ! NSA

##start of function
ospats<- function(data,
                  dRange,
                  nCycles,
                  dStart, # choose between kMeans (0) or CumrootSquare(1) or external (3)
                  ClusterStart = c() , # external for dStart == 3
                  dMaxrun,
                  dRSquare,
                  dStrata,
                  initialTemperature = 1,
                  coolingRate = 0.9999,                   debug=F,
                  verbose=T){

  ###################################################################################################################


  ###################EMBEDDED FUNCTION 2 ################################################################################################
  #distance matrix


  sourceCpp("src/fastPdist.cpp") #C++ function for calculating distances

  ###################E################################################################################################


  # Define structure for storing time series of criterion
  Eall<-NULL
  # set initial temperature
  Temp <- initialTemperature

  x<- data[,1]
  y<- data[,2]
  z<- data[,3]
  s2<-data[,4]
  s2s<- s2
  n = length(x)
  d20<- 0 # stroring distance in case of re-run

print("-----OSPATS INITIALISATION------")
# generate initial solution
if (dStart == 0)  #use kmeans of coordinates
{k1<- kmeans(x= cbind(x,y), centers= dStrata, nstart=10, iter.max = 500)
strat0<- as.matrix(k1$cluster)}
if (dStart == 1)  #use Cum-sqrt-f
{nclass<-  100
strat0<- cumsqfDel(z= z, nclass= nclass, ns= dStrata)[[1]]}
if (dStart == 3)  #use external solution
{ strat0<- as.matrix(ClusterStart) }


#                            INITIATION
# Vectorised calculation of distance matrix, without loop
if(d20==0){
  # calculate distance matrix
  xy = cbind(x,y)
  Lag<- fastPdist2(Ar = xy, Br = xy)


  # calculate prediction (z) difference (vectorized)
  z1 = array(z,dim=c(nrow(xy),1,1))
  z2 = array(z,dim=c(1,nrow(xy),1))
  xx<- matrix(rep(z1, length(x)), nrow=length(x),ncol= length(x))
  yy<-  matrix(rep(z2, each=length(x)), nrow=length(x),ncol= length(x))
  Z2<- (xx-yy)^2
  Z2<-  Z2/dRSquare        #compensation for leveling

  # calculate variance (s2)
  zs1 = array(s2,dim=c(nrow(xy),1,1))
  zs2 = array(s2,dim=c(1,nrow(xy),1))
  xx<- matrix(rep(zs1, length(x)), nrow=length(x),ncol= length(x))
  yy<- yy<- matrix(rep(zs2, each=length(x)), nrow=length(x),ncol= length(x))
  S2<- (xx+yy)


  # calculate matrix of maximum of covariance
  Cov_max <-  0.5 *S2
  Cov<-  Cov_max * exp(-3*Lag/dRange)
  d2<-  Z2 + S2 -2*Cov
  #d2 = S2 -2*Cov; % emulate equal z-pred's
  #d2 = Lag;    %use only Lag
  #d2 = Z2;     % use only Z (true or pred)
  d20<- d2} else {d2<- d20}

rm(Cov_max , Lag , yy , xx , zs1 , zs2 , Z2)
gc()

#format compact
#MinZ2 = min(min(Z2)), MeanZ2 = mean(mean(Z2)), MaxZ2 = max(max(Z2));
#MinS2 = min(min(S2)), MeanS2 = mean(mean(S2)), MaxS2 = max(max(S2));
#MinCov = min(min(Cov)), MeanCov = mean(mean(Cov)), MaxCov = max(max(Cov));

TOTd2<-  sum(colSums(d2))/2
ObjNoStr<-  sqrt(TOTd2)
ObarH1<-  ObjNoStr/n

# Calculate sums of d2's within strata (Sd2) and contributions of strata to Obj (cbObj).
Sd2<- matrix(0, nrow=1, ncol=dStrata)
for (strat in 1:dStrata){
  Sd2[1,strat]<-  0
  for (i in 1:(n-1)){
    if (strat0[i,1] == strat)
    {for (j in (i+1):n){
      if (strat0[j,1] == strat)
      {Sd2[1,strat]<-  Sd2[1,strat]+d2[i,j]}}}}}
gc()

Sd2_init = Sd2;

cbObj<- sqrt(Sd2)
Obj<- sum(cbObj)
ObarInit<- Obj/n

print("-----OSPATS START RUNS------")

#                          START RUNS
ObarBest<-  1000        #arbitrary large number to start with
StratBest<- matrix(0, nrow=n, ncol=1)
d1<- data

for (run in 1:dMaxrun){
  if (debug == TRUE){print("_________________________________________________________________________________________________________RUN::")}
  TotTransf<- 0
  Sd2<- Sd2_init;
  #                             ITERATIVE RE-ALLOCATION
  if (verbose==TRUE){print(paste("----------------------------------------Run:", run, sep=" "))}
  stratcy<-  strat0       # stores stratum number for each point
  change<-  1
  for (cycle in 1:nCycles){ # loop through cycles


    if (debug==TRUE){print("________________________________________________________________________________________________________CYCLES:::")}
    transfers<-  0

    u<- t(as.matrix(sample.int(n, size = n, replace = FALSE, prob = NULL))) # put grid points in random order

    for (gp in 1:ncol(u)){

      if (debug==TRUE){print("________________________________________________________GP:::")}

      if (debug==TRUE){print(paste("le gp est",gp)) }

      samp <- u[1,gp]
      Delta <-  0
      change <-  0         # indicator of transfer
      A <-  stratcy[samp]
      # remove t from A
      ij <-  which(stratcy == A)
      dA <- sum(d2[samp,ij])-d2[samp,samp]
      sumd2tinA <-  dA
      Sd2Amint <-  Sd2[1,A] - sumd2tinA

      if (Sd2Amint<0) break() # test
      cbObjA <- sqrt(Sd2Amint)

      for (stratnr in 1:dStrata){ #% idem between t and the points in different strata


        if (debug==TRUE) {print("___________________________________________________STRATNR:::")}

        Delta<-  0
        sumd2plus<-  0

        if (stratnr != A){


          if (debug==TRUE){print("stratnr != A")}
          # Add to B
          B <-  stratnr
          ij <-  which(stratcy == B)
          dB <- sum(d2[samp,ij])
          sumd2plus<- dB

          cbObjB<- sqrt(Sd2[1,B] + sumd2plus)
          Delta<- cbObjA + cbObjB - cbObj[1,A] - cbObj[1,B]


          if (Delta < -Obj*1e-10){
            # always accept improvements
            if (debug==TRUE){print("Delta < -Obj*1e-10")}
            if (debug==TRUE){print("Delta < -Obj*1e-10")}
            if (debug==TRUE) print(paste("le Delta est " , Delta))
            pr <- 1
          } else {
            pr <- exp(-abs(Delta) / Temp ) # use Boltzmann to judge if deteriorations should be accepted
            if (debug==TRUE){print("!!!!Delta > -Obj*1e-10 !!!!")}
            if (debug==TRUE) print(paste("le Delta est " , Delta))
            if (debug==TRUE) print(paste("la temp?rature est " , Temp))
            if (debug==TRUE) print(paste("le crit?re de Boltzman " , pr))
          } # end calculate prob

          pChange <- runif(n = 1) # draw uniform deviate

          if (pChange < pr) { # test  proposal

            change<- 1
            transfers <- transfers + 1
            stratcy[samp,1]<- B    # update stratification
            Sd2[1,A] <-  Sd2[1,A]- sumd2tinA
            Sd2[1,B] <- Sd2[1,B] + sumd2plus
            cbObj <-  sqrt(Sd2)
            Obj <-  sum(cbObj)
            Eall <- rbind( Eall , Obj )
            Delta <- 0
          }


        } # end if srtat!=A


        if (change == 1){break}

      } # End of For loop per stratum


    } # End of For loop per gp

    ###___________________________

    if(verbose==TRUE){print(paste(paste("Cycle Number=", cycle, sep=" "), paste("Transfers=", transfers, sep=" "), sep= "  :::---:::  "))}
    if(verbose==TRUE){print(paste('Temperature = ', Temp, sep=" "))}


    TotTransf<-  TotTransf + transfers
    # lower temperature
    Temp <- coolingRate * Temp

    if (transfers == 0) {break}
  } # End of cycles---------------------------------------------


  if(verbose==TRUE){print("-----END TRANSFER------")


    print(paste('Number of cycles= ', cycle, sep=" "))
    print(paste('Total number of transfers= ', TotTransf, sep=" "))}
  Obj<-  sum(cbObj)
  ObarFinal<- Obj/n
  print(paste('Objective function= ', ObarFinal, sep=" "))

  # Update if last ObarFinal is smaller than the smallest previous ObarFinal
  if (ObarFinal < ObarBest){
    ObarBest<- ObarFinal
    StratBest[,1]<- stratcy}}
d1<- cbind(data,StratBest)
names(d1)[ncol(d1)]<- "ospats_strata"
retval<- list(ObarBest, d1 , Eall)
return(retval)
}



