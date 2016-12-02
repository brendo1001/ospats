#For clustering the predictions based on quantiles of the data

##Inputs:
#z:
#nclass:
#ns:

cumsqfDel<- function(z, nclass, ns){ 
  # calculate the distribution based on nclass
  z2<- cbind(seq(1,length(z),by=1), z)
  d<- z2[order(z),2]
  io<- z2[order(z),1]
  nd<- length(d)
  dmin<- min(d)
  dmax<- max(d)
  step<- (dmax-dmin)/nclass
  h<- t(as.matrix(seq(dmin,dmax,by=step)))
  
  fx<- matrix(NA, nrow=nclass, ncol=2)
  fy<- matrix(NA, nrow=nclass, ncol=1)
  for (i in 1:nclass){
    ll<- h[1,i]
    ul<- h[1,i+1]
    ij=which(d>=ll & d<ul)
    if(i==nclass) {ij<- which(d>=ll & d<=ul)}
    nij<- length(ij)
    fx[i,1]<- ll
    fx[i,2]<- ul
    fy[i,]<- nij}
  
  sqf<- sqrt(fy)
  cumsqf<- cumsum(sqf)
  
  #set out boundaries based on ns (number of strata)
  cmax<- max(cumsqf)
  cmin<- min(cumsqf)
  cstep=(cmax-cmin) /ns
  cbound<- t(as.matrix(seq(cmin,cmax,by=cstep)))
  
  # find out corresponding class
  db<- matrix(NA, nrow=1, ncol= length(2:ns))
  for (i in 2:ns){
    cm<- abs(cumsqf-cbound[1,i])
    db[1,i-1]<- which(cm==min(cm))}
  
  idat<- matrix(0, nrow= nd, ncol= 1)
  xm<- matrix(0, nrow= ns, ncol=5)
  
  for (i in 1:ns){
    if(i==1) {ll<- 1} else {ll<- db[1,i-1]}
    if(i==ns) {ul<- nclass} else {ul<- db[1,i]}
    
    fll<- fx[ll,1]
    ful<- fx[ul,1]
    if(i==ns){ful<- fx[ul,2]}
    
    ij<- which(d>=fll & d<ful)
    if(i==ns) {ij<- which(d>=fll & d<=ful)}
    idat[ij,1]<- i
    dij<- d[ij]
    nij<- length(dij)
    
    xm[i,1]<- nij
    xm[i,2]<-mean(dij)
    xm[i,3]<-var(dij)
    xm[i,4]<-min(dij)
    xm[i,5]<-max(dij)}
  
  strat<- matrix(1,nrow=nd, ncol=1)
  strat[io,1]<- idat
  retval<- list(strat, xm)
  return(retval)}
###################################################################################################################