# Set up of data to be used for package

setwd("Z:/rdev/ospats/muddles/")

dat<- read.table("NowleyCStock2014_75m.txt", header = T, sep = ",")
dat<- dat[,c(1,2,4,5)]
names(dat)<- c("X", "Y", "Pred", "S2")

#data file for package
nowley_Cstock <- dat
save(nowley_Cstock, file= "Z:/rdev/ospats/data/nowley_Cstock.rda")
