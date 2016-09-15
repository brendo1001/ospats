
library(Ospats)
library(raster)
library(rasterVis)

#Get some data
data(nowley_Cstock)


##ospats parameters
tester_1<- ospatsF(data = nowley_Cstock, # input data (data frame). 4 columns [ X, Y, Pred, Var]
                   dRange = 4000, # spatial structure of prediction variance
                   nCycles = 400, # Number of allocation iterations 
                   dStart = 1, # choose between kMeans (0) or CumrootSquare(1) or external (3)
                   ClusterStart = c() , # external for dStart == 3 (Saby Input)
                   dMaxrun = 1, # Number of runs the algorithm will go through to find optimal allocation
                   dRSquare = 0.42, # Used for compensation of leveling
                   dStrata = 5, # Number of strata
                   initialTemperature = 1, #simulated annealing parameter
                   coolingRate = 0.95,  # simulated annealing parameter 
                   debug=F, # Useful during development for troubleshooting issues
                   verbose=T) # Prints messages during the running of function

str(tester_1)


#plot ospats stratification
r1<- rasterFromXYZ(tester_1[[2]][,c(1,2,5)])
map.c <- as.factor(r1)
rat <- levels(map.c)[[1]]
rat[["strata"]] <- c("HVT_001", "HVT_002", "HVT_003", "HVT_004", "HVT_005")
levels(map.c) <- rat
area_colors <- c("#FF0000", "#38A800", "#73DFFF", "#FFEBAF", "#800000")
levelplot(map.c, col.regions = area_colors, xlab = "", ylab = "", main= "Ospats Strata")
###
