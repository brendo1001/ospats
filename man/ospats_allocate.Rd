\name{ospats allocation}
\docType{methods}
\title{Allocation of ospats strata to further additional grid points.}

\description{These words are more-or-less verbatim from Jaap de Gruijter via personal communication. We assume that a grid of N points has been divided into H strata by Ospats \code{ospatsF}. We assume further that additional grid points (not included in the stratifcation process) are to be allocated to the existing Ospats strata. The same data are assumed to be available for the additional grid points as for the original points, i.e. predictions and prediction error variances. This function optimally allocates these additional data to the ospats strata. The principle of the proposed allocation method is to assign each additional grid point in turn to one of the strata while minimizing the same objective function O as was used to create the stratifcation. This means that, for the allocation of point xt, we have to fnd the stratum Sh of which the increase of its contribution to O by adding xt to Sh is minimal.}

\usage{ospats_allocate(grids = NULL, ospat.obj = NULL, cores = 1)}

\arguments{
  \item{grids}{\code{raster stack}; Map of predictions and map of associated error variance stacked together.}
  \item{ospat.obj}{\code{list}; The object that is created from running the \code{ospatsF} function.}     
  \item{cores}{\code{numeric}; The number compute cores the run the function on. Default is 1}
}

\value{This function returns raster stack of 2 maps. The first map is a map of the allocated strata. The second map is that is the objective function calculated at each pixel given the addtion of the new point.}

\note{At the moment the function works on input rasters. This may be a limitation if the addtional points are not in grid form and do not form a regular grid ie can not be converted to raster. }

\author{Brendan Malone using original script from Jaap de Gruijter and Budiman Minasny }
\references{
\itemize{
\item de Gruijter, J.J., Minasny, B., McBratney, A.B., (2015) \href{http://dx.doi.org/10.1093/jssam/smu024}{Optimizing Stratification and Allocation for Design-Based Estimation of Spatial Means Using Predictions with Error}. Journal of Survey Statistics and Methodology 3(1), 19-42. 
\item de Gruijter, J.J., McBratney, A.B., Minasny, B., Wheeler, I., Malone, B.P., Stockmann, U., (2016) \href{http://dx.doi.org/10.1016/j.geoderma.2015.11.010}{Farm-scale soil carbon auditing}. Geoderma 265, 120-130.
}}



\examples{
####NOT RUN
#library(raster)
#library(rasterVis)
#library(foreach)
#library(parallel)
#library(doParallel)

#Get some data
#data(nowley_Cstock)


##ospats parameters
#tester_1<- ospatsF(data = nowley_Cstock, # input data (data frame). 4 columns [ X, Y, Pred, Var]
#                   dRange = 4000, # spatial structure of prediction variance
#                   nCycles = 400, # Number of allocation iterations 
#                   dStart = 1, # choose between kMeans (0) or CumrootSquare(1) or external (3)
#                   ClusterStart = c() , # external for dStart == 3 (Saby Input)
#                   dMaxrun = 1, # Number of runs the algorithm will go through to find optimal allocation
#                   dRSquare = 0.42, # Used for compensation of leveling
#                   dStrata = 5, # Number of strata
#                   initialTemperature = 1, #simulated annealing parameter
#                   coolingRate = 0.95,  # simulated annealing parameter 
#                   debug=F, # Useful during development for troubleshooting issues
#                   verbose=T) # Prints messages during the running of function

#str(tester_1)


#plot ospats stratification
#r1<- rasterFromXYZ(tester_1[[2]][,c(1,2,5)])
#map.c <- as.factor(r1)
#rat <- levels(map.c)[[1]]
#rat[["strata"]] <- c("HVT_001", "HVT_002", "HVT_003", "HVT_004", "HVT_005")
#levels(map.c) <- rat
#area_colors <- c("#FF0000", "#38A800", "#73DFFF", "#FFEBAF", "#800000")
#levelplot(map.c, col.regions = area_colors, xlab = "", ylab = "", main= "Ospats Strata")


###Allocation
##Note that allocation would be used in situations in which additional grid points are avilable
## and need to be allocated to an exisiting ospats stratification. In the present example we are using
##the original data that was used for the stratification. We are just pretending they are additional data.

###
#rasterise data
#pred.r<- rasterFromXYZ(nowley_Cstock[,c(1,2,3)])
#s2.r<- rasterFromXYZ(nowley_Cstock[,c(1,2,4)])
#stacks<- stack(pred.r, s2.r)
#stacks
#####################################################

##run allocate function on 4 compute cores
#temp1<- ospats_allocate(grids = stacks , ospat.obj= tester_1, cores = 4)

}
\keyword{methods}
