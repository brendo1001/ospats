\name{ospatsF}
\docType{methods}
\title{Optimizing Stratification and Allocation for Design-Based Estimation of Spatial Means Using Predictions with Error.}

\description{Given some information (map) of a target variable with associated error variance, Ospats is able to produce stratifications that represent transitions between the ‘knowing nothing about the target variable’ (high mapping uncertainty) and ‘knowing a lot about the target variable’ (low mapping uncertainty) situations. The theoretical and mathematical explanation of Ospats is described in detail in de Gruijter et al (2015). Stratification is optimised by minimising the expected sampling variance. The objective function O is minimized by an iterative re-allocation algorithm similar to those for k-means, except that here mutual distances are taken between grid points instead of distances between grid points and centroids.}

\usage{ospatsF(data= NULL, dRange= NULL, nCycles= NULL, dStart= NULL, ClusterStart= NULL, dMaxrun=NULL, dRSquare=NULL, dStrata=NULL, initialTemperature = NULL, coolingRate = NULL, debug=FALSE, verbose=TRUE)}

\arguments{
  \item{data}{\code{data.frame}; Needs to be 4 columns. Column 1 and 2 are the coorodinates. Column 3 are the predictions. Column 4 are the prediction variances.}
  \item{dRange}{\code{numeric}; A distance value of the autocorrelation of prediction errors.}     
  \item{nCycles}{\code{numeric}; The number of allowable allocation cycles in order to minimise the objective function.}
  \item{dStart}{\code{numeric}; Starting stratifcation. 0 for strata based on clustering of coordinates, or 1 for clustering based on quantiles of mapped variable. }
  \item{ClusterStart}{\code{numeric}; Not sure how to specify paramter yet. It is included to allow some to be able to specify their own starting stratifcation instead of using ones specified for \code{dStart}.}
   \item{dMaxrun}{\code{numeric}; Number of iterations one want to perfrom the allocation cycles. The default is set to 1.}
  \item{dRSquare}{\code{numeric}; Goodness of fit of the model that created the input prediction map. Value is used for preferencing allocations.}
  \item{dStrata}{\code{numeric}; number of strata one wants to create.}
  \item{initialTemperature}{\code{numeric}; Annealing parameter.}
  \item{coolingRate}{\code{numeric}; Annealing parameter. Determine the rate of cooling.}
  \item{debug}{\code{logical}; function troubleshooting option.}
  \item{verbose}{\code{logical}; generate additional outputs while function is running.}
}

\value{This function returns a list containing a number of elements. 
\describe{
   \item{objective}{The is the value of the objective function. This is similar to a sampling variance of which ospats tries to minimise.}
   \item{stratFrame}{This is a dataframe of 5 columns: Every row is an instance of coordinates (x and Y), prediction, error variance and alocated strata number.}
   \item{annealOut}{SOme diognostic information regarding the annealling schedule.}  
   \item{covMat}{Covariance matrix of input data}  
   \item{sumSq}{Vector os strata sums of squares}  
   \item{ospatsParams}{Vector of the ospats auguments.}  
}}

\note{A near future update of this function will be the ability to use rasters as the input data. In theory the solution (stratification) should be the same irrespective of the starting stratification. It was noted in a very early development of this function that this would not always eventaute for every situation. The process of converging to a final solution was very much railroaded. The implementation of this function however is a little different in that simulated annealing allows alternate realisations to be accepted along the way to minimisaion of the objective function. This new features avoids the railroading issue. maybe in future function development, it may be optional to select the annealing routine or not. }

\author{Brendan Malone and Nicolas Saby}
\references{
\itemize{
\item de Gruijter, J.J., Minasny, B., McBratney, A.B., (2015) \href{http://dx.doi.org/10.1093/jssam/smu024}{Optimizing Stratification and Allocation for Design-Based Estimation of Spatial Means Using Predictions with Error}. Journal of Survey Statistics and Methodology 3(1), 19-42. 
\item de Gruijter, J.J., McBratney, A.B., Minasny, B., Wheeler, I., Malone, B.P., Stockmann, U., (2016) \href{http://dx.doi.org/10.1016/j.geoderma.2015.11.010}{Farm-scale soil carbon auditing}. Geoderma 265, 120-130.
}}



\examples{

####NOT RUN
#library(raster)
#library(rasterVis)

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
###

}
\keyword{methods}
