\name{ospatsF}
\docType{methods}
\alias{ospatsF}
\title{Optimizing Stratification and Allocation for Design-Based Estimation of Spatial Means Using Predictions with Error.}

\description{Need to populate this}

\usage{ospatsF(data= data, dRange= dRange, nCycles= nCycles, dStart= 1, ClusterStart= NULL, dMaxrun=1, dRSquare=0.5, dStrata=4, debug=FALSE, verbose=TRUE)}

\arguments{
  \item{data}{\code{data.frame}; Needs to be 4 columns. Column 1 and 2 are the coorodinates. Column 3 are the predictions. Column 4 are the prediction variances.}
  \item{dRange}{\code{numeric}; A distance value of the autocorrelation of prediction errors.}     
  \item{nCycles}{\code{numeric}; The number of allowable allocation cycles in order to minimise the objective function.}
  \item{dStart}{\code{numeric}; Starting stratifcation. 0 for strata based on clustering of coordinates, or 1 for clustering based on quantiles of mapped variable. }
  \item{ClusterStart}{\code{!!!}; Not sure how to specify paramter yet. It is included to allow some to be able to specify their own starting stratifcation instead of using ones specified for \code{dStart}.}
   \item{dMaxrun}{\code{numeric}; Number of iterations one want to perfrom the allocation cycles. The default is set to 1.}
  \item{dRSquare}{\code{numeric}; Goodness of fit of the model that created the input prediction map. Value is used for preferencing allocations.}
  \item{dStrata}{\code{numeric}; number of strata one wants to create.}
  \item{debug}{\code{logical}; function troubleshooting option.}
  \item{verbose}{\code{logical}; generate additional outputs while function is running.}
}

\value{NEED TO FILL THIS IN BETTER}. 
\describe{
   \item{\code{dsmartModels.Rdata}}{This \code{R} object contains in a list, the C5 model object from each realisation of the algorithm.}
   \item{\code{models}}{folder; This folder contains the summary output of each C5 model fitted from each realisation of the algorithm. }
   \item{\code{rasters}}{folder; This folder contains the rasters or soil class map created from each C5 model realisation and are labelled sequentially. This folder also contains the file \code{classLookupTable.txt} which is a lookup table to link unique numeric labels to soil class character names.}  
}}

\note{A near future update of this function will be the ability to use rasters as the input data. In theory the solution (stratification) should be the same irrespective of the starting stratification. It was noted in a very early development of this function that this would not always eventaute for every situation. The process of converging to a final solution was very much railroaded. The implementation of this function however is a little different in that simulated annealing allows alternate realisations to be accepted along the way to minimisaion of the objective function. This new features avoids the railroading issue. maybe in future function development, it may be optional to select the annealing routine or not. }

\author{Brendan Malone and Nicolas Saby}
\references{
\itemize{
\item de Gruijter, J.J., Minasny, B., McBratney, A.B., (2015) \href{http://dx.doi.org/10.1093/jssam/smu024}{Optimizing Stratification and Allocation for Design-Based Estimation of Spatial Means Using Predictions with Error}. Journal of Survey Statistics and Methodology 3(1), 19-42. 
\item de Gruijter, J.J., McBratney, A.B., Minasny, B., Wheeler, I., Malone, B.P., Stockmann, U., (2016) \href{http://dx.doi.org/10.1016/j.geoderma.2015.11.010}{Farm-scale soil carbon auditing}. Geoderma 265, 120-130.
}}



\examples{
#library(dsmart)
#library(raster)

#load the relevent data

#Covariates
#data(dsT_covariates)

#Polygons
#data(dsT_polygons)

#Map unit compositions
#data(dsT_composition)

#Run dsmart (with 15 samples per polygon, 5 C5 model realisations, using 1 compute node)
#dsmart(covariates = dsT_covariates, polygons = dsT_polygons, composition = dsT_composition, n=15, reals = 5, cpus=1)

}
\keyword{methods}
