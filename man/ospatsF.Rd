\name{ospatsF}
\docType{methods}
\alias{ospatsF}
\title{Optimizing Stratification and Allocation for Design-Based Estimation of Spatial Means Using Predictions with Error.}

\description{Need to populate this}

\usage{ospatsF(data= data, dRange= dRange, nCycles= nCycles, dStart= 1, ClusterStart= NULL, dMaxrun=1, dRSquare=0.5, dStrata=4, debug=FALSE, verbose=TRUE)}

\arguments{
  \item{covariates}{\code{RasterStack}; A stack of grids or rasters that are surrogates for environmental variables to both guide the C5 model fitting and provide the necessary environmental data for the spatial mapping of predicted soil classes.} 
  \item{polygons}{\code{SpatialPolygonsDataFrame}; See \code{data(dsT_polygons)} for an example of requirements. Each polygon needs to have the following necessary attribute data: A numeric unique map unit identifier, a character map unit name or code. And coupled columns of map unit classes and their respective proportions.}     
  \item{composition}{\code{data.frame}; See \code{data(dsT_composition)} for an example of requirements. First column indicates a numeric unique map unit identifier. Second column indicates character map unit name or code. Third column indicates map unit compositions. And fourth column indicates the percentage proportion of map unit compositions. All map unit codes and their subsequent compositions are row stacked together.}
  \item{n}{numeric; number of samples to take from each soil mapping polygon for C5 model fitting.}
  \item{reals}{numeric; number of C5 modeling fitting and mapping realisations to implement. }
  \item{cpus}{numeric; number of compute nodes to use. Default is 1.}
}

\value{This function returns a number of outputs, all of which are saved to the working directory in a folder called \code{dsmartOuts}. 
\describe{
   \item{\code{dsmartModels.Rdata}}{This \code{R} object contains in a list, the C5 model object from each realisation of the algorithm.}
   \item{\code{models}}{folder; This folder contains the summary output of each C5 model fitted from each realisation of the algorithm. }
   \item{\code{rasters}}{folder; This folder contains the rasters or soil class map created from each C5 model realisation and are labelled sequentially. This folder also contains the file \code{classLookupTable.txt} which is a lookup table to link unique numeric labels to soil class character names.}  
}}

\note{No notes at this stage}

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
