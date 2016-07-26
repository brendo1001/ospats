\name{nowley_Cstock}
\docType{data}
\alias{nowley_Cstock}
\title{Predictions of carbon stock with variance estimates}
\description{A \code{data.frame} with 4 columns which can be easily transformed into maps. Columns 1 and 2 contain spatial coordinates, while columns 3 and 4 contain predictions of carbon stocks and their prediction variance respectively. Each row more-or-less represents a pixel. This data is good for demonstration of using the ospats algorithm as described in de Gruijter et al. (2015).}

\usage{data(nowley_Cstock)}
\format{
\code{nowley_Cstock} is a 4 column, 4382 row \code{data.frame}.}
\details{This is an example data set in order to familiarize users to the working and inputs required for operationalising the \code{ospatsF} algorithm.}

\references{
\itemize{
\item de Gruijter, J.J., Minasny, B., McBratney, A.B., (2015) \href{http://dx.doi.org/10.1093/jssam/smu024}{Optimizing Stratification and Allocation for Design-Based Estimation of Spatial Means Using Predictions with Error}. Journal of Survey Statistics and Methodology 3(1), 19-42. 
\item de Gruijter, J.J., McBratney, A.B., Minasny, B., Wheeler, I., Malone, B.P., Stockmann, U., (2016) \href{http://dx.doi.org/10.1016/j.geoderma.2015.11.010}{Farm-scale soil carbon auditing}. Geoderma 265, 120-130.}}

\examples{
library(ospats)
data(nowley_Cstock)
str(nowley_Cstock)

###Not Run
#library(raster)
#prediction
#p.map<- rasterFromXYZ(nowley_Cstock[,c(1:3)])
#plot(p.map)

#variance
#u.map<- rasterFromXYZ(nowley_Cstock[,c(1,2,4)])
#plot(u.map)
###

}
\keyword{datasets}
