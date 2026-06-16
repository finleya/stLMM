\name{stunitco}
\alias{stunitco}
\docType{data}
\title{FIA state, unit, and county polygons}
\usage{
stunitco
}
\description{
An \code{sf} polygon data frame with FIA state, unit, and county codes and
names.
}
\format{
An object of class \code{sf} inheriting from \code{data.frame}, with 3,108 rows
and 8 columns:
\describe{
\item{\code{STATECD}}{FIA state code.}
\item{\code{STATENM}}{State name.}
\item{\code{UNITCD}}{FIA unit code.}
\item{\code{UNITNM}}{FIA unit name.}
\item{\code{COUNTYCD}}{County code.}
\item{\code{COUNTYNM}}{County name.}
\item{\code{COUNTYFIPS}}{County FIPS code.}
\item{\code{geom}}{Multipolygon geometry.}
}
The package copy is stored in NAD83 / Conus Albers, EPSG:5070.
}
\source{
Downloaded from the United States Census Bureau on 2019 November 3, originally
in Esri Shapefile format from
\url{https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html}.
The original projection was Geographic Coordinate System North American 1983,
EPSG:4269. The package copy is projected to EPSG:5070 and omits a small number
of very remote island polygons from the original source layer.
}
\examples{
data(stunitco)
nrow(stunitco)
names(stunitco)

\dontrun{
## This errors because the polygon layer has disconnected island components.
try(car_graph(stunitco, id = "COUNTYFIPS"))

## Add explicit nearest-neighbor bridge edges and inspect them.
g <- car_graph(stunitco, id = "COUNTYFIPS", island = "nearest")
head(g$island_added_edges)
}
}
\seealso{
\code{\link{car_graph}}
}
\keyword{datasets}
