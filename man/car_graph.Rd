\name{car_graph}
\alias{car_graph}
\alias{plot.stLMM_car_graph}
\alias{plot.stLMM_graph}
\title{Build a CAR adjacency graph}
\usage{
car_graph(graph, id = NULL, queen = TRUE, island = c("error", "nearest"),
  island_k = 1L)

\method{plot}{stLMM_car_graph}(x, show_ids = FALSE, show_nodes = TRUE,
  show_bridges = TRUE, color_by_degree = FALSE, ...)

\method{plot}{stLMM_graph}(x, show_ids = FALSE, show_nodes = TRUE,
  color_by_degree = FALSE, ...)
}
\arguments{
\item{graph}{An \code{sf} polygon object or a square symmetric adjacency
matrix.}
\item{x}{A graph object returned by \code{car_graph()} or an internal
\code{stLMM_graph} object stored in a fitted model backend.}
\item{id}{For \code{sf} inputs, the polygon ID column used to match
\code{car(area_id, graph = ...)}, \code{dagar(area_id, graph = ...)}, or
\code{dagar_time(area_id, time, graph = ...)}
values.}
\item{queen}{For \code{sf} inputs, whether to use queen-style boundary
touching. If \code{FALSE}, rook-style edge sharing is used.}
\item{island}{Policy for disconnected polygon components. The default
\code{"error"} reports disconnected or isolated areas. \code{"nearest"} adds
nearest-neighbor bridge edges between disconnected \code{sf} polygon
components.}
\item{island_k}{Positive integer number of nearest bridge edges to add from
the selected disconnected component at each connection step when
\code{island = "nearest"}.}
\item{show_ids}{Logical; if \code{TRUE}, draw graph IDs at node locations.}
\item{show_nodes}{Logical; if \code{TRUE}, draw representative node points.}
\item{show_bridges}{Logical; if \code{TRUE}, bridge edges added by
\code{island = "nearest"} are highlighted.}
\item{color_by_degree}{Logical; if \code{TRUE}, color polygons or nodes by
graph degree.}
\item{...}{Additional arguments passed to the internal graph plotting helper.}
}
\description{
Normalizes an areal adjacency structure for use in \code{car()},
\code{dagar()}, \code{car_time()}, and \code{dagar_time()} terms in
\code{\link{stLMM}} formulas.
}
\details{
The graph stores binary adjacency, not a row-standardized matrix. The CAR
implementation uses a degree-scaled proper CAR precision proportional to
\code{D - rho W}, where \code{D} is the diagonal degree matrix and \code{W} is
the symmetric binary adjacency matrix. This is equivalent to a CAR conditional
specification with row-standardized neighbor averaging in the conditional mean,
\eqn{E(w_i \mid w_{-i}) = rho d_i^{-1}\sum_{j \sim i} w_j}, with conditional
variance \eqn{sigma_sq / d_i}. The process variance \code{sigma_sq} is supplied
as the model term's variance parameter. Isolated areas currently error.
Disconnected components also error by default.

For matrix input, row and column names define the graph IDs. If names are not
present, integer IDs are used. Matrix input must be square, symmetric, and have
no isolated rows. Nonzero off-diagonal entries are treated as adjacency
indicators and are stored as binary edges. Nearest-neighbor island bridging is
not available for matrix input because no polygon geometry is available.

For \code{sf} input, \code{id} names the polygon identifier column. The helper
uses polygon boundary relationships to build a binary adjacency matrix. With
\code{queen = TRUE}, polygons that touch at a point or along an edge are
neighbors. With \code{queen = FALSE}, polygons must share an edge.

With \code{island = "nearest"}, disconnected components are connected by
artificial nearest-neighbor bridge edges based on polygon representative-point
distances. These edges are stored in \code{island_added_edges} so the graph
construction choice is auditable.

For teaching and auditability, a useful workflow is to first call
\code{car_graph()} with the default \code{island = "error"} and then rerun with
\code{island = "nearest"} only when disconnected island components are expected.

The \code{plot()} methods return \code{ggplot2} objects. For \code{sf}-based
graphs, polygons are drawn with graph edges overlaid. For matrix-based graphs,
a deterministic circular node layout is used. Internal DAGAR graph objects are
drawn with directed arrows from each parent node to its child node, making the
ordering-induced directed graph visible.
}
\value{
An object of class \code{stLMM_car_graph}. The object stores graph IDs, degree
information, and lower-triangle edge pairs used by \code{car()},
\code{dagar()}, \code{car_time()}, and \code{dagar_time()}. Island-bridging metadata are stored in
\code{island_policy}, \code{island_k}, \code{island_components_initial}, and
\code{island_added_edges}.
}
\examples{
adj <- matrix(
  c(
    0, 1, 0, 0,
    1, 0, 1, 0,
    0, 1, 0, 1,
    0, 0, 1, 0
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(paste0("a", 1:4), paste0("a", 1:4))
)
g <- car_graph(adj)

\dontrun{
data(stunitco)

## The bundled polygon layer has disconnected island components, so the
## conservative default reports the problem.
try(car_graph(stunitco, id = "COUNTYFIPS"))

## Make the graph-building rule explicit and inspect the artificial bridges.
g_stunitco <- car_graph(stunitco, id = "COUNTYFIPS", island = "nearest")
nrow(g_stunitco$island_added_edges)
head(g_stunitco$island_added_edges)

dat <- data.frame(
  y = rnorm(4),
  x = rnorm(4),
  area = paste0("a", 1:4)
)

fit <- stLMM(
  y ~ x + car(area, graph = g),
  data = dat,
  priors = list(
    resid = list(tau_sq = ig(2, 1)),
    car_1 = list(sigma_sq = ig(2, 1), rho = uniform(0.05, 0.95))
  ),
  n_samples = 1000
)
}
}
\seealso{
\code{\link{stunitco}}, \code{\link{stLMM}}, \code{\link{stLMM-terms}}
}
