#' @title Local Communication Efficeiency
#' @description The average inverse shortest path length is a measure known as
#'    the global efficiency (see Latora and Marchiori, 2001).
#'    We implement the local communication efficiency of node \eqn{i \in G} as
#'    the global communication efficency of the subgraph induced by \eqn{i},
#'    \code{\link{GCE}}.
#'
#' @param g a network.
#' @param directed default FALSE.
#' @param normalised logical, default TRUE.
#' @param weights edge weights, representing the appeal of each link.
#'    if weights is NULL (the default) and g has a weight edge
#'    attribute that will be passed to the igraph
#'    \code{\link[igraph]{distances}} function, as \eqn{w_{ij}^{-1}}.
#' @return local communication efficiency vector, for each node \eqn{i \in G}
#'    \eqn{LCE(G_i)} if \code{normalised = TRUE}, \eqn{E(G_i)} otherwise.
#' @keywords communication efficiency; integration; centrality;
#' @references Latora, V., & Marchiori, M. (2001). Efficient Behavior of Small-World Networks. \url{https://doi.org/10.1103/PhysRevLett.87.198701}
#' @examples
#' library(intsegration)
#' library(igraph)
#' karate <- make_graph("zachary")
#' LCE(karate)
#' @export
LCE <- function (g, directed = FALSE, normalised = TRUE, weights = NULL)
{
  if (is.numeric(weights)) {
    g <- igraph::set_edge_attr(g, name = "weight", value = weights)
  } else if (is.null(weights)) {
    # do nothing
  } else {
    warning("Check weight attribute! Uniform unit weights will be used.")
    igraph::E(g)$weight <- 1
  }
  # list of the neighbourhood g_i of each vertex i
  gi_list <- igraph::make_ego_graph(g, order = 1, nodes = igraph::V(g),
    mode = "all")
  sapply(gi_list, function(gi) {
    GCE(gi, directed = directed, weights = NULL, verbose = F)
  })
}
