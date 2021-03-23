#' @title Modularity
#'
#' @description Compute the modularity and number communities of graph g,
#'    w.r.t. its community structure, obtained through the chosen clustering
#'    method.
#'
#' @param g a graph. If the graph is not connected, the largest connected
#'    component is selected.
#' @param method default is \code{\link[igraph]{cluster_louvain}}, but others
#'    are possible: e.g. \code{\link[igraph]{cluster_fast_greedy}},
#'    \code{\link[igraph]{cluster_infomap}}.
#' @param ... additional parameters passed to the chosen igraph clustering
#'    method. The default \code{\link[igraph]{cluster_louvain}} has a
#'    \code{weights} parameter, ehich can be a positive weight vector, or the
#'    \code{E(g)$weight} attribute, or \code{NA} if the graph has a 'weight'
#'    edge attribute, but you want to ignore it.
#'    Larger edge weights correspond to stronger connections.
#'    Check out parameters for others methods.
#' @seealso [igraph::cluster_infomap()]
#' @return list(Q, N) where:
#' \itemize{
#'   \item{Q}{the modularity}
#'   \item{N}{the number of communities}
#' }
#' @keywords modularity; segregation; networks
#' @references Newman
#' @export
compute_modularity <- function(g, method = igraph::cluster_louvain, ...) {
    LCC_g <- LCC_subgraph(g)$LCC
    communities <- method(LCC_g, ...)
    res <- list('Q'= max(communities$modularity),
                'N'= max(communities$membership))
    return(res)
}

#' Isolate the Largest Connected Component
#'
#' @description
#' Compute modularity and number communities
#'
#' @param g a graph.
#' @return list(LCC, leftovers) where:
#'
#' * LCC is the largest connected components
#' * leftovers the nodes (labels or indeces depending on the network) that do
#' not belong to the LCC.
LCC_subgraph <- function(g) {
    clu <- igraph::components(g)
    lcc <- which(clu$membership == which.max(clu$csize))
    leftovers <- names(which(clu$membership != which.max(clu$csize)))
    g_LCC <- igraph::induced_subgraph(g, lcc)
    return(list(
      "LCC" = g_LCC,
      "leftovers" = leftovers
    ))
}
