####################################################
# Global Communication Efficiency: function for evaluation of GCE
#
# Version: 1.0
# Last update: May 2021
# Authors: Giulia Bertagnolli
#
# History:
#
# Dependencies:
#
# igraph, Rcpp
#
####################################################
#' @title Global Communication Efficiency
#' @description The average inverse shortest path length is a measure known as
#'    the global efficiency (see Latora and Marchiori, 2001).
#'    We implement the global communication efficiency (GCE) for weighted
#'    networks and propose a new normalisation method for the GCE.
#'    The global efficiency may be meaningfully computed on disconnected
#'    networks, as paths between disconnected nodes are defined to have
#'    infinite length and correspondingly zero efficiency.
#'    \strong{Assumption}: Edge weights are non-negative and represent the
#'    \emph{appeal} of links, hence we want shortest paths with heavy links.
#'
#'    \strong{Non-normalised GCE}:
#'    \deqn{E(g) = \frac{1}{N} \sum_{i \in V} \frac{\sum_{j \in V,i \neq j} d_{ij}}{N - 1},}
#'    where \eqn{N} is the number of vertices of graph \eqn{g} and \eqn{d_{ij}}
#'    is the shortest-path distance between \eqn{i} and \eqn{j}, computed by
#'    the Floyd-Warshall's algorithm with inverse edge weights as edge costs.
#'    See \url{https://en.wikipedia.org/wiki/Floyd-Warshall_algorithm}.
#'    \strong{Normalisation}: the previous efficiency should be divided by the
#'    efficency of an \emph{ideal} graph.
#'    \deqn{E(g_{ideal}) = \frac{1}{N} \sum_{i \in V} \frac{\sum_{j \in V,i \neq j} l_{ij}}{N - 1}.}
#'    The normalised communication efficiency is then given by
#'    \deqn{GCE(g)=\frac{E(g)}{E(g_{ideal})}}
#' @param g a network.
#' @param directed logical, if the directed network has to be considered.
#'    If FALSE (default) the network is taken as undirected.
#' @param normalised logical, default TRUE.
#' @param weights edge weights, representing the appeal of each link.
#'    if weights is NULL (the default) and g has a weight edge
#'    attribute that will be passed to the igraph
#'    \code{\link[igraph]{distances}} function, as \eqn{w_{ij}^{-1}}.
#'    If this is NA then no weights are used (even if the graph has a weight attribute).
#'    The normalisation of the GCE relises on the comparison with the efficiency of a
#'    \eqn{G_{ideal}}, characterised by the weighted adjacency matrix
#'    \eqn{G_{ideal} = \frac{W + \Phi}{2}}. See [2] for more details.
#' @param lengths Euclidean length of the edge, i.e., Euclidean distance between its
#'    two end-points, which are assumed to be non-negative reals. If this is NULL (the default)
#'    and \code{g} has an edge attribute \code{length} that will be used.
#'    If this is NA then no lengths are used (even if the graph edges have a length attribute).
#' @param verbose logical, default TRUE.
#' @return A list
#'   \itemize{
#'     \item non_normalised - communication efficiency \eqn{E(G)},
#'     \item normalised - normalised communication efficiency, \eqn{GCE(G)},
#'       or NULL if \code{normalised} is \code{FALSE}. If the network is topological (i.e. no flows)
#'       or if \code{weights = NA} then the normalised topological efficiency is equal to the
#'       non-normalised topological efficiency.
#'     \item normalised_lengths - normalised communication efficiency w.r.t. edge lengths,
#'       i.e. \eqn{GCE^{Eucl}(G)}: shortest-paths are computed using the flows (as in the usual
#'       normalied version), but now the path total cost and total strength are, respectively,
#'       the sum of Euclidean distances (lengths) and the sum of the reciprocals of the lengths of
#'       edges over the paths. NULL if there are non lengths.
#'       If \code{normalised = FALSE} but \code{lengths} are provided, then the \eqn{GCE^{Eucl}(G)}
#'       (that is \code{normalised_lengths}) will be computed anyway.
#'   }
#' @keywords communication efficiency; integration;
#' @references \describe{
#'     \item{[1]}{Latora, V. & Marchiori, M. (2001).
#'                Efficient Behavior of Small-World Networks.
#'                \url{https://doi.org/10.1103/PhysRevLett.87.198701}}
#'     \item{[2]}{Bertagnolli, G., Gallotti R. & De Domenico, M. (2020)
#'                Quantifying efficient information exchange in real flow
#'                networks. \url{arxiv}}
#'    }
#' @examples
#' library(intsegration)
#' library(igraph)
#' karate <- make_graph("zachary")
#' GCE(karate, directed = F)
#' @export
GCE <- function(g, directed = FALSE, normalised = TRUE, weights = NULL,
                lengths = NULL, verbose = TRUE) {
  # create return object
  res <- list(
    "non_normalised" = NULL,
    "normalised" = NULL,
    "normalised_lengths" = NULL
  )
  # aux variables
  topological <- FALSE
  use_lengths <- TRUE
  # directed or undirected graph
  if ((igraph::is.directed(g)) && (!directed)) g <- igraph::as.undirected(g)
  # set weight and inverse-of-weight attributes
  if (is.null(weights)) {
    # default case
    if ("weight" %in% igraph::edge_attr_names(g)) {
      # weight edge attribute
      g <- igraph::set_edge_attr(g, name = "weight_inv",
                                 value = (1. / igraph::E(g)$weight))
    } else {
      # no (numeric) weight edge attribute
      if (verbose) {
        cat("Unweighted graph, topological case.")
      }
      topological <- T
      g <- igraph::set_edge_attr(g, name = "weight", value = 1)
      g <- igraph::set_edge_attr(g, name = "weight_inv", value = 1)
    }
  } else if (length(weights) == 1 && is.na(weights)) {
    if (verbose) cat("Ignoring edge weights, topological case.")
    topological <- T
    g <- igraph::set_edge_attr(g, name = "weight", value = 1)
    g <- igraph::set_edge_attr(g, name = "weight_inv", value = 1)
  } else if (is.numeric(weights)) {
    # providing new weights
    if (verbose) {
      cat("Adding given weights and 1 / weights as edge attributes.\n")
    }
    igraph::E(g)$weight <- weights
    igraph::E(g)$weight_inv <- 1. / igraph::E(g)$weight
  } else {
    warning("Check edge weights, ignoring them. Topological efficiencies.")
    topological <- T
    g <- igraph::set_edge_attr(g, name = "weight", value = 1)
    g <- igraph::set_edge_attr(g, name = "weight_inv", value = 1)
  }
  # Euclidean distances
  if (is.null(lengths)) {
    # default case
    if ("length" %in% igraph::edge_attr_names(g)) {
      # using length edge attribute
      if (verbose) {
        cat("Using length edge attribute.")
      }
      # igraph::E(g)$length # length is an edge cost
      # 1 / igraph::E(g)$length # reciprocal of lengths are edge strengths
      g <- igraph::set_edge_attr(g, name = "length_inv",
                                 value = (1. / igraph::E(g)$length))
    } else {
      # no (numeric) length edge attribute
      if (verbose) {
        cat("No Euclidean lengths/distances provided.")
      }
      use_lengths <- FALSE
      # g <- igraph::set_edge_attr(g, name = "length", value = 1)
      # g <- igraph::set_edge_attr(g, name = "length_inv", value = 1)
    }
  } else if (length(lengths) == 1 && is.na(lengths)) {
    if (verbose) cat("Ignoring edge lengths.")
    use_lengths <- FALSE
    # g <- igraph::set_edge_attr(g, name = "length", value = 1)
    # g <- igraph::set_edge_attr(g, name = "length_inv", value = 1)
  } else if (is.numeric(lengths)) {
    # providing new edge lengths
    if (verbose) {
      cat("Adding given lengths as edge costs and 1 / lengths as edge strengths.\n")
    }
    igraph::E(g)$length <- lengths
    igraph::E(g)$length_inv <- 1. / igraph::E(g)$length
  } else {
    warning("Issues with edge lengths. They are currently ignored.")
    use_lengths <- FALSE
    # g <- igraph::set_edge_attr(g, name = "length", value = 1)
    # g <- igraph::set_edge_attr(g, name = "length_inv", value = 1)
  }
  # Check for multiedges and loops
  if (!igraph::is.simple(g)) {
    cat("graph is not simple. Aggrgating multiedges (sum of weights, min of lengths),
        removing self-loops.\n")
    min_na_rm <- function(x) {
      return(min(x, na.rm = TRUE))
    }
    g <- igraph::simplify(g, edge.attr.comb = list(
                               "weight" = "sum",
                               "length" = "min_na_rm")
                               )
    igraph::E(g)$weight_inv <- 1. / igraph::E(g)$weight
    igraph::E(g)$length_inv <- 1. / igraph::E(g)$length
  }
  # ---
  # Auxiliary objects
  # ---
  N <- igraph::gorder(g)
  A <- igraph::as_adjacency_matrix(g, attr = "weight")
  if (topological) {
    normalised <- FALSE
  }
  # ---
  # Start computing shortest-paths
  if (verbose) cat("Start computation of shortest-paths\n")
  x <- as.matrix(igraph::as_adjacency_matrix(g, attr = "weight_inv"))
  x[x == 0] <- .Machine$double.xmax         # almost as <- Inf
  if (!use_lengths) {
    # case without Euclidean lengths/distances
    g_all_shortest_paths <- rcpp_floyd_flow(x)
    rm(x)
  } else {
    # case with Euclidean lengths/distances
    l <- as.matrix(igraph::as_adjacency_matrix(g, attr = "length"))
    l[l == 0] <- .Machine$double.xmax         # almost as <- Inf
    g_all_shortest_paths <- rcpp_floyd_flow_length(x, l)
    rm(x)
    rm(l)
  }
  # split outputs
  # shortest-path distances
  D <- g_all_shortest_paths$D
  D[D == .Machine$double.xmax] <- NA
  diag(D) <- NA
  if (normalised) {
    # total flows along shortest-paths
    Phi <- g_all_shortest_paths$F
    diag(Phi) <- NA
  }
  if (use_lengths) {
    # total (sum of) Euclidean distance along (the same) shortest-paths
    D_eucl <- g_all_shortest_paths$D_eucl
    D_eucl[D_eucl == .Machine$double.xmax] <- NA
    diag(D_eucl) <- NA
    # sum of reciprocal Euclidean distances along (the same) shortest-paths
    Phi_eucl <- g_all_shortest_paths$F_eucl
    diag(Phi_eucl) <- NA
  }
  rm(g_all_shortest_paths)
  # ---
  # if (verbose) cat("Computing matrix Phi\n")
  # if (normalised) {
  #   if (topological) {
  #     Phi <- D
  #     Phi[is.na(Phi)] <- 0
  #     diag(Phi) <- NA
  #   } else {
  #     for (i in 1:(N - 1)) {
  #       if (!directed) {
  #         j0 <- i + 1
  #       } else {
  #         j0 <- 1
  #       }
  #       for (j in c(j0:N)) {
  #         Phi[i, j] <- sum(extract_paths_weights(g_all_shortest_paths, A, i, j))
  #       }
  #     }
  #     if (!directed) {
  #       Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(t(Phi))]
  #     }
  #   }
  # }
  if (verbose) cat("compute E(G)\n")
  E <- 1. / N / (N - 1) * sum(1. / D, na.rm = T)
  # save E(G) in res[["non_normalised"]]
  res$non_normalised <- E
  if (topological) {
    # if there are no flows, the topological efficiency is already normalised, hence
    # non_normalised and normalised efficencies are equal in the topological case.
    res$normalised <- E
  }
  if (normalised) {
    # also compute the normalised GCE
    W_ideal <- .5 * (A + Phi)
    E_ideal <- 1. / N / (N - 1) * sum(W_ideal, na.rm = T)
    # normalised GCE
    if (!identical(E_ideal, E)) {
      res$normalised <- E / E_ideal
    } else {
      if (verbose) {
        cat("G is not a true network, but a collaction of isolated nodes
            and pairs.\n")
      }
      res$normalised <- NA
    }
  }
  if (use_lengths) {
    E_eucl <- 1. / N / (N - 1) * sum(1. / D_eucl, na.rm = T)
    W_ideal <- .5 * (as.matrix(igraph::as_adjacency_matrix(g, attr = "length_inv")) + Phi_eucl)
    E_eucl_ideal <- 1. / N / (N - 1) * sum(W_ideal, na.rm = T)
    if (!identical(E_eucl_ideal, E_eucl)) {
      res$normalised_lengths <- E_eucl / E_eucl_ideal
    } else {
      if (verbose) {
        cat("G is not a true network, but a collaction of isolated nodes
            and pairs.\n")
      }
      res$normalised_lengths <- NA
    }
  }
  if (!igraph::is_connected(g)) {
    if (verbose) {
      warning("Network is not connected. Average among disconnected subgraphs.")
    }
    if (!topological) {
      # to-do:
      # check if weighted GCE is almost equal to topological
      # if TRUE return only topological with warning
      # condition on variance of weights?
    }
  }
  return(res)
}

# # auxiliary function
# extract_paths_weights <- function(obj, adj, start, end) {
#   sp <- e1071::extractPath(obj, start, end)
#   res <- numeric(0)
#   for (i in 2:length(sp)) {
#     res <- c(res, adj[sp[i - 1], sp[i]])
#   }
#   return(res)
# }
