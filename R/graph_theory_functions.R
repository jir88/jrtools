#' Trim Weighted Graph to k Nearest Neighbors
#'
#' Given an igraph object with weighted edges, this function will create a
#' trimmed igraph with only the k nearest neighbors of each vertex. In this
#' context, nearest neighbors are those with the largest (or smallest)
#' edge weights.
#'
#' @param graph the igraph object to be trimmed
#' @param k the number of neighbors to keep
#' @param method "max" to keep k largest edge weights, "min" to keep the k
#'   smallest edge weights
#' @param edge_mode "out" to consider only outbound edges, "in" to consider only
#'   inbound edges, or "all" to consider both edge types. Ignored for
#'   un-directed graphs.
#' @return a trimmed copy of the input graph
#' @details This is a really simplistic kludge. Use at your own risk!
#' @export
knn_trim <- function(graph, k = 3, method = "max", edge_mode = "all") {
  # check class of graph
  if(!("igraph" %in% class(graph))) {
    stop(simpleError(cat("Expected igraph object, got [", class(graph), "] instead!")))
  }
  # validate method parameter
  if(!(method %in% c("min", "max"))) {
    stop(simpleError("knn method must be either 'min' or 'max'!"))
  }
  # validate edge mode parameter
  if(!(edge_mode %in% c("all", "out", "in"))) {
    stop(simpleError("Edge direction must be 'in', 'out', or 'all'!"))
  }
  # vector for recording whether to keep an edge
  keep_edge <- vector(mode = "logical", length = length(igraph::E(graph)))

  # for each vertex in the graph
  for(i in seq_along(igraph::V(graph))) {
    v <- igraph::V(graph)[i]
    # get the edges connected to this vertex
    v_edges <- igraph::incident(graph, v, mode = edge_mode)
    # if we're keeping largest edge weights
    if(method == "max") {
      # rank edges smallest to largest
      edge_ranks <- rank(igraph::E(graph)[v_edges]$weight, ties.method = "max")
      # grab the highest k ranks (aka the largest weights)
      k_idx <- which(edge_ranks > (length(v_edges) - k))
    } else { # keeping smallest weights
      # rank edges smallest to largest
      edge_ranks <- rank(igraph::E(graph)[v_edges]$weight, ties.method = "min")
      # grab the lowest k ranks (aka the smallest weights)
      k_idx <- which(edge_ranks <= k)
    }
    # mark these edges to be kept
    keep_edge[v_edges[k_idx]] <- TRUE
  }

  # which edges are we NOT going to keep?
  drop_edges <- igraph::E(graph)[which(!keep_edge)]
  # delete these edges
  knn_graph <- igraph::delete_edges(graph, drop_edges)

  return(knn_graph)
}
