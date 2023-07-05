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

#' Find connected k-plexes in networks
#'
#' Given an igraph object, this function identifies connected components of
#' induced subgraphs where all N nodes in the component are connected to at
#' least N-k other nodes. Thus, a 1-plex is the same as a clique of size N,
#' because all nodes are connected to each other but not to themselves. A 2-plex
#' is a clique where some nodes are missing at most one edge to another node.
#'
#' WARNING: this method uses an untested algorithm! Use at your own risk! Only
#' works for connected k-plexes! May not find maximal k-plexes or all k-plexes.
#' Requires more formalized testing.
#'
#' @param g the igraph object to search for k-plexes
#' @param min_kc The minimum coreness to consider
#' @param kp the maximum number of missing edges per node in the plex
#' @return a list containing a list of the k-plex graphs, the minimum degrees of
#'   each plex, and the sizes of each plex
#' @details This is a really simplistic kludge. Use at your own risk! May not
#'   behave as intended!
#' @export
find_kplexes <- function(g, min_kc, kp) {
  sdf_coreness <- igraph::coreness(g)
  # V(g)$coreness <- sdf_coreness

  plex_list <- c()

  # for each observed k-core
  for(kc in seq.int(from = max(sdf_coreness), to = min_kc, by = -1)) {
    # get the k-core
    kcore <- igraph::induced_subgraph(g, igraph::V(g)[sdf_coreness >= kc])
    # get a list of the connected components
    kcore_membership <- igraph::components(kcore)$membership
    df <- rle(sort(kcore_membership))
    # a component is a k-plex if its size equals kc + kp
    plex_idx <- which(df$lengths == kc + kp)
    # if there are any detected plexes
    if(length(plex_idx) > 0) {
      # for each detected plex
      pl <- lapply(X = df$values[plex_idx], FUN = function(comp_id) {
        # grab the induced subgraph associated with the plex component's ID
        return(igraph::induced_subgraph(kcore, igraph::V(kcore)[comp_id == kcore_membership]))
      })
      plex_list <- append(plex_list, pl)
    }

  }

  # what is the minimum degree of the nodes in each plex?
  plex_min_deg <- sapply(X = plex_list, FUN = function(g) { return(min(igraph::degree(g)))})
  # how big is each plex?
  plex_size <- sapply(X = plex_list, FUN = function(g) { return(length(g))})

  return(list(plex_list = plex_list, plex_min_deg = plex_min_deg,
              plex_size = plex_size))
}
