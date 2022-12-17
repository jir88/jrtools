#' Change the order of clusters in existing K-means model
#'
#' Changes the order (and associated IDs) of the clusters defined by an
#' existing K-means model.
#'
#' @param x a fitted K-means model object
#' @param order a vector of the same length as the number of clusters, giving
#' the new order of the current cluster IDs
#' @param ... not used
#'
#' @return A reordered copy of the K-means model object
#'
#' @importFrom stats reorder
#' @export
reorder.kmeans <- function(x,
                           order,
                           ...) {
  # check that number of clusters is correct
  if(length(order) != length(x$size)) {
    stop("Parameter 'order' must have same length as number of clusters!")
  }
  # check that cluster IDs are in range
  if(!setequal(unique(x$cluster), unique(order))) {
    stop("Cluster ID values must match existing cluster IDs!")
  }

  # swap cluster assignments
  swapped_cluster <- x$cluster
  for(i in seq_along(order)) {
    cid <- order[i]
    swapped_cluster[x$cluster == cid] <- i
  }
  x$cluster <- swapped_cluster

  # reorder cluster centers
  x$centers <- x$centers[order, ]
  rownames(x$centers) <- sort(unique(order))

  # reorder other outputs
  x$withinss <- x$withinss[order]
  x$size <- x$size[order]

  return(x)
}
