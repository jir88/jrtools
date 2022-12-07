#' Assign new data to clusters in existing K-means model
#'
#' Assigns new samples to the clusters defined by an existing K-means model.
#'
#' @param object a fitted K-means model object
#' @param newdata a matrix-like object with one sample per row and the same
#' column order as the data used for the original K-means model
#' @param method whether to return the centers of the assigned clusters
#' ("centers") or the assigned cluster IDs ("classes")
#' @param ... not used
#'
#' @return An array containing either predicted cluster centers or cluster IDs
#'
#' @export
predict.kmeans <- function(object,
                           newdata,
                           method = c("centers", "classes"),
                           ...) {
  method <- match.arg(method)

  centers <- object$centers
  ss_by_center <- apply(centers, 1, function(x) {
    colSums((t(newdata) - x) ^ 2)
  })
  best_clusters <- apply(ss_by_center, 1, which.min)

  if (method == "centers") {
    return(centers[best_clusters, ])
  } else {
    return(best_clusters)
  }
}
