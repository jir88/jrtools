#' Memory-efficient search for similar numeric values
#'
#' Given a sorted numeric vector, this function finds all pairs of values that
#' are within a threshold difference of each other.
#'
#' @param x sorted numeric vector to be grouped
#' @param maxDiff maximum difference between pairs
#' @param ppm FALSE (default) to use absolute differences, or TRUE to use
#'   parts-per-million differences instead
#'
#' @return A two-column matrix giving the indices of all pairs.
#'
#' @export
find_sim_pairs <- function (x, maxDiff, ppm = FALSE)
{
  # if x hasn't been sorted, whine about it
  if (is.unsorted(x)) {
    stop("x must be sorted in ascending order!")
  }
  # if maxDiff is negative, whine about it
  if (maxDiff < 0) {
    stop("maxDiff must be non-negative!")
  }

  # Initialize an empty list to store the pairs
  pairs <- list()
  # pair list index
  p_idx <- 1

  # Iterate through the sorted vector
  for (i in seq_along(x)) {
    # Find the values within the threshold of the current value
    if(ppm) {
      within_threshold <- which(abs(x - x[i])/x[i]*1e6 < maxDiff)
    } else {
      within_threshold <- which(abs(x - x[i]) < maxDiff)
    }

    # Remove current and prior indices from the result
    within_threshold <- within_threshold[within_threshold > i]

    # If there are any values within the threshold, add them to the pairs list
    for(j in within_threshold) {
      pairs[[p_idx]] <- c(i, j)
      p_idx <- p_idx + 1
    }
  }
  pairs <- do.call(rbind, pairs)

  # Return the pairs
  return(pairs)
}
