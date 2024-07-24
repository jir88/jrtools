#' Fast value clustering with PPM threshold
#'
#' This is a lightly adapted version of \link[MsFeatures]{groupConsecutive}
#' from the MsFeatures package. Instead of an absolute difference threshold, it
#' uses a PPM threshold.
#'
#' @param x sorted numeric vector to be grouped
#' @param maxPPM threshold below which features are clustered into the same group
#'
#' @return An integer vector of group ID assignments for each value in x.
#'
#' @export
groupConsecutivePPM <- function (x, maxPPM = 5)
{
  if (is.unsorted(x)) {
    idx <- order(x)
    x <- x[idx]
  } else {
    idx <- integer()
  }
  x_len <- length(x)
  x_groups <- rep(NA_integer_, x_len)
  i <- 1
  group_id <- 1
  while (any(is.na(x_groups))) {
    grp <- which(abs(x - x[i])/x[i]*1e6 <= maxPPM)
    not_in_prev_grp <- is.na(x_groups[grp])
    in_prev_grp <- grp[!not_in_prev_grp]
    if (length(in_prev_grp)) {
      i_diff <- abs(x[in_prev_grp] - mean(x[grp]))
      prev_grp <- x_groups[in_prev_grp]
      to_rem <- rep(FALSE, length(in_prev_grp))
      for (j in unique(prev_grp)) {
        j_diff <- abs(x[in_prev_grp] - mean(x[which(x_groups ==
                                                      j)]))
        to_rem <- to_rem | j_diff < i_diff
      }
      grp <- c(in_prev_grp[!to_rem], grp[not_in_prev_grp])
    }
    x_groups[grp] <- group_id
    group_id <- group_id + 1
    i <- which.max(is.na(x_groups))
  }
  x_groups[idx] <- x_groups
  x_groups
}
