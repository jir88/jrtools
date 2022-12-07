#' Calculate Shannon Entropy
#'
#' Calculates the Shannon entropy of a set of data.
#'
#' @param x The data set in a format acceptable for the table function
#'
#' @return the number of bits of entropy (log base 2) in the data set
#'
#' @export
entropy <- function(x) {
  freq <- table(x)/length(x)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  return(-sum(vec * log2(vec)))
}
