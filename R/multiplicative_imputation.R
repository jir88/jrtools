#' Impute zeroes in compositional data
#'
#' Impute zeroes in compositional data using the multiplicative replacement
#' strategy described by Martin-Fernandez et al. (2003). Data will be closed
#' before zero replacement.
#'
#' @param X count or relative abundance matrix where rows are samples and columns are features
#' @param imp_factor small proportion to replace zero values with (default 1e-11)
#'
#' @references
#' Martin-Fernandez, J. A., Barcelo-Vidal, C., & Pawlowsky-Glahn, V. (2003). Dealing
#' with Zeros and Missing Values in Compositional Data Sets Using Nonparametric
#' Imputation. Mathematical Geology, 26.
#'
#' @return A matrix with data closed and imputed zeroes
#'
#' @export
multiplicative_imputation <- function(X, imp_factor = 1e-11){
  # close the data
  rs <- rowSums(X)
  if(sum(rs == 0) > 0) {
    warning("One or more rows have sum==0! This may generate unexpected behavior!")
  }
  nz <- rs != 0
  X[nz, ] <- X[nz, ]/rs[nz]

  if(is.null(imp_factor)){
    imp_factor <- min(X[X>0])/10
  }

  data_imputed <- apply(X = X, MARGIN = 1, FUN = function(row_data) {
    nz <- sum(row_data == 0)
    # row is all zeroes, undefined imputation since we don't actually know if
    # proportions of zero values are equal
    if(nz == length(row_data)) {
      return(rep_len(NA, nz))
    }
    imp_row <- row_data
    # impute zeroes
    imp_row[row_data==0] <- imp_factor
    # scale non-zero values
    imp_row[row_data!=0] <- imp_row[row_data!=0]*(1 - imp_factor*nz)
    return(imp_row)
  })
  # flip back to right orientation
  data_imputed <- t(data_imputed)
  return(data_imputed)
}
