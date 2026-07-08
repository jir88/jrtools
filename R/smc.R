
#' Significant Multivariate Correlation (SMC)
#'
#' Computes Significant Multivariate Correlation (SMC) F-statistics
#' for each variable in a dataset based on projection onto a regression
#' coefficient vector.
#'
#' @param b A numeric vector of regression coefficients (e.g., from PLS).
#' @param X A numeric matrix of predictors with dimensions \eqn{n \times p}
#'   (samples × variables).
#' @param alpha The desired significance level at which to calculate critical F-values.
#'
#' @details
#' The function projects the data matrix \code{X} onto the direction defined
#' by \code{b}, decomposing \code{X} into captured and residual components.
#' It then computes F-statistics to assess the significance of each variable's
#' contribution.
#'
#' Degrees of freedom:
#' \itemize{
#'   \item Captured variance: 1
#'   \item Residual variance: \eqn{n - 2}
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{smcF}{Numeric vector of SMC F-values for each variable.}
#'   \item{smcFcrit}{Critical F-value at significance level \eqn{\alpha = 0.05}.}
#'   \item{SSCapture}{Sum of squares captured for each variable.}
#'   \item{SSResidual}{Residual sum of squares for each variable.}
#' }
#'
#' @references
#' Thanh N. Tran (2014). Interpretation of variable importance in Partial Least
#' Squares with Significance Multivariate Correlation (sMC). DOI:10.1016/j.chemolab.2014.08.005
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(100 * 5), 100, 5)
#' b <- rnorm(5)
#' result <- smc(b, X)
#' result$smcF
#'
#' @export
smc <- function(b, X, alpha = 0.05) {

  n <- nrow(X)

  # Ensure b is a column vector
  b <- as.matrix(b)

  # Predicted response
  yhat <- X %*% b

  # Projection of X onto b
  Xhat <- (yhat %*% t(b)) / (norm(b, type = "F")^2)

  # Residual matrix
  Xresidual <- X - Xhat

  # Sum of squares
  SSCapture <- colSums(Xhat^2)
  SSResidual <- colSums(Xresidual^2)

  # Mean squares
  MSCapture <- SSCapture  # 1 degree of freedom
  MSResidual <- SSResidual / (n - 2)

  # F statistic
  smcF <- MSCapture / MSResidual

  # Critical F value
  smcFcrit <- stats::qf(1 - alpha, df1 = 1, df2 = n - 2)

  return(list(
    smcF = smcF,
    smcFcrit = smcFcrit,
    SSCapture = SSCapture,
    SSResidual = SSResidual
  ))
}
