#' Check if object is a date
#'
#' @keywords internal
#' @param x vector
is_date <- function(x) {
  inherits(x, c("POSIXt", "POSIXct", "POSIXlt", "Date"))
}

#' Correlation value plot
#'
#' Estimate correlation from the given data. If a color variable is supplied, the correlation will also be calculated per group.
#'
#' @param data data set using
#' @param mapping aesthetics being used
#' @param ... other arguments being supplied to \code{\link[ggplot2]{geom_text}()} for the title and groups
#' @param stars logical value which determines if the significance stars should be displayed.  Given the \code{\link[stats]{cor.test}} p-values, display \describe{
#'   \item{\code{"***"}}{if the p-value is \verb{< 0.001}}
#'   \item{\code{"**"}}{if the p-value is \verb{< 0.01}}
#'   \item{\code{"*"}}{if the p-value is \verb{< 0.05}}
#'   \item{\code{"."}}{if the p-value is \verb{< 0.10}}
#'   \item{\code{""}}{otherwise}
#' }
#' @param p_adjust_method method for adjusting correlation p-values if desired. Name passed to \code{\link[stats]{p.adjust}}.
#' @param p_adjust_n number of comparisons for p-value adjustment, because we can't access that info within this function.
#' @param method \code{method} supplied to cor function
#' @param use \code{use} supplied to \code{\link[stats]{cor}} function
#' @param display_grid if \code{TRUE}, display aligned panel grid lines. If \code{FALSE} (default), display a thin panel border.
#' @param digits number of digits to be displayed after the decimal point. See \code{\link[base]{formatC}} for how numbers are calculated.
#' @param title_args arguments being supplied to the title's \code{\link[ggplot2]{geom_text}()}
#' @param group_args arguments being supplied to the split-by-color group's \code{\link[ggplot2]{geom_text}()}
#' @param justify_labels \code{justify} argument supplied when \code{\link[base]{format}}ting the labels
#' @param align_percent relative align position of the text. When \code{justify_labels = 0.5}, this should not be needed to be set.
#' @param alignPercent,displayGrid deprecated. Please use their snake-case counterparts.
#' @param title title text to be displayed
#' @author Barret Schloerke
#' @seealso \code{\link[GGally]{ggally_statistic}}, \code{\link[GGally]{ggally_cor_v1_5}}, \code{\link[GGally]{ggally_cor}}
#' @export
#' @keywords hplot
ggally_adj_cor <- function(
    data,
    mapping,
    ...,
    stars = TRUE,
    p_adjust_method = "none",
    p_adjust_n = NA,
    method = "pearson",
    use = "complete.obs",
    display_grid = FALSE,
    digits = 3,
    title_args = list(...),
    group_args = list(...),
    justify_labels = "right",
    align_percent = 0.5,
    title = "Corr",
    alignPercent = warning("deprecated. Use `align_percent`"),
    displayGrid = warning("deprecated. Use `display_grid`")
) {
  if (!missing(alignPercent)) {
    warning("`alignPercent` is deprecated. Please use `align_percent` if alignment still needs to be adjusted")
    align_percent <- alignPercent
  }
  if (!missing(displayGrid)) {
    warning("`displayGrid` is deprecated. Please use `display_grid`")
    display_grid <- displayGrid
  }
  if(p_adjust_method != "none" & is.na(p_adjust_n)) {
    stop("Number of comparisons must be specified if adjusting p-values for multiple comparisons!")
  }

  na.rm <-
    if (missing(use)) {
      # display warnings
      NA
    } else {
      (use %in% c("complete.obs", "pairwise.complete.obs", "na.or.complete"))
    }

  GGally::ggally_statistic(
    data = data,
    mapping = mapping,
    na.rm = na.rm,
    align_percent = align_percent,
    display_grid = display_grid,
    title_args = title_args,
    group_args = group_args,
    justify_labels = justify_labels,
    justify_text = "left",
    sep = if ("colour" %in% names(mapping)) ": " else ":\n",
    title = title,
    text_fn = function(x, y) {
      if (is_date(x)) {
        x <- as.numeric(x)
      }
      if (is_date(y)) {
        y <- as.numeric(y)
      }

      corObj <- stats::cor.test(x, y, method = method, use = use)

      # make sure all values have X-many decimal places
      cor_est <- as.numeric(corObj$estimate)
      cor_txt <- formatC(cor_est, digits = digits, format = "f")

      # if we need to adjust the p-values
      if(p_adjust_method != "none") {
        corObj$p.value <- stats::p.adjust(p = corObj$p.value,
                                          method = p_adjust_method,
                                          n = p_adjust_n)
      }

      # if stars should be added
      if (isTRUE(stars)) {
        cor_txt <- stringr::str_c(
          cor_txt,
          GGally::signif_stars(corObj$p.value)
        )
      }

      cor_txt
    }
  )
}
