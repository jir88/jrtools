#' Generate data for a volcano plot
#'
#' Given a dataset with two labels, calculates the fold-change and p-values
#' for each feature in the dataset.
#'
#' @param feature_data Feature values for all samples, where rows are samples
#'   and columns are features.
#' @param sample_labels Factor with two levels specifying which group each.sample
#'   belongs to.
#' @param test The statistical test function to use for calculating p-values.
#' @param p_adjust The type of multiple comparisons correction to perform.
#'
#' @return Tibble with columns 'feature', 'fold_change', 'p_value', and 'p_adj'
#'
#' @importFrom rlang .data
#' @export
calculate_volcano <- function(feature_data, sample_labels, test = stats::t.test,
                              p_adjust = "none") {
  lab_levels <- levels(sample_labels)
  if(length(lab_levels) != 2) {
    stop(paste("Sample labels must have 2 levels, got",
               length(lab_levels)))
  }

  lev1_idx <- which(sample_labels == lab_levels[1])
  lev2_idx <- which(sample_labels == lab_levels[2])

  t_func <- function(x) test(x = x[lev2_idx], y = x[lev1_idx])$p.value
  fold_func <- function(x) (mean(x[lev2_idx])/mean(x[lev1_idx]))

  vol_data <- sapply(X = feature_data,
                        FUN = t_func) %>%
    tibble::as_tibble(rownames = "feature") %>%
    dplyr::rename(p_value = .data[["value"]]) %>%
    dplyr::mutate(p_adj = stats::p.adjust(.data[["p_value"]], method = p_adjust)) %>%
    dplyr::mutate(fold_change = sapply(X = feature_data,
                                FUN = fold_func))

  return(vol_data)
}

#' Create a volcano plot
#'
#' Given a set of fold-changes and p-values, makes and returns a volcano plot.
#'
#' @param volcano_data Matrix or data frame with at least columns named
#'   'fold_change' and 'p_adj'
#' @param fc_thresh The fold-change threshold to display on the plot, or NA for none
#' @param p_value_thresh The p-value threshold to display on the plot, or NA for none
#'
#' @return A ggplot object containing the graph
#'
#' @importFrom rlang .data
#' @export
plot_volcano <- function(volcano_data, fc_thresh = NA, p_value_thresh = NA) {
  g <- ggplot2::ggplot(volcano_data, ggplot2::aes(x = log2(.data[["fold_change"]]),
                                         y = -log10(.data[["p_adj"]])))
  if(!is.na(fc_thresh)) {
    g <- g + ggplot2::geom_vline(xintercept = c(-1,1)*log2(fc_thresh),
                        color = "darkgrey", linetype = "dashed")
  }
  if(!is.na(p_value_thresh)) {
    g <- g + ggplot2::geom_hline(yintercept = -log10(p_value_thresh),
                        color = "darkgrey", linetype = "dashed")
  }
  g <- g +
    ggplot2::geom_point(alpha = 0.25) +
    ggplot2::labs(x = "log2(fold change)", y = "-log10(p-value)")
  return(g)
}
