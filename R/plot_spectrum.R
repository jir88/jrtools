#' Plot a single mass spectrum
#'
#' Plots one mass spectrum with fragments labeled above a certain threshold
#'
#' @param spec Spectrum to plot. Must have columns named
#'   'mz' and 'intensity'.
#' @param label_thresh Minimum intensity for labeled fragments. Set to NULL for
#'   no labels.
#'
#' @return The ggplot object for rendering or further modification.
#'
#' @importFrom rlang .data
#' @export
plot_spectrum <- function(spec, label_thresh = NULL) {
  # select only the m/z and intensity columns to avoid issues with other columns
  spec <- tibble::as_tibble(spec[, c("mz", "intensity")])

  # peak labels
  if(is.numeric(label_thresh)) {
    # only peaks above threshold
    feat_labs <- dplyr::filter(spec, .data[["intensity"]] > label_thresh)
    feat_labs <- dplyr::mutate(feat_labs, Label = as.character(round(.data[["mz"]], 4)))
  }

  # return the plot
  g <- ggplot2::ggplot(data = spec,
                       mapping = ggplot2::aes(x = .data[["mz"]],
                                              xend = .data[["mz"]],
                                              y = .data[["intensity"]],
                                              yend = 0)) +
    ggplot2::geom_hline(yintercept = 0, color = "darkgrey") +
    ggplot2::geom_segment(color = "black") +
    ggplot2::xlab("m/z")

  # if we're adding labels
  if(is.numeric(label_thresh)) {
    g <- g + ggplot2::geom_text(data = feat_labs,
                                mapping = ggplot2::aes(label = .data[["Label"]]),
                                check_overlap = TRUE,
                                color = "darkgrey", vjust = 0, size = 3)
  }
  return(g)
}
