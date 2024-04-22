#' Plot two mass spectra
#'
#' Plots two mass spectra back-to-back for comparison. Fragments within a
#' given PPM tolerance of each other are marked.
#'
#' @param spec1 Spectrum to plot on top half of graph. Must have columns named
#'   'mz' and 'intensity'.
#' @param spec2 Spectrum to plot on bottom half of graph. Must have columns named
#'   'mz' and 'intensity'.
#' @param match_tol_ppm Fragment matching tolerance in PPM
#' @param label_thresh Minimum intensity for labeled fragments. Set to NULL for
#'   no labels.
#'
#' @return The ggplot object for rendering or further modification.
#'
#' @importFrom rlang .data
#' @export
plot_spectral_comparison <- function(spec1, spec2, match_tol_ppm = 5,
                                     label_thresh = NULL) {
  # select only the m/z and intensity columns to avoid issues with other columns
  spec1 <- as.matrix(spec1[, c("mz", "intensity"), drop = FALSE])
  spec2 <- as.matrix(spec2[, c("mz", "intensity"), drop = FALSE])

  # look for matching peaks
  tol_low <- 1.0 - match_tol_ppm/1e6
  tol_high <- 1.0 + match_tol_ppm/1e6

  # sort spectrum 2 if it wasn't already
  spec2 <- spec2[order(spec2[, "mz"]), , drop = FALSE]

  # locate fragment indices all at once
  s2_idx_low <- findInterval(spec1[, "mz"]*tol_low,
                             spec2[, "mz"])
  s2_idx_high <- findInterval(spec1[, "mz"]*tol_high,
                              spec2[, "mz"])

  # which fragments' tolerance ranges in spec1 intersect 1+ fragments in spec2?
  s1_matches <- (s2_idx_high - s2_idx_low) > 0
  # indices of matching spec1 fragments
  s1_match_idx <- which(s1_matches)

  # which spec2 fragments match spec1 fragments?
  s2_matches <- logical(length = nrow(spec2))
  # for each spec1 fragment that matches 1+ spec2 fragment(s)
  for(i in s1_match_idx) {
    # get indices of spec2 fragments
    idx <- seq.int(from = s2_idx_low[i] + 1,
                   to = s2_idx_high[i])
    # mark matches
    s2_matches[idx] <- TRUE
  }

  # put spectrum 2 on the lower half
  spec2[, "intensity"] <- -spec2[, "intensity"]
  # combine spectra
  combo_spec <- tibble::as_tibble(rbind(spec1, spec2))
  # mark matching fragments
  combo_spec$Matches = c(s1_matches, s2_matches)

  # build the plot
  g <- ggplot2::ggplot(combo_spec, ggplot2::aes(x = .data[["mz"]],
                                                xend = .data[["mz"]],
                                                y = .data[["intensity"]],
                                                yend = 0,
                                                color = .data[["Matches"]])) +
    ggplot2::geom_hline(yintercept = 0, color = "darkgrey") +
    ggplot2::geom_segment() +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue"))

  # if we're adding labels
  if(is.numeric(label_thresh)) {
    # which peaks are above the threshold?
    feat_labs <- dplyr::filter(combo_spec, abs(.data[["intensity"]]) > label_thresh)
    # create labels and set positions based on which spectrum they are in
    feat_labs <- dplyr::mutate(feat_labs, Label = as.character(round(.data[["mz"]], 4)),
                               vjust = dplyr::if_else(.data[["intensity"]] > 0,
                                                      0, 1))

    g <- g + ggplot2::geom_text(data = feat_labs,
                                mapping = ggplot2::aes(label = .data[["Label"]],
                                                       vjust = .data[["vjust"]]),
                                check_overlap = TRUE,
                                color = "darkgrey", size = 3)
  }
  return(g)
}
