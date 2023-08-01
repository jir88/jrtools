#' Match two sets of mass and retention time pairs
#'
#' Compares two sets of mass/retention time pairs and finds any matches within
#' a specified range of mass and retention time differences. WARNING: make sure
#' that both sets are using the same measure of feature mass. This function will
#' happily compare m/z values to predicted neutral masses and return bogus
#' results, since there's no way to differentiate them.
#'
#' @param mz1 First set of mass or m/z values
#' @param rt1 First set of retention times
#' @param mz2 Second set of mass or m/z values
#' @param rt2 Second set of retention times
#' @param ppm_range A length-two numeric giving the minimum and maximum allowed
#' differences between feature masses in parts per million
#' @param drt_range A length-two numeric giving the minimum and maximum allowed
#' differences between feature retention times
#'
#' @return A tibble of feature matches
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
match_ms_features <- function(mz1, rt1, mz2, rt2, ppm_range = c(-5, 5),
                             drt_range = c(-0.5, 0.5)) {
  # find close masses (within 25 ppm)
  delta_mass <- outer(mz1, mz2,
                      FUN = function(x, y) { return((x - y)/y*1e6) })
  # colnames(delta_mass) <- stringr::str_c(mdro_acquirex_features$Mass, "@",
  #                                        mdro_acquirex_features$RT)
  colnames(delta_mass) <- stringr::str_c(mz2, "@", rt2)
  rownames(delta_mass) <- stringr::str_c(mz1, "@", rt1)

  feature_comparison <- tibble::as_tibble(delta_mass, rownames = "Feature1",
                                          .name_repair = "minimal") %>%
    tidyr::pivot_longer(cols = -.data$Feature1, names_to = "Feature2", values_to = "ppm") %>%
    dplyr::filter(dplyr::between(.data$ppm, ppm_range[1], ppm_range[2])) %>%
    dplyr::mutate(mz1 = mz1[match(.data$Feature1, rownames(delta_mass))],
                  rt1 = rt1[match(.data$Feature1, rownames(delta_mass))]) %>%
    dplyr::mutate(mz2 = mz2[match(.data$Feature2, colnames(delta_mass))],
                  rt2 = rt2[match(.data$Feature2, colnames(delta_mass))]) %>%
    dplyr::mutate(dRT = rt1 - rt2) %>%
    dplyr::filter(dplyr::between(.data$dRT, drt_range[1], drt_range[2])) %>%
    dplyr::relocate(.data$ppm, .before = .data$dRT)

  return(feature_comparison)
}
