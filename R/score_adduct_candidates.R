#' Score all possible adduct identities for an MS1 peak
#'
#' Given a peak from an MS1 spectrum, tests how well different adduct assignments
#' match the MS1 spectrum.
#'
#' @param mz The m/z ratio of interest
#' @param mz_polarity The polarity of the m/z ratio. Must be either "positive" or
#'  "negative"
#' @param ms1_spec The MS1 spectrum being examined. Must be tibble-alike with
#'  columns "polarity", "mz", and "intensity"
#' @param all_adducts The possible adducts to be considered. Must be tibble-alike
#'  with columns "polarity", "Ion name", "Power", "Mult", and "Mass".
#' @param ppm_thresh The maximum difference between observed m/z values and
#'  proposed adduct m/z values
#'
#' @return A tibble of adduct candidate match scores.
#'
#' @importFrom rlang .data
#' @export
score_adduct_candidates <- function(mz, mz_polarity, ms1_spec, all_adducts,
                                    ppm_thresh=15) {
  # get candidate adduct list matching ion polarity
  adduct_candidates <- dplyr::filter(all_adducts, .data$polarity == mz_polarity)
  # for each candidate adduct
  oap_match_data <- lapply(X = seq_along(adduct_candidates$`Ion name`), FUN = function(j) {
    # calculate proposed adduct molecular weight
    prop_mw <- (mz - adduct_candidates$Mass[j])/(adduct_candidates$Mult[j]^adduct_candidates$Power[j])
    # predict other adduct m/z ratios in both polarities
    other_adduct_mz <- prop_mw*all_adducts$Mult^all_adducts$Power + all_adducts$Mass
    # calculate number of matching peaks in MS1 spectrum
    # find close masses (within 25 ppm)
    other_adduct_ppm <- outer(other_adduct_mz, ms1_spec$mz,
                              FUN = function(x, y) { return((x - y)/y*1e6) })
    # where are the hits?
    df <- which(abs(other_adduct_ppm) < ppm_thresh, arr.ind = TRUE)
    # number of unique rows is the number of other adducts that match 1+ MS1 peaks
    oap_match_idx <- unique(df[, "row"])
    num_oap_matches <- length(oap_match_idx)
    # calculate total intensity associated with all matched MS1 peaks
    oap_match_intensity <- sum(ms1_spec$intensity[unique(df[, "col"])])
    return(list("adduct_idx" = oap_match_idx,
                "adduct_names" = all_adducts$`Ion name`[oap_match_idx],
                "adduct_polarities" = all_adducts$polarity[oap_match_idx],
                "adduct_mz" = other_adduct_mz[oap_match_idx],
                "num_matches" = num_oap_matches,
                "match_intensity" = oap_match_intensity))
  })
  # proposed adduct that explains the most other adducts is 'best'
  # break ties using total explained intensity
  df2 <- tibble::as_tibble(purrr::list_transpose(oap_match_data))
  df2$Candidate = adduct_candidates$`Ion name`
  df2 <- dplyr::arrange(df2, dplyr::desc(.data$num_matches), dplyr::desc(.data$match_intensity))
  # convert to tidy arrangement
  df2 <- dplyr::mutate(df2, candidate_idx = dplyr::row_number(), .before = .data$adduct_idx)
  df2 <- tidyr::unnest(df2, cols = c("adduct_idx", "adduct_names", "adduct_polarities", "adduct_mz"))
  df2 <- dplyr::arrange(df2, .data$candidate_idx, .data$adduct_polarities, .data$adduct_mz)
  return(df2)
}
