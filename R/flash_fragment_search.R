#' Quickly search spectral library for fragment pattern
#'
#' Given a library of compound spectra and a set of fragment ions, this function
#' identifies all compounds containing one or more of the fragment ions. By
#' converting the library and fragment ion m/z values to neutral losses, one can
#' also do neutral loss searching
#'
#' Current implementation does not consider the relative fragment intensities.
#'
#' @param fragment_library matrix of compound spectra sorted by fragment m/z.
#'   Must have columns "Spectrum", "mz", and "intensity"
#' @param query_ions matrix with query fragment ions. Must have columns "mz" and
#'   "intensity"
#' @param ms2_tol_ppm Fragment matching tolerance in PPM
#'
#' @return A tibble containing the number of query ions found in each compound
#'   spectrum, plus the total intensity of those ions.
#'
#' @export
flash_fragment_search <- function(fragment_library, query_ions, ms2_tol_ppm = 5.0) {
  ms2_tol_low <- 1.0 - ms2_tol_ppm/1e6
  ms2_tol_high <- 1.0 + ms2_tol_ppm/1e6

  frag_library_mz <- fragment_library[, "mz"]

  # pre-allocate table
  frag_tab_spec_id <- sort(unique(fragment_library[, "Spectrum"]))
  n_ids <- length(frag_tab_spec_id)
  frag_tab_sim <- integer(length = n_ids)
  frag_tab_match_int <- numeric(length = n_ids)

  # locate query ion indices all at once
  match_idx_low <- findInterval(query_ions[, "mz"]*ms2_tol_low,
                                frag_library_mz)
  match_idx_high <- findInterval(query_ions[, "mz"]*ms2_tol_high,
                                 frag_library_mz)

  # which query ions have library fragment matches?
  query_match_idx <- which((match_idx_high - match_idx_low) > 0)

  # for each query ion with matches
  for(i in query_match_idx) {
    # find matching compound spectrum fragments
    match_idx <- (match_idx_low[i] + 1):match_idx_high[i]
    # which spectra contain each matching fragment?
    match_spectra <- fragment_library[match_idx, "Spectrum"]
    match_int <- fragment_library[match_idx, "intensity"]

    # increment match counts
    sim_idx <- match(match_spectra, frag_tab_spec_id)
    frag_tab_sim[sim_idx] <- frag_tab_sim[sim_idx] + 1
    frag_tab_match_int[sim_idx] <- frag_tab_match_int[sim_idx] + match_int
  }

  similarity_table <- tibble::tibble(Spectrum = frag_tab_spec_id,
                                     Matches = frag_tab_sim,
                                     MatchIntensity = frag_tab_match_int)
  return(similarity_table)
}
