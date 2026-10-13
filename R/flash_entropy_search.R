#' Perform flash entropy search
#'
#' Calculates flash entropy similarity between a query spectrum and a library of
#' reference spectra. By converting the library and query m/z values to neutral
#' losses, one can also do neutral loss matching.
#'
#' Current implementation only searches one spectrum in 'open'/similarity mode.
#'
#' @param fragment_library matrix of library spectra sorted by fragment m/z. Must
#'   have columns "Spectrum", "mz", and "intensity". Individual spectra MUST have
#'   their total fragment intensity normalized to 0.5.
#' @param query_spectrum matrix with query spectrum. Must have columns "mz" and
#'   "intensity". Total intensity MUST be normalized to 0.5.
#' @param ms2_tol_ppm Fragment matching tolerance in PPM
#'
#' @return A tibble containing spectral entropy similarity scores between the
#'   query spectrum and any library spectra with matching peaks.
#'
#' @references
#' Li, Y., Kind, T., Folz, J. et al. Spectral entropy outperforms MS/MS dot product
#' similarity for small-molecule compound identification. Nat Methods 18,
#' 1524–1531 (2021). https://doi.org/10.1038/s41592-021-01331-z
#'
#' Li, Y., Fiehn, O. Flash entropy search to query all mass spectral libraries in
#' real time. Nat Methods 20, 1475–1478 (2023). https://doi.org/10.1038/s41592-023-02012-9
#' @export
flash_entropy_search <- function(fragment_library, query_spectrum, ms2_tol_ppm = 5.0) {
  # make sure query spectrum has been normalized
  query_intensity <- query_spectrum[, "intensity"]
  query_sum <- sum(query_intensity)
  if(query_sum != 0.5) {
    query_spectrum[, "intensity"] <- 0.5*query_intensity/query_sum
  }

  ms2_tol_low <- 1.0 - ms2_tol_ppm/1e6
  ms2_tol_high <- 1.0 + ms2_tol_ppm/1e6

  # pull query spectrum m/z ratios once, not multiple times
  query_mz <- query_spectrum[, "mz"]
  # grab sub-library for faster searching
  min_mz <- min(query_mz)*ms2_tol_low
  max_mz <- max(query_mz)*ms2_tol_high
  idx <- findInterval(c(min_mz, max_mz), fragment_library[, "mz"])
  sub_library <- fragment_library[(idx[1] + 1):idx[2], , drop = FALSE]
  sub_library_mz <- sub_library[, "mz"]

  # pre-allocate table
  sim_tab_spec_id <- unique(sub_library[, "Spectrum"])
  sim_tab_sim <- numeric(length = length(sim_tab_spec_id))

  # locate fragment indices all at once
  match_idx_low <- findInterval(query_mz*ms2_tol_low,
                                sub_library_mz)
  match_idx_high <- findInterval(query_mz*ms2_tol_high,
                                sub_library_mz)

  # which query spectrum fragments have library fragment matches?
  query_match_idx <- which((match_idx_high - match_idx_low) > 0)

  # for each query fragment with matches
  for(i in query_match_idx) {
    frag_int <- query_spectrum[i, "intensity"]
    # find matching library fragments
    match_idx <- (match_idx_low[i] + 1):match_idx_high[i]
    # which spectra contain each matching fragment?
    match_spectra <- sub_library[match_idx, "Spectrum"]

    # calculate similarity contributions
    match_int <- sub_library[match_idx, "intensity"]
    combo_int <- match_int + frag_int
    sim_contrib <- combo_int*log2(combo_int) -
      match_int*log2(match_int) -
      frag_int*log2(frag_int)

    # add to similarity table
    sim_idx <- match(match_spectra, sim_tab_spec_id)
    sim_tab_sim[sim_idx] <- sim_tab_sim[sim_idx] + sim_contrib
  }

  similarity_table <- tibble::tibble(Spectrum = sim_tab_spec_id,
                             Similarity = sim_tab_sim)
  return(similarity_table)
}
