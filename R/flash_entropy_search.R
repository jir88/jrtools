#' Perform flash entropy search
#'
#' Calculates flash entropy similarity between a query spectrum and a library of
#' reference spectra. By converting the library and query m/z values to neutral
#' losses, one can also do neutral loss matching.
#'
#' Current implementation only searches one spectrum in 'open'/similarity mode.
#'
#' @param fragment_library matrix of library spectra sorted by fragment m/z. Must
#'   have columns "Spectrum", "mz", and "intensity"
#' @param query_spectrum matrix with query spectrum. Must have columns "mz" and
#'   "intensity"
#' @param ms2_tol_ppm Fragment matching tolerance in PPM
#'
#' @return A tibble containing spectral entropy similarity scores between the
#'   query spectrum and any library spectra with matching peaks.
#'
#' @export
flash_entropy_search <- function(fragment_library, query_spectrum, ms2_tol_ppm = 5.0) {
  # TODO: allow identity search (match precursors)
  ms2_tol_low <- 1.0 - ms2_tol_ppm/1e6
  ms2_tol_high <- 1.0 + ms2_tol_ppm/1e6

  # grab sub-library for faster searching
  min_mz <- min(query_spectrum[, "mz"])*ms2_tol_low
  max_mz <- max(query_spectrum[, "mz"])*ms2_tol_high
  idx <- findInterval(c(min_mz, max_mz), fragment_library[, "mz"]) + 1
  sub_library <- fragment_library[idx[1]:idx[2], ]
  sub_library_mz <- sub_library[, "mz"]

  # pre-allocate table
  sim_tab_spec_id <- sort(unique(sub_library[, "Spectrum"]))
  sim_tab_sim <- rlang::rep_along(sim_tab_spec_id, 0)

  # locate fragment indices all at once
  match_idx_low <- findInterval(query_spectrum[, "mz"]*ms2_tol_low,
                                sub_library_mz)
  match_idx_high <- findInterval(query_spectrum[, "mz"]*ms2_tol_high,
                                sub_library_mz)
  # for each query fragment
  for(i in 1:nrow(query_spectrum)) {
    frag_mz <- query_spectrum[i, "mz"]
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
