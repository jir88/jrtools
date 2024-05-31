#' Extract peak quality blobs in a Compound Discoverer alignment
#'
#' Compound Discoverer cdResult alignment databases store peak quality ratings
#' as SQLite blobs containing the data in binary form. This function extracts the
#' stored data from such blobs. The returned ratings are for the reference adduct
#' ion for each compound. Other areas must be retrieved manually from the relevant
#' tables.
#'
#' @note
#' Compound Discoverer software is produced by Thermo Fisher Scientific. This
#' package is not affiliated with Thermo Fisher Scientific in any way. For an
#' official Python interface to Compound Discoverer alignment files, see
#' \url{https://github.com/thermofisherlsms/pyeds}
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get ratings for, or NULL to get all peak ratings
#'
#' @return A tibble containing the peak ratings
#'
#' @export
extract_peak_ratings <- function(msa, ids = NULL) {
  # TODO: allow extracting sub-scores as well as combined rating

  # blobs contain alternating 64-bit floats and a binary byte specifying whether
  # the area was calculated. Presumably this is to differentiate between true
  # zeroes and NA values?

  # use all IDs
  if(is.null(ids)) {
    ids <- msa$unknown_compound_items$ID
  }
  # usually index and ID are the same, but we can't be sure
  idx <- match(ids, msa$unknown_compound_items$ID)
  # get only study file IDs that were actually integrated (Sample, QC, and Blank)
  sfids <- msa$input_files$StudyFileID[which(msa$input_files$SampleType != "Identification Only")]

  rating_blobs <- msa$unknown_compound_items$PeakRating[idx]

  rating_data <- lapply(X = rating_blobs,
                        FUN = function(blb) {
                          rating_values <- blb[-seq.int(from = 9, to = length(blb), by = 9)]
                          rating_values <- readBin(rating_values, what = numeric(), n = length(rating_values)/8)

                          # was this peak actually integrated in this file?
                          rating_flags <- blb[seq.int(from = 9, to = length(blb), by = 9)]
                          rating_flags <- readBin(rating_flags, what = logical(),
                                                size = 1, n = length(rating_flags))
                          return(list(rating = rating_values,
                                      flag = rating_flags))
                        })
  rating_data <- purrr::transpose(rating_data)

  ratings <- do.call(rbind, rating_data$rating)

  rownames(ratings) <- ids
  ratings <- tibble::as_tibble(t(ratings), .name_repair = "minimal")
  ratings <- dplyr::mutate(ratings, StudyFileID = sfids,
                         .before = 1)

  return(ratings)
}
