#' Extract gap fill status information for compounds in a Compound Discoverer alignment
#'
#' Compound Discoverer mass spec data alignment databases store gap filling data
#' for all consolidated compounds as binary blobs. This function extracts that
#' information and converts it into human readable factor data. Note that this
#' data is only for the reference/best ion in each compound. Information for
#' other ions must be extracted from the enormous MissingCompoundIonInstanceItems
#' table.
#'
#' The alignment offers two levels of status information. The simple one has
#' levels "No gap", "Missing ions", and "Full gap". The more detailed status
#' includes information on exactly how any gaps were filled. In the detailed
#' status information, 'unknown status' seems to correspond to features that
#' did not require gap filling.
#'
#' @note
#' Compound Discoverer software is produced by Thermo Fisher Scientific. This
#' package is not affiliated with Thermo Fisher Scientific in any way. For an
#' official Python interface to Compound Discoverer alignment files, see
#' \url{https://github.com/thermofisherlsms/pyeds}
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get gap status info for, or NULL to return all
#' @param detailed Should detailed gap filling status be returned? Defaults to
#'   simple status.
#'
#' @return A tibble containing the gap fill status data
#'
#' @importFrom rlang .data
#' @export
extract_gap_fill_status <- function(msa, ids = NULL, detailed = FALSE) {
  # blobs contain alternating 32-bit integer codes and a binary byte that seems
  # to always be '1', so we'll ignore it

  # use all IDs
  if(is.null(ids)) {
    ids <- msa$unknown_compound_items$ID
  }
  # usually index and ID are the same, but we can't be sure
  idx <- match(ids, msa$unknown_compound_items$ID)
  # get only study file IDs that were actually integrated (Sample, QC, and Blank)
  sfids <- msa$input_files$StudyFileID[which(msa$input_files$SampleType != "Identification Only")]

  # get the desired level of detail
  if(detailed) {
    gap_blobs <- msa$unknown_compound_items$GapFillStatus[idx]
    # sometimes compounds aren't gap filled, so deal with those
    found_status_idx <- which(!sapply(gap_blobs, is.null))
    gap_blobs <- gap_blobs[found_status_idx]
  } else {
    gap_blobs <- msa$unknown_compound_items$GapStatus[idx]
    # excluded compounds don't get integrated, so we need to deal with those
    found_status_idx <- which(!sapply(gap_blobs, is.null))
    gap_blobs <- gap_blobs[found_status_idx]
  }
  gap_data <- sapply(X = gap_blobs,
                      FUN = function(blb) {
                        status_values <- blb[-seq.int(from = 5, to = length(blb), by = 5)]
                        status_values <- readBin(status_values, what = integer(),
                                                 size = 4, # 32-bit int
                                                 n = length(status_values)/4)
                        return(status_values)
                      })
  # gap_data <- t(gap_data)


  # areas <- do.call(rbind, gap_data$area)

  # if there are any IDs with un-integrated peaks
  if(length(found_status_idx) < length(idx)) {
    df <- matrix(NA_real_, ncol = length(idx), nrow = nrow(gap_data))
    df[, found_status_idx] <- gap_data
    # assign expanded matrix
    gap_data <- df
  }

  colnames(gap_data) <- ids
  gap_data <- tibble::as_tibble(gap_data, .name_repair = "minimal")
  gap_data <- dplyr::mutate(gap_data, StudyFileID = sfids,
                         .before = 1)

  # convert numeric codes to factors
  if(detailed) {
    # 0: "Unknown status",
    #     1: "No gap to fill",
    #     2: "Unable to fill",
    #     4: "Filled by arbitrary value",
    #     8: "Filled by trace area",
    #     16: "Filled by simulated peak",
    #     32: "Filled by spectrum noise",
    #     64: "Filled by matching ion",
    #     128: "Filled by re-detected peak",
    #     256: "Imputed by low area value",
    #     512: "Imputed by group median",
    #     1024: "Imputed by Random Forest",
    #     2048: "Skipped"
    levs <- c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
    labs <- c("Unknown status",
              "No gap to fill",
              "Unable to fill",
              "Filled by arbitrary value",
              "Filled by trace area",
              "Filled by simulated peak",
              "Filled by spectrum noise",
              "Filled by matching ion",
              "Filled by re-detected peak",
              "Imputed by low area value",
              "Imputed by group median",
              "Imputed by Random Forest",
              "Skipped")
  } else {
    levs <- c(1, 2, 3)
    labs <- c("No gap", "Full gap", "Missing ions")
  }
  gap_data <- dplyr::mutate(gap_data,
                            dplyr::across(.cols = -.data$StudyFileID,
                                          .fns = ~ factor(.x,
                                                          levels = levs,
                                                          labels = labs)))

  return(gap_data)
}
