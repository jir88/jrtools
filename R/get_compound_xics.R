#' Retrieve extracted ion chromatograms (XICs) associated with compounds
#'
#' Queries a Compound Discoverer mass spec alignment database to get all the
#' extracted ion chromatograms associated with some or all compounds. This
#' includes all identified adducts and isotopologues. Note that this only gets
#' XICs where a peak was detected. For gap-filled XICs, see \code{\link{get_gap_xics}}.
#'
#' @note
#' Compound Discoverer software is produced by Thermo Fisher Scientific. This
#' package is not affiliated with Thermo Fisher Scientific in any way. For an
#' official Python interface to Compound Discoverer alignment files, see
#' \url{https://github.com/thermofisherlsms/pyeds}
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get XICs for, or NULL to get XICs associated with
#'   all compounds
#'
#' @return A tibble with the XIC data
#'
#' @importFrom rlang .data
#' @export
get_compound_xics <- function(msa, ids = NULL) {
  if(is.null(ids)) {
    ids <- msa$unknown_compound_items$ID
  }

  # key of unconsolidated compound IDs associated with all compounds
  ucids <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsUnknownCompoundInstanceItems")
  # table of unconsolidated unknown compounds
  unk_comp_items <- dplyr::tbl(msa$db_connection, "UnknownCompoundInstanceItems")
  # key of ions associated with unknown compounds
  iids <- dplyr::tbl(msa$db_connection, "UnknownCompoundInstanceItemsUnknownCompoundIonInstanceItems")
  # table of unknown compound ions
  unk_comp_ion_items <- dplyr::tbl(msa$db_connection, "UnknownCompoundIonInstanceItems")
  # key of XIC traces associated with unknown compound ions
  xic_ids <- dplyr::tbl(msa$db_connection, "UnknownCompoundIonInstanceItemsXicTraceItems")
  # table of XIC traces
  xic_traces <- dplyr::tbl(msa$db_connection, "XicTraceItems")

  # assemble ID keys only for the compounds we're interested in
  keys <- dplyr::filter(ucids, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
  keys <- dplyr::left_join(keys, iids,
                           by = c("UnknownCompoundInstanceItemsWorkflowID",
                                  "UnknownCompoundInstanceItemsID"))
  keys <- dplyr::left_join(keys, xic_ids,
                           by = c("UnknownCompoundIonInstanceItemsWorkflowID",
                                  "UnknownCompoundIonInstanceItemsID"))

  # now add in the actual data
  # add features in each file
  xics <- dplyr::left_join(keys, unk_comp_items,
                           by = c("UnknownCompoundInstanceItemsWorkflowID" = "WorkflowID",
                                  "UnknownCompoundInstanceItemsID" = "ID"))
  # add ions associated with each feature
  xics <- dplyr::left_join(xics, unk_comp_ion_items,
                           by = c("UnknownCompoundIonInstanceItemsWorkflowID" = "WorkflowID",
                                  "UnknownCompoundIonInstanceItemsID" = "ID",
                                  "FileID", "StudyFileID"),
                           suffix = c(".compound", ".ion"))
  # add XIC traces for each ion
  xics <- dplyr::left_join(xics, xic_traces,
                           by = c("XicTraceItemsWorkflowID" = "WorkflowID",
                                  "XicTraceItemsID" = "ID", "FileID"),
                           suffix = c(".compound", ".xic"))
  # get local copy of the data
  xics <- dplyr::collect(xics)

  return(xics)
}
