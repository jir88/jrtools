#' Retrieve extracted ion chromatograms (XICs) associated with gap filling
#'
#' Queries a Compound Discoverer mass spec alignment database to get all the
#' extracted ion chromatograms associated with gap filling for some or all
#' compounds. This includes all identified adducts and isotopologues. Note that
#' this only gets XICs where a peak NOT was detected and gap filling was
#' performed. For XICs of detected peaks, see \code{\link{get_compound_xics}}.
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
get_gap_xics <- function(msa, ids = NULL) {
  if(is.null(ids)) {
    ids <- msa$unknown_compound_items$ID
  }

  # key of missing ions associated with all compounds
  uc_mi_ids <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsMissingCompoundIonInstanceItems")
  # table of missing ions
  missing_ions <- dplyr::tbl(msa$db_connection, "MissingCompoundIonInstanceItems")
  # key of XIC traces associated with missing ions
  xic_ids <- dplyr::tbl(msa$db_connection, "MissingCompoundIonInstanceItemsXicTraceItems")
  # table of XIC traces
  xic_traces <- dplyr::tbl(msa$db_connection, "XicTraceItems")

  # assemble ID keys only for the compounds we're interested in
  keys <- dplyr::filter(uc_mi_ids, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
  keys <- dplyr::left_join(keys, xic_ids,
                           by = c("MissingCompoundIonInstanceItemsWorkflowID",
                                  "MissingCompoundIonInstanceItemsID"))

  # now add in the actual data
  # add ions in each file where gap filling was required
  xics <- dplyr::left_join(keys, missing_ions,
                           by = c("MissingCompoundIonInstanceItemsWorkflowID" = "WorkflowID",
                                  "MissingCompoundIonInstanceItemsID" = "ID"))
  # add XIC traces for each ion
  xics <- dplyr::left_join(xics, xic_traces,
                           by = c("XicTraceItemsWorkflowID" = "WorkflowID",
                                  "XicTraceItemsID" = "ID", "FileID"),
                           suffix = c(".ion", ".xic"))
  # get local copy of the data
  xics <- dplyr::collect(xics)

  # turn fill status column into a factor
  fs_levels <- jrtools::get_cdresult_enum(msa = msa,
                                          table = "MissingCompoundIonInstanceItems",
                                          enum = "FillStatus")
  xics$FillStatus <- factor(xics$FillStatus, levels = fs_levels$Value,
                            labels = fs_levels$DisplayName)

  return(xics)
}
