#' Retrieve extracted ion chromatograms (XICs) associated with compounds
#'
#' Queries a mass spec alignment database to get all the extracted ion chromatograms
#' associated with some or all compounds. This includes all identified adducts
#' and isotopologues.
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get XICs for, or NULL to get all XICs
#'
#' @return A tibble with the XIC data
#'
#' @importFrom rlang .data
#' @export
get_compound_xics <- function(msa, ids = NULL) {
  if(is.null(ids)) {
    ids <- msa$unknown_compound_items$ID
  }

  # have to go to the unconsolidated unknown compounds table to get all XICs
  # get the unconsolidated compound IDs associated with all compounds
  ucids <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsUnknownCompoundInstanceItems")
  # get the table of unconsolidated unknown compounds
  unk_comp_items <- dplyr::tbl(msa$db_connection, "UnknownCompoundInstanceItems")
  # get the ions associated with unknown compounds
  iids <- dplyr::tbl(msa$db_connection, "UnknownCompoundInstanceItemsUnknownCompoundIonInstanceItems")
  # get the chromatographic peak IDs associated with unknown compounds
  xic_ids <- dplyr::tbl(msa$db_connection, "UnknownCompoundIonInstanceItemsXicTraceItems")
  # get the table of chromatographic peaks
  xic_traces <- dplyr::tbl(msa$db_connection, "XicTraceItems")
  # get only the peaks associated with the compounds we're interested in
  xics <- dplyr::filter(ucids, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
  # merge tables
  xics <- dplyr::left_join(xics, unk_comp_items,
                           by = c("UnknownCompoundInstanceItemsWorkflowID" = "WorkflowID",
                                  "UnknownCompoundInstanceItemsID" = "ID"))
  xics <- dplyr::left_join(xics, iids,
                           by = c("UnknownCompoundInstanceItemsWorkflowID",
                                  "UnknownCompoundInstanceItemsID"))
  xics <- dplyr::left_join(xics, xic_ids,
                           by = c("UnknownCompoundIonInstanceItemsWorkflowID",
                                  "UnknownCompoundIonInstanceItemsID"))
  # some chromatographic peaks aren't associated with a compound
  xics <- dplyr::left_join(xics, xic_traces,
                           by = c("XicTraceItemsWorkflowID" = "WorkflowID",
                                  "XicTraceItemsID" = "ID"),
                           suffix = c(".compound", ".xic"))
  # get local copy of the data
  xics <- dplyr::collect(xics)

  return(xics)
}
