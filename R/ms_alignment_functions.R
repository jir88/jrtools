#' Retrieve predicted elemental compositions
#'
#' Queries a mass spec alignment database to get the predicted elemental
#' compositions for some or all compounds.
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get compositions for, or NULL to get all predictions
#'
#' @return A tibble with the predicted composition data
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
get_elemental_compositions <- function(msa, ids = NULL) {
  if(is.null(ids)) {
    # get the composition IDs associated with all compounds
    df <- DBI::dbReadTable(msa$db_connection, "ConsolidatedUnknownCompoundItemsPredictedCompositionItem")
    # read predicted elemental compositions of all items
    predicted_compositions <- DBI::dbReadTable(msa$db_connection, "PredictedCompositionItem")
    # merge tables
    predicted_compositions <- dplyr::full_join(x = df, y = predicted_compositions,
                                               by = c("PredictedCompositionItemWorkflowID" = "WorkflowID",
                                                      "PredictedCompositionItemID" = "ID"))
  } else {
    # get the composition IDs associated with compounds of interest
    cid_tbl <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsPredictedCompositionItem")
    query_cids <- dplyr::filter(cid_tbl, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
    # read predicted elemental compositions of those items
    pc_tbl <- dplyr::tbl(msa$db_connection, "PredictedCompositionItem")

    predicted_compositions <- dplyr::left_join(x = query_cids, y = pc_tbl,
                                        by = c("PredictedCompositionItemWorkflowID" = "WorkflowID",
                                               "PredictedCompositionItemID" = "ID"))
    # get local copy of the data
    predicted_compositions <- dplyr::collect(predicted_compositions)
  }

  return(predicted_compositions)
}

#' Retrieve spectra associated with compounds
#'
#' Queries a mass spec alignment database to get the mass spectra associated with
#' some or all compounds. The actual spectra are stored as XML blobs which can
#' be extracted using extract_spectral_blob. Note that this method does NOT return
#' all the spectra in the alignment. Some spectra may not be associated with a
#' particular compound.
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get spectra for, or NULL to get all spectra
#'
#' @return A tibble with the spectra
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
get_compound_spectra <- function(msa, ids = NULL) {
  if(is.null(ids)) {
    # get the best hit ion IDs associated with all compounds
    bhi_ids <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsBestHitIonInstanceItems")
    # get the spectrum IDs associated with the best hit ions
    bhi_spectrum_ids <- dplyr::tbl(msa$db_connection, "BestHitIonInstanceItemsMassSpectrumItems")
    df <- dplyr::left_join(x = bhi_ids, y = bhi_spectrum_ids,
                           by = c("BestHitIonInstanceItemsWorkflowID", "BestHitIonInstanceItemsID"))
    # read predicted elemental compositions of all items
    all_spectra <- dplyr::tbl(msa$db_connection, "MassSpectrumItems")
    # merge tables
    bhi_spectra <- dplyr::left_join(x = df, y = all_spectra,
                                               by = c("MassSpectrumItemsWorkflowID" = "WorkflowID",
                                                      "MassSpectrumItemsID" = "ID"))
    # pull down the results
    bhi_spectra <- dplyr::collect(bhi_spectra)
  } else {
    # get the best hit ion IDs associated with compounds of interest
    bhi_ids <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsBestHitIonInstanceItems")
    query_bhi_ids <- dplyr::filter(bhi_ids, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
    # get the spectrum IDs associated with the best hit ions
    bhi_spectrum_ids <- dplyr::tbl(msa$db_connection, "BestHitIonInstanceItemsMassSpectrumItems")
    df <- dplyr::left_join(x = query_bhi_ids, y = bhi_spectrum_ids,
                           by = c("BestHitIonInstanceItemsWorkflowID", "BestHitIonInstanceItemsID"))
    # get spectra associated with best hit ions for those items
    all_spectra <- dplyr::tbl(msa$db_connection, "MassSpectrumItems")

    query_bhi_spectra <- dplyr::left_join(x = df, y = all_spectra,
                                        by = c("MassSpectrumItemsWorkflowID" = "WorkflowID",
                                               "MassSpectrumItemsID" = "ID"))
    # get local copy of the data
    bhi_spectra <- dplyr::collect(query_bhi_spectra)
  }

  return(bhi_spectra)
}
