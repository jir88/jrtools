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

#' Retrieve cloud spectral library matches
#'
#' Queries a mass spec alignment database to get the cloud-based spectral libary
#' matches for some or all compounds.
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get matches for, or NULL to get all matches
#'
#' @return A tibble with the cloud library matches
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
get_spectral_cloud_matches <- function(msa, ids = NULL) {
  if(is.null(ids)) {
    # get the mzCloud hit IDs associated with all compounds
    df <- DBI::dbReadTable(msa$db_connection, "ConsolidatedUnknownCompoundItemsMzCloudHitItems")
    # read predicted elemental compositions of all items
    mzcloud_hits <- DBI::dbReadTable(msa$db_connection, "MzCloudHitItems")
    # merge tables
    mzcloud_hits <- dplyr::full_join(x = df, y = mzcloud_hits,
                                               by = c("MzCloudHitItemsID" = "ID"))
  } else {
    # get the composition IDs associated with compounds of interest
    mzcloud_id_tbl <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsMzCloudHitItems")
    query_mzcloud_ids <- dplyr::filter(mzcloud_id_tbl, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
    # read predicted elemental compositions of those items
    mzcloud_hit_tbl <- dplyr::tbl(msa$db_connection, "MzCloudHitItems")

    mzcloud_hits <- dplyr::left_join(x = query_mzcloud_ids, y = mzcloud_hit_tbl,
                                               by = c("MzCloudHitItemsID" = "ID"))
    # get local copy of the data
    mzcloud_hits <- dplyr::collect(mzcloud_hits)
  }

  return(mzcloud_hits)
}

#' Extract predicted isotopologue pattern data from an SQLite blob
#'
#' Certain mass spec data alignment formats store predicted isotopologue spectra
#' as SQLite blobs containing a compressed XML file. This function extracts the
#' stored data from such blobs.
#'
#' @param blb A raw vector containing an SQLite blob of isotopologue data
#' @param zip_dir Directory where XML files should be unzipped
#'
#' @return A tibble containing the isotopologue pattern
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
extract_isotope_pattern_blob <- function(blb, zip_dir = tempdir()) {
  # spectrum blobs are zipped XML files
  # format looks custom, but is human readable
  # contains scan info and centroided spectra

  # read an XML spectrum
  xml_dir <- read_zip_blob(zb = blb, blob_path = zip_dir)
  # should only be the one
  xml_file <- list.files(xml_dir, full.names = TRUE)
  # load the XML document
  spec_xml <- xml2::read_xml(xml_file)

  # pull out the pattern node and shove it into a tidy tibble
  pattern_peak_data <- xml2::xml_child(spec_xml, search = "PatternPeaks") %>%
    # get peak nodes
    xml2::xml_children() %>%
    # pull attribute data
    xml2::xml_attrs() %>%
    # mash into a tibble
    dplyr::bind_rows() %>%
    # convert text to numbers
    dplyr::mutate(mz = as.numeric(.data$X),
                  intensity = as.numeric(.data$Y),
                  charge = as.numeric(.data$Z),
                  resolution = as.numeric(.data$R),
                  signalNoiseRatio = as.numeric(.data$SN),
                  .keep = "none")

  return(pattern_peak_data)
}

# unzip a blob from the database
# zb: the blob
# path: directory to put unzipped blob in
# returns: path to directory with blob's contents
read_zip_blob <- function(zb, blob_path) {
  # write blob to file
  tmp_blob_file <- tempfile(tmpdir = blob_path, fileext = ".zip")
  writeBin(zb, tmp_blob_file)
  # create subdir for this blob's contents
  blob_content_dir <- tempfile(tmpdir = blob_path)
  # extract to directory
  zip::unzip(tmp_blob_file, exdir = blob_content_dir)
  # return directory path so we can do stuff with the contents
  return(blob_content_dir)
}
