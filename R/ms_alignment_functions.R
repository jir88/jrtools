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

#' Retrieve chromatogram peaks associated with compounds
#'
#' Queries a mass spec alignment database to get all the chromatographic peaks
#' associated with some or all compounds. This includes all identified adducts.
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get peaks for, or NULL to get all peaks
#'
#' @return A tibble with the chromatographic peaks
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
get_chromatogram_peaks <- function(msa, ids = NULL) {
  if(is.null(ids)) {
    # have to go to the unconsolidated unknown compounds table to get all peaks
    # get the unconsolidated compound IDs associated with all compounds
    ucids <- DBI::dbReadTable(msa$db_connection, "ConsolidatedUnknownCompoundItemsUnknownCompoundInstanceItems")
    # get the table of unconsolidated unknown compounds
    unk_comp_items <- DBI::dbReadTable(msa$db_connection, "UnknownCompoundInstanceItems")
    # get the ions associated with unknown compounds
    iids <- DBI::dbReadTable(msa$db_connection, "UnknownCompoundInstanceItemsUnknownCompoundIonInstanceItems")
    # get the chromatographic peak IDs associated with unknown compounds
    cpids <- DBI::dbReadTable(msa$db_connection, "UnknownCompoundIonInstanceItemsChromatogramPeakItems")
    # get the table of chromatographic peaks
    chrom_peaks <- DBI::dbReadTable(msa$db_connection, "ChromatogramPeakItems")
    # merge tables
    all_chrom_peaks <- dplyr::full_join(ucids, unk_comp_items,
                                        by = c("UnknownCompoundInstanceItemsWorkflowID" = "WorkflowID",
                                               "UnknownCompoundInstanceItemsID" = "ID"))
    all_chrom_peaks <- dplyr::full_join(all_chrom_peaks, iids,
                                        by = c("UnknownCompoundInstanceItemsWorkflowID",
                                               "UnknownCompoundInstanceItemsID"))
    all_chrom_peaks <- dplyr::full_join(all_chrom_peaks, cpids,
                                        by = c("UnknownCompoundIonInstanceItemsWorkflowID",
                                               "UnknownCompoundIonInstanceItemsID"))
    # some chromatographic peaks aren't associated with a compound
    all_chrom_peaks <- dplyr::left_join(all_chrom_peaks, chrom_peaks,
                                        by = c("ChromatogramPeakItemsWorkflowID" = "WorkflowID",
                                               "ChromatogramPeakItemsID" = "ID"),
                                        suffix = c(".compound", ".peak"))
  } else {
    # have to go to the unconsolidated unknown compounds table to get all peaks
    # get the unconsolidated compound IDs associated with all compounds
    ucids <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsUnknownCompoundInstanceItems")
    # get the table of unconsolidated unknown compounds
    unk_comp_items <- dplyr::tbl(msa$db_connection, "UnknownCompoundInstanceItems")
    # get the ions associated with unknown compounds
    iids <- dplyr::tbl(msa$db_connection, "UnknownCompoundInstanceItemsUnknownCompoundIonInstanceItems")
    # get the chromatographic peak IDs associated with unknown compounds
    cpids <- dplyr::tbl(msa$db_connection, "UnknownCompoundIonInstanceItemsChromatogramPeakItems")
    # get the table of chromatographic peaks
    chrom_peaks <- dplyr::tbl(msa$db_connection, "ChromatogramPeakItems")
    # get only the peaks associated with the compounds we're interested in
    all_chrom_peaks <- dplyr::filter(ucids, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
    # merge tables
    all_chrom_peaks <- dplyr::left_join(all_chrom_peaks, unk_comp_items,
                                        by = c("UnknownCompoundInstanceItemsWorkflowID" = "WorkflowID",
                                               "UnknownCompoundInstanceItemsID" = "ID"))
    all_chrom_peaks <- dplyr::left_join(all_chrom_peaks, iids,
                                        by = c("UnknownCompoundInstanceItemsWorkflowID",
                                               "UnknownCompoundInstanceItemsID"))
    all_chrom_peaks <- dplyr::left_join(all_chrom_peaks, cpids,
                                        by = c("UnknownCompoundIonInstanceItemsWorkflowID",
                                               "UnknownCompoundIonInstanceItemsID"))
    # some chromatographic peaks aren't associated with a compound
    all_chrom_peaks <- dplyr::left_join(all_chrom_peaks, chrom_peaks,
                                        by = c("ChromatogramPeakItemsWorkflowID" = "WorkflowID",
                                               "ChromatogramPeakItemsID" = "ID"),
                                        suffix = c(".compound", ".peak"))
    # get local copy of the data
    all_chrom_peaks <- dplyr::collect(all_chrom_peaks)
  }

  return(all_chrom_peaks)
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
    # get mzCloud search results associated with all compounds
    unk_comp_search_res <- DBI::dbReadTable(msa$db_connection,
                                            "ConsolidatedUnknownCompoundItemsMzCloudSearchResultItems")
    # get all mzCloud hits -- less detail, but points to relevant spectra
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
#' @importFrom magrittr %>%
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

#' Extract peak areas stored as blobs in a mass spec alignment
#'
#' Certain mass spec data alignment formats store peak areas as SQLite blobs
#' containing the data in binary form. This function extracts the
#' stored data from such blobs. The returned areas are for the reference adduct
#' ion for each compound. Other areas must be retrieved manually from the relevant
#' tables.
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get areas for, or NULL to get all peak areas
#'
#' @return A tibble containing the peak areas and another containing boolean
#'  flag reporting whether the peak exists or not
#'
#' @importFrom rlang .data
#' @export
extract_peak_areas <- function(msa, ids = NULL) {
  # blobs contain alternating 64-bit floats and a binary byte specifying whether
  # the area was calculated. Presumably this is to differentiate between true
  # zeroes and NA values?

  # use all IDs
  if(is.null(ids)) {
    ids <- msa$unknown_compound_items$ID
  }
  # usually index and ID are the same, but we can't be sure
  idx <- match(ids, msa$unknown_compound_items$ID)

  area_data <- lapply(X = msa$unknown_compound_items$Area[idx],
                      FUN = function(blb) {
                        area_values <- blb[-seq.int(from = 9, to = length(blb), by = 9)]
                        area_values <- readBin(area_values, what = numeric(), n = length(area_values)/8)

                        area_flags <- blb[seq.int(from = 9, to = length(blb), by = 9)]
                        area_flags <- readBin(area_flags, what = logical(),
                                              size = 1, n = length(area_flags))
                        return(list(area = area_values,
                                    flag = area_flags))
                      })
  area_data <- purrr::transpose(area_data)

  areas <- do.call(rbind, area_data$area)
  rownames(areas) <- ids
  areas <- tibble::as_tibble(t(areas), .name_repair = "minimal")
  areas <- dplyr::mutate(areas, StudyFileID = msa$input_files$StudyFileID,
                         .before = 1)

  flags <- do.call(rbind, area_data$flag)
  rownames(flags) <- ids
  flags <- tibble::as_tibble(t(flags), .name_repair = "minimal")
  flags <- dplyr::mutate(flags, StudyFileID = msa$input_files$StudyFileID,
                         .before = 1)

  return(list("areas" = areas, "flags" = flags))
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
