#' Get descriptive info for tables and columns in the alignment
#'
#' Queries a mass spec alignment database to get annotations for all top-level
#' data tables and their columns. Also gets level annotations associated with
#' categorical ("Enum") column types.
#'
#' @param msa An ms_alignment object to query
#'
#' @return A list of two tibbles, one with the table annotation data and the other
#'  with the categorical data types info
#'
#' @importFrom rlang .data
#' @export
get_data_table_info <- function(msa) {
  # get table listing all top-level data tables
  table_info <- dplyr::tbl(msa$db_connection, "DataTypes")
  # get the table of column metadata
  column_info <- dplyr::tbl(msa$db_connection, "DataTypesColumns")
  # fix mixed column to be all strings
  column_info <- dplyr::mutate(column_info, DefaultValue = as.character(.data$DefaultValue))
  # combine tables
  table_info <- dplyr::full_join(table_info, column_info,
                                 by = "DataTypeID",
                                 suffix = c(".table", ".column"))
  table_info <- dplyr::collect(table_info)

  # get table listing all enum types
  enum_info <- dplyr::tbl(msa$db_connection, "EnumDataTypes")
  # get table listing all enum type levels
  enum_level_info <- dplyr::tbl(msa$db_connection, "EnumDataTypeValues")
  # combine tables
  enum_info <- dplyr::full_join(enum_info, enum_level_info,
                                 by = "EnumID",
                                 suffix = c(".enum", ".values"))
  enum_info <- dplyr::collect(enum_info)

  return(list(table_info = table_info, enum_info = enum_info))
}

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

#' Get all ions associated with particular compounds
#'
#' Queries a mass spec alignment database to get all ions associated with each
#' compound on a per file basis.
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get ions for, or NULL to get all best ions
#'
#' @return A tibble with the ion data
#'
#' @importFrom rlang .data
#' @export
get_compound_ions <- function(msa, ids = NULL) {
  if(is.null(ids)) {
    ids <- msa$unknown_compound_items$ID
  }

  # have to go to the unconsolidated unknown compounds table to get all ions
  # get the unconsolidated compound IDs associated with all compounds
  ucids <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsUnknownCompoundInstanceItems")
  # get the table of unconsolidated unknown compounds
  unk_comp_items <- dplyr::tbl(msa$db_connection, "UnknownCompoundInstanceItems")
  # get the ion IDs associated with unknown compounds
  iids <- dplyr::tbl(msa$db_connection, "UnknownCompoundInstanceItemsUnknownCompoundIonInstanceItems")
  # get the actual ion table
  cpids <- dplyr::tbl(msa$db_connection, "UnknownCompoundIonInstanceItems")

  # get only the ions associated with the compounds we're interested in
  all_ions <- dplyr::filter(ucids, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
  # merge tables
  all_ions <- dplyr::left_join(all_ions, unk_comp_items,
                                      by = c("UnknownCompoundInstanceItemsWorkflowID" = "WorkflowID",
                                             "UnknownCompoundInstanceItemsID" = "ID"))
  all_ions <- dplyr::left_join(all_ions, iids,
                                      by = c("UnknownCompoundInstanceItemsWorkflowID",
                                             "UnknownCompoundInstanceItemsID"))
  all_ions <- dplyr::left_join(all_ions, cpids,
                                      by = c("UnknownCompoundIonInstanceItemsWorkflowID" = "WorkflowID",
                                             "UnknownCompoundIonInstanceItemsID" = "ID",
                                             "FileID", "StudyFileID", "IdentifyingNodeNumber"),
                               suffix = c(".compound", ".ion"))
  # get local copy of the data
  all_ions <- dplyr::collect(all_ions)

  return(all_ions)
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

#' Retrieve all spectra in alignment
#'
#' Queries a mass spec alignment database to get all the mass spectra in it. The
#' actual spectra are stored as XML blobs which can be extracted using
#' extract_spectral_blob.
#'
#' @param msa An ms_alignment object to query
#'
#' @return A tibble with the spectra
#'
#' @importFrom rlang .data
#' @export
get_all_spectra <- function(msa) {
  all_spectra <- dplyr::tbl(msa$db_connection, "MassSpectrumItems")
  return(dplyr::collect(all_spectra))
}

#' Retrieve cloud spectral library matches
#'
#' Queries a mass spec alignment database to get the cloud-based spectral library
#' matches for some or all compounds.
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get matches for, or NULL to get all matches
#'
#' @return A tibble with the cloud library matches, or NA if this alignment
#'  database does not contain such data. An empty tibble indicates that cloud
#'  library matches exist, but none of the requested compounds have such matches.
#'
#' @importFrom rlang .data
#' @export
get_spectral_cloud_matches <- function(msa, ids = NULL) {
  # check to see if there are any hits at all
  if(!DBI::dbExistsTable(msa$db_connection, "MzCloudSearchResultItems")) {
    simpleWarning("This alignment database does not have cloud library matches!")
    return(NA)
  }
  # if no IDs are specified, return hits for all consolidated compounds
  if(is.null(ids)) {
    ids <- msa$unknown_compound_items$ID
  }

  # old database format lacks the cloud "hits", so we start with "search results"
  # get the mzCloud search result IDs associated with compounds of interest
  mzcloud_id_tbl <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsMzCloudSearchResultItems")
  mzcloud_search_res_tbl <- dplyr::tbl(msa$db_connection, "MzCloudSearchResultItems")

  # get just the search results we're interested in
  cloud_matches <- dplyr::filter(mzcloud_id_tbl, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
  cloud_matches <- dplyr::left_join(x = cloud_matches, y = mzcloud_search_res_tbl,
                                    by = c("MzCloudSearchResultItemsID" = "ID"))

  # if this database has the "hits", we'll put them in too
  if(DBI::dbExistsTable(msa$db_connection, "MzCloudHitItems")) {
    # read the mzCloud hit info for all items
    mzcloud_hit_tbl <- dplyr::tbl(msa$db_connection, "MzCloudHitItems")
    # get the search to hit mapping
    mzcloud_hit_ids <- dplyr::tbl(msa$db_connection, "MzCloudHitItemsMzCloudSearchResultItems")
    # add hit info to search results
    cloud_matches <- dplyr::left_join(x = cloud_matches, y = mzcloud_hit_ids,
                                      by = c("MzCloudSearchResultItemsID"))
    cloud_matches <- dplyr::left_join(x = cloud_matches, y = mzcloud_hit_tbl,
                                      by = c("MzCloudHitItemsID" = "ID",
                                             "Name"),
                                      suffix = c(".result", ".hit"))
    # can also pull in associated spectral blobs
    library_spectra_tbl <- dplyr::tbl(msa$db_connection, "LibrarySpectrumItems")
    cloud_matches <- dplyr::left_join(x = cloud_matches, y = library_spectra_tbl,
                                      by = c("LibrarySpectrumId" = "ID"))
  }

  # get local copy of the data
  cloud_matches <- dplyr::collect(cloud_matches)

  return(cloud_matches)
}

#' Get ChemSpider matches associated with particular compounds
#'
#' Queries a mass spec alignment database to get any ChemSpider matches associated
#' with some or all compounds.
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get matches for, or NULL to get all matches
#'
#' @return A tibble with the ChemSpider matches
#'
#' @importFrom rlang .data
#' @export
get_chemspider_matches <- function(msa, ids = NULL) {
  if(is.null(ids)) {
    ids <- msa$unknown_compound_items$ID
  }

  # get the match IDs associated with all compounds
  ucids <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsChemSpiderResultItems")
  # get the table of all ChemSpider matches
  cs_results <- dplyr::tbl(msa$db_connection, "ChemSpiderResultItems")

  # get only the matches associated with the compounds we're interested in
  all_matches <- dplyr::filter(ucids, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
  # merge tables
  all_matches <- dplyr::left_join(all_matches, cs_results,
                               by = c("ChemSpiderResultItemsChemSpiderID" = "ChemSpiderID"))
  # get local copy of the data
  all_matches <- dplyr::collect(all_matches)

  return(all_matches)
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

#' Get best ions associated with compounds
#'
#' Queries a mass spec alignment database to get what the software considers to
#' be the "best representative compound data" for each compound.
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get best ions for, or NULL to get all best ions
#'
#' @details
#' The BestHitType column can be 0 (unknown), 1 (best MS1), 2 (best MS2), or
#' 4 (best deconvoluted MS).
#'
#' @return A tibble with the best ion data
#'
#' @importFrom rlang .data
#' @export
get_compound_best_ions <- function(msa, ids = NULL) {
  if(is.null(ids)) {
    ids <- msa$unknown_compound_items$ID
  }

  # get table matching CIDs to best hit ions
  ucids <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsBestHitIonInstanceItems")
  # get the table of best hit ions
  all_best_hit_ions <- dplyr::tbl(msa$db_connection, "BestHitIonInstanceItems")
  # get only the ions associated with the compounds we're interested in
  best_hit_ions <- dplyr::filter(ucids, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
  # merge tables
  best_hit_ions <- dplyr::left_join(best_hit_ions, all_best_hit_ions,
                           by = c("BestHitIonInstanceItemsWorkflowID" = "WorkflowID",
                                  "BestHitIonInstanceItemsID" = "ID"))
  # get local copy of the data
  best_hit_ions <- dplyr::collect(best_hit_ions)

  return(best_hit_ions)
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
#' @importFrom dplyr %>%
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
#' @param adj Should QC-adjusted peak areas be returned instead of raw peak
#'  areas? Defaults to FALSE.
#'
#' @return A tibble containing the peak areas and another containing boolean
#'  flag reporting whether the peak exists or not
#'
#' @importFrom rlang .data
#' @export
extract_peak_areas <- function(msa, ids = NULL, adj = FALSE) {
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

  # get the desired type of area blobs
  if(adj) {
    area_blobs <- msa$unknown_compound_items$NormArea[idx]
    # many compounds are not normalized, so we need to deal with those
    found_area_idx <- which(!sapply(area_blobs, is.null))
    area_blobs <- area_blobs[!sapply(area_blobs, is.null)]
  } else {
    area_blobs <- msa$unknown_compound_items$Area[idx]
    # excluded compounds don't get integrated, so we need to deal with those
    found_area_idx <- which(!sapply(area_blobs, is.null))
    area_blobs <- area_blobs[!sapply(area_blobs, is.null)]
  }
  area_data <- lapply(X = area_blobs,
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

  # if there are any IDs with un-integrated peaks
  if(length(found_area_idx) < length(idx)) {
    df <- matrix(NA_real_, nrow = length(idx), ncol = ncol(areas))
    df[found_area_idx, ] <- areas
    # assign expanded matrix
    areas <- df
  }

  rownames(areas) <- ids
  areas <- tibble::as_tibble(t(areas), .name_repair = "minimal")
  areas <- dplyr::mutate(areas, StudyFileID = sfids,
                         .before = 1)

  flags <- do.call(rbind, area_data$flag)

  # if there are any IDs with un-integrated peaks
  if(length(found_area_idx) < length(idx)) {
    # sensible value here is false, not NA
    df <- matrix(FALSE, nrow = length(idx), ncol = ncol(flags))
    df[found_area_idx, ] <- flags
    # assign expanded matrix
    flags <- df
  }

  rownames(flags) <- ids
  flags <- tibble::as_tibble(t(flags), .name_repair = "minimal")
  flags <- dplyr::mutate(flags, StudyFileID = sfids,
                         .before = 1)

  return(list("areas" = areas, "flags" = flags))
}

#' Get retention time corrections for scans in a mass spec alignment
#'
#' Certain mass spec data alignment formats store retention time correction data
#' as SQLite blobs containing the data in binary form. This function extracts the
#' stored data from such blobs. The returned retention times are presumably for
#' each scan in the raw data file.
#'
#' @param msa An ms_alignment object to query
#' @param file_id ID numbers for the file retention time corrections to return,
#'  or NULL to return all of them
#'
#' @return A tibble containing the original and corrected retention times for
#' each file in minutes.
#'
#' @importFrom rlang .data
#' @export
get_file_rt_corrections <- function(msa, file_id = NULL) {
  if(is.null(file_id)) {
    file_id <- msa$input_files$FileID
  }

  fa_cor_tbl <- dplyr::tbl(msa$db_connection, "FileAlignmentCorrectionItems")
  fa_cor_tbl <- dplyr::filter(fa_cor_tbl, .data$FileID %in% file_id)
  fa_cor_tbl <- dplyr::collect(fa_cor_tbl)

  # read retention time blobs
  rt_corrections <- lapply(X = seq.int(from = 1, to = nrow(fa_cor_tbl), by = 1), FUN = function(r) {
    # original retention times
    rt_blob <- fa_cor_tbl$OriginalRT[[r]]
    # grab length of blob
    blob_len <- readBin(rt_blob[1:4], n = 1, what = "integer", size = 4, endian = "little")
    # this is actually a blob of 64-bit doubles, NOT integers, idiot
    orig_rt <- readBin(rt_blob[-1:-4], n = blob_len, what = "double", size = 8,
                       endian = "little")
    # corrected retention times
    rt_blob <- fa_cor_tbl$CorrectedRT[[r]]
    # grab length of blob
    blob_len <- readBin(rt_blob[1:4], n = 1, what = "integer", size = 4, endian = "little")
    # this is actually a blob of 64-bit doubles, NOT integers, idiot
    corr_rt <- readBin(rt_blob[-1:-4], n = blob_len, what = "double", size = 8,
                       endian = "little")
    return(tibble::tibble(OriginalRT = orig_rt, CorrectedRT = corr_rt))
  })
  names(rt_corrections) <- fa_cor_tbl$StudyFileID
  rt_corrections <- dplyr::bind_rows(rt_corrections, .id = "StudyFileID")
  rt_corrections <- dplyr::full_join(x = dplyr::select(fa_cor_tbl, "WorkflowId":"StudyFileID"),
                                     y = rt_corrections, by = "StudyFileID")
  return(rt_corrections)
}

extract_rt_raster_trace <- function(msa, file_id = NULL) {
  if(is.null(file_id)) {
    file_id <- msa$input_files$FileID
  }

  rt_raster_tbl <- dplyr::tbl(msa$db_connection, "RetentionTimeRasterItem")
  rt_raster_tbl <- dplyr::filter(rt_raster_tbl, .data$FileID %in% file_id)
  rt_raster_tbl <- dplyr::collect(rt_raster_tbl)

  data_start <- 23
  traces <- lapply(X = rt_raster_tbl$Trace, FUN = function(compr_data) {
    # gunzip the data
    trace_data <- memDecompress(compr_data, type = "gzip")
    # how long are the chunks?
    chunk_len <- readBin(trace_data[18:21], what = "integer", n = 1, size = 4)
    # read the scan numbers (?)
    scan_index_end <- data_start + 4*chunk_len
    scan_index <- readBin(trace_data[data_start:(scan_index_end)],
                          what = "integer", n = chunk_len, size = 4)
    # next chunk
    chunk2_end <- scan_index_end + 2 + 4*chunk_len
    chunk2 <- readBin(trace_data[(scan_index_end + 1):(chunk2_end)],
                      what = "integer", n = chunk_len, size = 4)
    # next chunk, possibly alignment graph
    chunk3_end <- chunk2_end + 2 + 4*chunk_len
    chunk3 <- readBin(trace_data[(chunk2_end + 1):(chunk3_end)],
                      what = "integer", n = chunk_len, size = 4)
    return(dplyr::tibble(Chunk1 = scan_index,
                  Chunk2 = chunk2,
                  Chunk3 = chunk3))
  })
  names(traces) <- rt_raster_tbl$ID
  traces <- dplyr::bind_rows(traces, .id = "ID")
  traces <- dplyr::mutate(traces, ID = as.numeric(.data$ID))
  df <- dplyr::full_join(dplyr::select(rt_raster_tbl, -"Trace"), traces, by = "ID")

  return(df)
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


