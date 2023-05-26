#' Import the data in a mass spec data alignment
#'
#' Certain mass spec data alignment results are stored as SQLite databases. This
#' function imports relevant result tables for further analysis.
#'
#' @param f Path to file containing an alignment project in SQLite format.
#'
#' @return A list containing the alignment data tables.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
import_ms_alignment <- function(f) {
  # make a temp directory for unzipping blobs
  # tmp_zip_dir <- tempfile(pattern = "MSA_dir", tmpdir = tempdir())
  # dir.create(tmp_zip_dir)

  # connect to the alignment database
  con <- DBI::dbConnect(RSQLite::SQLite(), f)

  # which tables are available?
  all_tables <- DBI::dbListTables(con)
  # TODO: import more tables?
  # TODO: set up separate methods for importing certain table(s)?

  # what are the Compound Table headers?
  # DBI::dbListFields(con, "MassSpectrumItems")

  # read mass spectrum items
  spectrum_items <- DBI::dbReadTable(con, "MassSpectrumItems")
  # this table has all the spectra associated with a compound
  # actual centroided spectra are stored as blobs -- use extract_spectral_blob
  # to get individual spectra

  # read unknown compound items
  consolidated_unk_comp_items <- DBI::dbReadTable(con, "ConsolidatedUnknownCompoundItems")
  # this table is the main compounds table seen in CD
  # some of the info flags are encoded as binary blobs

  # read mzCloud search result items
  consolidated_unk_comp_items_mzcloud_hits <- DBI::dbReadTable(con, "ConsolidatedUnknownCompoundItemsMzCloudSearchResultItems")
  # this is the sub-table for mzCloud hits for each compound that has any
  # ID values point to associated compounds, library spectra, etc.

  # read library spectrum items
  library_spectrum_items <- DBI::dbReadTable(con, "LibrarySpectrumItems")
  # contains centroided spectra associated with mzCloud matches
  # same blob format as compound spectra, so can extract the same way

  # disconnect cleanly
  DBI::dbDisconnect(con)

  return(list(AvailableTableNames = all_tables,
              UnknownCompounds = consolidated_unk_comp_items,
              UnknownCompoundSpectra = spectrum_items,
              mzCloudHits = consolidated_unk_comp_items_mzcloud_hits,
              mzCloudHitSpectra = library_spectrum_items))

  # pull some spectrum blobs
  # blob_mass <- library_spectrum_items$Spectrum[[2]]
  # # these are zipped XML files in same format as feature spectra
  #
  # writeBin(blob_mass, "test.bin")
  # df <- readBin(blob_mass, numeric(), n = length(blob_mass)/8)
  #
  # # parse the spectral blobs
  # spec_masses <- lapply(spec_table$blobMass, FUN = function(b) {
  #   return(readBin(b, numeric(), n = length(b)/8))
  # })
  # spec_intensities <- lapply(spec_table$blobIntensity, FUN = function(b) {
  #   return(readBin(b, numeric(), n = length(b)/8))
  # })
  # spec_accuracies <- lapply(spec_table$blobAccuracy, FUN = function(b) {
  #   # some have missing accuracy values
  #   if(length(b) == 0) {
  #     return(c(0))
  #   }
  #   return(readBin(b, numeric(), n = length(b)/8))
  # })
  # spec_noises <- lapply(spec_table$blobNoises, FUN = function(b) {
  #   # some have missing noise values
  #   if(length(b) == 0) {
  #     return(c(0))
  #   }
  #   return(readBin(b, numeric(), n = length(b)/8))
  # })
  # not gonna bother with the top 10 peaks blob

  # # have to use a wrapper to make .data work right...
  # fn_rd <- function(mass, intensity, accuracy, noise) {
  #   return(tibble::tibble(Mass = mass, Intensity = intensity,
  #                 Accuracy = accuracy, Noise = noise))
  # }
  #
  # spec_data <- tibble::tibble(SpectrumId = spec_table$SpectrumId,
  #                             Mass = spec_masses,
  #                             Intensity = spec_intensities,
  #                             Accuracy = spec_accuracies,
  #                             Noise = spec_noises) %>%
  #   dplyr::rowwise() %>%
  #   dplyr::mutate(Spectrum = list(fn_rd(.data$Mass, .data$Intensity,
  #                                       .data$Accuracy, .data$Noise))) %>%
  #   dplyr::select(.data$SpectrumId, .data$Spectrum) %>%
  #   dplyr::ungroup()
  #
  # # now merge spectrum tibble column back into original table
  # df <- dplyr::full_join(x = spec_table, y = spec_data, by = "SpectrumId") %>%
  #   dplyr::select(-dplyr::starts_with("blob"))
  #
  # # disconnect cleanly
  # DBI::dbDisconnect(con)
  #
  # return(list(HeaderTable = tibble::as_tibble(hdr_table),
  #             MaintenanceTable = tibble::as_tibble(maint_table),
  #             CompoundTable = tibble::as_tibble(compound_table),
  #             SpectrumTable = df))

  # df <- spec_data %>%
  #   tidyr::unnest_longer(col = c(Mass, Intensity, Accuracy, Noise))
  # # read the first spectrum from blobs
  # blob_mass <- spec_table$blobMass[[1]]
  # spec_mass <- readBin(blob_mass, numeric(), n = length(blob_mass)/8)
  # blob_intensity <- spec_table$blobIntensity[[1]]
  # spec_intensity <- readBin(blob_intensity, numeric(), n = length(blob_intensity)/8)
  # blob_accuracy <- spec_table$blobAccuracy[[1]]
  # spec_accuracy <- readBin(blob_accuracy, numeric(), n = length(blob_accuracy)/8)
  # blob_flags <- spec_table$blobFlags[[1]]
  # spec_flags <- readBin(blob_flags, logical(), n = length(blob_flags)/1,
  #                       size = 1)
  # blob_top_peaks <- spec_table$blobTopPeaks[[1]]
  # spec_top_peaks <- readBin(blob_top_peaks, numeric(), n = length(blob_top_peaks)/8)
  # blob_noises <- spec_table$blobNoises[[1]] # blub blargh glurb
  # spec_noises <- readBin(blob_noises, numeric(), n = length(blob_noises)/8)
  #
  # # pull some spectrum blobs
  # blob_mass <- spec_table$blobNoises[[1]]
  # writeBin(blob_mass, "test.bin")
  # df <- readBin(blob_mass, numeric(), n = length(blob_mass)/8)
}
