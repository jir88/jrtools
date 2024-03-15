#' Extract spectral data from Compound Discoverer alignment
#'
#' Compound Discoverer cdResult alignment databases store feature spectra as SQLite
#' blobs containing a compressed XML file. This function extracts the stored
#' data from such blobs.
#'
#' Note: specifying a temp directory will speed up extraction of multiple blobs.
#'
#' @note
#' Compound Discoverer software is produced by Thermo Fisher Scientific. This
#' package is not affiliated with Thermo Fisher Scientific in any way. For an
#' official Python interface to Compound Discoverer alignment files, see
#' \url{https://github.com/thermofisherlsms/pyeds}
#'
#' @param blb A raw vector containing an SQLite blob of spectral data
#' @param meta Should spectrum metadata be returned, or just the tibble of
#'   mass spectral data? Default is FALSE.
#' @param zip_dir Directory where XML files should be unzipped
#'
#' @return A tibble containing the mass spectrum. If meta is TRUE, returns a list
#'   with the spectrum metadata and the mass spectrum tibble.
#'
#' @export
extract_spectral_blob <- function(blb, meta = FALSE, zip_dir = tempdir()) {
  # spectrum blobs are zipped XML files
  # format looks custom, but is human readable
  # contains scan info and centroided spectra

  # read an XML spectrum
  xml_dir <- read_zip_blob(zb = blb, blob_path = zip_dir)
  # should only be the one
  xml_file <- list.files(xml_dir, full.names = TRUE)
  # load the XML document
  spec_xml <- xml2::read_xml(xml_file)
  # turn it into a list
  # actually, puts lots of data into R attributes, which would be awful to extract
  # spec_list <- xml2::as_list(spec_xml)
  # grab the centroided peak data
  pks_node <- xml2::xml_child(spec_xml, search = "PeakCentroids")
  pk_centroids <- xml2::xml_children(pks_node)
  pk_centroid_data <- xml2::xml_attrs(pk_centroids)
  # shove it into a matrix
  pk_centroid_data <- do.call(rbind, pk_centroid_data)
  # convert to numeric
  class(pk_centroid_data) <- "numeric"
  # fix up column names without assuming consistent order
  cn <- colnames(pk_centroid_data)
  cn <- c("mz", "intensity", "Z", "resolution", "signalNoiseRatio")[match(cn, c("X", "Y", "Z", "R", "SN"))]
  colnames(pk_centroid_data) <- cn
  # stuff the data into a tidy tibble
  pk_centroid_data <- tibble::as_tibble(pk_centroid_data, .name_repair = "minimal")

  if(meta) {
    # currently we just rip the metadata into nested lists
    hdr <- xml2::as_list(xml2::xml_child(spec_xml, search = "Header"))
    scan_event <- xml2::as_list(xml2::xml_child(spec_xml, search = "ScanEvent"))
    prec_info <- xml2::as_list(xml2::xml_child(spec_xml, search = "PrecursorInfo"))

    return(list("Spectrum" = pk_centroid_data,
                "Header" = hdr,
                "ScanEvent" = scan_event,
                "PrecursorInfo" = prec_info))
  } else {
    return(pk_centroid_data)
  }
}
