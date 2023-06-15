#' Extract spectral data from an SQLite blob
#'
#' Certain mass spec data alignment formats store feature spectra as SQLite
#' blobs containing a compressed XML file. This function extracts the stored
#' data from such blobs.
#'
#' The XML file contains a wealth of metadata about the spectrum, which is not
#' currently extracted by this function.
#'
#' @param blb A raw vector containing an SQLite blob of spectral data
#' @param zip_dir Directory where XML files should be unzipped
#'
#' @return A list containing the library data tables.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
extract_spectral_blob <- function(blb, zip_dir = tempdir()) {
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
  # stuff the data into a tidy tibble
  pk_centroid_data <- dplyr::bind_rows(pk_centroid_data)
  # convert text to numbers
  pk_centroid_data <- dplyr::mutate(pk_centroid_data,
                                    X = as.numeric(.data$X),
                                    Y = as.numeric(.data$Y),
                                    Z = as.numeric(.data$Z),
                                    R = as.numeric(.data$R),
                                    SN = as.numeric(.data$SN))
  pk_centroid_data <- dplyr::rename(pk_centroid_data,
                                    mz = .data$X,
                                    intensity = .data$Y,
                                    #Z = .data$Z,
                                    resolution = .data$R,
                                    signalNoiseRatio = .data$SN)

  return(pk_centroid_data)
}
