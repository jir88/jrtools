#' Extract fitted peak model parameters from Compound Discoverer alignment
#'
#' Compound Discoverer cdResult alignment databases store fitted peak models as
#' SQLite blobs containing a compressed XML file. This function extracts the
#' stored data from such blobs.
#'
#' Note: specifying a temp directory will speed up extraction of multiple blobs.
#'
#' @note
#' Compound Discoverer software is produced by Thermo Fisher Scientific. This
#' package is not affiliated with Thermo Fisher Scientific in any way. For an
#' official Python interface to Compound Discoverer alignment files, see
#' \url{https://github.com/thermofisherlsms/pyeds}
#'
#' @param blb A raw vector containing an SQLite blob of a peak model
#' @param zip_dir Directory where XML files should be unzipped
#'
#' @return The XML structure describing the peak model.
#'
#' @export
extract_peak_model_blob <- function(blb, zip_dir = tempdir()) {
  # peak blob is a zipped XML file
  # containing a base64 encoded string
  # containing a zip file
  # containing another XML file with the actual model in it

  # read the outer zipped XML blob
  xml_dir <- read_zip_blob(zb = blb, blob_path = zip_dir)
  # should only be the one
  xml_file <- list.files(xml_dir, full.names = TRUE)
  # load the XML document
  spec_xml <- xml2::read_xml(xml_file)
  # model data is base64 encoded...
  pk_xml <- xml2::xml_child(spec_xml, search = "Data")
  pk_xml <- xml2::xml_text(pk_xml)
  pk_xml <- base64enc::base64decode(pk_xml)
  # the data is another zipped blob
  pk_xml <- read_zip_blob(zb = pk_xml, blob_path = zip_dir)
  # with another XML file in it
  pk_xml <- xml2::read_xml(list.files(pk_xml, full.names = TRUE))
  # mercifully this is the last layer
  return(pk_xml)
}
