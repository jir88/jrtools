% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_spectral_blob.R
\name{extract_spectral_blob}
\alias{extract_spectral_blob}
\title{Extract spectral data from Compound Discoverer alignment}
\usage{
extract_spectral_blob(blb, meta = FALSE, zip_dir = tempdir())
}
\arguments{
\item{blb}{A raw vector containing an SQLite blob of spectral data}

\item{meta}{Should spectrum metadata be returned, or just the tibble of
mass spectral data? Default is FALSE.}

\item{zip_dir}{Directory where XML files should be unzipped}
}
\value{
A tibble containing the mass spectrum. If meta is TRUE, returns a list
  with the spectrum metadata and the mass spectrum tibble.
}
\description{
Compound Discoverer cdResult alignment databases store feature spectra as SQLite
blobs containing a compressed XML file. This function extracts the stored
data from such blobs.
}
\details{
Note: specifying a temp directory will speed up extraction of multiple blobs.
}
\note{
Compound Discoverer software is produced by Thermo Fisher Scientific. This
package is not affiliated with Thermo Fisher Scientific in any way. For an
official Python interface to Compound Discoverer alignment files, see
\url{https://github.com/thermofisherlsms/pyeds}
}
