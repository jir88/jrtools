% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idx_functions.R
\name{export_idx_spectrum}
\alias{export_idx_spectrum}
\title{Export a Spectrum2 object in FreeStyle CSV format}
\usage{
export_idx_spectrum(spectrum, file)
}
\arguments{
\item{spectrum}{A Spectrum2 object to be exported}

\item{file}{File to write the spectrum to}
}
\value{
The input spectrum, invisibly
}
\description{
Uses the data in a [MSnbase::Spectrum2] object to write a CSV file in the
same format used by FreeStyle to export spectra. These files can be used to
search the mzCloud online library.
}
