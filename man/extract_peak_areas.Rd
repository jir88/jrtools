% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms_alignment_functions.R
\name{extract_peak_areas}
\alias{extract_peak_areas}
\title{Extract peak areas stored as blobs in a mass spec alignment}
\usage{
extract_peak_areas(msa, ids = NULL, adj = FALSE)
}
\arguments{
\item{msa}{An ms_alignment object to query}

\item{ids}{Compound IDs to get areas for, or NULL to get all peak areas}

\item{adj}{Should QC-adjusted peak areas be returned instead of raw peak
areas? Defaults to FALSE.}
}
\value{
A tibble containing the peak areas and another containing boolean
 flag reporting whether the peak exists or not
}
\description{
Certain mass spec data alignment formats store peak areas as SQLite blobs
containing the data in binary form. This function extracts the
stored data from such blobs. The returned areas are for the reference adduct
ion for each compound. Other areas must be retrieved manually from the relevant
tables.
}
