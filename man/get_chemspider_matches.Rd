% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms_alignment_functions.R
\name{get_chemspider_matches}
\alias{get_chemspider_matches}
\title{Get ChemSpider matches associated with particular compounds}
\usage{
get_chemspider_matches(msa, ids = NULL)
}
\arguments{
\item{msa}{An ms_alignment object to query}

\item{ids}{Compound IDs to get matches for, or NULL to get all matches}
}
\value{
A tibble with the ChemSpider matches
}
\description{
Queries a mass spec alignment database to get any ChemSpider matches associated
with some or all compounds.
}
