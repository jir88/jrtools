% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_spectral_comparison.R
\name{plot_spectral_comparison}
\alias{plot_spectral_comparison}
\title{Plot two mass spectra}
\usage{
plot_spectral_comparison(spec1, spec2, match_tol_ppm = 5)
}
\arguments{
\item{spec1}{Spectrum to plot on top half of graph. Must have columns named
'mz' and 'intensity'.}

\item{spec2}{Spectrum to plot on bottom half of graph. Must have columns named
'mz' and 'intensity'.}

\item{match_tol_ppm}{Fragment matching tolerance in PPM}
}
\value{
The ggplot object for rendering or further modification.
}
\description{
Plots two mass spectra back-to-back for comparison. Fragments within a
given PPM tolerance of each other are marked.
}
