% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphing_functions.R
\name{fancy_scientific}
\alias{fancy_scientific}
\title{Generate plotmath-formatted scientific notation expressions}
\usage{
fancy_scientific(l)
}
\arguments{
\item{l}{a single numeric value to be formatted}
}
\value{
an expression in plotmath format for displaying this value
}
\description{
Converts a number into scientific notation using plotmath formatting. The
resulting expression can be used in ggplot2 labels.
}
\details{
There isn't really any error-checking here. Put garbage in, you'll
  get garbage out.
}
