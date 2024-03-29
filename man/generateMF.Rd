% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/molecular_formulas.R
\name{generateMF}
\alias{generateMF}
\title{Molecular formula generation}
\usage{
generateMF(
  mass,
  ppm = 1,
  charge = 0,
  element_ranges = suitableElementRanges(mass)
)
}
\arguments{
\item{mass}{molecular mass}

\item{ppm}{parts per million error tolerance}

\item{charge}{total charge}

\item{element_ranges}{named list of element ranges}
}
\value{
A \code{tibble} containing the generated molecular formulas, their theoretical mass and their PPM error.
}
\description{
Exhaustive molecular formula generation for a specified molecular mass.
}
\examples{
generateMF(342.11621,
           element_ranges = list(C = c(0,12),
                                 H = c(0,22),
                                 O = c(0,11)))
}
