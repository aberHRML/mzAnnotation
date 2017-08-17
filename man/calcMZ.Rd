% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcMZ.R
\name{calcMZ}
\alias{calcMZ}
\title{calcM}
\usage{
calcMZ(M, adduct = "[M+H]1+", isotope = NA, transformation = NA,
  adducts = mzAnnotation::Adducts, isotopes = mzAnnotation::Isotopes,
  transformations = mzAnnotation::Transformations)
}
\arguments{
\item{M}{M for which to calculate an m/z}

\item{adduct}{adduct to apply}

\item{isotope}{isotope to apply}

\item{transformation}{transformation to apply}

\item{adducts}{adduct table containing available adduct rules. Defaults to \code{Adducts}.}

\item{isotopes}{isotope table containing available isotope rules. Defaults to \code{Isotopes}.}

\item{transformations}{transformations table containing available transformations rules. Defaults to \code{Transformations}.}
}
\description{
calculate an m/z for a given M, adduct, isotope and transformation
}
\examples{
calcMZ(116.05182,adduct = '[M+H]1+',isotope = 'C13',transformation = 'M - [O] + NH2]')
}