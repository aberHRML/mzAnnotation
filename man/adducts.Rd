% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rules.R
\name{adduct_rules}
\alias{adduct_rules}
\alias{adduct_names}
\title{Adduct formation rules}
\format{
A tibble containing 58 rows and 6 columns
\describe{
\item{Name}{Adduct name}
\item{Charge}{Charge added}
\item{xM}{Number of M needed}
\item{Add}{Mass value change}
\item{Nelec}{Electron change}
\item{AddAt}{pseudo formula of the atoms to be added to the molecular formula of one M}
\item{RemAt}{pseudo formula of the atoms to be removed to the molecular formula of one M}
\item{AddEx}{pseudo formula of the atoms to be added to the final ionisation product molecular formula}
\item{RemEx}{pseudo formula of the atoms to be removed to the final ionisation product molecular formula}
\item{Rule}{Structural rule for formation}
}
}
\usage{
adduct_rules()

adduct_names()
}
\description{
Formation rules for electrospray ionisation adducts
}
\examples{
adduct_names()

adduct_rules()
}
