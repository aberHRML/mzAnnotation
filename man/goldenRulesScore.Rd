% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/golden_rules.R
\name{goldenRulesScore}
\alias{goldenRulesScore}
\title{Molecular formula plausibility scores}
\usage{
goldenRulesScore(golden_rules)
}
\arguments{
\item{golden_rules}{a tibble containing golden rule heuristic checks results as from \code{goldenRules()}}
}
\value{
A tibble containing golden rules plausibility percentage scores.
}
\description{
Percentage scores for molecular formula plausibility based on Rules 2, 4, 5 and 6 of Kind et al 2007.
}
\examples{
c('H2O','C12H22O11') \%>\% 
  goldenRules() \%>\% 
  goldenRulesScore()
}
\references{
Kind, T. and Fiehn, O., 2007. Seven Golden Rules for heuristic filtering of molecular formulas obtained by
accurate mass spectrometry. \emph{BMC bioinformatics}, \emph{8}(1), pp.1-20.
}
