% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/golden_rules.R
\name{lewis}
\alias{lewis}
\alias{senior}
\title{LEWIS and SENIOR checks}
\usage{
lewis(
  element_frequencies,
  valences = list(C = 4, H = 1, N = 3, O = 2, P = 3, S = 4)
)

senior(
  element_frequencies,
  valences = list(C = 4, H = 1, N = 3, O = 2, P = 3, S = 4)
)
}
\arguments{
\item{element_frequencies}{table of element frequencies for a set of molecular formulas as returned by \code{elementFrequencies()}.}

\item{valences}{named list of element valences}
}
\value{
Boolean vector of check results for each molecular formula.
}
\description{
LEWIS molecular formula valence test and SENIOR test for the existence of molecular graphs.
}
\examples{
element_frequencies <- elementFrequencies(c('C12H22O11','C12H22NO11'))
lewis(element_frequencies)
senior(element_frequencies)
}
