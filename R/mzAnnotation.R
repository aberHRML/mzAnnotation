#' mzAnnotation
#' @name mzAnnotation
#' @useDynLib mzAnnotation
#' @importFrom Rcpp evalCpp
#' @name mzAnnotation

globalVariables(c('Name','ID','MF','Accurate Mass','Smile 1','Adduct',
                  'PPM Error','RelativeAbundance','Element','AtomicMass',
                  'Frequency','Isotope','MF Change','Probability',
                  'Relative Abundance','C','S','Elements','Error',
                  'Theoretical m/z','True','desc','m/z','Adducts',
                  'Negative Charge','Positive Charge','SMILE','Total Charge',
                  'NHH','COO','ACCESSION_ID'
                  ))