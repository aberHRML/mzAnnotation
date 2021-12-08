#' @useDynLib mzAnnotation
#' @importFrom Rcpp evalCpp

globalVariables(c('Name','ID','MF','Accurate Mass','Smile 1','Adduct',
                  'PPM error','RelativeAbundance','Element','AtomicMass',
                  'Frequency','Isotope','MF Change','Probability',
                  'Relative Abundance','C','S','Elements','Error',
                  'Theoretical m/z','True','desc','m/z','Adducts',
                  'Negative_Charge','Positive_Charge','SMILES','Total_Charge',
                  'NHH','COO','ACCESSION_ID','Accurate_Mass','Measured m/z',
                  'HBA1','TPSA','.','Possible','Rule','x','Transformation','ID1',
                  'ID2','M','M1','M2','m/z1','m/z2','fill','Mass','Theoretical M',
                  'Measured M','Score','RDBE','remainder','element','valence','frequency',
                  'total_valence','sum_valence','odd_valence_total','twice_maximum_valence',
                  'twice_atoms_minus_1'
))