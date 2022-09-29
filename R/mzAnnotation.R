#' @useDynLib mzAnnotation
#' @importFrom Rcpp evalCpp

globalVariables(c('element','count','MF','CHO','Name','.','RelativeAbundance','Element',
                  'AtomicMass','Frequency','Isotope','MF Change','name','probability_check',
                  'heuristic_check','result','check','label','ratio','operator','threshold',
                  'Mass','PPM error','CHO proportion','rule','score','LEWIS and SENIOR',
                  'Element ratios','Element counts','Plausibility (%)','Theoretical M',
                  'Measured m/z','Theoretical m/z','isotope_possible','Adduct','Measured M',
                  'Relative Abundance','m/z','possible','RDBE','remainder','valence','frequency',
                  'value','Abundance','Transformation','Probability','ID1','ID2','ID','M','M1',
                  'M2','Error','m/z1','m/z2','fill','Transformation1','Transformation2','total_valence',
                  'sum_valence','odd_valence_total','twice_maximum_valence','twice_atoms_minus_1',
                  'C','S'))