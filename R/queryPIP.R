#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select

queryPIP <- function(mz, ppm, add, iso = NA, adducts = mzAnnotation::Adducts, isotopes = mzAnnotation::Isotopes, DB = mzAnnotation::MZedDB, filter = T){
  # Retrieve database entries
  elementFreq <- DB$ELementFrequencies
  metaboliteRules <- DB$Rules
  DB <- DB$DB
  
  # Calculate the mz range over which to search based on ppm error
  mass <- ppmRange(mz,ppm)
  
  # Calculate M
  mass <- map(mass,calcM,adduct = add,isotope = iso)
  
  # Filter database entries based mz range
  DB <- filter(DB, `Accurate Mass` > mass$lower & `Accurate Mass` < mass$upper)
  
  # If an isotope has been selected; filter the database entries based on the isotope rules
  if (!is.na(iso)) {
    elementFreq <- elementFreq[elementFreq$ID %in% DB$ID,]
    C <- elementFreq[,'C']
    O <- elementFreq[,'O']
    Cl <- elementFreq[,'Cl']
    K <- elementFreq[,'K']
    S <- elementFreq[,'S']
    
    isotopes <- filter(isotopes,Isotope == iso)
    DB <- filter(DB,eval(parse(text = isotopes$Rule)))
  }
  
  # Filter the database based on the adduct formation rules
  metaboliteRules <- filter(metaboliteRules, ID %in% DB$ID)
  Nch <- metaboliteRules[,"Nch"]
  Nacc <- metaboliteRules[,"Nacc"]
  Ndon <- metaboliteRules[,"Ndon"]
  Nnhh <- metaboliteRules[,"Nnhh"]
  Noh <-  metaboliteRules[,"Noh"]
  Ncooh <- metaboliteRules[,"Ncooh"]
  Ncoo <-  metaboliteRules[,"Ncoo"]
  
  adducts <- filter(adducts,Name == add)
  
  DB <- filter(DB,eval(parse(text = adducts$Rule))) %>% 
    filterPIP() %>%
    mutate(Adduct = add, 
           Isotope = iso, 
           `Theoretical m/z` = calcMZ(`Accurate Mass`,add,iso),
           `PPM Error` = ppmError(mz,`Theoretical m/z`)
    ) %>%
    select(ID,Name,MF,`Accurate Mass`,`Smile 1`,Adduct:`PPM Error`)
  
  DB$`PPM Error` <- round(DB$`PPM Error`,5)
  
  colnames(DB)[5] <- 'Smile'
  
  return(DB)
}
