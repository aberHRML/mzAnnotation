#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select rename

queryPIP <- function(mz, ppm, add, iso = NA, adducts = mzAnnotation::Adducts, isotopes = mzAnnotation::Isotopes, DB = mzAnnotation::MZedDB, filter = T){
  # Retrieve database entries
  accessions <- DB@accessions
  descriptors <- DB@descriptors
  
  # elementFreq <- DB$ELementFrequencies
  
  # Calculate the mz range over which to search based on ppm error
  mass <- ppmRange(mz,ppm)
  
  # Calculate M
  mass <- map(mass,calcM,adduct = add,isotope = iso)
  
  # Filter database entries based mz range
  hits <- descriptors %>%
    filter(Accurate_Mass > mass$lower & Accurate_Mass < mass$upper)
  
  # If an isotope has been selected; filter the database entries based on the isotope rules
  if (!is.na(iso)) {
    elementFreq <- elementFreq[elementFreq$ID %in% DataBase$ID,]
    C <- elementFreq[,'C']
    O <- elementFreq[,'O']
    Cl <- elementFreq[,'Cl']
    K <- elementFreq[,'K']
    S <- elementFreq[,'S']
    
    isotopes <- filter(isotopes,Isotope == iso)
    DataBase <- DataBase %>%
      filter(eval(parse(text = isotopes$Rule)))
  }
  
  # Filter the database based on the adduct formation rules
  if (nrow(hits) > 0) {
    hits <- hits %>%
      rowwise() %>%
      mutate(True = applyIonisationRule(ID,add,adducts = adducts,DB = DB)) %>%
      filter(True == T)
    
    DataBase <- DataBase %>%
      filter(ID %in% metaboliteRules$ID) %>% 
      mutate(Adduct = add, 
             Isotope = iso, 
             `Theoretical m/z` = calcMZ(`Accurate Mass`,add,iso),
             `PPM Error` = ppmError(mz,`Theoretical m/z`)
      )
  } else {
    DataBase <- DataBase %>%
      mutate(Adduct = character(),
             Isotope = character(),
             `Theoretical m/z` = numeric(),
             `PPM Error` = numeric()
             )
  }
  
  return(DataBase)
}
