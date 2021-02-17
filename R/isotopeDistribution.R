#' Isotopic distribution calculator
#' @param MF the molecular formular to generate the isotope distribution
#' @param charge the charge of the molecular formula
#' @param limit the relative abundance threshold
#' @param elementTable the table containing elemental information. Defaults to \code{Elements}.
#' @importFrom stringr str_split str_extract str_replace_all
#' @importFrom dplyr bind_cols group_by summarise right_join rename
#' @export
#' @author Jasen Finch
#' @examples 
#' res <- isotopeDistribution('C4H5O5',charge = -1)

isotopeDistribution <- function(MF,charge, limit = 0.00009 , elementTable = elements()) {
  elementTable <- elementTable %>%
    mutate(Name = str_c(round(AtomicMass),Element))
  atomFrequencies <- tibble(Element = names(count.elements(MF)),Frequency = count.elements(MF))
  
  Isotopes <- filter(elementTable,Element %in% atomFrequencies$Element, RelativeAbundance != 1) %>%
    mutate(Isotope = str_c(round(AtomicMass),Element))
  
  Isotopes <- map(unique(Isotopes$Element),~{
    element <- .
    numberEle <- atomFrequencies$Frequency[atomFrequencies$Element == element]
    availableIsotopes <- Isotopes %>%
      filter(Element == element) 
    t <- map(availableIsotopes$Isotope, ~{
      abund <- availableIsotopes$RelativeAbundance[availableIsotopes$Isotope == .]
      relativeAbundance <- 1
      for(i in 2:(numberEle + 1)) {
        relativeAbundance[i] <- relativeAbundance[i - 1] * abund * (numberEle - (i - 1) + 1) / (i - 1)
      }
      relativeAbundance <- relativeAbundance[-1]
      iso <- tibble(Element = rep(element,length(relativeAbundance)),`Relative Abundance` = relativeAbundance,Isotope = rep(.,length(relativeAbundance)), Frequency = 1:length(relativeAbundance)) %>%
        mutate(Name = str_c(Isotope,Frequency,sep = ' '))
    }) %>%
      bind_rows()
    return(t) 
  }) %>%
    bind_rows() %>%
    filter(`Relative Abundance` > limit)
  
  combinations <- map(2:length(unique(Isotopes$Isotope)),~{
    numberCombinations <- .
    com <- Isotopes$Name %>% 
      combn(m = numberCombinations) %>%
      t() %>%
      {
        suppressMessages(as_tibble(.,.name_repair = 'unique'))
      }
    differentCombinations <- map(com,~{
      unlist(map(str_split(.,' '),~{
        .[[1]][1]
      }))
    }) %>%
      bind_cols() %>%
      apply(1,function(x){
        T %in% (table(x) > 1)
      })
    com <- com %>%
      filter(differentCombinations == FALSE) %>%
      apply(1,str_c,collapse = '; ')
    return(com)
  }) %>%
    unlist() %>%
    map(~{
      isos <- unlist(str_split(.,'; '))
      isos <- Isotopes %>%
        filter(Name %in% isos)
      tibble(Name = ., `Relative Abundance` = prod(isos$`Relative Abundance`))
    }) %>%
    bind_rows() %>%
    filter(`Relative Abundance` > limit)
  
  Isotopes <- Isotopes %>%
    select(Name, `Relative Abundance`) %>%
    bind_rows(combinations) %>%
    rename(Isotope = Name) %>% 
    bind_rows(tibble(Isotope = NA, `Relative Abundance` = 1)) %>%
    arrange(desc(`Relative Abundance`)) %>%
    mutate(Probability = `Relative Abundance`/sum(`Relative Abundance`))
  
  mzs <- map(Isotopes$Isotope,~{
    if(is.na(.)){
      calcAccurateMass(MF,charge = charge)
    } else {
      Name <- unlist(map(str_split(., '; '),~{unlist(map(str_split(.,' '),~{.[1]}))}))
      Frequency <- as.integer(unlist(map(str_split(., '; '),~{unlist(map(str_split(.,' '),~{.[2]}))})))
      isoDat <- tibble(Name = Name,Frequency = Frequency) %>%
        mutate(Element = str_replace_all(Name,'[:digit:]','')) %>%
        left_join(select(elementTable,Name,Element,AtomicMass),by = c('Name' = 'Name', 'Element' = 'Element')) %>%
        mutate(AtomicMass = AtomicMass * Frequency)
      isoMass <- sum(isoDat$AtomicMass)
      newMF <- isoDat %>%
        group_by(Element) %>%
        summarise(Frequency = sum(Frequency)) %>%
        right_join(atomFrequencies,by = c('Element' = 'Element'),suffix = c('1','2'))
      newMF$Frequency1[is.na(newMF$Frequency1)] <- 0
      newMF <- newMF %>%
        mutate(Frequency = Frequency2 - Frequency1,MF = str_c(Element,Frequency)) %>%
        select(MF) %>% 
        unlist() %>%
        str_c(collapse = '') %>%
        calcAccurateMass(charge = charge)
      mz <- isoMass + newMF
      
      return(mz)
    }
  }) %>%
    unlist() %>%
    round(5)
  
  Isotopes <- Isotopes %>%
    mutate(`m/z` = mzs) %>%
    select(Isotope,`m/z`,`Relative Abundance`,Probability)
  return(Isotopes)
}
