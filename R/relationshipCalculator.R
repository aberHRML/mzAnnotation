#' m/z relationship calculation
#' @description adduct, isotope and biotransfromation calculation.
#' @param mz numeric \code{vector} of accurate m/z
#' @param limit limit of deviation for thresholding associations. Defaults to 0.001
#' @param adducts character vector of adducts names to use
#' @param isotopes character vector of isotopes to use.
#' @param transformations character vector of transformations to use
#' @param adductTable table containing adduct formation rules. Defaults to \code{\link{adducts}()}.
#' @param isotopeTable table containing isotope rules. Defaults to \code{\link{isotopes}()}.
#' @param transformationTable table containing transformation rules. Defaults to \code{\link{transformations}()}.
#' @examples 
#' relationshipCalculator(c(132.03023,168.00691))
#' @author Jasen Finch
#' @export
#' @importFrom utils combn
#' @importFrom dplyr left_join contains inner_join
#' @importFrom stringr str_c
#' @importFrom tidyr expand_grid spread drop_na
#' @importFrom tibble rowid_to_column tibble
#' @importFrom stats setNames

relationshipCalculator <- function(mz, limit = 0.001, adducts = c("[M-H]1-","[M+Cl]1-","[M+K-2H]1-"), isotopes = NA, transformations = NA, adductTable = adducts(), isotopeTable = isotopes(), transformationTable = transformations()){
  
  Ms <- expand_grid(`m/z` = mz,
                    Adduct = adducts,
                    Isotope = isotopes,
                    Transformation = transformations) %>% 
    mutate(M = calculateMs(`m/z`,Adduct,Isotope,Transformation)) %>% 
    rowid_to_column(var = 'ID')
  
  relationships <- expand_grid(ID1 = Ms$ID,
                               ID2 = Ms$ID) %>%
    filter(ID1 != ID2) %>% 
  left_join(Ms %>% 
              select(ID1 = ID,M1 = M),
            by = 'ID1') %>% 
    left_join(Ms %>% 
                select(ID2 = ID,M2 = M),
              by = 'ID2') %>% 
    mutate(Error = abs(M1 - M2)) %>% 
    filter(Error <= limit) %>% 
    left_join(Ms %>% 
                setNames(str_c(names(.),'1')),
              by = c("ID1", "M1")) %>% 
    left_join(Ms %>% 
                setNames(str_c(names(.),'2')),
              by = c("ID2", "M2")) %>%
    filter(`m/z1` != `m/z2`) %>% 
    select(contains('ID'),
           contains('m/z'),
           contains('Adduct'),
           contains('Isotope'),
           contains('Transformation'),Error)
  
  unique_rel <- relationships %>% 
    select(contains('ID')) %>% 
    mutate(fill = 1) %>%
    spread(ID2,fill) %>%
    {
      id1 <- select(., ID1)
      . <- select(.,-ID1)
      
      .[lower.tri(.)] <- NA
      . <- bind_cols(.,id1)
      .
    } %>%
    gather(ID2,fill,-ID1) %>%
    drop_na() %>%
    select(-fill) %>% 
    mutate(ID2 = as.numeric(ID2))
  
  relationships <- relationships %>% 
    inner_join(unique_rel,
                     by = c("ID1", "ID2")) %>% 
    select(-contains('ID'))
  
  return(relationships)
}

calculateMs <- function(mzs, adducts, isotopes, transformations, adductTable = adducts(), isotopeTable = isotopes(), transformationTable = transformations()){
  
  if (!identical(length(mzs),length(adducts),length(isotopes),length(transformations))){
    stop('Arguments mzs, adducts, isotopes, transformations should be vectors of the same length',
         call. = FALSE) 
  }
  
  addRules <- adducts %>% 
    tibble(Adduct = .) %>% 
    left_join(adductTable ,
              by = c('Adduct' = 'Name'))
  
  isoRules <- isotopes %>% 
    tibble(Isotope = .) %>% 
    left_join(isotopeTable, by = "Isotope") %>% 
    {
      .$`Mass Difference`[is.na(.$`Mass Difference`)] <- 0
      .
    }
  
  transRules <- transformations %>% 
    tibble(Transformation  = .) %>% 
    left_join(transformationTable, by = c("Transformation" = "MF Change")) %>% 
    {
      .$Difference[is.na(.$Difference)] <- 0
      .
    }
  
  M <- ((mzs - addRules$Add) * addRules$Charge) 
  M <- (M - isoRules$`Mass Difference`)
  M <- M / addRules$xM
  M <- M - transRules$Difference
  
  return(round(M,5))
}