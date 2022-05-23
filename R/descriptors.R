#' Chemical descriptors
#' @param SMILES a character vector of valid SMILES
#' @importFrom dplyr mutate
#' @importFrom purrr map
#' @importFrom furrr future_map_dbl future_map_int future_map_chr furrr_options
#' @export
#' @examples
#' data(amino_acids)
#' chemicalDescriptors(amino_acids$SMILES)

chemicalDescriptors <- function(SMILES){
  desc <- c('HBA1',
            'HBA2',
            'HBD',
            'logP',
            'TPSA'
  )
  
  descs <- desc %>%
    map(~{
      descType <- .
      future_map_dbl(SMILES,~{
        m <- .
        m %>%
          descriptor(descType)
      },
      .options = furrr_options(seed = TRUE))
    })
  names(descs) <- desc
  
  descs <- descs %>%
    bind_cols()
  
  Fgroups <- tibble(Name = c('Negative_Charge',
                             'Positive_Charge',
                             'NHH',
                             'OH',
                             'COOH',
                             'COO'),
                    String = c('[-]',
                               '[+]',
                               "[NX3;H2]",
                               "[OX2H]",
                               "[CX3](=O)[OX2H1]",
                               "[CX3](=O)[OX1H0-]")
  )
  
  
  groups <- Fgroups %>%
    split(1:nrow(Fgroups)) %>%
    map(~{
      string <- .$String
      g <- future_map_int(SMILES,~{
        s <- .
        s %>%
          smartsSearch(string)
      },
      .options = furrr_options(seed = TRUE)) %>%
        as_tibble()
      names(g) <- .$Name
      return(g)
    }) %>%
    bind_cols()
  
  desc <- bind_cols(SMILES = SMILES,descs,groups) %>%
    mutate(Total_Charge = -Negative_Charge + Positive_Charge,
           MF = future_map_chr(SMILES,smilesToMF,
                               .options = furrr_options(seed = TRUE)),
           `Accurate_Mass` = future_map_dbl(SMILES,smilesToAccurateMass,
                                            .options = furrr_options(seed = TRUE)) %>% round(5)
           ) %>%
    
    select(SMILES,MF,Accurate_Mass,Negative_Charge,Positive_Charge,Total_Charge,HBA1:TPSA,NHH:COO)
  
  return(desc)
}
