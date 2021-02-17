#' Molecular descriptors
#' @param SMILES a character vector of valid SMILES
#' @importFrom dplyr mutate
#' @importFrom parallel makeCluster parLapply stopCluster
#' @importFrom purrr map_dbl map_int map_chr map
#' @export
#' @examples
#' data(aminoAcids)
#' descriptors(aminoAcids$SMILES)

descriptors <- function(SMILES){
  desc <- c('HBA1',
            'HBA2',
            'HBD',
            'logP',
            'TPSA'
  )
  
  descs <- desc %>%
    map(~{
      descType <- .
      map_dbl(SMILES,~{
        m <- .
        m %>%
          descriptor(descType)
      })
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
      g <- map_int(SMILES,~{
        s <- .
        s %>%
          smartsSearch(string)
      }) %>%
        as_tibble()
      names(g) <- .$Name
      return(g)
    }) %>%
    bind_cols()
  
  desc <- bind_cols(SMILES = SMILES,descs,groups) %>%
    mutate(Total_Charge = -Negative_Charge + Positive_Charge,
           MF = map_chr(SMILES,smileToMF),
           `Accurate_Mass` = map_dbl(SMILES,smileToAccurateMass) %>% round(5)
           ) %>%
    
    select(SMILES,MF,Accurate_Mass,Negative_Charge,Positive_Charge,Total_Charge,HBA1:TPSA,NHH:COO)
  
  return(desc)
}
