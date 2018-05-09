#' convert SMILES to a series of molecular descriptors
#'
#' @param smiles vector of valid SMILEs
#' @importFrom parallel makeCluster parLapply stopCluster
#' @importFrom purrr map_dbl map_int
#' @export
#' @examples
#' data(aminoAcids)
#' descriptors(aminoAcids$SMILE)

descriptors <- function(smiles){
 
  desc <- c('HBA1',
            'HBA2',
            'HBD',
            'logP',
            'TPSA'
            )

 descs <- desc %>%
   map(~{
     descType <- .
     map_dbl(smiles,~{
       m <- .
       m %>%
         descriptor(descType)
     })
   })
 names(descs) <- desc

 descs <- descs %>%
   bind_cols()

  Fgroups <- tibble(Name = c('Negative Charge',
                             'Positive Charge',
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
      g <- map_int(smiles,~{
        s <- .
        s %>%
          smartsSearch(string)
      }) %>%
        as_tibble()
      names(g) <- .$Name
      return(g)
    }) %>%
    bind_cols()

  desc <- bind_cols(SMILE = smiles,descs,groups) %>%
    mutate(`Total Charge` = -`Negative Charge` + `Positive Charge`,
           MF = map_chr(SMILE,smileToMF),
           `Accurate Mass` = map_dbl(MF,calcAccurateMass)) %>%
    select(SMILE,MF,`Accurate Mass`,`Negative Charge`,`Positive Charge`,`Total Charge`,NHH:COO)

  return(desc)
}
