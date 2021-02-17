#' Ionisation products
#' @description Calculate ionisation products for a given SMILES.
#' @param SMILES a valid SMILES string
#' @param adductTable table of adduct rules. Defaults to adducts()
#' @examples 
#' ionisationProducts(aminoAcids$SMILES[1])
#' @importFrom dplyr ungroup
#' @export

ionisationProducts <- function(SMILES,adductTable = adducts()){
  desc <- descriptors(SMILES)
  
  adductTable %>%
    select(Name,Rule) %>%
    bind_cols(desc) %>%
    rowwise() %>%
    mutate(Possible = eval(parse(text = Rule)),
           `m/z` = calcMZ(Accurate_Mass,Name,adductTable = adductTable),
           MF = adductTransformMF(MF,Name,Adducts = adductTable)) %>%
    select(Adduct = Name,`m/z`,MF,Possible) %>%
    ungroup()
}
