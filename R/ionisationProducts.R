#' Ionisation products
#' @description Calculate ionisation products for a given SMILES.
#' @param SMILES a valid SMILES string
#' @param adduct_rules_table table of adduct rules. Defaults to adducts()
#' @examples 
#' ionisationProducts(amino_acids$SMILES[1])
#' @importFrom dplyr ungroup rowwise
#' @export

ionisationProducts <- function(SMILES,adduct_rules_table = adduct_rules()){
  
  desc <- chemicalDescriptors(SMILES)
  
  adduct_rules_table %>%
    select(Name,Rule) %>%
    bind_cols(desc) %>%
    rowwise() %>%
    mutate(Possible = eval(parse(text = Rule)),
           `m/z` = calcMZ(Accurate_Mass,
                          Name,
                          adduct_rules_table = adduct_rules_table),
           MF = adductTransformMF(MF,
                                  Name,
                                  adduct_rules_table = adduct_rules())) %>%
    select(Adduct = Name,`m/z`,MF,Possible) %>%
    ungroup()
}
