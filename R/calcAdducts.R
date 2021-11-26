#' Calculate adducts
#' @rdname calcAdducts
#' @description Calculate ionisation products for given database accession.
#' @param db object of class \code{MetaboliteDatabase}
#' @param id accession id
#' @param adduct_rules_table table of adduct rules. Defaults to adducts()
#' @examples 
#' db <- metaboliteDB(amino_acids,descriptors(amino_acids$SMILES))
#' add <- calcAdducts(db,1)
#' @importFrom purrr map_df
#' @importFrom tibble deframe
#' @export

setGeneric('calcAdducts',function(db,id,adduct_rules_table = adduct_rules())
  standardGeneric('calcAdducts'))

#' @rdname calcAdducts

setMethod('calcAdducts',signature = 'MetaboliteDatabase',
          function(db,id,adduct_rules_table = adduct_rules()){
            
            smiles <- db %>%
              getDescriptors() %>%
              filter(ID == id) %>%
              select(SMILES) %>%
              deframe()
            
            smiles %>%
              ionisationProducts(adduct_rules_table = adduct_rules_table)
          })
