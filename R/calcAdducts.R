#' Calculate adducts
#' @rdname calcAdducts
#' @description Calculate ionisation products for given database accession.
#' @param db object of class \code{MetaboliteDatabase}
#' @param id accession id
#' @param adductTable table of adduct rules. Defaults to adducts()
#' @examples 
#' db <- metaboliteDB(aminoAcids,descriptors(aminoAcids$SMILES))
#' add <- calcAdducts(db,1)
#' @importFrom purrr map_df
#' @importFrom tibble deframe
#' @export

setMethod('calcAdducts',signature = 'MetaboliteDatabase',
          function(db,id,adductTable = adducts()){
            
            smiles <- db %>%
              getDescriptors() %>%
              filter(ID == id) %>%
              select(SMILES) %>%
              deframe()
            
            smiles %>%
              ionisationProducts(adductTable = adductTable)
          })