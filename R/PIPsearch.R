#' Putative ionisation product search
#' @rdname PIPsearch
#' @param db object of class `MetaboliteDatabase`.
#' @param mz the accurate m/z to search.
#' @param ppm the parts per million threshold to search.
#' @param adduct the adduct name to search.
#' @param isotope the isotope name to search. Defaults to NA for non-isotopic searches.
#' @param adduct_rules_table adduct formation rules tabl. Defaults to table returned by `adducts_rules()`.
#' @param isotope_rules_table isotope table containing available isotope rules. Defaults to table returned by `isotope_rules()`.
#' @export
#' @author  Jasen Finch
#' @importFrom dplyr bind_rows select filter
#' @importFrom magrittr %>%
#' @examples
#' res <- PIPsearch(metaboliteDB(amino_acids,descriptors(amino_acids$SMILES)),132.03023,5,'[M-H]1-')

setGeneric('PIPsearch',function(db,
                                mz,
                                ppm,
                                adduct,
                                isotope = NA, 
                                adduct_rules_table = adduct_rules(),
                                isotope_rules_table = isotope_rules()) 
  standardGeneric('PIPsearch')
)

#' @rdname PIPsearch

setMethod('PIPsearch',signature = 'MetaboliteDatabase',
          function(db,
                   mz,
                   ppm,
                   adduct,
                   isotope = NA, 
                   adduct_rules_table = adduct_rules(),
                   isotope_rules_table = isotope_rules()){
  M <- calcM(mz,
             adduct = adduct,
             isotope = isotope,
             adduct_rules_table = adduct_rules_table,
             isotope_rules_table =  isotope_rules_table)
  mr <- ppmRange(M,ppm)
  
  res <- db %>%
    filterMR(mr$lower,mr$upper)
  
  if (!is.na(isotope) & nrow(res@accessions[[1]]) > 0) {
    isoRule <- isotope_rules_table$Rule[isotope_rules_table$Isotope == isotope]
    res <- res %>%
      filterER(isoRule)
  }
  
  addRule <- adduct_rules_table$Rule[adduct_rules_table$Name == adduct]
  
  res <- res %>%
    filterIP(addRule)
  
  res <- res %>%
  {left_join(.@accessions[[1]],.@descriptors[[1]],by = c("ID", "SMILES"))} %>%
    select(ID:Accurate_Mass) %>%
    mutate(Isotope = isotope,
           Adduct = adduct,
           `Measured m/z` = mz,
           `Theoretical m/z` = calcMZ(Accurate_Mass,adduct,isotope),
           `PPM Error` = ppmError(`Measured m/z`,`Theoretical m/z`)
    ) 
  return(res)
})
