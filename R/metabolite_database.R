
setClass(
  'MetaboliteDatabase',
  slots = list(
    entries = 'tbl_df',
    descriptors = 'tbl_df'
  ),
  prototype = list(
    entries = tibble(),
    descriptors = tibble()
  ))

#' Create a metabolite database
#' @description Build a metabolite database ready for use.
#' @param entries tibble containing accession information. If \code{type = 'remote'} this should be the name of the table containing the accession information within the SQL database.
#' @examples 
#' db <- metaboliteDB(amino_acids,descriptors(amino_acids$SMILES))
#' @importFrom dplyr tbl
#' @importFrom purrr map_chr
#' @importFrom methods new
#' @export

metaboliteDB <- function(entries){
  
  db <- new(
    'MetaboliteDatabase',
    entries = as_tibble(entries)
  )
  
  descriptors(db) <- chemicalDescriptors(entries$SMILES)
  
  return(db)
}

#' Retrieve entries
#' @rdname accessors
#' @description return accession data table from MetaboliteDatabase object
#' @param db MetaboliteDatabase
#' @export

setGeneric('entries',function(db) {
  standardGeneric('entries')
})

#' @rdname accessors

setMethod('entries',signature = 'MetaboliteDatabase',
          function(db){
            db@entries
          }
)

#' Retrieve descriptors
#' @rdname getDescriptors
#' @description return descriptor data table from MetaboliteDatabase object
#' @param db MetaboliteDatabase
#' @export

setMethod('getDescriptors',signature = 'MetaboliteDatabase',
          function(db){
            db@descriptors[[1]]
          }
)

#' Filter a mass range
#' @rdname filterMR
#' @description Filter a MetaboliteDatabase for a given mass range.
#' @param db S4 object of class MetaboliteDatabase
#' @param lower lower mass boundary
#' @param upper upper mass boundary
#' @examples 
#' db <- metaboliteDB(amino_acids,descriptors(amino_acids$SMILES))
#' db <- filterMR(db,100,120)
#' @export

setMethod('filterMR',signature = 'MetaboliteDatabase',
          function(db,lower,upper){
            desc <- db@descriptors[[1]]
            desc <- desc %>%
              filter(Accurate_Mass > lower & Accurate_Mass < upper)
            acc <- db@accessions[[1]]
            acc <- acc %>%
              filter(SMILES %in% desc$SMILES)
            db@descriptors <- list(desc)
            db@accessions <- list(acc)
            return(db)
          }
)

#' Filter using an elemental rule
#' @rdname filterER
#' @description Filter a MetaboliteDatabase based on an elemental rule.
#' @param db S4 object of class MetaboliteDatabase
#' @param rule elemental rule given as a string
#' @examples 
#' db <- metaboliteDB(amino_acids,descriptors(amino_acids$SMILES))
#' db <- filterER(db,'S>0')
#' @export

setMethod('filterER',signature = 'MetaboliteDatabase',
          function(db,rule){
            ef <- elementFreq(db)
            if (str_extract(rule,'[:alpha:]') %in% colnames(ef)) {
              ef <- ef %>%
                filter(eval(parse(text = rule)))   
            } else {
              ef[0,]
            }
            db@descriptors[[1]] <- db@descriptors[[1]] %>%
              filter(ID %in% ef$ID)
            db@accessions[[1]] <- db@accessions[[1]] %>%
              filter(ID %in% ef$ID)
            return(db)
          }
)

#' Filter using an ionisation rule
#' @rdname filterIP
#' @description Filter MetaboliteDatabase based on an ionisation rule
#' @param db S4 object of class MetaboliteDatabase
#' @param rule Character containing ionisation rule.
#' @examples 
#' rule <- adduct_rules()$Rule[52]
#' db <- metaboliteDB(amino_acids,descriptors(amino_acids$SMILES))
#' db <- filterIP(db,rule)
#' @export

setMethod('filterIP',signature = 'MetaboliteDatabase',
          function(db,rule){
            desc <- db@descriptors[[1]]
            desc <- desc %>%
              filter(eval(parse(text = rule)))
            acc <- db@accessions[[1]]
            acc <- acc %>%
              filter(SMILES %in% desc$SMILES)
            db@descriptors <- list(desc)
            db@accessions <- list(acc)
            return(db)
          }
)

#' Filter accessions
#' @rdname filterACCESSIONS
#' @description Filter a MetaboliteDatabase based on given accession IDs.
#' @param db S4 object of class MetaboliteDatabase
#' @param ids vector of accession IDs
#' @examples 
#' db <- metaboliteDB(amino_acids,descriptors(amino_acids$SMILES))
#' db <- filterACCESSIONS(db,c(1,2))
#' @export

setMethod('filterACCESSIONS',signature = 'MetaboliteDatabase',
          function(db,ids){
            db@accessions[[1]] <- db %>%
              entries() %>%
              filter(ID %in% ids)
            db@descriptors[[1]] <- db %>%
              getDescriptors() %>%
              filter(ID %in% ids)
            return(db)
          }
)

#' Filter a molecular formula
#' @rdname filterMF
#' @description Filter a MetaboliteDatabase based on given molecular formulas
#' @param db S4 object of class MetaboliteDatabase
#' @param mf character vector of molecular formulas
#' @examples 
#' db <- metaboliteDB(amino_acids,descriptors(amino_acids$SMILES))
#' db <- filterMF(db,c('C3H7NO2','C5H10N2O3'))
#' @export

setMethod('filterMF',signature = 'MetaboliteDatabase',
          function(db,mf){
            db@descriptors[[1]] <- db@descriptors[[1]] %>%
              filter(MF %in% mf)
            
            db@accessions[[1]] <- db@accessions[[1]] %>%
              filter(ID %in% db@descriptors[[1]]$ID)
            
            return(db)
          }
)

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

setGeneric("elementFreq", function(db) {
  standardGeneric("elementFreq")
})

#' @importFrom dplyr everything

setMethod('elementFreq',signature = 'MetaboliteDatabase',
          function(db){
            MFs <- db %>%
              getDescriptors() %>%
              .$MF %>%
              unique() %>%
              map(~{
                mf <- .
                mf %>%
                  count.elements() %>%
                  as.list() %>%
                  as_tibble()
              })
            names(MFs) <- db %>%
              getDescriptors() %>%
              .$MF %>% 
              unique()
            MFs <- MFs %>% 
              bind_rows(.id = 'MF') %>%
              right_join(db %>%
                           getDescriptors() %>%
                           select(ID,MF), by = "MF") %>%
              select(ID,everything())
            return(MFs)
          })
