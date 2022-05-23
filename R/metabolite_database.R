
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

setMethod('show',signature = 'MetaboliteDatabase',
          function(object){
            type <- object@type
            accessions <- object %>%
              entries() %>%
              nrow()
            cat('\n',type,' MetaboliteDatabase object containing ',accessions,' accessions\n\n',sep = '')
          }
)

#' Create a metabolite database
#' @description Build a metabolite database ready for use.
#' @param entries tibble containing accession information. If \code{type = 'remote'} this should be the name of the table containing the accession information within the SQL database.
#' @examples 
#' db <- metaboliteDB(amino_acids)
#' @importFrom dplyr tbl
#' @importFrom purrr map_chr
#' @importFrom methods new
#' @export

metaboliteDB <- function(entries){
  
  checkEntries(entries)
  
  metabolite_descriptors <- chemicalDescriptors(entries$SMILES)
  
  db <- new(
    'MetaboliteDatabase',
    entries = as_tibble(entries),
    descrtors = as_tibble(metabolite_descriptors)
  )
  
  return(db)
}

#' Retrieve database entries
#' @rdname accessors
#' @description Get or set tables in a `MetaboliteDatabase` class object.
#' @param db object of S4 class `MetaboliteDatabase`
#' 
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

#' @rdname accessors
#' @export

setGeneric('descriptors',function(db) {
  standardGeneric('descriptors')
})

#' @rdname accessors

setMethod('descriptors',signature = 'MetaboliteDatabase',
          function(db){
            db@descriptors
          }
)

#' Filter a mass range
#' @rdname utilities
#' @description Filter a MetaboliteDatabase for a given mass range.
#' @param db S4 object of class MetaboliteDatabase
#' @param lower lower mass boundary
#' @param upper upper mass boundary
#' @examples 
#' db <- metaboliteDB(amino_acids,descriptors(amino_acids$SMILES))
#' db <- filterMR(db,100,120)
#' @export

setGeneric("filterMR", function(db,lower,upper) {
  standardGeneric("filterMR")
})

#' @rdname utilities

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

#' @rdname utilities
#' @export

setGeneric("filterER", function(db,rule) {
  standardGeneric("filterER")
})

#' @rdname utilities

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

#' @rdname utilities
#' @export

setGeneric("filterIP", function(db,rule) {
  standardGeneric("filterIP")
})

#' @rdname utilities

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

#' @rdname utilities
#' @export

setGeneric("filterACCESSIONS", function(db,ids) {
  standardGeneric("filterACCESSIONS")
})

#' @rdname utilities

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

#' @rdname utilities
#' @export

setGeneric('filterMF', function(db,mf){
  standardGeneric('filterMF')
})

#' @rdname utilities

setMethod('filterMF',signature = 'MetaboliteDatabase',
          function(db,mf){
            db@descriptors[[1]] <- db@descriptors[[1]] %>%
              filter(MF %in% mf)
            
            db@accessions[[1]] <- db@accessions[[1]] %>%
              filter(ID %in% db@descriptors[[1]]$ID)
            
            return(db)
          }
)

#' @rdname utilities
#' @export

setGeneric('calcAdducts',function(db,id,adduct_rules_table = adduct_rules())
  standardGeneric('calcAdducts'))

#' @rdname utilities

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

#' @rdname utilities
#' @importFrom dplyr bind_rows select filter
#' @export

setGeneric('PIPsearch',function(db,
                                mz,
                                ppm,
                                adduct,
                                isotope = NA, 
                                adduct_rules_table = adduct_rules(),
                                isotope_rules_table = isotope_rules()) 
  standardGeneric('PIPsearch')
)

#' @rdname utilities

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

#' @rdname utilities
#' @export

setGeneric("elementFreq", function(db) {
  standardGeneric("elementFreq")
})

#' @rdname utilities
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
