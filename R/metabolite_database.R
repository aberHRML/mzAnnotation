
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
            entries <- object %>%
              entries() %>%
              nrow()
            cat('\nMetaboliteDatabase object containing ',entries,' entries\n\n',sep = '')
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
  
  entries <- as_tibble(entries)
  
  db <- new(
    'MetaboliteDatabase',
    entries = entries,
    descriptors = bind_cols(
      select(entries,ID),
      as_tibble(metabolite_descriptors)
    )
  )
  
  return(db)
}

#' Retrieve database entries
#' @rdname accessors
#' @description Get or set tables in a `MetaboliteDatabase` class object.
#' @param db object of S4 class `MetaboliteDatabase`
#' 
#' @export

setGeneric('entries',function(db)
  standardGeneric('entries'))

#' @rdname accessors

setMethod('entries',signature = 'MetaboliteDatabase',
          function(db){
            db@entries
          }
)

setGeneric('entries<-',function(db,value)
  standardGeneric('entries<-'))

setMethod('entries<-',signature = 'MetaboliteDatabase',
          function(db,value){
            db@entries <- value
            return(db)
          }
)

#' @rdname accessors
#' @export

setGeneric('descriptors',function(db) 
  standardGeneric('descriptors'))

#' @rdname accessors

setMethod('descriptors',signature = 'MetaboliteDatabase',
          function(db){
            db@descriptors
          }
)

setGeneric('descriptors<-',function(db,value)
  standardGeneric('descriptors<-'))

setMethod('descriptors<-',signature = 'MetaboliteDatabase',
          function(db,value){
            db@descriptors <- value
            return(db)
          }
)

#' Metabolite database utilities
#' @rdname utilities
#' @description Utilities for working with metabolite databases.
#' @param db S4 object of class `MetaboliteDatabase`
#' @param IDs a numeric vector of entry IDs
#' @param mf a molecular formula to filter
#' @param rule a filtering expression
#' @param lower lower mass boundary
#' @param upper upper mass boundary
#' @return An S4 object of class `MetaboliteDatabase`.
#' @examples 
#' metabolite_database <- metaboliteDB(amino_acids)
#' 
#' ## Return the number of database entries
#' nEntries(metabolite_database)
#' 
#' ## Filter database entries
#' filterEntries(metabolite_database,c(1:5))
#' 
#' ## Filter database using a mass range
#' filterMR(metabolite_database,100,120)
#' 
#' ## Filter the database by an element frequency rule
#' filterER(metabolite_database,C > 2)
#' 
#' ## Filter the database by an ionisation product rule
#' filterIP(metabolite_database,HBA2>0 & Total_Charge==0)
#' 
#' ## Filter a database by a molecular formula
#' filterMF(metabolite_database,"C3H7NO2")
#' @export

setGeneric("nEntries", function(db) {
  standardGeneric("nEntries")
})

#' @rdname utilities

setMethod('nEntries',signature = 'MetaboliteDatabase',
          function(db){
          nrow(entries(db))
          }
)

#' @rdname utilities
#' @export

setGeneric("filterMR", function(db,lower,upper) {
  standardGeneric("filterMR")
})

#' @rdname utilities

setMethod('filterMR',signature = 'MetaboliteDatabase',
          function(db,lower,upper){
            desc <- descriptors(db) %>%
              filter(Accurate_Mass > lower,
                     Accurate_Mass < upper)
            
            acc <- entries(db) %>%
              filter(SMILES %in% desc$SMILES)
            
            descriptors(db) <- desc
            entries(db) <- acc
            return(db)
          }
)

#' @rdname utilities
#' @export

setGeneric("filterER", function(db,rule) {
  standardGeneric("filterER")
})

#' @rdname utilities
#' @importFrom rlang enexpr expr_text
#' @importFrom stringr str_extract_all

setMethod('filterER',signature = 'MetaboliteDatabase',
          function(db,rule){
            
            rule <- enexpr(rule)
            
            ef <- elementFreq(db)
            
            if (all(str_extract_all(expr_text(rule),'[:alpha:]')[[1]] %in% 
                    colnames(ef))) {
              ef <- ef %>%
                filter(!!rule)   
            } else {
              ef[0,]
            }
           
            db <- filterEntries(db,
                                IDs = ef)
            
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
            
            rule <- enexpr(rule)
            
            desc <- descriptors(db) %>%
              filter(!!rule)
            
            db <- filterEntries(
              db,
              desc$ID)
            
            return(db)
          }
)

#' @rdname utilities
#' @export

setGeneric("filterEntries", function(db,IDs) {
  standardGeneric("filterEntries")
})

#' @rdname utilities

setMethod('filterEntries',signature = 'MetaboliteDatabase',
          function(db,IDs){
            entries(db) <- db %>%
              entries() %>%
              filter(ID %in% IDs)
            
            descriptors(db) <- db %>%
              descriptors() %>%
              filter(ID %in% IDs)
            
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
            mfs <- descriptors(db) %>% 
              filter(MF == mf)
            
            db <- filterEntries(
              db,
              mfs$ID
            )
            
            return(db)
          }
)

#' Metabolite database ionisation product searches
#' @rdname products
#' @description Methods for metabolite database ionisation product searches
#' @param db S4 object of class `MetaboliteDatabase`
#' @param id Database entry ID
#' @param mz search mass to charge ratio
#' @param adduct search adduct
#' @param ppm search ppm error
#' @param isotope search isotope
#' @param adduct_rules_table table containing adduct formation rules. Defaults to `adduct_rules()`.
#' @param isotope_rules_table table containing isotope rules. Defaults to `isotope_rules()`.
#' @return 
#' A tibble containing the relevant product search results depending on the method used
#' @examples 
#' metabolite_database <- metaboliteDB(amino_acids)
#' 
#' ## Perform a putative ionisation product search
#' PIPsearch(
#'   metabolite_database,
#'   mz = 133.03358,
#'   adduct = '[M-H]1-',
#'   isotope = '13C')
#' 
#' ## Calculate adduct m/z for a database entry
#' calcAdducts(metabolite_database,1)
#' 
#' @export

setGeneric('calcAdducts',function(db,id,adduct_rules_table = adduct_rules())
  standardGeneric('calcAdducts'))

#' @rdname products

setMethod('calcAdducts',signature = 'MetaboliteDatabase',
          function(db,id,adduct_rules_table = adduct_rules()){
            
            smiles <- db %>%
              descriptors() %>%
              filter(ID == id) %>%
              select(SMILES) %>%
              deframe()
            
            smiles %>%
              ionisationProducts(adduct_rules_table = adduct_rules_table)
          })

#' @rdname products
#' @importFrom dplyr bind_rows select filter
#' @export

setGeneric('PIPsearch',function(db,
                                mz,
                                adduct,
                                ppm = 6,
                                isotope = NA, 
                                adduct_rules_table = adduct_rules(),
                                isotope_rules_table = isotope_rules()) 
  standardGeneric('PIPsearch')
)

#' @rdname products
#' @importFrom rlang expr

setMethod('PIPsearch',signature = 'MetaboliteDatabase',
          function(db,
                   mz,
                   adduct,
                   ppm = 6,
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
            
            if (!is.na(isotope) & nEntries(res) > 0) {
              isoRule <- isotope_rules_table$Rule[isotope_rules_table$Isotope == isotope] %>% 
                parse_expr()
              
              res <- res %>%
                filterER(!!isoRule)
            }
            
            addRule <- adduct_rules_table$Rule[adduct_rules_table$Name == adduct] %>% 
              parse_expr()
            
            res <- res %>%
              filterIP(!!addRule)
            
            res <- res %>%
              {left_join(entries(.),
                         descriptors(.),
                         by = c("ID", "SMILES"))} %>%
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
              descriptors() %>%
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
              descriptors() %>%
              .$MF %>% 
              unique()
            MFs <- MFs %>% 
              bind_rows(.id = 'MF') %>%
              right_join(db %>%
                           descriptors() %>%
                           select(ID,MF), by = "MF") %>%
              select(ID,everything())
            return(MFs)
          })
