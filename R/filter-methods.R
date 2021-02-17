#' Filter a mass range
#' @rdname filterMR
#' @description Filter a MetaboliteDatabase for a given mass range.
#' @param db S4 object of class MetaboliteDatabase
#' @param lower lower mass boundary
#' @param upper upper mass boundary
#' @examples 
#' db <- metaboliteDB(aminoAcids,descriptors(aminoAcids$SMILES))
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
#' db <- metaboliteDB(aminoAcids,descriptors(aminoAcids$SMILES))
#' db <- filterER(db,'S>0')
#' @export

setMethod('filterER',signature = 'MetaboliteDatabase',
          function(db,rule){
            ef <- elementFrequencies(db)
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
#' rule <- adducts()$Rule[52]
#' db <- metaboliteDB(aminoAcids,descriptors(aminoAcids$SMILES))
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
#' db <- metaboliteDB(aminoAcids,descriptors(aminoAcids$SMILES))
#' db <- filterACCESSIONS(db,c(1,2))
#' @export

setMethod('filterACCESSIONS',signature = 'MetaboliteDatabase',
          function(db,ids){
            db@accessions[[1]] <- db %>%
              getAccessions() %>%
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
#' db <- metaboliteDB(aminoAcids,descriptors(aminoAcids$SMILES))
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