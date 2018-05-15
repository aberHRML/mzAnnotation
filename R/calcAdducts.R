#' calcAdducts
#' @description calculate ionisation products for given database accession
#' @param id accession id.
#' @param db object of class \code{MetaboliteDatabase}.
#' @param adducts table of adduct rules. Defaults to Adducts.
#' @examples 
#' db <- metaboliteDB(aminoAcids,descriptors(aminoAcids))
#' add <- calcAdducts(1,db)
#' @importFrom purrr map_df
#' @export

calcAdducts <- function(id,db,adducts = mzAnnotation::Adducts){
  
  desc <- db@descriptors[[1]] %>%
    filter(ACCESSION_ID == id) %>%
    map_df(rep,times = nrow(adducts))
  
  add <- adducts %>%
    select(Name,Rule) %>%
    bind_cols(desc) %>%
    rowwise() %>%
    mutate(Possible = eval(parse(text = Rule)),`m/z` = calcMZ(Accurate_Mass,Name,adducts = adducts)) %>%
    select(Name,Possible,`m/z`)
  return(add)
}