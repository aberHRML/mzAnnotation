#' Extract element frequencies
#' @description Extract the element frequencies from a molecular formula.
#' @param MF a vector of molecular formulas
#' @return A tibble containing element frequencies.
#' @examples elementFrequencies(c('H2O','C12H22O11'))
#' @export

elementFrequencies <- function(MF){
  MF %>% 
    map_dfr(~{
      .x %>% 
        count.elements() %>% 
        tibble(MF = .x,
               element = names(.),
               frequency = .) %>% 
        spread(element,
               frequency)
    }) 
}


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