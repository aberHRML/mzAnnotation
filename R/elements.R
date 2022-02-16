#' Extract element frequencies
#' @description Extract the element frequencies from a molecular formula.
#' @param MF a vector of molecular formulas
#' @return A tibble containing element frequencies.
#' @examples elementFrequencies(c('H2O','C12H22O11'))
#' @importFrom purrr map_dfr
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

#' Calculate element ratios
#' @description Calculate ratios between all element frequency combinations.
#' @param element_frequencies a tibble containing element frequencies as returned by `elementFrequencies()`
#' @return A tibble containing element frequency ratios.
#' @examples 
#' elementFrequencies(c('H2O','C12H22O11')) %>% 
#'   elementRatios()
#' @importFrom purrr map_dfc
#' @export

elementRatios <- function(element_frequencies,
                          elements = c('C','H','N','O','P','S')){
  elements %>% 
    expand_grid(x = .,y = .) %>% 
    filter(x != y) %>% 
    split(1:nrow(.)) %>% 
    map_dfc(~{
      ratio <- paste0(.x$x,
                      '/',
                      .x$y)
      
      element_ratio <- tibble(
        !!ratio := element_frequencies[[.x$x]] / element_frequencies[[.x$y]]
      )
      
      if (nrow(element_ratio) > 0) return(element_ratio)
      else return(NULL)
    }) %>% 
    bind_cols(element_frequencies %>% 
                select(MF)) %>% 
    select(MF,everything())
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