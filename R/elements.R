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
#' @param elements a character vector of elements for which to calculate ratios
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
    filter(.data$x != .data$y) %>% 
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
