#' @importFrom dplyr everything

elementFrequencies <- function(db){
  MFs <- db@descriptors[[1]]$MF %>%
    unique() %>%
    map(~{
      mf <- .
      mf %>%
        CHNOSZ::makeup() %>%
        as.list() %>%
        as_tibble()
    })
  names(MFs) <- db@descriptors[[1]]$MF %>% unique()
  MFs <- MFs %>% 
    bind_rows(.id = 'MF') %>%
    right_join(db@descriptors[[1]] %>% select(ACCESSION_ID,MF), by = "MF") %>%
    select(ACCESSION_ID,everything())
  return(MFs)
}