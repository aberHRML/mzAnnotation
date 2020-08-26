#' @importFrom dplyr everything

setMethod('elementFrequencies',signature = 'MetaboliteDatabase',
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