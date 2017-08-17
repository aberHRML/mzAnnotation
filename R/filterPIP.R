#' @importFrom stringr str_replace_all coll
#' @importFrom tibble tribble

filterPIP <- function(DB){
  patterns <- tribble(
    ~pattern, ~replacement,
    "ic acid", 'ate',
    "keto", 'oxo',
    "<i>", '',
    "</i>", '',
    "(R)-", '',
    "(S)-", '',
    "-L-", '',
    "L-", '',
    "(+)-", '',
    "cis,", '',
    "cis-", '',
    "trans-", '',
    "&beta;-", '',
    "-D-", '',
    "D-", '',
    "D.", '',
    "&gamma;-", '',
    "-n-", '',
    "N-", '',
    "&alpha;-", '',
    "&alpha;", '',
    "&beta;", ''
  )
  
  for (i in 1:nrow(patterns)) {
    DB$Name <- str_replace_all(DB$Name,coll(patterns$pattern[i]),patterns$replacement[i])
  }
  DB$Name <- tolower(DB$Name)
  DB <- filter(DB,!duplicated(DB$Name))
  return(DB)
}