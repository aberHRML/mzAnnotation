#' Safe Convert
#' @description safely convert between SMILES and Inchi and to InchiKey using callr to catch potential segfaults
#' @param input a valid SMILE or Inchi
#' @param inputType either "smiles" or "inchi", denoting the input type
#' @param outputType either "smiles", "inchi" or "inchikey", denoting the output type
#' @examples
#' safeConvert(aminoAcids$SMILE[1],'smiles','inchi')
#' @export

safeConvert <- function(input, inputType, outputType) {
  output <- tryCatch(
    callr::r(
      function(x)
        convert(x, inputType, outputType),
      args = list(x),
      error = 'stack'
    ),
    error = function(e)
      return(e$error)
  )
  
  return(output)
}
