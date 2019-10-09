#' convert
#' @description convert between SMILES and Inchi and to InchiKey
#' @param input a valid SMILE or Inchi
#' @param inputType either "smiles" or "inchi", denoting the input type
#' @param outputType either "smiles", "inchi" or "inchikey", denoting the output type
#' @examples
#' convert(aminoAcids$SMILE[1],'smiles','inchi')
#' @export

convert <- function(input, inputType, outputType) {
  output <- tryCatch(
    callr::r(
      function(x)
        cnvrt(x, inputType, outputType),
      args = list(x),
      error = 'stack'
    ),
    error = function(e)
      return(e$error)
  )
  
  return(output)
}
