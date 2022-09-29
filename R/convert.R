#' Convert chemical notation
#' @description convert between SMILES and Inchi and to InchiKey
#' @param input a valid SMILE or Inchi
#' @param inputType either "smiles" or "inchi", denoting the input type
#' @param outputType either "smiles", "inchi" or "inchikey", denoting the output type
#' @examples
#' convert(amino_acids$SMILES[1],'smiles','inchi')
#' @importFrom callr r
#' @importFrom utils getFromNamespace
#' @export

convert <- function(input, inputType, outputType) {
  output <- r(function(input,inputType,outputType){
    cnvrt <- getFromNamespace('cnvrt','mzAnnotation')
    cnvrt(input,inputType,outputType)
  },
      args = list(input = input,inputType = inputType,outputType = outputType),
      error = 'stack'
    )
  
  if (output == 'NA') output <- NA
  
  return(output)
}
