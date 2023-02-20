#' Adduct formation rules
#' @rdname adducts
#' @description Formation rules for electrospray ionisation adducts
#' @format A tibble containing 58 rows and 6 columns
#' \describe{
#' \item{Name}{Adduct name}
#' \item{Charge}{Charge added}
#' \item{xM}{Number of M needed}
#' \item{Add}{Mass value change}
#' \item{Nelec}{Electron change}
#' \item{AddAt}{pseudo formula of the atoms to be added to the molecular formula of one M}
#' \item{RemAt}{pseudo formula of the atoms to be removed to the molecular formula of one M}
#' \item{AddEx}{pseudo formula of the atoms to be added to the final ionisation product molecular formula}
#' \item{RemEx}{pseudo formula of the atoms to be removed to the final ionisation product molecular formula}
#' \item{Rule}{Structural rule for formation}
#' }
#' @examples 
#' adduct_names()
#' 
#' adduct_rules()
#' @export

adduct_rules <- function(){
  return(Adducts)
}

#' @rdname adducts
#' @export

adduct_names <- function(){
  return(Adducts$Name)
}

#' Isotopic rules
#' @rdname isotopes
#' @description Isotope rules
#' @format A tibble containing 8 rows and 3 columns.
#' \describe{
#' \item{Isotope}{Isotope name}
#' \item{Mass Difference}{Isotopic mass difference}
#' \item{Rule}{Elemental composition rules}
#' }
#' @examples 
#' isotope_names()
#' 
#' isotope_rules()
#' @export

isotope_rules <- function(){
  return(Isotopes)
}

#' @rdname isotopes
#' @export

isotope_names <- function(){
  return(Isotopes$Isotope)
}

#' Transformation rules
#' @rdname transformations
#' @description Transformation rules
#' @format A tibble containing 12 rows and 9 columns.
#' \describe{
#' \item{Name}{Transformation name}
#' \item{MF Change}{Molecular formula change}
#' \item{Difference}{Mass difference}
#' \item{C}{Change in carbon atoms}
#' \item{H}{Change in hydrogen atoms}
#' \item{O}{Change in oxygen atoms}
#' \item{N}{Change in nitrogen atoms}
#' \item{P}{Change in phosphorus atoms}
#' \item{S}{Change in sulphur atoms}
#' }
#' @examples 
#' transformation_names()
#' 
#' transformation_rules()
#' @export

transformation_rules <- function(){
  return(Transformations)
}

#' @rdname transformations
#' @export

transformation_names <- function(){
  return(Transformations$`MF Change`)
}