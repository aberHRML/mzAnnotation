#' applyIonisationRule
#' @description apply an adduct formation rule to a database entry
#' @param ID database ID of the metabolite  forr which to apply formation rule
#' @param adduct adduct of which to apply rule
#' @param adducts adduct table containing available adduct rules
#' @param DB Database to use. Defaults to \code{MZedDB}
#' @examples 
#' #applyIonisationRule('D27880',"[M-H]1-")
#' @export

applyIonisationRule <- function(ID, adduct, adducts = mzAnnotation::Adducts, DB = mzAnnotation::MZedDB) {
  metaboliteID <- ID
  metaboliteRule <- filter(DB$Rules,ID == metaboliteID)
  
  Nch <- unlist(metaboliteRule[1,"Charge"],use.names = F)
  Nacc <- unlist(metaboliteRule[1,"Nacc"],use.names = F)
  Ndon <- unlist(metaboliteRule[1,"Ndon"],use.names = F)
  Nnhh <- unlist(metaboliteRule[1,"Nnhh"],use.names = F)
  Noh <-  unlist(metaboliteRule[1,"Noh"],use.names = F)
  Ncooh <- unlist(metaboliteRule[1,"Ncooh"],use.names = F)
  Ncoo <-  unlist(metaboliteRule[1,"Ncoo"],use.names = F)
  
  adductRule <- filter(adducts, Name == adduct)
  eval(parse(text = adductRule$Rule))
}