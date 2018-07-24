#' generateMF
#' @description Molecular formula generation
#' @param mass accurate mass
#' @param ppm ppm tolerance
#' @param charge charge
#' @param validation \code{boolean}, apply validation rules
#' @param composition numeric \code{vector} of maximum elemental composition
#' @details Multiple molecular formulas are available:
#' \describe{
#' \item{HR2}{the HR2 generator available at \url{http://maltese.dbs.aber.ac.uk:8888/hrmet/supp/rhrmet.html}.}
#' \item{CDK}{the generator availble in the \code{rcdk} package using \code{\link[rcdk]{generate.formula}}.}
#' \item{Rdisop}{the generator available in the \code{Rdisop} package using \code{\link[Rdisop]{decomposeMass}}.}
#' }
#' @author Jasen Finch
#' @importFrom CHNOSZ makeup
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange
#' @export
#' @return A \code{tibble} containing the generated MFs, their theoretical mass and their PPM error.
#' @examples
#' res <- generateMF(342.11621,
#'                   composition=c(C = 12,H = 22,N = 0,
#'                                 O = 11,P = 0,S = 0))

generateMF <- function(mass, ppm = 1, charge = 0, validation = TRUE, composition = c(C = 12,H = 22,N = 0,O = 11,P = 0,S = 0)){

    comp = c(C = 0,iC = 0,H = 0,iH = 0,N = 0,iN = 0,O = 0,iO = 0,F = 0 ,Na = 0,
             Si = 0,P = 0,S = 0,Cl = 0,iCl = 0,Br = 0,iBr = 0,K = 0,iK = 0)
    comp[names(composition)] <- composition
    
    mmu <- ppm/10^6*mass*1000
    res <- HR2(mass,comp,rep(0,19),mmu,charge,validation)	
    
    if (length(res) == 0) {
      res <- data.frame(matrix(nrow = 0,ncol = 3))
      colnames(res) <- c("MF", "Mass", "PPM Error")
    } else{
      res <- data.frame(matrix(unlist(res), nrow = length(res),byrow = T),stringsAsFactors = F)
      res <- res[,-c(1,3,5)]
      colnames(res) <- c("MF", "Mass","PPM Error")
      res$Mass <- round(as.numeric(res$Mass),5)
      res$`PPM Error` <- sapply(res$Mass,ppmError,measured = mass) %>%
        round(5)
      res$MF <- sapply(res$MF,function(x){
        mf <- makeup(x)
        if (1 %in% mf) {
          mf <- sapply(names(mf),function(y,freq){
            freq <- freq[y]
            if (freq == 1) {
              y
            } else {
              paste(y,freq,sep = '')
            }
          },freq = mf)
          mf <- paste(mf,collapse = '')
          return(mf)
        } else {
          return(x)
        }
      })
    } 
    res <- res %>%
      as_tibble() %>%
      arrange(`PPM Error`)
 
  return(res)
}
