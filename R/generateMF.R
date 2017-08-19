#' Molecular formula generation
#' @param mass accurate mass for MF generation
#' @param ppm ppm tolerance for MF generation
#' @param charge charge to apply to MF generation
#' @param validation \code{boolean} apply validation rules
#' @param composition numeric \code{vector} of maximum elemental composition
#' @details Uses the HR2 molecular formula generator available at \url{http://maltese.dbs.aber.ac.uk:8888/hrmet/supp/rhrmet.html}.
#' @author Jasen Finch
#' @importFrom CHNOSZ makeup
#' @importFrom tibble as_tibble
#' @importFrom rcdk generate.formula
#' @importFrom dplyr arrange
#' @importFrom Rdisop decomposeMass initializeCHNOPS
#' @export
#' @return A \code{data.frame} containing the generated MFs, their theoretical mass and their PPM error.
#' @examples
#' res <- generateMF(342.11621,
#'                   composition=c(C = 12,H = 22,N = 0,
#'                                 O = 11,P = 0,S = 0))

generateMF <- function(mass, ppm = 1, charge = 0, validation = TRUE, composition = c(C = 12,H = 22,N = 0,O = 11,P = 0,S = 0), generator = 'HR2'){
  generators <- c('HR2','CDK','Rdisop')
  
  if (!(generator %in% generators)) {
    stop('Unrecognised MF generator')
  }
  
  if (generator == 'HR2') {
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
  }
  
  if (generator == 'CDK') {
    amu <- ppm/10^6*mass
    composition <- composition[composition != 0]
    composition = map(1:length(composition),~{
      c <- composition[.]
      e <- names(c)
      names(c) <- NULL
      return(c(e,0,c))
    })
    res <- suppressWarnings(try(generate.formula(mass,window = amu,validation = validation,elements = composition,charge = charge),silent = T))
    if (is.list(res)) {
      MF <- sapply(res,function(x){return(x@string)})
      m <- sapply(res,function(x){return(x@mass)}) %>%
        round(5)
      Error <- sapply(m,ppmError,measured = mass) %>% 
        round(5)
      res <- tibble(MF = MF, Mass = m, `PPM Error` = Error)
    } else {
      res <- tibble(MF = character(), Mass = numeric(), `PPM Error` = numeric())
    }
  }
  
  if (generator == 'Rdisop') {
    composition <- composition[composition != 0]
    minElements <- str_c(str_c(names(composition),0),collapse = '')
    maxElements <- str_c(str_c(names(composition),composition),collapse = '')
    elements <- initializeCHNOPS()
    elements <- map(elements,~{
      if(.$name %in% names(composition)) {
        return(.)
      }
    })
    elements <- elements[!sapply(elements,is.null)]
    res <- decomposeMass(mass,ppm = ppm,filter = validation,elements = elements,
                         minElements = minElements,maxElements = maxElements,mzabs = 0)
    if (length(res) > 0) {
      MF <- res$formula
      m <- res$exactmass %>%
        round(5)
      Error <- sapply(m,ppmError,measured = mass) %>% 
        round(5)
      res <- tibble(MF = MF, Mass = m, `PPM Error` = Error)
    } else {
      res <- tibble(MF = character(), Mass = numeric(), `PPM Error` = numeric())
    }
   
  }
  
  return(res)
}
