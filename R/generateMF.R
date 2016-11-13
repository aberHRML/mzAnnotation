#' Molecular formula generation
#' @param mz accurate m/z for MF generation
#' @param ppm ppm tolerance for MF generation
#' @param charge charge to apply to MF generation
#' @param applygr \code{boolean} denoting whether to apply the 7 golden rules
#' @param composition numeric \code{vector} of maximum elemental composition
#' @details Uses the HR2 molecular formula generator available at \url{http://maltese.dbs.aber.ac.uk:8888/hrmet/supp/rhrmet.html}.
#' @author Jasen Finch
#' @export
#' @return A \code{data.frame} containing the generated MFs, their theoretical m/z and their PPM error.
#' @examples
#' res <- generateMF(341.10894,ppm = 5,charge = -1, applygr = TRUE,
#'              composition=c(C = 12,iC = 0,H = 22,iH = 0,N = 0,iN = 0,O = 11,iO = 0,F = 0 ,Na = 0,
#'                    Si = 0,P = 0,S = 0,Cl = 0,iCl = 0,Br = 0,iBr = 0,K = 0,iK = 0))

generateMF <- function(mz, ppm = 5, charge = 0, applygr = TRUE, composition = c(C = 12,iC = 0,H = 22,iH = 0,N = 0,iN = 0,O = 11,iO = 0,
                                                         F = 0 ,Na = 0,Si = 0,P = 0,S = 0,Cl = 0,iCl = 0,Br = 0,
                                                         iBr = 0,K = 0,iK = 0)){
    mmu <- ppm/10^6*mz*1000
    res <- mfGen(mz,composition,rep(0,19),mmu,charge,applygr)	
    
    if (length(res) == 0) {
      res <- data.frame(matrix(nrow = 0,ncol = 3))
      colnames(res) <- c("MF", "m/z", "PPM Error")
    } else{
      res <- data.frame(matrix(unlist(res), nrow = length(res),byrow = T),stringsAsFactors = F)
      res <- res[,-c(1,3,5)]
      colnames(res) <- c("MF", "m/z","PPM Error")
      res[,3] <- sapply(as.numeric(as.character(res[,2])),function(x,mass){x <- as.numeric(x);x <- (x - mass)/mass*10^6; return(round(x,5))},mass = mz)
    } 
    return(res)
  }
