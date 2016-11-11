#' Correlation Analysis for mz data matirx
#' @description Calculates correlations with a p-value above a given threshold between all variables in supplied data. M/z relationships between correlated vairables are given for a supplied variable list.
#' @param data \code{matrix} containing the spectrally binned data to correlate
#' @param varlist vector of explanatory mz to extract from correlation matrix
#' @param mode either 'p' or 'n' denoting positive or negative acquisition modes respectively
#' @param pval the p value threshold for returning significant correlations
#' @author Jasen Finch \email{jsf9@@aber.ac.uk}
#' @export
#' @importFrom Hmisc rcorr
#' @importFrom stats p.adjust

corAnalysis <-
function(data,varlist,mode,pval=0.05){
  cors <- rcorr(as.matrix(data))
  cors$P <- apply(cors$P,1,p.adjust,method = "bonferroni")
  cors$r[cors$P > pval] <- 0
  cors$r[cors$r == 1] <- 0
  cors$r[cors$r < 0] <- 0
  cors <- cors$r
  cors <- cors[,colnames(cors) %in% varlist]
  s <- apply(cors,1,sum)
  cors <- cors[s > 0,]
  cors.lists <- corLists(cors)
  cors.lists <- lapply(cors.lists,calcDiff,mo = mode)
  names(cors.lists) <- varlist
  return(cors.lists)
}