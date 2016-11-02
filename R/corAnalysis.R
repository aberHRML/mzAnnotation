#' Correlation Analysis for mz data matirx
#' @description Calculates correlations with a p-value above a given threshold between all variables in supplied data. M/z relationships between correlated vairables are given for a supplied variable list.
#' @usage corAnalysis(data,varlist,pval=0.05)
#' @param data
#' @param varlist
#' @param pval
#' @author Jasen Finch \email{jsf9@@aber.ac.uk}
#' @seealso corLists, calcDiff
#' @export
#' @import Hmisc


corAnalysis <-
function(data,varlist,mode,pval=0.05){
  cors <- rcorr(as.matrix(data))
  cors$P <- apply(cors$P,1,p.adjust,method="bonferroni")
  cors$r[cors$P>pval] <- 0
  cors$r[cors$r==1] <- 0
  cors$r[cors$r < 0] <- 0
  cors <- cors$r
  cors <- cors[,colnames(cors) %in% varlist]
  s <- apply(cors,1,sum)
  cors <- cors[s>0,]
  cors.lists <- corLists(cors)
  cors.lists <- lapply(cors.lists,calcDiff,mo=mode)
  names(cors.lists) <- varlist
  return(cors.lists)
}