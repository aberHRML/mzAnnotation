#' Correlation Analysis for mz data matirx
#' @description Calculates correlations with a p-value above a given threshold between all variables in supplied data. M/z relationships between correlated vairables are given for a supplied variable list.
#' @usage corAnalysis(data,varlist,pval=0.05)
#' @param data
#' @param varlist
#' @param pval
#' @author Jasen Finch \email{jsf9@@aber.ac.uk}
#' @seealso corLists, calcDiff


corAnalysis <-
function(data,varlist,mode,pval=0.05){
  cors <- rcorr(as.matrix(data))
  cors$P <- apply(cors$P,1,p.adjust,method="bonferroni")
  cors$r[cors$P>pval] <- 0
  cors$r[cors$r==1] <- 0
  cors$r[cors$r < 0] <- 0
  s <- apply(cors$r,2,sum)
  cors$r <- cors$r[s>0,s>0]
  cors.lists <- corLists(cors$r)
  cors <- NULL
  seq.1 <- seq(1,ncol(cors.lists))
  col <- colnames(cors.lists) %in% varlist
  seq.1 <- seq.1[col]
  seq.2 <- seq.1 + 1
  cors.lists <- cors.lists[,sort(c(seq.1,seq.2))]
  cors.lists[cors.lists==0] <- NA
  add_pred <- NULL
  for(i in 1:(ncol(cors.lists)/2)){
    add_pred[i] <- list(calcDiff(cors.lists[,(i*2-1):(i*2)],mode))    
  }
  names(add_pred) <- colnames(cors.lists)[seq(1,ncol(cors.lists),2)] 
  return(add_pred)
}