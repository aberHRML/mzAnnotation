#' adduct relationship prediction
#' @param mz numeric \code{vector} of accurate m/z
#' @param mode string of either 'p' or 'n; denoting the acquisition mode
#' @param limit limit of deviation for thresholding associations. Defaults to 0.001
#' @author Jasen Finch
#' @export
#' @importFrom utils combn
#' @examples 
#' res <- relationshipPredictor(c(65.51148,132.03023,168.00691),'n')

relationshipPredictor <- function(mz,mode,limit=0.001){
  adducts <- MZedDB$ADDUCT_FORMATION_RULES
  
  if (mode == 'p') {
    adducts <- adducts[adducts$Nelec < 0,]  
  }
  if (mode == 'n') {
    adducts <- adducts[adducts$Nelec > 0,]
  }
  
  M  <- lapply(mz,function(x,add,charge,xM,name){
    m <- ((x - add)*charge)/xM
    names(m) <- name
    return(m)
  },add = adducts$Add,charge = adducts$Charge,xM = adducts$xM,name = adducts$Name)
  
  names(M) <- mz
  
  combinations <- combn(mz,2)
  
  diffs <- apply(combinations, MARGIN = 2,function(x,M,limit){
    add1 <- M[[as.character(x[1])]]
    add2 <- M[[as.character(x[2])]]
    diffs <- lapply(add1,function(x,a){
      sqrt((x - a)^2)
    },a = add2)
    diffs <- as.data.frame(diffs)
    colnames(diffs) <- rownames(diffs)
    res <- data.frame(Adduct1 = colnames(diffs)[which(diffs < limit,arr.ind = T)[,2]],
                      Adduct2 = rownames(diffs)[which(diffs < limit,arr.ind = T)[,1]],
                      stringsAsFactors = F)
    if (nrow(res) > 0) {
      errors <- apply(res,1,function(x,diffs){diffs[x[2],x[1]]},diffs = diffs)
      res <- data.frame(res,Error = errors,stringsAsFactors = F)
    }
    return(res)
  },M = M,limit = limit)
  
  names(diffs) <- apply(combinations,2,function(x){paste(x,collapse = '~')})
  
  diffs <- ldply(diffs,stringsAsFactors = F)
  mzS <- strsplit(diffs$.id,'~')
  mzS <- ldply(mzS,stringsAsFactors = F)
  mzS <- apply(mzS,2,as.numeric)
  colnames(mzS) <- c('mz1','mz2')
  diffs <- diffs[,-1]
  diffs <- data.frame(mzS,diffs,stringsAsFactors = F)
  return(diffs)
}

