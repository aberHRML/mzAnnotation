#' m/z relationship prediction
#' @description adduct, isotope and biotransfromation prediction.
#' @param mz numeric \code{vector} of accurate m/z
#' @param mode string of either 'p' or 'n; denoting the acquisition mode
#' @param limit limit of deviation for thresholding associations. Defaults to 0.001
#' @param add character vector of adducts to use. If \code{NULL} all available adducts will be used
#' @author Jasen Finch
#' @export
#' @importFrom utils combn
#' @examples 
#' res <- relationshipPredictor(c(132.03023,133.01425,133.03359,168.00691),'n')

relationshipPredictor <- function(mz,mode,limit=0.001, add = NULL){
  options(digits = 15)
  adducts <- MZedDB$ADDUCT_FORMATION_RULES
  isotopes <- MZedDB$ISOTOPE_RULES
  isotopes <- isotopes[!(isotopes$Isotope %in% c('Cl37','K41')),]
  isotopes <- data.frame(Isotope = c(NA,isotopes$Isotope),Difference = c(0,isotopes$Mass.Difference),stringsAsFactors = F)
  transformations <- MZedDB$BIOTRANSFORMATION_RULES
  transformations <- data.frame(Transformation = c(NA,transformations$MF.Change),Difference = c(0,transformations$Difference),stringsAsFactors = F)
  if (length(mode) == 1) {
    if (mode == 'p') {
      adducts <- adducts[adducts$Nelec < 0,]  
    }
    if (mode == 'n') {
      adducts <- adducts[adducts$Nelec > 0,]
    } 
  } else {
    adducts <- adducts[adducts$Nelec > 0 | adducts$Nelec < 0,]
  }
  
  if (!is.null(add)) {
    adducts <- adducts[adducts$Name %in% add,]
    if (nrow(adducts) == 0) {
      stop('no adducts selected')
    }
  }
  
  M  <- lapply(mz,function(x,add,charge,xM,name){
    m <- ((x - add)*charge)/xM
    names(m) <- name
    return(m)
  },add = adducts$Add,charge = adducts$Charge,xM = adducts$xM,name = adducts$Name)
  
  names(M) <- mz
  
  combinations <- combn(mz,2)
  
  diffs <- apply(combinations, MARGIN = 2,function(x,M,limit,isotopes,transformations){
    add1 <- M[[as.character(x[1])]]
    add2 <- M[[as.character(x[2])]]
    diffs <- lapply(add1,function(x,a){
      abs(x - a)
    },a = add2)
    diffs <- as.data.frame(diffs)
    colnames(diffs) <- rownames(diffs)
    res <- apply(isotopes,1,function(y,diffs,limit,transformations){
      diffs <- abs(diffs - as.numeric(y[2]))
      res <- apply(transformations,1,function(z,diffs,limit){
        diffs <- abs(diffs - abs(as.numeric(z[2])))
        res <- data.frame(Adduct1 = colnames(diffs)[which(diffs < limit,arr.ind = T)[,2]],
                          Adduct2 = rownames(diffs)[which(diffs < limit,arr.ind = T)[,1]],
                          stringsAsFactors = F)
        if (nrow(res) > 0) {
          errors <- apply(res,1,function(a,diffs){diffs[a[2],a[1]]},diffs = diffs)
          res <- data.frame(res,Error = errors,stringsAsFactors = F)
        }
        return(res)
      },diffs = diffs,limit = limit)
      names(res) <- transformations$Transformation
      res <- ldply(res,.id = 'Transformation')
      return(res)
    },diffs = diffs,limit = limit,transformations = transformations)
    names(res) <- isotopes$Isotope
    res <- ldply(res,.id = 'Isotope')
    return(res)
  },M = M,limit = limit,isotopes = isotopes,transformations = transformations)
  
  names(diffs) <- apply(combinations,2,function(x){paste(x,collapse = '~')})
  
  diffs <- ldply(diffs,stringsAsFactors = F)
  mzS <- strsplit(diffs$.id,'~')
  mzS <- ldply(mzS,stringsAsFactors = F)
  mzS <- data.frame(matrix(apply(mzS,2,as.numeric),ncol = 2))
  colnames(mzS) <- c('mz1','mz2')
  diffs <- diffs[,-1]
  diffs <- data.frame(mzS,diffs,stringsAsFactors = F)
  if (nrow(diffs) > 0) {
    adducts <- MZedDB$ADDUCT_FORMATION_RULES
    adducts[,c(2,3,4,5)] <- apply(adducts[,c(2,3,4,5)],2,as.numeric)
    diffs <- apply(diffs,1,function(x,adducts,isotopes,trans){
      mz1 <- ((as.numeric(x[1]) - adducts$xM[which(adducts$Name == x[5])]) * adducts$xM[which(adducts$Name == x[5])])/adducts$xM[which(adducts$Name == x[5])]
      mz2 <- ((as.numeric(x[2]) - adducts$xM[which(adducts$Name == x[6])]) * adducts$xM[which(adducts$Name == x[6])])/adducts$xM[which(adducts$Name == x[6])]
      Isotope <- rep(NA,2)
      Transformation <- rep(NA,2)
      if (!is.na(x[3])) {
        i <- isotopes$Mass.Difference[which(isotopes$Isotope == x[3])]
        if (!is.na(x[4])) {
          t <- trans$Difference[which(trans$MF.Change == x[4])]
          iso <- c(c(abs((mz1 - t) - mz2),abs(mz1 - (mz2 - t))) - i)
          Isotope[which(iso == min(iso))] <- x[3]
        } else {
          Isotope[which(c(mz1,mz2) == max(c(mz1,mz2)))] <- x[3]
        }
      }
      if (!is.na(x[4])) {
        t <- trans$Difference[which(trans$MF.Change == x[4])]
        if (!is.na(x[3])) {
          i <- isotopes$Mass.Difference[which(isotopes$Isotope == x[3])]
          tran <- c(c(abs((mz1 - i) - mz2),abs(mz1 - (mz2 - i))) - t)
          Transformation[which(tran == min(tran))] <- x[4]
        } else {
          if (t < 0) {
            Transformation[which(c(mz1,mz2) == min(c(mz1,mz2)))] <- x[4]
          } else {
            Transformation[which(c(mz1,mz2) == max(c(mz1,mz2)))] <- x[4]
          }
        }
      }
      return(c(x[1:2], Isotope1 = Isotope[1],Isotope2 = Isotope[2],Transformation1 = Transformation[1],Transformation2 = Transformation[2],x[5:7]))
    },adducts = adducts,isotopes = MZedDB$ISOTOPE_RULES, trans = MZedDB$BIOTRANSFORMATION_RULES)
    diffs <- data.frame(t(diffs),stringsAsFactors = F)
    diffs[,1:2] <- apply(diffs[,1:2],2,as.numeric)
  } else {
    diffs <- data.frame(diffs[,1:2],Isotope1 = character(),Isotope2 = character(),Transformation1 = character(),Transformation2 = character(),diffs[,5:6])
  }
  return(diffs)
}

