#' corAnalysisLC
#' @export

corAnalysisLC <- function(data,varlist,mode,pval=0.05){
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
  cor.names <- sapply(cors.lists,function(x){return(colnames(x)[1])})
  names(cors.lists) <- cor.names
  cors.lists <- lapply(cors.lists,function(x){
    rt <- strsplit(as.character(x[,1]),'T')
    rt <- as.numeric(sapply(rt,function(x){return(x[2])}))
    mz <- strsplit(as.character(x[,1]),'T')
    mz <- sapply(mz,function(x){return(x[1])})
    mz <- as.numeric(gsub(paste(mode,'M',sep=''),'',mz))
    return(data.frame(x,mz=mz,rt=rt))
  })
  cors.diff <- lapply(names(cors.lists),function(x,l,m){
    l <- l[[x]]
    l[,1] <- l$mz
    l <- l[,1:2]
    mz <- strsplit(colnames(l)[1],'T')
    mz <- sapply(mz,function(y){return(y[1])})
    mz <- gsub('M','',mz)
    colnames(l)[1] <- mz
    return(l)
  },l=cors.lists,m=mode)
  cors.diff <- lapply(cors.diff,calcDiff,mo=mode)
  cors.lists <- lapply(1:length(cors.lists),function(x,l,d){
    l <- l[[x]]
    d <- d[[x]]
    return(data.frame(l,d[,3:5]))
  },l=cors.lists,d <- cors.diff)
  names(cors.lists) <- cor.names
  return(cors.lists)
}