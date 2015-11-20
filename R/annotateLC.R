
annotateLC <- function(features,data,xset,mode,ppm){
  peakTab <- peakTable(xset)
  mz.id <- round(peakTab$mz, digits=2)
  rt.id <- round(peakTab$rt/60, digits=2)
  f.id <- paste(mode,"M",mz.id,"T",rt.id,sep="")
  feat.ind <- which(f.id %in% features)
  featTab <- data.frame(ID=f.id[feat.ind],peakTab[feat.ind,])
  cat('\nRunning Correlation Analysis\n')
  corrAnalysis <- corAnalysisLC(data,features,mode)
  cat('\nRunning PIP searches\n\n')
  featPIP <- lapply(featTab$mz,getPIP,mode=mode,ppm=ppm)
  names(featPIP) <- featTab$ID
  return(list(featPIP=featPIP,corAnalysis=corrAnalysis,peakTable=data.frame(ID=f.id,peakTab[,1:6])))
}
