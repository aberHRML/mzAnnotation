
annotateLC <- function(features,xset,mode,ppm){
  data("CAMERArules")
  peakTab <- peakTable(xset)
  mz.id <- round(peakTab$mz, digits=2)
  rt.id <- round(peakTab$rt/60, digits=2)
  f.id <- paste(mode,"M",mz.id,"T",rt.id,sep="")
  feat.ind <- which(f.id %in% features)
  featTab <- data.frame(ID=f.id[feat.ind],peakTab[feat.ind,])
  cat('\nRunning PIP searches\n\n')
  featPIP <- lapply(featTab$mz,getPIP,mode=mode,ppm=ppm)
  names(featPIP) <- featTab$ID
  xsa <- xsAnnotate(xset)
  anF <- groupFWHM(xsa, perfwhm = 0.6)
  cat('\n')
  anI <- findIsotopes(anF, mzabs = 0.01)
  cat('\n')
  anIC <- groupCorr(anI, cor_eic_th = 0.75)
  cat('\n')
  anFA <- findAdducts(anIC, polarity="negative",rules = CAMERArules[[mode]])
  cat('\n')
  peaklist <- getPeaklist(anFA)
  peaklist <- data.frame(ID=f.id,peaklist)
  featPC <- unique(peaklist[feat.ind,'pcgroup'])
  peaklist <- peaklist[which(peaklist$pcgroup %in% featPC),]
  return(list(featPIP=featPIP,featRel=peaklist,peakTable=data.frame(ID=f.id,peakTab[,1:6])))
}
