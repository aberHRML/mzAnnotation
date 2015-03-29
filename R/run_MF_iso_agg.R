run_MF_iso_agg <-
function(x,explan,peak.mat,mode.1,ppm.1,add_pred.1,MZedDB.1=MZedDB,lat.1=lat,srce.1=srce){
  source(srce.1)
  if (x %in% names(add_pred.1)){
    annot <- MF_iso_agg(explan[x,],peak.mat,mode.1,ppm.1,add_pred.1[as.character(x)],MZedDB.1,lmax=lat.1)
  } else {
    annot <- MF_iso_agg(explan[x,],peak.mat,mode.1,ppm.1,"No Correlations",MZedDB.1,lmax=lat.1)
  }
  return(annot)
}
