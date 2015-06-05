#' Annotation of FIE-HRMS m/z

annotateFIE <-
function(explan_mass,data,peaks,mode,ppm,MZedDB,gencors=T,genmf=T,geniso=T,mzsearch=T,lat=list(c(C=0,iC=0,H=0,iH=0,N=0,iN=0,O=0,iO=0,F=0,Na=0,Si=0,P=0,S=0,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0),c(C=90,iC=1,H=180,iH=0,N=10,iN=0,O=40,F=0,Na=1,Si=0,P=6,S=5,Cl=0,iCl=0,Br=0,iBr=0,K=1,iK=1,iO=1))){
  annot_all <- list()
  if(gencors==T){
    annot_all[["Correlation_Analysis"]] <- corAnalysis(data,rownames(explan_mass),mode=mode)
  }
  if(genmf==T){
    annot_all[["Molecular_Formulas"]] <- apply(explan_mass,1,getMF,mode=mode,lat=lat,ppm=5)
  }
  if(geniso==T){
    if(genmf==F){
      cat("Unable to calculate theoretical isotope distributions. Molecular FOrmulas need to be generated first! Set genmf to TRUE")
    } else {
      annot_all["Theoretical Isotope Distributions"] <- lapply(annot_all$Molecular_Formulas,getIsoDist,mode=mode)
    }
  }
  if(mzsearch==T){
    annot_all["MZedDB Hits"] <- apply(explan_mass)
  }
  
  names(info.annot) <- bin.list
  return(info.annot)
}
