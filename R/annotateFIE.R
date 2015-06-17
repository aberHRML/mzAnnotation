#' Annotation of FIE-HRMS m/z
#'
#' @param explan_mass 
#' @param data 
#' @param peaks 
#' @param mode 
#' @param ppm 
#' @param gencors 
#' @param genmf 
#' @param genpip 
#' @param lat 

annotateFIE <-
function(explan_mass,data,peaks,mode,ppm,gencors=T,genmf=T,genpip=T,lat=list(c(C=0,iC=0,H=0,iH=0,N=0,iN=0,O=0,iO=0,F=0,Na=0,Si=0,P=0,S=0,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0),c(C=90,iC=1,H=180,iH=0,N=10,iN=0,O=40,F=0,Na=1,Si=0,P=6,S=5,Cl=0,iCl=0,Br=0,iBr=0,K=1,iK=1,iO=1))){
  annot_all <- list()
  annot_all[["Accurate m/z"]] <- cbind(Bin=rownames(explan_mass),explan_mass)
  if(gencors==T){
    annot_all[["Correlation Analysis"]] <- corAnalysis(data,rownames(explan_mass),mode=mode)
  }
  if(genmf==T){
    annot_all[["Molecular Formulas"]] <- apply(explan_mass,1,getMF,mode=mode,lat=lat,ppm=ppm)
    annot_all["Theoretical Isotope Distributions"] <- list(lapply(annot_all$'Molecular Formulas',getIsoDist,mode=mode))
  }
  if(genpip==T){
    annot_all["Putative Ionisation Products"] <- list(lapply(explan_mass[,1],getPIP,mode=mode,ppm=ppm))
  }
  return(annot_all)
}
