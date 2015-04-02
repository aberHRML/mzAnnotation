#' Annotation of FIE-HRMS m/z

annotateFIE <-
function(explan.mass,peak.mat,dat.all,mode,ppm,MZedDB,Path,lat=c(C=90,iC=1,H=180,iH=0,N=10,iN=0,O=40,F=0,Na=1,Si=0,P=6,S=5,Cl=0,iCl=0,Br=0,iBr=0,K=1,iK=1,iO=1)){
  cat("\n","-",mode,sep=""); flush.console()
  explan <- explan.mass[[y]]
  explan.1 <- data.frame(rownames(explan),explan[,1])
  if(mode=="p"){
    dat <- as.matrix(dat.all[[1]])
    peak.mat.1 <- peak.mat[[1]]
    lat <- latmp
  }
  if(mode=="n"){
    dat <- as.matrix(dat.all[[2]])
    peak.mat.1 <- peak.mat[[2]]
    lat <- latmn
    } 
   info.annot <-lapply(bin.list,MFIsoAgg,explan=explan.1,peak.mat=peak.mat.1,mode.1=mode,ppm.1=ppm,add_pred.1=add_pred,MZedDB.1=MZedDB,lat.1=lat,srce.1=srce)  
  names(info.annot) <- bin.list
  cat("\n",'...  done in ',FIEmspro:::timer_end(time1)$dt,sep="")
  return(info.annot)
}
