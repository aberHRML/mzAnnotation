mf_gen <-
function(measured_mass,tolerance,chrge,applygr,mini=c(C=0,iC=0,H=0,iH=0,N=0,iN=0,O=0,iO=0,F=0,Na=0,Si=0,P=0,S=0,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0),maxi=c(C=41,iC=0,H=72,iH=0,N=15,iN=0,O=30,iO=0,F=0,Na=0,Si=0,P=2,S=2,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0)){
  res <- do_calculations(measured_mass,maxi,mini,tolerance,chrge,applygr)
  if (length(res)>0){
    res <- matrix(unlist(res), nrow = length(res),byrow=T)
    colnames(res) <- c("Clean_MF", "MF", "RDB", "m/z", "nadd","niso","Error")
  }
  return(res)
}
