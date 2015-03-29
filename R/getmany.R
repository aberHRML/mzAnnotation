getmany <-
function(host="localhost",port=8888,isrules=TRUE,nohalide=TRUE,precmass=0.7,
 nommass=NULL,accmass=NULL,massacc=1,add="[M+H]1+", sources=NULL,printlong=FALSE,tmpf=TRUE){
tab<-res<-NULL
if(!is.null(nommass))
  lmass=nommass
if(!is.null(accmass))
  lmass=accmass

tmpfile=paste("tmp",floor(abs(rnorm(1)*1000000)),".csv",sep="");
if(tmpf)
 print(paste("Temporary file:",tmpfile))

n<-0;
for(cmass in lmass){
n<-n+1
#for(iadd in add){
#print(paste("Retrieving mass:", cmass, "- adduct:",iadd))
print(paste(n,". Retrieving mass:", cmass))
if(!is.null(nommass))
 res=getone(host,port,isrules,nohalide,precmass,cmass,NULL,massacc,add,sources,printlong=printlong)

if(!is.null(accmass))
 res=getone(host,port,isrules,nohalide,precmass,NULL,cmass,massacc,add,sources,printlong=printlong)

if(!is.null(res)){
  tab=rbind(tab,cbind(rep(cmass, nrow(res)),res))
  if(tmpf)
    write.csv(tab,file=tmpfile)
}


}
tabnames<-c('mz','Formula','adduct','mass','cpdid','smile','smileCanon')

if(printlong)
tabnames<-c('mz','Formula','adduct','mass','cpdid','name',
                 'sourceDB','sourceID','smile','smileCanon')
if (length(tab[,1])>0) {
  colnames(tab)<-tabnames
}
return(tab)
}
