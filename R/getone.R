getone <-
function(host="localhost",port=8888,isrules=TRUE,nohalide=TRUE,
   precmass=0.7,nommass=NULL,accmass=NULL,massacc=1,add="[M+H]1+",sources=NULL,printlong=FALSE){

url=""
if(printlong)
 url="printlong=TRUE"
if(isrules)
 url=paste(url,"&applyrule=TRUE",sep="")
if(nohalide)
 url=paste(url,"&nohalide=TRUE",sep="")
if(precmass)
 url=paste(url,"&precmass=",precmass,sep="")
if(massacc)
 url=paste(url,"&massacc=",massacc,sep="")

iadd=0
for(addstr in add){
  add0=gsub('\\+','%2b',addstr)
  #url=paste(url,"&adductsel[0]=",add0,sep="")
  url=paste(url,"&adductsel[",iadd,"]=",add0,sep="")
  iadd <- iadd + 1
}

if (!is.null(sources)) {
  for(i in c(1:length(sources)))
   url = paste(url,"&rulesel[",(i-1),"]=",sources[i],sep="")
}

if(!is.null(nommass))
 url=paste(url,"&nommass=",nommass,sep="")
if(!is.null(accmass))
 url=paste(url,"&accmass=",accmass,sep="")

tab<-NULL
tmp=simplePostToHost(host,"/hrmet/search/addsearch.php","",url,port)
if(length(grep('No results',tmp))==0){
	reg=regexpr('tmpres([0-9]+)\\.txt',tmp)
	fi=paste(strsplit(tmp,'')[[1]][reg+1:attr(reg,"match.length")-1],collapse="")
	download.file(paste("http://",host,":",port,"/hrmet/search/tmp/",fi,sep=""),"tmphrmet",quiet=T)
	tmpf=scan("tmphrmet",what="raw",sep="\n")
	tab=do.call("rbind",lapply(tmpf,function(x) strsplit(x,'\t')[[1]]))
}
return(tab)
}
