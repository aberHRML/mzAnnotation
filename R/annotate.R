annotate <-
function(explan.mass,peak.mat,dat.all,ppm,MZedDB,Path,DF,dn,srce,par.mode=T,nCores=1,j.mem = "10g",plots=T,latmp=NULL,latmn=NULL){
  if(!file.exists(paste(Path,DF,paste(DF,"Annotations",sep="_"),sep="/"))){
    dir.create(paste(Path,DF,paste(DF,"Annotations",sep="_"),sep="/"))
  }
  options(java.parameters= paste("\"-Xmx",j.mem,"\"",sep=""))
  time1 <- FIEmspro:::timer_start()
  suppressPackageStartupMessages(library(xlsx))
  wb <- createWorkbook()
  gc()
  info.annot.1 <- NULL
  modes <- c("pos","neg")
  for (y in 1:length(explan.mass)){
    cat("\n","-",dn[y],sep=""); flush.console()
    explan <- explan.mass[[y]]
    explan.1 <- data.frame(rownames(explan),explan[,1])
    mode <- modes[y]
    if(mode=="P"|mode=="p"|mode=="pos"){
      mode <- "p"
      dat <- as.matrix(dat.all[[1]])
      peak.mat.1 <- peak.mat[[1]]
      lat <- latmp
    }
    if(mode=="N"|mode=="n"|mode=="neg"){
      mode <- "n"
      dat <- as.matrix(dat.all[[2]])
      peak.mat.1 <- peak.mat[[2]]
      lat <- latmn
    }
    ## Calculate correlations
    cat("\n","Calculating Correlations",sep="");flush.console()
    suppressPackageStartupMessages(library(Hmisc))
    cors <- rcorr(dat)
    cors$P <- apply(cors$P,1,p.adjust,method="bonferroni")
    cors$r[cors$P>0.05] <- 0
    cors$r[cors$r==1] <- 0
    cors$r[cors$r < 0] <- 0
    s <- apply(cors$r,2,sum)
    cors$r <- cors$r[s>0,s>0]
    cors.lists <- cor.lists(cors$r)
    cors <- NULL
    gc()
    seq.1 <- seq(1,ncol(cors.lists))
    bin.list <- as.character(explan.1[,1])
    col <- colnames(cors.lists) %in% bin.list
    seq.1 <- seq.1[col]
    seq.2 <- seq.1 + 1
    cors.lists <- cors.lists[,sort(c(seq.1,seq.2))]
    cors.lists[cors.lists==0] <- NA
    add_pred <- NULL
    for(i in 1:(ncol(cors.lists)/2)){
      add_pred[i] <- list(calc_diff(cors.lists[,(i*2-1):(i*2)],mode))  	
    }
    names(add_pred) <- colnames(cors.lists)[seq(1,ncol(cors.lists),2)] 
    cat("\n","Collecting MFs, Isotope Distributions and MZedDB hits","\n",sep="");flush.console()
    if(par.mode==T){
      library(parallel)
      clust = makeCluster(nCores, type="PSOCK")
      info.annot <- clusterApplyLB(clust, bin.list, fun=run_MF_iso_agg,explan=explan.1,peak.mat=peak.mat.1,mode.1=mode,ppm.1=ppm,add_pred.1=add_pred,MZedDB.1=MZedDB,lat.1=lat,srce.1=srce)
      stopCluster(clust)
    } else {
      info.annot <-lapply(bin.list,run_MF_iso_agg,explan=explan.1,peak.mat=peak.mat.1,mode.1=mode,ppm.1=ppm,add_pred.1=add_pred,MZedDB.1=MZedDB,lat.1=lat,srce.1=srce)  
    }
    names(info.annot) <- bin.list
    annot_all <- annot_sheet.1(info.annot)
    wb <- xlsx.annot(wb,annot_all,dn[y],Path,DF,plots=plots)
    #write.table(annot_all,file=paste(Path,paste(DF,dn[y],"annotation_sheet.csv",sep="_"),sep="/"),sep=",",row.names=F,col.names=F)
    info.annot.1[y] <- list(info.annot)
  }
  names(info.annot.1) <- dn
  save(info.annot.1,file=paste(Path,DF,paste(DF,"Annotations",sep="_"),paste(DF,"annotation.RData",sep="_"),sep="/"))
  saveWorkbook(wb,paste(Path,DF,paste(DF,"Annotations",sep="_"),paste(DF,"annotation_sheet.xlsx",sep="_"),sep="/"))
  cat("\n",'...  done in ',FIEmspro:::timer_end(time1)$dt,sep="")
  cat("\n",'...  done in ',FIEmspro:::timer_end(time1)$dt,sep="",file=paste(Path,DF,paste(DF,"log-file.txt",sep="_"),sep="/"),append=T)
  return(info.annot.1)
}
