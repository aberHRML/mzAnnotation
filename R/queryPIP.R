
queryPIP <-
function(add,mass_low,mass_high,MZedDB,mz,filter=T){
	rules <- MZedDB$ADDUCT_FORMATION_RULES
  all <- MZedDB$MZedDB_ALL
  metrules <- MZedDB$MZedDB_METRULES
  add <- rules[which(rules[,1]==add),]
  mass_low <- ((mass_low-add[1,"Add"])*add[1,"Charge"])/add[1,"xM"]
  mass_high <- ((mass_high-add[1,"Add"])*add[1,"Charge"])/add[1,"xM"]
  rows <- which(all["Accurate.Mass"]>mass_low & all["Accurate.Mass"] <mass_high)
  metrules.1 <- metrules[rows,]
  Nch <- metrules.1[,"Nch"]
  Nacc <- metrules.1[,"Nacc"]
  Ndon <- metrules.1[,"Ndon"]
  Nnhh <- metrules.1[,"Nnhh"]
  Noh <-  metrules.1[,"Noh"]
  Ncooh <- metrules.1[,"Ncooh"]
  Ncoo <-  metrules.1[,"Ncoo"]
  res <- all[rows[which(eval(parse(text=as.character(add["Rule"]))))],]
  addmz <- sapply(res["Accurate.Mass"],function(x,rules){y <- (x*rules[1,"xM"])/rules[1,"Charge"]+rules[1,"Add"];return(y)},rules=add)
	ppmerr <- sapply(addmz,function(x,mz){y <- (x-mz)/mz*10^6;return(y)},mz=mz)
 	adduct <- rep(add[1,1],nrow(res))
  res <- data.frame(res,Adduct=adduct,Adduct_MZ=addmz,PPMErr=ppmerr)
  res <- res[,-c(3,6,7,8,10)]
  if(filter==T & nrow(res)>0){
   res <- filterPIP(res)
  }
  return(res)
}
