iso.distr <-
function(mf, chrg=0, limri=0.00000001, limitfin=0.009, rv=10000){
## -------------------------------------------------------------------------

mzv=cbind("1",	"12C",	"C",	"100",	"0.9889",	"12")
mzv=rbind(mzv,cbind("2",	"13C",	"iC",	"1.12245929821013",	"0.0111",	"13.003354"))
mzv=rbind(mzv,cbind("3",	"1H",	"H",	"100",	"0.9999",	"1.007825"))
mzv=rbind(mzv,cbind("4",	"2H",	"iH",	"0.01000100010001",	"0.0001",	"2.014102"))
mzv=rbind(mzv,cbind("5",	"79Br",	"Br",	"100",	"0.5054",	"78.918348"))
mzv=rbind(mzv,cbind("6",	"81Br",	"iBr",	"97.8630787495053",	"0.4946",	"80.916344"))
mzv=rbind(mzv,cbind("7",	"35Cl",	"Cl",	"100",	"0.7553",	"34.968854"))
mzv=rbind(mzv,cbind("8",	"37Cl",	"iCl",	"32.3977227591685",	"0.2447",	"36.965896"))
mzv=rbind(mzv,cbind("9",	"19F",	"F",	"100",	"1",	"18.998405"))
mzv=rbind(mzv,cbind("10",	"127I",	"I",	"100",	"1",	"126.904352"))
mzv=rbind(mzv,cbind("11",	"39K",	"K",	"100",	"0.93258",	"38.9637069"))
mzv=rbind(mzv,cbind("12",	"41K",	"iK",	"7.21653906367282",	"0.0673",	"40.9618259"))
mzv=rbind(mzv,cbind("13",	"40K",	"i2K",	"0.12545840571318",	"0.00117",	"39.9639986"))
mzv=rbind(mzv,cbind("14",	"24Mg",	"Mg",	"100",	"0.7899",	"23.9850419"))
mzv=rbind(mzv,cbind("15",	"26Mg",	"iMg",	"13.9384732244588",	"0.1101",	"25.98259304"))
mzv=rbind(mzv,cbind("16",	"25Mg",	"i2Mg",	"12.6598303582732",	"0.1",	"24.98583702"))
mzv=rbind(mzv,cbind("17",	"14N",	"N",	"100",	"0.9963",	"14.003074"))
mzv=rbind(mzv,cbind("18",	"15N",	"iN",	"0.37137408411121",	"0.0037",	"15.000108"))
mzv=rbind(mzv,cbind("19",	"23Na",	"Na",	"100",	"1",	"22.989773"))
mzv=rbind(mzv,cbind("20",	"16O",	"O",	"100",	"0.9976",	"15.994915"))
mzv=rbind(mzv,cbind("21",	"18O",	"iO",	"0.20048115477145",	"0.002",	"17.99916"))
mzv=rbind(mzv,cbind("22",	"17O",	"i2O",	"0.04009623095429",	"0.0004",	"16.999133"))
mzv=rbind(mzv,cbind("23",	"31P",	"P",	"100",	"1",	"30.973763"))
mzv=rbind(mzv,cbind("24",	"32S",	"S",	"100",	"0.9502",	"31.972074"))
mzv=rbind(mzv,cbind("25",	"34S",	"iS",	"4.44117027994107",	"0.0422",	"33.967864"))
mzv=rbind(mzv,cbind("26",	"33S",	"i2S",	"0.79983161439697",	"0.0076",	"32.97146"))
mzv=rbind(mzv,cbind("27",	"36S",	"i3S",	"0.01052410018943",	"0.0001",	"35.967091"))
mzv=rbind(mzv,cbind("28",	"28Si",	"Si",	"100",	"0.9221",	"27.976927"))
mzv=rbind(mzv,cbind("29",	"29Si",	"iSi",	"5.09706105628457",	"0.047",	"28.976491"))
mzv=rbind(mzv,cbind("30",	"30Si",	"i2Si",	"3.35104652423815",	"0.0309",	"29.973761"))
mzv=rbind(mzv,cbind("31",	"e",	"e",	"100",	"1",	"0.0005484"))

mzv <- mzv[,-1]
colnames(mzv)=c("ICP.Isotope",	"Isotope",	"PercAbundance",	"Abundance",	"MW")
mzv <- data.frame(mzv)
mzv[,3] = as.numeric(as.character(mzv[,3]))
mzv[,4] = as.numeric(as.character(mzv[,4]))
mzv[,5] = as.numeric(as.character(mzv[,5]))
## -------------------------------------------------------------------------
exptn=c("O","K","Si","Mg","S")                  
# =============================
# R2.6:
#greg <- gregexpr("[[:upper:]]{1}[[:lower:]]?[[:digit:]]*", mf, perl = TRUE)
#ugreg <- c(unlist(greg), nchar(mf)+1)
#c <- sapply(1:(length(ugreg)-1), function(z) substr(mf, ugreg[z], ugreg[z+1]-1) )
#at <- sub("([0-9]*)$", "", sapply(1:(length(ugreg)-1),   # blank out the digits
#                 function(z) substr(mf, ugreg[z], ugreg[z+1]-1) ) )
#nu <- sub("^$", "1", sub("([a-z]*)", "",    # just [a-z] works,  see above
#                     sapply(1:(length(ugreg)-1),
#                           function(z) substr(mf, ugreg[z], ugreg[z+1]-1) ) ) )
# =============================
# R3.0:
greg <- gregexpr("[A-Z]{1}[a-z]?[0-9]*", mf)  # lower case does not work
ugreg <- c(unlist(greg), nchar(mf)+1)
c <- sapply(1:(length(ugreg)-1), function(z) substr(mf, ugreg[z], ugreg[z+1]-1) )
at <- sub("(\\d*)$", "", sapply(1:(length(ugreg)-1),   # blank out the digits
                 function(z) substr(mf, ugreg[z], ugreg[z+1]-1) ) )
nu <- sub("^$", "1", sub("([A-Za-z]*)", "",    # subst "1" for empty strings
                     sapply(1:(length(ugreg)-1),
                           function(z) substr(mf, ugreg[z], ugreg[z+1]-1) ) ) )
# =============================
nu = as.numeric(nu)
at2=at
at2[] <- "0"
at3=at2
at4=at2
mw1=nu
mw1[] <- 0
mw2=mw1
mw3=mw1
mw4=mw1
dmw2=mw1
dmw3=mw1
dmw4=mw1
ap=mw1
ap2=mw1
ap3=mw1
ap4=mw1
for (i in 1:length(at)) {
  tmp <- which(mzv[,2]==at[i])
  mw1[i] <- as.numeric(mzv[tmp,5])
  mw2[i] <- as.numeric(mzv[tmp+1,5])
  dmw2[i] <- as.numeric(mzv[tmp+1,5])-as.numeric(mzv[tmp,5])
  at2[i] <- as.character(mzv[tmp+1,1])
  if (at[i]=="O") {
    mw3[i] <- as.numeric(mzv[tmp+2,5])
    dmw3[i] <- as.numeric(mzv[tmp+2,5])-as.numeric(mzv[tmp,5])
    at3[i] <- as.character(mzv[tmp+2,1])
    ap3[i] <- mzv[tmp+2,3]
  }
  if (at[i]=="S") {
    mw3[i] <- as.numeric(mzv[tmp+2,5])
    dmw3[i] <- as.numeric(mzv[tmp+2,5])-as.numeric(mzv[tmp,5])
    at3[i] <- as.character(mzv[tmp+2,1])
    ap3[i] <- mzv[tmp+2,3]
    # 36S
    mw4[i] <- as.numeric(mzv[tmp+3,5])
    dmw4[i] <- as.numeric(mzv[tmp+3,5])-as.numeric(mzv[tmp,5])
    at4[i] <- as.character(mzv[tmp+3,1])
    ap4[i] <- mzv[tmp+3,3]
  }
  if (at[i]=="K") {
    mw3[i] <- as.numeric(mzv[tmp+2,5])
    dmw3[i] <- as.numeric(mzv[tmp+2,5])-as.numeric(mzv[tmp,5])
    at3[i] <- as.character(mzv[tmp+2,1])
    ap3[i] <- mzv[tmp+2,3]
  }
  if (at[i]=="Si") {
    mw3[i] <- as.numeric(mzv[tmp+2,5])
    dmw3[i] <- as.numeric(mzv[tmp+2,5])-as.numeric(mzv[tmp,5])
    at3[i] <- as.character(mzv[tmp+2,1])
    ap3[i] <- mzv[tmp+2,3]
  }
  ap[i] <-  mzv[tmp,3]
  ap2[i] <- mzv[tmp+1,3]
}
if (length(grep("F",at)>0)) {
   ap2[grep("F",at)]=0
   mw2[grep("F",at)]=0
   dmw2[grep("F",at)]=0
   at2[grep("F",at)]="0"
}
if (length(grep("I",at)>0)) {
   ap2[grep("I",at)]=0
   mw2[grep("I",at)]=0
   dmw2[grep("I",at)]=0
   at2[grep("I",at)]="0"
}
if (length(grep("Na",at)>0)) {
   ap2[grep("Na",at)]=0
   mw2[grep("Na",at)]=0
   dmw2[grep("Na",at)]=0
   at2[grep("Na",at)]="0"
}
if (length(grep("P",at)>0)) {
   ap2[grep("P",at)]=0
   mw2[grep("P",at)]=0
   dmw2[grep("P",at)]=0
   at2[grep("P",at)]="0"
}
# =============================
res <- 1:(rv*4) ; dim(res) <- c(rv,4)
res[] <- "0"
colnames(res) <- c("m/z", "Prob", "rel Int", "Isotope")
# =============================
#------------------------------
# mono isotopic mz: (0)
i=1
merki=0
elmt=1
loop=0
res[i,4] <- as.matrix(mzv[i,1])
# rel Int
res[i,3] <- ap[elmt] #100   
# m/z
for (j in 1:length(at)) {   
  res[i,1] <- as.numeric(res[i,1]) + as.numeric(nu[j]) * as.numeric(mw1[j])
}
moiso = as.numeric(res[i,1])
#------------------------------
# 13C mz and single: (1)
for (elmt in 1:length(at)) {
  for (loop in 1: min(nu[elmt],10)) {
    if (nu[elmt]-loop>=0) {
      # rel Int
      i=i+1
      if (loop==1) {        
        merki = i
        res[i,3] <- nu[elmt] * ap2[elmt]
      } else {
        res[i,3] <- as.numeric(res[merki+loop-2,3])/ap[elmt]*ap2[elmt]*(nu[elmt]-loop+1)/loop
      }
      if (as.numeric(res[i,3])<limri) {
        res[i,3] <- 0
        i=i-1
      } else {
        # m/z
        res[i,1] <- moiso + dmw2[elmt]*loop      
        # isotope
        res[i,4] <- as.matrix(paste(at2[elmt],loop))  
      }
    }
  }
  if (length(grep(at[elmt],exptn))>0) {  # extra for exeptions 
    for (loop in 1: min(nu[elmt],6)) {
      if (nu[elmt]-loop>=0) {
        # rel Int
        i=i+1
        if (loop==1) {     
          merki = i
          res[i,3] <- nu[elmt] * ap3[elmt]
        } else {
          res[i,3] <- as.numeric(res[merki+loop-2,3])/ap[elmt]*ap3[elmt]*(nu[elmt]-loop+1)/loop
        }
        if (as.numeric(res[i,3])<limri) {
          res[i,3] <- 0
          i=i-1
        } else {
          # m/z
          res[i,1] <- moiso + dmw3[elmt]*loop     
          # isotope
          res[i,4] <- as.matrix(paste(at3[elmt],loop))  
        }
      }
    }
    if (at[elmt]=="S") {         # extra for 36S only
      for (loop in 1: min(nu[elmt],6)) {
        if (nu[elmt]-loop>=0) {
          # rel Int
          i=i+1
          if (loop==1) {  
            merki = i
            res[i,3] <- nu[elmt] * ap4[elmt]
          } else {
            res[i,3] <- as.numeric(res[merki+loop-2,3])/ap[elmt]*ap4[elmt]*(nu[elmt]-loop+1)/loop
          }
          if (as.numeric(res[i,3])<limri) {
            res[i,3] <- 0
            i=i-1
          } else {
            # m/z
            res[i,1] <- moiso + dmw4[elmt]*loop   
            # isotope
            res[i,4] <- as.matrix(paste(at4[elmt],loop))  
          }
        }
      }
    }
  }
}
if (length(grep("O",at))>0) {
  dmw2 <- c(dmw2,dmw3[grep("O",at)])
  at2 <- c(at2,at3[grep("O",at)])
}
if (length(grep("K",at))>0) {
  dmw2 <- c(dmw2,dmw3[grep("K",at)])
  at2 <- c(at2,at3[grep("K",at)])
}
if (length(grep("S",at))>0) {
  dmw2 <- c(dmw2,dmw3[grep("S",at)])
  at2 <- c(at2,at3[grep("S",at)])
  # 36S
  dmw2 <- c(dmw2,dmw4[grep("S",at)])
  at2 <- c(at2,at4[grep("S",at)])
}
if (length(grep("Si",at))>0) {
  dmw2 <- c(dmw2,dmw3[grep("Si",at)])
  at2 <- c(at2,at3[grep("Si",at)])
}
if (length(grep("Mg",at))>0) {
  dmw2 <- c(dmw2,dmw3[grep("Mg",at)])
  at2 <- c(at2,at3[grep("Mg",at)])
}
#------------------------------
# 13C mz and double crosses: (2)
merki=i
merkiorg=merki
for (elmt in 2:(merki-1)) {
  for (loop in (elmt+1):merki) {
    tmp=strsplit(res[elmt,4]," ")
    if (length(grep(tmp[[1]][1] ,res[loop,4]))==0) {
      i=i+1
      # rel Int
      res[i,3] <- as.numeric(res[elmt,3])*as.numeric(res[loop,3])/100  
      if (as.numeric(res[i,3])<limri) {
        res[i,3] <- 0
        i=i-1
      } else {
        # m/z
        tmp2=strsplit(res[loop,4]," ")   
        res[i,1] <- moiso + dmw2[grep(tmp[[1]][1],at2)] * as.numeric(tmp[[1]][2])+
                            dmw2[grep(tmp2[[1]][1],at2)] * as.numeric(tmp2[[1]][2]) 
        # isotope
        res[i,4] <- paste(res[elmt,4],res[loop,4])  
      }
    }
  }
}
#------------------------------
# 13C mz and other crosses: (3) 
if (length(at)>3) {           
  merki=merkiorg
  for (tr in 1:(length(at)-3)) {
    merki2=i
    for (elmt in 2:merkiorg) {
      for (loop in (merki+1):merki2) {
        tmp=strsplit(res[elmt,4]," ")
        if (length(grep(tmp[[1]][1] ,res[loop,4]))==0) {
          i=i+1
          # rel Int
          res[i,3] <- as.numeric(res[elmt,3])*as.numeric(res[loop,3])/100  
          if (as.numeric(res[i,3])<limri) {
            res[i,3] <- 0
            i=i-1
          } else {
            # m/z
            tmp2=strsplit(res[loop,4]," ") 
            ltp2=length(tmp2[[1]])/2
            res[i,1] <- moiso + dmw2[grep(tmp[[1]][1],at2)] * as.numeric(tmp[[1]][2])
            for (mii in 2*(1:ltp2)) {
              res[i,1] <- as.numeric(res[i,1])+
                        dmw2[grep(tmp2[[1]][mii-1],at2)] * as.numeric(tmp2[[1]][mii])
            }
            # isotope
            res[i,4] <- paste(res[elmt,4],res[loop,4])
          }
        } else {
          merki=merki+1
        }
      }
      merki=merki-1
    }
    merki=merki2
  }
}
# =============================
# format res:
if (1==1) {
  # remove empty rows:
  res <- res[-which(res[,1]==0),]
  # sort mz:
  res <- res[order(as.numeric(res[,1])),]
  # remove mz doubles:
  ci=0
  for (loop in 2:length(res[,1])) {
    if (res[loop,1]==res[loop-1,1]) {
      ci=c(ci,loop)
    }
  }
  ci <- ci[-1]
  if (length(ci)>0) res <- res[-ci,]
  # highest to 100% rel Int:
  tmp <- max(as.numeric(res[,3]))
  res[,3] <- 100/tmp*as.numeric(res[,3])
  # calc probability: 
  sri <- sum(as.numeric(res[,3]))
  res[,2] <- (as.numeric(res[,3]))/sri
  # remove low  rel Int:
  if (limitfin>limri) res <- res[-which(as.numeric(res[,3])<max(limri,limitfin)),]
  # apply charge:
  tmp <- which(mzv[,2]=="e")
  if (chrg>0) {
    res[,1] <- (as.numeric(res[,1]))/chrg - as.numeric(mzv[tmp,5])
  } else {
    res[,1] <- (as.numeric(res[,1]))/abs(chrg) + as.numeric(mzv[tmp,5])
  }
} # if 
## -------------------------------------------------------------------------
#View(res)
	res <- res[order(as.numeric(res[,3]),decreasing=TRUE),]
  res[,1] <- round(as.numeric(res[,1]),digits =5)
	res[,2:3] <- round(as.numeric(res[,2:3]),digits=3)
  return(res)
}
