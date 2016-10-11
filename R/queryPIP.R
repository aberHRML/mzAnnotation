
queryPIP <-
function(add,MZedDB,mz,ppm,filter=T,iso=NULL,bio=NULL){
  # Retrieve database entries
  all <- MZedDB$MZedDB_ALL
  metrules <- MZedDB$MZedDB_METRULES
  
  # Retrieve adduct formation rules
  add_rules <- MZedDB$ADDUCT_FORMATION_RULES
  add <- add_rules[which(add_rules[,1]==add),]
  
  # If an isotope has been selected; retieve isotope formation rules and rules for individual metabolites
  if(!is.null(bio)){
    bio_rules <- MZedDB$BIOTRANSFORMATION_RULES
    bio <- bio_rules[which(bio_rules[,1]==bio),]
  }
  
  # If a biotransformation has been selected; retieve isotope formation rules and rules for individual metabolites
  if(!is.null(iso)){
    iso_rules <- MZedDB$ISOTOPE_RULES
    iso <- iso_rules[which(iso_rules[,1]==iso),]
    elementFreq <- MZedDB$MZedDB_ELEMENT_FREQUENCIES
  }
  
  # Calculate the mz range over which to search based on ppm error
  mass_low <- mz - (mz * ppm * 10^-6)
  mass_high <- mz + (mz * ppm * 10^-6)
  
  # If a biotransformation has been selected; account for this in the mz range
  if(!is.null(bio)){
    mass_low <- mass_low + bio$Mass.Difference
    mass_high <- mass_high + bio$Mass.Difference
  }
  
  # Account for adduct transformation in mz range
  mass_low <- ((mass_low-add[1,"Add"])*add[1,"Charge"])/add[1,"xM"]
  mass_high <- ((mass_high-add[1,"Add"])*add[1,"Charge"])/add[1,"xM"]
  
  # If an isotope has been selected; account for this in the mz range
  if(!is.null(iso)){
    mass_low <- mass_low - iso$Mass.Difference
    mass_high <- mass_high - iso$Mass.Difference
  }
  
  # Filter database entries based mz range
  all <- all[which(all[,"Accurate.Mass"]>mass_low & all[,"Accurate.Mass"] <mass_high),]
  
  # If an isotope has been selected; filter the database entries based on the isotope rules
  if(!is.null(iso)){
    elementFreq <- elementFreq[elementFreq$ID %in% all$ID,]
    C <- elementFreq[,'C']
    O <- elementFreq[,'O']
    Cl <- elementFreq[,'Cl']
    K <- elementFreq[,'K']
    S <- elementFreq[,'S']
  
    all <- all[which(eval(parse(text=as.character(iso["Rule"])))),]
  }
  
  # Filter the database based on the adduct formation rules
  metrules <- metrules[metrules$ID %in% all$ID,]
  Nch <- metrules[,"Nch"]
  Nacc <- metrules[,"Nacc"]
  Ndon <- metrules[,"Ndon"]
  Nnhh <- metrules[,"Nnhh"]
  Noh <-  metrules[,"Noh"]
  Ncooh <- metrules[,"Ncooh"]
  Ncoo <-  metrules[,"Ncoo"]
  
  res <- all[which(eval(parse(text=as.character(add["Rule"])))),]
  
  # Calculate ppm errors of results
  if(!is.null(bio)){
    addmz <- sapply(res[,"Accurate.Mass"],function(x,rules){y <- x - rules[1,'Mass.Difference'];return(y)},rules=bio)
    addmz <- sapply(res[,"Accurate.Mass"],function(x,rules){y <- x + rules[1,'Mass.Difference'];return(y)},rules=iso)
    addmz <- sapply(addmz,function(x,rules){y <- (x*rules[1,"xM"])/rules[1,"Charge"]+rules[1,"Add"];return(y)},rules=add)
  } else {
    if(!is.null(iso)){
      addmz <- sapply(res[,"Accurate.Mass"],function(x,rules){y <- x + rules[1,'Mass.Difference'];return(y)},rules=iso)
      addmz <- sapply(addmz,function(x,rules){y <- (x*rules[1,"xM"])/rules[1,"Charge"]+rules[1,"Add"];return(y)},rules=add)
    } else {
      addmz <- sapply(res[,"Accurate.Mass"],function(x,rules){y <- (x*rules[1,"xM"])/rules[1,"Charge"]+rules[1,"Add"];return(y)},rules=add)
    }
  }
	ppmerr <- sapply(addmz,function(x,mz){y <- (x-mz)/mz*10^6;return(y)},mz=mz)
 	
	# Add adduct and isotope information
	adduct <- rep(add[1,1],nrow(res))
 	if(!is.null(iso)){
 	  isotope <- rep(as.character(iso[1,1]),nrow(res))
 	} else {
 	  isotope <- rep('',nrow(res))
 	}
	if(!is.null(bio)){
	  biotransformation <- rep(as.character(bio[2,1]),nrow(res))
	} else {
	  biotransformation <- rep('',nrow(res))
	}
 	
	# Format results
  res <- data.frame(res,Adduct=adduct,Isotope=isotope,Biotransformation=biotransformation,Adduct_MZ=addmz,PPMErr=ppmerr)
  res <- res[,-c(3,6,7,8,10)]
  if(filter==T & nrow(res)>0){
   res <- filterPIP(res)
  }
  
  return(res)
}
