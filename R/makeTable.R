makeTable <-
  function(annot_info){
    # combine iso and mf tables if present
    MF.iso <- list()
    if("Molecular Formulas" %in% names(annot_info) & "Theoretical Isotope Distributions" %in% names(annot_info)){
      if (nrow(annot_info[["Molecular Formulas"]])>0 & !is.null(annot_info[["Theoretical Isotope Distributions"]])){
        MF.tab <- annot_info[["Molecular Formulas"]]
        MF.tab <- matrix(apply(MF.tab,2,as.character),ncol=ncol(MF.tab))
        colnames(MF.tab) <- colnames(annot_info[["Molecular Formulas"]])
        Iso.tab <- annot_info[["Theoretical Isotope Distributions"]]
        for (i in 1:nrow(MF.tab)){
          if(MF.tab[i,3] %in% names(Iso.tab)){
            iso.cur <- Iso.tab[[which(names(Iso.tab)==MF.tab[i,3])]]
          } else {
            iso.cur <- matrix(ncol=4,nrow=1)
          }
          colnames(iso.cur) <- c("m.z.1","Prob","rel.Int","Isotope")
          MF.tab.cur <- matrix(ncol=ncol(MF.tab),nrow=nrow(iso.cur)-1)
          colnames(MF.tab.cur) <- colnames(MF.tab)
          MF.tab.cur <- rbind(as.character(MF.tab[i,]),MF.tab.cur)
          MF.iso.tab <- cbind(MF.tab.cur,iso.cur)
          MF.iso[[i]] <- MF.iso.tab
        }
        MF.iso <- ldply(MF.iso,data.frame,stringsAsFactors=F)
        MF.iso <- apply(MF.iso,2,as.character)
        annot_info$`Molecular Formulas` <- MF.iso
        annot_info <- annot_info[-which(names(annot_info) == "Theoretical Isotope Distributions")]
      } else {
        if (nrow(annot_info[["Molecular Formulas"]])>0){
          MF.iso <- cbind(annot_info[["Molecular Formulas"]],matrix(ncol=4,nrow=nrow(annot_info[["Molecular Formulas"]])))
        } else {
          MF.iso <- matrix(ncol=10,nrow=1)
        }
          annot_info$`Molecular Formulas` <- MF.iso
          annot_info <- annot_info[-which(names(annot_info) == "Theoretical Isotope Distributions")]
      }
    }
    # rename columns for each table
    for (i in names(annot_info)){
      if (i == "Accurate m/z"){
        colnames(annot_info[[i]]) <- c("Bin","Accurate Mass","Intensity")
      }
      if (i == "Correlation Analysis"){
        colnames(annot_info[[i]]) <- c("m/z","Coefficient","Difference","Relation_to","Relation_from")
      }
      if (i == "Molecular Formulas"){
        colnames(annot_info[[i]]) <- c( "m/z","Original MF","Clean MF","Theoretical m/z","RDB",
                                 "PPM Error","m/z","Probability","Relative Intensity","Isotope")
      }
      if (i == "Putative Ionisation Product"){
        colnames(annot_info[[i]]) <- c("m/z","Theoretical m/z","MF","Name","Adduct","PPM Error")
      }
    }
    # calculate number of rows needed for table
    row <- list()
    for (i in 1:length(annot_info)){
        row[[i]] <- dim(annot_info[[i]])[1]
    }
    names(row) <- names(annot_info)
    row <- max(unlist(row))
    # expand each table to the max row number
    annot_info <- lapply(annot_info, function(x,row){
      row.add <- row - nrow(x)
      new.tab <- matrix(ncol=ncol(x),nrow=row.add)
      colnames(new.tab) <- colnames(x)
      new.tab <- rbind(x,new.tab)
      return(new.tab)
    },row=row)
    # combine each table
    annot_table <- do.call("cbind", annot_info)
    return(annot_table)
  }