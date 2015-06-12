xlsxAnnot <-
function(wb,annot_all,dn,Path,DF,plots=F){
	sheet <- createSheet(wb, sheetName = dn)
	bins <- unique(annot_all[,1])
	bins[bins==""] <- NA
	bins <- bins[!(bins=="Bin" | is.na(bins))]
	pos.bin <- NULL
	for (i in 1:length(bins)){
		pos.bin[i] <- which(annot_all[,1]==bins[i])
	}
	pos.bin <- unlist(pos.bin)
	if (plots==T){
		annot_all <- cbind(annot_all[,1:2],matrix("",nrow=nrow(annot_all),ncol=14),annot_all[,3:ncol(annot_all)])
		annot_all[1:2,3:9] <- "Boxplots"
		annot_all[1:2,10:16] <- "Bin Plots"
		border.col <- c(1,2,9,16,21,27,31,37,43)
		fills <- list(1,2,3:9,10:16,17:21,22:27,28:31,32:37,38:43)
		colours <- c("lightblue","lightsalmon","Red","Blue","yellow","bisque1","skyblue","lightgreen","darkorange") 
		pos.bin.1 <- c(pos.bin,nrow(annot_all))
		if (pos.bin.1[length(pos.bin.1)-1]==pos.bin.1[length(pos.bin.1)]){
			annot_all <- rbind(annot_all,matrix("",nrow=2,ncol=ncol(annot_all)))
			pos.bin.1 <- c(pos.bin,nrow(annot_all))
		}
		pos.diff <- c()
		for (i in 1:(length(pos.bin.1)-1)){
			pos.diff[i] <-  pos.bin.1[i+1] - pos.bin.1[i]
		}
		for (i in 1:length(pos.bin)){
			if (pos.diff[i] < 15){
				annot_all <- rbind(annot_all[1:(pos.bin.1[i+1]-1),],matrix("",nrow=15-pos.diff[i],ncol=ncol(annot_all)),annot_all[(pos.bin.1[i+1]):nrow(annot_all),])
				for (i in 1:length(bins)){
					pos.bin[i] <- which(annot_all[,1]==bins[i])
				}
				pos.bin <- unlist(pos.bin)
				pos.bin.1 <- c(pos.bin,nrow(annot_all))
			#	if (pos.bin.1[length(pos.bin.1)-1]==pos.bin.1[length(pos.bin.1)]){
			#		pos.bin.1[length(pos.bin.1)] <- pos.bin.1[length(pos.bin.1)] +2
		  #	}
			}
		}
		for (i in 1:length(bins)){
			pos.bin[i] <- which(annot_all[,1]==bins[i])
		}
		pos.bin <- unlist(pos.bin)
		for (i in 1:length(bins)){
			if (file.exists(paste(Path,DF,paste(DF,"Boxplots",sep="_"),paste(paste(DF,bins[i],sep="_"),".jpeg",sep=""),sep="/"))){
				addPicture(paste(Path,DF,paste(DF,"Boxplots",sep="_"),paste(paste(DF,bins[i],sep="_"),".jpeg",sep=""),sep="/"),sheet,startRow=pos.bin[i]+1,startColumn=4,scale=0.4)
			}
			if (file.exists(paste(Path,DF,paste(DF,"Bin_Plots",sep="_"),paste(paste(DF,"Peaks",bins[i],sep="_"),".jpeg",sep=""),sep="/"))){
				addPicture(paste(Path,DF,paste(DF,"Bin_Plots",sep="_"),paste(paste(DF,"Peaks",bins[i],sep="_"),".jpeg",sep=""),sep="/"),sheet,startRow=pos.bin[i]+1,startColumn=11,scale=0.2)
			}
		}
	} else {
		border.col <- c(1,2,7,13,17,23,29)
		fills <- list(1,2,3:7,8:13,14:17,18:23,24:29)
		colours <- c("lightblue","lightsalmon","yellow","bisque1","skyblue","lightgreen","darkorange") 
	}
	cb <- CellBlock(sheet, 1, 1, nrow(annot_all), ncol(annot_all))
	alignment <- Alignment(horizontal="ALIGN_CENTER")
	cell.style <- CellStyle(wb, alignment=alignment)
	CB.setMatrixData(cb, annot_all, 1, 1,,cellStyle=cell.style)
	border.r <- Border(position="RIGHT")
	border.b <- Border(position="TOP")
	for (i in pos.bin){
		CB.setBorder(cb, border.b, i, 1:(ncol(annot_all)))
	}
	for (i in border.col){
		CB.setBorder(cb, border.r, 1:(nrow(annot_all)), i)
	}
	for (i in 1:length(colours)){
		fill <- Fill(foregroundColor = colours[i], backgroundColor=colours[i])
		fills.1 <- fills[[i]]
		for (x in 1:length(fills.1)){
			CB.setFill(cb, fill, 1:2, fills.1[x])
		}	  
	}
	font <-  Font(wb, heightInPoints=8)
	for (i in 1:ncol(annot_all)){
		CB.setFont(cb, font,1:(nrow(annot_all)),i)
	}
	return(wb)
}
