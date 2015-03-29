annot_sheet <-
function(x) {
  annot_all <- matrix(ncol=29)
  bins <- names(x)
  for (i in 1:length(x)){
  	print(i)
    current.mass <- x[[i]]
    acc.mass <- current.mass[[1]]
    cors <- current.mass[[2]]
    cors <- cors[[1]]
    mfs <- current.mass[[3]]
    isodists <- current.mass[[4]]
    sum_meas_iso <- current.mass[[6]]
    mzeddb <- current.mass[[7]]
	if (class(mfs)=="character"){
    	mfs <- matrix(mfs,ncol=6)
    }
    if(length(isodists)>0 & !is.null(names(isodists))){
      isod.rows <- sapply(isodists,nrow)
      isod.len <- sum(isod.rows)+nrow(mfs)-length(isodists) 
      iso.mat <- matrix(ncol=10)
      for (y in 1:nrow(mfs)){
        if(mfs[y,3] %in% names(isodists)){
          posi <- grep(mfs[y,3],names(isodists))
          iso.mat.1 <-matrix(ncol=10,nrow=nrow(isodists[[posi]]))
          iso.mat.1[1,1:6] <- as.matrix(mfs[y,])
          iso.mat.1[1:nrow(isodists[[posi]]),7:10] <- isodists[[posi]]
          iso.mat <- rbind(iso.mat,iso.mat.1)
        } else {
          mat <- matrix(ncol=10)
          mat[1,1:6] <- as.matrix(mfs[y,])
          iso.mat <- rbind(iso.mat,mat)
        }        
      }
      iso.mat <- iso.mat[-1,]
    } else {
      iso.len <- 1
      isod.len <- 0
    }
    num.rows <- c(nrow(cors),nrow(mfs),nrow(mzeddb))
    max.rows <- max(num.rows)
    if (isod.len > max.rows){
    	n.rows <- isod.len
    } else {
    	n.rows <- max.rows
    }
    curr.mat <- matrix(ncol=29,nrow=n.rows)
  	if(nrow(curr.mat)==0){
  		curr.mat <- rbind(matrix(nrow=1,ncol=29))
  	}
    curr.mat[1,1] <- bins[i] 
    curr.mat[1,2] <- acc.mass
    if(length(cors)>1){
    	cors[,2] <- round(as.numeric(as.character(cors[,2])),3)
    	curr.mat[1:nrow(cors),3:7] <- as.matrix(cors)
    }
    if(length(isodists)>0 & !is.null(names(isodists))){
    	curr.mat[1:nrow(iso.mat),8:17] <- iso.mat
    }
    ###### Add isodists here ######
    if(length(sum_meas_iso)>1){
    	sum_meas_iso[,3:6] <- round(as.numeric(as.character(sum_meas_iso[,3:6])),3)
    	curr.mat[1,19:23] <- as.matrix(sum_meas_iso[,2:6])
    	curr.mat[1,18] <- "Carbon"
    }
  	if (class(mzeddb)=="character"){
    	mzeddb <- matrix(mzeddb,ncol=6)
    }
    if(!is.null(nrow(mzeddb))){
    	mzeddb[,6] <- round(as.numeric(as.character(mzeddb[,6])),2)
    	curr.mat[1:nrow(mzeddb),24:29] <- as.matrix(mzeddb)
    }
    annot_all <- rbind(annot_all,curr.mat)
  }
  annot_all[1,] <- c("Bin","Accurate Mass","m/z","Correlation","Difference","Relation_to","Relation_from",
  									 "m/z","Original MF","Clean MF","Theoretical m/z","No. Isotopes",
  									 "PPM Error","m/z","Probability","Relative Intensity","Isotope","Element","Number","Mean","Standard Deviation",
  									 "Standard Error","95% Conf Int","m/z","Theoretical m/z","MF","Name","Adduct","PPM Error")
  title <- c("Bin","Accurate Mass",rep("Correlations",5),rep("MF's",6),rep("Theoretical Isotope Distributions",4),rep("Measured Isotope Distirbutions",6),rep("MZedDB Hits",6))
  annot_all <- rbind(title,annot_all)
  annot_all[is.na(annot_all)] <- ""
  return(annot_all)
}
