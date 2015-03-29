acc.mass <-
function(master_mat,cls,masses){		# collects accurate masses(5dp) for a given list of explanatory masses(2dp) based on the mosted intense signal
		mat.acc <- NULL
		masses.p <- as.numeric(gsub("p","",masses[[1]]))
		masses.n <- as.numeric(gsub("n","",masses[[2]]))
		mass <- list(masses.p,masses.n)
		for (i in 1:length(master_mat)){
			acc.mat <- matrix(0,nrow=length(mass[[i]]),ncol=2*length(unique(cls)))
			mat <- master_mat[[i]]
			masses.1 <- mass[[i]]
			for (x in 1:length(masses.1)){
				mat.mass <- mat[which(mat[,3]==masses.1[x]),1:ncol(mat)]
				if(class(mat.mass)=="numeric"){
					mat.mass <- matrix(mat.mass,ncol=length(mat.mass))
				}
				mass.acc <- matrix(0,nrow=1,ncol=2*length(unique(cls)))
				s.1 <- seq(1,2*length(unique(cls)),2)
				for (y in 1:(ncol(mat.mass)-4)){
					m <- mat.mass[grep(max(mat.mass[,y+4]),mat.mass[,y+4]),4]
					int <- max(mat.mass[,y+4])
					int <- round(int,0)
					mi <- c(m[1],int)
					if (length(m)>1){
						mi <- c(m[1],int)
					}else{
						mi <- c(m,int)
					}
					mass.acc[,s.1[y]:(s.1[y]+1)] <- mi
				}
				acc.mat[x,] <- mass.acc
				rownames(acc.mat) <- masses[[i]]
				s.2 <- seq(1,2*length(unique(cls)),2)
				nam <- matrix(0,nrow=1,ncol=2*length(unique(cls)))
				nam[,s.2] <- as.character(unique(cls))
				nam[,s.2+1] <- "Intensity"
				names <- nam[1,]
				colnames(acc.mat) <- names
		}	
		mat.acc[i] <- list(acc.mat)
	}
	return(mat.acc)
}
