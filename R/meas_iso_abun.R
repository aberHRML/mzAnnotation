meas_iso_abun <-
function(m,acc_mat,mo){
	colnames(acc_mat) <- gsub(mo,"",colnames(acc_mat))
	i.u.rang <- (m + 1.0033548) + ((m + 1.0033548) * 1 * 10^-6)
	i.l.rang <- (m + 1.0033548) - ((m + 1.0033548) * 1 * 10^-6)
	p.u.rang <- m +  ((m + 1.0033548) * 1 * 10^-6)
	p.l.rang <- m -  ((m + 1.0033548) * 1 * 10^-6)
	m.nam.i <- as.numeric(colnames(acc_mat))
	m.nam.i <- m.nam.i[m.nam.i>i.l.rang & m.nam.i<i.u.rang]
	m.nam.p <- as.numeric(colnames(acc_mat))
	m.nam.p <- m.nam.p[m.nam.p>p.l.rang & m.nam.p<p.u.rang]
	iso.mass <- acc_mat[,as.character(m.nam.i)]
	par.mass <- acc_mat[,as.character(m.nam.p)]
	if (length(iso.mass)>0){
		if(class(iso.mass)=="numeric"){
			nam <- colnames(acc_mat)[grep(as.character(m.nam.i),colnames(acc_mat))]
			iso.max <- col.max(iso.mass,name=nam)
		} else {
			iso.max <- col.max(iso.mass)
		}
		par.max <- col.max(par.mass)
		if(ncol(par.max)==1){
			par.max <- cbind(rep(m.nam.p,nrow(par.max)),par.max)
			colnames(par.max) <- c("mz","intensity")
		}
		rat <- as.numeric(iso.max[,"intensity"])/as.numeric(par.max[,"intensity"])*100
		num <- rat/1.07
		all.max <- data.frame(Parent_mz=par.max[,"mz"],Parent_intensity=par.max[,"intensity"],Iso_mz=iso.max,Iso_intensity=iso.max[,"intensity"],Relative_Intensity=rat,No_Element=num)
	} else {
		all.max <- "No masses found in isotope range"
	}
	return(all.max)
}
