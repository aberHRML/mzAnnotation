calc_diff <-
function(mat,mo){
	mat <- na.omit(mat)
	par <- colnames(mat)[1]
	par <- as.numeric(gsub(mo,"",par))
	if(mo=="n"){
		add <- c("iC13","iO18, iCl37, iS34","[M+Cl]1-","[M-2H]2-","[2M-H]1-","Same Peak?",
						 "-[H]+[CH2]","-[H]+[OH]","[M+K-2H]1-","-[NO]+[H]","-[PO3H2]+[H]",
						 "-[CO2]+[H]","+[H2]-[O2]","-[H2O]","-[H2]","-[O]+[NH2]","-[OH]+[NH2]",
						 "+[H2]","-[H2]+[O]","-[H]+[NH2]","+[H2O]","-2[H]+2[CH2]","-[H3]+[HO2]",
							"-2[H]+2[OH]","-[H]+[C2H2O]","-[H] +[CHO2]","-[H]+[SO3]","-[H]+[PO3H2]",
						 "-2[H]+[SO3]+[CH2]","-2[H]+2[SO3]","-2[H]+[C2H2O]+[C3H7NO2S]","-3[H]+3[SO3]",
						 "-[H]+[C6H8O6]","-2[H]+[C6H8O6]+[CH2]","-[H]+[C6H10O5]","-2[H]+[C6H10O5]+[SO3]",
						 "-2[H]+[C6H8O6]+[SO3]","-[H]+[C10H17N3O6S]","-2[H]+2[C6H8O6]","-3[H]+2[C6H8O6]+[CH2]",
						 "[M+Na-2H]1-")
		diff <- c(1.00,2.00, 35.98,-(par-1.01)/2,(par+1.01)*2,0.01,14.02,15.99,37.96,
							-44.99,-79.97,-43.99,-29.97,-18.01,-2.02,-0.02,0.98,2.02,13.98,15.01
							,18.01,28.03,29.97,31.99,42.01,43.99,79.96,79.97,93.97,159.91,163.03,
							239.87,176.03,190.05,162.05,242.01,255.99,307.08,352.06,366.08,21.98)
	}
	if(mo=="p"){
		add <- c("iC13","iO18, iK41, iS34","[M+Na]1+","[M+K]1+","K-Na",
						 "[M+H-H2O]1+, -[H2O]","[M+H-FA]1+","[M+NH4]1+","[2M+H]1+","[2M+Na]1+",
						 "[2M+K]1+","Same Peak?","-[H]+[CH2]","-[H]+[OH]","-[NO]+[H]",
						 "-[PO3H2]+[H]","-[CO2]+[H]","+[H2]-[O2]","-[H2]","-[O]+[NH2]","-[OH]+[NH2]",
						 "+[H2]","-[H2]+[O]","-[H]+[NH2]","+[H2O]","-2[H]+2[CH2]","-[H3]+[HO2]",
						 "-2[H]+2[OH]","-[H]+[C2H2O]","-[H] +[CHO2]","-[H]+[SO3]","-[H]+[PO3H2]",
						 "-2[H]+[SO3]+[CH2]","-2[H]+2[SO3]","-2[H]+[C2H2O]+[C3H7NO2S]","-3[H]+3[SO3]",
						 "-[H]+[C6H8O6]","-2[H]+[C6H8O6]+[CH2]","-[H]+[C6H10O5]","-2[H]+[C6H10O5]+[SO3]",
						 "-2[H]+[C6H8O6]+[SO3]","-[H]+[C10H17N3O6S]","[M+H-NH3]1+","[M+2Na-H]1+","[M+2K-H]1+")
		diff <- c(1.00,2.00,21.98,37.96,15.97,-18.01,-46.01,17.03,par*2-1.01,par*2-1.01+21.98
							,par*2-1.01+37.96,0.01,14.02,15.99,-44.99,-79.97,-43.99,-29.97,-2.02,-0.02,
							0.98,2.02,13.98,15.01,18.01,28.03,29.97,31.99,42.01,43.99,79.96,79.97,93.97,
							159.91,163.03,239.87,176.03,190.05,162.05,242.01,255.99,307.08,352.06,366.08,
							-17.03,43.96,75.91)
	}
	cors <- mat[,1]
	cors <- as.numeric(gsub(mo,"",cors))
	cors.1 <- round(cors-par,2)
	mat.1 <- cbind(mat,cors.1)
	colnames(mat.1)[3] <- "Difference"
	posi.add <-sapply(mat.1[,3],function(x,dif,ad){matc <- x==dif;if (TRUE %in% matc) {add <- ad[matc]} else {add <- NA};return(add)},dif=diff,ad=add)
	add.itself <- sapply(mat.1[,3],function(x,dif,ad){matc <- x==dif;if (TRUE %in% matc) {add <- ad[matc]} else {add <- NA};return(add)},dif=-diff,ad=add)
	mat.1 <- cbind(mat.1,posi.add,add.itself)
	rownames(mat.1) <- NULL
	return(data.frame(mat.1))
}
