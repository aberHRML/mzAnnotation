#' filterPIP
#' @export

filterPIP <- 
  function(res){
		patterns <- c("ic acid","keto","<i>","</i>","(R)-","(S)-","-L-","L-","(+)-",
                  "cis,","cis-","trans-","&beta;-","-D-",
                  "D-","D.","&gamma;-","-n-","N-","&alpha;-","&alpha;","&beta;")
		replacements <- c("ate","oxo","","","","","","","","","","","","","","","","","","","","")
		for(i in 1:length(patterns)){
			res[,2] <- sapply(res[,2],gsub,pattern=patterns[i],replacement=replacements[i],fixed=TRUE)
		}
		res[,2] <- tolower(res[,2])
		res <- res[!duplicated(res[,2]),]
		return(res)
  }