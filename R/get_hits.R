get_hits <-
function(add_cur,limit_low.1,limit_high.1,add_masses,mzeddb.1){		
		query_masses <- which(add_masses[,add_cur] > limit_low.1 & add_masses[,add_cur] < limit_high.1)
		query_hits <- as.matrix(mzeddb.1[query_masses,])
		if(nrow(query_hits)>0){
			hits <- matrix(ncol=(ncol(query_hits)+2),nrow=nrow(query_hits))
			hits[1:nrow(query_hits),1:ncol(query_hits)] <- query_hits
			hits[1:nrow(query_hits),(ncol(query_hits)+1)] <- rep(add_cur,nrow(query_hits))
			hits[1:nrow(query_hits),(ncol(query_hits)+2)] <- add_masses[query_masses,add_cur]
			hits <- hits[,-c(3,6,7,8)]
			if(class(hits)=="character"){
				hits <- matrix(hits,nrow=1)
			}
		}
	  if (exists("hits")){
			return(hits)
	  }
}
