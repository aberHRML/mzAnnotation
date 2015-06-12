annotSheet <-
function(annot_all,genpip=T) {
  # Caluclate number of columns in sheet
  bins <- length(annot_all[[1]])
  # Build table for each bin
  for (i in 1:bins){
    annot_info <- makeTable(lapply(annot_all,function(x,y){return(x[[y]])},y=i))
  } 
  annot_info <- ldply(annot_info,data.frame)
  return(annot_info)
}
