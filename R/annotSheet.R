annotSheet <-
function(annot_all) {
  # put accurate m/z into list
  annot_all$`Accurate m/z` <- lapply(apply(annot_all$`Accurate m/z`,1, function(x){return(list(x))}),function(x){return(matrix(x[[1]],ncol=3))})
  # Caluclate number of columns in sheet
  bins <- length(annot_all[[1]])
  # Build table for each bin
  annot_info <- list()
  for (i in 1:bins){
    annot_info[[i]] <- makeTable(lapply(annot_all,function(x,y){return(x[[y]])},y=i))
  } 
  annot_info <- do.call(rbind,annot_info)
  return(annot_info)
}
