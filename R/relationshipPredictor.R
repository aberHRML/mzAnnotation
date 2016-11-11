
mz <- c(66.00349,133.01425,134.0176,168.99093,267.03577)

relationshipPredictor <- function(mz,mode,limit=0.001){

      diffMatrix <- data.frame(as.matrix(dist(mz)))
      colnames(diffMatrix) <- mz
      rownames(diffMatrix) <- mz
  
}
