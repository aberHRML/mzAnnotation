
context('isotopic distributions')

test_that('isotopeDistribution works',{
  res <- isotopeDistribution('C4H5O5',charge = -1)
  
  expect_false(F %in% (class(res) == c("tbl_df","tbl","data.frame")))
  expect_true(nrow(res) == 8)
  expect_true(ncol(res) == 4)
  expect_false(F %in% (colnames(res) == c("Isotope","m/z","Relative Abundance","Probability")))
  expect_false(F %in% (res$Isotope == c(NA,"13C 1","18O 1","17O 1","13C 2","2H 1","13C 1; 18O 1","13C 1; 17O 1")))
  expect_false(F %in% (round(res$`Relative Abundance`,5) == c(1.00000,0.04490,0.01002,
                                                              0.00200,0.00076,0.00050,
                                                              0.00045,0.00009)))
})