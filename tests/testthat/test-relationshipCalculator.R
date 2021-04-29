
context('relationship calculation')

test_that('relationshipCalculator works',{
  rel <- relationshipCalculator(c(132.03023,172.00067),
                                adducts = c("[M-H]1-","[M+Cl]1-","[M+K-2H]1-",'[M+H]1+','[M+K]1+','[M+Na]1+'))
  
  expect_false(F %in% (class(rel) == c("tbl_df","tbl","data.frame")))
  expect_true(nrow(rel) == 1)
  expect_true(ncol(rel) == 9)
  expect_false(F %in% (colnames(rel) == c("m/z1","m/z2","Adduct1",
                                          "Adduct2","Isotope1","Isotope2",
                                          "Transformation1","Transformation2",
                                          "Error" )))
  expect_true(rel$Adduct1[1] == '[M-H]1-')
  expect_true(rel$Adduct2[1] == '[M+K]1+')
  expect_true(rel$Error[1] == 0)
})