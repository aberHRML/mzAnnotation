
context('molecular formula generation')

test_that('generateMF works',{
  res <- generateMF(342.11621,
                    element_max = c(C = 12,H = 22,N = 0,
                                  O = 11,P = 0,S = 0))
  expect_false(F %in% (class(res) == c("tbl_df","tbl","data.frame")))
  expect_true(nrow(res) == 1)
  expect_true(ncol(res) == 3)
  expect_false(F %in% (colnames(res) == c( "MF","Mass","PPM Error")))
  expect_true(res$MF == 'C12H22O11')
  expect_true(res$Mass == 342.11621)
  expect_true(res$`PPM Error` == 0)
})