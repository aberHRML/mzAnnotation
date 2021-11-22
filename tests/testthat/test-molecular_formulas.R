
context('molecular formula generation')

test_that('generateMF works',{
  res <- generateMF(342.11621)
  
  expect_s3_class(res,"tbl_df")
  expect_equal(nrow(res),4)
  expect_equal(ncol(res),3)
  expect_identical(colnames(res),c( "MF","Mass","PPM Error"))
})

test_that('generateMF returns correctly when no MFs are generated',{
  expect_s3_class(generateMF(1),"tbl_df")
})

test_that('ipMF works',{
  res <- ipMF(118.08626,adduct = '[M+H]1+')
  
  expect_s3_class(res,"tbl_df")
})

test_that('ipMF works for mass below 100',{
  res <- ipMF(99,adduct = '[M+H]1+')
  
  expect_s3_class(res,"tbl_df")
})

test_that('ipMF works for mass over 200',{
  res <- ipMF(202,adduct = '[M+H]1+')
  
  expect_s3_class(res,"tbl_df")
})