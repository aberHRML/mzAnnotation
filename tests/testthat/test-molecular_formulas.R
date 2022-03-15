
context('molecular formula generation')

test_that('generateMF works',{
  res <- generateMF(342.11621)
  
  expect_s3_class(res,"tbl_df")
  expect_equal(nrow(res),25)
  expect_equal(ncol(res),3)
  expect_identical(colnames(res),c( "MF","Mass","PPM error"))
})

test_that('generateMF returns correctly when no MFs are generated',{
  expect_s3_class(generateMF(111),"tbl_df")
})

test_that('ipMF works',{
  res <- ipMF(118.08626,adduct = '[M+H]1+')
  
  expect_s3_class(res,"tbl_df")
})

test_that('ipMF returns correctly when no molecular formulas are generated',{
  res <- ipMF(118.08630,
              adduct = '[M+H]1+',
              ppm = 0.1)
  
  expect_equal(nrow(res),0)
})


test_that('ipMF works for mass below 100',{
  res <- ipMF(99,adduct = '[M+H]1+')
  
  expect_s3_class(res,"tbl_df")
})

test_that('ipMF works for mass over 200',{
  res <- ipMF(202,
              adduct = '[M+H]1+',
              isotope = '13C')
  
  expect_s3_class(res,"tbl_df")
})

test_that('isotopePossible throws error if incorrect isotope specified',{
  expect_error(isotopePossible('H2O',isotope = 'incorrect'))
})
