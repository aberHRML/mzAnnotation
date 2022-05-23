
metabolite_database <- metaboliteDB(amino_acids)

test_that('PIPsearch works',{
  res <- PIPsearch(
    metabolite_database,
    mz = 133.03358,
    adduct = '[M-H]1-',
    isotope = '13C')
  
  expect_s3_class(res,'tbl_df')
})

test_that('Adducts can be calculated',{
  res <- calcAdducts(
    metabolite_database,
    1
  )
  
  expect_equal(nrow(res),58)
})