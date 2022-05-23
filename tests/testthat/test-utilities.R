
metabilite_database <- metaboliteDB(amino_acids)

test_that('PIPsearch works',{
  res <- PIPsearch(metabilite_database,133.03358,5,'[M-H]1-','13C')
  
  expect_s3_class(res,'tbl_df')
})

