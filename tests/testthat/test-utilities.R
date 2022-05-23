
metabolite_database <- metaboliteDB(amino_acids)

test_that('Metabolite database show method works',{
  expect_output(print(metabolite_database),
                'MetaboliteDatabase')
})

test_that('Metabolite database can be filtered by molecular formula',{
  res <- filterMF(
    metabolite_database,
    "C3H7NO2"
  )
  
  expect_equal(nEntries(res),1)
})


test_that('filterER returns correctly when no matches found',{
  res <- filterER(
    metabolite_database,
    K > 1
  )
  
  expect_equal(nEntries(res),0)
})
