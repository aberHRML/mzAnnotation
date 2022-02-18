test_that("Heteroatom ratios can be checked", {
  ratio_checks <- elementFrequencies(c('H2O','C12H22O11')) %>% 
    elementRatios() %>% 
    heteroatomRatioCheck()
  
  expect_s3_class(ratio_checks,'tbl_df')
})
