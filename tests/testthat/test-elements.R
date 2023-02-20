test_that("MF element frequencies can be calculated", {
  element_frequencies <- elementFrequencies(c('H2O','C12H22O11'))
  
  expect_s3_class(element_frequencies,'tbl_df')
})

test_that('MF element frequency ratios can be calculated',{
  element_ratios <- elementFrequencies(c('H2O','C12H22O11')) %>% 
    elementRatios()
  
  expect_s3_class(element_ratios,'tbl_df')
})
