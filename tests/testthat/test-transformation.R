
context('transformation')

test_that('adductTransfromMF works',{
  mf <-  adductTransformMF('C4H5O5','[M-H]1-')
  expect_true(mf == "C4H4O5")
})

test_that('transformMF works',{
  transformed_MF <- transformMF('C4H5O5')
  
  expect_equal(transformed_MF, "C4H7NO4")
})