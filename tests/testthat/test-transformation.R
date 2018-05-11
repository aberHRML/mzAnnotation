
context('transformation')

test_that('adductTransfromMF works',{
  mf <-  adductTransformMF('C4H5O5','[M-H]1-')
  expect_true(mf == "C4H4O5")
})