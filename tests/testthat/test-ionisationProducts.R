
context('ionisation products')

test_that('ionisationProducts',{
  ips <- ionisationProducts(aminoAcids$SMILES[1])
  expect_equal(nrow(ips),58)
  expect_equal(ncol(ips),4)
  expect_equal(ips$`m/z`[1],52.67178)
})