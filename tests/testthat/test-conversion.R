
context('conversion')

test_that('conversion works',{
  testthat::skip_on_os('windows')
  
  inchi <- convert(amino_acids$SMILES[1],'smiles','inchi')
  expect_true(identical(inchi,"InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1"))
})
