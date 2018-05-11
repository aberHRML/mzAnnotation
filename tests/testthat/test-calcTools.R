
context('calculation tools')

test_that('calcAccurateMass works',{
  am <- calcAccurateMass('C4H5O5',charge = 0)  
  expect_true(am == 133.0137)
})

test_that('calcM works',{
  M <- calcM(118.08626,adduct = '[M+H]1+',isotope = '13C',transformation = 'M - [O] + NH2]')
  expect_true(M == 116.05182)
})

test_that('calcMZ works',{
  mz <- calcMZ(116.05182,adduct = '[M+H]1+',isotope = '13C',transformation = 'M - [O] + NH2]')
  expect_true(mz == 118.08626)
})