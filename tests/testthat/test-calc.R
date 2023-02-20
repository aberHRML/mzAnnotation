
context('calculation tools')

test_that('calcAccurateMass works',{
  am <- calcAccurateMass('C4H5O5',charge = 0)  
  expect_equal(am,133.0137)
})

test_that('calculate accurate mass with charge != 0',{
  am <- calcAccurateMass('C4H5O5',charge = 1)  
  expect_equal(am,133.01315)
})

test_that('calcM works',{
  M <- calcM(
    118.08626,
    adduct = '[M+H]1+',
    isotope = '13C',
    transformation = 'M - [O] + [NH2]')
  expect_equal(M,116.05182)
})

test_that('calcM errors when incorrect adduct specified',{
  expect_error(calcM(
    118.08626,
    adduct = 'wrong',
    isotope = '13C',
    transformation = 'M - [O] + [NH2]'))
})

test_that('calcM errors when incorrect isotope specified',{
  expect_error(calcM(
    118.08626,
    adduct = '[M+H]1+',
    isotope = 'wrong',
    transformation = 'M - [O] + [NH2]'))
})

test_that('calcM errors when incorrect transformation specified',{
  expect_error(calcM(
    118.08626,
    adduct = '[M+H]1+',
    isotope = '13C',
    transformation = 'wrong'))
})

test_that('calcMZ works',{
  mz <- calcMZ(
    116.05182,
    adduct = '[M+H]1+',
    isotope = '13C',
    transformation = 'M - [O] + [NH2]')
  expect_equal(mz,118.08626)
})

test_that('calcMZ errors when incorrect adduct specified',{
  expect_error(calcMZ(
    116.05182,
    adduct = 'wrong',
    isotope = '13C',
    transformation = 'M - [O] + [NH2]'))
})

test_that('calcMZ errors when incorrect isotope specified',{
  expect_error(calcMZ(
    116.05182,
    adduct = '[M+H]1+',
    isotope = 'wrong',
    transformation = 'M - [O] + [NH2]'))
})

test_that('calcMZ errors when incorrect transformation specified',{
  expect_error(calcMZ(
    116.05182,
    adduct = '[M+H]1+',
    isotope = '13C',
    transformation = 'wrong'))
})
