
test_that("error thrown if incorrect adduct table specified", {
  expect_error(checkAdductTable(tibble()))
})

test_that("error thrown if incorrect adduct table specified", {
  expect_error(checkIsotopeTable(tibble()))
})

test_that("error thrown if incorrect adduct table specified", {
  expect_error(checkTransformationTable(tibble()))
})

test_that('error thrown if incorrect transformation specified',{
  expect_error(checkTransformation('incorrect'))
})

test_that('error trhown if incorrect entries specified',{
  expect_error(checkEntries(tibble()))
})
