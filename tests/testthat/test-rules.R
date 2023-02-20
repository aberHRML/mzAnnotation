test_that("adduct names returned", {
  expect_identical(adduct_names(),
                   adduct_rules()$Name)
})

test_that("isotope names returned", {
  expect_identical(isotope_names(),
                   isotope_rules()$Isotope)
})

test_that("adduct names returned", {
  expect_identical(transformation_names(),
                   transformation_rules()$`MF Change`)
})
