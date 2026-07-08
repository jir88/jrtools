test_that("multiplicative_imputation preserves rows without zeros", {
  X <- matrix(c(2, 3), nrow = 1)

  result <- multiplicative_imputation(X)

  expect_equal(result, matrix(c(0.4, 0.6), nrow = 1))
})

test_that("multiplicative_imputation imputes zeros and preserves row sums", {
  X <- matrix(c(1, 0, 3), nrow = 1)

  result <- multiplicative_imputation(X, imp_factor = 0.01)

  expected <- matrix(c(
    0.25 * (1 - 0.01),
    0.01,
    0.75 * (1 - 0.01)
  ), nrow = 1)

  expect_equal(result, expected)
  expect_equal(sum(result), 1)
})

test_that("multiplicative_imputation returns NA for all-zero rows", {
  X <- matrix(c(0, 0, 0), nrow = 1)

  expect_warning(
    result <- multiplicative_imputation(X)
  )

  expect_true(all(is.na(result)))
})

test_that("multiplicative_imputation computes imp_factor automatically", {
  X <- matrix(c(0.2, 0, 0.8), nrow = 1)

  result <- multiplicative_imputation(X, imp_factor = NULL)

  expected_imp <- min(c(0.2, 0.8)) / 10

  expect_equal(result[1, 2], expected_imp)
  expect_equal(sum(result), 1)
})
