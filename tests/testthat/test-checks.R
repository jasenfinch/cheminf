
test_that('error trhown if incorrect entries specified',{
  expect_error(checkEntries(tibble()))
})
