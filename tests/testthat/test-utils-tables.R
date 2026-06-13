test_that("collapse and coercion helpers are deterministic", {
  expect_equal(forgeKI:::hdr_collapse_nonempty(c("B", "A", "", NA, "A")), "A;B")
  expect_true(is.na(forgeKI:::hdr_collapse_nonempty(c("", NA))))
  expect_equal(forgeKI:::hdr_chr0(c("x", NA)), c("x", ""))
  expect_equal(forgeKI:::hdr_na_chr(c("x", "")), c("x", NA_character_))
  expect_equal(forgeKI:::hdr_num(c("1", "bad")), c(1, NA))
  expect_equal(forgeKI:::hdr_bool(c("yes", "", NA), default = FALSE), c(TRUE, FALSE, FALSE))
})

test_that("table helpers preserve existing columns and safe joins", {
  df <- tibble::tibble(id = c("a", "b"), x = 1:2)
  out <- forgeKI:::hdr_add_missing_columns(df, list(x = 9, y = NA_character_))
  expect_equal(out$x, 1:2)
  expect_true("y" %in% names(out))
  expect_equal(forgeKI:::hdr_first_existing_col(out, c("z", "y")), "y")

  y <- tibble::tibble(id = c("a", "b"), score = c(10, 20))
  joined <- forgeKI:::hdr_left_join_safely(df, y, by = "id", y_cols = "score")
  expect_equal(joined$score, c(10, 20))
  expect_equal(forgeKI:::hdr_left_join_safely(df, y, by = "missing"), df)
})
