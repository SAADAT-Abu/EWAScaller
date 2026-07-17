test_that(".is_valid_probe_id recognises Illumina-style identifiers", {
  expect_true(all(EWAScaller:::.is_valid_probe_id(c("cg05575921", "cg11903855", "rs123456"))))
  expect_false(EWAScaller:::.is_valid_probe_id("not_a_probe"))
  expect_equal(
    EWAScaller:::.is_valid_probe_id(c("cg05575921", "AHRR", "cg00240195")),
    c(TRUE, FALSE, TRUE)
  )
})

test_that(".ewas_job_id is order-independent and background-sensitive", {
  id1 <- EWAScaller:::.ewas_job_id(c("cg1", "cg2", "cg3"), "450K")
  id2 <- EWAScaller:::.ewas_job_id(c("cg3", "cg1", "cg2"), "450K")
  id3 <- EWAScaller:::.ewas_job_id(c("cg1", "cg2", "cg3"), "850K")
  expect_identical(id1, id2)
  expect_false(identical(id1, id3))
  expect_match(id1, "^[0-9a-f]{32}$")
})

test_that(".collapse_unique drops missing values and de-duplicates", {
  expect_identical(EWAScaller:::.collapse_unique(c("A", "A", "B", NA, "")), "A; B")
  expect_identical(EWAScaller:::.collapse_unique(c(NA, "")), NA_character_)
})

test_that(".run_parallel preserves input order sequentially and in parallel", {
  out_seq <- EWAScaller:::.run_parallel(as.list(1:5), function(x) x * 2, workers = 1, delay = 0, progress = FALSE)
  expect_equal(unlist(out_seq), c(2, 4, 6, 8, 10))

  out_par <- EWAScaller:::.run_parallel(as.list(1:5), function(x) x * 2, workers = 2, delay = 0, progress = FALSE)
  expect_equal(unlist(out_par), c(2, 4, 6, 8, 10))
})

test_that(".coerce_numeric_cols only converts fully-numeric character columns", {
  df <- data.frame(
    p = c("0.01", "0.02"),
    trait = c("smoking", "aging"),
    mixed = c("1", "not_numeric"),
    stringsAsFactors = FALSE
  )
  out <- EWAScaller:::.coerce_numeric_cols(df)
  expect_type(out$p, "double")
  expect_type(out$trait, "character")
  expect_type(out$mixed, "character")
})
