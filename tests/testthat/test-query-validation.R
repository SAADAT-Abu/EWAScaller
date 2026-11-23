test_that("query_cpg rejects empty input", {
  expect_error(query_cpg(character(0)), "non-missing probe identifier")
  expect_error(query_cpg(c(NA, "")), "non-missing probe identifier")
})

test_that("query_cpg warns on and drops malformed identifiers", {
  expect_error(
    expect_warning(query_cpg("not_a_probe"), "do not look like probe identifiers"),
    "No valid probe identifiers"
  )
})

test_that("query_gene rejects empty input", {
  expect_error(query_gene(character(0)), "non-missing gene symbol")
})

test_that("query_region requires start/end or a data frame", {
  expect_error(query_region("1"), "required when")
  expect_error(
    query_region(data.frame(chr = "1", start = 100, end = 50)),
    "greater than or equal to"
  )
  expect_error(query_region("1", start = NA, end = 50), "must not be missing")
  expect_error(
    query_region(data.frame(chr = NA, start = 1, end = 50)),
    "must not be missing"
  )
})

test_that("ewas_enrichment enforces probe-count bounds", {
  expect_error(ewas_enrichment(paste0("cg", 1:5)), "at least 20 valid")
  expect_error(ewas_enrichment(paste0("cg", seq_len(5001))), "at most 5000")
})

test_that("ewas_enrichment validates custom background size", {
  probes <- sprintf("cg%08d", seq_len(25))
  expect_error(ewas_enrichment(probes, background = c("cg00000001", "cg00000002")), "at least 20 probe")
})
