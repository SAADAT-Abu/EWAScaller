make_fake_result <- function(query, ok = TRUE) {
  if (!ok) {
    return(list(query = query, ok = FALSE, message = "not found", probes = NULL, associations = NULL))
  }
  probes <- data.frame(
    probeId = paste0("cg_", query), chrHg19 = "1", posHg19 = 100,
    cpgIsland = "Island", gene = "GENE1", stringsAsFactors = FALSE
  )
  associations <- data.frame(
    probeId = paste0("cg_", query), studyId = "ES00001", trait = "smoking",
    correlation = "neg", rank = 1L, pmid = 12345L, stringsAsFactors = FALSE
  )
  list(query = query, ok = TRUE, probes = probes, associations = associations)
}

test_that(".assemble_ewas_result builds associations, probes and failed tables", {
  results <- list(make_fake_result("a"), make_fake_result("b"), make_fake_result("c", ok = FALSE))
  res <- EWAScaller:::.assemble_ewas_result(results, query_type = "cpg")

  expect_s3_class(res, "ewas_result")
  expect_equal(res$query_type, "cpg")
  expect_equal(res$n_queries, 3)
  expect_equal(nrow(res$probes), 2)
  expect_equal(nrow(res$associations), 2)
  expect_equal(nrow(res$failed), 1)
  expect_equal(res$failed$query, "c")
})

test_that("ewas_result methods work on an assembled object", {
  results <- list(make_fake_result("a"), make_fake_result("b"))
  res <- EWAScaller:::.assemble_ewas_result(results, query_type = "cpg")

  expect_output(print(res), "ewas_result")
  expect_output(summary(res), "Top 1 trait")

  top <- top_traits(res, n = 5)
  expect_equal(top$trait, "smoking")
  expect_equal(top$n, 2)

  df <- as.data.frame(res)
  expect_identical(df, res$associations)
})

test_that("empty results still produce well-formed tables", {
  res <- EWAScaller:::.assemble_ewas_result(list(), query_type = "gene")
  expect_equal(nrow(res$associations), 0)
  expect_equal(nrow(res$probes), 0)
  expect_equal(names(res$associations), c("probeId", "studyId", "trait", "correlation", "rank", "pmid", "query"))
})
