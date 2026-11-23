test_that("probe ID validation accepts real 450K/EPIC formats, including ch probes", {
  ok <- c("cg05575921", "ch.4.37949012R", "ch.X.97651759F", "ch.17.4249285F", "rs1234567")
  expect_true(all(EWAScaller:::.is_valid_probe_id(ok)))
  expect_false(any(EWAScaller:::.is_valid_probe_id(c("AHRR", "ch123", "nv_1", "cg1_BC21"))))
})

test_that("EPIC v2 IDs are converted to v1 names with a message", {
  ids <- c("cg05575921_BC21", "cg05575921_TC22", "ch.4.37949012R_TO11", "cg00240195")
  expect_message(out <- EWAScaller:::.normalize_probe_ids(ids), "converted 3 EPIC v2.*1 replicate")
  expect_identical(out, c("cg05575921", "cg05575921", "ch.4.37949012R", "cg00240195"))
  expect_silent(EWAScaller:::.normalize_probe_ids(c("cg05575921", NA)))
  expect_silent(EWAScaller:::.normalize_probe_ids(ids, verbose = FALSE))
})

test_that("query_cpg queries EPIC v2 probes under their v1 names, once per merged replicate", {
  seen <- character()
  local_mocked_bindings(.query_one_probe = function(probe_id) {
    seen <<- c(seen, probe_id)
    list(query = probe_id, ok = FALSE, message = "mock", probes = NULL, associations = NULL)
  })
  expect_message(
    query_cpg(c("cg05575921_BC21", "cg05575921_TC22", "cg00240195"), workers = 1, delay = 0, progress = FALSE),
    "EPIC v2"
  )
  expect_setequal(seen, c("cg05575921", "cg00240195"))
})

test_that("background = 'EPICv2' runs against 850K and says so", {
  bk_sent <- NULL
  local_mocked_bindings(.ewas_post = function(path, body = list(), ...) {
    switch(path,
      "/checkJob" = list(ok = TRUE, value = "not exist"),
      "/checkInput" = { bk_sent <<- body$bk; list(ok = TRUE, value = "success") },
      "/analysis" = list(ok = TRUE, value = "ok"),
      "/jobStatus" = list(ok = TRUE, value = "finished"),
      "/getFile" = list(ok = TRUE, value = list(code = 0L, data = data.frame(trait = "x", p = "0.1"))),
      stop("unexpected ", path)
    )
  })
  v2 <- paste0(sprintf("cg%08d", 1:20), "_BC11")
  msgs <- character()
  enr <- withCallingHandlers(
    ewas_enrichment(v2, background = "EPICv2", progress = FALSE),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )
  expect_identical(bk_sent, "850K")
  expect_match(enr$background, "EPICv2")
  expect_identical(enr$probes, sprintf("cg%08d", 1:20))
  expect_true(any(grepl("no EPIC v2 background", msgs)))
  expect_true(any(grepl("converted 20 EPIC v2", msgs)))
})

test_that("word clouds with few words get a larger font", {
  few <- EWAScaller:::.wordcloud_plot(c("A", "B", "C"), 1:3, NULL)
  many <- EWAScaller:::.wordcloud_plot(paste0("w", 1:80), 1:80, NULL)
  expect_gt(ggplot2::ggplot_build(few)$plot$scales$get_scales("size")$palette(1),
            ggplot2::ggplot_build(many)$plot$scales$get_scales("size")$palette(1))
  fixed <- EWAScaller:::.wordcloud_plot(c("A", "B"), 1:2, NULL, max_size = 10)
  expect_equal(fixed$scales$get_scales("size")$palette(1), 10)
})
