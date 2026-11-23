probes20 <- sprintf("cg%08d", seq_len(20))

test_that("ewas_enrichment uses the 450K preset when background is left at its default", {
  seen_bk <- NULL
  local_mocked_bindings(
    .ewas_post = function(path, body = list(), times = 3L, as = "text") {
      if (path == "/checkJob") return(list(ok = TRUE, value = "not exist"))
      if (path == "/checkInput") {
        seen_bk <<- body$bk
        return(list(ok = TRUE, value = "success"))
      }
      if (path == "/analysis") return(list(ok = TRUE, value = "ok"))
      if (path == "/jobStatus") return(list(ok = TRUE, value = "finished"))
      if (path == "/getFile") {
        return(list(ok = TRUE, value = list(code = 0L, data = data.frame(
          trait = "smoking", p = "0.001", OR = "2.5", DMC = "3", stringsAsFactors = FALSE
        ))))
      }
      stop("unexpected path ", path)
    }
  )
  enr <- ewas_enrichment(probes20, progress = FALSE)
  expect_identical(seen_bk, "450K")
  expect_identical(enr$background, "450K")
  expect_identical(enr$job, EWAScaller:::.ewas_job_id(probes20, "450K"))
  expect_type(enr$results$trait$p, "double")
})

test_that("a server-side job failure is reported without waiting for the timeout", {
  local_mocked_bindings(
    .ewas_post = function(path, body = list(), times = 3L, as = "text") {
      switch(path,
        "/checkJob" = list(ok = TRUE, value = "exist"),
        "/analysis" = list(ok = TRUE, value = "ok"),
        "/jobStatus" = list(ok = TRUE, value = "error"),
        stop("unexpected path ", path)
      )
    }
  )
  elapsed <- system.time(
    enr <- ewas_enrichment(probes20, background = "850K", poll_interval = 5,
                           timeout = 60, progress = FALSE)
  )[["elapsed"]]
  expect_lt(elapsed, 5)
  expect_match(enr$failed$trait, "failed on the server")
})

test_that("plot_wordcloud does not change the user's random seed", {
  set.seed(42)
  before <- .Random.seed
  res <- EWAScaller:::.assemble_ewas_result(list(list(
    query = "a", ok = TRUE,
    probes = data.frame(probeId = "cg1", chrHg19 = "1", posHg19 = 1, cpgIsland = "Island",
                        gene = "G", stringsAsFactors = FALSE),
    associations = data.frame(probeId = "cg1", studyId = "ES1", trait = c("smoking", "aging"),
                              correlation = "neg", rank = 1L, pmid = 1L, stringsAsFactors = FALSE)
  )), query_type = "cpg")
  p <- plot_wordcloud(res)
  expect_s3_class(p, "ggplot")
  expect_identical(.Random.seed, before)
})
