probes20 <- sprintf("cg%08d", seq_len(20))

fake_pairs <- data.frame(
  p = c("1e-10", "0.001", "0.02", "1e-5"),
  gene = c("ENSG01", "ENSG02", "ENSG01", "ENSG03"),
  cor = c(-0.4, 0.3, 0.2, 0.5),
  gene_symbol = c("AHRR", "GFI1", "AHRR", "MYO1G"),
  location = "body",
  tissue = c("brain", "colon", "colon", "kidney"),
  source = "TCGA",
  probe = c("cg00000001", "cg00000002", "cg00000001", "cg00000003"),
  stringsAsFactors = FALSE
)

mock_toolkit <- function(data = fake_pairs) {
  function(path, body = list(), times = 3L, as = "text") {
    switch(path,
      "/checkJob" = list(ok = TRUE, value = "exist"),
      "/analysis" = list(ok = TRUE, value = "ok"),
      "/jobStatus" = list(ok = TRUE, value = "finished"),
      "/getFile" = {
        stopifnot(identical(body$file, "expression_correlation.tsv"))
        list(ok = TRUE, value = list(code = 0L, data = data))
      },
      stop("unexpected path ", path)
    )
  }
}

fake_eqtm <- function() {
  EWAScaller:::new_ewas_eqtm(EWAScaller:::.tidy_eqtm_pairs(fake_pairs),
                             job = "abc", probes = probes20, background = "850K")
}

test_that("ewas_eqtm returns tidy, typed, p-sorted pairs", {
  local_mocked_bindings(.ewas_post = mock_toolkit())
  eq <- ewas_eqtm(probes20, background = "850K", progress = FALSE)
  expect_s3_class(eq, "ewas_eqtm")
  expect_named(eq$pairs, c("probeId", "gene_symbol", "gene_id", "location", "cor", "p", "tissue", "source"))
  expect_type(eq$pairs$p, "double")
  expect_false(is.unsorted(eq$pairs$p))
  expect_identical(eq$pairs$gene_id[1], "ENSG01")
  expect_identical(as.data.frame(eq), eq$pairs)
  expect_output(print(eq), "4 \\(3 gene\\(s\\), 3 tissue\\(s\\)\\)")
  expect_output(s <- summary(eq), "Pairs per tissue")
  expect_equal(s$by_tissue$tissue[1], "colon")
})

test_that("ewas_eqtm handles an empty result and reports failures", {
  local_mocked_bindings(.ewas_post = mock_toolkit(data = list()))
  eq <- ewas_eqtm(probes20, background = "850K", progress = FALSE)
  expect_equal(nrow(eq$pairs), 0)
  expect_error(autoplot(eq), "No eQTM pairs")

  local_mocked_bindings(.ewas_post = function(path, ...) {
    if (path == "/jobStatus") return(list(ok = TRUE, value = "error"))
    list(ok = TRUE, value = "exist")
  })
  expect_error(ewas_eqtm(probes20, background = "850K", progress = FALSE), "eQTM analysis failed")
})

test_that("ewas_eqtm enforces the toolkit's probe-count limits", {
  expect_error(ewas_eqtm(probes20[1:5]), "at least 20 valid")
})

test_that("eqtm_samples validates input and parses the sample matrix", {
  expect_error(eqtm_samples(c("cg1", "cg2"), "ENSG01", "brain"), "single non-empty string")
  local_mocked_bindings(.ewas_toolkit_get = function(path, query, ...) {
    expect_identical(path, "/expression")
    list(ok = TRUE, value = matrix(c(0.1, 0.5, 10, 12), ncol = 2))
  })
  s <- eqtm_samples("cg00000001", "ENSG01", "brain")
  expect_equal(s, data.frame(methylation = c(0.1, 0.5), expression = c(10, 12)))
})

test_that("every eQTM plot type builds", {
  eq <- fake_eqtm()
  for (type in c("heatmap", "volcano", "tissue")) {
    p <- autoplot(eq, type = type)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
  expect_s3_class(plot_wordcloud(eq), "ggplot")

  local_mocked_bindings(.ewas_toolkit_get = function(...) {
    list(ok = TRUE, value = matrix(c(0.1, 0.5, 0.9, 10, 12, 11), ncol = 2))
  })
  p <- autoplot(eq, type = "scatter", probe = "cg00000001", gene = "AHRR", tissue = "brain")
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle, "r = -0.400")
})

test_that("scatter plots explain how to pick a single pair", {
  eq <- fake_eqtm()
  expect_error(autoplot(eq, type = "scatter"), "`probe` and `gene` are required")
  expect_error(autoplot(eq, type = "scatter", probe = "cg00000001", gene = "AHRR"),
               "several tissues \\(brain, colon\\)")
  expect_error(autoplot(eq, type = "scatter", probe = "cg00000001", gene = "NOPE"), "No eQTM pair")
})

test_that("tissue filtering limits what is plotted", {
  p <- autoplot(fake_eqtm(), type = "volcano", tissue = "colon", label_top = 0)
  expect_setequal(unique(p$data$tissue), "colon")
  expect_error(autoplot(fake_eqtm(), tissue = "lung"), "No eQTM pairs")
})
