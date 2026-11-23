fake_result <- function() {
  EWAScaller:::.assemble_ewas_result(list(list(
    query = "a", ok = TRUE,
    probes = data.frame(probeId = "cg1", chrHg19 = "1", posHg19 = 1, cpgIsland = "Island",
                        gene = "G", stringsAsFactors = FALSE),
    associations = data.frame(probeId = "cg1", studyId = "ES1",
                              trait = c("smoking", "a very long trait name that should be wrapped onto lines"),
                              correlation = c("neg", "pos"), rank = 1L, pmid = 1L,
                              stringsAsFactors = FALSE)
  )), query_type = "cpg")
}

fake_enrichment <- function() {
  EWAScaller:::new_ewas_enrichment(
    results = list(trait = data.frame(trait = c("smoking", "aging"), p = c(1e-5, 0.01),
                                      OR = c(2, 1.5), DMC = c(4, 2))),
    job = "abc", probes = "cg1", background = "450K", failed = list()
  )
}

test_that("title: NULL keeps the default, a string replaces it, NA removes it", {
  res <- fake_result()
  expect_identical(autoplot(res)$labels$title, "Top EWAS Atlas trait terms")
  expect_identical(autoplot(res, title = "Mine", subtitle = "Sub")$labels$title, "Mine")
  expect_identical(autoplot(res, subtitle = "Sub")$labels$subtitle, "Sub")
  expect_null(autoplot(res, title = NA)$labels$title)
  expect_null(autoplot(res, title = "")$labels$title)
})

test_that("theme accepts a function (with base_size) or a theme object, and rejects others", {
  res <- fake_result()
  p <- autoplot(res, theme = ggplot2::theme_bw, base_size = 18)
  expect_equal(p$theme$text$size, 18)
  p2 <- autoplot(res, theme = ggplot2::theme_classic(base_size = 9))
  expect_equal(p2$theme$text$size, 9)
  expect_error(autoplot(res, theme = "bw"), "must be a ggplot2 theme")
})

test_that("plot-specific theme tweaks survive a custom theme", {
  p <- autoplot(fake_result(), type = "direction", theme = ggplot2::theme_gray)
  expect_identical(p$theme$legend.position, "none")
})

test_that("colors are applied", {
  b <- ggplot2::ggplot_build(autoplot(fake_result(), colors = "#123456"))
  expect_true(all(b$data[[1]]$fill == "#123456"))
  b <- ggplot2::ggplot_build(autoplot(fake_result(), type = "direction", colors = c(pos = "red", neg = "blue")))
  expect_setequal(b$data[[1]]$fill, c("red", "blue"))
  b <- ggplot2::ggplot_build(autoplot(fake_enrichment(), style = "bar", colors = "purple"))
  expect_true(all(b$data[[1]]$fill == "purple"))
})

test_that("wrap controls label wrapping", {
  wrapped <- autoplot(fake_result(), wrap = 20)$data$label
  expect_true(any(grepl("\n", wrapped)))
  unwrapped <- autoplot(fake_result(), wrap = NULL)$data$label
  expect_false(any(grepl("\n", unwrapped)))
})

test_that("enrichment plots support dot and bar styles", {
  enr <- fake_enrichment()
  expect_s3_class(ggplot2::ggplot_build(autoplot(enr)), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(autoplot(enr, style = "bar")), "ggplot_built")
  expect_error(autoplot(enr, style = "pie"))
})

test_that("word clouds accept title, max_size and seed", {
  p <- plot_wordcloud(fake_result(), title = "Cloud", max_size = 30, seed = 7)
  expect_identical(p$labels$title, "Cloud")
})
