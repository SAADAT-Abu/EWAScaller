# These tests contact the live EWAS Atlas web service. They are skipped on
# CRAN and whenever the service is unreachable (e.g. no internet access),
# so they never cause a CRAN check failure.

skip_if_offline <- function() {
  testthat::skip_on_cran()
  ok <- tryCatch({
    con <- url("https://ngdc.cncb.ac.cn/ewas/rest/probe?probeId=cg05575921", "rb")
    on.exit(close(con), add = TRUE)
    readBin(con, "raw", n = 1)
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(ok, "EWAS Atlas is not reachable")
}

test_that("query_cpg retrieves a known probe", {
  skip_if_offline()
  res <- query_cpg("cg05575921", workers = 1, delay = 0, progress = FALSE)
  expect_s3_class(res, "ewas_result")
  expect_true("cg05575921" %in% res$probes$probeId)
  expect_true(nrow(res$associations) > 0)
  expect_true("smoking" %in% res$associations$trait)
})

test_that("query_gene retrieves probes for a known gene", {
  skip_if_offline()
  res <- query_gene("AHRR", workers = 1, delay = 0, progress = FALSE)
  expect_s3_class(res, "ewas_result")
  expect_true(nrow(res$probes) > 0)
})

test_that("query_region retrieves probes in a known interval", {
  skip_if_offline()
  res <- query_region(chr = "5", start = 373000, end = 374000, workers = 1, delay = 0, progress = FALSE)
  expect_s3_class(res, "ewas_result")
  expect_true(nrow(res$probes) > 0)
})

test_that("ewas_enrichment runs trait enrichment end to end", {
  skip_if_offline()
  probes <- c(
    "cg05575921", "cg00240195", "cg16419584", "cg11903855", "cg21566642",
    "cg23771366", "cg19859270", "cg00382138", "cg25648203", "cg01940273",
    "cg23576855", "cg00401753", "cg00300637", "cg12803068", "cg09935388",
    "cg01899089", "cg04180046", "cg00050873", "cg02228185", "cg18147296"
  )
  enr <- ewas_enrichment(probes, background = "850K", types = "trait",
                          poll_interval = 3, timeout = 120, progress = FALSE)
  expect_s3_class(enr, "ewas_enrichment")
  expect_true(nrow(enr$results$trait) > 0)
  expect_true(all(c("trait", "p", "OR", "DMC") %in% names(enr$results$trait)))
})

test_that("ewas_eqtm and eqtm_samples return live CpG-gene expression data", {
  skip_if_offline()
  probes <- c(
    "cg05575921", "cg00240195", "cg16419584", "cg11903855", "cg21566642",
    "cg23771366", "cg19859270", "cg00382138", "cg25648203", "cg01940273",
    "cg23576855", "cg00401753", "cg00300637", "cg12803068", "cg09935388",
    "cg01899089", "cg04180046", "cg00050873", "cg02228185", "cg18147296"
  )
  eq <- ewas_eqtm(probes, background = "850K", poll_interval = 3, timeout = 300, progress = FALSE)
  expect_s3_class(eq, "ewas_eqtm")
  expect_true(nrow(eq$pairs) > 0)
  expect_type(eq$pairs$cor, "double")
  first <- eq$pairs[1, ]
  s <- eqtm_samples(first$probeId, first$gene_id, first$tissue)
  expect_true(nrow(s) > 10)
  expect_true(all(s$methylation >= 0 & s$methylation <= 1, na.rm = TRUE))
})
