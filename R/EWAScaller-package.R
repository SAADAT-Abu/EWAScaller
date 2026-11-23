#' EWAScaller: Query and Analyse the EWAS Atlas Database
#'
#' Tools to query the EWAS Atlas (<https://ngdc.cncb.ac.cn/ewas/>) by CpG
#' probe, gene symbol, or genomic region, and to run enrichment analyses on
#' probe sets. Query functions dispatch requests concurrently across a small
#' worker pool with configurable pacing to stay within reasonable request
#' rates. Results are returned as tidy data frames with summary and plotting
#' helpers, including trait-term word clouds.
#'
#' @section Main functions:
#' \itemize{
#'   \item [query_cpg()], [query_gene()], [query_region()] -- retrieve
#'     probe-trait associations from the EWAS Atlas.
#'   \item [ewas_enrichment()] -- run trait, gene ontology, KEGG pathway, or
#'     genomic location enrichment on a probe set.
#'   \item [ewas_eqtm()], [eqtm_samples()] -- find CpG-gene expression
#'     quantitative trait methylation (eQTM) pairs and their per-sample
#'     data.
#'   \item [plot_wordcloud()], [autoplot.ewas_result()],
#'     [autoplot.ewas_enrichment()], [autoplot.ewas_eqtm()] -- visualise
#'     results.
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom ggplot2 .data
#' @importFrom stats reorder setNames
#' @importFrom utils head
## usethis namespace: end
NULL

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
