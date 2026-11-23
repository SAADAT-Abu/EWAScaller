#' Example CpG probe list
#'
#' A small example dataset of Illumina CpG probe identifiers, useful for
#' trying out [query_cpg()] without needing to supply your own probe list.
#' It contains only 3 distinct probes (with repeats), so it is too small
#' for [ewas_enrichment()], which requires 20-5000 distinct probes.
#'
#' @format A data frame with 30 rows and 1 variable:
#' \describe{
#'   \item{query_cpg}{Character. CpG probe identifier.}
#' }
#' @source Illustrative subset of probes with known EWAS Atlas trait
#'   associations.
#' @examples
#' data(example_cpgs)
#' head(example_cpgs)
"example_cpgs"
