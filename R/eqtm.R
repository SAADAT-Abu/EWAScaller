#' Find expression quantitative trait methylation (eQTM) pairs for CpG probes
#'
#' Runs the EWAS Atlas toolkit's "Expression Regulation" analysis, which
#' reports CpG-gene pairs whose DNA methylation is correlated with gene
#' expression (expression quantitative trait methylation, eQTM). The
#' correlations were computed by EWAS Atlas from matched methylation and
#' RNA-seq data from The Cancer Genome Atlas (TCGA) across several tissues.
#'
#' Like [ewas_enrichment()], the analysis runs as an asynchronous toolkit
#' job and requires between 20 and 5000 probes. A probe set that has
#' already been submitted (for example by an earlier [ewas_enrichment()]
#' call with the same background) reuses the existing job.
#'
#' @inheritParams ewas_enrichment
#' @param probes Character vector of CpG probe identifiers (20-5000).
#'
#' @return An [ewas_eqtm-class] object whose `pairs` element is a data frame
#'   with one row per CpG-gene-tissue correlation and columns `probeId`,
#'   `gene_symbol`, `gene_id` (Ensembl), `location` (position of the CpG
#'   relative to the gene), `cor` (Pearson correlation), `p`, `tissue`, and
#'   `source` (TCGA project(s)).
#'
#' @seealso [eqtm_samples()] for the per-sample data behind one pair, and
#'   [autoplot.ewas_eqtm()] for plots.
#'
#' @examples
#' \donttest{
#' probes <- c(
#'   "cg05575921", "cg00240195", "cg16419584", "cg11903855", "cg21566642",
#'   "cg23771366", "cg19859270", "cg00382138", "cg25648203", "cg01940273",
#'   "cg23576855", "cg00401753", "cg00300637", "cg12803068", "cg09935388",
#'   "cg01899089", "cg04180046", "cg00050873", "cg02228185", "cg18147296"
#' )
#' eq <- ewas_eqtm(probes, background = "850K")
#' summary(eq)
#' head(as.data.frame(eq))
#' }
#'
#' @export
ewas_eqtm <- function(probes, background = c("450K", "850K"),
                      poll_interval = 5, timeout = 300, progress = interactive()) {
  probes <- .prepare_toolkit_probes(probes)
  bk <- .resolve_background(background)
  job <- .submit_toolkit_job(probes, bk$string)

  if (isTRUE(progress)) message("EWAScaller: running eQTM (expression regulation) analysis...")
  out <- .run_toolkit_analysis(job, "expression", "expression_correlation.tsv",
                               poll_interval = poll_interval, timeout = timeout, progress = progress)
  if (!isTRUE(out$ok)) {
    stop(sprintf("eQTM analysis failed: %s", out$message), call. = FALSE)
  }

  new_ewas_eqtm(.tidy_eqtm_pairs(out$data), job = job, probes = probes, background = bk$label)
}

#' @keywords internal
#' @noRd
.empty_eqtm_pairs <- function() {
  data.frame(
    probeId = character(), gene_symbol = character(), gene_id = character(),
    location = character(), cor = numeric(), p = numeric(),
    tissue = character(), source = character(), stringsAsFactors = FALSE
  )
}

#' Rename and order the columns of the toolkit's expression_correlation table
#' @keywords internal
#' @noRd
.tidy_eqtm_pairs <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(.empty_eqtm_pairs())
  col <- function(name, fallback) if (name %in% names(df)) df[[name]] else rep(fallback, nrow(df))
  out <- data.frame(
    probeId = as.character(col("probe", NA_character_)),
    gene_symbol = as.character(col("gene_symbol", NA_character_)),
    gene_id = as.character(col("gene", NA_character_)),
    location = as.character(col("location", NA_character_)),
    cor = suppressWarnings(as.numeric(col("cor", NA_real_))),
    p = suppressWarnings(as.numeric(col("p", NA_real_))),
    tissue = as.character(col("tissue", NA_character_)),
    source = as.character(col("source", NA_character_)),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$p, -abs(out$cor)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' EWAS Atlas eQTM result
#'
#' An `ewas_eqtm` holds the CpG-gene expression correlations returned by
#' [ewas_eqtm()].
#'
#' @section Fields:
#' \describe{
#'   \item{pairs}{Data frame of CpG-gene-tissue correlations, sorted by
#'     `p`; see [ewas_eqtm()] for its columns.}
#'   \item{job}{The EWAS Atlas toolkit job identifier.}
#'   \item{probes}{The probe identifiers submitted.}
#'   \item{background}{The background set used (`"450K"`, `"850K"`, or
#'     `"other"`).}
#' }
#'
#' @return An object of class `ewas_eqtm`: a named list with the `pairs`,
#'   `job`, `probes`, and `background` elements described above.
#'   `summary.ewas_eqtm()` invisibly returns a list of summary tables
#'   (after printing them), and `as.data.frame.ewas_eqtm()` returns the
#'   `pairs` data frame.
#' @name ewas_eqtm-class
#' @keywords internal
NULL

#' @keywords internal
#' @noRd
new_ewas_eqtm <- function(pairs, job, probes, background) {
  structure(
    list(pairs = pairs, job = job, probes = probes, background = background),
    class = "ewas_eqtm"
  )
}

#' @export
print.ewas_eqtm <- function(x, ...) {
  pr <- x$pairs
  cat("<ewas_eqtm>\n")
  cat(sprintf("  Job:        %s\n", x$job))
  cat(sprintf("  Background: %s\n", x$background))
  cat(sprintf("  Probes:     %d submitted, %d with an eQTM\n",
              length(x$probes), length(unique(pr$probeId))))
  cat(sprintf("  Pairs:      %d (%d gene(s), %d tissue(s))\n",
              nrow(pr), length(unique(pr$gene_id)), length(unique(pr$tissue))))
  invisible(x)
}

#' @param object An [ewas_eqtm-class] object.
#' @param top_n Number of strongest pairs to list.
#' @param ... Unused.
#' @rdname ewas_eqtm-class
#' @export
summary.ewas_eqtm <- function(object, top_n = 10L, ...) {
  pr <- object$pairs
  cat("EWAS Atlas eQTM summary\n")
  cat(sprintf("  Background:        %s\n", object$background))
  cat(sprintf("  Probes submitted:  %d\n", length(object$probes)))
  cat(sprintf("  Probes with eQTM:  %d\n", length(unique(pr$probeId))))
  cat(sprintf("  Genes:             %d\n", length(unique(pr$gene_id))))
  cat(sprintf("  Pairs:             %d CpG-gene-tissue (%d negative, %d positive)\n",
              nrow(pr), sum(pr$cor < 0, na.rm = TRUE), sum(pr$cor > 0, na.rm = TRUE)))

  by_tissue <- if (nrow(pr) > 0) {
    tab <- as.data.frame(table(tissue = pr$tissue), stringsAsFactors = FALSE)
    names(tab) <- c("tissue", "n_pairs")
    tab[order(-tab$n_pairs, tab$tissue), , drop = FALSE]
  } else {
    data.frame(tissue = character(), n_pairs = integer())
  }
  rownames(by_tissue) <- NULL
  top <- utils::head(pr[, c("probeId", "gene_symbol", "tissue", "cor", "p"), drop = FALSE], top_n)

  if (nrow(by_tissue) > 0) {
    cat("\nPairs per tissue:\n")
    print(by_tissue, row.names = FALSE)
  }
  if (nrow(top) > 0) {
    cat(sprintf("\nTop %d pair(s) by p-value:\n", nrow(top)))
    print(top, row.names = FALSE)
  }
  invisible(list(by_tissue = by_tissue, top_pairs = top))
}

#' @param x An [ewas_eqtm-class] object.
#' @param row.names,optional Unused; present for S3 consistency.
#' @rdname ewas_eqtm-class
#' @export
as.data.frame.ewas_eqtm <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$pairs
}

#' Per-sample methylation and expression values for one eQTM pair
#'
#' Downloads the matched DNA methylation (beta value) and gene expression
#' (TPM) measurements behind one CpG-gene-tissue correlation reported by
#' [ewas_eqtm()].
#'
#' @param probe CpG probe identifier, e.g. `"cg05575921"` (EPIC v2 names
#'   such as `"cg05575921_BC21"` are accepted).
#' @param gene Ensembl gene identifier, e.g. `"ENSG00000063438"` (the
#'   `gene_id` column of an [ewas_eqtm()] result).
#' @param tissue Tissue name as reported in the `tissue` column, e.g.
#'   `"colon"`.
#'
#' @return A data frame with one row per sample and columns `methylation`
#'   and `expression`.
#'
#' @examples
#' \donttest{
#' s <- eqtm_samples("cg05575921", "ENSG00000063438", "brain")
#' head(s)
#' }
#'
#' @export
eqtm_samples <- function(probe, gene, tissue) {
  for (arg in c("probe", "gene", "tissue")) {
    v <- get(arg)
    if (!is.character(v) || length(v) != 1L || is.na(v) || !nzchar(v)) {
      stop(sprintf("`%s` must be a single non-empty string.", arg), call. = FALSE)
    }
  }
  probe <- .normalize_probe_ids(probe, verbose = FALSE)
  res <- .ewas_toolkit_get("/expression", list(probe = probe, gene = gene, tissue = tissue))
  if (!isTRUE(res$ok)) {
    stop(sprintf("Could not retrieve eQTM samples: %s", res$message), call. = FALSE)
  }
  .tidy_eqtm_samples(res$value)
}

#' @keywords internal
#' @noRd
.tidy_eqtm_samples <- function(value) {
  empty <- data.frame(methylation = numeric(), expression = numeric())
  if (is.null(value) || length(value) == 0) return(empty)
  m <- if (is.matrix(value)) value else do.call(rbind, lapply(value, unlist))
  if (is.null(m) || ncol(m) < 2) return(empty)
  data.frame(
    methylation = suppressWarnings(as.numeric(m[, 1])),
    expression = suppressWarnings(as.numeric(m[, 2]))
  )
}
