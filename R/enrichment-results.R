#' EWAS Atlas enrichment result
#'
#' An `ewas_enrichment` bundles the analyses run by [ewas_enrichment()].
#'
#' @section Fields:
#' \describe{
#'   \item{results}{Named list of data frames, one per requested analysis
#'     type (e.g. `results$trait`).}
#'   \item{job}{The EWAS Atlas toolkit job identifier.}
#'   \item{probes}{The probe identifiers submitted.}
#'   \item{background}{The background set used (`"450K"`, `"850K"`,
#'     `"EPICv2 (approximated by 850K)"`, or `"other"`).}
#'   \item{failed}{Named list of analysis types that could not be
#'     retrieved, with an explanatory message.}
#'   \item{excluded}{Input probes dropped because they are not in the
#'     background (see the `background` argument of [ewas_enrichment()]).}
#' }
#'
#' @return An object of class `ewas_enrichment`: a named list with the
#'   `results`, `job`, `probes`, `background`, `failed`, and `excluded` elements
#'   described above. `summary.ewas_enrichment()` invisibly returns the
#'   object it was called on (after printing a summary), and
#'   `as.data.frame.ewas_enrichment()` returns the requested analysis's
#'   result data frame on its own.
#' @name ewas_enrichment-class
#' @keywords internal
NULL

#' @keywords internal
#' @noRd
new_ewas_enrichment <- function(results, job, probes, background, failed, excluded = character()) {
  structure(
    list(results = results, job = job, probes = probes, background = background,
         failed = failed, excluded = excluded),
    class = "ewas_enrichment"
  )
}

#' @export
print.ewas_enrichment <- function(x, ...) {
  cat("<ewas_enrichment>\n")
  cat(sprintf("  Job:        %s\n", x$job))
  cat(sprintf("  Background: %s\n", x$background))
  cat(sprintf("  Probes:     %d\n", length(x$probes)))
  if (length(x$excluded) > 0) {
    cat(sprintf("  Excluded:   %d probe(s) not in the background (see $excluded)\n", length(x$excluded)))
  }
  for (nm in names(x$results)) {
    cat(sprintf("  %-18s %d row(s)\n", nm, nrow(x$results[[nm]])))
  }
  if (length(x$failed) > 0) {
    cat(sprintf("  Failed:     %s\n", paste(names(x$failed), collapse = ", ")))
  }
  invisible(x)
}

#' @param object An [ewas_enrichment-class] object.
#' @param top_n Number of top rows to report per analysis type.
#' @param ... Unused.
#' @rdname ewas_enrichment-class
#' @export
summary.ewas_enrichment <- function(object, top_n = 10L, ...) {
  cat("EWAS Atlas enrichment summary\n")
  cat(sprintf("  Job:        %s\n", object$job))
  cat(sprintf("  Background: %s\n", object$background))
  cat(sprintf("  Probes:     %d\n", length(object$probes)))
  if ("trait" %in% names(object$results)) {
    cat(sprintf("\nTop %d enriched trait(s):\n", top_n))
    print(top_traits(object, n = top_n))
  }
  invisible(object)
}

#' @rdname top_traits
#' @export
top_traits.ewas_enrichment <- function(x, n = 10L, ...) {
  df <- x$results[["trait"]]
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(trait = character(), p = numeric(), OR = numeric(), DMC = integer()))
  }
  df <- df[order(df$p), , drop = FALSE]
  cols <- intersect(c("trait", "p", "OR", "DMC", "background"), names(df))
  rownames(df) <- NULL
  utils::head(df[, cols, drop = FALSE], n)
}

#' @param x An [ewas_enrichment-class] object.
#' @param row.names,optional Unused; present for S3 consistency.
#' @param type Which analysis result to extract (default `"trait"`).
#' @rdname ewas_enrichment-class
#' @export
as.data.frame.ewas_enrichment <- function(x, row.names = NULL, optional = FALSE, type = "trait", ...) {
  x$results[[type]]
}
