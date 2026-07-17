#' EWAS Atlas query result
#'
#' An `ewas_result` bundles the tables returned by [query_cpg()],
#' [query_gene()], or [query_region()].
#'
#' @section Fields:
#' \describe{
#'   \item{associations}{One row per probe-trait association, with columns
#'     `probeId`, `studyId`, `trait`, `correlation`, `rank`, `pmid`, and
#'     `query` (the input term that returned this row).}
#'   \item{probes}{One row per unique probe, with columns `probeId`,
#'     `chrHg19`, `posHg19`, `cpgIsland`, `gene`, and `query`.}
#'   \item{failed}{Input terms that could not be resolved, with an error
#'     `message` column.}
#'   \item{query_type}{One of `"cpg"`, `"gene"`, `"region"`.}
#'   \item{n_queries}{Number of input terms submitted.}
#' }
#'
#' @name ewas_result
#' @keywords internal
NULL

#' @keywords internal
#' @noRd
new_ewas_result <- function(associations, probes, failed, query_type, n_queries) {
  structure(
    list(
      associations = associations,
      probes = probes,
      failed = failed,
      query_type = query_type,
      n_queries = n_queries
    ),
    class = "ewas_result"
  )
}

#' @export
print.ewas_result <- function(x, ...) {
  cat(sprintf("<ewas_result> %s query\n", x$query_type))
  cat(sprintf(
    "  %d/%d input term(s) resolved | %d probe(s) | %d association(s) | %d unique trait(s)\n",
    x$n_queries - nrow(x$failed), x$n_queries, nrow(x$probes), nrow(x$associations),
    length(unique(x$associations$trait))
  ))
  if (nrow(x$failed) > 0) {
    cat(sprintf("  %d input(s) failed (see $failed)\n", nrow(x$failed)))
  }
  invisible(x)
}

#' Summarise an EWAS Atlas query result
#'
#' @param object An [ewas_result] object.
#' @param top_n Number of top traits to report.
#' @param ... Unused.
#' @return Invisibly, a list with summary statistics; also printed.
#' @export
summary.ewas_result <- function(object, top_n = 10L, ...) {
  assoc <- object$associations
  probes <- object$probes
  top <- top_traits(object, n = top_n)

  cat(sprintf("EWAS Atlas %s query summary\n", object$query_type))
  cat(sprintf("  Input terms:        %d\n", object$n_queries))
  cat(sprintf("  Resolved terms:     %d\n", object$n_queries - nrow(object$failed)))
  cat(sprintf("  Unique probes:      %d\n", nrow(probes)))
  cat(sprintf("  Associations:       %d\n", nrow(assoc)))
  cat(sprintf("  Unique traits:      %d\n", length(unique(assoc$trait))))
  cat(sprintf("  Unique studies:     %d\n", length(unique(stats::na.omit(assoc$studyId)))))
  if (nrow(top) > 0) {
    cat(sprintf("\nTop %d trait(s) by association count:\n", nrow(top)))
    for (i in seq_len(nrow(top))) {
      cat(sprintf("  %-45s %d\n", top$trait[i], top$n[i]))
    }
  }
  invisible(list(
    n_input = object$n_queries,
    n_resolved = object$n_queries - nrow(object$failed),
    n_probes = nrow(probes),
    n_associations = nrow(assoc),
    n_traits = length(unique(assoc$trait)),
    top_traits = top
  ))
}

#' @param row.names,optional Passed on for S3 consistency; unused.
#' @rdname ewas_result
#' @export
as.data.frame.ewas_result <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$associations
}

#' Rank traits by association frequency
#'
#' @param x An [ewas_result] or [ewas_enrichment] object.
#' @param n Number of top traits to return.
#' @param ... Passed to methods.
#' @return A data frame with columns `trait` and `n` (or `p`/`OR` for
#'   enrichment results), ordered from most to least significant/frequent.
#' @export
top_traits <- function(x, n = 10L, ...) UseMethod("top_traits")

#' @rdname top_traits
#' @export
top_traits.ewas_result <- function(x, n = 10L, ...) {
  assoc <- x$associations
  if (nrow(assoc) == 0) return(data.frame(trait = character(), n = integer()))
  tab <- as.data.frame(table(trait = assoc$trait), stringsAsFactors = FALSE)
  names(tab) <- c("trait", "n")
  tab <- tab[order(-tab$n, tab$trait), , drop = FALSE]
  rownames(tab) <- NULL
  utils::head(tab, n)
}
