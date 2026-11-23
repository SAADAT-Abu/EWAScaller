#' Query the EWAS Atlas by CpG probe identifier
#'
#' Retrieves trait associations and annotation for one or more Illumina CpG
#' probe identifiers (e.g. `"cg05575921"`) from the EWAS Atlas.
#'
#' @param cpgs Character vector of CpG probe identifiers. Duplicates and
#'   missing values are dropped before querying.
#' @param workers Integer; maximum number of concurrent requests. Defaults
#'   to `2`, a conservative value intended to avoid overloading the remote
#'   service.
#' @param delay Numeric; seconds to pause between batches of `workers`
#'   requests.
#' @param progress Logical; print a simple progress message as batches
#'   complete.
#'
#' @return An [ewas_result] object.
#'
#' @examples
#' \donttest{
#' res <- query_cpg(c("cg05575921", "cg11903855"), workers = 1, delay = 1)
#' summary(res)
#' }
#'
#' @export
query_cpg <- function(cpgs, workers = 2L, delay = 1, progress = interactive()) {
  cpgs <- as.character(cpgs)
  cpgs <- unique(cpgs[!is.na(cpgs) & nzchar(cpgs)])
  if (length(cpgs) == 0) {
    stop("`cpgs` must contain at least one non-missing probe identifier.", call. = FALSE)
  }
  invalid <- cpgs[!.is_valid_probe_id(cpgs)]
  if (length(invalid) > 0) {
    warning(sprintf(
      "Skipping %d input(s) that do not look like probe identifiers: %s",
      length(invalid), paste(utils::head(invalid, 5), collapse = ", ")
    ), call. = FALSE)
    cpgs <- setdiff(cpgs, invalid)
  }
  if (length(cpgs) == 0) stop("No valid probe identifiers to query.", call. = FALSE)

  results <- .run_parallel(as.list(cpgs), .query_one_probe, workers = workers, delay = delay, progress = progress)
  .assemble_ewas_result(results, query_type = "cpg")
}

#' @keywords internal
#' @noRd
.query_one_probe <- function(probe_id) {
  res <- .ewas_get("/probe", list(probeId = probe_id))
  if (!isTRUE(res$ok)) {
    return(list(query = probe_id, ok = FALSE, message = res$message, probes = NULL, associations = NULL))
  }
  parsed <- .rows_from_probe_data(res$value)
  list(query = probe_id, ok = TRUE, probes = parsed$probes, associations = parsed$associations)
}

#' Query the EWAS Atlas by gene symbol
#'
#' Retrieves probes and trait associations linked to one or more gene
#' symbols (e.g. `"AHRR"`) from the EWAS Atlas.
#'
#' @param genes Character vector of gene symbols.
#' @inheritParams query_cpg
#'
#' @return An [ewas_result] object.
#'
#' @examples
#' \donttest{
#' res <- query_gene("AHRR", workers = 1, delay = 1)
#' summary(res)
#' }
#'
#' @export
query_gene <- function(genes, workers = 2L, delay = 1, progress = interactive()) {
  genes <- as.character(genes)
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) == 0) {
    stop("`genes` must contain at least one non-missing gene symbol.", call. = FALSE)
  }

  results <- .run_parallel(as.list(genes), .query_one_gene, workers = workers, delay = delay, progress = progress)
  .assemble_ewas_result(results, query_type = "gene")
}

#' @keywords internal
#' @noRd
.query_one_gene <- function(gene_symbol) {
  res <- .ewas_get("/gene", list(geneSymbol = gene_symbol))
  if (!isTRUE(res$ok)) {
    return(list(query = gene_symbol, ok = FALSE, message = res$message, probes = NULL, associations = NULL))
  }
  data <- res$value
  probe_list <- if (is.list(data) && !is.data.frame(data) && !is.null(data$probeList)) data$probeList else data
  parsed <- .rows_from_probe_data(probe_list)
  list(query = gene_symbol, ok = TRUE, probes = parsed$probes, associations = parsed$associations)
}

#' Query the EWAS Atlas by genomic region
#'
#' Retrieves probes and trait associations located within one or more
#' genomic intervals (hg19 coordinates) from the EWAS Atlas.
#'
#' @param chr Chromosome (e.g. `"1"`, `"X"`), or a data frame with columns
#'   `chr`, `start`, `end` giving one or more regions, in which case `start`
#'   and `end` are ignored.
#' @param start Start position (hg19, 1-based).
#' @param end End position (hg19, 1-based).
#' @inheritParams query_cpg
#'
#' @return An [ewas_result] object.
#'
#' @examples
#' \donttest{
#' res <- query_region(chr = "5", start = 373000, end = 374000, workers = 1)
#' summary(res)
#' }
#'
#' @export
query_region <- function(chr, start = NULL, end = NULL, workers = 2L, delay = 1, progress = interactive()) {
  if (is.data.frame(chr)) {
    df <- chr
    if (!all(c("chr", "start", "end") %in% names(df))) {
      stop("`chr` data frame must have columns 'chr', 'start', 'end'.", call. = FALSE)
    }
  } else {
    if (is.null(start) || is.null(end)) {
      stop("`start` and `end` are required when `chr` is not a data frame.", call. = FALSE)
    }
    df <- data.frame(chr = chr, start = start, end = end, stringsAsFactors = FALSE)
  }
  df$chr <- gsub("^chr", "", as.character(df$chr), ignore.case = TRUE)
  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)
  if (anyNA(df$chr) || anyNA(df$start) || anyNA(df$end) || any(!nzchar(df$chr))) {
    stop("`chr`, `start`, and `end` must not be missing.", call. = FALSE)
  }
  if (any(df$end < df$start)) stop("`end` must be greater than or equal to `start` for every region.", call. = FALSE)

  regions <- lapply(split(df, seq_len(nrow(df))), as.list)
  results <- .run_parallel(regions, .query_one_region, workers = workers, delay = delay, progress = progress)
  .assemble_ewas_result(results, query_type = "region")
}

#' @keywords internal
#' @noRd
.query_one_region <- function(region) {
  label <- sprintf("chr%s:%s-%s", region$chr, region$start, region$end)
  res <- .ewas_get("/pos", list(chr = region$chr, start = region$start, end = region$end))
  if (!isTRUE(res$ok)) {
    return(list(query = label, ok = FALSE, message = res$message, probes = NULL, associations = NULL))
  }
  parsed <- .rows_from_probe_data(res$value)
  list(query = label, ok = TRUE, probes = parsed$probes, associations = parsed$associations)
}

#' @keywords internal
#' @noRd
.assemble_ewas_result <- function(results, query_type) {
  ok_results <- Filter(function(r) isTRUE(r$ok), results)
  failed_results <- Filter(function(r) !isTRUE(r$ok), results)

  probes_tagged <- lapply(ok_results, function(r) {
    if (is.null(r$probes)) return(NULL)
    d <- r$probes
    d$query <- r$query
    d
  })
  assoc_tagged <- lapply(ok_results, function(r) {
    if (is.null(r$associations)) return(NULL)
    d <- r$associations
    d$query <- r$query
    d
  })

  probes_df <- do.call(rbind, Filter(Negate(is.null), probes_tagged))
  assoc_df <- do.call(rbind, Filter(Negate(is.null), assoc_tagged))

  if (!is.null(probes_df)) {
    probes_df <- unique(probes_df)
    rownames(probes_df) <- NULL
  }
  if (!is.null(assoc_df)) rownames(assoc_df) <- NULL

  failed_df <- if (length(failed_results) > 0) {
    data.frame(
      query = vapply(failed_results, function(r) as.character(r$query), character(1)),
      message = vapply(failed_results, function(r) as.character(r$message %||% "no results"), character(1)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(query = character(), message = character(), stringsAsFactors = FALSE)
  }

  new_ewas_result(
    associations = if (is.null(assoc_df)) .empty_associations() else assoc_df,
    probes = if (is.null(probes_df)) .empty_probes() else probes_df,
    failed = failed_df,
    query_type = query_type,
    n_queries = length(results)
  )
}
