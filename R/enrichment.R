.enrichment_files <- c(
  trait = "trait.tsv",
  genomic_location = "genomic_location.tsv",
  gene_ontology = "gene_ontology.tsv",
  kegg = "kegg.tsv"
)

#' Run EWAS Atlas enrichment analysis on a set of CpG probes
#'
#' Submits a probe set to the EWAS Atlas toolkit and retrieves one or more
#' enrichment analyses run against it: EWAS trait-term enrichment, genomic
#' location enrichment, Gene Ontology enrichment, KEGG pathway enrichment,
#' tissue-specificity enrichment, and gene-expression correlation.
#'
#' The EWAS Atlas toolkit requires between 20 and 5000 probes per submission
#' and computes results asynchronously; this function submits the job and
#' polls for completion.
#'
#' @param probes Character vector of CpG probe identifiers (20-5000).
#' @param background One of `"450K"`, `"850K"`, or a character vector of
#'   probe identifiers to use as a custom background set.
#' @param types Character vector of analyses to run: any of `"trait"`
#'   (EWAS Atlas trait-term enrichment), `"genomic_location"` (enrichment
#'   across gene/CpG-island-relative locations), `"gene_ontology"` (GO term
#'   enrichment), or `"kegg"` (KEGG pathway enrichment). Defaults to
#'   `"trait"`.
#' @param poll_interval Seconds between job-status polls.
#' @param timeout Maximum seconds to wait for each analysis to finish.
#' @param progress Logical; print progress messages while polling.
#'
#' @return An [ewas_enrichment] object.
#'
#' @examples
#' \donttest{
#' probes <- c(
#'   "cg05575921", "cg00240195", "cg16419584", "cg11903855", "cg21566642",
#'   "cg23771366", "cg19859270", "cg00382138", "cg25648203", "cg01940273",
#'   "cg23576855", "cg00401753", "cg00300637", "cg12803068", "cg09935388",
#'   "cg01899089", "cg04180046", "cg00050873", "cg02228185", "cg18147296"
#' )
#' enr <- ewas_enrichment(probes, background = "850K", types = "trait")
#' top_traits(enr)
#' }
#'
#' @export
ewas_enrichment <- function(probes, background = c("450K", "850K"), types = "trait",
                             poll_interval = 5, timeout = 300, progress = interactive()) {
  probes <- unique(as.character(probes))
  probes <- probes[!is.na(probes) & nzchar(probes)]
  invalid <- probes[!.is_valid_probe_id(probes)]
  if (length(invalid) > 0) {
    warning(sprintf(
      "Skipping %d input(s) that do not look like probe identifiers.", length(invalid)
    ), call. = FALSE)
    probes <- setdiff(probes, invalid)
  }
  if (length(probes) < 20) {
    stop("EWAS Atlas enrichment requires at least 20 valid probe identifiers.", call. = FALSE)
  }
  if (length(probes) > 5000) {
    stop("EWAS Atlas enrichment accepts at most 5000 probe identifiers.", call. = FALSE)
  }

  types <- match.arg(types, choices = names(.enrichment_files), several.ok = TRUE)

  if (length(background) > 1 || !background[1] %in% c("450K", "850K")) {
    custom <- unique(as.character(background))
    custom <- custom[!is.na(custom) & nzchar(custom)]
    if (length(custom) < 20) {
      stop("A custom `background` must contain at least 20 probe identifiers.", call. = FALSE)
    }
    bk_string <- paste(custom, collapse = "\n")
    bk_label <- "other"
  } else {
    bk_string <- background[1]
    bk_label <- background[1]
  }

  job <- .ewas_job_id(probes, bk_string)

  exist <- .ewas_post("/checkJob", list(job = job, date = .now_string()))
  if (!isTRUE(exist$ok)) stop("Could not reach the EWAS Atlas toolkit service.", call. = FALSE)

  if (!identical(trimws(exist$value), "exist")) {
    check <- .ewas_post("/checkInput", list(
      job = job, input = paste(probes, collapse = "\n"), bk = bk_string, date = .now_string()
    ))
    if (!isTRUE(check$ok)) stop("Could not reach the EWAS Atlas toolkit service.", call. = FALSE)
    if (!identical(trimws(check$value), "success")) {
      stop(sprintf("EWAS Atlas rejected the input: %s", trimws(check$value)), call. = FALSE)
    }
  }

  results <- list()
  failed <- list()
  for (type in types) {
    if (isTRUE(progress)) message(sprintf("EWAScaller: running '%s' enrichment...", type))
    trigger <- .ewas_post("/analysis", list(job = job, type = type))
    if (!isTRUE(trigger$ok)) {
      failed[[type]] <- "could not start analysis"
      next
    }
    ok <- .poll_job_status(job, type, poll_interval = poll_interval, timeout = timeout, progress = progress)
    if (!ok) {
      failed[[type]] <- "timed out waiting for results"
      next
    }
    file <- .ewas_post("/getFile", list(file = .enrichment_files[[type]], job = job), as = "json")
    if (!isTRUE(file$ok) || is.null(file$value) || !identical(file$value$code, 0L)) {
      failed[[type]] <- "could not retrieve results"
      next
    }
    df <- file$value$data
    if (is.null(df) || (is.data.frame(df) && nrow(df) == 0)) {
      results[[type]] <- data.frame()
    } else {
      results[[type]] <- .coerce_numeric_cols(as.data.frame(df, stringsAsFactors = FALSE))
    }
  }

  new_ewas_enrichment(
    results = results, job = job, probes = probes,
    background = bk_label, failed = failed
  )
}

#' @keywords internal
#' @noRd
.poll_job_status <- function(job, type, poll_interval, timeout, progress) {
  deadline <- Sys.time() + timeout
  repeat {
    status <- .ewas_post("/jobStatus", list(job = job, type = type))
    if (isTRUE(status$ok) && identical(trimws(status$value), "finished")) return(TRUE)
    if (Sys.time() > deadline) return(FALSE)
    Sys.sleep(poll_interval)
  }
}

#' @keywords internal
#' @noRd
.now_string <- function() format(Sys.time(), "%a %b %d %Y %H:%M:%S GMT+0000 (Coordinated Universal Time)")

#' @keywords internal
#' @noRd
.coerce_numeric_cols <- function(df) {
  for (col in names(df)) {
    x <- df[[col]]
    if (is.character(x)) {
      num <- suppressWarnings(as.numeric(x))
      if (!anyNA(num[!is.na(x) & nzchar(x)]) || all(is.na(x) | !nzchar(x))) {
        # only coerce if every non-empty value parsed successfully
        if (!any(is.na(num) & !is.na(x) & nzchar(x))) df[[col]] <- num
      }
    }
  }
  df
}
