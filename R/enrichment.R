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
#' @param probes Character vector of probe identifiers (20-5000) from the
#'   450K, EPIC (850K), or EPIC v2 array. EPIC v2 names such as
#'   `"cg05575921_BC21"` are converted to their 450K/EPIC v1 names, and
#'   replicate probes are merged.
#' @param background The array the probes were selected from: `"450K"`
#'   (the default), `"850K"` (EPIC v1), `"EPICv2"`, or a character vector of
#'   probe identifiers to use as a custom background set. EWAS Atlas has no
#'   EPIC v2 background, so `"EPICv2"` uses the EPIC/850K background; this
#'   is an approximation because that background also contains EPIC v1
#'   probes that are not on EPIC v2.
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
  probes <- .prepare_toolkit_probes(probes)
  types <- match.arg(types, choices = names(.enrichment_files), several.ok = TRUE)
  bk <- .resolve_background(background)
  job <- .submit_toolkit_job(probes, bk$string)

  results <- list()
  failed <- list()
  for (type in types) {
    if (isTRUE(progress)) message(sprintf("EWAScaller: running '%s' enrichment...", type))
    out <- .run_toolkit_analysis(job, type, .enrichment_files[[type]],
                                 poll_interval = poll_interval, timeout = timeout, progress = progress)
    if (!isTRUE(out$ok)) {
      failed[[type]] <- out$message
    } else {
      results[[type]] <- out$data
    }
  }

  new_ewas_enrichment(
    results = results, job = job, probes = probes,
    background = bk$label, failed = failed
  )
}

#' Validate a probe set for submission to the EWAS Atlas toolkit
#' @keywords internal
#' @noRd
.prepare_toolkit_probes <- function(probes) {
  probes <- unique(.normalize_probe_ids(probes))
  probes <- probes[!is.na(probes) & nzchar(probes)]
  invalid <- probes[!.is_valid_probe_id(probes)]
  if (length(invalid) > 0) {
    warning(sprintf(
      "Skipping %d input(s) that do not look like probe identifiers.", length(invalid)
    ), call. = FALSE)
    probes <- setdiff(probes, invalid)
  }
  if (length(probes) < 20) {
    stop("The EWAS Atlas toolkit requires at least 20 valid probe identifiers.", call. = FALSE)
  }
  if (length(probes) > 5000) {
    stop("The EWAS Atlas toolkit accepts at most 5000 probe identifiers.", call. = FALSE)
  }
  probes
}

#' Turn a `background` argument into the toolkit's string and a label
#' @keywords internal
#' @noRd
.resolve_background <- function(background) {
  # An untouched default means "use the first preset" (450K), not a custom set.
  if (identical(background, c("450K", "850K"))) background <- "450K"
  if (length(background) == 1L && toupper(background) %in% c("EPICV2", "EPIC_V2", "935K")) {
    message(
      "EWAScaller: EWAS Atlas has no EPIC v2 background, so the EPIC/850K ",
      "background is used instead. It includes EPIC v1 probes that are not on ",
      "EPIC v2, so p-values are approximate."
    )
    return(list(string = "850K", label = "EPICv2 (approximated by 850K)"))
  }
  if (length(background) > 1 || !background[1] %in% c("450K", "850K")) {
    custom <- unique(.normalize_probe_ids(background, verbose = FALSE))
    custom <- custom[!is.na(custom) & nzchar(custom)]
    if (length(custom) < 20) {
      stop("A custom `background` must contain at least 20 probe identifiers.", call. = FALSE)
    }
    return(list(string = paste(custom, collapse = "\n"), label = "other"))
  }
  list(string = background[1], label = background[1])
}

#' Register a probe set with the toolkit (if not already known) and return the job id
#' @keywords internal
#' @noRd
.submit_toolkit_job <- function(probes, bk_string) {
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
  job
}

#' Trigger one toolkit analysis, wait for it, and download its result table
#' @return `list(ok = TRUE, data = <data.frame>)` or `list(ok = FALSE, message = <chr>)`.
#' @keywords internal
#' @noRd
.run_toolkit_analysis <- function(job, type, file, poll_interval, timeout, progress) {
  trigger <- .ewas_post("/analysis", list(job = job, type = type))
  if (!isTRUE(trigger$ok)) return(list(ok = FALSE, message = "could not start analysis"))
  status <- .poll_job_status(job, type, poll_interval = poll_interval, timeout = timeout, progress = progress)
  if (!isTRUE(status)) return(list(ok = FALSE, message = status))
  res <- .ewas_post("/getFile", list(file = file, job = job), as = "json")
  if (!isTRUE(res$ok) || is.null(res$value) || !identical(res$value$code, 0L)) {
    return(list(ok = FALSE, message = "could not retrieve results"))
  }
  df <- res$value$data
  if (is.null(df) || (is.data.frame(df) && nrow(df) == 0)) {
    return(list(ok = TRUE, data = data.frame()))
  }
  list(ok = TRUE, data = .coerce_numeric_cols(as.data.frame(df, stringsAsFactors = FALSE)))
}

#' @keywords internal
#' @noRd
.poll_job_status <- function(job, type, poll_interval, timeout, progress) {
  deadline <- Sys.time() + timeout
  repeat {
    status <- .ewas_post("/jobStatus", list(job = job, type = type))
    if (isTRUE(status$ok)) {
      value <- trimws(status$value)
      if (identical(value, "finished")) return(TRUE)
      # Stop early if the server reports the job failed rather than
      # waiting out the full timeout.
      if (grepl("error|fail", value, ignore.case = TRUE)) {
        return(sprintf("analysis failed on the server: %s", value))
      }
    }
    if (Sys.time() > deadline) return("timed out waiting for results")
    Sys.sleep(poll_interval)
  }
}

#' @keywords internal
#' @noRd
.now_string <- function() {
  # Mimic a browser's Date.toString() in UTC, with English day/month names
  # regardless of the user's locale.
  old <- Sys.getlocale("LC_TIME")
  on.exit(invisible(suppressWarnings(Sys.setlocale("LC_TIME", old))), add = TRUE)
  invisible(suppressWarnings(Sys.setlocale("LC_TIME", "C")))
  format(Sys.time(), "%a %b %d %Y %H:%M:%S GMT+0000 (Coordinated Universal Time)", tz = "UTC")
}

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
