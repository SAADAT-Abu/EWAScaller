.ewas_rest_base <- "https://ngdc.cncb.ac.cn/ewas/rest"
.ewas_toolkit_base <- "https://ngdc.cncb.ac.cn/ewas/toolkit"
.ewas_user_agent <- "EWAScaller (https://github.com/SAADAT-Abu/EWAScaller)"

#' Perform a GET request against the EWAS Atlas REST API
#'
#' @param path API path appended to the REST base URL.
#' @param query Named list of query parameters.
#' @param times Number of attempts before giving up.
#' @return The parsed `data` element of the JSON response, or an
#'   `ewas_api_error` condition object on failure.
#' @keywords internal
#' @noRd
.ewas_get <- function(path, query = list(), times = 3L) {
  url <- paste0(.ewas_rest_base, path)
  resp <- tryCatch(
    httr::RETRY(
      "GET", url,
      query = query,
      httr::user_agent(.ewas_user_agent),
      times = times,
      pause_base = 1,
      pause_cap = 8,
      terminate_on = c(400, 401, 403, 404),
      quiet = TRUE
    ),
    error = function(e) e
  )
  .parse_ewas_response(resp, context = query)
}

#' Perform a POST request against the EWAS Atlas toolkit API
#'
#' @inheritParams .ewas_get
#' @param body Named list encoded as form data.
#' @param as Return type passed to [httr::content()].
#' @keywords internal
#' @noRd
.ewas_post <- function(path, body = list(), times = 3L, as = "text") {
  url <- paste0(.ewas_toolkit_base, path)
  resp <- tryCatch(
    httr::RETRY(
      "POST", url,
      body = body, encode = "form",
      httr::user_agent(.ewas_user_agent),
      times = times,
      pause_base = 1,
      pause_cap = 8,
      terminate_on = c(400, 401, 403, 404),
      quiet = TRUE
    ),
    error = function(e) e
  )
  if (inherits(resp, "error")) {
    return(structure(list(ok = FALSE, message = conditionMessage(resp)), class = "ewas_api_error"))
  }
  status <- httr::status_code(resp)
  if (status >= 400) {
    return(structure(list(ok = FALSE, message = paste("HTTP", status)), class = "ewas_api_error"))
  }
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  if (identical(as, "json")) {
    parsed <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = TRUE), error = function(e) NULL)
    return(list(ok = TRUE, value = parsed, raw = txt))
  }
  list(ok = TRUE, value = txt, raw = txt)
}

#' @keywords internal
#' @noRd
.parse_ewas_response <- function(resp, context = NULL) {
  if (inherits(resp, "error")) {
    return(structure(
      list(ok = FALSE, message = conditionMessage(resp), context = context),
      class = "ewas_api_error"
    ))
  }
  status <- httr::status_code(resp)
  if (status >= 400) {
    return(structure(
      list(ok = FALSE, message = paste("HTTP", status), context = context),
      class = "ewas_api_error"
    ))
  }
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  parsed <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(parsed)) {
    return(structure(
      list(ok = FALSE, message = "Could not parse API response", context = context),
      class = "ewas_api_error"
    ))
  }
  if (!identical(parsed$code, 0L) && !identical(parsed$code, 0)) {
    msg <- if (!is.null(parsed$msg)) parsed$msg else "EWAS Atlas API returned an error"
    return(structure(
      list(ok = FALSE, message = msg, context = context),
      class = "ewas_api_error"
    ))
  }
  list(ok = TRUE, value = parsed$data)
}

#' Validate CpG probe identifiers
#'
#' @param x Character vector of probe identifiers.
#' @return Logical vector, `TRUE` where `x` matches the Illumina probe ID
#'   pattern (e.g. `"cg05575921"`).
#' @keywords internal
#' @noRd
.is_valid_probe_id <- function(x) {
  grepl("^(cg|ch|rs)[0-9]+$", x, ignore.case = TRUE)
}

#' Dispatch a function over a vector of items with bounded concurrency
#'
#' Splits `items` into batches of size `workers` and evaluates each batch in
#' parallel, pausing `delay` seconds between batches. This bounds the number
#' of simultaneous requests sent to the remote server and paces successive
#' batches to avoid excessive call rates.
#'
#' @param items List or vector of inputs, one per call to `fun`.
#' @param fun Function taking a single element of `items` (plus `...`).
#' @param workers Maximum number of concurrent workers.
#' @param delay Seconds to pause between batches.
#' @param progress Print a simple progress indicator to the console.
#' @param ... Additional arguments passed to `fun`.
#' @return A list of results, one per element of `items`, in the original
#'   order.
#' @keywords internal
#' @noRd
.run_parallel <- function(items, fun, workers = 2L, delay = 1, progress = TRUE, ...) {
  n <- length(items)
  out <- vector("list", n)
  if (n == 0L) return(out)

  workers <- max(1L, as.integer(workers))
  done <- 0L

  report <- function(k) {
    if (isTRUE(progress)) {
      message(sprintf("\rEWAScaller: %d/%d queries complete", k, n), appendLF = (k == n))
    }
  }

  if (workers <= 1L || n == 1L) {
    for (i in seq_len(n)) {
      out[[i]] <- fun(items[[i]], ...)
      done <- done + 1L
      report(done)
      if (i < n && delay > 0) Sys.sleep(delay)
    }
    return(out)
  }

  cl_type <- if (.Platform$OS.type == "windows") "PSOCK" else "FORK"
  cl <- parallel::makeCluster(min(workers, n), type = cl_type)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  chunks <- split(seq_len(n), ceiling(seq_len(n) / workers))
  for (idx in chunks) {
    res <- parallel::parLapply(cl, items[idx], fun, ...)
    out[idx] <- res
    done <- done + length(idx)
    report(done)
    if (done < n && delay > 0) Sys.sleep(delay)
  }
  out
}

#' Compute the EWAS Atlas toolkit job identifier for a probe set
#'
#' The toolkit web service identifies an analysis job by the MD5 digest of
#' the sorted, newline-joined probe list concatenated with the background
#' identifier.
#'
#' @param probes Character vector of probe identifiers.
#' @param background Background identifier string (e.g. `"450K"`, `"850K"`,
#'   or a custom probe list string when background type is `"other"`).
#' @keywords internal
#' @noRd
.ewas_job_id <- function(probes, background) {
  key <- paste0(paste(sort(probes), collapse = "\n"), background)
  digest::digest(key, algo = "md5", serialize = FALSE)
}

#' @keywords internal
#' @noRd
.collapse_unique <- function(x, sep = "; ") {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = sep)
}
