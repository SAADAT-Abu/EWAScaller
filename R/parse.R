`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

#' Normalise a probe API payload into probe and association tables
#'
#' Handles both the single-object payload returned by `/probe` and the
#' data-frame-of-probes payload returned by `/gene` and `/pos`.
#'
#' @param data Parsed `data` element from an EWAS Atlas API response.
#' @return A list with `probes` and `associations` data frames (either may
#'   be `NULL` if no records were present).
#' @keywords internal
#' @noRd
.rows_from_probe_data <- function(data) {
  if (is.null(data) || length(data) == 0) return(list(probes = NULL, associations = NULL))

  if (is.data.frame(data)) {
    n <- nrow(data)
    if (n == 0) return(list(probes = NULL, associations = NULL))
    get_field <- function(name, i) if (name %in% names(data)) data[[name]][[i]] else NULL
    ids <- if ("probeId" %in% names(data)) data$probeId else rep(NA_character_, n)
  } else {
    n <- 1L
    get_field <- function(name, i) data[[name]]
    ids <- data$probeId %||% NA_character_
  }

  probe_rows <- vector("list", n)
  assoc_rows <- vector("list", n)

  for (i in seq_len(n)) {
    pid <- if (length(ids) >= i) as.character(ids[[i]]) else NA_character_

    chr <- get_field("chrHg19", i)
    chr <- if (length(chr) == 0) NA_character_ else as.character(chr)

    pos <- get_field("posHg19", i)
    pos <- if (length(pos) == 0) NA_real_ else as.numeric(pos)

    island <- get_field("cpgIsland", i)
    island <- if (length(island) == 0) NA_character_ else as.character(island)

    rt <- get_field("relatedTranscription", i)
    gene <- if (is.data.frame(rt) && nrow(rt) > 0 && "geneName" %in% names(rt)) {
      .collapse_unique(as.character(rt$geneName))
    } else {
      NA_character_
    }

    probe_rows[[i]] <- data.frame(
      probeId = pid, chrHg19 = chr, posHg19 = pos, cpgIsland = island, gene = gene,
      stringsAsFactors = FALSE
    )

    al <- get_field("associationList", i)
    if (is.data.frame(al) && nrow(al) > 0) {
      assoc_rows[[i]] <- data.frame(
        probeId = pid,
        studyId = if ("studyId" %in% names(al)) as.character(al$studyId) else NA_character_,
        trait = if ("trait" %in% names(al)) as.character(al$trait) else NA_character_,
        correlation = if ("correlation" %in% names(al)) as.character(al$correlation) else NA_character_,
        rank = if ("rank" %in% names(al)) suppressWarnings(as.integer(al$rank)) else NA_integer_,
        pmid = if ("pmid" %in% names(al)) suppressWarnings(as.integer(al$pmid)) else NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  }

  probes <- do.call(rbind, probe_rows)
  assoc_rows <- Filter(Negate(is.null), assoc_rows)
  associations <- if (length(assoc_rows) > 0) do.call(rbind, assoc_rows) else NULL

  list(probes = probes, associations = associations)
}

#' @keywords internal
#' @noRd
.empty_probes <- function() {
  data.frame(
    probeId = character(), chrHg19 = character(), posHg19 = numeric(),
    cpgIsland = character(), gene = character(), query = character(),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.empty_associations <- function() {
  data.frame(
    probeId = character(), studyId = character(), trait = character(),
    correlation = character(), rank = integer(), pmid = integer(),
    query = character(), stringsAsFactors = FALSE
  )
}
