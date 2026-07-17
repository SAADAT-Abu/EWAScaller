#' Word cloud of EWAS Atlas trait terms
#'
#' Draws a word cloud sized by association frequency (for [ewas_result]
#' objects) or by enrichment significance (for [ewas_enrichment] objects).
#'
#' @param x An [ewas_result] or [ewas_enrichment] object.
#' @param max_words Maximum number of trait terms to display.
#' @param colors A vector of colors to cycle through.
#' @param ... Passed to methods.
#'
#' @return A `ggplot` object.
#' @export
plot_wordcloud <- function(x, max_words = 75L, colors = NULL, ...) UseMethod("plot_wordcloud")

#' @keywords internal
#' @noRd
.wordcloud_plot <- function(words, weights, colors, seed = 1) {
  if (is.null(colors)) {
    colors <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")
  }
  set.seed(seed)
  df <- data.frame(word = words, weight = weights, stringsAsFactors = FALSE)
  df$color <- sample(colors, nrow(df), replace = TRUE)
  ggplot2::ggplot(df, ggplot2::aes(
    label = .data$word, size = .data$weight, color = .data$color
  )) +
    ggwordcloud::geom_text_wordcloud_area(rm_outside = TRUE) +
    ggplot2::scale_size_area(max_size = 16) +
    ggplot2::scale_color_identity() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "EWAS Atlas trait terms")
}

#' @rdname plot_wordcloud
#' @export
plot_wordcloud.ewas_result <- function(x, max_words = 75L, colors = NULL, ...) {
  top <- top_traits(x, n = max_words)
  if (nrow(top) == 0) stop("No trait associations available to plot.", call. = FALSE)
  .wordcloud_plot(top$trait, top$n, colors)
}

#' @param type Which enrichment analysis to plot the word cloud for
#'   (default `"trait"`).
#' @param weight Weighting statistic: `"neglog10p"` (default, significance)
#'   or `"count"` (differentially methylated CpG/gene count for the term).
#' @rdname plot_wordcloud
#' @export
plot_wordcloud.ewas_enrichment <- function(x, max_words = 75L, colors = NULL, type = "trait",
                                            weight = c("neglog10p", "count"), ...) {
  weight <- match.arg(weight)
  df <- x$results[[type]]
  if (is.null(df) || nrow(df) == 0) stop("No results available for that analysis type.", call. = FALSE)
  label_col <- intersect(c("trait", "term", "pathway", "location", "gene"), names(df))[1]
  if (is.na(label_col)) stop("Could not find a label column to plot.", call. = FALSE)
  df <- df[order(df$p), , drop = FALSE]
  df <- utils::head(df, max_words)
  count_col <- intersect(c("DE", "DMC", "DMG"), names(df))[1]
  w <- if (weight == "count" && !is.na(count_col)) df[[count_col]] else -log10(pmax(df$p, 1e-300))
  .wordcloud_plot(df[[label_col]], w, colors)
}

#' Plot summary statistics for an EWAS Atlas query result
#'
#' @param object An [ewas_result] object.
#' @param type One of `"traits"` (top trait frequency), `"chromosome"`
#'   (probes per chromosome), or `"direction"` (hyper-/hypomethylation
#'   balance).
#' @param top_n Number of categories to show for `type = "traits"` or
#'   `type = "chromosome"`.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
autoplot.ewas_result <- function(object, type = c("traits", "chromosome", "direction"), top_n = 15L, ...) {
  type <- match.arg(type)
  if (type == "traits") {
    top <- top_traits(object, n = top_n)
    if (nrow(top) == 0) stop("No trait associations available to plot.", call. = FALSE)
    ggplot2::ggplot(top, ggplot2::aes(x = stats::reorder(.data$trait, .data$n), y = .data$n)) +
      ggplot2::geom_col(fill = "#2a9d8f") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Number of associations", title = "Top EWAS Atlas trait terms") +
      ggplot2::theme_minimal()
  } else if (type == "chromosome") {
    probes <- object$probes
    if (nrow(probes) == 0) stop("No probes available to plot.", call. = FALSE)
    tab <- as.data.frame(table(chr = probes$chrHg19), stringsAsFactors = FALSE)
    names(tab) <- c("chr", "n")
    chr_order <- c(as.character(1:22), "X", "Y", "M")
    tab$chr <- factor(tab$chr, levels = intersect(chr_order, unique(tab$chr)))
    ggplot2::ggplot(tab, ggplot2::aes(x = .data$chr, y = .data$n)) +
      ggplot2::geom_col(fill = "#264653") +
      ggplot2::labs(x = "Chromosome (hg19)", y = "Number of probes", title = "Probe distribution by chromosome") +
      ggplot2::theme_minimal()
  } else {
    assoc <- object$associations
    if (nrow(assoc) == 0) stop("No trait associations available to plot.", call. = FALSE)
    tab <- as.data.frame(table(correlation = assoc$correlation), stringsAsFactors = FALSE)
    names(tab) <- c("correlation", "n")
    ggplot2::ggplot(tab, ggplot2::aes(x = .data$correlation, y = .data$n, fill = .data$correlation)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = c(pos = "#e76f51", neg = "#2a9d8f", NA_ = "grey70"), na.value = "grey70") +
      ggplot2::labs(x = "Correlation direction", y = "Number of associations",
                    title = "Direction of methylation-trait correlation") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
  }
}

#' Plot an EWAS Atlas enrichment result
#'
#' Draws a dot plot of enrichment significance, mirroring the EWAS Atlas
#' toolkit's trait-enrichment chart: one row per term, position given by
#' `-log10(p)` and point size by the differentially-methylated CpG count.
#'
#' @param object An [ewas_enrichment] object.
#' @param type Which analysis result to plot (default `"trait"`).
#' @param top_n Number of top terms to display.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
autoplot.ewas_enrichment <- function(object, type = "trait", top_n = 20L, ...) {
  df <- object$results[[type]]
  if (is.null(df) || nrow(df) == 0) stop("No results available for that analysis type.", call. = FALSE)
  label_col <- intersect(c("trait", "term", "pathway", "location", "gene"), names(df))[1]
  if (is.na(label_col)) stop("Could not find a label column to plot.", call. = FALSE)
  df <- df[order(df$p), , drop = FALSE]
  df <- utils::head(df, top_n)
  df$label <- df[[label_col]]
  df$neglog10p <- -log10(pmax(df$p, 1e-300))
  size_col <- intersect(c("DE", "DMC", "DMG"), names(df))[1]
  size_col <- if (length(size_col) == 0) NA else size_col
  df$size_val <- if (!is.na(size_col)) df[[size_col]] else 1

  ggplot2::ggplot(df, ggplot2::aes(
    x = .data$neglog10p, y = stats::reorder(.data$label, .data$neglog10p), size = .data$size_val
  )) +
    ggplot2::geom_point(color = "#2a9d8f") +
    ggplot2::labs(
      x = expression(-log[10](italic(p))), y = NULL, size = if (!is.na(size_col)) size_col else NULL,
      title = sprintf("EWAS Atlas %s enrichment", type)
    ) +
    ggplot2::theme_minimal()
}
