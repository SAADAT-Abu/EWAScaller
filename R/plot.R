#' Word cloud of EWAS Atlas terms
#'
#' Draws a word cloud of trait terms sized by association frequency (for
#' [ewas_result] objects) or by enrichment significance (for
#' [ewas_enrichment] objects), or of eQTM genes sized by their strongest
#' association (for [ewas_eqtm-class] objects).
#'
#' @param x An [ewas_result], [ewas_enrichment], or [ewas_eqtm-class]
#'   object.
#' @param max_words Maximum number of terms to display.
#' @param colors A vector of colours to cycle through.
#' @param ... Passed to methods. All methods also accept `title`,
#'   `subtitle`, `base_size`, `theme` (see [autoplot.ewas_result()]),
#'   `max_size` (largest font size; by default 16, scaled up to 40 when
#'   there are only a few words), and `seed` (word layout seed, default 1).
#'
#' @return A `ggplot` object.
#' @export
plot_wordcloud <- function(x, max_words = 75L, colors = NULL, ...) UseMethod("plot_wordcloud")

#' @keywords internal
#' @noRd
.wordcloud_plot <- function(words, weights, colors, default_title = "EWAS Atlas trait terms",
                            title = NULL, subtitle = NULL, base_size = 11,
                            theme = ggplot2::theme_minimal, max_size = NULL, seed = 1, ...) {
  # With few words a fixed font size leaves the cloud tiny in an empty
  # panel, so scale the largest word up as the word count goes down
  # (16 for ~60+ words, up to 40 for a handful).
  if (is.null(max_size)) max_size <- max(16, min(40, 120 / sqrt(length(words))))
  colors <- .plot_colors(colors, c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51"))
  # Assign colours deterministically and pass the seed to the geom so the
  # layout is reproducible without touching the user's global RNG state.
  df <- data.frame(word = words, weight = weights, stringsAsFactors = FALSE)
  df$color <- rep_len(colors, nrow(df))
  p <- ggplot2::ggplot(df, ggplot2::aes(
    label = .data$word, size = .data$weight, color = .data$color
  )) +
    ggwordcloud::geom_text_wordcloud_area(rm_outside = TRUE, seed = seed) +
    ggplot2::scale_size_area(max_size = max_size) +
    ggplot2::scale_color_identity()
  .finish_plot(p, title, subtitle, default_title, base_size, theme)
}

#' @rdname plot_wordcloud
#' @export
plot_wordcloud.ewas_result <- function(x, max_words = 75L, colors = NULL, ...) {
  top <- top_traits(x, n = max_words)
  if (nrow(top) == 0) stop("No trait associations available to plot.", call. = FALSE)
  .wordcloud_plot(top$trait, top$n, colors, ...)
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
  label_col <- .enrichment_label_col(df)
  df <- df[order(df$p), , drop = FALSE]
  df <- utils::head(df, max_words)
  count_col <- intersect(c("DE", "DMC", "DMG"), names(df))[1]
  w <- if (weight == "count" && !is.na(count_col)) df[[count_col]] else -log10(pmax(df$p, 1e-300))
  .wordcloud_plot(df[[label_col]], w, colors,
                  default_title = sprintf("EWAS Atlas %s terms", gsub("_", " ", type)), ...)
}

#' @rdname plot_wordcloud
#' @export
plot_wordcloud.ewas_eqtm <- function(x, max_words = 75L, colors = NULL, ...) {
  pr <- x$pairs
  pr <- pr[!is.na(pr$gene_symbol) & nzchar(pr$gene_symbol), , drop = FALSE]
  if (nrow(pr) == 0) stop("No eQTM pairs available to plot.", call. = FALSE)
  best <- tapply(pr$p, pr$gene_symbol, min)
  best <- utils::head(sort(best), max_words)
  .wordcloud_plot(names(best), -log10(pmax(as.numeric(best), 1e-300)), colors,
                  default_title = "eQTM genes", ...)
}

#' @keywords internal
#' @noRd
.enrichment_label_col <- function(df) {
  label_col <- intersect(c("trait", "term", "pathway", "location", "gene"), names(df))[1]
  if (is.na(label_col)) stop("Could not find a label column to plot.", call. = FALSE)
  label_col
}

#' Plot summary statistics for an EWAS Atlas query result
#'
#' All EWAScaller plots return ordinary `ggplot` objects, so they can be
#' extended further with `+`. The arguments below cover the common
#' adjustments and are shared by every `autoplot()` method in the package.
#'
#' @param object An [ewas_result] object.
#' @param type One of `"traits"` (top trait frequency), `"chromosome"`
#'   (probes per chromosome), or `"direction"` (hyper-/hypomethylation
#'   balance).
#' @param top_n Number of categories to show for `type = "traits"`.
#' @param colors Colours to use. For single-colour plots the first colour
#'   is used. For `type = "direction"`, a vector named `pos` and `neg` (or
#'   two unnamed colours in that order).
#' @param title,subtitle Plot title and subtitle. `NULL` (default) uses the
#'   built-in title; `NA` or `""` removes it.
#' @param base_size Base font size passed to `theme`.
#' @param theme A ggplot2 theme function such as [ggplot2::theme_bw()]
#'   (called with `base_size`), or a complete theme object.
#' @param wrap Wrap long category labels at this many characters; `NULL`
#'   to disable wrapping.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
autoplot.ewas_result <- function(object, type = c("traits", "chromosome", "direction"), top_n = 15L,
                                 colors = NULL, title = NULL, subtitle = NULL, base_size = 11,
                                 theme = ggplot2::theme_minimal, wrap = 40, ...) {
  type <- match.arg(type)
  if (type == "traits") {
    top <- top_traits(object, n = top_n)
    if (nrow(top) == 0) stop("No trait associations available to plot.", call. = FALSE)
    top$label <- .wrap_labels(top$trait, wrap)
    p <- ggplot2::ggplot(top, ggplot2::aes(x = stats::reorder(.data$label, .data$n), y = .data$n)) +
      ggplot2::geom_col(fill = .plot_colors(colors, "#2a9d8f")[1]) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Number of associations")
    .finish_plot(p, title, subtitle, "Top EWAS Atlas trait terms", base_size, theme)
  } else if (type == "chromosome") {
    probes <- object$probes
    if (nrow(probes) == 0) stop("No probes available to plot.", call. = FALSE)
    tab <- as.data.frame(table(chr = probes$chrHg19), stringsAsFactors = FALSE)
    names(tab) <- c("chr", "n")
    chr_order <- c(as.character(1:22), "X", "Y", "M")
    tab$chr <- factor(tab$chr, levels = intersect(chr_order, unique(tab$chr)))
    p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$chr, y = .data$n)) +
      ggplot2::geom_col(fill = .plot_colors(colors, "#264653")[1]) +
      ggplot2::labs(x = "Chromosome (hg19)", y = "Number of probes")
    .finish_plot(p, title, subtitle, "Probe distribution by chromosome", base_size, theme)
  } else {
    assoc <- object$associations
    if (nrow(assoc) == 0) stop("No trait associations available to plot.", call. = FALSE)
    tab <- as.data.frame(table(correlation = assoc$correlation), stringsAsFactors = FALSE)
    names(tab) <- c("correlation", "n")
    fills <- .plot_colors(colors, c(pos = "#e76f51", neg = "#2a9d8f"))
    if (is.null(names(fills))) names(fills) <- c("pos", "neg")[seq_along(fills)]
    p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$correlation, y = .data$n, fill = .data$correlation)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = fills, na.value = "grey70") +
      ggplot2::labs(x = "Correlation direction", y = "Number of associations")
    .finish_plot(p, title, subtitle, "Direction of methylation-trait correlation", base_size, theme,
                 extra = ggplot2::theme(legend.position = "none"))
  }
}

#' Plot an EWAS Atlas enrichment result
#'
#' Draws the top enriched terms of one analysis, either as a dot plot
#' (mirroring the EWAS Atlas toolkit chart: position given by `-log10(p)`,
#' point size by the differentially-methylated CpG/gene count) or as a bar
#' chart of `-log10(p)`.
#'
#' @param object An [ewas_enrichment] object.
#' @param type Which analysis result to plot (default `"trait"`).
#' @param top_n Number of top terms to display.
#' @param style `"dot"` (default) or `"bar"`.
#' @param colors Colour of the points or bars (first element used).
#' @inheritParams autoplot.ewas_result
#' @return A `ggplot` object.
#' @export
autoplot.ewas_enrichment <- function(object, type = "trait", top_n = 20L, style = c("dot", "bar"),
                                     colors = NULL, title = NULL, subtitle = NULL, base_size = 11,
                                     theme = ggplot2::theme_minimal, wrap = 40, ...) {
  style <- match.arg(style)
  df <- object$results[[type]]
  if (is.null(df) || nrow(df) == 0) stop("No results available for that analysis type.", call. = FALSE)
  label_col <- .enrichment_label_col(df)
  df <- df[order(df$p), , drop = FALSE]
  df <- utils::head(df, top_n)
  df$label <- .wrap_labels(df[[label_col]], wrap)
  df$neglog10p <- -log10(pmax(df$p, 1e-300))
  col <- .plot_colors(colors, "#2a9d8f")[1]
  default_title <- sprintf("EWAS Atlas %s enrichment", gsub("_", " ", type))
  xlab <- expression(-log[10](italic(p)))

  if (style == "bar") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(.data$label, .data$neglog10p), y = .data$neglog10p)) +
      ggplot2::geom_col(fill = col) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = xlab)
    return(.finish_plot(p, title, subtitle, default_title, base_size, theme))
  }

  size_col <- intersect(c("DE", "DMC", "DMG"), names(df))[1]
  df$size_val <- if (!is.na(size_col)) df[[size_col]] else 1
  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$neglog10p, y = stats::reorder(.data$label, .data$neglog10p), size = .data$size_val
  )) +
    ggplot2::geom_point(color = col) +
    ggplot2::labs(x = xlab, y = NULL, size = if (!is.na(size_col)) size_col else NULL)
  .finish_plot(p, title, subtitle, default_title, base_size, theme)
}

#' Plot an EWAS Atlas eQTM result
#'
#' Visualises CpG-gene expression correlations from [ewas_eqtm()]:
#' \describe{
#'   \item{`"heatmap"`}{CpG probes against genes, coloured by correlation,
#'     one panel per tissue.}
#'   \item{`"volcano"`}{Correlation against `-log10(p)` for every pair,
#'     coloured by tissue, with the strongest pairs labelled.}
#'   \item{`"tissue"`}{Number of pairs per tissue, split into negative and
#'     positive correlations.}
#'   \item{`"scatter"`}{Per-sample methylation against expression for one
#'     pair (chosen with `probe`, `gene`, and `tissue`). This downloads the
#'     sample data with [eqtm_samples()], so it needs internet access.}
#' }
#'
#' @param object An [ewas_eqtm-class] object.
#' @param type One of `"heatmap"`, `"volcano"`, `"tissue"`, `"scatter"`.
#' @param tissue Optional character vector of tissues to keep. Required for
#'   `type = "scatter"` when the pair was found in more than one tissue.
#' @param top_n For `"heatmap"`, the number of strongest pairs (by p-value)
#'   to show.
#' @param label_top For `"volcano"`, how many of the strongest pairs to
#'   label (`0` for none).
#' @param probe,gene For `"scatter"`, the probe and gene of the pair to
#'   plot. `gene` may be a gene symbol or an Ensembl identifier.
#' @param colors Colours to use. `"heatmap"`: low, mid, and high colours of
#'   the correlation scale (default teal-white-orange). `"volcano"`: one
#'   colour per tissue (optionally named by tissue). `"tissue"`: a vector
#'   named `negative` and `positive` (or two unnamed colours in that
#'   order). `"scatter"`: point colour, then trend-line colour.
#' @inheritParams autoplot.ewas_result
#' @return A `ggplot` object.
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
#' autoplot(eq, type = "heatmap")
#' autoplot(eq, type = "volcano", theme = ggplot2::theme_bw, base_size = 14)
#' autoplot(eq, type = "scatter", probe = "cg05575921", gene = "AHRR", tissue = "brain")
#' }
#' @export
autoplot.ewas_eqtm <- function(object, type = c("heatmap", "volcano", "tissue", "scatter"),
                               tissue = NULL, top_n = 30L, label_top = 5L, probe = NULL, gene = NULL,
                               colors = NULL, title = NULL, subtitle = NULL, base_size = 11,
                               theme = ggplot2::theme_minimal, wrap = 40, ...) {
  type <- match.arg(type)
  pr <- object$pairs
  if (!is.null(tissue)) pr <- pr[pr$tissue %in% tissue, , drop = FALSE]
  if (nrow(pr) == 0) stop("No eQTM pairs available to plot.", call. = FALSE)
  pr$neglog10p <- -log10(pmax(pr$p, 1e-300))
  pr$direction <- ifelse(pr$cor < 0, "negative", "positive")

  switch(type,
    heatmap = .eqtm_heatmap(pr, top_n, colors, title, subtitle, base_size, theme, wrap),
    volcano = .eqtm_volcano(pr, label_top, colors, title, subtitle, base_size, theme),
    tissue = .eqtm_tissue(pr, colors, title, subtitle, base_size, theme),
    scatter = .eqtm_scatter(pr, probe, gene, colors, title, subtitle, base_size, theme)
  )
}

#' @keywords internal
#' @noRd
.eqtm_heatmap <- function(pr, top_n, colors, title, subtitle, base_size, theme, wrap) {
  pr <- utils::head(pr[order(pr$p), , drop = FALSE], top_n)
  pr$gene_label <- .wrap_labels(ifelse(is.na(pr$gene_symbol) | !nzchar(pr$gene_symbol),
                                       pr$gene_id, pr$gene_symbol), wrap)
  cols <- .plot_colors(colors, c("#2a9d8f", "white", "#e76f51"))
  if (length(cols) == 2) cols <- c(cols[1], "white", cols[2])
  lim <- max(abs(pr$cor), na.rm = TRUE)
  p <- ggplot2::ggplot(pr, ggplot2::aes(x = .data$probeId, y = .data$gene_label, fill = .data$cor)) +
    ggplot2::geom_tile(color = "grey85") +
    ggplot2::scale_fill_gradient2(low = cols[1], mid = cols[2], high = cols[3], midpoint = 0,
                                  limits = c(-lim, lim), name = "Pearson r") +
    ggplot2::labs(x = "CpG probe", y = "Gene")
  if (length(unique(pr$tissue)) > 1) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$tissue), scales = "free")
  }
  .finish_plot(p, title, subtitle, "eQTM: methylation-expression correlation", base_size, theme,
               extra = ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
}

#' @keywords internal
#' @noRd
.eqtm_volcano <- function(pr, label_top, colors, title, subtitle, base_size, theme) {
  tissues <- sort(unique(pr$tissue))
  p <- ggplot2::ggplot(pr, ggplot2::aes(x = .data$cor, y = .data$neglog10p, color = .data$tissue)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey70", linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::scale_color_manual(values = .category_colors(tissues, colors), name = "Tissue") +
    ggplot2::labs(x = "Pearson correlation (methylation vs expression)",
                  y = expression(-log[10](italic(p))))
  if (label_top > 0) {
    lab <- utils::head(pr[order(pr$p), , drop = FALSE], label_top)
    lab$label <- paste(lab$probeId, lab$gene_symbol, sep = " / ")
    # Put labels on the inner side of each point so they stay inside the
    # panel, and draw them in a neutral colour so pale tissue colours stay
    # readable.
    lab$hjust <- ifelse(lab$cor > 0, 1, 0)
    p <- p + ggplot2::geom_text(
      data = lab,
      mapping = ggplot2::aes(x = .data$cor, y = .data$neglog10p, label = .data$label, hjust = .data$hjust),
      inherit.aes = FALSE, color = "grey20", vjust = -0.8, size = 3, check_overlap = TRUE
    ) +
      ggplot2::coord_cartesian(clip = "off")
  }
  .finish_plot(p, title, subtitle, "eQTM volcano plot", base_size, theme)
}

#' @keywords internal
#' @noRd
.eqtm_tissue <- function(pr, colors, title, subtitle, base_size, theme) {
  tab <- as.data.frame(table(tissue = pr$tissue, direction = pr$direction), stringsAsFactors = FALSE)
  names(tab) <- c("tissue", "direction", "n")
  totals <- tapply(tab$n, tab$tissue, sum)
  tab$tissue <- factor(tab$tissue, levels = names(sort(totals)))
  fills <- .plot_colors(colors, c(negative = "#2a9d8f", positive = "#e76f51"))
  if (is.null(names(fills))) names(fills) <- c("negative", "positive")[seq_along(fills)]
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$tissue, y = .data$n, fill = .data$direction)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = fills, name = "Correlation") +
    ggplot2::labs(x = NULL, y = "Number of CpG-gene pairs")
  .finish_plot(p, title, subtitle, "eQTM pairs per tissue", base_size, theme)
}

#' @keywords internal
#' @noRd
.eqtm_scatter <- function(pr, probe, gene, colors, title, subtitle, base_size, theme) {
  if (is.null(probe) || is.null(gene)) {
    stop("`probe` and `gene` are required for type = \"scatter\".", call. = FALSE)
  }
  probe <- .normalize_probe_ids(probe, verbose = FALSE)
  hit <- pr[pr$probeId == probe & (pr$gene_symbol == gene | pr$gene_id == gene), , drop = FALSE]
  if (nrow(hit) == 0) {
    stop(sprintf("No eQTM pair for probe '%s' and gene '%s'%s.", probe, gene,
                 " (check `tissue` if you set it)"), call. = FALSE)
  }
  if (nrow(hit) > 1) {
    stop(sprintf("Probe '%s' and gene '%s' are paired in several tissues (%s); choose one with `tissue`.",
                 probe, gene, paste(unique(hit$tissue), collapse = ", ")), call. = FALSE)
  }
  samples <- eqtm_samples(hit$probeId, hit$gene_id, hit$tissue)
  if (nrow(samples) == 0) stop("No sample data returned for this pair.", call. = FALSE)

  cols <- .plot_colors(colors, c("#264653", "#e76f51"))
  if (length(cols) < 2) cols <- c(cols, "#e76f51")
  default_sub <- sprintf("r = %.3f, p = %s | %s (%s)", hit$cor, format(signif(hit$p, 3)),
                         hit$tissue, hit$source)
  p <- ggplot2::ggplot(samples, ggplot2::aes(x = .data$methylation, y = .data$expression)) +
    ggplot2::geom_point(color = cols[1], alpha = 0.6) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = cols[2], se = TRUE) +
    ggplot2::labs(x = sprintf("%s methylation (beta)", hit$probeId),
                  y = sprintf("%s expression (TPM)", hit$gene_symbol))
  p <- .finish_plot(p, title, NA, sprintf("%s vs %s", hit$probeId, hit$gene_symbol), base_size, theme)
  p + ggplot2::labs(subtitle = if (is.null(subtitle)) default_sub else .plot_title(subtitle, NULL))
}
