# Shared helpers that give every EWAScaller plot the same customisation
# arguments: colors, title, subtitle, base_size, theme, and wrap.

.default_palette <- c(
  "#264653", "#2a9d8f", "#e9c46a", "#f4a261",
  "#e76f51", "#8ab17d", "#6d597a", "#b56576"
)

#' Resolve a user-supplied title: NULL means "use the default", while NA
#' or "" means "no title".
#' @keywords internal
#' @noRd
.plot_title <- function(title, default) {
  if (is.null(title)) return(default)
  if (is.atomic(title) && length(title) == 1L && (is.na(title) || identical(title, ""))) return(NULL)
  title
}

#' Build the base theme from a theme function or a theme object
#' @keywords internal
#' @noRd
.plot_theme <- function(theme, base_size) {
  if (is.function(theme)) return(theme(base_size = base_size))
  if (inherits(theme, "theme")) return(theme)
  stop("`theme` must be a ggplot2 theme function (e.g. ggplot2::theme_bw) or a theme object.",
       call. = FALSE)
}

#' Wrap long axis labels onto several lines
#' @keywords internal
#' @noRd
.wrap_labels <- function(x, width) {
  x <- as.character(x)
  if (is.null(width) || !is.finite(width)) return(x)
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    paste(strwrap(s, width = width), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

#' Use the user's colours if supplied, otherwise the defaults
#' @keywords internal
#' @noRd
.plot_colors <- function(colors, default) {
  if (is.null(colors) || length(colors) == 0) default else colors
}

#' Map a vector of category values to colours, honouring named colour vectors
#' @keywords internal
#' @noRd
.category_colors <- function(levels, colors) {
  colors <- .plot_colors(colors, .default_palette)
  if (!is.null(names(colors)) && all(levels %in% names(colors))) return(colors[levels])
  stats::setNames(rep_len(unname(colors), length(levels)), levels)
}

#' Apply title, subtitle, and theme to a finished plot
#'
#' `extra` holds plot-specific theme tweaks (e.g. hiding a legend); they are
#' added after the base theme so the base theme cannot undo them.
#' @keywords internal
#' @noRd
.finish_plot <- function(p, title, subtitle, default_title, base_size, theme, extra = NULL) {
  p <- p +
    .plot_theme(theme, base_size) +
    ggplot2::labs(
      title = .plot_title(title, default_title),
      subtitle = .plot_title(subtitle, NULL)
    )
  if (!is.null(extra)) p <- p + extra
  p
}
