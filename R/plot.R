#' Plot marginal means from fe_mm results
#'
#' A ggplot2 wrapper for the standard conjoint MM plot: coloured point-ranges
#' with a matching vertical reference line at `mean_outcome`, percentages on
#' the x-axis, and markdown-rendered axis labels.
#'
#' @param data            Results returned by [fe_mm()] or [fe_mm_by_group()].
#' @param colors          Character vector of colors passed to 
#'                        `scale_color_manual`. Default `"black"`.
#' @param x_limits        Numeric vector of length 2 for x-axis limits 
#'                        or `NULL` for automatic limits (default `NULL`).
#' @param x_label         Label for the x axis (default `"Marginal Mean"`).
#' @param facet_var       Name of a column in `data` to facet on, or `NULL`
#'                        (default `NULL`).
#' @param facet_ncol      Number of columns in the facet layout (default `NULL`).
#' @param legend_position Position of the legend passed to `theme()` (default
#'                        `"none"`).
#' @param legend_nrow     Number of rows in the color legend guide (default `1`).
#'
#' @return A `ggplot` object.
#' @export
fe_mm_plot <- function(data,
                       colors          = "black",
                       x_limits        = NULL,
                       x_label         = "Marginal Mean",
                       facet_var       = NULL,
                       facet_ncol      = NULL,
                       legend_position = "bottom",
                       legend_nrow     = 1) {

  p <- ggplot2::ggplot(data, ggplot2::aes(x     = .data$estimate,
                                          y     = .data$term_nice,
                                          color = .data$group)) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = .data$mean_outcome, color = .data$group)
    ) +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = .data$boot_lower, xmax = .data$boot_upper),
      position = ggplot2::position_dodge2(width = 1, reverse = TRUE)
    ) +
    ggplot2::scale_y_discrete() +
    ggplot2::scale_color_manual(values = colors, na.translate = FALSE) +
    ggplot2::labs(x = x_label, y = NULL, color = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x     = ggtext::element_markdown(),
      axis.text.y     = ggtext::element_markdown(),
      axis.title.x    = ggtext::element_markdown(),
      strip.text.x    = ggplot2::element_text(),
      legend.position = legend_position,
      legend.margin   = ggplot2::margin(-10, 0, 10, 0)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = legend_nrow))

  if (!is.null(x_limits))
    p <- p + ggplot2::coord_cartesian(xlim = x_limits)

  if (!is.null(facet_var))
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)),
                                 ncol = facet_ncol)

  p
}
