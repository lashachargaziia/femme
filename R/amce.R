#' Estimate AMCEs for one subgroup
#'
#' Fits a linear model with respondent fixed effects and returns a tidy
#' coefficient table annotated with user-defined labels.
#'
#' @param data           Conjoint data frame.
#' @param group_var      Name of the factor variable used to split the sample
#'                       (e.g. `"gender"`).
#' @param group_val      Value of `group_var` for the subgroup of interest.
#' @param feature_vars   Character vector of conjoint attribute columns.
#' @param outcome_var    Name of the outcome column (default `"choice"`). Can
#'                       be both continuous and binary.
#' @param respondent_id  Name of the unique observation identifier
#'                       (default `"id"`). Set to `NULL` to estimate without
#'                       respondent fixed effects (plain OLS with HC-robust SEs).
#' @param cluster_var    Optional column name used for clustering standard
#'                       errors. When `NULL` (default), clustering uses
#'                       `respondent_id`. Set this when you want to cluster by
#'                       respondent but estimate without individual fixed effects
#'                       (`respondent_id = NULL`).
#' @param variable_lookup Optional data frame with two columns - `variable` and
#'                        `variable_nice` - for user-defined variable labels.
#'                        Use `NULL` to skip labelling.
#' @param reverse        Reverse the `term_nice` factor order so the first
#'                       attribute plots at the top of a horizontal
#'                       `y = term_nice` plot (default `TRUE`). `term_nice` is
#'                       returned as an ordered factor - attributes in
#'                       `feature_vars` order, levels in their data factor order,
#'                       with one bold-faced header row per attribute - so it can
#'                       be plotted directly without any manual re-levelling.
#' @param header         Insert one bold-faced header row per attribute
#'                       (default `TRUE`). Set `FALSE` to omit the header rows;
#'                       `term_nice` is still returned as the same ordered factor,
#'                       just without the header levels.
#'
#' @return A tidy data frame with one bold-faced header row per attribute and
#'   `term_nice` as an ordered factor ready for plotting.
#' @export
#'
#' @examples
#' \dontrun{
#' fe_amce(
#'   data            = icecream,
#'   group_var       = "gender",
#'   group_val       = "Female",
#'   feature_vars    = c("color", "taste", "cost"),
#'   variable_lookup = my_lookup
#' )
#' }
fe_amce <- function(data,
                    group_var       = NULL,
                    group_val       = NULL,
                    feature_vars,
                    outcome_var     = "choice",
                    respondent_id   = "id",
                    cluster_var     = NULL,
                    variable_lookup = NULL,
                    reverse         = TRUE,
                    header          = TRUE) suppressWarnings(suppressMessages({

  if (!is.null(group_var) && !is.null(group_val)) {
    var_expr      <- rlang::sym(group_var)
    filtered_data <- dplyr::filter(data, !!var_expr == group_val)
  } else {
    filtered_data <- data
  }

  rhs <- paste(paste0("`", feature_vars, "`"), collapse = " + ")
  if (!is.null(respondent_id)) {
    fml <- stats::as.formula(
      paste0("`", outcome_var, "` ~ ", rhs, " | `", respondent_id, "`")
    )
  } else {
    fml <- stats::as.formula(paste0("`", outcome_var, "` ~ ", rhs))
  }

  cl <- cluster_var %||% respondent_id
  cluster_fml <- if (!is.null(cl)) stats::as.formula(paste0("~", cl)) else NULL

  model <- fixest::feols(fml, data = filtered_data, cluster = cluster_fml)

  tidy_model <- model |>
    broom.helpers::tidy_and_attach(include = "nobs") |>
    broom.helpers::tidy_add_reference_rows() |>
    broom.helpers::tidy_add_estimate_to_reference_rows() |>
    broom.helpers::tidy_add_variable_labels() |>
    dplyr::filter(!.data$variable %in% c("(Intercept)", respondent_id,
                                         paste0("factor(", respondent_id, ")"))) |>
    dplyr::rename(se = "std.error") |>
    dplyr::mutate(term_nice = stringr::str_remove(.data$term, .data$variable))

  if (!is.null(variable_lookup)) {
    tidy_model <- tidy_model |>
      dplyr::left_join(variable_lookup, by = dplyr::join_by("variable")) |>
      dplyr::mutate(variable_nice = forcats::fct_inorder(.data$variable_nice))
  } else {
    tidy_model <- tidy_model |>
      dplyr::mutate(variable_nice = forcats::fct_inorder(.data$variable))
  }

  tidy_model$nobs  <- stats::nobs(model)
  tidy_model$group <- as.factor(rep(group_val %||% "Pooled", nrow(tidy_model)))

  # Keep only the meaningful columns (dropping broom.helpers bookkeeping such as
  # `original_term`, `instrumental`, `var_label`, `var_class`, `contrasts`) and
  # name them consistently with fe_mm()/fe_dmm(): `feature`, `is_reference`,
  # `p_value`.
  tidy_model <- dplyr::select(tidy_model, dplyr::any_of(c(
    "variable_nice", "term_nice", feature = "variable",
    "estimate", is_reference = "reference_row",
    "se", "statistic", p_value = "p.value", "conf.low", "conf.high",
    "group", "nobs"
  )))

  term_order <- unlist(lapply(feature_vars,
                              function(v) levels(as.factor(filtered_data[[v]]))),
                       use.names = FALSE)
  .add_header_rows(tidy_model, level_order = term_order, reverse = reverse,
                   header = header)
}))

#' Estimate AMCEs for all levels of a grouping variable
#'
#' A wrapper around [fe_amce()] that iterates over every factor level of
#' `group_var` and binds the results.
#'
#' @inheritParams fe_amce
#' @param level_order Optional character vector giving the desired within-attribute
#'   order of the attribute levels. If `NULL` the data's factor-level order is
#'   used.
#' @param reverse Reverse the `term_nice` factor order so the first attribute
#'   plots at the top of a horizontal `y = term_nice` plot (default `TRUE`).
#'   `term_nice` is returned as an ordered factor - attributes in `feature_vars`
#'   order, levels in their data factor order, with one bold-faced header row per
#'   attribute - so it can be plotted directly without any manual re-levelling.
#' @param header Insert one bold-faced header row per attribute (default `TRUE`).
#'   Set `FALSE` to omit the header rows; `term_nice` is still returned as the
#'   same ordered factor, just without the header levels.
#'
#' @return A stacked data frame (one block per group level) with the same
#'   columns as [fe_amce()], plus header rows that bold-face the attribute
#'   name for convenient plotting.
#' @export
#'
#' @examples
#' \dontrun{
#' fe_amce_by_group(
#'   data            = icecream,
#'   group_var       = "gender",
#'   feature_vars    = c("color", "taste", "cost"),
#'   variable_lookup = my_lookup
#' )
#' }
fe_amce_by_group <- function(data,
                             group_var,
                             feature_vars,
                             outcome_var     = "choice",
                             respondent_id   = "id",
                             cluster_var     = NULL,
                             variable_lookup = NULL,
                             level_order     = NULL,
                             reverse         = TRUE,
                             header          = TRUE) suppressWarnings(suppressMessages({

  if (!group_var %in% names(data))
    stop("`group_var` must be a column in `data`.")

  data[[group_var]] <- as.factor(data[[group_var]])
  group_names       <- levels(data[[group_var]])

  results <- lapply(group_names, function(gn) {
    fe_amce(data,
            group_var       = group_var,
            group_val       = gn,
            feature_vars    = feature_vars,
            outcome_var     = outcome_var,
            respondent_id   = respondent_id,
            cluster_var     = cluster_var,
            variable_lookup = variable_lookup,
            reverse         = reverse,
            header          = header)
  })
  combined <- do.call(rbind, results)

  combined$variable_nice <- factor(combined$variable_nice,
                                   levels = levels(results[[1]]$variable_nice))
  combined$estimand      <- "AMCE"

  term_order <- level_order %||%
    unlist(lapply(feature_vars, function(v) levels(as.factor(data[[v]]))),
           use.names = FALSE)

  .add_header_rows(combined, level_order = term_order, reverse = reverse,
                   header = header)
}))
