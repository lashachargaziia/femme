#' Estimate Differences-in-Marginal-Means for one subgroup (with bootstrap CIs)
#'
#' Computes the difference-in-marginal-means (dMM) for every non-reference level
#' of every attribute in `feature_vars`. For each attribute a *univariate* model
#' `outcome ~ attribute` is fit - optionally with respondent fixed effects - so
#' that each level estimate is the difference in (within-respondent) marginal
#' means relative to the attribute's reference level. Stratified percentile
#' bootstrap confidence intervals and two-sided bootstrap p-values are returned.
#'
#' Unlike [fe_amce()], which fits all attributes jointly in a single regression,
#' `fe_dmm` estimates each attribute in isolation; the two coincide under a fully
#' randomised, orthogonal design but can differ otherwise.
#'
#' @param data           Conjoint data frame.
#' @param feature_vars   Character vector of conjoint attribute columns.
#' @param outcome_var    Name of the outcome column (default `"choice"`). Can
#'                       be both continuous and binary.
#' @param respondent_id  Name of the unique observation identifier
#'                       (default `"id"`). Set to `NULL` to estimate without
#'                       respondent fixed effects (plain OLS with HC-robust SEs
#'                       and unstratified bootstrap).
#' @param weight_var     Optional column name of post-stratification weights.
#'                       When supplied the model is weighted (default `NULL`).
#' @param group_var      Name of the factor variable used to split the sample
#'                       (e.g. `"gender"`).
#' @param group_val      Value of `group_var` for the subgroup of interest.
#' @param boot_reps      Bootstrap replications (default `999`).
#' @param conf_level     Confidence level for intervals (default `0.95`).
#' @param variable_lookup Optional data frame with two columns - `variable` and
#'                        `variable_nice` - for user-defined variable labels.
#'                        Use `NULL` to skip labelling.
#' @param cluster_var    Optional column name used for clustering standard
#'                       errors and stratifying the bootstrap. When `NULL`
#'                       (default), clustering uses `respondent_id`. Set this
#'                       when you want to cluster by respondent but estimate
#'                       without individual fixed effects
#'                       (`respondent_id = NULL`).
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
#' @param seed           Random seed for reproducibility (default `123`).
#' @param parallel       Parallelization for `boot::boot`: `"multicore"`
#'                       (Unix/Mac), `"snow"` (Windows), or `"no"`
#'                       (default `"no"`).
#' @param ncpus          Number of CPUs when `parallel != "no"`
#'                       (default `parallel::detectCores() - 1`).
#'
#' @return A tidy data frame with one bold-faced header row per attribute and
#'   `term_nice` as an ordered factor ready for plotting.
#' @export
#'
#' @examples
#' \dontrun{
#' fe_dmm(
#'   data            = icecream,
#'   feature_vars    = c("color", "taste", "cost"),
#'   variable_lookup = my_lookup
#' )
#' }
fe_dmm <- function(data,
                   feature_vars,
                   outcome_var     = "choice",
                   respondent_id   = "id",
                   weight_var      = NULL,
                   group_var       = NULL,
                   group_val       = NULL,
                   boot_reps       = 999,
                   conf_level      = 0.95,
                   variable_lookup = NULL,
                   cluster_var     = NULL,
                   reverse         = TRUE,
                   header          = TRUE,
                   seed            = 123,
                   parallel        = "no",
                   ncpus           = parallel::detectCores() - 1) suppressWarnings(suppressMessages({

  if (!is.null(group_var) && !is.null(group_val)) {
    var_expr      <- rlang::sym(group_var)
    filtered_data <- dplyr::filter(data, !!var_expr == group_val)
  } else {
    filtered_data <- data
  }
  if (!is.null(weight_var))
    filtered_data <- dplyr::filter(filtered_data, !is.na(.data[[weight_var]]))

  for (v in feature_vars) filtered_data[[v]] <- as.factor(filtered_data[[v]])

  strata_var <- cluster_var %||% respondent_id

  boot_fn <- function(bd, indices) {
    .dmm_point(bd[indices, ],
               outcome_var   = outcome_var,
               respondent_id = respondent_id,
               feature_vars  = feature_vars,
               weight_var    = weight_var,
               cluster_var   = cluster_var)$estimate
  }

  boot_res <- .run_boot(seed, parallel, function() boot::boot(
    data      = filtered_data,
    statistic = boot_fn,
    R         = boot_reps,
    strata    = if (!is.null(strata_var))
                  as.numeric(as.factor(filtered_data[[strata_var]]))
                else rep(1L, nrow(filtered_data)),
    parallel  = parallel,
    ncpus     = if (parallel != "no") ncpus else 1L
  ))

  dmm_final <- .dmm_point(filtered_data,
                          outcome_var   = outcome_var,
                          respondent_id = respondent_id,
                          feature_vars  = feature_vars,
                          weight_var    = weight_var,
                          cluster_var   = cluster_var)
  nobs_val <- attr(dmm_final, "nobs")

  boot_estimates <- colMeans(boot_res$t, na.rm = TRUE)
  ci_list <- lapply(seq_along(boot_estimates), function(i) {
    if (dmm_final$is_reference[i]) return(c(NA_real_, NA_real_))
    ci <- boot::boot.ci(boot_res, type = "perc", index = i, conf = conf_level)
    if (!is.null(ci$perc)) ci$perc[4:5] else c(NA_real_, NA_real_)
  })
  ci_bounds <- do.call(rbind, ci_list)

  p_values <- vapply(seq_along(boot_estimates), function(i) {
    if (dmm_final$is_reference[i]) return(NA_real_)
    draws <- boot_res$t[, i]
    min(2 * p_positive(draws), 2 * p_negative(draws))
  }, numeric(1))

  results_df <- data.frame(
    level         = dmm_final$level,
    boot_estimate = boot_estimates,
    boot_lower    = ci_bounds[, 1],
    boot_upper    = ci_bounds[, 2],
    p_value       = p_values
  )
  dmm_final <- dmm_final |> dplyr::right_join(results_df, by = "level")

  if (!is.null(variable_lookup)) {
    dmm_final <- dmm_final |>
      dplyr::left_join(variable_lookup,
                       by = dplyr::join_by("feature" == "variable")) |>
      dplyr::rename(term_nice = "level") |>
      dplyr::mutate(variable_nice = forcats::fct_inorder(.data$variable_nice))
  } else {
    dmm_final <- dmm_final |>
      dplyr::rename(term_nice = "level") |>
      dplyr::mutate(variable_nice = forcats::fct_inorder(.data$feature))
  }

  dmm_final$group <- as.factor(rep(group_val %||% "Pooled", nrow(dmm_final)))
  dmm_final$nobs  <- nobs_val

  term_order <- unlist(lapply(feature_vars,
                              function(v) levels(as.factor(filtered_data[[v]]))),
                       use.names = FALSE)
  .add_header_rows(dmm_final, level_order = term_order, reverse = reverse,
                   header = header)
}))


#' Estimate Differences-in-Marginal-Means for all group levels (with bootstrap CIs)
#'
#' Iterates [fe_dmm()] over every level of `group_var` and binds the results.
#'
#' @inheritParams fe_dmm
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
#' @return A tidy data frame.
#' @export
fe_dmm_by_group <- function(data,
                            group_var,
                            feature_vars,
                            outcome_var     = "choice",
                            respondent_id   = "id",
                            weight_var      = NULL,
                            boot_reps       = 999,
                            conf_level      = 0.95,
                            variable_lookup = NULL,
                            level_order     = NULL,
                            cluster_var     = NULL,
                            reverse         = TRUE,
                            header          = TRUE,
                            seed            = 123,
                            parallel        = "no",
                            ncpus           = parallel::detectCores() - 1) suppressWarnings(suppressMessages({

  if (!group_var %in% names(data))
    stop("`group_var` must be a column in `data`.")

  data[[group_var]] <- as.factor(data[[group_var]])
  group_names       <- levels(data[[group_var]])

  results <- lapply(group_names, function(gn) {
    fe_dmm(data,
           group_var       = group_var,
           group_val       = gn,
           feature_vars    = feature_vars,
           outcome_var     = outcome_var,
           respondent_id   = respondent_id,
           weight_var      = weight_var,
           boot_reps       = boot_reps,
           conf_level      = conf_level,
           variable_lookup = variable_lookup,
           cluster_var     = cluster_var,
           reverse         = reverse,
           header          = header,
           seed            = seed,
           parallel        = parallel,
           ncpus           = ncpus)
  })
  combined <- do.call(rbind, results)

  combined$variable_nice <- factor(combined$variable_nice,
                                   levels = levels(results[[1]]$variable_nice))
  combined$estimand      <- "dMM"

  term_order <- level_order %||%
    unlist(lapply(feature_vars, function(v) levels(as.factor(data[[v]]))),
           use.names = FALSE)

  .add_header_rows(combined, level_order = term_order, reverse = reverse,
                   header = header)
}))
