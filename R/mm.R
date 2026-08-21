#' Estimate Marginal Means for one subgroup (with bootstrap CIs)
#'
#' Computes marginal means adjusted for individual fixed-effects for every level
#' of every attribute in `feature_vars`, and produces stratified percentile 
#' bootstrap confidence intervals.
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
#'                       When supplied the model is weighted and the mean fixed
#'                       effect is computed as a weighted mean (default `NULL`).
#' @param group_var      Name of the factor variable used to split the sample
#'                       (e.g. `"gender"`).
#' @param group_val      Value of `group_var` for the subgroup of interest.
#' @param offset         Adjust for within-subgroup mean (default `"yes"`).
#' @param boot_reps      Bootstrap replications (default `999`).
#' @param conf_level     Confidence level for intervals (default `0.95`).
#' @param variable_lookup Optional data frame with two columns - `variable` and
#'                        `variable_nice` - for user-defined variable labels.
#'                        Use `NULL` to skip labelling.
#' @param seed           Random seed for reproducibility (default `123`).
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
#' @param parallel       Parallelization for `boot::boot`: `"multicore"`
#'                       (Unix/Mac), `"snow"` (Windows), or `"no"`
#'                       (default `"no"`).
#' @param ncpus          Number of CPUs when `parallel != "no"`
#'                       (default `parallel::detectCores() - 1`).
#'
#' @return A tidy data frame. Alongside the per-level `estimate`, `boot_lower`
#'   and `boot_upper`, it carries the grand mean of the outcome in
#'   `mean_outcome` (the usual `geom_vline` reference) together with an analytic
#'   Wald interval in `mean_outcome_lower` / `mean_outcome_upper` (from an
#'   intercept-only model, cluster-robust by respondent and weight-aware), so the
#'   reference line can be drawn with an uncertainty band.
#' @export
fe_mm <- function(data,
                  feature_vars,
                  outcome_var     = "choice",
                  respondent_id   = "id",
                  weight_var      = NULL,
                  group_var       = NULL,
                  group_val       = NULL,
                  offset          = "yes",
                  boot_reps       = 999,
                  conf_level      = 0.95,
                  variable_lookup = NULL,
                  cluster_var     = NULL,
                  reverse         = TRUE,
                  header          = TRUE,
                  seed            = 123,
                  parallel        = "no",
                  ncpus           = parallel::detectCores() - 1) suppressWarnings(suppressMessages({

  offset <- match.arg(offset, c("yes", "no"))

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
    .mm_point(bd[indices, ],
              outcome_var   = outcome_var,
              respondent_id = respondent_id,
              feature_vars  = feature_vars,
              offset        = offset,
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

  boot_estimates <- colMeans(boot_res$t, na.rm = TRUE)
  ci_list        <- lapply(seq_along(boot_estimates), function(i) {
    ci <- boot::boot.ci(boot_res, type = "perc", index = i, conf = conf_level)
    if (!is.null(ci$perc)) ci$perc[4:5] else c(NA_real_, NA_real_)
  })
  ci_bounds <- do.call(rbind, ci_list)

  mm_final <- .mm_point(filtered_data,
                        outcome_var   = outcome_var,
                        respondent_id = respondent_id,
                        feature_vars  = feature_vars,
                        offset        = offset,
                        weight_var    = weight_var,
                        cluster_var   = cluster_var)

  nobs_val <- attr(mm_final, "nobs")
  obs_used <- attr(mm_final, "obs_used")
  mo_data  <- if (!is.null(obs_used))
    filtered_data[obs_used, , drop = FALSE] else filtered_data

  results_df <- data.frame(
    level         = mm_final$level,
    boot_estimate = boot_estimates,
    boot_lower    = ci_bounds[, 1],
    boot_upper    = ci_bounds[, 2]
  )
  mm_final <- mm_final |> dplyr::right_join(results_df, by = "level")

  if (!is.null(variable_lookup)) {
    mm_final <- mm_final |>
      dplyr::left_join(variable_lookup,
                       by = dplyr::join_by("feature" == "variable")) |>
      dplyr::rename(term_nice = "level") |>
      dplyr::mutate(variable_nice = forcats::fct_inorder(.data$variable_nice))
  } else {
    mm_final <- mm_final |>
      dplyr::rename(term_nice = "level") |>
      dplyr::mutate(variable_nice = forcats::fct_inorder(.data$feature))
  }

  mm_final$group        <- as.factor(rep(group_val %||% "Pooled", nrow(mm_final)))
  mm_final$mean_outcome <- if (!is.null(weight_var)) {
    stats::weighted.mean(mo_data[[outcome_var]], mo_data[[weight_var]], na.rm = TRUE)
  } else {
    mean(mo_data[[outcome_var]], na.rm = TRUE)
  }

  cl_mo  <- cluster_var %||% respondent_id
  mo_mod <- fixest::feols(
    stats::as.formula(paste0("`", outcome_var, "` ~ 1")),
    data    = mo_data,
    weights = if (!is.null(weight_var)) stats::as.formula(paste0("~", weight_var)) else NULL,
    cluster = if (!is.null(cl_mo)) stats::as.formula(paste0("~", cl_mo)) else NULL
  )
  mo_se <- sqrt(stats::vcov(mo_mod)[1, 1])
  mo_z  <- stats::qnorm(1 - (1 - conf_level) / 2)
  mm_final$mean_outcome_lower <- mm_final$mean_outcome - mo_z * mo_se
  mm_final$mean_outcome_upper <- mm_final$mean_outcome + mo_z * mo_se
  mm_final$offset <- offset
  mm_final$nobs   <- nobs_val

  term_order <- unlist(lapply(feature_vars,
                              function(v) levels(as.factor(filtered_data[[v]]))),
                       use.names = FALSE)
  .add_header_rows(mm_final, level_order = term_order, reverse = reverse,
                   header = header)
}))


#' Estimate Marginal Means for all group levels (with bootstrap CIs)
#'
#' Iterates [fe_mm()] over every level of `group_var` and binds the results.
#'
#' @inheritParams fe_mm
#' @param level_order Optional character vector giving the desired within-attribute
#'   order of the attribute levels. If `NULL` the data's factor-level order is
#'   used.
#'
#' @return A tidy data frame.
#' @export
fe_mm_by_group <- function(data,
                           group_var,
                           feature_vars,
                           outcome_var     = "choice",
                           respondent_id   = "id",
                           weight_var      = NULL,
                           offset          = "yes",
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

  offset <- match.arg(offset, c("yes", "no"))

  data[[group_var]] <- as.factor(data[[group_var]])
  group_names       <- levels(data[[group_var]])

  results <- lapply(group_names, function(gn) {
    fe_mm(data,
          group_var       = group_var,
          group_val       = gn,
          feature_vars    = feature_vars,
          outcome_var     = outcome_var,
          respondent_id   = respondent_id,
          weight_var      = weight_var,
          offset          = offset,
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

  term_order <- level_order %||%
    unlist(lapply(feature_vars, function(v) levels(as.factor(data[[v]]))),
           use.names = FALSE)

  .add_header_rows(combined, level_order = term_order, reverse = reverse,
                   header = header) |>
    dplyr::mutate(
      boot_se        = (.data$boot_upper - .data$boot_lower) / (2 * stats::qnorm(1 - (1 - conf_level) / 2)),
      diff_from_mean = .data$estimate - .data$mean_outcome,
      z_from_mean    = .data$diff_from_mean / .data$boot_se,
      p_from_mean    = 2 * stats::pnorm(-abs(.data$z_from_mean))
    )
}))
