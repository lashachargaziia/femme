#' Within-subgroup pairwise MM differences with bootstrap inference
#'
#' For each level of `group_var`, computes all pairwise differences between
#' attribute levels within that subgroup, together with two-sided p-values.
#'
#' @param data           Conjoint data frame.
#' @param group_var      Name of the factor variable used to split the sample
#'                       (e.g. `"gender"`).
#' @param feature_vars   Character vector of conjoint attribute columns.
#' @param outcome_var    Name of the outcome column (default `"choice"`). Can
#'                       be both continuous and binary.
#' @param respondent_id  Name of the unique observation identifier
#'                       (default `"id"`). Set to `NULL` to estimate without
#'                       respondent fixed effects (plain OLS with HC-robust SEs
#'                       and unstratified bootstrap).
#' @param boot_reps      Bootstrap replications (default `999`).
#' @param conf_level     Confidence level (default `0.95`).
#' @param offset         Adjust for within-subgroup mean (default `"yes"`).
#' @param pretty         Collapse each difference into a single readable string
#'                       column (e.g. `"15.2% (CI [13.6; 16.8]), p-value =
#'                       0.000"`) instead of separate `estimate`/`ci_*`/`p`
#'                       columns (default `FALSE`, which keeps the numeric
#'                       columns). In [fe_mm_pairwise_across()] one such column
#'                       is produced per subgroup pair, named after the pair;
#'                       in [fe_mm_pairwise_within()] the string lands in a
#'                       single `difference` column.
#' @param pretty_scale   Multiplier applied to the estimate and CI bounds when
#'                       `pretty = TRUE` (default `1`; use `100` for a 0-1
#'                       outcome).
#' @param pretty_suffix  String appended to each formatted number when
#'                       `pretty = TRUE` (default `"%"`).
#' @param pretty_digits  Decimal places for the estimate and CI when
#'                       `pretty = TRUE` (default `1`).
#' @param cluster_var    Optional column name used for clustering standard
#'                       errors and stratifying the bootstrap. When `NULL`
#'                       (default), clustering uses `respondent_id`. Set this
#'                       when you want to cluster by respondent but estimate
#'                       without individual fixed effects (`respondent_id = NULL`).
#' @param seed           Random seed (default `123`) used to reproduce results.
#' @param parallel       Parallelization for `boot::boot`: `"multicore"`
#'                       (Unix/Mac), `"snow"` (Windows), or `"no"`
#'                       (default `"no"`).
#' @param ncpus          Number of CPUs when `parallel != "no"`
#'                       (default `parallel::detectCores() - 1`).
#'
#' @return A tidy data frame.
#' @export
#'
#' @examples
#' \dontrun{
#' fe_mm_pairwise_within(
#'   data          = icecream,
#'   group_var     = "gender",
#'   feature_vars  = c("color", "taste", "cost"),
#'   boot_reps     = 100
#' )
#' }
fe_mm_pairwise_within <- function(data,
                                  group_var,
                                  feature_vars,
                                  outcome_var   = "choice",
                                  respondent_id = "id",
                                  cluster_var   = NULL,
                                  boot_reps     = 999,
                                  conf_level    = 0.95,
                                  offset        = "yes",
                                  pretty        = FALSE,
                                  pretty_scale  = 1,
                                  pretty_suffix = "%",
                                  pretty_digits = 1,
                                  seed          = 123,
                                  parallel      = "no",
                                  ncpus         = parallel::detectCores() - 1) suppressWarnings(suppressMessages({

  offset <- match.arg(offset, c("yes", "no"))

  for (v in feature_vars) data[[v]] <- as.factor(data[[v]])

  var_expr  <- rlang::sym(group_var)
  subgroups <- unique(data[[group_var]])
  out_list  <- vector("list", length(subgroups))
  names(out_list) <- as.character(subgroups)

  strata_var <- cluster_var %||% respondent_id

  boot_fn <- function(df, indices) {
    .mm_point(df[indices, ],
              outcome_var   = outcome_var,
              respondent_id = respondent_id,
              feature_vars  = feature_vars,
              offset        = offset,
              cluster_var   = cluster_var)$estimate
  }

  for (sg in subgroups) {
    df_sg <- dplyr::filter(data, !!var_expr == sg)

    bt <- .run_boot(seed, parallel, function() boot::boot(
      data      = df_sg,
      statistic = boot_fn,
      R         = boot_reps,
      strata    = if (!is.null(strata_var))
                    as.numeric(as.factor(df_sg[[strata_var]]))
                  else rep(1L, nrow(df_sg)),
      parallel  = parallel,
      ncpus     = if (parallel != "no") ncpus else 1L
    ))

    mm_df    <- .mm_point(df_sg,
                          outcome_var   = outcome_var,
                          respondent_id = respondent_id,
                          feature_vars  = feature_vars,
                          offset        = offset,
                          cluster_var   = cluster_var)
    bootframe        <- as.data.frame(t(bt$t))
    names(bootframe) <- paste0("V", seq_len(ncol(bootframe)))

    res_sg <- .within_pairwise(cbind(mm_df, bootframe), conf_level) |>
      dplyr::mutate(!!group_var := sg, .before = 1)

    out_list[[as.character(sg)]] <- res_sg
  }

  out <- dplyr::bind_rows(out_list)

  if (pretty && all(c("estimate", "p", "ci_low", "ci_high") %in% names(out))) {
    out <- out |>
      dplyr::mutate(difference = .fmt_pairwise(
        .data$estimate, .data$ci_low, .data$ci_high, .data$p,
        scale = pretty_scale, suffix = pretty_suffix, digits = pretty_digits)) |>
      dplyr::select(-dplyr::any_of(c("estimate", "p_pos", "p_neg", "p",
                                     "ci_low", "ci_high")))
  }

  out
}))


#' Cross-subgroup pairwise MM differences with bootstrap inference
#'
#' Computes all pairwise differences between subgroups for every attribute
#' level, with two-sided p-values.
#'
#' @inheritParams fe_mm_pairwise_within
#'
#' @return A tidy data frame.
#' @export
#'
#' @examples
#' \dontrun{
#' fe_mm_pairwise_across(
#'   data          = icecream,
#'   group_var     = "gender",
#'   feature_vars  = c("color", "taste", "cost"),
#'   boot_reps     = 100
#' )
#' }
fe_mm_pairwise_across <- function(data,
                                  group_var,
                                  feature_vars,
                                  outcome_var   = "choice",
                                  respondent_id = "id",
                                  cluster_var   = NULL,
                                  boot_reps     = 999,
                                  conf_level    = 0.95,
                                  offset        = "yes",
                                  pretty        = FALSE,
                                  pretty_scale  = 1,
                                  pretty_suffix = "%",
                                  pretty_digits = 1,
                                  seed          = 123,
                                  parallel      = "no",
                                  ncpus         = parallel::detectCores() - 1) suppressWarnings(suppressMessages({

  offset <- match.arg(offset, c("yes", "no"))
  for (v in feature_vars) data[[v]] <- as.factor(data[[v]])

  var_expr   <- rlang::sym(group_var)
  subgroups  <- unique(data[[group_var]])
  alpha      <- 1 - conf_level
  strata_var <- cluster_var %||% respondent_id

  boot_fn <- function(bd, indices) {
    .mm_point(bd[indices, ],
              outcome_var   = outcome_var,
              respondent_id = respondent_id,
              feature_vars  = feature_vars,
              offset        = offset,
              cluster_var   = cluster_var) |>
      dplyr::arrange(.data$feature, .data$level) |>
      dplyr::pull(.data$estimate)
  }

  boot_list <- lapply(stats::setNames(subgroups, subgroups), function(sg) {
    df_sg <- dplyr::filter(data, !!var_expr == sg)
    .run_boot(seed, parallel, function() boot::boot(
               data      = df_sg,
               statistic = boot_fn,
               R         = boot_reps,
               strata    = if (!is.null(strata_var))
                             as.numeric(as.factor(df_sg[[strata_var]]))
                  else rep(1L, nrow(df_sg)),
               parallel  = parallel,
               ncpus     = if (parallel != "no") ncpus else 1L))
  })

  sg_points <- lapply(stats::setNames(subgroups, subgroups), function(sg) {
    df_sg <- dplyr::filter(data, !!var_expr == sg)
    pt    <- .mm_point(df_sg,
                       outcome_var   = outcome_var,
                       respondent_id = respondent_id,
                       feature_vars  = feature_vars,
                       offset        = offset,
                       cluster_var   = cluster_var) |>
      dplyr::arrange(.data$feature, .data$level)
    nm <- paste0(as.character(sg), "_estimate")
    pt |> dplyr::rename_with(~ nm, .cols = "estimate")
  })

  sg_df <- Reduce(function(x, y) {
    y_only <- y[, setdiff(names(y), names(x)), drop = FALSE]
    cbind(x, y_only)
  }, sg_points)

  bootframes <- do.call(cbind, lapply(names(boot_list), function(nm) {
    bf        <- as.data.frame(t(boot_list[[nm]]$t))
    names(bf) <- paste0(nm, "_", paste0("V", seq_len(ncol(bf))))
    bf
  }))
  mm_full <- cbind(sg_df, bootframes)

  prefixes <- as.character(subgroups)
  pref_re  <- paste0("^(", paste(prefixes, collapse = "|"), ")_")

  res <- mm_full |>
    dplyr::select("feature", "level", dplyr::matches(pref_re)) |>
    dplyr::inner_join(
      mm_full |>
        dplyr::select("feature", "level", dplyr::matches(pref_re)),
      by = "feature", suffix = c("_1", "_2")
    ) |>
    dplyr::filter(as.character(.data$level_1) < as.character(.data$level_2)) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::matches(paste0(pref_re, ".*_1$")),
        ~ . - get(sub("_1$", "_2", dplyr::cur_column())),
        .names = "{sub('_1$', '', .col)}"
      )
    ) |>
    dplyr::select("feature", l1_level = "level_1", l2_level = "level_2",
                  dplyr::matches(pref_re)) |>
    dplyr::select(-dplyr::ends_with("_1"), -dplyr::ends_with("_2"))

  res2 <- .make_pairwise_cols(res)

  pairs <- names(res2) |>
    stringr::str_match("^(.*)_(V\\d+)$") |>
    (\(m) m[, 2])() |>
    unique() |>
    stats::na.omit() |>
    as.character()

  res3 <- res2
  for (pair in pairs) {
    vcols <- names(res3)[stringr::str_detect(
      names(res3), paste0("^", stringr::fixed(pair), "_V\\d+$"))]
    if (!length(vcols)) next

    boot_mat <- as.matrix(res3[, vcols])
    res3[[paste0(pair, "_p_pos")]]   <- 2 * rowMeans(boot_mat > 0, na.rm = TRUE)
    res3[[paste0(pair, "_p_neg")]]   <- 2 * rowMeans(boot_mat < 0, na.rm = TRUE)
    res3[[paste0(pair, "_p")]]       <- pmin(res3[[paste0(pair, "_p_pos")]],
                                             res3[[paste0(pair, "_p_neg")]])
    res3[[paste0(pair, "_ci_low")]]  <- apply(boot_mat, 1, stats::quantile,
                                              probs = alpha / 2, na.rm = TRUE, type = 7)
    res3[[paste0(pair, "_ci_high")]] <- apply(boot_mat, 1, stats::quantile,
                                              probs = 1 - alpha / 2, na.rm = TRUE, type = 7)
  }

  if (pretty) {
    for (pair in pairs) {
      res3[[pair]] <- .fmt_pairwise(
        res3[[paste0(pair, "_estimate")]],
        res3[[paste0(pair, "_ci_low")]],
        res3[[paste0(pair, "_ci_high")]],
        res3[[paste0(pair, "_p")]],
        scale = pretty_scale, suffix = pretty_suffix, digits = pretty_digits)
    }
    return(res3[, c("feature", "l1_level", "l2_level", pairs), drop = FALSE])
  }

  res3[, !grepl("_V\\d+$", names(res3))]
}))
