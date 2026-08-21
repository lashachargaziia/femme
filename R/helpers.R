#' @keywords internal

p_positive <- function(x) mean(x > 0, na.rm = TRUE)
p_negative <- function(x) mean(x < 0, na.rm = TRUE)

#' @keywords internal
.fmt_pairwise <- function(estimate, ci_low, ci_high, p,
                          scale = 1, suffix = "%", digits = 1, p_digits = 3) {
  num        <- paste0("%.", digits, "f")
  suffix_esc <- gsub("%", "%%", suffix)
  fmt <- paste0(num, suffix_esc,
                " (CI [", num, "; ", num, "]), p-value = %.", p_digits, "f")
  sprintf(fmt, estimate * scale, ci_low * scale, ci_high * scale, p)
}

#' @keywords internal
.run_boot <- function(seed, parallel, boot_call) {
  if (!identical(parallel, "no")) {
    old_kind <- RNGkind("L'Ecuyer-CMRG")
    on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)
  }
  set.seed(seed)
  boot_call()
}

#' @keywords internal
.model_obs <- function(model, n) {
  used <- seq_len(n)
  for (s in model$obs_selection) used <- used[s]
  used
}

#' @keywords internal
.mm_point <- function(df, outcome_var, respondent_id, feature_vars,
                      offset = "yes", weight_var = NULL, cluster_var = NULL) {
  mm_list <- vector("list", length(feature_vars))
  names(mm_list) <- feature_vars

  w_fml <- if (!is.null(weight_var)) stats::as.formula(paste0("~", weight_var)) else NULL

  cl <- cluster_var %||% respondent_id
  cluster_fml <- if (!is.null(cl)) stats::as.formula(paste0("~", cl)) else NULL

  nobs_val <- NULL
  obs_used <- NULL

  for (var in feature_vars) {
    df[[var]] <- as.factor(df[[var]])
    levs      <- levels(df[[var]])

    if (!is.null(respondent_id)) {
      fml <- stats::as.formula(
        paste0("`", outcome_var, "` ~ `", var, "` | `", respondent_id, "`")
      )
    } else {
      fml <- stats::as.formula(
        paste0("`", outcome_var, "` ~ `", var, "`")
      )
    }

    model <- fixest::feols(fml, data = df,
                           cluster = cluster_fml,
                           weights = w_fml)

    if (is.null(nobs_val)) {
      nobs_val <- stats::nobs(model)
      obs_used <- .model_obs(model, nrow(df))
    }

    if (!is.null(respondent_id)) {
      fe_id    <- fixest::fixef(model)[[respondent_id]]
      mean_FEs <- if (!is.null(weight_var)) {
        wts <- df[[weight_var]][match(names(fe_id), df[[respondent_id]])]
        stats::weighted.mean(fe_id, wts, na.rm = TRUE)
      } else {
        mean(fe_id, na.rm = TRUE)
      }
      coefs <- stats::coef(model)
    } else {
      all_coefs <- stats::coef(model)
      mean_FEs  <- all_coefs[["(Intercept)"]]
      coefs     <- all_coefs[names(all_coefs) != "(Intercept)"]
    }

    expected <- paste0(var, levs[-1])
    level_coefs <- vapply(expected, function(nm)
      if (nm %in% names(coefs)) as.numeric(coefs[[nm]]) else 0,
      numeric(1), USE.NAMES = FALSE)
    abs_mm   <- c(0, level_coefs) + mean_FEs
    allcoefs <- if (offset == "yes") abs_mm else abs_mm - mean(abs_mm)

    mm_list[[var]] <- data.frame(
      feature  = var,
      level    = levs,
      estimate = allcoefs,
      stringsAsFactors = FALSE
    )
  }

  out <- dplyr::bind_rows(mm_list)   # keep `feature_vars` order
  attr(out, "nobs")     <- nobs_val
  attr(out, "obs_used") <- obs_used
  out
}

#' @keywords internal
.dmm_point <- function(df, outcome_var, respondent_id, feature_vars,
                       weight_var = NULL, cluster_var = NULL) {
  dmm_list <- vector("list", length(feature_vars))
  names(dmm_list) <- feature_vars

  w_fml <- if (!is.null(weight_var)) stats::as.formula(paste0("~", weight_var)) else NULL

  cl <- cluster_var %||% respondent_id
  cluster_fml <- if (!is.null(cl)) stats::as.formula(paste0("~", cl)) else NULL

  nobs_val <- NULL

  for (var in feature_vars) {
    df[[var]] <- as.factor(df[[var]])
    levs      <- levels(df[[var]])

    if (!is.null(respondent_id)) {
      fml <- stats::as.formula(
        paste0("`", outcome_var, "` ~ `", var, "` | `", respondent_id, "`")
      )
    } else {
      fml <- stats::as.formula(
        paste0("`", outcome_var, "` ~ `", var, "`")
      )
    }

    model <- fixest::feols(fml, data = df,
                           cluster = cluster_fml,
                           weights = w_fml)

    if (is.null(nobs_val)) nobs_val <- stats::nobs(model)

    coefs    <- stats::coef(model)
    coefs    <- coefs[names(coefs) != "(Intercept)"]
    expected <- paste0(var, levs[-1])
    level_coefs <- vapply(expected, function(nm)
      if (nm %in% names(coefs)) as.numeric(coefs[[nm]]) else 0,
      numeric(1), USE.NAMES = FALSE)

    dmm_list[[var]] <- data.frame(
      feature      = var,
      level        = levs,
      estimate     = c(0, level_coefs),
      is_reference = c(TRUE, rep(FALSE, length(level_coefs))),
      stringsAsFactors = FALSE
    )
  }

  out <- dplyr::bind_rows(dmm_list)   # keep `feature_vars` order
  attr(out, "nobs") <- nobs_val
  out
}

#' @keywords internal
.within_pairwise <- function(mm_df_w_boot, conf_level = 0.95) {
  alpha  <- 1 - conf_level
  vcols  <- grep("^V\\d+$", names(mm_df_w_boot), value = TRUE)

  joined <- mm_df_w_boot |>
    dplyr::select(dplyr::all_of(c("feature", "level", "estimate", vcols))) |>
    dplyr::inner_join(
      mm_df_w_boot |>
        dplyr::select(dplyr::all_of(c("feature", "level", "estimate", vcols))),
      by = "feature", suffix = c("_1", "_2")
    ) |>
    dplyr::filter(as.character(.data$level_1) < as.character(.data$level_2))

  joined$estimate <- joined$estimate_2 - joined$estimate_1

  if (length(vcols)) {
    v1   <- paste0(vcols, "_1")
    v2   <- paste0(vcols, "_2")
    ok   <- v1 %in% names(joined) & v2 %in% names(joined)
    v1   <- v1[ok]; v2 <- v2[ok]
    newn <- vcols[ok]
    joined[newn] <- joined[v2] - joined[v1]
  }

  res <- joined |>
    dplyr::select(
      feature  = "feature",
      l1_level = "level_1",
      l2_level = "level_2",
      estimate = "estimate",
      dplyr::all_of(vcols)
    )

  if (length(vcols)) {
    boot_mat    <- as.matrix(res[, vcols])
    res$p_pos   <- 2 * rowMeans(boot_mat > 0, na.rm = TRUE)
    res$p_neg   <- 2 * rowMeans(boot_mat < 0, na.rm = TRUE)
    res$p       <- pmin(res$p_pos, res$p_neg)
    res$ci_low  <- apply(boot_mat, 1, stats::quantile,
                         probs = alpha / 2, na.rm = TRUE, type = 7)
    res$ci_high <- apply(boot_mat, 1, stats::quantile,
                         probs = 1 - alpha / 2, na.rm = TRUE, type = 7)
    res <- res[, !(names(res) %in% vcols)]
  }

  res
}

#' @keywords internal
.make_pairwise_cols <- function(df, combine = function(a, b) b - a) {
  stopifnot(is.data.frame(df))
  nm        <- names(df)
  is_target <- endsWith(nm, "_estimate") | grepl("_V\\d+$", nm)
  cols      <- nm[is_target]
  if (!length(cols)) return(df)

  suffix <- ifelse(endsWith(cols, "_estimate"), "estimate",
                   sub(".*_(V\\d+)$", "\\1", cols))
  prefix <- sub("_(estimate|V\\d+)$", "", cols)

  by_suffix_idx <- split(seq_along(cols), suffix)
  new_data      <- list()

  for (s in names(by_suffix_idx)) {
    idx          <- by_suffix_idx[[s]]
    per_pref_col <- stats::setNames(cols[idx], prefix[idx])
    per_pref_col <- per_pref_col[!duplicated(names(per_pref_col))]
    prefs        <- names(per_pref_col)
    if (length(prefs) < 2) next
    pairs <- utils::combn(prefs, 2, simplify = FALSE)
    for (pr in pairs) {
      p1 <- pr[[1]]; p2 <- pr[[2]]
      new_name             <- paste(p1, p2, s, sep = "_")
      new_data[[new_name]] <- combine(df[[per_pref_col[[p1]]]],
                                      df[[per_pref_col[[p2]]]])
    }
  }

  out           <- cbind(df[, !is_target, drop = FALSE],
                         as.data.frame(new_data, check.names = FALSE))
  rownames(out) <- NULL
  out
}

#' @keywords internal
.add_header_rows <- function(combined, level_order = NULL, reverse = TRUE,
                             header = TRUE) {
  header_of <- function(v) paste0("**", v, "**")

  var_chr  <- as.character(combined$variable_nice)
  term_chr <- as.character(combined$term_nice)

  keep     <- term_chr != header_of(var_chr)
  combined <- combined[keep, , drop = FALSE]
  var_chr  <- var_chr[keep]
  term_chr <- term_chr[keep]

  var_order <- if (is.factor(combined$variable_nice))
    levels(droplevels(combined$variable_nice)) else unique(var_chr)

  term_levels <- unlist(lapply(var_order, function(v) {
    lv <- unique(term_chr[var_chr == v])
    if (!is.null(level_order)) lv <- lv[order(match(lv, level_order))]
    if (header) c(header_of(v), lv) else lv
  }), use.names = FALSE)
  term_levels <- term_levels[!duplicated(term_levels)]

  out <- if (header) {
    headers <- combined |>
      dplyr::distinct(.data$variable_nice) |>
      dplyr::mutate(term_nice = header_of(as.character(.data$variable_nice)))
    dplyr::bind_rows(headers, combined)
  } else {
    combined
  }

  out$term_nice <- factor(as.character(out$term_nice),
                          levels = if (reverse) rev(term_levels) else term_levels)
  dplyr::arrange(out, .data$term_nice)
}
