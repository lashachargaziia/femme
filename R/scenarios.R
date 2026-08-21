#' Predict scenario-level choice probabilities.
#'
#' @param data          Conjoint data frame with. Requires a pre-computed weight column.
#' @param scenarios     Named list of scenario vectors.
#' @param group_label   String label added to every output row.
#' @param weight_var    Name of the pre-computed weight column (default `"wt"`).
#' @param outcome_var   Name of the outcome column (default `"choice"`). Can
#'                       be both continuous and binary.
#' @param respondent_id Name of the unique observation identifier
#'                       (default `"id"`). Set to `NULL` to estimate without
#'                       respondent fixed effects (plain weighted OLS with
#'                       unstratified bootstrap).
#' @param feature_vars  Character vector of conjoint attribute columns.
#' @param boot_reps     Bootstrap replications (default `999`).
#' @param order         Interaction order (default `1`). `1` fits main effects
#'                      only; `2` adds all pairwise interactions, and so on.
#' @param shrinkage     `"LASSO"` or `"none"` (default `"none"`).
#' @param cluster_var   Optional column name used for clustering standard
#'                      errors and stratifying the bootstrap. When `NULL`
#'                      (default), clustering uses `respondent_id`. Set this
#'                      when you want to cluster by respondent but estimate
#'                      without individual fixed effects (`respondent_id = NULL`).
#' @param conf_level    Confidence level (default `0.95`).
#' @param lasso_seed    Seed for LASSO CV (default `123`).
#' @param seed          Seed for bootstrap reproducibility (default `123`).
#' @param parallel      Parallelization for `boot::boot`:
#'                      `"multicore"` (Unix/Mac), `"snow"` (Windows),
#'                      or `"no"` (default `"no"`).
#' @param ncpus         Number of CPUs when `parallel != "no"`
#'                      (default `parallel::detectCores() - 1`).
#' @return Data frame of scenario estimates and pairwise differences.
#' @export
predict_preference <- function(data,
                              scenarios,
                              group_label,
                              weight_var      = "wt",
                              outcome_var     = "choice",
                              respondent_id   = "id",
                              cluster_var     = NULL,
                              feature_vars,
                              boot_reps       = 999,
                              order           = 1,
                              shrinkage       = "none",
                              conf_level      = 0.95,
                              lasso_seed      = 123,
                              seed            = 123,
                              parallel        = "no",
                              ncpus           = parallel::detectCores() - 1) suppressWarnings(suppressMessages({

  .fast_estimate <- function(model, weights, intercept) {
    coefs  <- stats::coef(model)
    names(coefs) <- gsub("`", "", names(coefs))
    shared <- intersect(names(weights), names(coefs))
    sum(weights[shared] * coefs[shared]) + intercept
  }

  stopifnot(weight_var %in% names(data))
  data <- dplyr::filter(data, !is.na(.data[[weight_var]]))

  for (v in feature_vars) if (!is.factor(data[[v]])) data[[v]] <- as.factor(data[[v]])

  ord_keys <- c(respondent_id, outcome_var, feature_vars)
  ord_keys <- ord_keys[!vapply(ord_keys, is.null, logical(1))]
  data <- data[do.call("order", lapply(ord_keys, function(k) data[[k]])), , drop = FALSE]
  rownames(data) <- NULL

  weight_fml  <- stats::as.formula(paste0("~`", weight_var, "`"))
  cl <- cluster_var %||% respondent_id
  cluster_fml <- if (!is.null(cl))
    stats::as.formula(paste0("~`", cl, "`")) else NULL
  strata_var  <- cluster_var %||% respondent_id

  if (shrinkage == "LASSO") {


    if (!requireNamespace("Matrix", quietly = TRUE) ||
        !requireNamespace("glmnet", quietly = TRUE))
      stop('Packages "Matrix" and "glmnet" are required for shrinkage = "LASSO".\n',
           'Install them with: install.packages(c("Matrix", "glmnet"))',
           call. = FALSE)

    .to_weights <- function(scenario_vec, coefs, sel_vars, safe_names) {
      stopifnot(!is.null(names(scenario_vec)), length(scenario_vec) > 0)

      split_vl <- function(nm) {
        hit <- feature_vars[startsWith(nm, feature_vars)]
        if (length(hit) != 1) return(NULL)
        list(var = hit, level = sub(paste0("^", hit), "", nm))
      }

      baseline   <- vapply(feature_vars, function(v) levels(data[[v]])[1], character(1))
      scen       <- as.list(baseline)
      main_mult  <- list()
      inter_mult <- list()

      for (nm in names(scenario_vec)) {
        if (grepl(":", nm, fixed = TRUE)) {
          inter_mult[[nm]] <- as.numeric(scenario_vec[[nm]])
        } else {
          pl <- split_vl(nm)
          if (!is.null(pl)) {
            if (!pl$level %in% levels(data[[pl$var]]))
              stop("Level '", pl$level, "' not found for '", pl$var, "'.")
            scen[[pl$var]]  <- pl$level
            main_mult[[nm]] <- as.numeric(scenario_vec[[nm]])
          }
        }
      }

      scen <- as.data.frame(scen, stringsAsFactors = FALSE)
      for (v in feature_vars) scen[[v]] <- factor(scen[[v]], levels = levels(data[[v]]))

      f <- if (order == 1) {
        stats::as.formula(paste0("~ (", paste(feature_vars, collapse = " + "), ") - 1"))
      } else {
        stats::as.formula(paste0("~ (", paste(feature_vars, collapse = " + "), ")^", order, " - 1"))
      }

      Xrow <- stats::model.matrix(f, data = scen)
      Xsel <- Xrow[, sel_vars, drop = FALSE]

      for (nm in names(main_mult))
        if (nm %in% colnames(Xsel)) Xsel[, nm] <- Xsel[, nm] * main_mult[[nm]]

      for (ic in colnames(Xsel)[grepl(":", colnames(Xsel), fixed = TRUE)]) {
        mult <- if (!is.null(inter_mult[[ic]])) {
          inter_mult[[ic]]
        } else {
          parts <- strsplit(ic, ":", fixed = TRUE)[[1]]
          prod(sapply(parts, function(p) if (!is.null(main_mult[[p]])) main_mult[[p]] else 1))
        }
        Xsel[, ic] <- Xsel[, ic] * mult
      }

      colnames(Xsel) <- safe_names
      keep <- colnames(Xsel) %in% coefs
      if (!any(keep)) return(stats::setNames(0, coefs[1]))
      w <- as.numeric(Xsel[1, keep, drop = FALSE])
      names(w) <- colnames(Xsel)[keep]
      w
    }

    f <- stats::as.formula(
      if (order == 1)
        paste0("~ (", paste(feature_vars, collapse = " + "), ") - 1")
      else
        paste0("~ (", paste(feature_vars, collapse = " + "), ")^", order, " - 1")
    )
    X <- stats::model.matrix(f, data = data)

    y <- data[[outcome_var]]
    w <- as.numeric(data[[weight_var]])

    if (!is.null(respondent_id)) {
      S   <- Matrix::sparse.model.matrix(
        stats::as.formula(paste0("~ 0 + ", respondent_id)), data = data)
      gw  <- as.numeric(Matrix::t(S) %*% w)
      iGw <- 1 / gw
      gWy <- as.numeric(Matrix::t(S) %*% (w * y))
      gWX <- Matrix::t(S) %*% (Matrix::Diagonal(x = w) %*% X)

      set.seed(lasso_seed)
      uid    <- sort(unique(as.character(data[[respondent_id]])))
      u_fold <- sample(rep(seq_len(10), length.out = length(uid)))
      foldid <- u_fold[match(as.character(data[[respondent_id]]), uid)]
      cvfit <- glmnet::cv.glmnet(
        x           = X - S %*% Matrix::Diagonal(x = iGw) %*% gWX,
        y           = y - as.numeric(S %*% (iGw * gWy)),
        family      = "gaussian", alpha = 1, foldid = foldid,
        standardize = TRUE, intercept = FALSE
      )
    } else {
      set.seed(lasso_seed)
      cvfit <- glmnet::cv.glmnet(
        x           = X,
        y           = y,
        weights     = w,
        family      = "gaussian", alpha = 1, nfolds = 10,
        standardize = TRUE, intercept = TRUE
      )
    }
    b        <- as.matrix(stats::coef(cvfit, s = "lambda.1se"))
    sel_vars <- colnames(X)[which(b[-1, 1] != 0)]
    stopifnot(length(sel_vars) > 0)

    n_lasso_vars <- length(sel_vars)

    X_sel      <- as.data.frame(as.matrix(X[, sel_vars, drop = FALSE]))
    safe_names <- make.names(colnames(X_sel), unique = TRUE)
    colnames(X_sel) <- safe_names
    data_model <- cbind(data, X_sel)

    fml <- if (!is.null(respondent_id)) {
      stats::as.formula(
        paste0(outcome_var, " ~ ", paste(safe_names, collapse = " + "), " | ", respondent_id)
      )
    } else {
      stats::as.formula(
        paste0(outcome_var, " ~ ", paste(safe_names, collapse = " + "))
      )
    }

    pilot_coefs <- names(stats::coef(fixest::feols(
      fml, data = data_model, cluster = cluster_fml,
      weights = weight_fml, fixef.rm = "none"
    )))

    precomputed_weights <- lapply(scenarios, .to_weights,
                                  coefs      = pilot_coefs,
                                  sel_vars   = sel_vars,
                                  safe_names = safe_names)

  } else {

    n_lasso_vars <- NA_integer_

    fml_rhs <- paste(paste0("`", feature_vars, "`"), collapse = " + ")
    fml_str <- if (order == 1) {
      paste0("`", outcome_var, "` ~ ", fml_rhs)
    } else {
      paste0("`", outcome_var, "` ~ (", fml_rhs, ")^", order)
    }
    if (!is.null(respondent_id)) fml_str <- paste0(fml_str, " | `", respondent_id, "`")
    fml <- stats::as.formula(fml_str)

    pilot_coefs <- gsub("`", "", names(stats::coef(fixest::feols(
      fml, data = data, cluster = cluster_fml,
      weights = weight_fml, fixef.rm = "none"
    ))))

    precomputed_weights <- lapply(scenarios, function(sv) {
      main_mult  <- list()
      inter_mult <- list()
      for (nm in names(sv)) {
        if (grepl(":", nm, fixed = TRUE)) {
          inter_mult[[nm]] <- as.numeric(sv[[nm]])
        } else {
          main_mult[[nm]] <- as.numeric(sv[[nm]])
        }
      }
      w <- numeric(0)
      for (cn in pilot_coefs) {
        if (grepl(":", cn, fixed = TRUE)) {
          if (!is.null(inter_mult[[cn]])) {
            w[cn] <- inter_mult[[cn]]
          } else {
            parts <- strsplit(cn, ":", fixed = TRUE)[[1]]
            if (all(parts %in% names(main_mult)))
              w[cn] <- prod(vapply(parts, function(p) main_mult[[p]], numeric(1)))
          }
        } else if (cn %in% names(main_mult)) {
          w[cn] <- main_mult[[cn]]
        }
      }
      w
    })

    data_model <- data
  }

  boot_fn <- function(d, indices) {
    rd <- d[indices, ]
    m  <- fixest::feols(
      fml, data = rd, cluster = cluster_fml,
      weights = weight_fml, fixef.rm = "none", keep_model = TRUE
    )

    FE_mean <- if (!is.null(respondent_id)) {
      fe_vals     <- fixest::fixef(m)[[respondent_id]]
      resp_counts <- table(rd[[respondent_id]])
      stats::weighted.mean(fe_vals,
                           as.numeric(resp_counts[names(fe_vals)]),
                           na.rm = TRUE)
    } else {
      stats::coef(m)[["(Intercept)"]]
    }

    vapply(precomputed_weights, .fast_estimate, numeric(1),
           model = m, intercept = FE_mean)
  }

  bootstrap_res <- .run_boot(seed, parallel, function() boot::boot(
    data      = data_model,
    statistic = boot_fn,
    R         = boot_reps,
    strata    = if (!is.null(strata_var))
                  as.numeric(as.factor(data_model[[strata_var]]))
                else rep(1L, nrow(data_model)),
    parallel  = parallel,
    ncpus     = if (parallel != "no") ncpus else 1L
  ))

  alpha          <- 1 - conf_level
  K              <- length(scenarios)
  scenario_names <- names(scenarios)
  boot_estimates <- colMeans(bootstrap_res$t, na.rm = TRUE)

  ci_bounds <- do.call(rbind, lapply(seq_len(K), function(i) {
    boot::boot.ci(bootstrap_res, type = "perc", index = i)$perc[4:5]
  }))

  scenario_tbl <- data.frame(
    scenario   = factor(scenario_names, levels = scenario_names),
    estimate   = boot_estimates,
    lower_CI   = ci_bounds[, 1],
    upper_CI   = ci_bounds[, 2],
    p          = NA_real_,
    comparison = NA_character_,
    base       = 1L,
    stringsAsFactors = FALSE
  )

  diff_tbl <- dplyr::bind_rows(lapply(
    utils::combn(seq_len(K), 2, simplify = FALSE),
    function(idx) {
      i <- idx[[1]]; j <- idx[[2]]
      draws <- bootstrap_res$t[, j] - bootstrap_res$t[, i]
      ci    <- stats::quantile(draws, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, type = 7)
      data.frame(
        scenario   = factor(paste0(scenario_names[j], " - ", scenario_names[i])),
        estimate   = mean(draws, na.rm = TRUE),
        lower_CI   = as.numeric(ci[1]),
        upper_CI   = as.numeric(ci[2]),
        p          = min(2 * p_positive(draws), 2 * p_negative(draws)),
        comparison = paste0(scenario_names[j], " minus ", scenario_names[i]),
        base       = 0L
      )
    }
  ))

  dplyr::bind_rows(scenario_tbl, diff_tbl) |>
    dplyr::mutate(group = group_label, n_lasso_vars = n_lasso_vars)
}))
