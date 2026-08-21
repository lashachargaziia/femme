#' femme: Marginal Means and AMCE Estimation for Conjoint Experiments
#'
#' @description
#' `femme` provides a complete workflow for analysing conjoint survey
#' experiments:
#'
#' * **AMCE estimation** – [fe_amce()] / [fe_amce_by_group()]
#' * **Marginal Means** – [fe_mm()] / [fe_mm_by_group()] (pass `weight_var` for weighted estimation)
#' * **Differences-in-Marginal-Means** – [fe_dmm()] / [fe_dmm_by_group()]
#' * **Within-subgroup pairwise inference** – [fe_mm_pairwise_within()]
#' * **Cross-subgroup pairwise inference** – [fe_mm_pairwise_across()]
#' * **Scenario predictions** – [predict_preference()]
#'
#' All estimators use `fixest::feols` for fast fixed-effect estimation and 
#' stratified bootstrapping for inference.
#'
#' @section Data format:
#' All functions expect data in long format: one row per profile shown,
#' with columns for each attribute, an outcome, and a respondent
#' identifier.
#'
#' @section Variable lookup:
#' Several functions accept a `variable_lookup` argument — a two-column
#' `data.frame` with `variable` (the column name in the data) and
#' `variable_nice` (the user-defined label).
#'
#' @keywords internal
"_PACKAGE"

## suppress R CMD CHECK notes for non-standard column references
utils::globalVariables(c(
  ".", ".data", ":=", "feature", "level", "estimate", "variable", "variable_nice",
  "term_nice", "term", "std.error", "p.value", "nobs",
  "p_pos", "p_neg", "p", "ci_low", "ci_high",
  "boot_estimate", "boot_lower", "boot_upper", "is_reference", "p_value",
  "boot_se", "diff_from_mean", "z_from_mean", "p_from_mean",
  "mean_outcome", "group", "offset",
  "l1_level", "l2_level", "level_1", "level_2",
  "estimate_1", "estimate_2"
))
