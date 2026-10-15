## power_core.R: shared functions for CTA power simulations
## Defines functions only: no runs, no seeds, no file output.

#' Beta parameters for simulated propensity scores
#' 
#' Propensity scores are simulated from a Beta distribution because it is
#' bounded between 0 and 1. Given the mean propensity score (expected proportion
#' treated) and its coefficient of variation (CV), this function returns the
#' matching Beta(alpha, beta) parameters and checks that the specification is 
#' possible. 
#' 
#' @param p0   Mean propensity score (expected proportion treated), in (0, 1).
#' @param p.cv Coefficient of variation of the propensity score, >0.  
#'             Must be below sqrt((1-p0) / p0).
#' @returns Named numeric vector c(alpha = , beta = ). 
#' @examples 
#' # 1/3 treated, moderate variation in the propensity score
#' get_beta_params(1/3, 0.5)
get_beta_params = function(p0, p.cv){
  # input checks
  if (p0 <= 0 || p0 >= 1) stop("p0 must be strictly between 0 and 1")
  if (p.cv <= 0) stop("p.cv must be > 0")
  
  max_cv = sqrt((1-p0)/p0)
  if (p.cv >= max_cv)
    stop(sprintf("p.cv = %.2f is impossible for p0 = %.2f (must be < %.2f)", 
                 p.cv, p0, max_cv))

  s2 = (p.cv * p0)^2 # target variance
  alpha = p0*(p0*(1-p0)/s2 - 1)
  beta = alpha*(1-p0)/p0
  
  if (alpha < 1 || beta < 1)
    warning(sprintf(paste0("alpha = %.2f, beta = %.2f: propensity scores pile up ",
                           "near 0 or 1, so expect extreme weights"), alpha, beta))
  
  c(alpha = alpha, beta = beta)
}

#' Simulate a cohort with propensity scores, treatment, and IPT weights
#' 
#' Draw propensity score from Beta(alpha, beta), assigns treatment as Bernoulli(e),
#' and computes inverse probability of treatment weights for the chosen estimand
#' 
#' @param n        Cohort size.
#' @param alpha    First shape parameter of the propensity score Beta distribution.
#' @param beta     Second shape parameter.
#' @param estimand "ATE" (average treatment effect) or "ATT" (effect on the treated).
#' @returns Data frame with columns e (propensity score), X (0/1 treatment),
#' and W (weight).
#' @examples
#' p = get_beta_params(1/3, 0.5)
#' co = simulate_cohort(100, p["alpha"], p["beta"], estimand = "ATE")
simulate_cohort = function(n, alpha, beta, estimand = c("ATE", "ATT")){
  estimand = match.arg(estimand)

  e_raw = rbeta(n, alpha, beta)
  n_clamped = sum(e_raw < 1e-6 | e_raw > 1 - 1e-6)
  if (n_clamped > 0)
    warning(sprintf("%d propensity score(s) within 1e-6 of 0 or 1 were clamped", n_clamped))
  e = pmin(pmax(e_raw, 1e-6), 1 - 1e-6)  # avoid infinite weights
  X = rbinom(n, 1, e)
  
  if (estimand == "ATE"){
    W = X/e + (1-X)/(1-e)
  } else {
    W = X + (1-X)*e/(1-e)
  }
  
  data.frame(e = e, X = X, W = W)
}

#' Standardized propensity score
#'
#' Puts the propensity score on an SD scale (mean 0, SD 1), so confounding
#' strength means the same thing whatever p0 and p.cv are.
#'
#' @param e    Propensity scores.
#' @param p0   Mean propensity score (same value given to get_beta_params()).
#' @param p.cv Propensity score CV (same value given to get_beta_params()).
#' @returns Numeric vector z = (e - p0) / (p.cv * p0).
ps_z = function(e, p0, p.cv) {
  (e - p0) / (p.cv * p0)
}

#' Combine potential outcomes into observed outcomes
#' 
#' For every pair of columns <name>0 / <name>1 in po, creates the observed 
#' column <name> = <name>1 if treated (X = 1), <name>0 otherwise. Columns in po
#' that are not part of a pair (e.g., baseline Yb, offset log_fu) are carried
#' along unchanged.
#' 
#' @param cohort  Data frame from simulate_cohort, with 0/1 column X.
#' @param po      Data frame of potential outcomes, same number of rows as cohort.
#' @param keep_po If true, also keep the potential outcome columns (useful for
#'                checking the true effect, e.g., mean(Yf1 - Yf0)).
#' @returns cohort with observed outcome columns (and other po columns) added.
#' @examples 
#' p = get_beta_params(1/3, 0.5)
#' co = simulate_cohort(100, p["alpha"], p["beta"], estimand = "ATE")
#' po = data.frame(Yb = rnorm(100), Yf0 = rnorm(100), Yf1 = rnorm(100))
#' d = combine_outcomes(co, po)
combine_outcomes = function(cohort, po, keep_po = FALSE) {
  # input checks
  if (nrow(cohort) != nrow(po)){
    stop("cohort and po must have the same number of rows")
  }
  if (!all(cohort$X %in% c(0,1))){
    stop("cohort$X must be coded 0/1")
  }
  
  # find complete pairs: <name>0 and <name>1 both present -> potential outcomes
  nm0 = grep("0$", names(po), value = TRUE)
  base = sub("0$", "", nm0)
  base = base[paste0(base, "1") %in% names(po)]
  if (length(base) == 0){
    stop("po has no <name>0 / <name>1 column pairs")
  }
  
  pair_cols = c(paste0(base, "0"), paste0(base, "1"))
  other = setdiff(names(po), pair_cols)
  
  # don't silently overwrite cohort columns
  new_cols = c(base, other, if (keep_po) pair_cols)
  clash = intersect(new_cols, names(cohort))
  if (length(clash) > 0){
    stop("column name(s) already in cohort: ", paste(clash, collapse = ", "))
  }
  
  # observed outcome = potential outcome under the treatment actually received
  obs = lapply(base, function(b)
    ifelse(cohort$X == 1, po[[paste0(b, "1")]], po[[paste0(b, "0")]]))
  names(obs) = base
  
  out = cbind(cohort, po[other], as.data.frame(obs))
  if (keep_po){
    out = cbind(out, po[pair_cols])
  }
  out
}
