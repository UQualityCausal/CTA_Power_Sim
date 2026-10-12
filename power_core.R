#' Beta parameters for simulated propensity scores
#' 
#' Propensity scores are simulated from a Beta distribution because it is
#' bounded between 0 and 1. Given the mean propensity score (expected proportion
#' treated) and its coefficient of variation (CV), this function returns the
#' matching Beta(alpha, beta) parameters and checks that the specification is 
#' possible. 
#' 
#' @param p0 Mean propensity score (expected proportion treated), in (0, 1).
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
  if (p.cv > max_cv)
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

#' Simulate a cohort with propensity scores, treatment, and IPT wights
#' 
#' Draw propensity score from Beta(alpha, beta), assigns treatment as Bernoulli(e),
#' and computes inverse probability of treatment weights for the chosen estimand
#' 
#' @param n Cohort size.
#' @param alpha First shape parameter of the propensity score Beta distribution.
#' @param beta Second shape parameter.
#' @param estimand "ATE" (average treatment effect) or "ATT" (effect on the treated).
#' @returns Data frame with with columns e (propensity score), X (0/1 treatment),
#' and W (weight).
#' @examples
#' p = get_beta_params(1/3, 0.5)
#' d = simulate_cohort(100, p["alpha"], p["beta"], estimand = "ATE")
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
