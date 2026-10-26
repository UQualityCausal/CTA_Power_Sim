## itt_ate_count.R: power for the ITT ATE (rate ratio) on a count
## outcome, Poisson or negative binomial data, no clustering. IPTW Poisson
## regression with a robust (HC0) SE, valid for NB data (saturated model).
source("power_core.R")

#' Per-setting constants for the count-outcome simulation
#'
#' Computed once per setting, before the simulation loop: Beta parameters for
#' the PS; c_conf, which centers the PS multiplier so E[M] = 1 (so mu0 and mu1
#' stay the population means); EM2 = E[M^2]; and, for NB, the dispersion
#' kw0 / kw1 that makes each arm's total SD equal sd0 / sd1. c_conf and EM2 are
#' exact integrals over the Beta distribution (no random draws).
#'
#' @param p0     Mean propensity score (expected proportion treated).
#' @param p.cv   Propensity score CV.
#' @param family "poisson" or "nb" (negative binomial, for overdispersed counts).
#' @param mu0    Population mean of Y0 (count if untreated), outcome units.
#' @param mu1    Population mean of Y1 (count if treated), outcome units.
#'               mu1 / mu0 is the true rate ratio.
#' @param sd0    Population SD of Y0 (NB only; needs sd0^2 > mu0).
#' @param sd1    Population SD of Y1 (NB only; needs sd1^2 > mu1).
#' @param HRy    Confounding: rate ratio per 1-CV increase in the PS (see
#'               ps_cv()). > 1: higher PS -> larger counts; < 1: smaller.
#'               1 = no confounding.
#' @returns List with the inputs plus alpha, beta (Beta parameters), c_conf,
#'          EM2, kw0, kw1 (NA for Poisson) and RR_true = mu1 / mu0.
#' @examples
#' s = setup_count(p0 = 0.193, p.cv = 0.5, family = "nb", mu0 = 5.95,
#'                 mu1 = 3.71, sd0 = 3.01, sd1 = 1.97, HRy = 1 / 1.3)
setup_count = function(p0, p.cv, family = c("poisson", "nb"), mu0, mu1,
                       sd0 = NA, sd1 = NA, HRy = 1) {
  family = match.arg(family)
  if (mu0 <= 0 || mu1 <= 0) stop("mu0 and mu1 must be > 0")
  if (HRy <= 0) stop("HRy must be > 0")

  # 1. Beta parameters for the PS
  bp = get_beta_params(p0, p.cv)

  # 2. Centering constants: c_conf so E[M] = 1, EM2 = E[M^2]
  lp  = function(x) log(HRy) * ps_cv(x, p0, p.cv)
  Ebeta = function(f) integrate(function(x) f(x) * dbeta(x, bp["alpha"], bp["beta"]),
                                0, 1)$value
  c_conf = log(Ebeta(function(x) exp(lp(x))))
  EM2    = Ebeta(function(x) exp(2 * (lp(x) - c_conf)))

  # 3. NB within-person dispersion so total SD matches sd0 / sd1
  #    Var(Y) = mu + mu^2 * [(1 + kw) * EM2 - 1]
  if (family == "nb") {
    if (is.na(sd0) || is.na(sd1)) stop("family = 'nb' needs sd0 and sd1")
    kappa0 = (sd0^2 - mu0) / mu0^2               # total overdispersion, arm 0
    kappa1 = (sd1^2 - mu1) / mu1^2               # total overdispersion, arm 1
    if (kappa0 <= 0 || kappa1 <= 0) stop("sd^2 <= mean: not overdispersed; use poisson")
    kw0 = (1 + kappa0) / EM2 - 1
    kw1 = (1 + kappa1) / EM2 - 1
    if (kw0 <= 0 || kw1 <= 0)
      stop("PS effect (HRy) too large for the given SDs")
  } else {
    kw0 = kw1 = NA
  }

  list(p0 = p0, p.cv = p.cv, family = family, mu0 = mu0, mu1 = mu1,
       sd0 = sd0, sd1 = sd1, HRy = HRy,
       alpha = unname(bp["alpha"]), beta = unname(bp["beta"]),
       c_conf = c_conf, EM2 = EM2, kw0 = kw0, kw1 = kw1,
       RR_true = mu1 / mu0)
}

#' Simulate potential count outcomes (Poisson or negative binomial)
#'
#' Each person's mean is mu0 x M (untreated) or mu1 x M (treated), where
#' M = exp(log(HRy) * ps_cv(e) - c_conf) is the PS multiplier (mean 1). The
#' same M is used for Y0 and Y1; the count draws are separate.
#'
#' @param e Propensity scores, from simulate_cohort()$e.
#' @param s Setting list from setup_count().
#' @returns data.frame(Y0, Y1): count if untreated and if treated.
sim_po_count = function(e, s) {
  n = length(e)
  M = exp(log(s$HRy) * ps_cv(e, s$p0, s$p.cv) - s$c_conf)   # confounding, mean 1

  if (s$family == "poisson") {
    Y0 = rpois(n, lambda = s$mu0 * M)
    Y1 = rpois(n, lambda = s$mu1 * M)
  } else {
    Y0 = rnbinom(n, size = 1 / s$kw0, mu = s$mu0 * M)
    Y1 = rnbinom(n, size = 1 / s$kw1, mu = s$mu1 * M)
  }

  data.frame(Y0 = Y0, Y1 = Y1)
}

#' IPTW estimate and robust SE for the log rate ratio
#'
#' Weighted Poisson regression Y ~ X with IPT weights and a robust (sandwich,
#' HC0) SE. Valid for NB data too: Y ~ X is saturated, so the estimate is the
#' ratio of weighted means either way, and the robust SE handles
#' overdispersion. Same estimate and SE as survey::svyglm(family =
#' quasipoisson()) (up to sqrt(n/(n-1))), but faster on large data.
#'
#' @param d Data frame from combine_outcomes(), with Y (count), X (0/1
#'          treatment) and W (IPT weight).
#' @returns c(est = , se = ): estimated log rate ratio and its robust SE.
#'          exp(est) is the rate ratio.
est_itt_ate_count = function(d) {
  fit = suppressWarnings(glm(Y ~ X, family = poisson, data = d, weights = W))
  # (warnings are from non-integer weighted counts; harmless here)

  est = coef(fit)["X"]                                          # log RR
  se  = sqrt(sandwich::vcovHC(fit, type = "HC0")["X", "X"])     # robust SE

  c(est = unname(est), se = unname(se))
}

#' Simulate nsim datasets and compute power for the IPTW rate ratio
#'
#' @param n         Sample size.
#' @param p0        Expected proportion treated (allocation; 1:3 -> 1/4).
#' @param p.cv      Propensity score CV.
#' @param family    "poisson" or "nb": how counts are generated.
#' @param mu0, mu1  Population means of Y0 / Y1 (outcome units).
#' @param sd0, sd1  Population SDs of Y0 / Y1 (NB only).
#' @param HRy       Confounding: rate ratio per 1-CV increase in the PS.
#'                  1 = no confounding.
#' @param nsim      Number of simulated datasets.
#' @param alpha_adj Two-sided significance level after multiplicity adjustment.
#' @returns One-row data.frame:
#'   RR_true   mu1 / mu0
#'   mean_est  average log RR estimate
#'   marg_RR   exp(mean_est): close to RR_true = no bias
#'   crude_RR  exp(mean unweighted log ratio of arm means): shows confounding
#'   emp_sd    SD of the log RR estimates (the true SE)
#'   mean_se   average robust SE (close to emp_sd = SE is accurate)
#'   power     share of datasets that reject H0: RR = 1
#' See the examples at the bottom of this file.
power_itt_ate_count = function(n, p0, p.cv, family = c("poisson", "nb"),
                               mu0, mu1, sd0 = NA, sd1 = NA, HRy = 1,
                               nsim = 1000, alpha_adj = 0.05) {
  family = match.arg(family)
  s    = setup_count(p0, p.cv, family, mu0, mu1, sd0, sd1, HRy)  # once per setting
  crit = qnorm(1 - alpha_adj / 2)

  res = t(replicate(nsim, {
    co = simulate_cohort(n, s$alpha, s$beta, estimand = "ATE")
    po = sim_po_count(co$e, s)
    d  = combine_outcomes(co, po)
    c(est_itt_ate_count(d),
      crude = log(mean(d$Y[d$X == 1]) / mean(d$Y[d$X == 0])))
  }))

  data.frame(RR_true  = s$RR_true,
             mean_est = mean(res[, "est"]),
             marg_RR  = exp(mean(res[, "est"])),
             crude_RR = exp(mean(res[, "crude"])),
             emp_sd   = sd(res[, "est"]),
             mean_se  = mean(res[, "se"]),
             power    = mean(abs(res[, "est"] / res[, "se"]) > crit))
}

# Examples
if (FALSE){
## 1. No effect: power should be about alpha_adj
set.seed(1)
power_itt_ate_count(n = 2000, p0 = 1/4, p.cv = 0.25, family = "poisson",
                    mu0 = 2, mu1 = 2, HRy = 1.3, nsim = 1000)
#   marg_RR ~1   power ~0.05 (0.035-0.065 is fine)

## 2. Poisson, no confounding: mean_se should match the textbook SE
##    sqrt(1/(n*p0*mu1) + 1/(n*(1-p0)*mu0)) = 0.0389
set.seed(2)
power_itt_ate_count(n = 2000, p0 = 1/4, p.cv = 0.25, family = "poisson",
                    mu0 = 2, mu1 = 1.7, HRy = 1, nsim = 500)
#   marg_RR ~0.85   emp_sd ~ mean_se ~0.040   power ~0.98

## 3. NB setting (antibiotics Aim 3), confounding: higher PS -> fewer days
set.seed(3)
power_itt_ate_count(n = 3908, p0 = 0.193, p.cv = 0.5, family = "nb",
                    mu0 = 5.95, mu1 = 3.71, sd0 = 3.01, sd1 = 1.97,
                    HRy = 1 / 1.3, nsim = 500)
#   marg_RR ~0.62 (= RR_true)   crude_RR ~0.60 (confounded)   power ~1

## 4. Power curve over the rate ratio (NB, smaller effects)
rrs   = c(0.90, 0.93, 0.95, 0.97)
curve = do.call(rbind, lapply(rrs, function(r) {
  set.seed(4)
  power_itt_ate_count(n = 3908, p0 = 0.193, p.cv = 0.5, family = "nb",
                      mu0 = 5.95, mu1 = 5.95 * r, sd0 = 3.01, sd1 = 3.01 * r,
                      HRy = 1 / 1.3, nsim = 200)
}))
curve[, c("RR_true", "marg_RR", "mean_se", "power")]
#   power ~0.99, 0.85, 0.59, 0.26
}
