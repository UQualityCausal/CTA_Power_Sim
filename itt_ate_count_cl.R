## itt_ate_count_cl.R: power for the ITT ATE (rate ratio) on a count
## outcome clustered by hospital. IPTW Poisson or NB regression with a
## hospital-clustered robust SE; superiority, non-inferiority, inferiority
## or equivalence tests via decide().
source("power_core.R")

#' Per-setting constants for the count-outcome simulation
#'
#' Computes, once per setting: Beta parameters for the PS; the PS-multiplier
#' constants c_conf (so E[M_ps] = 1) and EM2 = E[M_ps^2]; the hospital variance
#' tau2 for the target ICC; and, for NB, the within dispersion kw0/kw1 so each
#' arm's total SD matches sd0/sd1.
#'
#' @param p0, p.cv Propensity score mean and CV.
#' @param family   "poisson" or "nb".
#' @param mu0, mu1 Population means of Y0/Y1.
#' @param sd0, sd1 SD of the outcome by arm (NB only; need sd^2 > mean).
#' @param HRy      Size of PS effect on the outcome: ratio per 1-CV of PS.
#' @param sgn      Direction: -1 = higher PS -> smaller outcome; +1 = larger.
#' @param icc      Hospital ICC on the count scale, calibrated on the untreated
#'                 arm (0 = no clustering).
#' @param ndraw    Monte Carlo draws used for c_conf and EM2.
#' @returns List with the inputs (p0, p.cv, family, mu0, mu1, HRy, sgn, icc)
#'          and alpha, beta, c_conf, EM2, tau2, kw0, kw1 (NA for Poisson),
#'          RR_true = mu1 / mu0.
setup_count_cl = function(p0, p.cv, family = c("poisson", "nb"), mu0, mu1,
                       sd0 = NA, sd1 = NA, HRy = 1, sgn = 1, icc = 0,
                       ndraw = 1e6) {
  family = match.arg(family)
  
  # 1. Beta parameters for the PS
  bp = get_beta_params(p0, p.cv)
  if (bp["alpha"] <= 1 || bp["beta"] <= 1)
    stop("IPW weights have infinite mean: reduce p.cv")
  
  # 2. PS multiplier constants: c_conf so E[M_ps] = 1; EM2 = E[M_ps^2]
  e_big  = rbeta(ndraw, bp["alpha"], bp["beta"])
  lp_big = sgn * log(HRy) * ps_cv(e_big, p0, p.cv)
  c_conf = log(mean(exp(lp_big)))
  EM2    = mean(exp(2 * (lp_big - c_conf)))
  
  # 3. Hospital variance tau2 and NB within-dispersion kw
  if (family == "nb") {
    if (is.na(sd0) || is.na(sd1)) stop("family = 'nb' needs sd0 and sd1")
    kappa0 = (sd0^2 - mu0) / mu0^2                 # total dispersion, arm 0
    kappa1 = (sd1^2 - mu1) / mu1^2                 # total dispersion, arm 1
    if (kappa0 <= 0 || kappa1 <= 0) stop("sd^2 <= mean: not overdispersed")
    tau2 = log(1 + icc * (1/mu0 + kappa0))         # ICC calibrated on arm 0
    kw0  = (1 + kappa0) / (EM2 * exp(tau2)) - 1
    kw1  = (1 + kappa1) / (EM2 * exp(tau2)) - 1
    if (kw0 <= 0 || kw1 <= 0)
      stop("ICC or PS effect too large for the observed dispersion")
  } else {
    if (icc * EM2 >= 1) stop("ICC too large")
    tau2 = log((1 + icc * (1/mu0 - 1)) / (1 - icc * EM2))
    kw0 = kw1 = NA
  }
  
  list(p0 = p0, p.cv = p.cv, family = family, mu0 = mu0, mu1 = mu1,
       HRy = HRy, sgn = sgn, icc = icc,
       alpha = unname(bp["alpha"]), beta = unname(bp["beta"]),
       c_conf = c_conf, EM2 = EM2, tau2 = tau2, kw0 = kw0, kw1 = kw1,
       RR_true = mu1 / mu0)
}

## ==== Decision rules ====

#' Decision rule for one or more simulated datasets
#'
#' Compares one-sided confidence bounds for est with thresholds. All
#' quantities are on the estimate's scale (log RR / log HR for ratios; raw
#' difference for means). Works on vectors.
#'
#' @param est          Estimate(s) (e.g. log RR).
#' @param se           Standard error(s).
#' @param test         "superiority2" (two-sided, H1: effect != 0),
#'                     "superiority" (one-sided, H1: better),
#'                     "noninferiority" (H1: not worse by more than margin),
#'                     "inferiority" (H1: worse; by more than margin if margin > 0),
#'                     "equivalence" (H1: within +/- margin; TOST).
#' @param alpha        Two-sided level for "superiority2"; one-sided level per
#'                     bound for the others (e.g. 0.025).
#' @param margin       Margin on the estimate's scale (e.g. log(1.1)); >= 0.
#'                     Ignored for "superiority2" and "superiority".
#' @param lower_better TRUE if smaller outcomes are better.
#' @returns Logical vector: TRUE where the test's alternative is declared.
decide = function(est, se, test = c("superiority2", "superiority", "noninferiority",
                                    "inferiority", "equivalence"),
                  alpha = 0.05, margin = 0, lower_better = TRUE) {
  test = match.arg(test)
  b = if (lower_better) est else -est           # b > 0 means "worse"
  
  if (test == "superiority2") return(abs(b / se) > qnorm(1 - alpha / 2))
  
  z = qnorm(1 - alpha)
  upper = b + z * se
  lower = b - z * se
  
  switch(test,
         superiority    = upper < 0,
         noninferiority = upper < margin,
         inferiority    = lower > margin,
         equivalence    = lower > -margin & upper < margin)
}

#' Simulate hospital membership and shared hospital multipliers
#'
#' Each person is assigned to one of H hospitals at random. Each hospital gets
#' one multiplier M_h = exp(u), u ~ N(-tau2/2, tau2), so E[M_h] = 1. Everyone
#' in a hospital shares the same M_h, which creates within-hospital correlation.
#'
#' @param n    Number of people.
#' @param H    Number of hospitals (whole number >= 1).
#' @param tau2 Variance of the log hospital effect, e.g. setup_count_cl()$tau2.
#'             0 = no clustering (all M_h = 1).
#' @returns data.frame(hosp, M_h): hospital ID (1..H) and that hospital's
#'          multiplier, one row per person.
sim_cluster = function(n, H, tau2 = 0) {
  # input checks
  if (H < 1 || H != round(H)) stop("H must be a whole number >= 1")
  if (tau2 < 0) stop("tau2 must be >= 0")
  
  hosp = sample(1:H, n, replace = TRUE)               # hospital for each person
  u    = rnorm(H, mean = -tau2 / 2, sd = sqrt(tau2))  # one effect per hospital
  
  data.frame(hosp = hosp, M_h = exp(u[hosp]))         # look up own hospital's effect
}

#' Simulate potential count outcomes (Poisson or negative binomial)
#'
#' Each person's mean is mu x M, with M = PS multiplier x hospital multiplier
#' (both mean 1). The same M is used for Y0 and Y1 (same person, same
#' hospital); the count draws are separate.
#'
#' @param e   Propensity scores, from simulate_cohort()$e.
#' @param M_h Hospital multiplier per person, from sim_cluster()$M_h.
#' @param s   Setting list from setup_count_cl() (family, mu0, mu1, HRy, sgn, p0,
#'            p.cv, c_conf, kw0, kw1).
#' @returns data.frame(Y0, Y1): count if X = 0 and if X = 1.
sim_po_count_cl = function(e, M_h, s) {
  if (length(M_h) != length(e)) stop("e and M_h must have the same length")
  n = length(e)
  
  # PS multiplier (confounding) x hospital multiplier (clustering); both mean 1
  M = exp(s$sgn * log(s$HRy) * ps_cv(e, s$p0, s$p.cv) - s$c_conf) * M_h
  
  if (s$family == "poisson") {
    Y0 = rpois(n, lambda = s$mu0 * M)
    Y1 = rpois(n, lambda = s$mu1 * M)
  } else {
    Y0 = rnbinom(n, size = 1 / s$kw0, mu = s$mu0 * M)
    Y1 = rnbinom(n, size = 1 / s$kw1, mu = s$mu1 * M)
  }
  
  data.frame(Y0 = Y0, Y1 = Y1)
}

#' IPTW estimate and hospital-clustered robust SE for the log rate ratio
#'
#' Weighted Poisson or negative binomial regression of Y on X with IPT weights.
#' SE is the sandwich estimator with scores summed within hospital
#' (sandwich::vcovCL, HC0).
#'
#' @param d          Data frame from combine_outcomes() with Y (count), X (0/1),
#'                   W (IPT weight) and hosp (hospital ID).
#' @param fit_family "poisson" or "nb": the model used for the analysis.
#' @returns c(est = , se = ): estimated log rate ratio (X = 1 vs X = 0) and its
#'          hospital-clustered robust SE. NA for both if the fit fails.
est_itt_ate_count_cl = function(d, fit_family = c("poisson", "nb")) {
  fit_family = match.arg(fit_family)
  
  fit = tryCatch(
    if (fit_family == "poisson") {
      suppressWarnings(glm(Y ~ X, family = poisson, data = d, weights = W))
    } else {
      suppressWarnings(MASS::glm.nb(Y ~ X, data = d, weights = W))
    },
    error = function(err) NULL)
  if (is.null(fit)) return(c(est = NA, se = NA))
  
  est = coef(fit)["X"]
  se  = sqrt(sandwich::vcovCL(fit, cluster = d$hosp, type = "HC0")["X", "X"])
  
  c(est = unname(est), se = unname(se))
}

#' Simulate nsim datasets and compute power for the IPTW rate ratio
#'
#' Per dataset: cohort (PS, treatment, weights) -> hospitals -> potential
#' counts -> observed count -> IPTW log RR with hospital-clustered robust SE
#' -> decision via decide(). Random draws follow the order of the
#' original scripts: e, X, hosp, u, Y0, Y1.
#'
#' @param n            Sample size.
#' @param s            Setting list from setup_count_cl().
#' @param H            Number of hospitals.
#' @param test         "noninferiority", "superiority2", "superiority",
#'                     "inferiority" or "equivalence" (see decide()).
#' @param alpha        One-sided level per bound (two-sided for "superiority2").
#' @param margin       Margin on the RATIO scale (e.g. 1.1); converted to log.
#'                     Ignored for superiority tests.
#' @param lower_better TRUE if smaller counts are better (e.g. fewer days).
#' @param fit_family   Analysis model, "poisson" (fast) or "nb" (same est/SE).
#' @param nsim         Number of simulated datasets.
#' @returns One-row data.frame:
#'   RR_true   mu1 / mu0
#'   mean_RR   exp(mean log RR estimate): close to RR_true = no bias
#'   crude_RR  exp(mean unweighted log ratio): shows confounding bias
#'   emp_sd    SD of log RR estimates (the true SE)
#'   mean_se   average robust SE (close to emp_sd = SE is accurate)
#'   power     share of datasets declaring the test's alternative
#'   n_fail    number of failed fits (dropped)
power_itt_ate_count_cl = function(n, s, H, test = "noninferiority", alpha = 0.025,
                           margin = 1, lower_better = TRUE,
                           fit_family = c("poisson", "nb"), nsim = 1000) {
  fit_family = match.arg(fit_family)
  
  res = t(replicate(nsim, {
    co = simulate_cohort(n, s$alpha, s$beta)                 # e, X, W
    cl = sim_cluster(n, H = H, tau2 = s$tau2)                # hosp, M_h
    po = sim_po_count_cl(co$e, cl$M_h, s)                       # Y0, Y1
    d  = combine_outcomes(co, po)                            # observed Y
    d$hosp = cl$hosp
    c(est_itt_ate_count_cl(d, fit_family),
      crude = log(mean(d$Y[d$X == 1]) / mean(d$Y[d$X == 0])))
  }))
  
  ok = !is.na(res[, "est"])
  data.frame(RR_true  = s$RR_true,
             mean_RR  = exp(mean(res[ok, "est"])),
             crude_RR = exp(mean(res[ok, "crude"])),
             emp_sd   = sd(res[ok, "est"]),
             mean_se  = mean(res[ok, "se"]),
             power    = mean(decide(res[ok, "est"], res[ok, "se"], test = test,
                                    alpha = alpha, margin = log(margin),
                                    lower_better = lower_better)),
             n_fail   = sum(!ok))
}

## ==== Examples (run by hand) ====
if (FALSE) {
  
  ## ---- Setting from the NB script (Aim 3: discontinue vs continue antibiotics)
  set.seed(2025)
  s_null = setup_count_cl(p0 = 0.193, p.cv = 0.5, family = "nb",
                       mu0 = 5.95, mu1 = 3.71,
                       sd0 = 3.01, sd1 = 1.97,
                       HRy = 1.3, sgn = -1, icc = 0.01)
  set.seed(1)
  power_itt_ate_count_cl(n = 3908, s = s_null, H = 69, test = "noninferiority",
                  alpha = 0.025, margin = 1.1, nsim = 1000)
  
  ## ---- Setting from the Poisson script
  set.seed(2025)
  s_po = setup_count_cl(p0 = 0.161, p.cv = 0.5, family = "poisson",
                     mu0 = 1.84, mu1 = 1.65, HRy = 1.3, sgn = -1, icc = 0.01)
  
  set.seed(123)
  power_itt_ate_count_cl(n = 2696, s = s_po, H = 69, test = "noninferiority",
                  alpha = 0.025, margin = 1.1, nsim = 100)
  #   RR_true 0.897   mean_RR ~0.90   compare power with the Poisson script
}