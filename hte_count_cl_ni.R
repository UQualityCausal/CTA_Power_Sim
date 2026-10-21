## =============================================================================
## hte_count_cl_ni.R -- Power for effect modification (HTE) on a count outcome
##
## Assumptions:
##   - EF is independent of the propensity score (it is not a confounder), so
##     the PS model does not need EF.
##   - mu0, mu1 (and sd0, sd1 for NB) given to count_setup() describe the
##     EF = 0 group. The EF = 1 group's means are scaled by ef_ratio (both
##     arms) and MOR (treated arm only).
##   - NB: the within-person dispersion (kw0, kw1) is the same in both EF
##     groups, so SD in the EF = 1 group follows from its mean.
## =============================================================================

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
#' @param mu0, mu1 Population mean of Y0/Y1.
#' @param sd0, sd1 SD of the outcome by arm (NB only; need sd^2 > mean).
#' @param HRy      Size of PS effect on the outcome: ratio per 1-CV of PS.
#' @param sgn      Direction: -1 = higher PS -> smaller outcome; +1 = larger.
#' @param icc      Hospital ICC on the count scale, calibrated on the untreated
#'                 arm (0 = no clustering).
#' @param ndraw    Monte Carlo draws used for c_conf and EM2.
#' @returns List with the inputs (p0, p.cv, family, mu0, mu1, HRy, sgn, icc)
#'          and alpha, beta, c_conf, EM2, tau2, kw0, kw1 (NA for Poisson),
#'          RR_true = mu1 / mu0.
count_setup = function(p0, p.cv, family = c("poisson", "nb"), mu0, mu1,
                       sd0 = NA, sd1 = NA, HRy = 1, sgn = 1, icc = 0,
                       ndraw = 1e6) {
  family = match.arg(family)
  
  # 1. Beta parameters for the PS
  bp = get_beta_params(p0, p.cv)
  if (bp["alpha"] <= 1 || bp["beta"] <= 1)
    stop("IPW weights have infinite mean: reduce p.cv")
  
  # 2. PS multiplier constants: c_conf so E[M_ps] = 1; EM2 = E[M_ps^2]
  e_big  = rbeta(ndraw, bp["alpha"], bp["beta"])
  lp_big = sgn * log(HRy) * ps_cv(e, p0, p.cv)
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
#' @param tau2 Variance of the log hospital effect, e.g. count_setup()$tau2.
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

#' Simulate potential count outcomes with a binary effect modifier
#'
#' Same as sim_po_count(), plus EF. Each person's mean is
#'   X = 0:  mu0 x ef_ratio^EF          x M
#'   X = 1:  mu1 x ef_ratio^EF x MOR^EF x M
#' so the rate ratio is RR0 = mu1/mu0 when EF = 0 and RR0 x MOR when EF = 1.
#' M = PS multiplier x hospital multiplier (both mean 1), as in sim_po_count().
#'
#' @param e        Propensity scores, from simulate_cohort()$e.
#' @param M_h      Hospital multiplier per person, from sim_cluster()$M_h.
#' @param EF       0/1 effect modifier per person.
#' @param s        Setting list from count_setup() (EF = 0 group means).
#' @param MOR      Ratio of rate ratios, EF = 1 vs EF = 0. 1 = no HTE.
#' @param ef_ratio Main effect of EF: control-arm mean in EF = 1 divided by
#'                 control-arm mean in EF = 0. 1 = EF does not shift baseline.
#' @returns data.frame(Y0, Y1): count if X = 0 and if X = 1.
sim_po_count_hte = function(e, M_h, EF, s, MOR = 1, ef_ratio = 1) {
  n = length(e)
  if (length(M_h) != n || length(EF) != n)
    stop("e, M_h and EF must have the same length")

  # PS multiplier x hospital multiplier x EF main effect
  M  = exp(s$sgn * log(s$HRy) * ps_cv(e_big, s$p0, s$p.cv) - s$c_conf) * M_h *
       ef_ratio^EF
  m0 = s$mu0 * M                    # mean if untreated
  m1 = s$mu1 * MOR^EF * M           # mean if treated (HTE enters here only)

  if (s$family == "poisson") {
    Y0 = rpois(n, lambda = m0)
    Y1 = rpois(n, lambda = m1)
  } else {
    Y0 = rnbinom(n, size = 1 / s$kw0, mu = m0)
    Y1 = rnbinom(n, size = 1 / s$kw1, mu = m1)
  }

  data.frame(Y0 = Y0, Y1 = Y1)
}


#' IPTW estimate and hospital-clustered robust SE for the log MOR
#'
#' Weighted Poisson or NB regression Y ~ X + EF + X:EF with IPT weights.
#' SE from sandwich::vcovCL (HC0), scores summed within hospital.
#'
#' @param d          Data frame with Y, X (0/1 numeric), EF (0/1), W, hosp.
#' @param fit_family "poisson" or "nb": the analysis model.
#' @returns c(est, se, est_x): log MOR (X:EF coefficient) and its robust SE,
#'          plus the log RR in the EF = 0 group (X coefficient). NA if the
#'          fit fails.
est_iptw_count_hte = function(d, fit_family = c("poisson", "nb")) {
  fit_family = match.arg(fit_family)

  fit = tryCatch(
    if (fit_family == "poisson") {
      suppressWarnings(glm(Y ~ X * EF, family = poisson, data = d, weights = W))
    } else {
      suppressWarnings(MASS::glm.nb(Y ~ X * EF, data = d, weights = W))
    },
    error = function(err) NULL)
  if (is.null(fit)) return(c(est = NA, se = NA, est_x = NA))

  b = coef(fit)
  V = sandwich::vcovCL(fit, cluster = d$hosp, type = "HC0")

  c(est   = unname(b["X:EF"]),
    se    = sqrt(V["X:EF", "X:EF"]),
    est_x = unname(b["X"]))
}


#' Simulate nsim datasets and compute power for the IPTW ratio of rate ratios
#'
#' Per dataset: cohort (e, X, W) -> hospitals -> EF -> potential counts ->
#' observed count -> IPTW log MOR with hospital-clustered SE -> decide().
#' Random draws follow the order: e, X, hosp, u, EF, Y0, Y1.
#'
#' @param n            Sample size.
#' @param s            Setting list from count_setup() (EF = 0 group).
#' @param H            Number of hospitals.
#' @param pi_ef        Pr(EF = 1), e.g. 0.2 for a 20% subgroup.
#' @param MOR          True ratio of rate ratios (EF = 1 vs EF = 0).
#' @param ef_ratio     Main effect of EF on the control-arm mean (1 = none).
#' @param test         Passed to decide(); HTE is usually "superiority2".
#' @param alpha        Level passed to decide() (two-sided for
#'                     "superiority2"). For several modifiers, use
#'                     0.05 / number of modifiers.
#' @param margin       Margin on the ratio scale; ignored for superiority.
#' @param lower_better TRUE if smaller counts are better.
#' @param fit_family   Analysis model, "poisson" or "nb".
#' @param nsim         Number of simulated datasets.
#' @returns One-row data.frame:
#'   MOR_true   input MOR
#'   mean_MOR   exp(mean log MOR estimate): close to MOR_true = no bias
#'   crude_MOR  exp(mean unweighted ratio of cell-mean RRs)
#'   RR_ef0     exp(mean log RR, EF = 0): should be close to mu1 / mu0
#'   RR_ef1     exp(mean log RR, EF = 1): should be close to mu1 / mu0 x MOR
#'   emp_sd     SD of log MOR estimates (the true SE)
#'   mean_se    average robust SE (close to emp_sd = SE is accurate)
#'   power      share of datasets declaring the test's alternative
#'   n_fail     number of failed fits (dropped)
sim_power_count_hte = function(n, s, H, pi_ef, MOR, ef_ratio = 1,
                               test = "superiority2", alpha = 0.05,
                               margin = 1, lower_better = TRUE,
                               fit_family = c("poisson", "nb"), nsim = 1000) {
  fit_family = match.arg(fit_family)

  res = t(replicate(nsim, {
    co = simulate_cohort(n, s$alpha, s$beta)                  # e, X, W
    cl = sim_cluster(n, H = H, tau2 = s$tau2)                 # hosp, M_h
    EF = sim_ef(n, pi_ef)                                     # effect modifier (core)
    po = sim_po_count_hte(co$e, cl$M_h, EF, s, MOR, ef_ratio) # Y0, Y1
    d  = combine_outcomes(co, po)                             # observed Y
    d$hosp = cl$hosp
    d$EF   = EF

    # crude (unweighted) log MOR from the four cell means
    cm = tapply(d$Y, list(d$X, d$EF), mean)                   # rows X, cols EF
    crude = log((cm["1", "1"] / cm["0", "1"]) / (cm["1", "0"] / cm["0", "0"]))

    c(est_iptw_count_hte(d, fit_family), crude = crude)
  }))

  ok = !is.na(res[, "est"])
  data.frame(MOR_true  = MOR,
             mean_MOR  = exp(mean(res[ok, "est"])),
             crude_MOR = exp(mean(res[ok, "crude"])),
             RR_ef0    = exp(mean(res[ok, "est_x"])),
             RR_ef1    = exp(mean(res[ok, "est_x"] + res[ok, "est"])),
             emp_sd    = sd(res[ok, "est"]),
             mean_se   = mean(res[ok, "se"]),
             power     = mean(decide(res[ok, "est"], res[ok, "se"],
                                     test = test, alpha = alpha,
                                     margin = log(margin),
                                     lower_better = lower_better)),
             n_fail    = sum(!ok))
}

## ==== Examples (run by hand) ====
if (FALSE) {

  ## ---- 1. Type I error: MOR = 1 -> power should be close to alpha
  set.seed(2025)
  s_nb = count_setup(p0 = 0.193, p.cv = 0.5, family = "nb",
                     mu0 = 5.95, mu1 = 3.71, sd0 = 3.01, sd1 = 1.97,
                     HRy = 1.3, sgn = -1, icc = 0.01)
  set.seed(1)
  sim_power_count_hte(n = 3908, s = s_nb, H = 69, pi_ef = 0.3, MOR = 1,
                      nsim = 1000)

  ## ---- 2. Power across a grid of MORs (NB setting from ate_count.R)
  MORs = c(0.80, 0.85, 0.90)
  pow = do.call(rbind, lapply(MORs, function(m) {
    set.seed(1)
    sim_power_count_hte(n = 3908, s = s_nb, H = 69, pi_ef = 0.3, MOR = m,
                        nsim = 1000)
  }))
  pow
}
