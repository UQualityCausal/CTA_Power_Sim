## itt_hte_count_cl.R: power for ITT HTE on a count outcome clustered by
## hospital: ratio of rate ratios between levels of a binary effect modifier
## (X:EF interaction). IPTW Poisson or NB regression with a hospital-clustered
## robust SE; tests via decide(). mu0, mu1 describe the EF = 0 group.
source("itt_ate_count_cl.R")

#' Simulate potential count outcomes with a binary effect modifier
#'
#' Same as sim_po_count_cl(), plus EF. Each person's mean is
#'   X = 0:  mu0 x ef_ratio^EF          x M
#'   X = 1:  mu1 x ef_ratio^EF x MOR^EF x M
#' so the rate ratio is RR0 = mu1/mu0 when EF = 0 and RR0 x MOR when EF = 1.
#' M = PS multiplier x hospital multiplier (both mean 1), as in sim_po_count_cl().
#'
#' @param e        Propensity scores, from simulate_cohort()$e.
#' @param M_h      Hospital multiplier per person, from sim_cluster()$M_h.
#' @param EF       0/1 effect modifier per person.
#' @param s        Setting list from setup_count_cl() (EF = 0 group means).
#' @param MOR      Ratio of rate ratios, EF = 1 vs EF = 0. 1 = no HTE.
#' @param ef_ratio Main effect of EF: control-arm mean in EF = 1 divided by
#'                 control-arm mean in EF = 0. 1 = EF does not shift baseline.
#' @returns data.frame(Y0, Y1): count if X = 0 and if X = 1.
sim_po_count_hte_cl = function(e, M_h, EF, s, MOR = 1, ef_ratio = 1) {
  n = length(e)
  if (length(M_h) != n || length(EF) != n)
    stop("e, M_h and EF must have the same length")

  # PS multiplier x hospital multiplier x EF main effect
  M  = exp(s$sgn * log(s$HRy) * ps_cv(e, s$p0, s$p.cv) - s$c_conf) * M_h *
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
est_itt_hte_count_cl = function(d, fit_family = c("poisson", "nb")) {
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
#' @param s            Setting list from setup_count_cl() (EF = 0 group).
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
power_itt_hte_count_cl = function(n, s, H, pi_ef, MOR, ef_ratio = 1,
                               test = "superiority2", alpha = 0.05,
                               margin = 1, lower_better = TRUE,
                               fit_family = c("poisson", "nb"), nsim = 1000) {
  fit_family = match.arg(fit_family)

  res = t(replicate(nsim, {
    co = simulate_cohort(n, s$alpha, s$beta)                  # e, X, W
    cl = sim_cluster(n, H = H, tau2 = s$tau2)                 # hosp, M_h
    EF = sim_ef(n, pi_ef)                                     # effect modifier (core)
    po = sim_po_count_hte_cl(co$e, cl$M_h, EF, s, MOR, ef_ratio) # Y0, Y1
    d  = combine_outcomes(co, po)                             # observed Y
    d$hosp = cl$hosp
    d$EF   = EF

    # crude (unweighted) log MOR from the four cell means
    cm = tapply(d$Y, list(d$X, d$EF), mean)                   # rows X, cols EF
    crude = log((cm["1", "1"] / cm["0", "1"]) / (cm["1", "0"] / cm["0", "0"]))

    c(est_itt_hte_count_cl(d, fit_family), crude = crude)
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
  s_nb = setup_count_cl(p0 = 0.193, p.cv = 0.5, family = "nb",
                     mu0 = 5.95, mu1 = 3.71, sd0 = 3.01, sd1 = 1.97,
                     HRy = 1.3, sgn = -1, icc = 0.01)
  set.seed(1)
  power_itt_hte_count_cl(n = 3908, s = s_nb, H = 69, pi_ef = 0.3, MOR = 1,
                      nsim = 1000)

  ## ---- 2. Power across a grid of MORs (NB setting from itt_ate_count_cl.R)
  MORs = c(0.80, 0.85, 0.90)
  pow = do.call(rbind, lapply(MORs, function(m) {
    set.seed(1)
    power_itt_hte_count_cl(n = 3908, s = s_nb, H = 69, pi_ef = 0.3, MOR = m,
                        nsim = 1000)
  }))
  pow
}
