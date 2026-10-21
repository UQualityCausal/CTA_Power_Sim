source("itt_ate_count.R")

#' Simulate potential count outcomes with a binary effect modifier
#'
#' Same as sim_po_count(), plus EF. Each person's mean is
#'   untreated: mu0 x ef_ratio^EF          x M
#'   treated:   mu1 x ef_ratio^EF x MOR^EF x M
#' so the rate ratio is mu1/mu0 when EF = 0 and mu1/mu0 x MOR when EF = 1.
#'
#' @param e        Propensity scores, from simulate_cohort()$e.
#' @param EF       Effect modifier, 0/1 per person, from sim_ef().
#' @param s        Setting list from setup_count() (EF = 0 group).
#' @param MOR      Ratio of rate ratios, EF = 1 vs EF = 0 (the tested quantity).
#'                 1 = no HTE.
#' @param ef_ratio EF's own effect: untreated mean in EF = 1 divided by
#'                 untreated mean in EF = 0. 1 = none. Changes power (more
#'                 events in the EF = 1 group -> smaller SE).
#' @returns data.frame(Y0, Y1): count if untreated and if treated.
sim_po_count_hte = function(e, EF, s, MOR = 1, ef_ratio = 1) {
  if (length(EF) != length(e)) stop("EF and e must have the same length")

  n  = length(e)
  M  = exp(log(s$HRy) * ps_cv(e, s$p0, s$p.cv) - s$c_conf) * ef_ratio^EF
  m0 = s$mu0 * M                  # untreated mean
  m1 = s$mu1 * MOR^EF * M         # treated mean: effect depends on EF

  if (s$family == "poisson") {
    Y0 = rpois(n, lambda = m0)
    Y1 = rpois(n, lambda = m1)
  } else {
    Y0 = rnbinom(n, size = 1 / s$kw0, mu = m0)
    Y1 = rnbinom(n, size = 1 / s$kw1, mu = m1)
  }

  data.frame(Y0 = Y0, Y1 = Y1)
}

#' IPTW estimate and robust SE for the ratio of rate ratios (HTE)
#'
#' Weighted Poisson regression Y ~ X + EF + X:EF with IPT weights and a robust
#' (sandwich, HC0) SE. Valid for NB data too (saturated model; see
#' est_itt_ate_count()).
#'
#' @param d Data frame from combine_outcomes(), with Y (count), X (0/1
#'          treatment), EF (0/1 effect modifier) and W (IPT weight).
#' @returns c(est = , se = , est_x = ): estimated log MOR (X:EF) and its
#'          robust SE, plus the log RR in the EF = 0 group (X).
est_itt_hte_count = function(d) {
  fit = suppressWarnings(glm(Y ~ X * EF, family = poisson, data = d, weights = W))
  V   = sandwich::vcovHC(fit, type = "HC0")

  c(est   = unname(coef(fit)["X:EF"]),
    se    = sqrt(V["X:EF", "X:EF"]),
    est_x = unname(coef(fit)["X"]))
}

#' Simulate nsim datasets and compute power for the ratio of rate ratios
#'
#' @param n         Sample size.
#' @param pi_ef     Proportion of people with EF = 1, strictly between 0 and 1.
#' @param p0        Expected proportion treated (allocation; 1:3 -> 1/4).
#' @param p.cv      Propensity score CV.
#' @param family    "poisson" or "nb": how counts are generated.
#' @param mu0, mu1  Population means of Y0 / Y1 in the EF = 0 group.
#' @param sd0, sd1  Population SDs of Y0 / Y1 in the EF = 0 group (NB only).
#' @param HRy       Confounding: rate ratio per 1-CV increase in the PS.
#'                  1 = no confounding.
#' @param MOR       Ratio of rate ratios, EF = 1 vs EF = 0 (tested).
#' @param ef_ratio  EF's own effect on the untreated mean. 1 = none.
#' @param nsim      Number of simulated datasets.
#' @param alpha_adj Two-sided significance level after multiplicity adjustment.
#' @returns One-row data.frame:
#'   MOR_true  input MOR
#'   mean_est  average log MOR estimate
#'   marg_MOR  exp(mean_est): close to MOR_true = no bias
#'   RR_ef0    exp(mean log RR, EF = 0): close to mu1 / mu0
#'   RR_ef1    exp(mean log RR, EF = 1): close to mu1 / mu0 x MOR
#'   emp_sd    SD of the log MOR estimates (the true SE)
#'   mean_se   average robust SE (close to emp_sd = SE is accurate)
#'   power     share of datasets that reject H0: MOR = 1
#' See the examples at the bottom of this file.
power_itt_hte_count = function(n, pi_ef, p0, p.cv, family = c("poisson", "nb"),
                               mu0, mu1, sd0 = NA, sd1 = NA, HRy = 1,
                               MOR, ef_ratio = 1, nsim = 1000, alpha_adj = 0.05) {
  family = match.arg(family)
  s    = setup_count(p0, p.cv, family, mu0, mu1, sd0, sd1, HRy)  # once per setting
  crit = qnorm(1 - alpha_adj / 2)

  res = t(replicate(nsim, {
    co    = simulate_cohort(n, s$alpha, s$beta, estimand = "ATE")
    co$EF = sim_ef(n, pi_ef)
    po    = sim_po_count_hte(co$e, co$EF, s, MOR = MOR, ef_ratio = ef_ratio)
    d     = combine_outcomes(co, po)
    est_itt_hte_count(d)
  }))

  data.frame(MOR_true = MOR,
             mean_est = mean(res[, "est"]),
             marg_MOR = exp(mean(res[, "est"])),
             RR_ef0   = exp(mean(res[, "est_x"])),
             RR_ef1   = exp(mean(res[, "est_x"] + res[, "est"])),
             emp_sd   = sd(res[, "est"]),
             mean_se  = mean(res[, "se"]),
             power    = mean(abs(res[, "est"] / res[, "se"]) > crit))
}

# Examples
if (FALSE){
## 1. No effect modification: power should be about alpha_adj
set.seed(1)
power_itt_hte_count(n = 3908, pi_ef = 0.3, p0 = 0.193, p.cv = 0.5,
                    family = "nb", mu0 = 5.95, mu1 = 3.71, sd0 = 3.01, sd1 = 1.97,
                    HRy = 1 / 1.3, MOR = 1, nsim = 1000)
#   marg_MOR ~1   RR_ef0 ~ RR_ef1 ~0.62   power ~0.05

## 2. Basic call: MOR = 0.85, pi = 0.3
set.seed(2)
power_itt_hte_count(n = 3908, pi_ef = 0.3, p0 = 0.193, p.cv = 0.5,
                    family = "nb", mu0 = 5.95, mu1 = 3.71, sd0 = 3.01, sd1 = 1.97,
                    HRy = 1 / 1.3, MOR = 0.85, nsim = 500)
#   marg_MOR ~0.85   RR_ef1 ~0.53   emp_sd ~ mean_se ~0.056   power ~0.84

## 3. ef_ratio changes power (compare with example 2)
set.seed(2)
power_itt_hte_count(n = 3908, pi_ef = 0.3, p0 = 0.193, p.cv = 0.5,
                    family = "nb", mu0 = 5.95, mu1 = 3.71, sd0 = 3.01, sd1 = 1.97,
                    HRy = 1 / 1.3, MOR = 0.85, ef_ratio = 2, nsim = 500)
#   mean_se ~0.046   power ~0.95

## 4. Power curve over MOR
mors  = c(0.80, 0.85, 0.90, 0.95)
curve = do.call(rbind, lapply(mors, function(m) {
  set.seed(4)
  power_itt_hte_count(n = 3908, pi_ef = 0.3, p0 = 0.193, p.cv = 0.5,
                      family = "nb", mu0 = 5.95, mu1 = 3.71, sd0 = 3.01,
                      sd1 = 1.97, HRy = 1 / 1.3, MOR = m, nsim = 200)
}))
curve[, c("MOR_true", "marg_MOR", "mean_se", "power")]
}
