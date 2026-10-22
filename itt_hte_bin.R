## Power for HTE on a binary outcome: ratio of treatment ORs between levels of
## a binary effect modifier (X:EF interaction), IPTW logistic regression.
source("itt_ate_bin.R")

#' Simulate potential binary outcomes with effect modification
#'
#' Logistic model, conditional on the PS:
#'   untreated: logit P(Y0 = 1) = b0 + b1 * ps_cv(e, p0, p.cv) + log(OR_ef) * EF
#'   treated:   logit P(Y1 = 1) = same + log(TOR) + log(MOR) * EF
#' so the conditional treatment OR is TOR when EF = 0 and TOR * MOR when EF = 1.
#'
#' @param e     Propensity scores, from simulate_cohort()$e.
#' @param EF    Effect modifier, 0/1 per person, from sim_ef().
#' @param py0   Risk if untreated at the average PS, EF = 0 group.
#' @param HRy   Confounding: OR per 1-CV increase in the PS (see ps_cv()).
#'              1 = no confounding.
#' @param TOR   Conditional treatment OR in the EF = 0 group.
#' @param MOR   Ratio of conditional treatment ORs, EF = 1 vs EF = 0 (tested).
#' @param OR_ef EF's own effect on the odds, untreated. 1 = none. Changes power.
#' @param p0    Mean PS (same as in get_beta_params()).
#' @param p.cv  PS CV (same as in get_beta_params()).
#' @returns data.frame(Y0, Y1): 0/1 outcome if untreated and if treated.
sim_po_bin_hte = function(e, EF, py0, HRy = 1, TOR, MOR, OR_ef = 1, p0, p.cv) {
  if (length(EF) != length(e)) stop("EF and e must have the same length")
  if (py0 <= 0 || py0 >= 1) stop("py0 should be strictly between 0 and 1")
  
  n  = length(e)
  b0 = qlogis(py0)
  b1 = log(HRy)
  
  lp0 = b0 + b1 * ps_cv(e, p0, p.cv) + log(OR_ef) * EF   # untreated log-odds
  lp1 = lp0 + log(TOR) + log(MOR) * EF                   # treated: effect depends on EF
  
  Y0 = rbinom(n, 1, plogis(lp0))
  Y1 = rbinom(n, 1, plogis(lp1))
  
  data.frame(Y0 = Y0, Y1 = Y1)
}

#' IPTW estimate and robust SE for the ratio of ORs (HTE)
#'
#' Weighted logistic regression Y ~ X + EF + X:EF with IPT weights and a
#' robust (HC0) SE.
#'
#' @param d Data frame from combine_outcomes(), with Y (0/1), X (0/1 treatment),
#'          EF (0/1 effect modifier) and W (IPT weight).
#' @returns c(est = , se = , est_x = ): log ratio of marginal ORs (X:EF) and
#'          its robust SE, plus the log OR in the EF = 0 group (X).
est_itt_hte_bin = function(d) {
  fit = suppressWarnings(glm(Y ~ X * EF, family = binomial, data = d, weights = W))
  V   = sandwich::vcovHC(fit, type = "HC0")
  
  c(est   = unname(coef(fit)["X:EF"]),
    se    = sqrt(V["X:EF", "X:EF"]),
    est_x = unname(coef(fit)["X"]))
}

#' Simulate nsim datasets and compute empirical power for the ratio of ORs
#'
#' @param n         Sample size.
#' @param pi_ef     Proportion of people with EF = 1, strictly between 0 and 1.
#' @param p0        Expected proportion treated (allocation; 1:3 -> 1/4).
#' @param p.cv      PS CV.
#' @param py0       Risk if untreated at the average PS, EF = 0 group.
#' @param HRy       Confounding OR per 1-CV increase in the PS. 1 = none.
#' @param TOR       Conditional treatment OR in the EF = 0 group.
#' @param MOR       Ratio of conditional treatment ORs, EF = 1 vs EF = 0 (tested).
#' @param OR_ef     EF's own effect on the odds, untreated. 1 = none.
#' @param nsim      Number of simulated datasets.
#' @param alpha_adj Two-sided significance level after multiplicity adjustment.
#' @returns One-row data.frame:
#'   MOR            input (conditional) ratio of ORs
#'   mean_est       average log ratio-of-ORs estimate
#'   marg_MOR       exp(mean_est): the marginal ratio of ORs = the detectable
#'                  effect to report (= MOR when HRy = 1)
#'   OR_ef0         exp(mean log OR, EF = 0): close to TOR
#'   OR_ef1         exp(mean log OR, EF = 1): close to TOR * MOR
#'   emp_sd         SD of the estimates (the true SE)
#'   mean_se        average robust SE (close to emp_sd = SE is accurate)
#'   mean_events11  average Y = 1 among treated with EF = 1 (smallest cell)
#'   power          share of datasets that reject H0: MOR = 1
power_itt_hte_bin = function(n, pi_ef, p0, p.cv, py0, HRy = 1, TOR, MOR,
                             OR_ef = 1, nsim = 1000, alpha_adj = 0.05) {
  bp   = get_beta_params(p0, p.cv)
  crit = qnorm(1 - alpha_adj / 2)
  
  res = t(replicate(nsim, {
    co    = simulate_cohort(n, bp["alpha"], bp["beta"], estimand = "ATE")
    co$EF = sim_ef(n, pi_ef)
    po    = sim_po_bin_hte(e = co$e, EF = co$EF, py0 = py0, HRy = HRy,
                           TOR = TOR, MOR = MOR, OR_ef = OR_ef,
                           p0 = p0, p.cv = p.cv)
    d     = combine_outcomes(co, po)
    c(est_itt_hte_bin(d),
      events11 = sum(d$Y[d$X == 1 & d$EF == 1]))   # events among X = 1 & EF = 1
  }))
  
  data.frame(MOR           = MOR,
             mean_est      = mean(res[, "est"]),
             marg_MOR      = exp(mean(res[, "est"])),
             OR_ef0        = exp(mean(res[, "est_x"])),
             OR_ef1        = exp(mean(res[, "est_x"] + res[, "est"])),
             emp_sd        = sd(res[, "est"]),
             mean_se       = mean(res[, "se"]),
             mean_events11 = mean(res[, "events11"]),
             power         = mean(abs(res[, "est"] / res[, "se"]) > crit))
}

# Examples
if (FALSE){
  ## Shared settings
  n   = 4000
  p0  = 1/4
  py0 = 0.3
  
  ## 1. No effect modification: power should be about alpha_adj
  set.seed(1)
  power_itt_hte_bin(n, pi_ef = 0.3, p0, p.cv = 0.25, py0, HRy = 2,
                    TOR = 0.8, MOR = 1, nsim = 1000)
  #   marg_MOR ~1   OR_ef0 ~ OR_ef1 ~0.80   mean_events11 ~80   power ~0.05
  
  ## 2. Basic call: MOR = 0.6, pi = 0.3
  set.seed(2)
  power_itt_hte_bin(n, pi_ef = 0.3, p0, p.cv = 0.25, py0, HRy = 2,
                    TOR = 0.8, MOR = 0.6, nsim = 500)
  #   marg_MOR ~0.60   OR_ef0 ~0.80   OR_ef1 ~0.48 (= 0.8 x 0.6)
  #   emp_sd ~ mean_se ~0.20   mean_events11 ~53   power ~0.76
  
  ## 3. Effect of pi: a balanced modifier gives more power
  set.seed(3)
  power_itt_hte_bin(n, pi_ef = 0.5, p0, p.cv = 0.25, py0, HRy = 2,
                    TOR = 0.8, MOR = 0.6, nsim = 500)
  #   mean_se ~0.18   mean_events11 ~89   power ~0.82
  
  ## 4. OR_ef changes power: more events in EF = 1 (compare with example 2)
  set.seed(2)
  power_itt_hte_bin(n, pi_ef = 0.3, p0, p.cv = 0.25, py0, HRy = 2,
                    TOR = 0.8, MOR = 0.6, OR_ef = 2, nsim = 500)
  #   mean_se ~0.18   mean_events11 ~91   power ~0.83
  
  ## 5. With vs without confounding (same everything else)
  set.seed(5)
  power_itt_hte_bin(n, pi_ef = 0.3, p0, p.cv = 0.25, py0, HRy = 1,
                    TOR = 0.8, MOR = 0.6, nsim = 500)
  #   marg_MOR ~0.60   mean_se ~0.20   power ~0.72
  set.seed(5)
  power_itt_hte_bin(n, pi_ef = 0.3, p0, p.cv = 0.25, py0, HRy = 4,
                    TOR = 0.8, MOR = 0.6, nsim = 500)
  #   marg_MOR ~0.60 (IPTW removes confounding)   power ~0.73 
  
  ## 6. Power curve over MOR (to find the MOR with 80% power)
  mors  = c(0.50, 0.55, 0.60, 0.65, 0.70)
  curve = do.call(rbind, lapply(mors, function(m) {
    set.seed(6)
    power_itt_hte_bin(n, pi_ef = 0.3, p0, p.cv = 0.25, py0, HRy = 2,
                      TOR = 0.8, MOR = m, nsim = 300)
  }))
  curve[, c("MOR", "marg_MOR", "mean_events11", "power")]
  #   power ~0.93, 0.85, 0.75, 0.61, 0.44  -> 80% power at MOR ~0.57
}