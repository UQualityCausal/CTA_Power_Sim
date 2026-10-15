## hte_cts_change.R
## Power for the HTE on change in continuous variable (e.g., blood biomarker),
## baseline to one follow-up
source("power_core.R")

# copy of sd_change() in ate_cts_change.R — keep identical
#' SD of the change score (follow-up - baseline)
#'
#' @param sd_b SD of the biomarker at baseline (original units).
#' @param sd_f SD of the biomarker at follow-up, untreated (original units).
#' @param rho  Correlation between a person's baseline and follow-up values.
#' @returns SD of the change, original units.
#' @examples
#' sd_change(1.25, 1.25, 0.5)   # 1.25
sd_change = function(sd_b, sd_f, rho){
  # input checks
  if (sd_b <= 0 || sd_f <= 0){
    stop("sd must be positive")
  }
  if (abs(rho) > 1){
    stop("rho must be between -1 and 1")
  }
  
  sqrt(sd_b^2 + sd_f^2 - 2 * rho * sd_b * sd_f)
}

#' Simulate potential outcomes for the change score, with effect modification
#'
#' @param e         Propensity scores, from simulate_cohort()$e.
#' @param EF        Effect modifier, 0/1 per person, from sim_ef().
#' @param delta     Treatment effect in the EF = 0 group (original units).
#' @param delta_mod Extra treatment effect in the EF = 1 group: their effect is
#'                  delta + delta_mod (original units). This is what HTE tests.
#' @param b_ef      Extra change without treatment for EF = 1 (EF's own effect on
#'                  the outcome; original units). Does not affect power. Default 0.
#' @param sd_chg    SD of the change (original units), from sd_change().
#' @param conf      Confounding: extra untreated change for a 1-SD higher PS
#'                  (original units). 0 means none.
#' @param p0        Mean propensity score (same as in get_beta_params()).
#' @param p.cv      Propensity score CV (same as in get_beta_params()).
#' @returns data.frame(chg0, chg1): change without and with treatment.
sim_po_chg_hte = function(e, EF, delta, delta_mod, b_ef = 0, sd_chg, conf = 0, p0, p.cv){
  if (length(EF) != length(e)) stop("EF and e must have the same length")
  
  n = length(e)
  z = ps_z(e, p0, p.cv) # Standardized PS: mean 0, SD 1
  
  chg0 = rnorm(n, b_ef * EF + conf * z, sd_chg)     # untreated: main effect + confounding
  chg1 = chg0 + delta + delta_mod * EF              # treated: effect depends on EF
  
  data.frame(chg0 = chg0, chg1 = chg1)
}

#' IPTW estimate and robust SE for the difference in treatment effects (HTE)
#'
#' Fits the weighted regression chg ~ X + EF + X:EF with IPT weights and returns
#' the X:EF interaction: the treatment effect in EF = 1 minus the treatment
#' effect in EF = 0, with its robust (sandwich) SE.
#'
#' @param d Data frame from combine_outcomes(), with chg (observed change,
#'          original units), X (0/1 treatment), EF (0/1 effect modifier) and
#'          W (IPT weight).
#' @returns c(est = , se = ): estimated difference in treatment effects
#'          (EF = 1 minus EF = 0) and its robust SE, original units.
est_iptw_hte = function(d) {
  des = survey::svydesign(ids = ~1, weights = ~W, data = d)
  fit = survey::svyglm(chg ~ X * EF, design = des) # X + EF + X:EF
  
  cf  = summary(fit)$coefficients
  est = cf["X:EF", "Estimate"]
  se  = cf["X:EF", "Std. Error"]
  
  c(est = unname(est), se = unname(se))
}

#' Simulate nsim datasets and compute power for HTE (X:EF interaction)
#' 
#' @param n         Sample size.
#' @param pi_ef.    Proportion of people with EF = 1, strictly between 0 and 1.
#' @param p0        Expected proportion treated (allocation; 1:2 -> 1/3).
#' @param p.cv      Propensity score CV.
#' @param sd_chg    SD of the change (original units), from sd_change().
#' @param conf      Confounding: extra change for a 1-SD higher PS, untreated
#'                  (original units).
#' @param delta     True treatment effect on the change (original units).
#' @param delta_mod Extra treatment effect in the EF = 1 group (original units):
#'                  the difference in effects being tested.
#' @param b_ef      Extra untreated change for EF = 1 (original units).
#'                  Does not affect power. Default 0.
#' @param nsim      Number of simulated datasets.
#' @param alpha_adj Two-sided significance level after multiplicity adjustment.
#' @returns One-row data.frame:
#'   mean_est  average X:EF estimate (close to delta_mod = no bias)
#'   emp_sd    SD of the estimates across datasets (the true SE)
#'   mean_se   average robust SE (close to emp_sd = SE is accurate)
#'   power     share of datasets that reject H0: no difference in effects
#' See the examples at the bottom of this file.
sim_power_hte = function(n, pi_ef, p0, p.cv, sd_chg, conf, delta, delta_mod,
                         b_ef = 0, nsim = 1000, alpha_adj = 0.05) {
  bp   = get_beta_params(p0, p.cv)      # once per setting
  crit = qnorm(1 - alpha_adj / 2)  # critical value from alpha_adj
  
  res = t(replicate(nsim, {
    co = simulate_cohort(n, bp["alpha"], bp["beta"]) #cohort
    co$EF = sim_ef(n, pi_ef)
    po = sim_po_chg_hte(e = co$e, EF = co$EF, delta = delta, 
                        delta_mod = delta_mod, b_ef = b_ef, 
                        sd_chg = sd_chg, conf = conf, 
                        p0 = p0, p.cv = p.cv) # potential outcomes
    d  = combine_outcomes(co, po) # observed data
    est_iptw_hte(d) # c(est, se)
  }))
  
  data.frame(mean_est = mean(res[, "est"]),
             emp_sd   = sd(res[, "est"]), # the spread of the estimates across datasets
             mean_se  = mean(res[, "se"]), # the average of the SEs
             power    = mean(abs(res[, "est"] / res[, "se"]) > crit))
}

# Examples
## 1. No effect modification: power should be about alpha_adj
set.seed(1)
sim_power_hte(n = 1000, pi_ef = 0.2, p0 = 1/3, p.cv = 0.25, sd_chg = 1.25,
              conf = 0.1, delta = 0.2, delta_mod = 0,
              nsim = 1000, alpha_adj = 0.05)
#   mean_est ~0   power ~0.05

## 2. Basic call: difference in effects of 0.5 units, pi = 0.2
set.seed(2)
sim_power_hte(1000, pi_ef = 0.2, p0 = 1/3, p.cv = 0.25, sd_chg = 1.25,
              conf = 0.1, delta = 0.2, delta_mod = 0.5,
              nsim = 200, alpha_adj = 0.05)
#   mean_est ~0.5   emp_sd ~ mean_se ~0.21   power ~0.66

## 3. Effect of pi: a more balanced modifier gives more power
set.seed(3)
sim_power_hte(1000, pi_ef = 0.5, p0 = 1/3, p.cv = 0.25, sd_chg = 1.25,
              conf = 0.1, delta = 0.2, delta_mod = 0.5,
              nsim = 200, alpha_adj = 0.05)
#   mean_se ~0.13   power ~0.84

## 4. delta and b_ef don't change HTE power (compare with example 2)
set.seed(2)
sim_power_hte(1000, pi_ef = 0.2, p0 = 1/3, p.cv = 0.25, sd_chg = 1.25,
              conf = 0.1, delta = 1, delta_mod = 0.5, b_ef = 0.5,
              nsim = 200, alpha_adj = 0.05)
#   mean_se and power about the same as example 2