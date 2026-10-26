## itt_ate_chg.R: power for the ITT ATE on change in a continuous variable
## (e.g. a blood biomarker), baseline to one follow-up.
## IPTW linear regression (survey::svyglm) with a design-based SE.
source("power_core.R")

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
#' Simulate potential outcomes
#' 
#' @param e         Propensity scores, one per person, from simulate_cohort()$e.
#' @param delta     True treatment effect on the change, original units.
#' @param sd_chg    The sd of the change between baseline and follow-up.
#' @param conf      Confounding: extra change, in both arms, for someone 
#'                  whose propensity score is 1 CV above average. 0 means none.
#' @param p0        Mean propensity score.
#' @param p.cv      Propensity score CV.
#' @returns         Data frame with potential outcomes for treatment and control
#' @examples
#' po = sim_po_chg(co$e, 0.5, sd_change(1.25, 1.25, 0.5), 0.01, 1/3, 0.5)
sim_po_chg = function(e, delta, sd_chg, conf = 0, p0, p.cv){
  n = length(e)
  z = ps_cv(e, p0, p.cv) # PS in CV units
  
  chg0 = rnorm(n, conf * z, sd_chg)
  chg1 = chg0 + delta
  
  data.frame(chg0 = chg0, chg1 = chg1)
}

#' IPTW estimate and robust SE for the change score
#' 
#' @param d Data frame from combine_outcomes(), with chg (observed change, 
#'          original units), X (0/1 treatment) and W (IPT weight).
#' @returns c(est = , se = ): weighted difference in mean change (treated minus
#'          control) and its robust SE, original units.
#' @examples 
#' res = est_itt_ate_chg(d)
est_itt_ate_chg=function(d){
  des = survey::svydesign(ids = ~1, weights = ~W, data =d)
  fit = survey::svyglm(chg~X, design = des)
  
  cf = summary(fit)$coefficients
  est = cf["X", "Estimate"]
  se  = cf["X", "Std. Error"]                  
  
  c(est = unname(est), se = unname(se))
  
}

#' Simulate nsim datasets and compute power
#' 
#' @param n         Sample size.
#' @param p0        Expected proportion treated (allocation; 1:2 -> 1/3).
#' @param p.cv      Propensity score CV.
#' @param sd_chg    SD of the change (original units), from sd_change().
#' @param conf      Confounding: extra change for a 1-CV higher PS, in both arms
#'                  (original units).
#' @param delta     True treatment effect on the change (original units).
#' @param nsim      Number of simulated datasets.
#' @param alpha_adj Two-sided significance level after multiplicity adjustment.
#' @returns One-row data.frame:
#'   mean_est  average estimate (close to delta = no bias)
#'   emp_sd    SD of the estimates across datasets (the true SE)
#'   mean_se   average robust SE (close to emp_sd = SE is accurate)
#'   power     share of datasets that reject H0: effect = 0
#' See the examples at the bottom of this file.
power_itt_ate_chg = function(n, p0, p.cv, sd_chg, conf, delta, nsim = 1000,
                     alpha_adj = 0.05) {
  bp   = get_beta_params(p0, p.cv)      # once per setting
  crit = qnorm(1 - alpha_adj / 2)  # critical value from alpha_adj
  
  res = t(replicate(nsim, {
    co = simulate_cohort(n, bp["alpha"], bp["beta"]) #cohort
    po = sim_po_chg(co$e, delta = delta, sd_chg = sd_chg, conf = conf, 
                    p0 = p0, p.cv = p.cv) # potential outcomes
    d  = combine_outcomes(co, po) # observed data
    est_itt_ate_chg(d) # c(est, se)
  }))
  
  data.frame(mean_est = mean(res[, "est"]),
             emp_sd   = sd(res[, "est"]), # the spread of the estimates across datasets
             mean_se  = mean(res[, "se"]), # the average of the SEs
             power    = mean(abs(res[, "est"] / res[, "se"]) > crit))
}

if (FALSE){
# Examples
## 1. Basic call: pTau217, rho = 0.5, low PS variation, small confounding
set.seed(1)
power_itt_ate_chg(n = 1000, p0 = 1/3, p.cv = 0.25,
          sd_chg = sd_change(1.25, 1.25, 0.5),     # = 1.25
          conf = 0.1, delta = 0.3,
          nsim = 200, alpha_adj = 0.05 / 4)
#   mean_est ~0.30   emp_sd ~0.085   mean_se ~0.086   power ~0.85


## 2. No effect: power should be about alpha_adj (type I error check)
set.seed(2)
power_itt_ate_chg(1000, 1/3, 0.25, sd_chg = 1.25, conf = 0.1, delta = 0,
          nsim = 1000, alpha_adj = 0.0125)
#   mean_est ~0   power ~0.0125


## 3. Effect of propensity score variation: more variation, less power
set.seed(3)
power_itt_ate_chg(1000, 1/3, p.cv = 0.25, sd_chg = 1.25, conf = 0, delta = 0.2,
          nsim = 200, alpha_adj = 0.0125)          # power ~0.42
power_itt_ate_chg(1000, 1/3, p.cv = 0.50, sd_chg = 1.25, conf = 0, delta = 0.2,
          nsim = 200, alpha_adj = 0.0125)          # power ~0.33


## 4. Effect of rho: higher correlation, smaller change SD, more power
for (r in c(0.3, 0.5, 0.7)) {
  set.seed(4)
  out = power_itt_ate_chg(1000, 1/3, 0.25, sd_chg = sd_change(1.25, 1.25, r),
                  conf = 0, delta = 0.2, nsim = 200, alpha_adj = 0.0125)
  cat("rho =", r, " power =", out$power, "\n")
}


## 5. With vs. without confounding (same everything else)
set.seed(5)
power_itt_ate_chg(1000, 1/3, 0.25, sd_chg = 1.25, conf = 0,   delta = 0.2,
          nsim = 200, alpha_adj = 0.0125)
power_itt_ate_chg(1000, 1/3, 0.25, sd_chg = 1.25, conf = 0.3, delta = 0.2,
          nsim = 200, alpha_adj = 0.0125)
# mean_est stays ~0.2 in both (IPTW removes the bias); power drops a little


## 6. Power curve over effect sizes
deltas = c(0.1, 0.15, 0.2, 0.25, 0.3)
curve = do.call(rbind, lapply(deltas, function(dl) {
  set.seed(6)
  cbind(delta = dl,
        power_itt_ate_chg(1000, 1/3, 0.25, sd_chg = 1.25, conf = 0.1, delta = dl,
                  nsim = 200, alpha_adj = 0.0125))
}))
curve
}