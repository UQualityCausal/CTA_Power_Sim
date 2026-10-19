## hte_surv.R
## Power for HTE on a time-to-event outcome: ratio of treatment hazard ratios
## between levels of a binary effect modifier (X:EF interaction), IPTW Cox.
source("power_core.R")

#' Simulate potential survival outcomes with effect modification
#'
#' Exponential event times with log hazard
#'   untreated: mu0 = b0 + b1 * (e - p0)/p.cv + log(HR_ef) * EF
#'   treated:   mu1 = mu0 + log(TOR) + log(MOR) * EF
#' so the treatment HR is TOR when EF = 0 and TOR * MOR when EF = 1.
#' Times are censored at Ymax.
#'
#' @param e     Propensity scores, from simulate_cohort()$e.
#' @param EF    Effect modifier, 0/1 per person, from sim_ef().
#' @param py0   Untreated event rate by Ymax, at PS = p0 and EF = 0.
#' @param Ymax  Follow-up length = administrative censoring time (years).
#' @param HRy   Confounding: HR per one-unit increase in (e - p0)/p.cv,
#'              untreated. 1 = no confounding.
#' @param TOR   Treatment hazard ratio in the EF = 0 group.
#' @param MOR   Ratio of treatment HRs, EF = 1 vs EF = 0 (the tested quantity).
#' @param HR_ef EF's own effect on the hazard, untreated. 1 = none.
#'              Does not affect power.
#' @param p0    Mean propensity score (same as in get_beta_params()).
#' @param p.cv  Propensity score CV (same as in get_beta_params()).
#' @returns data.frame(Y0, D0, Y1, D1): time (capped at Ymax) and event
#'          indicator, without and with treatment.
sim_po_surv_hte = function(e, EF, py0, Ymax, HRy = 1, TOR, MOR, HR_ef = 1,
                           p0, p.cv) {
  if (length(EF) != length(e)) stop("EF and e must have the same length")
  
  n  = length(e)
  b0 = log(py0 / Ymax)
  b1 = log(HRy)
  
  # untreated
  mu0 = b0 + b1*(e-p0)/p.cv + log(HR_ef) * EF   # baseline + confounding + EF main effect
  Y0  = -log(runif(n)) / exp(mu0)
  D0  = as.integer(Y0 < Ymax)
  
  # treated
  mu1 = mu0 + log(TOR) + log(MOR) * EF   # treated: effect depends on EF # mu0 + treatment effect that depends on EF
  Y1  = -log(runif(n)) / exp(mu1)
  D1  = as.integer(Y1 < Ymax)
  
  data.frame(Y0 = pmin(Y0, Ymax), D0 = D0,
             Y1 = pmin(Y1, Ymax), D1 = D1)
}

#' IPTW estimate and robust SE for the ratio of hazard ratios (HTE)
#'
#' Weighted Cox model Surv(Y, D) ~ X + EF + X:EF with IPT weights, fit with
#' survival::coxph() and a robust (sandwich) SE. Same estimate and robust SE as
#' survey::svycoxph() (up to sqrt(n/(n-1))), but faster on large data.
#'
#' @param d Data frame from combine_outcomes(), with Y (follow-up time),
#'          D (event indicator), X (0/1 treatment), EF (0/1 effect modifier)
#'          and W (IPT weight).
#' @returns c(est = , se = ): estimated log ratio of HRs and its robust SE.
#'          exp(est) is the ratio of HRs (MOR).
est_iptw_surv_hte = function(d) {
  fit = survival::coxph(survival::Surv(Y, D) ~ X * EF, data = d,
                        weights = W, robust = TRUE)
  
  est = coef(fit)["X:EF"]                     # log MOR
  se  = sqrt(diag(vcov(fit)))["X:EF"]         # robust SE
  
  c(est = unname(est), se = unname(se))
}

#' Simulate nsim datasets and compute power for the ratio of HRs (HTE)
#'
#' See the examples at the bottom of this file.
#'
#' @param n         Sample size.
#' @param pi_ef     Proportion of people with EF = 1, strictly between 0 and 1.
#' @param p0        Expected proportion treated (allocation; 1:3 -> 1/4).
#' @param p.cv      Propensity score CV.
#' @param py0       Untreated event rate by Ymax, at PS = p0 and EF = 0.
#' @param Ymax      Follow-up length = administrative censoring time (years).
#' @param HRy       Confounding: HR per one-unit increase in (e - p0)/p.cv,
#'                  untreated. 1 = no confounding.
#' @param TOR       Treatment hazard ratio in the EF = 0 group.
#' @param MOR       Ratio of treatment HRs, EF = 1 vs EF = 0 (tested).
#' @param HR_ef     EF's own effect on the hazard, untreated. 1 = none.
#' @param nsim      Number of simulated datasets.
#' @param alpha_adj Two-sided significance level after multiplicity adjustment.
#' @returns One-row data.frame:
#'   mean_est       average log(MOR) estimate
#'   marg_MOR       exp(mean_est)
#'   emp_sd         SD of the log(MOR) estimates (the true SE)
#'   mean_se        average robust SE (close to emp_sd = SE is accurate)
#'   mean_events    average total events per dataset
#'   mean_events11  average events among treated with EF = 1 (smallest cell)
#'   power          share of datasets that reject H0: MOR = 1
sim_power_surv_hte = function(n, pi_ef, p0, p.cv, py0, Ymax, HRy = 1, TOR, MOR,
                              HR_ef = 1, nsim = 1000, alpha_adj = 0.05) {
  bp   = get_beta_params(p0, p.cv)
  crit = qnorm(1 - alpha_adj / 2)
  
  res = t(replicate(nsim, {
    co    = simulate_cohort(n, bp["alpha"], bp["beta"], estimand = "ATE")
    co$EF = sim_ef(n, pi_ef)
    po    = sim_po_surv_hte(e = co$e, EF = co$EF, py0 = py0, Ymax = Ymax,
                            HRy = HRy, TOR = TOR, MOR = MOR, HR_ef = HR_ef,
                            p0 = p0, p.cv = p.cv)
    d     = combine_outcomes(co, po)
    c(est_iptw_surv_hte(d),
      events   = sum(d$D),                            # all events
      events11 = sum(d[d$X == 1 & d$EF == 1,"D"]))    # events among X = 1 & EF = 1
  }))
  
  data.frame(mean_est      = mean(res[, "est"]),
             marg_MOR      = exp(mean(res[, "est"])),
             emp_sd        = sd(res[, "est"]),
             mean_se       = mean(res[, "se"]),
             mean_events   = mean(res[, "events"]),
             mean_events11 = mean(res[, "events11"]),
             power         = mean(abs(res[, "est"] / res[, "se"]) > crit))
}

if (FALSE){
# Examples
## Shared settings (KP)
n    = 8797
p0   = 1/4
p.cv = 0.25
py0  = 0.13
Ymax = 7

## 1. No effect modification: type I error should be about alpha_adj
set.seed(1)
sim_power_surv_hte(n = n, pi_ef = 0.2, p0 = p0, p.cv = p.cv, py0 = py0,
                   Ymax = Ymax, HRy = 1, TOR = 0.85, MOR = 1,
                   nsim = 1000, alpha_adj = 0.05)
#   marg_MOR ~1   mean_events ~1035   mean_events11 ~46   power ~0.05

## 2. Basic call: MOR = 1.5, pi = 0.2
set.seed(2)
out = sim_power_surv_hte(n = n, pi_ef = 0.2, p0 = p0, p.cv = p.cv, py0 = py0,
                         Ymax = Ymax, HRy = 1, TOR = 0.85, MOR = 1.5,
                         nsim = 500, alpha_adj = 0.05)
out
#   marg_MOR ~1.5   mean_events11 ~67   emp_sd ~ mean_se ~0.17   power ~0.65

## 3. Effect of pi: a balanced modifier gives more power
set.seed(3)
sim_power_surv_hte(n = n, pi_ef = 0.5, p0 = p0, p.cv = p.cv, py0 = py0,
                   Ymax = Ymax, HRy = 1, TOR = 0.85, MOR = 1.5,
                   nsim = 500, alpha_adj = 0.05)
#   mean_events11 ~168   mean_se ~0.14   power ~0.80

## 5. With confounding (KP's old HRy = 1.3)
set.seed(4)
sim_power_surv_hte(n = n, pi_ef = 0.2, p0 = p0, p.cv = p.cv, py0 = py0,
                   Ymax = Ymax, HRy = 1.3, TOR = 0.85, MOR = 1.5,
                   nsim = 500, alpha_adj = 0.05)
#   marg_MOR ~1.5   power about the same as example 2

## 6. Power curve over MOR (to find the MOR with 80% power)
mors  = c(1.3, 1.4, 1.5, 1.6, 1.7)
curve = do.call(rbind, lapply(mors, function(m) {
  set.seed(5)
  cbind(MOR = m,
        sim_power_surv_hte(n = n, pi_ef = 0.2, p0 = p0, p.cv = p.cv,
                           py0 = py0, Ymax = Ymax, HRy = 1, TOR = 0.85,
                           MOR = m, nsim = 200, alpha_adj = 0.05))
}))
curve[, c("MOR", "marg_MOR", "mean_events11", "power")]
}