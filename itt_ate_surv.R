## itt_ate_surv.R: power for the ITT ATE (marginal hazard ratio) on a
## time-to-event outcome. IPTW Cox model with a robust SE. Exponential event
## times; the only censoring is administrative, at Ymax.
source("power_core.R")

#' Hazard rate from cumulative risk (exponential model)
#'
#' Converts "a share p_event0 of untreated people have the event by year tau"
#' into the constant (exponential) hazard rate lambda, by solving
#' 1 - exp(-lambda * tau) = p_event0.
#'
#' @param p_event0 Cumulative incidence by tau among untreated people with an
#'                 average propensity score; strictly between 0 and 1.
#' @param tau      Length of follow-up (years); also the administrative
#'                 censoring time. Must be > 0.
#' @returns Hazard rate lambda (events per person-year).
#' @examples
#' rate_from_risk(251/2022, 10)   # 0.01325
rate_from_risk = function(p_event0, tau){
  # input checks
  if (p_event0 <= 0 || p_event0 >=1) stop("p_event0 should be strictly between 0 and 1")
  if (tau <= 0 ) stop("tau should be greater than 0")
  -log(1 - p_event0) / tau
}

#' Simulate potential survival outcomes
#'
#' Exponential event times with log hazard
#'   untreated: mu0 = b0 + b1 * ps_cv(e, p0, p.cv)
#'   treated:   mu1 = mu0 + log(TOR)
#' where b0 = log(rate_from_risk(py0, Ymax)) and b1 = log(HRy). TOR is a
#' conditional HR: the same for everyone at a given PS. Times are
#' administratively censored at Ymax (no other censoring).
#'
#' @param e    Propensity scores, from simulate_cohort()$e.
#' @param py0  Cumulative risk of the event by Ymax if untreated, for someone
#'             at the average PS (e = p0). Close to, but not exactly, the
#'             population risk when HRy != 1.
#' @param Ymax Follow-up length = administrative censoring time (years).
#' @param HRy  Confounding: HR per 1-CV increase in the PS. 1 = no confounding.
#' @param TOR  Conditional treatment hazard ratio.
#' @param p0   Mean propensity score (same as in get_beta_params()).
#' @param p.cv Propensity score CV (same as in get_beta_params()).
#' @returns data.frame(Y0, D0, Y1, D1): time (capped at Ymax) and event
#'          indicator, without and with treatment.
#' @examples 
#' po = sim_po_surv(co$e, 0.13, 7, 1.3, 0.85, 1/3, 0.5)
sim_po_surv = function(e, py0, Ymax, HRy = 1, TOR, p0, p.cv) {
  n  = length(e)
  b0 = log(rate_from_risk(py0, Ymax))  # log baseline hazard (at PS = p0)
  b1 = log(HRy)                        # confounding
  
  # untreated 
  mu0 = b0 + b1 * ps_cv(e, p0, p.cv)
  Y0  = -log(runif(n)) / exp(mu0)
  D0  = as.integer(Y0 < Ymax)
  
  # treated 
  mu1 = mu0 + log(TOR)
  Y1  = -log(runif(n)) / exp(mu1)
  D1  = as.integer(Y1 < Ymax)
  
  data.frame(Y0 = pmin(Y0, Ymax), D0 = D0,
             Y1 = pmin(Y1, Ymax), D1 = D1)
}

#' IPTW estimate and robust SE for the log hazard ratio
#'
#' Weighted Cox model with IPT weights, fit with survival::coxph() and a
#' robust (sandwich) SE. Gives the same estimate and robust SE as
#' survey::svycoxph() with svydesign(ids = ~1, weights = ~W) (up to a factor
#' sqrt(n/(n-1))), but is faster on large data.
#'
#' @param d Data frame from combine_outcomes(), with Y (follow-up time),
#'          D (event indicator: 1 = event, 0 = censored), X (0/1 treatment)
#'          and W (IPT weight).
#' @returns c(est = , se = ): estimated log hazard ratio and its robust SE.
#'          exp(est) is the hazard ratio.
est_itt_ate_surv = function(d) {
  fit = survival::coxph(survival::Surv(Y, D) ~ X, data = d,
                        weights = W, robust = TRUE)
  
  est = coef(fit)["X"]                        # log HR
  se  = sqrt(diag(vcov(fit)))["X"]            # robust SE (vcov is robust when robust = TRUE)
  
  c(est = unname(est), se = unname(se))
}

#' Simulate nsim datasets and compute power for the IPTW Cox hazard ratio
#'
#' @param n         Sample size.
#' @param p0        Expected proportion treated (allocation; 1:3 -> 1/4).
#' @param p.cv      Propensity score CV.
#' @param py0       Cumulative risk by Ymax if untreated, at the average PS.
#' @param Ymax      Follow-up length = administrative censoring time (years).
#' @param HRy       Confounding: HR per 1-CV increase in the PS, untreated.
#'                  1 = no confounding.
#' @param TOR       Conditional treatment HR used to generate the data.
#' @param nsim      Number of simulated datasets.
#' @param alpha_adj Two-sided significance level after multiplicity adjustment.
#' @returns One-row data.frame:
#'   mean_est     average log HR estimate
#'   marg_HR      exp(mean_est): the marginal HR = the detectable effect to
#'                report (~ TOR when HRy = 1; can differ slightly otherwise,
#'                because the HR is non-collapsible)
#'   emp_sd       SD of the log HR estimates (the true SE)
#'   mean_se      average robust SE (close to emp_sd = SE is accurate)
#'   mean_events  average number of events (both arms) per dataset
#'   power        share of datasets that reject H0: HR = 1
#'   
#' See the examples at the bottom of this file.
power_itt_ate_surv = function(n, p0, p.cv, py0, Ymax, HRy = 1, TOR,
                          nsim = 1000, alpha_adj = 0.05) {
  bp   = get_beta_params(p0, p.cv)
  crit = qnorm(1 - alpha_adj / 2)
  
  res = t(replicate(nsim, {
    co = simulate_cohort(n, bp["alpha"], bp["beta"], estimand = "ATE")
    po = sim_po_surv(e = co$e, py0 = py0, Ymax = Ymax, HRy = HRy, TOR = TOR,
                     p0 = p0, p.cv = p.cv)
    d  = combine_outcomes(co, po)
    c(est_itt_ate_surv(d), events = sum(d$D))
  }))
  
  data.frame(mean_est    = mean(res[, "est"]),
             marg_HR     = exp(mean(res[, "est"])),
             emp_sd      = sd(res[, "est"]),
             mean_se     = mean(res[, "se"]),
             mean_events = mean(res[, "events"]),
             power       = mean(abs(res[, "est"] / res[, "se"]) > crit))
}

# Examples
if (FALSE){
## 1. No effect: power should be about alpha_adj
set.seed(1)
power_itt_ate_surv(n = 2022, p0 = 1/4, p.cv = 0.25, py0 = 251/2022, Ymax = 10,
               HRy = 1, TOR = 1, nsim = 1000, alpha_adj = 0.05)
#   marg_HR ~1   mean_events ~251   power ~0.05 (0.035-0.065 is fine)

## 2. Basic call, no confounding: marg_HR should match TOR
set.seed(2)
out = power_itt_ate_surv(2022, p0 = 1/4, p.cv = 0.25, py0 = 251/2022, Ymax = 10,
                     HRy = 1, TOR = 0.62, nsim = 200, alpha_adj = 0.05)
out
#   marg_HR ~0.61   emp_sd ~ mean_se ~0.18   mean_events ~229   power ~0.82

## 3. With confounding: IPTW still recovers the HR; power about the same
set.seed(3)
power_itt_ate_surv(2022, p0 = 1/4, p.cv = 0.25, py0 = 251/2022, Ymax = 10,
               HRy = 2, TOR = 0.62, nsim = 1000, alpha_adj = 0.05)
#   marg_HR ~0.61   emp_sd ~ mean_se ~0.18   power ~0.83
}