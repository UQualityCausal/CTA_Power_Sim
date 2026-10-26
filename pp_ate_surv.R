## pp_ate_surv.R: power for the per-protocol ATE (hazard ratio) on a
## time-to-event outcome: artificial censoring at non-adherence, stabilized
## IP weights for censoring and adherence (times baseline IPTW), weighted
## pooled logistic regression with a person-clustered SE.

source("itt_ate_surv.R")   # sim_po_surv(), rate_from_risk(), power_core.R
source("pp_core.R")        # sim_pp_surv(), wt_pp()

#' Weighted pooled logistic estimate of the per-protocol effect
#'
#' Among rows still under protocol in the interval (uncensored, alive,
#' adherent), fits Y_k ~ ns(time) + X with weights sw. The X coefficient is
#' the log PP odds ratio per interval (~ log HR when interval risk is small).
#' SE: sandwich, clustered by person (a person's rows are not independent).
#'
#' @param pp      Data from wt_pp(), with id, time, X, sw, C_k, D_k, A_k, Y_k.
#' @param df_time Spline df for the baseline hazard in time.
#' @returns c(est = , se = ): log PP hazard (odds) ratio and its clustered SE.
#' 
est_pp_ate_surv = function(pp, df_time = 4) {
  # keep the rows still under protocol in the interval (event rows included)
  dat = pp[pp$C_k == 0 & pp$D_k == 0 & pp$A_k == 1, ]
  
  # weighted pooled logistic model (warning = non-integer weights; harmless)
  fit = suppressWarnings(glm(Y_k ~ splines::ns(time, df = df_time) + X,
                             family = binomial, data = dat, weights = sw))
  
  est = coef(fit)["X"]                                                  # log PP HR
  se  = sqrt(sandwich::vcovCL(fit, cluster = dat$id, type = "HC0")["X", "X"])  # clustered by person
  
  c(est = unname(est), se = unname(se))
}

#' Simulate nsim datasets and compute power for the per-protocol HR
#'
#' @param n         Sample size simulated per dataset.
#' @param p0        Expected proportion treated.
#' @param p.cv      Propensity score CV.
#' @param py0       Cumulative outcome risk by the end of follow-up (K * dt_len)
#'                  if untreated, at the average PS, ignoring death.
#' @param HRy       Confounding: HR per 1-CV increase in the PS. 1 = none.
#' @param TOR       Conditional per-protocol HR for the outcome.
#' @param K         Number of intervals.
#' @param dt_len    Interval length in years.
#' @param py_death  Cumulative death risk by the end of follow-up, untreated.
#' @param pp_death  Per-protocol HR for death.
#' @param Lproc     Covariate process, passed to sim_pp_surv().
#' @param cens      Censoring coefficients, passed to sim_pp_surv(); NULL = none.
#' @param adh       Adherence coefficients, passed to sim_pp_surv().
#' @param df_time   Spline df, passed to wt_pp() and est_pp_ate_surv().
#' @param trunc     Weight truncation quantile, passed to wt_pp().
#' @param nsim      Number of simulated datasets.
#' @param alpha_adj Two-sided significance level after multiplicity adjustment.
#' @returns One-row data.frame:
#'   TOR           input (conditional) PP hazard ratio
#'   mean_est      average log PP HR estimate
#'   marg_HR       exp(mean_est): the detectable effect to report
#'   emp_sd        SD of the estimates (at the simulated n)
#'   mean_se       average clustered SE
#'   mean_events   average number of events counted under protocol
#'   mean_adh_end  average share of people adherent through all K intervals
#'   mean_max_sw   average largest weight (before truncation)
#'   power         share of datasets that reject H0: HR = 1
power_pp_ate_surv = function(n, p0, p.cv, py0, HRy = 1, TOR,
                             K = 10, dt_len = 0.5, py_death, pp_death = 1,
                             Lproc = list(mean0 = 183, sd0 = 43, a = 2.5, rho = 0.9,
                                          bX = -0.5, sd = 1),
                             cens  = list(b0 = 0.3, bL = 0.06, bX = -0.1),
                             adh   = list(b0 = 0.1, bL = 0.015, bX = -0.07),
                             df_time = 4, trunc = 0.99,
                             nsim = 1000, alpha_adj = 0.05) {
  bp   = get_beta_params(p0, p.cv)
  Ymax = K * dt_len
  crit = qnorm(1 - alpha_adj / 2)
  
  res = t(replicate(nsim, {
    co = simulate_cohort(n, bp["alpha"], bp["beta"], estimand = "ATE")
    po = sim_po_surv(co$e, py0 = py0, Ymax = Ymax, HRy = HRy, TOR = TOR,
                     p0 = p0, p.cv = p.cv)                          # Y0, D0, Y1, D1
    pd = sim_po_surv(co$e, py0 = py_death, Ymax = Ymax, HRy = 1, TOR = pp_death,
                     p0 = p0, p.cv = p.cv)                          # death, no confounding
    names(pd) = c("dt0", "dD0", "dt1", "dD1")
    d = combine_outcomes(co, cbind(po, pd))                         # Y, D, dt, dD
    d$id = seq_len(n)
    
    pp = sim_pp_surv(d, K, dt_len, Lproc, cens, adh)
    pp = wt_pp(pp, df_time = df_time, trunc = trunc, cens = !is.null(cens))
    
    last = pp[!duplicated(pp$id, fromLast = TRUE), ]                # each person's last row
    c(est_pp_ate_surv(pp, df_time = df_time),
      events  = sum(pp$Y_k),
      adh_end = mean(last$time == K & last$A_k == 1 & last$Y_k == 0),
      max_sw  = attr(pp, "max_sw"))
  }))
  
  data.frame(TOR          = TOR,
             mean_est     = mean(res[, "est"]),
             marg_HR      = exp(mean(res[, "est"])),
             emp_sd       = sd(res[, "est"]),
             mean_se      = mean(res[, "se"]),
             mean_events  = mean(res[, "events"]),
             mean_adh_end = mean(res[, "adh_end"]),
             mean_max_sw  = mean(res[, "max_sw"]),
             power        = mean(abs(res[, "est"] / res[, "se"]) > crit))
}

# Examples
if (FALSE){
  ## Shared settings: 20% treated, 20% 5-year CVD risk untreated,
  ## 12% 5-year death risk; original script's covariate, censoring, adherence
  n = 5000; p0 = 0.2; p.cv = 0.25; py0 = 0.2; HRy = 2; py_death = 0.12
  
  ## 1. No effect: power should be about alpha_adj
  set.seed(1)
  power_pp_ate_surv(n, p0, p.cv, py0, HRy, TOR = 1, py_death = py_death, nsim = 300)
  #   marg_HR ~1.01   emp_sd ~ mean_se ~0.11   mean_adh_end ~0.18   power ~0.04
  #   (glm warnings from the near-zero censoring model are expected)
  
  ## 2. Basic call
  set.seed(2)
  power_pp_ate_surv(n, p0, p.cv, py0, HRy, TOR = 0.7, py_death = py_death, nsim = 200)
  #   marg_HR ~0.69   emp_sd ~ mean_se ~0.13   mean_events ~540   power ~0.86
  
  ## 3. Full adherence, no censoring: PP should match ITT
  set.seed(3)
  power_pp_ate_surv(n, p0, p.cv, py0, HRy, TOR = 0.7, py_death = py_death,
                    cens = NULL, adh = list(b0 = 20, bL = 0, bX = 0), nsim = 200)
  #   marg_HR ~0.70   mean_se ~0.095   mean_events ~890   mean_adh_end ~0.71   power ~0.99
  set.seed(3)
  power_itt_ate_surv(n, p0, p.cv, py0, Ymax = 5, HRy = HRy, TOR = 0.7, nsim = 200)
  #   marg_HR ~0.71   mean_se ~0.091   mean_events ~954   power ~0.97
  #   (PP has fewer events: death ends follow-up there, not in ITT)
  
  ## 4. Better adherence -> more events kept -> more power (compare with example 2)
  set.seed(4)
  power_pp_ate_surv(n, p0, p.cv, py0, HRy, TOR = 0.7, py_death = py_death,
                    adh = list(b0 = 2, bL = 0.015, bX = -0.07), nsim = 200)
  #   mean_adh_end ~0.55   mean_events ~800   mean_se ~0.10   power ~0.95
  
  ## 5. Power curve over TOR (to find the HR with 80% power)
  tors  = c(0.65, 0.70, 0.75, 0.80)
  curve = do.call(rbind, lapply(tors, function(t) {
    set.seed(6)
    power_pp_ate_surv(n, p0, p.cv, py0, HRy, TOR = t, py_death = py_death, nsim = 150)
  }))
  curve[, c("TOR", "marg_HR", "mean_se", "power")]
  #   power ~0.94, 0.83, 0.67, 0.45  -> 80% power at HR ~0.71
}