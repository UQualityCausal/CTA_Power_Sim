## Power for per-protocol HTE on a time-to-event outcome: ratio of PP hazard
## ratios between levels of a binary effect modifier (X:EF interaction).
## Artificial censoring at non-adherence, stabilized IP weights, weighted
## pooled logistic regression with a person-clustered SE.
source("itt_hte_surv.R")   # sim_po_surv_hte(); via itt_ate_surv.R: sim_po_surv(), power_core.R
source("pp_core.R")        # sim_pp_surv(), wt_pp()

#' Weighted pooled logistic estimate of per-protocol effect modification
#'
#' Among rows still under protocol in the interval (uncensored, alive,
#' adherent), fits Y_k ~ ns(time) + X + EF + X:EF with weights sw.
#' X:EF = log ratio of PP hazard (odds) ratios, EF = 1 vs EF = 0.
#' SE: sandwich, clustered by person.
#'
#' @param pp      Data from wt_pp(), with id, time, X, EF, sw, C_k, D_k, A_k, Y_k.
#' @param df_time Spline df for the baseline hazard in time.
#' @returns c(est = , se = , est_x = ): log MOR (X:EF) and its clustered SE,
#'          plus the log PP HR in the EF = 0 group (X).
est_pp_hte_surv = function(pp, df_time = 4) {
  # keep the rows still under protocol in the interval (event rows included)
  dat = pp[pp$C_k == 0 & pp$D_k == 0 & pp$A_k == 1, ]
  
  # weighted pooled logistic model (warning = non-integer weights; harmless)
  fit = suppressWarnings(glm(Y_k ~ splines::ns(time, df = df_time) + X*EF,
                             family = binomial, data = dat, weights = sw))
  
  b  = coef(fit)
  V  = sandwich::vcovCL(fit, cluster = dat$id, type = "HC0")
  c(est   = unname(b["X:EF"]),          # log MOR
    se    = sqrt(V["X:EF", "X:EF"]),    # clustered by person
    est_x = unname(b["X"]))             # log PP HR, EF = 0
}

#' Simulate nsim datasets and compute power for per-protocol HTE
#'
#' @param n         Sample size simulated per dataset.
#' @param pi_ef     Proportion with EF = 1, strictly between 0 and 1.
#' @param p0        Expected proportion treated.
#' @param p.cv      Propensity score CV.
#' @param py0       Cumulative outcome risk by the end of follow-up if untreated,
#'                  at the average PS, EF = 0 group, ignoring death.
#' @param HRy       Confounding: HR per 1-CV increase in the PS. 1 = none.
#' @param TOR       Conditional PP HR in the EF = 0 group.
#' @param MOR       Ratio of conditional PP HRs, EF = 1 vs EF = 0 (tested).
#' @param HR_ef     EF's own effect on the hazard, untreated. 1 = none.
#'                  Changes power.
#' @param K, dt_len, py_death, pp_death, Lproc, cens, adh, df_time, trunc
#'                  As in power_pp_ate_surv().
#' @param nsim      Number of simulated datasets.
#' @param alpha_adj Two-sided significance level after multiplicity adjustment.
#' @returns One-row data.frame:
#'   MOR            input (conditional) ratio of PP HRs
#'   mean_est       average log MOR estimate
#'   marg_MOR       exp(mean_est): the detectable effect to report
#'   HR_ef0         exp(mean log HR, EF = 0): close to TOR
#'   HR_ef1         exp(mean log HR, EF = 1): close to TOR * MOR
#'   emp_sd         SD of the log MOR estimates
#'   mean_se        average clustered SE (close to emp_sd = SE is accurate)
#'   mean_events    average events counted under protocol
#'   mean_events11  average events among treated with EF = 1 (smallest cell)
#'   mean_adh_end   average share adherent through all K intervals
#'   mean_max_sw    average largest weight (before truncation)
#'   power          share of datasets that reject H0: MOR = 1
power_pp_hte_surv = function(n, pi_ef, p0, p.cv, py0, HRy = 1, TOR, MOR, HR_ef = 1,
                             K = 10, dt_len = 0.5, py_death, pp_death = 1,
                             Lproc = list(mean0 = 183, sd0 = 43, a = 2.5, rho = 0.9,
                                          bX = -0.5, sd = 1),
                             cens  = list(b0 = 0.3, bL = 0.06, bX = -0.1),
                             adh   = list(b0 = 0.1, bL = 0.015, bX = -0.07),
                             df_time = 4, trunc = 0.99,
                             nsim = 1000, alpha_adj = 0.05) {
  bp = get_beta_params(p0, p.cv)
  Ymax = K * dt_len
  crit = qnorm(1 - alpha_adj / 2)
  
  res = t(replicate(nsim, {
    co = simulate_cohort(n, bp["alpha"], bp["beta"], estimand = "ATE")
    co$EF = sim_ef(n, pi_ef)
    co$id = seq_len(n)
    po = sim_po_surv_hte(e = co$e, EF = co$EF, py0 = py0, Ymax = Ymax,
                         HRy = HRy, TOR = TOR, MOR = MOR, HR_ef = HR_ef,
                         p0 = p0, p.cv = p.cv) # Y0, D0, Y1, D1
    pd = sim_po_surv(co$e, py0 = py_death, Ymax = Ymax, HRy = 1, TOR = pp_death,
                     p0 = p0, p.cv = p.cv) # death, no confounding
    names(pd) = c("dt0", "dD0", "dt1", "dD1")
    d = combine_outcomes(co, cbind(po, pd))                         # Y, D, dt, dD
    d$id = seq_len(n)
    
    # person-period data, then attach EF by person, then weights
    pp = sim_pp_surv(d, K, dt_len, Lproc, cens, adh)
    pp$EF = d$EF[pp$id]          # attach each row's person-level EF
    pp = wt_pp(pp, df_time = df_time, trunc = trunc, cens = !is.null(cens))
    
    last = pp[!duplicated(pp$id, fromLast = TRUE), ]                # each person's last row
    c(est_pp_hte_surv(pp, df_time = df_time),
      events  = sum(pp$Y_k),
      adh_end = mean(last$time == K & last$A_k == 1 & last$Y_k == 0),
      max_sw  = attr(pp, "max_sw"),
      events11 = sum(pp$Y_k[pp$X == 1 & pp$EF == 1]))
  }))
  
  data.frame(MOR           = MOR,
             mean_est      = mean(res[, "est"]),
             marg_MOR      = exp(mean(res[, "est"])),
             HR_ef0        = exp(mean(res[, "est_x"])),
             HR_ef1        = exp(mean(res[, "est_x"] + res[, "est"])),
             emp_sd        = sd(res[, "est"]),
             mean_se       = mean(res[, "se"]),
             mean_events   = mean(res[, "events"]),
             mean_events11 = mean(res[, "events11"]),
             mean_adh_end  = mean(res[, "adh_end"]),
             mean_max_sw   = mean(res[, "max_sw"]),
             power         = mean(abs(res[, "est"] / res[, "se"]) > crit))
}

# Examples
if (FALSE){
  # Examples
  if (FALSE){
    ## Shared settings: 2,000 people, 20% treated, 30% with EF = 1,
    ## 20% 5-year outcome risk untreated, 12% 5-year death risk;
    ## original script's covariate, censoring, adherence (~18% adherent at 5 y)
    n = 2000; pi_ef = 0.3; p0 = 0.2; p.cv = 0.25; py0 = 0.2; HRy = 2; py_death = 0.12
    
    ## 1. No effect modification: power should be about alpha_adj
    set.seed(1)
    power_pp_hte_surv(n, pi_ef, p0, p.cv, py0, HRy, TOR = 0.8, MOR = 1,
                      py_death = py_death, nsim = 300)
    #   marg_MOR ~1.00   HR_ef0 ~ HR_ef1 ~0.80   emp_sd ~ mean_se ~0.42
    #   mean_events11 ~11   power ~0.05  
    
    ## 2. Basic call: MOR = 0.5 (HR 0.8 in EF = 0, 0.4 in EF = 1)
    set.seed(2)
    power_pp_hte_surv(n, pi_ef, p0, p.cv, py0, HRy, TOR = 0.8, MOR = 0.5,
                      py_death = py_death, nsim = 300)
    #   marg_MOR ~0.47   HR_ef0 ~0.78   HR_ef1 ~0.37   mean_se ~0.54
    #   mean_events11 ~6   power ~0.25
    
    ## 3. A balanced modifier gives more power (compare with example 2)
    set.seed(3)
    power_pp_hte_surv(n, pi_ef = 0.5, p0, p.cv, py0, HRy, TOR = 0.8, MOR = 0.5,
                      py_death = py_death, nsim = 300)
    #   mean_events11 ~9   mean_se ~0.46   power ~0.39
    
    ## 4. HR_ef changes power: more events in EF = 1 (compare with example 2)
    set.seed(4)
    power_pp_hte_surv(n, pi_ef, p0, p.cv, py0, HRy, TOR = 0.8, MOR = 0.5, HR_ef = 2,
                      py_death = py_death, nsim = 300)
    #   mean_events11 ~11   mean_se ~0.41   power ~0.41
    
    ## 5. Full adherence, no censoring: PP should match ITT
    set.seed(5)
    power_pp_hte_surv(n, pi_ef, p0, p.cv, py0, HRy, TOR = 0.8, MOR = 0.5,
                      py_death = py_death, cens = NULL,
                      adh = list(b0 = 20, bL = 0, bX = 0), nsim = 300)
    #   marg_MOR ~0.48   mean_se ~0.40   mean_events ~354   power ~0.47
    set.seed(5)
    power_itt_hte_surv(n, pi_ef, p0, p.cv, py0, Ymax = 5, HRy = HRy,
                       TOR = 0.8, MOR = 0.5, nsim = 300)
    #   marg_MOR ~0.48   mean_se ~0.38   mean_events ~378   power ~0.47
    #   (PP has slightly fewer events: death ends follow-up there)
    
    ## 6. Power curve over MOR
    mors  = c(0.3, 0.4, 0.5, 0.6)
    curve = do.call(rbind, lapply(mors, function(m) {
      set.seed(6)
      power_pp_hte_surv(n, pi_ef, p0, p.cv, py0, HRy, TOR = 0.8, MOR = m,
                        py_death = py_death, nsim = 200)
    }))
    curve[, c("MOR", "marg_MOR", "mean_events11", "mean_se", "power")]
    #   power ~0.50, 0.32, 0.25, 0.15
    #   Note marg_MOR ~0.13 at MOR = 0.3: with only ~3-7 events in the treated
    #   EF = 1 cell, strong effects are estimated with large small-sample bias.
    #   n = 2,000 is too small for PP HTE with this adherence; increase n until
    #   mean_events11 is at least ~50 before using results.
  }
}