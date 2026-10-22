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
power_pp_ate_surv = function(n, p0, p.cv, py0, HRy = 1, TOR,
                             K = 10, dt_len = 0.5, py_death, pp_death = 1,
                             Lproc = list(mean0 = 183, sd0 = 43, a = 2.5, rho = 0.9,
                                          bX = -0.5, sd = 1),
                             cens  = list(b0 = 0.3, bL = 0.06, bX = -0.1),
                             adh   = list(b0 = 0.1, bL = 0.015, bX = -0.07),
                             df_time = 4, trunc = 0.99, fold = 1,
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
  
  se = res[, "se"] / sqrt(fold)                                     # scale to n * fold
  data.frame(TOR          = TOR,
             mean_est     = mean(res[, "est"]),
             marg_HR      = exp(mean(res[, "est"])),
             emp_sd       = sd(res[, "est"]),
             mean_se      = mean(se),
             mean_events  = mean(res[, "events"]),
             mean_adh_end = mean(res[, "adh_end"]),
             mean_max_sw  = mean(res[, "max_sw"]),
             power        = mean(abs(res[, "est"] / se) > crit))
}
