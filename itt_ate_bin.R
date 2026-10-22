source("power_core.R")

#' Simulate potential binary outcomes
#'
#' Logistic model, conditional on the PS:
#'   untreated: logit P(Y0 = 1) = b0 + b1 * ps_cv(e, p0, p.cv)
#'   treated:   logit P(Y1 = 1) = same + log(TOR)
#' where b0 = qlogis(py0) and b1 = log(HRy).
#'
#' @param e    Propensity scores, from simulate_cohort()$e.
#' @param py0  Risk if untreated for someone at the average PS (e = p0).
#'             Close to, but not exactly, the population risk when HRy != 1.
#' @param HRy  Confounding: OR per 1-CV increase in the PS (see ps_cv()).
#'             1 = no confounding.
#' @param TOR  Conditional treatment OR (same for everyone at a given PS).
#' @param p0   Mean PS (same as in get_beta_params()).
#' @param p.cv PS CV (same as in get_beta_params()).
#' @returns data.frame(Y0, Y1): 0/1 outcome if untreated and if treated.
sim_po_bin = function(e, py0, HRy = 1, TOR, p0, p.cv) {
  # input check: py0 strictly in (0, 1)
  if (py0 <= 0 || py0 >=1) stop("py0 should be strictly between 0 and 1")
  
  n   = length(e)
  b0 = qlogis(py0)
  b1 = log(HRy)
  lp0 = b0 + b1 * ps_cv(e, p0, p.cv)   # untreated log-odds
  lp1 = lp0 + log(TOR)                 # treated log-odds
  
  Y0 = rbinom(n, 1, plogis(lp0))
  Y1 = rbinom(n, 1, plogis(lp1))
  
  data.frame(Y0 = Y0, Y1 = Y1)
}

#' IPTW estimate and robust SE for the log marginal OR
#'
#' Weighted logistic regression Y ~ X with IPT weights and a robust (HC0) SE.
#'
#' @param d Data frame from combine_outcomes(), with Y (0/1), X (0/1), W.
#' @returns c(est = , se = ): log marginal OR and its robust SE.
est_itt_ate_bin = function(d) {
  fit = suppressWarnings(glm(Y ~ X, family = binomial, data = d, weights = W))

  est = coef(fit)["X"]
  se  = sqrt(sandwich::vcovHC(fit, type = "HC0")["X", "X"])
  c(est = unname(est), se = unname(se))
}

#' Simulate nsim datasets and compute empirical power for the IPTW OR
#'
#' @param n         Sample size.
#' @param p0        Expected proportion treated.
#' @param p.cv      PS CV.
#' @param py0       Risk if untreated at the average PS.
#' @param HRy       Confounding OR per 1-CV increase in the PS. 1 = none.
#' @param TOR       Conditional treatment OR used to generate the data.
#' @param nsim      Number of simulated datasets.
#' @param alpha_adj Two-sided significance level after multiplicity adjustment.
#' @returns One-row data.frame:
#'   TOR          input (conditional) OR
#'   mean_est     average log OR estimate
#'   marg_OR      exp(mean_est): the marginal OR = the detectable effect to report
#'                (= TOR when HRy = 1; closer to 1 otherwise)
#'   crude_OR     exp(mean unweighted log OR): shows confounding
#'   emp_sd       SD of the log OR estimates (the true SE)
#'   mean_se      average robust SE (close to emp_sd = SE is accurate)
#'   mean_events  average number of Y = 1 per dataset
#'   power        share of datasets that reject H0: OR = 1
power_itt_ate_bin = function(n, p0, p.cv, py0, HRy = 1, TOR,
                             nsim = 1000, alpha_adj = 0.05) {
  bp   = get_beta_params(p0, p.cv)
  crit = qnorm(1 - alpha_adj / 2)
  
  res = t(replicate(nsim, {
    co = simulate_cohort(n, bp["alpha"], bp["beta"], estimand = "ATE")
    po = sim_po_bin(co$e, py0 = py0, HRy = HRy, TOR = TOR, p0 = p0, p.cv = p.cv)
    d  = combine_outcomes(co, po)
    r1 = mean(d$Y[d$X == 1]); r0 = mean(d$Y[d$X == 0])        # crude arm risks
    c(est_itt_ate_bin(d),
      crude  = qlogis(r1) - qlogis(r0),
      events = sum(d$Y))
  }))
  
  data.frame(TOR         = TOR,
             mean_est    = mean(res[, "est"]),
             marg_OR     = exp(mean(res[, "est"])),
             crude_OR    = exp(mean(res[, "crude"])),
             emp_sd      = sd(res[, "est"]),
             mean_se     = mean(res[, "se"]),
             mean_events = mean(res[, "events"]),
             power       = mean(abs(res[, "est"] / res[, "se"]) > crit))
}

# Examples
if (FALSE){
  ## Shared settings: 1:3 allocation, 30% risk if untreated at the average PS
  n   = 2000
  p0  = 1/4
  py0 = 0.3
  
  ## 1. No effect: power should be about alpha_adj (type I error check)
  set.seed(1)
  power_itt_ate_bin(n = n, p0 = p0, p.cv = 0.25, py0 = py0, HRy = 2, TOR = 1,
                    nsim = 1000, alpha_adj = 0.05)
  #   marg_OR ~1   crude_OR ~1.06 (confounded)   power ~0.05 (0.035-0.065 is fine)
  
  ## 2. Basic call, no confounding: marg_OR should match TOR, and mean_se the
  ##    textbook SE sqrt(1/(n*p0*p1*(1-p1)) + 1/(n*(1-p0)*py0*(1-py0))) = 0.120,
  ##    with p1 = plogis(qlogis(py0) + log(TOR)) = 0.231
  set.seed(2)
  power_itt_ate_bin(n, p0, p.cv = 0.25, py0, HRy = 1, TOR = 0.7, nsim = 500)
  #   marg_OR ~0.70   emp_sd ~ mean_se ~0.12   mean_events ~566   power ~0.83
  
  ## 3. With confounding: IPTW removes the bias, the crude OR does not
  set.seed(3)
  power_itt_ate_bin(n, p0, p.cv = 0.25, py0, HRy = 2, TOR = 0.7, nsim = 500)
  #   marg_OR ~0.70 (non-collapsibility is tiny on the ps_cv scale)
  #   crude_OR ~0.74 (biased toward 1)   power ~0.86
  
  ## 4. Effect of propensity score variation: more variation, larger SE, less power
  set.seed(4)
  power_itt_ate_bin(n, p0, p.cv = 0.25, py0, HRy = 2, TOR = 0.7, nsim = 500)
  #   mean_se ~0.122   power ~0.87
  set.seed(4)
  power_itt_ate_bin(n, p0, p.cv = 0.50, py0, HRy = 2, TOR = 0.7, nsim = 500)
  #   mean_se ~0.135   power ~0.76
  
  ## 5. Power curve over TOR (to find the OR with 80% power)
  tors  = c(0.65, 0.70, 0.75, 0.80)
  curve = do.call(rbind, lapply(tors, function(t) {
    set.seed(5)
    power_itt_ate_bin(n, p0, p.cv = 0.25, py0, HRy = 2, TOR = t, nsim = 300)
  }))
  curve[, c("TOR", "marg_OR", "mean_se", "power")]
  #   power ~0.96, 0.85, 0.69, 0.45  -> report marg_OR at 80% power (~0.72)
}