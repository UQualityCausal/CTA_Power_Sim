## pp_core.R: shared functions for per-protocol (PP) power
## simulations: person-period data (sim_pp_surv) and stabilized censoring /
## adherence weights (wt_pp). Defines functions only: no runs, no seeds,
## no file output.
## Naming: X = treatment, W = baseline IPTW, L = time-varying covariate.

#' Build person-period data for a per-protocol analysis
#'
#' Within each interval t, in this order: censoring -> death -> adherence ->
#' outcome event. A person leaves at the first of these (non-adherence =
#' artificial censoring). L is then updated:
#'   L_t = a + rho * L_{t-1} + bX * X + N(0, sd^2).
#' Censoring and adherence in interval t depend on L_{t-1} (the value
#' entering the interval) and X.
#'
#' @param d      One row per person from combine_outcomes(), with id, X (0/1),
#'               W (baseline IPTW), Y, D (outcome time / event indicator),
#'               dt, dD (death time / indicator). Times capped at K * dt_len.
#' @param K      Number of intervals.
#' @param dt_len Interval length in years (follow-up = K * dt_len).
#' @param Lproc  Covariate process: list(mean0, sd0, a, rho, bX, sd).
#'               L0 ~ N(mean0, sd0^2); a, rho, bX, sd as in the formula above.
#' @param cens   Censoring: list(b0, bL, bX), coefficients on the logit of
#'               being UNcensored in an interval. NULL = no censoring.
#' @param adh    Adherence: list(b0, bL, bX), coefficients on the logit of
#'               being adherent in an interval.
#' @returns Long data.frame, one row per person per interval at risk:
#'   id, time (1..K), X, W, L0, Lprev (L entering the interval),
#'   C_k, D_k, A_k, Y_k (0/1 indicators for this interval).
#' @examples
#' pp = sim_pp_surv(d, K = 10, dt_len = 0.5,
#'                  Lproc = list(mean0 = 183, sd0 = 43, a = 2.5, rho = 0.9, bX = -0.5, sd = 1),
#'                  cens  = list(b0 = 0.3, bL = 0.06, bX = -0.1),
#'                  adh   = list(b0 = 0.1, bL = 0.015, bX = -0.07))
#' head(pp)
sim_pp_surv = function(d, K, dt_len, Lproc, cens = NULL, adh) {
  n = nrow(d)
  
  # 1. Covariate trajectories
  # draw L0 for everyone
  L0 = rnorm(n, mean = Lproc$mean0, sd = Lproc$sd0)
  # simulate L for intervals 1..K from the update formula, for everyone (whether or not they are still followed)
  Lprev = matrix(NA, nrow = n, ncol = K)
  Lprev[,1] = L0
  for (i in seq_len(K)[-1]){
    Lprev[,i] =  Lproc$a +Lproc$rho * Lprev[,i-1] + Lproc$bX * d$X + rnorm(n, 0, Lproc$sd)
  }
  
  # 2. Per-interval draws (n x K)
  # censoring: cens models the logit of being UNcensored
  if (is.null(cens)){
    Cm = matrix(0L, nrow = n, ncol = K)
  } else {
    pC = 1 - plogis(cens$b0 + cens$bL * Lprev + cens$bX * d$X)
    Cm = matrix(rbinom(n*K, 1, pC), nrow = n, ncol = K)
  }
  
  # adherence
  pA = plogis(adh$b0 + adh$bL * Lprev + adh$bX * d$X)
  Am = matrix(rbinom(n*K, 1, pA), nrow = n, ncol = K)
  
  # dead / outcome by the end of each interval
  tend = seq_len(K) * dt_len
  Dm = outer(d$dt, tend, "<=") & (d$dD == 1)
  Ym = outer(d$Y,  tend, "<=") & (d$D  == 1)
  
  # 3. Exit interval
  # leaves in interval t for any reason (only the first TRUE row matters)
  stop = (Cm == 1) | Dm | (Am == 0) | Ym
  # first flagged interval; K if never flagged
  exit = ifelse(rowSums(stop) > 0, max.col(stop, ties.method = "first"), K)
  
  # 4. Long format: one row per person for intervals 1..exit (id, time)
  id = rep(seq_len(n), exit)
  time = sequence(exit)
  ix = cbind(id, time) # (person, interval) for each row
  # indicators, in the order censoring -> death -> adherence -> event
  C_k = Cm[ix] == 1
  D_k = !C_k & Dm[ix]
  A_k = !C_k & !D_k & (Am[ix] == 1)
  Y_k = A_k & Ym[ix]
  
  data.frame(id = d$id[id], time = time,
             X = d$X[id], W = d$W[id], L0 = L0[id], Lprev = Lprev[ix],
             C_k = as.integer(C_k), D_k = as.integer(D_k),
             A_k = as.integer(A_k), Y_k = as.integer(Y_k))
}

#' Stabilized inverse probability weights for censoring and non-adherence
#'
#' Pooled logistic models, per interval:
#'   censoring (all rows; everyone is uncensored entering an interval)
#'     denominator: C_k ~ ns(time) + L0 + Lprev + X ; numerator: C_k ~ ns(time)
#'   adherence (rows uncensored and alive in the interval)
#'     denominator: A_k ~ ns(time) + L0 + Lprev + X ; numerator: A_k ~ ns(time)
#' Ratios num/den are multiplied over intervals within person, times the
#' baseline IPTW W, then truncated.
#'
#' @param pp      Long data from sim_pp_surv(), sorted by id and time.
#' @param df_time Degrees of freedom for the natural spline in time.
#' @param trunc   Lower and upper truncation quantiles for the final weight,
#'                e.g. c(0.01, 0.99); c(0, 1) = none. A single number (e.g.
#'                0.99) is read as c(0, 0.99), upper truncation only.
#' @param num     Extra numerator terms, as a formula string (e.g. "+ X").
#'                Only variables that are also in the outcome model. "" = time only.
wt_pp = function(pp, df_time = 4, trunc = c(0.01, 0.99), cens = TRUE, num = "") {
  if (length(trunc) == 1) trunc = c(0, trunc)          # old calls: upper only
  tf = function(rhs) as.formula(paste("~ splines::ns(time, df =", df_time, ")", rhs))
  r_c = rep(1, nrow(pp))
  if (cens){
    den_c = glm(update(tf("+ L0 + Lprev + X"), C_k ~ .), binomial, data = pp)
    num_c = glm(update(tf(num), C_k ~ .), binomial, data = pp)
    r_c = (1 - fitted(num_c)) / (1 - fitted(den_c)) # P(uncensored): model predicts P(censored)
  }
  sub = pp$C_k == 0 & pp$D_k == 0
  den_a = glm(update(tf("+ L0 + Lprev + X"), A_k ~ .), binomial, data = pp[sub,])
  num_a = glm(update(tf(num), A_k ~ .), binomial, data = pp[sub,])
  r_a = rep(1, nrow(pp))
  r_a[sub] = fitted(num_a) / fitted(den_a) # P(adherent): model predicts P(adherent)
  
  # final weight = W x cumulative censoring ratio x cumulative adherence ratio
  sw = pp$W * ave(r_c, pp$id, FUN = cumprod) * ave(r_a, pp$id, FUN = cumprod)
  attr(pp, "max_sw") = max(sw)          # largest weight before truncation
  # truncate at the trunc quantiles
  lim = quantile(sw, trunc)
  pp$sw = pmin(pmax(sw, lim[1]), lim[2])
  pp
}

# Examples
if (FALSE){
  source("itt_ate_surv.R")   # sim_po_surv(), power_core.R
  
  ## Settings (defaults from the original SMBPT script)
  n = 10000; K = 10; dt_len = 0.5; p0 = 0.2; p.cv = 0.25
  Lproc = list(mean0 = 183, sd0 = 43, a = 2.5, rho = 0.9, bX = -0.5, sd = 1)
  cens  = list(b0 = 0.3, bL = 0.06, bX = -0.1)
  adh   = list(b0 = 0.1, bL = 0.015, bX = -0.07)
  
  ## One person per row: cohort + outcome and death potential outcomes
  set.seed(1)
  bp = get_beta_params(p0, p.cv)
  co = simulate_cohort(n, bp["alpha"], bp["beta"])
  po = sim_po_surv(co$e, py0 = 0.2, Ymax = K * dt_len, HRy = 2, TOR = 0.97,
                   p0 = p0, p.cv = p.cv)                        # Y0, D0, Y1, D1
  pd = sim_po_surv(co$e, py0 = 0.12, Ymax = K * dt_len, HRy = 1, TOR = 0.95,
                   p0 = p0, p.cv = p.cv)                        # death
  names(pd) = c("dt0", "dD0", "dt1", "dD1")
  d = combine_outcomes(co, cbind(po, pd))                       # Y, D, dt, dD
  d$id = seq_len(n)
  
  ## 1. Build person-period data
  pp = sim_pp_surv(d, K, dt_len, Lproc, cens, adh)
  head(pp, 12)
  #   person 1: adherent through all 10 intervals, Lprev falling from ~285
  #   person 2: A_k = 0 in interval 2 (stops adhering), no rows after
  nrow(pp)        # ~58,500 person-intervals for 10,000 people
  
  ## 2. Why each person left (from their last row)
  last   = pp[!duplicated(pp$id, fromLast = TRUE), ]
  reason = with(last, ifelse(C_k == 1, "censored",
                             ifelse(D_k == 1, "died",
                                    ifelse(A_k == 0, "non-adherent",
                                           ifelse(Y_k == 1, "event", "completed")))))
  round(prop.table(table(reason)), 3)
  #   censored 0.007   completed 0.184   died 0.078   event 0.112
  #   non-adherent 0.619
  
  ## 3. Sanity: no exit reason before a person's last row
  early = pp[duplicated(pp$id, fromLast = TRUE), ]
  sum(early$C_k | early$D_k | early$A_k == 0 | early$Y_k)   # 0
  
  ## 4. Full adherence and no censoring: everyone stays until event, death, or K
  pp_full = sim_pp_surv(d, K, dt_len, Lproc, cens = NULL,
                        adh = list(b0 = 20, bL = 0, bX = 0))
  mean(pp_full$A_k == 0)                                     # 0
}