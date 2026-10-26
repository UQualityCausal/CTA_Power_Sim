# CTA_Power_Sim

Simulation-based power analysis for propensity-score-weighted comparative treatment analyses (CTA).

Each script simulates cohorts with propensity scores, treatment, and outcomes; estimates the treatment effect with inverse probability of treatment weighting (IPTW) and a robust standard error; and reports power as the share of simulated datasets that reject the null (empirical power). Inputs and results are in the outcome's original units or ratio scale.

## Scenarios

Scripts are organized by **effect** (ITT = intention-to-treat, the effect of initiating a strategy; PP = per-protocol, the effect of adhering to it), **estimand** (ATE = average treatment effect; HTE = difference in treatment effects between levels of a binary effect modifier, the X × EF interaction), and **outcome**.

**Study design** (randomized or observational) is set by arguments, not by file: every script handles both.

| Study design | Set | Propensity score | Confounding |
|------------------|------------------|------------------|------------------|
| Observational | `p.cv > 0` | varies around `p0` | `HRy` (or `conf`) applies |
| Randomized trial, or cloned strategies | `p.cv = 0` | equals `p0` for everyone (allocation: 0.5 for 1:1 or cloning, 1/3 for 1:2) | none; `HRy` and `conf` are ignored |

In a randomized design the weights are constant, so IPTW gives the same estimate as the unweighted analysis. Per-protocol scripts still apply adherence and censoring weights, since deviation after randomization is not randomized.

| Outcome | Effect measure | ATE | HTE |
|------------------|------------------|------------------|------------------|
| Continuous change (baseline → one follow-up) | Mean difference | `itt_ate_chg.R` | `itt_hte_chg.R` |
| Time to event | Hazard ratio (Cox) | `itt_ate_surv.R` | `itt_hte_surv.R` |
| Count | Rate ratio (Poisson / negative binomial) | `itt_ate_count.R` | `itt_hte_count.R` |
| Count, clustered by hospital | Rate ratio; non-inferiority and other tests | `itt_ate_count_cl.R` | `itt_hte_count_cl.R` |
| Binary | Odds ratio (logistic) | `itt_ate_bin.R` | `itt_hte_bin.R` |
| Time to event, per-protocol (adherence and censoring weights) | Hazard ratio (weighted pooled logistic) | `pp_ate_surv.R` | `pp_hte_surv.R` |

**Sequential trial emulation**. (pending)

## Files

| File | Contents |
|------------------------------------|------------------------------------|
| `power_core.R` | Shared functions: propensity score parameters (`get_beta_params()`), cohort simulation (`simulate_cohort()`), PS in CV units (`ps_cv()`), effect modifier generation (`sim_ef()`), combining potential outcomes (`combine_outcomes()`). Defines functions only. |
| `itt_ate_<outcome>.R` | ATE scenario. Sources `power_core.R`. Holds any outcome-specific helpers (e.g. `sd_change()`, `rate_from_risk()`, `setup_count()`). |
| `itt_hte_<outcome>.R` | HTE scenario. Sources the matching ATE file (and through it `power_core.R`), so it reuses that file's helpers. |
| `pp_core.R` | Shared per-protocol functions: person-period builder (`sim_pp_surv()`) and stabilized censoring / adherence weights (`wt_pp()`). Defines functions only. |
| `pp_<estimand>_surv.R` | Per-protocol scenario. Sources the matching ITT survival file (for `sim_po_surv()` / `sim_po_surv_hte()`) and `pp_core.R`. |

Every scenario file defines three functions, plus an examples block:

| Function | Does |
|------------------------------------|------------------------------------|
| `sim_po_<outcome>()` / `sim_po_<outcome>_hte()` | Simulates potential outcomes (`Y0`/`Y1`, or `chg0`/`chg1`) for one dataset |
| `est_itt_<estimand>_<outcome>()` | Fits the IPTW model to one dataset; returns the estimate and robust SE |
| `power_itt_<estimand>_<outcome>()` | Runs `nsim` datasets; returns bias, SE and power checks |

Per-protocol scripts follow the same pattern (`est_pp_*`, `power_pp_*`) and add two shared steps between generation and estimation: `sim_pp_surv()` (person-period data) and `wt_pp()` (weights). They reuse the ITT data generators: the same baseline potential outcomes feed both designs.

## How it works

For each simulated dataset:

1.  **Propensity score and treatment.** Each person's propensity score (PS) is drawn from a Beta distribution with mean `p0` (the expected proportion treated) and coefficient of variation `p.cv`. Treatment is drawn from that probability, and ATE weights are computed: 1/e if treated, 1/(1−e) if not.

2.  **Potential outcomes.** Each person's outcome is simulated with and without treatment. Confounding links the outcome to the PS (see *Confounding scale*). The treatment effect is added to the treated outcome; in HTE scripts it depends on EF.

3.  **Observed outcome.** Each person's observed outcome is the potential outcome for the treatment they received.

4.  **Estimate.** IPTW regression with a robust (sandwich) SE:

    | Outcome | Model | Robust SE |
    |------------------------|------------------------|------------------------|
    | Change | `survey::svyglm()`, Gaussian | design-based |
    | Time to event | `survival::coxph(..., robust = TRUE)` | equals `survey::svycoxph()` up to √(n/(n−1)) |
    | Count | `glm()` Poisson | `sandwich::vcovHC()` HC0; valid for NB data (saturated model) |
    | Count, clustered | `glm()` Poisson or `MASS::glm.nb()` | `sandwich::vcovCL()` by hospital |
    | Binary | `glm()` logistic | `sandwich::vcovHC()` HC0; equals `svyglm(quasibinomial)` |
    | Time to event, per-protocol | `glm()` pooled logistic, weighted by stabilized IP weights | `sandwich::vcovCL()` clustered by person |

5.  **Power.** Repeat `nsim` times. Power is the share of datasets with \|estimate / SE\| above the two-sided critical value at `alpha_adj` (the multiplicity-adjusted α). The clustered count scripts instead use `decide()`, which also supports one-sided, non-inferiority, inferiority and equivalence tests.

Each `power_*()` function returns, per setting: the mean estimate (bias check), the empirical SD of the estimates vs. the mean robust SE (SE check: they should be close), event counts where relevant, and power.

## Per-protocol scripts

Per-protocol (PP) effects compare strategies as if everyone adhered for the whole follow-up. For each simulated dataset:

1.  **Baseline.** Cohort, PS and treatment as in ITT; outcome times from `sim_po_surv()` (or `sim_po_surv_hte()`), and death times from `sim_po_surv()` without confounding.
2.  **Person-period data** (`sim_pp_surv()`). Follow-up is split into `K` intervals of `dt_len` years (follow-up = `K * dt_len`, which also sets the time scale of `py0` and `py_death`). In each interval, in order: censoring → death → adherence → outcome event. A person leaves at the first of these; non-adherence is artificial censoring. A time-varying covariate `L` evolves as $L_t = a + rho·L_{t−1} + bX·X + noise$.
3.  **Weights** (`wt_pp()`). Pooled logistic models for remaining uncensored and for remaining adherent. Denominator: time spline + `L0` + `Lprev` + `X`; numerator: time spline only. The num/den ratios are multiplied over intervals within person, times the baseline IPTW, then truncated (default: 99th percentile).
4.  **Outcome model.** Weighted pooled logistic `Y_k ~ ns(time) + X` (`+ X:EF` for HTE) on the intervals still under protocol, with an SE clustered by person.

**Assumptions to know:**

- **Baseline confounding is real; deviation after baseline is not informative.** With `p.cv > 0` and `HRy` ≠ 1, the PS drives both treatment and outcome, and the baseline IPTW removes that bias. After baseline, `L` drives adherence and censoring but not the outcome, so leaving the protocol is unrelated to outcome risk and creates no bias. The adherence and censoring weights therefore have nothing to correct here, but they add the variance a real per-protocol analysis would have.
- Censoring and adherence in interval t depend on `Lprev`, the value of `L` entering the interval.
- Death is treated as censoring, without its own weights.

**Diagnostics** returned: events kept under protocol (`mean_events`; `mean_events11` for the treated EF = 1 cell in HTE), share adherent through follow-up (`mean_adh_end`), and the largest weight before truncation (`mean_max_sw`). For HTE, aim for `mean_events11` of at least \~50: with fewer, estimates of strong effects are biased and power is unreliable.

`glm.fit` warnings ("fitted probabilities numerically 0 or 1") come from the censoring model when censoring is near zero. They are expected and do not affect the results.

## Key conventions

### Confounding scale (`ps_cv()`)

Confounding strength (`HRy` for ratio outcomes, `conf` for change scores) is the effect **per 1-CV increase in the PS**: per increase of `p.cv` in the PS on the absolute probability scale, `(e − p0) / p.cv`. This keeps the original scripts' scaling. Note that one such unit equals 1/`p0` SDs of the PS, so the same `HRy` means weaker confounding per SD when `p0` is small (per SD, the HR is `HRy^p0`).

For count outcomes, the direction comes from `HRy` itself: `HRy` \> 1 means higher PS → larger outcome; `HRy` \< 1 means smaller. The clustered count scripts keep the original `sgn` argument instead.

### Effect scale: report the marginal effect

IPTW `Y ~ X` estimates the **marginal** effect (population if everyone treated vs. no one treated).

- **Change, count:** the effect you enter is the marginal effect (mean differences and rate ratios are collapsible). For counts, `mu0`/`mu1` are population means of Y0/Y1.
- **Time to event, binary:** `TOR` (and `MOR`) are **conditional** HRs/ORs (the same for everyone at a given PS), and `py0` is the risk at the average PS. HRs and ORs are non-collapsible, so the marginal effect (`marg_HR`, `marg_OR`, `marg_MOR`) can differ slightly from `TOR` when `HRy` ≠ 1. On the `ps_cv` scale the difference is small. **Report the `marg_*` value as the minimum detectable effect.**
- **Per-protocol:** `TOR` / `MOR` are conditional HRs, as in ITT survival. The pooled logistic OR per interval approximates the HR when per-interval risk is small. Report `marg_HR` / `marg_MOR`.

### Survival baseline

`py0` is the cumulative risk by `Ymax` (by `K * dt_len` in PP scripts); it is converted to a hazard with `rate_from_risk()`: λ = −log(1 − `py0`) / `Ymax`.

### Naming

- Files: `<design>_<estimand>_<outcome>[_cl].R`, e.g. `itt_hte_surv.R`; `_cl` = clustered.
- Functions: `<step>_<design>_<estimand>_<outcome>`, where step is `est` or `power`; data generators are `sim_po_<outcome>[_hte][_cl]` (they do not depend on the design, so PP scripts reuse them).
- PP functions: `sim_pp_surv()`, `wt_pp()`, `est_pp_<estimand>_surv()`, `power_pp_<estimand>_surv()`. Within PP code, `X` = treatment, `W` = baseline IPTW, `L` = time-varying covariate. (The original scripts used `W` for treatment and `X` for the covariate.)
- Every function name is unique across the library: `source()` silently overwrites a function defined twice.
- Potential outcomes are named `<name>0` / `<name>1`; `combine_outcomes()` creates `<name>`. No other column name should end in `0` or `1`.

## Examples

Examples inside function documentation (`@examples`) are short and may rely on objects created by earlier examples or functions in `power_core.R`. Full, self-contained examples, with expected results in comments, are in an `if (FALSE) {}` block at the end of each scenario script, so sourcing a script never runs a simulation. Each block includes a type I error check (no effect → power ≈ `alpha_adj`).

## Requirements

- R (≥ 4.1)
- [`survey`](https://cran.r-project.org/package=survey) (change score)
- [`survival`](https://cran.r-project.org/package=survival) (time to event)
- [`sandwich`](https://cran.r-project.org/package=sandwich) (count, binary, per-protocol)
- [`MASS`](https://cran.r-project.org/package=MASS) (clustered count, NB fits)

Packages are called with `pkg::function()` rather than `library()`, so sourcing a script does not attach anything.

## Data

Never commit real data. Only summary values (SDs, correlations, proportions, event rates) belong in settings.
