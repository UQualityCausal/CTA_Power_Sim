# CTA_Power_Sim

Simulation-based power analysis for propensity-score-weighted comparative treatment analyses (CTA).

Each script simulates cohorts with propensity scores, treatment, and outcomes; estimates the treatment effect with inverse probability of treatment weighting (IPTW) and a robust standard error; and reports power as the share of simulated datasets that reject the null (empirical power). Inputs and results are in the outcome's original units or ratio scale.

## Scenarios

Scripts are organized by **design** (ITT = intention-to-treat; PP = per-protocol), **estimand** (ATE = average treatment effect; HTE = difference in treatment effects between levels of a binary effect modifier, the X × EF interaction), and **outcome**.

| Outcome | Effect measure | ATE | HTE |
|----|----|----|----|
| Continuous change (baseline → one follow-up) | Mean difference | `itt_ate_chg.R` | `itt_hte_chg.R` |
| Time to event | Hazard ratio (Cox) | `itt_ate_surv.R` | `itt_hte_surv.R` |
| Count | Rate ratio (Poisson / negative binomial) | `itt_ate_count.R` | `itt_hte_count.R` |
| Count, clustered by hospital | Rate ratio; non-inferiority and other tests | `itt_ate_count_cl.R` | `itt_hte_count_cl.R` |
| Binary | Odds ratio (logistic) | `itt_ate_bin.R` | `itt_hte_bin.R` |
| Per-protocol (adherence and censoring weights) | — | planned | planned |

**Sequential trial emulation** needs no separate script. Following the original approach, each person contributes to one trial only, so it is a point-treatment analysis run with conservative inputs (first-trial sample size, half the follow-up, a conservative `p0`). To allow for adherence weights, inflate the SE afterwards (e.g. multiply `mean_se` by 1.1 and recompute power), or enter n / 1.1² as the sample size.

## Files

| File | Contents |
|----|----|
| `power_core.R` | Shared functions: propensity score parameters (`get_beta_params()`), cohort simulation (`simulate_cohort()`), PS in CV units (`ps_cv()`), effect modifier (`sim_ef()`), combining potential outcomes (`combine_outcomes()`). Defines functions only. |
| `itt_ate_<outcome>.R` | ATE scenario. Sources `power_core.R`. Holds any outcome-specific helpers (e.g. `sd_change()`, `rate_from_risk()`, `setup_count()`). |
| `itt_hte_<outcome>.R` | HTE scenario. Sources the matching ATE file (and through it `power_core.R`), so it reuses that file's helpers. |

Every scenario file defines three functions, plus an examples block:

| Function | Does |
|----|----|
| `sim_po_<outcome>()` / `sim_po_<outcome>_hte()` | Simulates potential outcomes (`Y0`/`Y1`, or `chg0`/`chg1`) for one dataset |
| `est_itt_<estimand>_<outcome>()` | Fits the IPTW model to one dataset; returns the estimate and robust SE |
| `power_itt_<estimand>_<outcome>()` | Runs `nsim` datasets; returns bias, SE and power checks |

## How it works

For each simulated dataset:

1.  **Propensity score and treatment.** Each person's propensity score (PS) is drawn from a Beta distribution with mean `p0` (the expected proportion treated) and coefficient of variation `p.cv`. Treatment is drawn from that probability, and ATE weights are computed: 1/e if treated, 1/(1−e) if not.

2.  **Potential outcomes.** Each person's outcome is simulated with and without treatment. Confounding links the outcome to the PS (see *Confounding scale*). The treatment effect is added to the treated outcome; in HTE scripts it depends on EF.

3.  **Observed outcome.** Each person's observed outcome is the potential outcome for the treatment they received.

4.  **Estimate.** IPTW regression with a robust (sandwich) SE:

    | Outcome | Model | Robust SE |
    |----|----|----|
    | Change | `survey::svyglm()`, Gaussian | design-based |
    | Time to event | `survival::coxph(..., robust = TRUE)` | equals `survey::svycoxph()` up to √(n/(n−1)) |
    | Count | `glm()` Poisson | `sandwich::vcovHC()` HC0; valid for NB data (saturated model) |
    | Count, clustered | `glm()` Poisson or `MASS::glm.nb()` | `sandwich::vcovCL()` by hospital |
    | Binary | `glm()` logistic | `sandwich::vcovHC()` HC0; equals `svyglm(quasibinomial)` |

5.  **Power.** Repeat `nsim` times. Power is the share of datasets with \|estimate / SE\| above the two-sided critical value at `alpha_adj` (the multiplicity-adjusted α). The clustered count scripts instead use `decide()`, which also supports one-sided, non-inferiority, inferiority and equivalence tests.

Each `power_*()` function returns, per setting: the mean estimate (bias check), the empirical SD of the estimates vs. the mean robust SE (SE check: they should be close), event counts where relevant, and power.

## Key conventions

### Confounding scale (`ps_cv()`)

Confounding strength (`HRy` for ratio outcomes, `conf` for change scores) is the effect **per 1-CV increase in the PS**: per increase of `p.cv` in the PS on the absolute probability scale, `(e − p0) / p.cv`. This keeps the original scripts' scaling. Note that one such unit equals 1/`p0` SDs of the PS, so the same `HRy` means weaker confounding per SD when `p0` is small (per SD, the HR is `HRy^p0`).

For count outcomes, the direction comes from `HRy` itself: `HRy` \> 1 means higher PS → larger outcome; `HRy` \< 1 means smaller. The clustered count scripts keep the original `sgn` argument instead.

### Effect scale: report the marginal effect

IPTW `Y ~ X` estimates the **marginal** effect (population if everyone treated vs. no one treated).

- **Change, count:** the effect you enter is the marginal effect (mean differences and rate ratios are collapsible). For counts, `mu0`/`mu1` are population means of Y0/Y1.
- **Time to event, binary:** `TOR` (and `MOR`) are **conditional** HRs/ORs (the same for everyone at a given PS), and `py0` is the risk at the average PS. HRs and ORs are non-collapsible, so the marginal effect (`marg_HR`, `marg_OR`, `marg_MOR`) can differ slightly from `TOR` when `HRy` ≠ 1. On the `ps_cv` scale the difference is small. **Report the `marg_*` value as the minimum detectable effect.**

### Survival baseline

`py0` is the cumulative risk by `Ymax`; it is converted to a hazard with `rate_from_risk()`: λ = −log(1 − `py0`) / `Ymax`. (The original scripts used `py0 / Ymax`, which treats `py0` as a cumulative hazard; results differ slightly from those.)

### Naming

- Files: `<design>_<estimand>_<outcome>[_cl].R`, e.g. `itt_hte_surv.R`; `_cl` = clustered.
- Functions: `<step>_<design>_<estimand>_<outcome>`, where step is `est` or `power`; data generators are `sim_po_<outcome>[_hte][_cl]` (they do not depend on the design, so PP scripts will reuse them).
- Every function name is unique across the library: `source()` silently overwrites a function defined twice.
- Potential outcomes are named `<name>0` / `<name>1`; `combine_outcomes()` creates `<name>`. No other column name should end in `0` or `1`.

## Usage

``` r
source("itt_ate_surv.R")
set.seed(1)
power_itt_ate_surv(n = 2022, p0 = 1/4, p.cv = 0.25, py0 = 251/2022, Ymax = 10,
                   HRy = 2, TOR = 0.62, nsim = 1000, alpha_adj = 0.05)
```

Run scripts from the library folder (paths in `source()` are relative), or use `here::here()`. To find the minimum detectable effect, run a grid of effect sizes and read off where power crosses 80% (or 90%); use nsim ≥ 1000 for final numbers. With nsim = 1000, the Monte Carlo SE of a power near 0.8 is about 0.013.

## Examples

Examples inside function documentation (`@examples`) are short and may rely on objects created by earlier examples or functions in `power_core.R`. Full, self-contained examples, with expected results in comments, are in an `if (FALSE) {}` block at the end of each scenario script, so sourcing a script never runs a simulation. Each block includes a type I error check (no effect → power ≈ `alpha_adj`).

## Requirements

- R (≥ 4.1)
- [`survey`](https://cran.r-project.org/package=survey) (change score)
- [`survival`](https://cran.r-project.org/package=survival) (time to event)
- [`sandwich`](https://cran.r-project.org/package=sandwich) (count, binary)
- [`MASS`](https://cran.r-project.org/package=MASS) (clustered count, NB fits)

Packages are called with `pkg::function()` rather than `library()`, so sourcing a script does not attach anything.

## Data

Never commit real data. Only summary values (SDs, correlations, proportions, event rates) belong in settings.
