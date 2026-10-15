# CTA_Power_Sim

Simulation-based power analysis for propensity-score-weighted comparative treatment analyses (CTA).

Each script simulates cohorts with propensity scores, treatment, and outcomes; estimates the treatment effect with inverse probability of treatment weighting (IPTW) and robust standard errors; and reports power as the share of simulated datasets that reject the null. All inputs and results are in the outcome's original units.

## Scenarios

- **Continuous change, ATE** (`ate_cts_change.R`): average treatment effect on the change from baseline to one follow-up
- **Continuous change, HTE** (`hte_cts_change.R`): difference in treatment effects between levels of a binary effect modifier (X × EF interaction)
- **Count outcome**: planned
- **Survival outcome**: planned

## Files

| File | Contents |
|------------------------------------|------------------------------------|
| `power_core.R` | Shared functions used by every scenario: propensity score parameters (`get_beta_params()`), cohort simulation (`simulate_cohort()`), standardized propensity score (`ps_z()`), effect modifier (`sim_ef()`), and combining potential outcomes (`combine_outcomes()`). Defines functions only. |
| `ate_*.R`, `hte_*.R` | One script per scenario: outcome-specific functions, settings, run, and save. Each script sources `power_core.R`. |

## How it works

For each simulated dataset:

1.  **Propensity score and treatment.** Each person's propensity score is drawn from a Beta distribution with mean `p0` (the expected proportion treated) and coefficient of variation `p.cv`. Treatment is drawn from that probability, and IPT weights are computed.
2.  **Potential outcomes.** Each person's outcome is simulated under no treatment and under treatment. Confounding (`conf`) links the untreated outcome to the propensity score; the treatment effect is added to the treated outcome.
3.  **Observed outcome.** Each person's observed outcome is the potential outcome for the treatment they received.
4.  **Estimate.** IPTW regression (`survey::svyglm`) with robust standard errors.
5.  **Power.** Repeat `nsim` times; power is the share of datasets where \|estimate / SE\| exceeds the critical value at the multiplicity-adjusted α.

Each scenario script reports, per setting: mean estimate (bias check), empirical SD of estimates vs. mean robust SE (SE check), power, and *the minimum detectable effect, alongside the same quantities for a randomized trial of equal size and allocation.*

## Examples

Examples inside function documentation (`@examples`) are short and may rely on objects or functions created by earlier examples in the same script or in `power_core.R`. Full, self-contained examples are given in a block at the end of each scenario script.

## Requirements

- R (≥ 4.1)
- [`survey`](https://cran.r-project.org/package=survey)

Packages are called with `pkg::function()` rather than `library()`, so sourcing a script does not attach anything.

## Conventions

- Potential outcomes are named `<name>0` / `<name>1` (e.g. `chg0`, `chg1`); `combine_outcomes()` creates `<name>`. No other column name should end in `0` or `1`.
- File names: `<estimand>_<outcome>_<setting>.R`.
- Never commit real data. Only summary values (SDs, correlations, proportions) belong in settings.
