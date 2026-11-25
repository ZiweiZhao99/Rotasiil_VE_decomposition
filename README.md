# Decomposing the time-varying efficacy of a rotavirus vaccine

Analysis code for Zhao Z, Cheung YB. *Decomposing the Time-Varying Efficacy of a Rotavirus Vaccine into the Protection Provided by Vaccine-Induced Immunity and Impact of Reduction in Naturally Acquired Immunity: A Re-analysis of a Randomized Placebo-Controlled Trial of the Rotasiil Vaccine in Niger.*

## Data availability

The Rotasiil trial data (NCT02145000) are owned by Epicentre and are **not** included here. Requests for access may be submitted to Epicentre. To rerun the trial analyses, place `Baseline.xlsx` and `Outcome.xlsx` from the Epicentre data extraction in a folder of your choice. Then point `ROSE_DATA_DIR` at that folder (default `./data`, which is git-ignored).

The simulation needs no external data. Simulated datasets are regenerated from fixed seeds.

## Software

- R 4.5.2: survival 3.8.3, dplyr, tidyverse, readxl, ggplot2, scales, parallel.
- Stata 15.1 with `strmcure` (Xu & Cheung, *Stata J* 2015;15:135–54). The version used is in `strmcure_ado/`; copy it to your Stata `PERSONAL` ado folder.

## Contents and what each script reproduces

| Folder / script | Reproduces | Run log / output |
|---|---|---|
| `1_data_preparation/prepare_data.R` | Analysis data (`ROSE_ph.csv`), data split for Stata (`ROSE_split_tt.csv`), **Table 1** | `prepare_data.log`, `table1_episode_counts.csv` |
| `2_main_analysis/ch2roseCF_tt.do` (Stata) | Composite effect (Andersen–Gill, `stcox`), primary effect and frailty variance (conditional frailty, `strmcure`): **Table 2** AG/CF coefficients, θ̂ | `ch2roseCF_tt.log` |
| `2_main_analysis/model_selection_timeform.R` | Comparison of log t, t, √t, 1/t (and constant VE) by information criteria | `model_selection_timeform.log`, `.csv` |
| `3_cases_averted/rose_NCA.R` | Episodes averted per 1,000 children (composite, primary, secondary) | `rose_NCA.log`, `rose_NCA_monthly.csv` |
| `4_bootstrap/rose_boot_optionC.R` | 95% CI of the AUC-scale secondary effect, 1,000 bootstrap replicates, θ fixed at 0.8236 | `rose_boot_C_out/` |
| `4_bootstrap/boot_lower.R`, `boot_upper.R` | Sensitivity to θ fixed at 0.7507 and 0.8965: **Table S1** | `rose_boot_out_sensitivity_*/` |
| `5_figures/manuscript_figures.R` | **Figures 1–2**, first-event model and pointwise VE in **Table 2**, pointwise secondary-effect band | `manuscript_figures_runlog.txt`, `Table2_pointwise_VE.csv`, `Figure2_secondary_band.csv` |
| `6_simulation/stata/` | Simulation (thinning algorithm; θ = 1 and 0.2; 100 replicates each), Cox/AG/CF fits in Stata: **Table 3**, Stata rows of **Table S2** | `results_theta*/summary_stats.csv`, `results_all_*.csv`, per-replicate logs in `results_theta*/logs/`, `rep_one_tt_full.log` |
| `6_simulation/R_vs_Stata/compare_R_vs_Stata_CF.R` | R (`coxph` + `frailty`) fits on the same simulated datasets: R rows of **Table S2**; recurrence proportions quoted in Results | `compare_R_vs_Stata_CF.log` (fits), `compare_R_vs_Stata_CF_summary.log`, `R_sim_fits.csv`, `R_vs_Stata_CF_bias.csv`, `simulation_recurrence_proportion.csv` |

## How to run

```sh
export ROSE_DATA_DIR=/path/to/rotasiil/data     # holds Baseline.xlsx, Outcome.xlsx

cd 1_data_preparation && Rscript prepare_data.R > prepare_data.log 2>&1 && cd ..
# Stata: edit $datadir in 2_main_analysis/ch2roseCF_tt.do, then  do ch2roseCF_tt.do
cd 2_main_analysis   && Rscript model_selection_timeform.R > model_selection_timeform.log 2>&1 && cd ..
cd 3_cases_averted   && Rscript rose_NCA.R > rose_NCA.log 2>&1 && cd ..
cd 4_bootstrap       && Rscript rose_boot_optionC.R && Rscript boot_lower.R && Rscript boot_upper.R && cd ..
cd 5_figures         && N_CORES=4 Rscript manuscript_figures.R && cd ..
```

Simulation:

1. In Stata, run `6_simulation/stata/run_parallel_tt_full.do`. It launches `rep_one_tt_full.do` for replicates 1–100 with seed 100000 + r. Set `theta_true` in `rep_one_tt_full.do` to 1 or 0.2, and set the Stata executable and output paths at the top of the driver. Then run `summarize_partial.do <outdir>`.
2. Run `run_save_datasets.do`, which saves the same simulated datasets (same seeds) as CSV.
3. Run `SIM_DATA_DIR=<folder containing out_tt_full_1/ and out_tt_full_0.2/> Rscript 6_simulation/R_vs_Stata/compare_R_vs_Stata_CF.R`. Add `REUSE_FITS=1` to summarize the saved `R_sim_fits.csv` without refitting.

## Notes

- Memory: each `coxph(... tt())` fit expands the data to one row per subject per event time (~1–2 GB). On an 8 GB machine use `N_CORES=2` for `manuscript_figures.R` and `compare_R_vs_Stata_CF.R`, and do not run them at the same time.

- Point estimates of the AG and CF models, and θ̂, are from Stata (`ch2roseCF_tt.do`). In R the bootstrap refits both models with `coxph` + `tt()`, with θ fixed at the Stata estimate (Supplementary Material 1). The first-event Cox model in Table 2 and Figure 1 is fitted in R (`manuscript_figures.R`).
- `manuscript_figures.R` and `rose_boot_optionC.R` use the same resampling scheme and seeds (`set.seed(20260101 + b)`), so their AUC-scale bootstrap draws are identical.
- Absolute paths in the Stata simulation drivers (`run_parallel_*.do`, `summarize_partial.do`, and the default in `rep_one_tt_full.do`) refer to the machine on which the simulation was run, and need editing before rerunning.

## Reproduction status (R 4.5.2, November 2026)

| Manuscript quantity | Reproduced value | Source |
|---|---|---|
| Table 1 (all cells) | identical | `table1_episode_counts.csv` |
| Analysis / split data | byte-identical to the files analysed | `prepare_data.R` |
| AUC secondary effect, 95% CI | −0.0335 (−0.0572, −0.0149) | `rose_boot_C_out/`, `manuscript_figures_runlog.txt` (identical draws) |
| Table S1 | −0.031 (−0.054, −0.013); −0.036 (−0.061, −0.016) | `rose_boot_out_sensitivity_*/` |
| Secondary effect at 1, 6, 12, 20 months | −0.0135, −0.0296, −0.0384, −0.0460 | `Table2_pointwise_VE.csv` |
| First-event model | b0 −0.532, b1 0.156, AUC 0.189 | `manuscript_figures_runlog.txt` |
| Episodes averted per 1,000 | 51.0, 55.7, −4.7 | `rose_NCA.log` |
| Table 3 | identical | `6_simulation/stata/results_theta*/summary_stats.csv` |
| Table S2; R θ̂ ≈ 0 in 55% / 88% of datasets | identical | `compare_R_vs_Stata_CF_summary.log` |
| Recurrence among children with ≥1 episode, simulated | 8.1% (θ = 1), 5.1% (θ = 0.2) | `simulation_recurrence_proportion.csv` |
