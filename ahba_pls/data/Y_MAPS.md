# Imaging Y maps

NOTES
- Built by `code/01_build_y_matrix.py` from the settled 7.0 runs (see README Design table).
- dCT/CT from `thickness_dsk_70_139406217085`; dT1T2/T1T2 from `t1t2_ratio_dsk_70_9a62dde44370`; slopePC1-3 from `docs/developmental_maps_noglobal.csv`.
- `dX` = fixed-effect `age_c` slope (per year), `X` = `(Intercept)` = fitted value at the run's mean age for females with random effects at zero. Age centres: CT run 12.7973, T1T2 run 12.8307 (sample-mean centring; config `age_centre: null`).
- dCT/CT agree with docs slope_total/baseline_thickness to <4.8e-06.
- `y_maps_68.csv`: 68 rows (34 lh_ then 34 rh_, ggseg DK labels). `y_maps_bilateral_34.csv`: LH/RH mean, indexed by the lh_ label.
- T1T2 run: rh_temporalpole fit is singular/non-converged (fits/diagnostics.parquet); its estimates are kept but that region's T1T2/dT1T2 values are less reliable.
- `results/y_map_lr_agreement.csv`: LH-vs-RH Pearson r and Spearman rho per map across the 34 regions.
- `results/y_map_correlations.csv`: bilateral map correlations, Pearson below the diagonal, Spearman above.

## LH vs RH agreement

| map | pearson_r | spearman_rho |
|---|---|---|
| dCT | 0.962 | 0.95 |
| CT | 0.992 | 0.985 |
| dT1T2 | 0.993 | 0.994 |
| T1T2 | 0.995 | 0.991 |
| slopePC1 | 0.974 | 0.968 |
| slopePC2 | 0.989 | 0.976 |
| slopePC3 | 0.983 | 0.98 |

## Bilateral map correlations (Pearson lower / Spearman upper)

| pearson_lower__spearman_upper | dCT | CT | dT1T2 | T1T2 | slopePC1 | slopePC2 | slopePC3 |
|---|---|---|---|---|---|---|---|
| dCT | 1.0 | 0.17 | -0.09 | 0.39 | 0.27 | -0.21 | -0.48 |
| CT | 0.23 | 1.0 | -0.66 | -0.26 | 0.12 | 0.54 | -0.04 |
| dT1T2 | -0.14 | -0.59 | 1.0 | 0.39 | -0.33 | -0.59 | 0.36 |
| T1T2 | 0.31 | -0.31 | 0.31 | 1.0 | -0.26 | -0.83 | -0.13 |
| slopePC1 | 0.27 | 0.17 | -0.36 | -0.33 | 1.0 | 0.37 | -0.53 |
| slopePC2 | -0.16 | 0.63 | -0.51 | -0.79 | 0.35 | 1.0 | 0.01 |
| slopePC3 | -0.43 | 0.02 | 0.38 | -0.07 | -0.54 | 0.02 | 1.0 |
