# Work Log: TvsQCD pass/fail stacks with alternate $U_T$ binnings (2022EE)

**Date:** 2026-08-04
**Input file:** `hadmonotop2022EE_0702.scaled`
**Base script:** `temp_plotting_2022EE.py`

## Background

`temp_plotting_2022EE.py` draws data/MC stack + Data/MC ratio plots per region,
summing over the `TvsQCD` axis (`{'TvsQCD': sum}`). This work splits that axis
into its two top-tag categories instead of summing over it, and repeats the
exercise across five candidate $U_T$ binnings.

Checked the `TvsQCD` axis of the input histograms directly:
```
Variable([0, 0.33, 1], underflow=False, overflow=False, name='TvsQCD')
```
So the axis has exactly two bins: index 0 = fail (score 0–0.33), index 1 = pass
(score 0.33–1). Same convention already used in `plot_stack_toptag.py` for the
2023 file.

## Scripts produced

### `temp_plotting_2022EE_toptag.py`
Direct adaptation of `temp_plotting_2022EE.py`: every `{'TvsQCD': sum, ...}`
selector is replaced with `{'TvsQCD': tag_idx, ...}`, looping `tag_idx` over
`{toptag_fail: 0, toptag_pass: 1}`. Uses the single fixed binning
`new_bin = [250, 400, 600, 1000]`. Everything else (stack composition per
region, colors, signal overlay, ratio panel, log-y, CMS label) is unchanged
from the original.

Output: `plots2022post_2022EE_0702_toptag_{fail,pass}/{region}_ut_2022EE_0702_bins250-400-600-1000_{fail,pass}.png`

### `temp_plotting_2022EE_toptag_binning.py`
Generalization of the above: the plotting body is wrapped in
`make_plots(case_name, new_bin)` and called once per binning case, so the
input file is only loaded/VV-merged once. Same pass/fail split, same stack
logic, same 8 regions (`sr, wmcr, tmcr, wecr, tecr, gcr, zmcr, zecr`).

Binning cases used:
| Case  | Edges |
|-------|-------|
| case1 | `[250, 300, ..., 1000]` (uniform 50 GeV steps) |
| case2 | `[250, 300, 350, 400, 450, 500, 550, 650, 1000]` |
| case3 | `[250, 300, 350, 400, 450, 500, 600, 700, 850, 1000]` |
| case4 | `[250, 350, 450, 550, 650, 1000]` |
| case5 | `[250, 400, 600, 1000]` (same as the fixed binning in the non-looping script) |

Output: `plots2022post_2022EE_0702_{case1..case5}_toptag_{fail,pass}/{region}_ut_2022EE_0702_bins{edges}_{fail,pass}.png`

## Run results

Both scripts were executed with the plotting conda env:
```bash
/Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python temp_plotting_2022EE_toptag.py
/Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python temp_plotting_2022EE_toptag_binning.py
```
- `temp_plotting_2022EE_toptag.py`: 16 plots (8 regions x pass/fail), all saved successfully.
- `temp_plotting_2022EE_toptag_binning.py`: 80 plots (5 binning cases x 8 regions x pass/fail), all saved successfully (exit code 0).

Known, pre-existing gaps (same behavior as the original `temp_plotting_2022EE.py`,
not introduced by this change): `zmcr` is missing `QCD Multijet`,
`Z ($\nu\nu$) + Jets`, `G + Jets`, and both signal samples; `zecr` is missing
`Z ($\nu\nu$) + Jets` and both signal samples; `gcr` is missing the
`sig_Mphi-200_Mchi-150` signal sample. These samples/regions are simply absent
from the input histograms and are skipped with a printed warning rather than
causing a failure.

## Output directories
```
plots2022post_2022EE_0702_toptag_fail/
plots2022post_2022EE_0702_toptag_pass/
plots2022post_2022EE_0702_case1_toptag_fail/
plots2022post_2022EE_0702_case1_toptag_pass/
plots2022post_2022EE_0702_case2_toptag_fail/
plots2022post_2022EE_0702_case2_toptag_pass/
plots2022post_2022EE_0702_case3_toptag_fail/
plots2022post_2022EE_0702_case3_toptag_pass/
plots2022post_2022EE_0702_case4_toptag_fail/
plots2022post_2022EE_0702_case4_toptag_pass/
plots2022post_2022EE_0702_case5_toptag_fail/
plots2022post_2022EE_0702_case5_toptag_pass/
```
