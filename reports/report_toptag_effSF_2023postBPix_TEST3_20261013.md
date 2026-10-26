# Top-tag α, MC efficiency, Data efficiency, SF and uncertainties — 2023postBPix (TEST3)

**Input:** `hadmonotop2023BPix_toptagTEST3.scaled` (updated 2026-10-13 13:12)
**Report date:** 2026-10-13
**Tagger / WP:** particleTransformer `TvsQCD` > 0.33, leading AK15 jet
**Luminosity:** 2023postBPix 9.68 fb⁻¹ (13.6 TeV)
**PDF version (with figures):** `reports/report_toptag_effSF_2023postBPix_TEST3_20261013.pdf`
**Companions:** `report_toptag_effSF_2023preBPix_TEST3_20260930.md`, `report_toptag_effSF_2022preEE_TEST3_20260930.md`, `report_toptag_effSF_2022postEE_TEST3_20261005.md` (same method)

Definitions follow the AN (AN-2020/061 Sec. 6.6): α and the Top-jet class yields use **all background MC**.

## Summary

- **α depends strongly on pT in both CRs.** γ CR: 0.38 / 0.33 / 0.56 / 0.96. tt̄ CR (μ+e): 0.67 / 0.61 / 0.73 / 0.94. In 2023preBPix both were roughly flat (γ 0.47–0.57, tt̄ 0.82–0.94).
- **Two separate normalization effects, both checked per fb⁻¹ against 2023preBPix:**
  - **Data deficit below 400 GeV.** The MC scales with luminosity (0.94–1.03) but the data does not: 0.69–0.76 in the γ CR and 0.72–0.75 in the tt̄ CR. The worst case is `tmcr` at 250–325 GeV (0.58). This is the same kind of data-side effect seen in 2022postEE. Possible causes: the lumi mask / lumi value, missing data files, or a BPix-specific veto applied differently to data and MC.
  - **G + Jets MC deficit at 600–2000 GeV.** G + Jets per fb⁻¹ is 0.62× the 2023preBPix value in that bin, against 0.99 below. The γ CR α therefore jumps to 0.96 there. A G + Jets sample (high HT or pT) is probably missing or incomplete in this file.
- **SF_QCD (γ CR):** 1.13 ± 0.47, 2.85 ± 0.48, 2.11 ± 0.29, 2.29 ± 0.38. The 325–2000 GeV bins agree with 2023preBPix (2.5–2.8) within 1σ. The 250–325 GeV bin has only 19 tagged data events (about half the per-fb⁻¹ rate of 2023preBPix), so its σ is 41%.
- **SF_Top (tt̄ CR, μ+e), SF_QCD-corrected (in the JSON):** 1.03 ± 0.38, 0.44 ± 0.20, 0.69 ± 0.10, 0.71 ± 0.16. With raw subtraction: 1.06, 0.70, 0.77, 0.85. The two upper bins agree with 2023preBPix (0.72, 0.74). The 325–400 GeV bin sits 0.8σ below it, and the correction lowers it by 0.26, the largest shift of any bin.
- **The subtraction is very large at low pT.** α·(N(B) + N(QCD)) removes 83% / 78% of the pre-tag data in the two lowest bins, leaving 39 / 101 events in the eff_Data denominator.
- **Uncertainties:** below 400 GeV, data statistics are the largest source (0.185 / 0.096 on Top-jet eff_Data). Above 400 GeV, SF_B = 1 ± 0.5 is the largest (0.067 / 0.101).
- **No high-weight MC outlier like the one in 2022postEE.** QCD Multijet in `tmcr` has weights ≈ 4 (21 ± 9 events at 325–400 GeV). It adds at most 0.044 to eff_Data.
- **eff_MC agrees with 2023preBPix** (Top-jet 0.68–0.80, QCD-jet 1.9–3.1%). The SF differences therefore come from data and from the normalization effects above.
- **JSON written:** `plots_effSF_toptagTEST3/sf_json/toptag_sf_2023postBPix.json`. It contains every bin as measured, including 250–325 and 325–400 GeV.

## 1. Steps

All steps run with the `plot_env` interpreter (Python 3.10) and take the `.scaled` file with `-f`. The `2023BPix` in the file name maps to era 2023postBPix (9.68 fb⁻¹). Outputs go to `plots_effSF_toptagTEST3/<kind>/2023postBPix/`.

```bash
cd /Users/hongjieun/working/toptag
PY=/Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python
F=hadmonotop2023BPix_toptagTEST3
$PY studies/plot_alpha_toptagTEST.py           -f $F   # Step 1: α
$PY studies/plot_toptag_eff_toptagTEST.py      -f $F   # Step 2: eff_MC
$PY studies/plot_toptag_eff_data_toptagTEST.py -f $F   # Step 3: eff_Data + SF (+ sf_json)
```

Full console output: `plots_effSF_toptagTEST3/run_log_2023postBPix.txt`. SF file: `plots_effSF_toptagTEST3/sf_json/toptag_sf_2023postBPix.json`. The composition, per-source uncertainty breakdown (Section 6), per-fb⁻¹ comparison (Section 3) and MC outlier scan (Section 5.3) come from a one-off script. It re-uses `measure()` and the yield functions from Step 3, and its totals reproduce `measure()` exactly.

| Step | Script | Reads | Computes | Output |
|---|---|---|---|---|
| 1 | `plot_alpha_toptagTEST.py` | `fjpt_bin` (bkg + data) | α per region and pT bin: γ CR, tt̄ CR μ+e, μ/e split | `alpha/2023postBPix/` |
| 2 | `plot_toptag_eff_toptagTEST.py` | `fjpt_cat` (bkg) | eff_MC: Top-jet (tmcr, tecr, μ+e, all bkg), QCD-jet (gcr, G + Jets) | `eff/2023postBPix/` |
| 3 | `plot_toptag_eff_data_toptagTEST.py` | both | eff_Data and SF: QCD-jet first, then Top-jet (raw and SF_QCD-corrected); correctionlib JSON | `effdata/2023postBPix/`, `sf_json/` |

### 1.1 Inputs

| Histogram | Axes | Filled for | Used for |
|---|---|---|---|
| `fjpt_bin` | region × pT (250, 325, 400, 600, 2000) × TvsQCD | bkg + data | N(Data), N(Data, tag), N(MC) → α |
| `fjpt_cat` | region × pT (same edges) × TvsQCD × n_qinfj × n_binfj | bkg only | class-split MC yields, eff_MC |

TvsQCD edges [0, 0.33, 1]: bin 1 = tagged. Data per region: `tmcr` → MET, `tecr` and `gcr` → EGamma. The `fjpt_bin`/`fjpt_cat` MC consistency check passed in every region.

Jet classes: **Top-jet** = classes 5–8, **B-jet** = 2, 4, **QCD-jet** = 1, 3.

### 1.2 Formulas

- **α** = N(Data) / N(MC), pre-tag, per region and pT bin; N(MC) = all background processes, all jet classes.
- **eff_MC** = N(MC, class, tag) / N(MC, class).
- **QCD-jet (γ CR, G + Jets):** eff_Data = [N(Data, tag) − α·SF_B·N(GJ, B, tag)] / [N(Data) − α·N(GJ, B)]
- **Top-jet (tt̄ CR, all bkg):** eff_Data = [N(Data, tag) − α·(SF_B·N(B, tag) + SF_QCD·N(QCD, tag))] / [N(Data) − α·(N(B) + N(QCD))]
  - *raw*: SF_QCD = 1; *SF_QCD-corrected*: SF_QCD from the γ CR, per pT bin.
- **SF** = eff_Data / eff_MC.

## 2. Sample composition

Fraction (%) of the pre-tag background MC (250–325 / 325–400 / 400–600 / 600–2000 GeV):

| Process | tt̄ CR (μ+e) | γ CR |
|---|---|---|
| TT | 34.6 / 35.5 / 43.3 / 40.1 | 0 / 0 / 0 / 0.1 |
| ST | 9.2 / 7.9 / 7.9 / 10.9 | 0 |
| W (ℓν) + Jets | 54.8 / 51.8 / 43.7 / 44.6 | 0.4 / 0.3 / 0.4 / 1.1 |
| QCD Multijet | 0.0 / 2.8 / 3.5 / 2.6 | 0.4 / 1.5 / 1.6 / 2.1 |
| G + Jets | 0 | 99.2 / 98.1 / 97.9 / 96.4 |
| other (VV, Z+jets) | 1.4 / 1.9 / 1.6 / 1.9 | < 0.2 |

TT + ST makes up 44–51% of the tt̄ CR MC, the same as in 2023preBPix (45–53%). The γ CR is ≥ 96% G + Jets.

Jet-class yields in the tt̄ CR (μ+e), all background MC:

| pT [GeV] | Top pre-tag | Top tag | B pre-tag | B tag | QCD pre-tag | QCD tag | MC tag | Data tag |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 250–325 | 58.2 | 39.5 | 40.1 | 14.9 | 236.7 | 9.5 | 63.8 | 44 |
| 325–400 | 166.1 | 123.2 | 78.3 | 30.1 | 509.7 | 17.0 | 170.2 | 81 |
| 400–600 | 878.3 | 700.2 | 300.7 | 117.8 | 1657.7 | 51.6 | 869.6 | 519 |
| 600–2000 | 252.0 | 200.9 | 150.6 | 50.6 | 690.4 | 20.3 | 271.9 | 226 |

- Top-jets are 17–31% of the pre-tag MC but 62–81% of the tagged MC.
- The subtracted tagged terms are 20–38% of the tagged MC (B-jet 14–23%, QCD-jet 6–15%). They are largest in the lowest pT bin.
- The subtracted pre-tag term α·(N(B) + N(QCD)) is 83 / 78 / 69 / 77% of N(Data). The eff_Data denominators that remain are **39 / 101 / 643 / 237** events. The two lowest bins are therefore fragile.

## 3. Step 1 — α = N(Data) / N(all bkg MC)

| Region | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| γ CR (gcr) | 0.385 ± 0.019 | 0.327 ± 0.009 | 0.556 ± 0.006 | 0.959 ± 0.021 |
| tt̄(μ) CR (tmcr) | 0.572 ± 0.055 | 0.631 ± 0.040 | 0.729 ± 0.023 | 0.936 ± 0.041 |
| tt̄(e) CR (tecr) | 0.807 ± 0.082 | 0.572 ± 0.046 | 0.736 ± 0.027 | 0.943 ± 0.046 |
| **tt̄ CR (μ+e)** | **0.666 ± 0.047** | **0.609 ± 0.030** | **0.732 ± 0.017** | **0.939 ± 0.031** |

N(Data) / N(MC), tt̄ CR μ+e: 223/335.0, 459/754.0, 2076/2836.7, 1026/1092.9; γ CR: 411/1068, 1538/4702.5, 8487/15263, 2110/2200.6.

- In both CRs α rises with pT, from about 0.35 (γ) or 0.6 (tt̄) to about 0.95 at 600–2000 GeV. In 2023preBPix it was flat: γ 0.47–0.57, tt̄ 0.82–0.94.
- μ and e agree above 400 GeV. In 2023preBPix the electron CR was about 0.1 lower in every bin; that offset has disappeared here. At 250–325 GeV the electron CR is 2.4σ higher (0.81 vs 0.57), driven by a low `tmcr` data count (Table below).

**Where the pT dependence comes from.** Yield per fb⁻¹, 2023postBPix / 2023preBPix (lumi ratio 9.68 / 17.96 = 0.54):

| | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| γ CR, data | 0.76 | 0.69 | 0.98 | 1.11 |
| γ CR, MC (all bkg) | 0.94 | 0.99 | 1.00 | **0.63** |
| γ CR, G + Jets MC | 0.99 | 0.99 | 1.00 | **0.62** |
| tt̄ CR, data | 0.72 | 0.75 | 0.87 | 1.05 |
| tt̄ CR, MC | 1.01 | 1.01 | 1.03 | 0.98 |
| tt̄(μ) CR, data | **0.58** | 0.74 | 0.83 | 0.94 |
| tt̄(e) CR, data | 0.97 | 0.75 | 0.95 | 1.22 |

- **The tt̄ CR MC scales with luminosity in every bin, but the data falls short below 600 GeV.** The deficit is 25–28% below 400 GeV in both CRs and is largest in `tmcr` at 250–325 GeV (MET data, 115 events). Because it depends on pT and is seen in both lepton channels and in the γ CR, a flat luminosity offset alone is not a sufficient explanation. Things to check: the 2023postBPix lumi mask and luminosity value (9.68 fb⁻¹), completeness of the MET and EGamma data files, and whether the BPix-hole veto is applied the same way to data and MC.
- **At 600–2000 GeV the γ CR MC is 37% low, and it is G + Jets alone.** The data in that bin scales normally (1.11). The tt̄ CR MC in the same bin is normal (0.98). This points to a missing or partly processed G + Jets sample covering high photon or jet pT. γ CR α = 0.96 in that bin is therefore an artifact, and the γ CR α would otherwise be about 0.5–0.6 at all pT, as in 2023preBPix.

Figures: `plots_effSF_toptagTEST3/alpha/2023postBPix/alpha_gcr_zoom.pdf`, `alpha_ttcr_split_zoom.pdf`.

## 4. Step 2 — MC efficiency (eff_MC)

### Top-jet tag efficiency (tt̄ CR, all bkg)

| Region | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| tmcr | 0.664 ± 0.019 | 0.744 ± 0.011 | 0.799 ± 0.004 | 0.799 ± 0.009 |
| tecr | 0.700 ± 0.024 | 0.737 ± 0.013 | 0.795 ± 0.005 | 0.795 ± 0.010 |
| **μ+e** | **0.678 ± 0.015** | **0.742 ± 0.008** | **0.797 ± 0.003** | **0.797 ± 0.006** |

### QCD-jet mis-tag rate (γ CR, G + Jets)

| 250–325 | 325–400 | 400–600 | 600–2000 |
|---:|---:|---:|---:|
| 0.0310 ± 0.0018 | 0.0209 ± 0.0007 | 0.0206 ± 0.0003 | 0.0189 ± 0.0004 |

- The Top-jet efficiency rises from 0.68 to 0.80 with pT and agrees with 2023preBPix (0.69–0.81) within 1σ, except at 600–2000 GeV (0.797 vs 0.811, 1.9σ). μ and e agree.
- The QCD-jet mis-tag rate is 1.9–3.1%. Above 325 GeV it matches 2023preBPix (1.9–2.2%). At 250–325 GeV it is 2.0σ higher (0.031 vs 0.027).
- The missing high-pT G + Jets MC (Section 3) does not visibly change the 600–2000 GeV mis-tag rate (0.0189 vs 0.0186).

Figures: `plots_effSF_toptagTEST3/eff/2023postBPix/eff_toptag_combined.pdf`, `eff_qcdjet_gcr_zoom.pdf`.

## 5. Step 3 — Data efficiency and SF

### 5.1 QCD-jet (γ CR, B-jet subtracted)

| pT [GeV] | α | N(Data) | N(Data, tag) | eff_MC | eff_Data | SF_QCD |
|---|---:|---:|---:|---:|---:|---:|
| 250–325 | 0.385 | 411 | 19 | 0.0310 | 0.0352 ± 0.0144 | **1.13 ± 0.47** |
| 325–400 | 0.327 | 1538 | 103 | 0.0209 | 0.0596 ± 0.0098 | **2.85 ± 0.48** |
| 400–600 | 0.556 | 8487 | 420 | 0.0206 | 0.0436 ± 0.0059 | **2.11 ± 0.29** |
| 600–2000 | 0.959 | 2110 | 103 | 0.0189 | 0.0433 ± 0.0072 | **2.29 ± 0.38** |

Tagged data / tagged MC in the γ CR = 0.38, 0.75, 1.00, 1.75 (2023preBPix: 0.96, 0.99, 1.16, 1.20).

- As in 2023preBPix, SF_QCD ≈ 2–3 comes mainly from the pre-tag normalization. In the three upper bins the data mis-tag rate (4.3–6.0%) is close to 2023preBPix (4.9–5.9%), so SF_QCD agrees with it (2.64, 2.54, 2.77) within 1σ. This holds even at 600–2000 GeV, where α and the tagged ratio both move because of the G + Jets deficit: the SF uses rates, and the G + Jets loss there does not change the MC rate.
- **The 250–325 GeV bin is low (1.13).** It has 19 tagged data events. Scaling 2023preBPix (70 events in 17.96 fb⁻¹) to 9.68 fb⁻¹ predicts about 38, and the pre-tag deficit (0.76) lowers that to about 29. 19 is therefore about a 2σ downward fluctuation on top of the data deficit. The bin is 1.8–2.6σ below the other bins. Its 41% error makes it effectively unconstrained.

### 5.2 Top-jet (tt̄ CR, μ+e)

| pT [GeV] | α | N(Data) | N(Data, tag) | eff_MC | eff_Data raw | SF_Top raw | eff_Data SF_QCD-corr. | SF_Top SF_QCD-corr. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 250–325 | 0.666 | 223 | 44 | 0.678 | 0.718 ± 0.248 | 1.06 ± 0.37 | 0.696 ± 0.259 | **1.03 ± 0.38** |
| 325–400 | 0.609 | 459 | 81 | 0.742 | 0.518 ± 0.142 | 0.70 ± 0.19 | 0.328 ± 0.146 | **0.44 ± 0.20** |
| 400–600 | 0.732 | 2076 | 519 | 0.797 | 0.615 ± 0.081 | 0.77 ± 0.10 | 0.549 ± 0.082 | **0.69 ± 0.10** |
| 600–2000 | 0.939 | 1026 | 226 | 0.797 | 0.674 ± 0.127 | 0.85 ± 0.16 | 0.569 ± 0.130 | **0.71 ± 0.16** |

Lepton-flavour split (SF_Top raw / SF_QCD-corrected):

| pT [GeV] | tmcr raw | tmcr corr. | tecr raw | tecr corr. |
|---|---:|---:|---:|---:|
| 250–325 | 1.30 ± 0.51 | 1.27 ± 0.52 | 0.80 ± 0.45 | 0.77 ± 0.46 |
| 325–400 | 0.71 ± 0.23 | 0.46 ± 0.23 | 0.69 ± 0.26 | 0.43 ± 0.26 |
| 400–600 | 0.73 ± 0.11 | 0.64 ± 0.11 | 0.83 ± 0.12 | 0.75 ± 0.12 |
| 600–2000 | 0.81 ± 0.18 | 0.66 ± 0.19 | 0.89 ± 0.19 | 0.78 ± 0.19 |

- μ and e agree within errors in every bin.
- The SF_QCD correction lowers SF_Top by 0.03 / 0.26 / 0.08 / 0.13. The 325–400 GeV shift is the largest because SF_QCD is largest there (2.85) and the denominator is small (101 events). The 250–325 GeV shift is small because SF_QCD there is only 1.13.
- Tagged data / tagged MC = 0.69, 0.48, 0.60, 0.83 vs α = 0.67, 0.61, 0.73, 0.94. As in every other era, SF_Top < 1 because tagged data falls further below MC than the pre-tag data does. At 250–325 GeV the two ratios are equal, so SF_Top ≈ 1.
- **400–2000 GeV agrees with 2023preBPix** (corrected 0.72 ± 0.09, 0.74 ± 0.14; raw 0.83, 0.93). The 325–400 GeV corrected value (0.44 ± 0.20) is 0.8σ below 2023preBPix (0.66 ± 0.17), using the combined error. It should not be over-interpreted while the data deficit in that bin is unexplained.

### 5.3 MC outlier scan

Unlike in 2022postEE, no process has a few high-weight events dominating a bin. The largest effective weights (Σw² / Σw) are:

| Region | Process | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---|---:|---:|---:|---:|
| tmcr | QCD Multijet | 0.1 ± 0.1 | 21.1 ± 9.4 | 99.5 ± 18.3 | 26.6 ± 3.1 |
| gcr | QCD Multijet | 4.2 ± 4.2 | 71.7 ± 17.4 | 237.6 ± 30.1 | 46.9 ± 4.5 |

Weights are ≈ 4 below 600 GeV (about 5 events in `tmcr` at 325–400 GeV). Every other process has weights ≤ 0.5. In the Top-jet eff_Data, the QCD-jet pre-tag MC term is at most 0.086 (250–325 GeV, raw) and 0.028–0.044 at 325–400 GeV (Section 6.3). That is below the data statistics in every bin.

Figures: `plots_effSF_toptagTEST3/effdata/2023postBPix/effdata_qcdjet_gcr.pdf`, `effdata_topjet_ttcr_combined_zoom.pdf`, `effdata_topjet_tmcr_zoom.pdf`, `effdata_topjet_tecr_zoom.pdf`.

## 6. Uncertainties

The horizontal bars in all plots are the bin half-widths, not uncertainties. Histogram variances are Σw² for MC and N (Poisson) for data.

### 6.1 How each error is computed

Same method as the other TEST3 reports:

| Quantity | Error formula | Where |
|---|---|---|
| α = D / M | σ_α/α = √(var_D/D² + var_M/M²); data and MC independent, regions summed before the ratio | `alpha()`, Step 1 |
| eff_MC = P / (P + F) | σ² = [(1−ε)²·var_P + ε²·var_F] / N², the weighted-event binomial error (pass and fail disjoint) | `efficiency()`, Step 2 |
| eff_Data | linearized propagation: shift each input by +1σ with the others fixed, add the shifts of eff_Data in quadrature (inputs uncorrelated) | `propagate()`, Step 3 |
| SF = eff_Data / eff_MC | σ_SF/SF = √((σ_Data/eff_Data)² + (σ_MC/eff_MC)²), the two efficiencies treated as independent | `measure()`, Step 3 |

α is not a separate input to eff_Data: it is recomputed inside the formula as N(Data)/N(MC), so its correlation with N(Data) is handled exactly.

| Input | σ | Type |
|---|---|---|
| N(Data, tag), N(Data) | √N | data stat |
| N(MC) (all bkg, for α) | √Σw² | MC stat |
| B-jet MC, tagged and pre-tag | √Σw²; γ CR adds a 50% normalization uncertainty (AN "BJetNorm") | MC stat (+ syst in γ CR) |
| QCD-jet MC, tagged and pre-tag (Top-jet only) | √Σw² | MC stat |
| SF_B | 1 ± 0.5 (no dedicated B-jet CR) | syst |
| SF_QCD (Top-jet, corrected variant) | γ CR value ± its total error, per bin | syst |

### 6.2 Breakdown — QCD-jet eff_Data (γ CR)

Contribution of each input to σ(eff_Data) (+1σ shift):

| Source | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| N(Data, tag) | 0.0110 | 0.0068 | 0.0025 | 0.0050 |
| N(Data) | 0.0022 | 0.0017 | 0.0005 | 0.0011 |
| N(MC) (α) | 0.0001 | < 0.0001 | < 0.0001 | < 0.0001 |
| B-jet tagged: 50% norm. + MC stat | 0.0065 | 0.0048 | 0.0038 | 0.0036 |
|   of which MC stat | 0.0012 | 0.0005 | 0.0002 | 0.0003 |
| B-jet pre-tag | 0.0006 | 0.0010 | 0.0007 | 0.0007 |
| SF_B = 1 ± 0.5 | 0.0063 | 0.0048 | 0.0038 | 0.0035 |
| **Total** | **0.0144** | **0.0098** | **0.0059** | **0.0072** |
| stat / syst | 0.0113 / 0.0090 | 0.0070 / 0.0068 | 0.0026 / 0.0054 | 0.0051 / 0.0051 |

Relative σ(SF_QCD): from eff_Data 41 / 16 / 14 / 17%, from eff_MC 5.7 / 3.5 / 1.5 / 2.1%.

- Data statistics dominate at 250–325 GeV and equal the systematics at 325–400 and 600–2000 GeV. The B-jet 50% normalization and SF_B ± 0.5 dominate only at 400–600 GeV. With 9.68 fb⁻¹, data statistics matter much more than in 2022postEE.
- The total does not cover the normalization questions of Section 3.

### 6.3 Breakdown — Top-jet eff_Data (tt̄ CR, μ+e, SF_QCD-corrected)

| Source | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| N(Data, tag) | 0.171 | 0.089 | 0.035 | 0.064 |
| N(Data) | 0.071 | 0.036 | 0.017 | 0.029 |
| N(MC) (α) | 0.051 | 0.011 | 0.008 | 0.014 |
| B-jet tagged | 0.015 | 0.007 | 0.006 | 0.006 |
| B-jet pre-tag | 0.018 | 0.005 | 0.005 | 0.006 |
| QCD-jet tagged | 0.015 | 0.020 | 0.005 | 0.017 |
| QCD-jet pre-tag | 0.083 | 0.028 | 0.015 | 0.021 |
| SF_B = 1 ± 0.5 | 0.128 | 0.091 | 0.067 | 0.101 |
| SF_QCD (γ CR) | 0.076 | 0.049 | 0.017 | 0.031 |
| **Total** | **0.259** | **0.146** | **0.082** | **0.130** |

Grouped:

| Group | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| Data stat | 0.185 | 0.096 | 0.040 | 0.070 |
| MC stat | 0.102 | 0.037 | 0.019 | 0.031 |
| SF_B | 0.128 | 0.091 | 0.067 | 0.101 |
| SF_QCD | 0.076 | 0.049 | 0.017 | 0.031 |

Relative σ(SF_Top): from eff_Data 37 / 44 / 15 / 23%, from eff_MC 2.2 / 1.1 / 0.4 / 0.8%.

- **Data statistics are the largest source below 400 GeV**, and SF_B ± 0.5 above. At 325–400 GeV the two are about equal (0.096 vs 0.091).
- MC statistics are 0.10 at 250–325 GeV. This comes mostly from the QCD-jet pre-tag W + Jets MC, which enters a denominator of only 39 events. Above 325 GeV MC statistics are ≤ 0.037.
- The SF_QCD term (0.017–0.076) is the third-largest source. It is larger here than in 2022postEE because SF_QCD has 17–41% errors in this era.
- In the raw variant the SF_QCD term is absent; the totals are 0.248 / 0.142 / 0.081 / 0.127, with groups: data stat 0.185 / 0.096 / 0.040 / 0.070, MC stat 0.105 / 0.052 / 0.021 / 0.033, SF_B 0.128 / 0.091 / 0.067 / 0.101.
- eff_MC contributes ≤ 2.2% to σ(SF_Top), so the SF error is almost entirely the eff_Data error.

### 6.4 Validity checks

**Linearity.** Repeating the propagation with −1σ shifts changes every total by ≤ 0.5% (Top-jet 250–325 GeV: 0.259 → 0.260 corrected, 0.248 → 0.250 raw). The linear error is valid in every bin. The 2022postEE outlier problem does not occur here.

**Correlation of tagged and total data.** N(Data, tag) is a subset of N(Data), but they are propagated as independent. Re-parametrizing with the independent counts N(Data, tag) and N(Data, fail) = N(Data) − N(Data, tag):

| | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| Top-jet corr.: data stat, current | 0.185 | 0.096 | 0.040 | 0.070 |
| Top-jet corr.: data stat, (tag, fail) | 0.148 | 0.079 | 0.030 | 0.055 |
| Top-jet corr.: total, current | 0.259 | 0.146 | 0.082 | 0.130 |
| Top-jet corr.: total, (tag, fail) | 0.233 | 0.135 | 0.078 | 0.123 |
| Top-jet raw: total, current | 0.248 | 0.142 | 0.081 | 0.127 |
| Top-jet raw: total, (tag, fail) | 0.222 | 0.131 | 0.077 | 0.119 |
| QCD-jet: total, current | 0.0144 | 0.0098 | 0.0059 | 0.0072 |
| QCD-jet: total, (tag, fail) | 0.0139 | 0.0095 | 0.0059 | 0.0070 |

The data-stat part shrinks by 18–23%. Because data statistics matter more in this era, the Top-jet total shrinks by 5–10%, against 2–4% in 2022postEE. The QCD-jet total shrinks by ≤ 3%. The current numbers, and the JSON, are therefore conservative. The recommendation to switch to the (tag, fail) parametrization stands.

**eff_MC–eff_Data correlation.** N(MC) in α includes the target class, so eff_MC and eff_Data share MC events. With eff_MC contributing ≤ 5.7% (QCD) and ≤ 2.2% (Top) of σ(SF), neglecting this is safe.

**Not included:** JES/JER, parton shower, top-pT reweighting, scale/PDF, pileup, lepton/photon ID and trigger SFs, the WP choice, and the normalization questions of Section 3.

## 7. Correctionlib file

`plots_effSF_toptagTEST3/sf_json/toptag_sf_2023postBPix.json`, same layout as the other eras (`studies/sf_json.py`). It loads with `correctionlib.CorrectionSet.from_file()`, and the evaluated values match the tables above.

| Correction | Inputs | Content (250–325 / 325–400 / 400–600 / 600–2000 GeV) |
|---|---|---|
| `toptag_sf`, `top` | `pt`, `jetclass`, `systematic` | 1.027 / 0.443 / 0.689 / 0.714 (tt̄ μ+e, SF_QCD-corrected); up/down = ±1σ total |
| `toptag_sf`, `qcd` | | 1.133 / 2.852 / 2.113 / 2.294 (γ CR); up/down = ±1σ total |
| `toptag_sf`, `b` | | 1 ± 0.5 |
| `toptag_efficiency` | `pt`, `jetclass` | Top 0.678 / 0.742 / 0.797 / 0.797; QCD 0.0310 / 0.0209 / 0.0206 / 0.0189; B 0.370 / 0.384 / 0.392 / 0.336 |

pT edges 250, 325, 400, 600, 2000 with `flow="clamp"`. No bin needed the non-finite fallback (1 ± 1). The file contains every bin as measured, with no merging or capping. Before it is used for event weights, the 250–325 GeV SF_QCD and 325–400 GeV SF_Top in particular should be reviewed against the open points below.

## 8. Conclusions and next steps

1. **Data normalization below 400 GeV:** the data per fb⁻¹ is 25–28% below 2023preBPix in both CRs while the MC scales normally (`tmcr` 250–325 GeV: 42% below). Check the postBPix lumi / lumi mask, data-file completeness (MET and EGamma), and the BPix-hole veto in data vs MC. This is the same kind of effect as in 2022postEE.
2. **G + Jets MC at 600–2000 GeV:** 37% below the expected per-fb⁻¹ yield. Check that every G + Jets sample (especially the high-HT / high-pT ones) was processed for 2023postBPix. SF_QCD in that bin is not affected much (it is a rate), but the γ CR α is.
3. **SF_QCD:** 2.1–2.9 above 325 GeV, consistent with 2023preBPix. The 250–325 GeV value (1.13 ± 0.47) has 19 tagged events. Consider combining it with 2023preBPix, or using a 2023-combined SF_QCD in that bin.
4. **SF_Top:** 0.69 / 0.71 (corrected) at 400–2000 GeV, consistent with 2023preBPix. The two lower bins are dominated by data statistics, and the subtraction removes about 80% of the pre-tag data there. Combining with 2023preBPix or merging 250–400 GeV would make them usable. The raw vs corrected choice (0.03–0.26 here) still needs a decision.
5. **Uncertainties:** switching to (tagged, failed) data counts lowers the Top-jet errors by 5–10% in this era. SF_B ± 0.5 is the dominant Top-jet source above 400 GeV.
6. Still missing: the systematics listed in 6.4, and a decision on whether to combine pre/postBPix (and pre/postEE) before deriving the final SFs.
