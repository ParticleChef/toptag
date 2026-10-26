# Top-tag α, MC efficiency, Data efficiency, SF and uncertainties — 2022preEE (TEST3)

**Input:** `hadmonotop2022_toptagTEST3.scaled` (updated 2026-09-30 17:07)
**Report date:** 2026-09-30
**Tagger / WP:** particleTransformer `TvsQCD` > 0.33, leading AK15 jet
**Luminosity:** 2022preEE 7.99 fb⁻¹ (13.6 TeV)
**PDF version (with figures):** `reports/report_toptag_effSF_2022preEE_TEST3_20260930.pdf`
**Companion:** `report_toptag_effSF_2023preBPix_TEST3_20260930.md` (same method, 2023preBPix)

Definitions follow the AN (AN-2020/061 Sec. 6.6): α and the Top-jet class yields use **all background MC**.

## Summary

- **The MC normalization is good in the updated 2022 file.** α = 0.90–1.05 in the γ CR and the combined tt̄ CR in every pT bin (0.87–1.12 in the single-lepton CRs). In the tt̄ CR, TT + ST is 44–54% of the MC.
- **Data statistics are now sufficient** for eff_Data / SF: 302–2392 events per bin in the tt̄ CR and 970–13672 in the γ CR. The earlier note that the 2022 TEST3 statistics were too low no longer applies.
- **SF_QCD (γ CR):** 1.32–1.63 (±16–21%). Unlike 2023 (SF_QCD ≈ 2.5, driven by α ≈ 0.5), the 2022 SF_QCD comes from a genuine excess of tagged data (tagged data / tagged MC = 1.17–1.35 with α ≈ 1).
- **SF_Top (tt̄ CR, μ+e):** 0.74–1.05 with raw subtraction, 0.69–0.87 with SF_QCD-corrected subtraction. The 400–600 GeV bin is the most precise (0.84 ± 0.10).
- **Uncertainties:** SF_Top ±12–38% per bin. The two largest sources are data statistics of the tagged yield and **SF_B = 1 ± 0.5**. The SF_B term alone is 0.065–0.135 on eff_Data, larger than data statistics in the 400–600 GeV bin.
- **The error bars are slightly conservative.** The propagation treats N(Data, tag) and N(Data) as independent although one is a subset of the other. Using independent (tagged, failed) counts instead lowers the Top-jet eff_Data error by 6–10% and the QCD-jet error by ≤ 4%.

## 1. Steps

All steps run with the `plot_env` interpreter (Python 3.10) and now take the `.scaled` file with `-f`. Outputs go to `plots_effSF_<version>/<kind>/<era>/`, with the version and era taken from the file name.

```bash
cd /Users/hongjieun/working/toptag
PY=/Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python
F=hadmonotop2022_toptagTEST3
$PY studies/plot_alpha_toptagTEST.py           -f $F   # Step 1: α
$PY studies/plot_toptag_eff_toptagTEST.py      -f $F   # Step 2: eff_MC
$PY studies/plot_toptag_eff_data_toptagTEST.py -f $F   # Step 3: eff_Data + SF
```

Full console output: `plots_effSF_toptagTEST3/run_log_2022preEE.txt`. The per-source uncertainty breakdown (Section 6) was computed with a one-off script that re-uses `measure()` from Step 3 and records each input's shift.

| Step | Script | Reads | Computes | Output |
|---|---|---|---|---|
| 1 | `plot_alpha_toptagTEST.py` | `fjpt_bin` (bkg + data) | α per region and pT bin: γ CR, tt̄ CR μ+e, μ/e split | `alpha/2022preEE/` |
| 2 | `plot_toptag_eff_toptagTEST.py` | `fjpt_cat` (bkg) | eff_MC: Top-jet (tmcr, tecr, μ+e, all bkg), QCD-jet (gcr, G + Jets) | `eff/2022preEE/` |
| 3 | `plot_toptag_eff_data_toptagTEST.py` | both | eff_Data and SF: QCD-jet first, then Top-jet (raw and SF_QCD-corrected) | `effdata/2022preEE/` |

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
| TT | 36.4 / 35.6 / 45.7 / 40.8 | 0 |
| ST | 9.1 / 8.5 / 8.2 / 11.1 | 0 |
| W (ℓν) + Jets | 51.3 / 52.9 / 41.9 / 43.4 | 0.2 / 0.2 / 0.3 / 0.5 |
| QCD Multijet | 1.7 / 1.4 / 2.4 / 3.0 | 0.6 / 1.2 / 1.7 / 1.3 |
| G + Jets | 0 | 99.2 / 98.6 / 97.9 / 98.0 |
| other (VV, Z+jets) | 1.5 / 1.6 / 1.7 / 1.8 | < 0.1 |

Jet-class yields in the tt̄ CR (μ+e), all background MC:

| pT [GeV] | Top pre-tag | Top tag | B pre-tag | B tag | QCD pre-tag | QCD tag | MC tag | Data tag |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 250–325 | 55.6 | 40.0 | 38.0 | 15.0 | 216.9 | 11.8 | 66.8 | 67 |
| 325–400 | 160.7 | 124.2 | 74.8 | 30.6 | 470.6 | 22.3 | 177.1 | 176 |
| 400–600 | 803.9 | 654.6 | 259.8 | 104.9 | 1395.0 | 64.2 | 823.6 | 726 |
| 600–2000 | 212.5 | 177.0 | 132.6 | 48.6 | 588.2 | 26.6 | 252.2 | 186 |

- Top-jets are 18–33% of the pre-tag MC but 60–80% of the tagged MC.
- The subtracted tagged terms are 20–40% of the tagged MC (B-jet 13–23%, QCD-jet 8–18%); they are largest in the lowest pT bin.
- The composition is nearly identical to 2023preBPix TEST3.

## 3. Step 1 — α = N(Data) / N(all bkg MC)

| Region | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| γ CR (gcr) | 0.936 ± 0.031 | 1.035 ± 0.016 | 0.949 ± 0.008 | 0.962 ± 0.017 |
| tt̄(μ) CR (tmcr) | 0.869 ± 0.071 | 1.014 ± 0.052 | 0.932 ± 0.026 | 0.883 ± 0.042 |
| tt̄(e) CR (tecr) | 1.124 ± 0.099 | 1.112 ± 0.065 | 1.032 ± 0.033 | 0.921 ± 0.049 |
| **tt̄ CR (μ+e)** | **0.973 ± 0.058** | **1.054 ± 0.041** | **0.973 ± 0.021** | **0.899 ± 0.032** |

N(Data) / N(MC), tt̄ CR μ+e: 302/310.5, 744/706.0, 2392/2458.7, 839/933.3; γ CR: 970/1036.2, 4907/4739.1, 13672/14414, 3141/3264.1.

- α is within 10% of 1 in the γ CR and the combined tt̄ CR; the single-lepton CRs spread from 0.87 to 1.12. The γ CR does **not** show the 2× G + Jets excess seen in 2023preBPix (α ≈ 0.5).
- The electron tt̄ CR is 0.04–0.25 higher than the muon one (2.4σ at 400–600 GeV). In 2023 the electron CR was lower instead, so this is not a consistent electron-side offset.
- α falls slightly in the highest bin in both CRs (0.90 in tt̄).

Figures: `plots_effSF_toptagTEST3/alpha/2022preEE/alpha_gcr_zoom.pdf`, `alpha_ttcr_split_zoom.pdf`.

## 4. Step 2 — MC efficiency (eff_MC)

### Top-jet tag efficiency (tt̄ CR, all bkg)

| Region | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| tmcr | 0.720 ± 0.012 | 0.773 ± 0.007 | 0.817 ± 0.003 | 0.835 ± 0.005 |
| tecr | 0.716 ± 0.014 | 0.773 ± 0.008 | 0.810 ± 0.003 | 0.829 ± 0.006 |
| **μ+e** | **0.719 ± 0.009** | **0.773 ± 0.005** | **0.814 ± 0.002** | **0.833 ± 0.004** |

### QCD-jet mis-tag rate (γ CR, G + Jets)

| 250–325 | 325–400 | 400–600 | 600–2000 |
|---:|---:|---:|---:|
| 0.0392 ± 0.0019 | 0.0284 ± 0.0008 | 0.0273 ± 0.0003 | 0.0262 ± 0.0003 |

- The Top-jet efficiency rises from 0.72 to 0.83 with pT, 0.02–0.03 higher than in 2023preBPix (0.69–0.81). μ and e agree.
- The QCD-jet mis-tag rate is 2.6–3.9%, higher than in 2023 (1.9–2.7%).

Figures: `plots_effSF_toptagTEST3/eff/2022preEE/eff_toptag_combined.pdf`, `eff_qcdjet_gcr_zoom.pdf`.

## 5. Step 3 — Data efficiency and SF

### QCD-jet (γ CR, B-jet subtracted)

| pT [GeV] | α | N(Data) | N(Data, tag) | eff_MC | eff_Data | SF_QCD |
|---|---:|---:|---:|---:|---:|---:|
| 250–325 | 0.936 | 970 | 72 | 0.0392 | 0.0639 ± 0.0133 | **1.63 ± 0.35** |
| 325–400 | 1.035 | 4907 | 225 | 0.0284 | 0.0376 ± 0.0077 | **1.32 ± 0.27** |
| 400–600 | 0.949 | 13672 | 630 | 0.0273 | 0.0393 ± 0.0062 | **1.44 ± 0.23** |
| 600–2000 | 0.962 | 3141 | 129 | 0.0262 | 0.0347 ± 0.0067 | **1.33 ± 0.26** |

Tagged data / tagged MC in the γ CR = 1.35, 1.27, 1.24, 1.17 at α ≈ 1. The data mis-tag rate is therefore genuinely ~30–60% above MC. In 2023 the tagged yields agreed and SF_QCD ≈ 2.5 came from the pre-tag normalization instead, so the two eras disagree in both the SF value and its origin.

### Top-jet (tt̄ CR, μ+e)

| pT [GeV] | α | N(Data) | N(Data, tag) | eff_MC | eff_Data raw | SF_Top raw | eff_Data SF_QCD-corr. | SF_Top SF_QCD-corr. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 250–325 | 0.973 | 302 | 67 | 0.719 | 0.756 ± 0.232 | **1.05 ± 0.32** | 0.623 ± 0.239 | **0.87 ± 0.33** |
| 325–400 | 1.054 | 744 | 176 | 0.773 | 0.710 ± 0.137 | **0.92 ± 0.18** | 0.665 ± 0.142 | **0.86 ± 0.18** |
| 400–600 | 0.973 | 2392 | 726 | 0.814 | 0.718 ± 0.078 | **0.88 ± 0.10** | 0.683 ± 0.080 | **0.84 ± 0.10** |
| 600–2000 | 0.899 | 839 | 186 | 0.833 | 0.619 ± 0.142 | **0.74 ± 0.17** | 0.578 ± 0.145 | **0.69 ± 0.17** |

Lepton-flavour split (SF_Top raw / SF_QCD-corrected):

| pT [GeV] | tmcr raw | tmcr corr. | tecr raw | tecr corr. |
|---|---:|---:|---:|---:|
| 250–325 | 0.97 ± 0.40 | 0.78 ± 0.40 | 1.13 ± 0.44 | 0.94 ± 0.44 |
| 325–400 | 0.93 ± 0.21 | 0.86 ± 0.22 | 0.91 ± 0.22 | 0.86 ± 0.23 |
| 400–600 | 0.91 ± 0.11 | 0.86 ± 0.11 | 0.85 ± 0.11 | 0.81 ± 0.11 |
| 600–2000 | 0.84 ± 0.19 | 0.78 ± 0.20 | 0.63 ± 0.20 | 0.59 ± 0.20 |

- μ and e agree within errors in every bin.
- Because SF_QCD ≈ 1.4 (not 2.5 as in 2023), the SF_QCD correction lowers SF_Top by only 0.04–0.19, well below 1σ.
- Tagged data / tagged MC = 1.00, 0.99, 0.88, 0.74 vs α = 0.97, 1.05, 0.97, 0.90. SF_Top < 1 above 400 GeV comes from this tagged-data deficit, as in 2023.
- SF_Top agrees with 2023preBPix within errors (2023 corrected: 0.66–0.78).

Figures: `plots_effSF_toptagTEST3/effdata/2022preEE/effdata_qcdjet_gcr.pdf`, `effdata_topjet_ttcr_combined_zoom.pdf`, `effdata_topjet_tmcr_zoom.pdf`, `effdata_topjet_tecr_zoom.pdf`.

## 6. Uncertainties

The horizontal bars in all plots are the bin half-widths, not uncertainties. Histogram variances are Σw² for MC and N (Poisson) for data.

### 6.1 How each error is computed

| Quantity | Error formula | Where |
|---|---|---|
| α = D / M | σ_α/α = √(var_D/D² + var_M/M²); data and MC independent, regions summed before the ratio | `alpha()`, Step 1 |
| eff_MC = P / (P + F) | σ² = [(1−ε)²·var_P + ε²·var_F] / N², the weighted-event binomial error (pass and fail disjoint) | `efficiency()`, Step 2 |
| eff_Data | linearized propagation: shift each input by +1σ with the others fixed, add the shifts of eff_Data in quadrature (inputs uncorrelated) | `propagate()`, Step 3 |
| SF = eff_Data / eff_MC | σ_SF/SF = √((σ_Data/eff_Data)² + (σ_MC/eff_MC)²), the two efficiencies treated as independent | `measure()`, Step 3 |

α is not a separate input to eff_Data: it is recomputed inside the formula as N(Data)/N(MC), so its correlation with N(Data) is handled exactly.

Inputs to the eff_Data propagation:

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
| N(Data, tag) | 0.0091 | 0.0032 | 0.0019 | 0.0037 |
| N(Data) | 0.0024 | 0.0007 | 0.0004 | 0.0007 |
| N(MC) (α) | 0.0001 | < 0.0001 | < 0.0001 | < 0.0001 |
| B-jet tagged: 50% norm. + MC stat | 0.0067 | 0.0049 | 0.0042 | 0.0039 |
|   of which MC stat | 0.0010 | 0.0004 | 0.0002 | 0.0002 |
| B-jet pre-tag | 0.0012 | 0.0007 | 0.0006 | 0.0006 |
| SF_B = 1 ± 0.5 | 0.0066 | 0.0049 | 0.0042 | 0.0039 |
| **Total** | **0.0133** | **0.0077** | **0.0062** | **0.0067** |
| stat / syst | 0.0094 / 0.0094 | 0.0033 / 0.0070 | 0.0019 / 0.0059 | 0.0038 / 0.0055 |

Relative σ(SF_QCD): from eff_Data 21 / 21 / 16 / 19%, from eff_MC 4.9 / 2.7 / 1.2 / 1.0%.

- Above 325 GeV the systematics dominate: the B-jet 50% normalization and SF_B ± 0.5 contribute about equally.
- Only the lowest bin is statistics-limited (72 tagged data events).

### 6.3 Breakdown — Top-jet eff_Data (tt̄ CR, μ+e, SF_QCD-corrected)

| Source | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| N(Data, tag) | 0.151 | 0.078 | 0.035 | 0.071 |
| N(Data) | 0.067 | 0.037 | 0.019 | 0.033 |
| N(MC) (α) | 0.034 | 0.022 | 0.007 | 0.012 |
| B-jet tagged | 0.012 | 0.006 | 0.002 | 0.009 |
| B-jet pre-tag | 0.016 | 0.008 | 0.003 | 0.007 |
| QCD-jet tagged | 0.023 | 0.011 | 0.005 | 0.010 |
| QCD-jet pre-tag | 0.060 | 0.036 | 0.011 | 0.019 |
| SF_B = 1 ± 0.5 | 0.135 | 0.095 | 0.065 | 0.114 |
| SF_QCD (γ CR) | 0.074 | 0.038 | 0.018 | 0.032 |
| **Total** | **0.239** | **0.142** | **0.080** | **0.145** |

Grouped:

| Group | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| Data stat | 0.166 | 0.087 | 0.039 | 0.078 |
| MC stat | 0.076 | 0.045 | 0.015 | 0.027 |
| SF_B | 0.135 | 0.095 | 0.065 | 0.114 |
| SF_QCD | 0.074 | 0.038 | 0.018 | 0.032 |

Relative σ(SF_Top): from eff_Data 38 / 21 / 12 / 25%, from eff_MC 1.2 / 0.6 / 0.3 / 0.5%.

- **SF_B ± 0.5 is the largest single systematic** and is larger than data statistics at 325–2000 GeV. A B-jet SF measurement, or even a tighter prior, would give the biggest improvement.
- The QCD-jet pre-tag MC statistics (0.060 at 250–325 GeV) are the largest MC-stat term: this class comes largely from W + jets, whose events carry large weights.
- In the raw variant the SF_QCD term is absent; the totals are 0.232 / 0.137 / 0.078 / 0.142.
- eff_MC contributes ≤ 1.2% to σ(SF_Top), so the SF error is entirely the eff_Data error.

### 6.4 Validity checks

**Linearity.** Repeating the propagation with −1σ shifts changes the totals by < 1% (e.g. 0.2320 → 0.2331 for the raw Top-jet 250–325 GeV bin). The symmetric linear error is adequate at the 1σ level. For the 250–325 GeV Top-jet bin (±38%) the true distribution is still likely non-Gaussian.

**Correlation of tagged and total data.** N(Data, tag) is a subset of N(Data), but they are propagated as independent. Re-parametrizing with the independent counts N(Data, tag) and N(Data, fail) = N(Data) − N(Data, tag):

| | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| Top-jet raw: data stat, current | 0.166 | 0.087 | 0.039 | 0.078 |
| Top-jet raw: data stat, (tag, fail) | 0.129 | 0.067 | 0.028 | 0.062 |
| Top-jet raw: total, current | 0.232 | 0.137 | 0.078 | 0.142 |
| Top-jet raw: total, (tag, fail) | 0.208 | 0.126 | 0.073 | 0.133 |
| QCD-jet: total, current | 0.0133 | 0.0077 | 0.0062 | 0.0067 |
| QCD-jet: total, (tag, fail) | 0.0128 | 0.0077 | 0.0062 | 0.0066 |

The current bars overestimate the Top-jet eff_Data error by 6–10% (data-stat part by 21–28%) and the QCD-jet error by ≤ 4%. The (tag, fail) parametrization should become the default.

**eff_MC–eff_Data correlation.** N(MC) in α includes the target class, so eff_MC and eff_Data share MC events. With eff_MC contributing ≤ 5% (QCD) and ≤ 1.2% (Top) of σ(SF), neglecting this is safe.

**Not included:** JES/JER, parton shower, top-pT reweighting, scale/PDF, pileup, lepton/photon ID and trigger SFs, and the WP choice.

## 7. Conclusions and next steps

1. The updated 2022preEE TEST3 file is usable: α ≈ 1 everywhere, and eff_Data / SF are physical in every bin.
2. SF_Top ≈ 0.85–0.9 below 600 GeV and ≈ 0.7 above; the raw vs SF_QCD-corrected choice matters much less than in 2023.
3. SF_QCD = 1.3–1.6 in 2022 vs 2.3–2.8 in 2023, with a different origin (tagged-data excess vs pre-tag normalization). The γ CR must be understood before combining eras.
4. SF_B = 1 ± 0.5 dominates the Top-jet systematic. Constraining it is the most effective improvement.
5. Switch the data-stat propagation to (tagged, failed) counts; the Top-jet bars shrink by 6–10%.
6. Still missing: the systematics listed in 6.4, the other eras (2022postEE, 2023postBPix), and the correctionlib SF file (Stage 3).
