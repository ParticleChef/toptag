# Top-tag α, MC efficiency, Data efficiency and SF — 2023preBPix (TEST3)

**Input:** `hadmonotop2023_toptagTEST3.scaled`
**Report date:** 2026-09-30
**Tagger / WP:** particleTransformer `TvsQCD` > 0.33, leading AK15 jet
**Luminosity:** 2023preBPix 17.96 fb⁻¹ (13.6 TeV)
**PDF version (with figures):** `reports/report_toptag_effSF_2023preBPix_TEST3_20260930.pdf`

This report replaces the revision-2 report (`report_toptag_effSF_alpha_20260928.md`) for 2023preBPix. It uses the two updated definitions:

1. **α uses all background MC** (AN-2020/061 Sec. 6.6.1): N(MC) = total yield of all 10 background processes in the region, all jet classes.
2. **Top-jet efficiency uses all background MC:** eff_MC and the subtracted B/QCD-jet terms in the tt̄ CRs are summed over every background process, not TT + ST only.

## Summary

- **The MC normalization problem is fixed in TEST3.** In the tt̄ CRs, TT + ST is 45–53% of the MC (previously 3–16%) and α = 0.82–0.94 (previously 100–170 with the TT+ST α, or ~20 with the all-bkg α).
- **Top-jet eff_Data is physical in every bin** (0.66–0.75 raw, 0.50–0.60 SF_QCD-corrected), with no over-subtraction.
- **SF_Top (tt̄ CR, μ+e):** 0.83–1.08 with raw subtraction, 0.66–0.78 with SF_QCD-corrected subtraction; ±10–30% per bin.
- **SF_QCD (γ CR):** 2.30–2.77 (±13–19%), consistent with the earlier 2023preBPix result (2.2–2.7).
- **Open point:** the choice between raw and SF_QCD-corrected subtraction moves SF_Top by 0.15–0.3, about 1σ per bin. It needs a decision.
- **Open point:** the γ CR MC is ~2× the data pre-tag (α ≈ 0.5), while the tagged MC ≈ tagged data. SF_QCD ≈ 2.5 comes mostly from this.

## 1. Steps

All steps run on the laptop with the `plot_env` interpreter (Python 3.10; the system Python 3.7 cannot read the `.scaled` pickles). Each script takes the era key as its only argument. The `2023preBPix_TEST3` key (added to `ERAS` in all three scripts) selects the input file, the 17.96 fb⁻¹ label and the output directory `plots_effSF_TEST3/<kind>/2023preBPix/`.

```bash
cd /Users/hongjieun/working/toptag
PY=/Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python
$PY studies/plot_alpha_toptagTEST.py           2023preBPix_TEST3   # Step 1: α
$PY studies/plot_toptag_eff_toptagTEST.py      2023preBPix_TEST3   # Step 2: eff_MC
$PY studies/plot_toptag_eff_data_toptagTEST.py 2023preBPix_TEST3   # Step 3: eff_Data + SF
```

Full console output: `plots_effSF_TEST3/run_log_2023preBPix.txt`.

**Plot style:** MC `#1f77b4` (legend "MC"), Data `#ff7f0e`, SF black. The Top-jet `effdata_topjet_*` plots show only MC and the SF_QCD-corrected Data; both subtraction variants are printed in the log and tabulated below.

| Step | Script | Reads | Computes | Output |
|---|---|---|---|---|
| 0 | (check) | `fjpt_bin`, `fjpt_cat` | process composition, `fjpt_bin` MC = Σ classes of `fjpt_cat` MC | — |
| 1 | `plot_alpha_toptagTEST.py` | `fjpt_bin` (bkg + data) | α per region and pT bin: γ CR, tt̄ CR μ+e, μ/e split | `alpha/2023preBPix/` |
| 2 | `plot_toptag_eff_toptagTEST.py` | `fjpt_cat` (bkg) | eff_MC: Top-jet (tmcr, tecr, μ+e, all bkg), QCD-jet (gcr, G + Jets) | `eff/2023preBPix/` |
| 3 | `plot_toptag_eff_data_toptagTEST.py` | both | eff_Data and SF: QCD-jet first, then Top-jet (raw and SF_QCD-corrected) | `effdata/2023preBPix/` |

Step 3 must derive SF_QCD before the Top-jet measurement, because the SF_QCD-corrected variant uses it.

### 1.1 Inputs

| Histogram | Axes | Filled for | Used for |
|---|---|---|---|
| `fjpt_bin` | region × pT (250, 325, 400, 600, 2000) × TvsQCD | bkg + data | N(Data), N(Data, tag), N(MC) → α |
| `fjpt_cat` | region × pT (same edges) × TvsQCD × n_qinfj × n_binfj | bkg only | class-split MC yields, eff_MC |

TvsQCD edges [0, 0.33, 1]: bin 1 = tagged. Data per region: `tmcr` → MET, `tecr` and `gcr` → EGamma. The file has all 10 background processes, so one file serves both calibrations. The `fjpt_bin`/`fjpt_cat` consistency check passed in every region (no warning).

Jet classes: **Top-jet** = classes 5–8 (≥ 2 top-decay quarks in the AK15 jet), **B-jet** = 2, 4, **QCD-jet** = 1, 3.

### 1.2 Formulas

- **α** = N(Data) / N(MC), pre-tag (TvsQCD summed), per region and pT bin; N(MC) = all background processes, all jet classes. Regions are summed before the ratio.
- **eff_MC** = N(MC, class, tag) / N(MC, class).
- **QCD-jet (γ CR, G + Jets):**
  eff_Data = [N(Data, tag) − α·SF_B·N(GJ, B, tag)] / [N(Data) − α·N(GJ, B)]
- **Top-jet (tt̄ CR, all bkg):**
  eff_Data = [N(Data, tag) − α·(SF_B·N(B, tag) + SF_QCD·N(QCD, tag))] / [N(Data) − α·(N(B) + N(QCD))]
  - *raw*: SF_QCD = 1; *SF_QCD-corrected*: SF_QCD from the γ CR, per pT bin.
- **SF** = eff_Data / eff_MC.

SF_B and SF_QCD correct the tag rate of the subtracted classes, so they multiply only the tagged terms. Because α normalizes the full MC to data, the Top-jet denominator equals α·N(MC, Top), and SF_Top = [N(Data, tag) − α·(subtracted tagged terms)] / [α·N(MC, Top, tag)].

**Uncertainties:** statistics (data and MC), SF_B = 1 ± 0.5, a 50% normalization uncertainty on the G + Jets B-jet MC (γ CR), and in the corrected variant the SF_QCD uncertainty. Inputs are treated as uncorrelated (one-at-a-time finite-difference propagation). No JES, parton-shower, top-pT or scale/PDF systematics yet.

## 2. Sample composition

Fraction (%) of the pre-tag background MC per pT bin (250–325 / 325–400 / 400–600 / 600–2000 GeV):

| Process | tt̄ CR (μ+e) | γ CR |
|---|---|---|
| TT | 35.5 / 35.9 / 44.2 / 39.7 | 0 |
| ST | 9.5 / 8.6 / 8.1 / 10.8 | 0 |
| W (ℓν) + Jets | 52.2 / 53.0 / 43.1 / 44.4 | 0.3 / 0.3 / 0.3 / 0.7 |
| QCD Multijet | 1.3 / 0.8 / 2.9 / 3.4 | 5.6 / 1.6 / 1.7 / 1.6 |
| G + Jets | 0 | 94.1 / 98.1 / 97.9 / 97.5 |
| other (VV, Z+jets) | 1.6 / 1.5 / 1.6 / 1.7 | < 0.1 |

Jet-class yields in the tt̄ CR (μ+e), all background MC:

| pT [GeV] | Top pre-tag | Top tag | B pre-tag | B tag | QCD pre-tag | QCD tag | MC tag | Data tag |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 250–325 | 110.3 | 76.0 | 70.7 | 27.5 | 433.1 | 17.1 | 120.7 | 119 |
| 325–400 | 310.2 | 232.4 | 149.3 | 56.4 | 927.1 | 32.7 | 321.5 | 244 |
| 400–600 | 1615.8 | 1282.7 | 538.9 | 205.4 | 2955.8 | 92.7 | 1580.8 | 1178 |
| 600–2000 | 473.9 | 384.4 | 277.8 | 90.2 | 1313.7 | 41.0 | 515.6 | 429 |

- W + jets is about half of the tt̄ CR MC, but it contributes mostly untagged QCD-jets. Top-jets are only 18–32% of the pre-tag MC, but 63–81% of the tagged MC.
- The subtracted terms are 19–37% of the tagged MC (B-jet 13–23%, QCD-jet 6–14%). This is why the SF_QCD correction has a visible effect.
- In the γ CR, G + Jets is 94–98% of the MC, so the AN's γ+jets-dominance approximation holds well.

## 3. Step 1 — α = N(Data) / N(all bkg MC)

| Region | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| γ CR (gcr) | 0.478 ± 0.030 | 0.467 ± 0.008 | 0.566 ± 0.005 | 0.545 ± 0.009 |
| tt̄(μ) CR (tmcr) | 1.020 ± 0.058 | 0.888 ± 0.035 | 0.907 ± 0.019 | 0.960 ± 0.030 |
| tt̄(e) CR (tecr) | 0.827 ± 0.061 | 0.724 ± 0.038 | 0.795 ± 0.020 | 0.773 ± 0.031 |
| **tt̄ CR (μ+e)** | **0.941 ± 0.042** | **0.822 ± 0.026** | **0.862 ± 0.014** | **0.880 ± 0.022** |

N(Data) / N(MC), tt̄ CR μ+e: 578/614, 1140/1387, 4406/5111, 1818/2065; γ CR: 1009/2113, 4122/8823, 16121/28472, 3520/6463.

- In the tt̄ CRs, α is close to 1 (MC 6–22% above data) and flat in pT. The electron CR is ~0.1 lower than the muon CR in every bin, a 3–4σ difference in the three upper bins. This points to an electron-side effect (trigger, ID or EGamma data) rather than the tagger.
- In the γ CR, the MC is ~2× the data (α ≈ 0.47–0.57). This is a G + Jets normalization issue, which α corrects pre-tag.

Figures: `plots_effSF_TEST3/alpha/2023preBPix/alpha_gcr.pdf`, `alpha_ttcr_split_zoom.pdf`.

## 4. Step 2 — MC efficiency (eff_MC)

### Top-jet tag efficiency (tt̄ CR, all bkg)

| Region | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| tmcr | 0.708 ± 0.012 | 0.742 ± 0.007 | 0.795 ± 0.003 | 0.811 ± 0.005 |
| tecr | 0.664 ± 0.015 | 0.759 ± 0.008 | 0.792 ± 0.003 | 0.812 ± 0.006 |
| **μ+e** | **0.689 ± 0.010** | **0.749 ± 0.005** | **0.794 ± 0.002** | **0.811 ± 0.004** |

### QCD-jet mis-tag rate (γ CR, G + Jets)

| 250–325 | 325–400 | 400–600 | 600–2000 |
|---:|---:|---:|---:|
| 0.0267 ± 0.0012 | 0.0223 ± 0.0005 | 0.0194 ± 0.0002 | 0.0186 ± 0.0002 |

- The Top-jet efficiency rises with pT from 0.69 to 0.81, as expected for more fully merged tops. μ and e agree within errors.
- The QCD-jet mis-tag rate is ~2%, falling slightly with pT (compare the ~1% background efficiency at the WP from the ROC study).

Figures: `plots_effSF_TEST3/eff/2023preBPix/eff_toptag_combined.pdf`, `eff_qcdjet_gcr_zoom.pdf`.

## 5. Step 3 — Data efficiency and SF

### QCD-jet (γ CR, B-jet subtracted)

| pT [GeV] | α | N(Data) | N(Data, tag) | eff_MC | eff_Data | SF_QCD |
|---|---:|---:|---:|---:|---:|---:|
| 250–325 | 0.478 | 1009 | 70 | 0.0267 | 0.0615 ± 0.0115 | **2.30 ± 0.44** |
| 325–400 | 0.467 | 4122 | 271 | 0.0223 | 0.0587 ± 0.0079 | **2.64 ± 0.36** |
| 400–600 | 0.566 | 16121 | 889 | 0.0194 | 0.0491 ± 0.0059 | **2.54 ± 0.31** |
| 600–2000 | 0.545 | 3520 | 200 | 0.0186 | 0.0517 ± 0.0066 | **2.77 ± 0.36** |

eff_Data uncertainty split: stat 0.009 / 0.004 / 0.002 / 0.004, syst (SF_B, B-jet norm.) 0.007 / 0.007 / 0.006 / 0.005.

The tagged yields agree (data/MC = 0.96, 0.99, 1.16, 1.20), but the pre-tag MC is 2× the data. After α scales the MC down by ~0.5, the data mis-tag rate is ~2.5× the MC rate. SF_QCD ≈ 2.5 is therefore driven by the pre-tag G + Jets normalization, not by an excess of tagged data. It should be cross-checked, for example against the G + Jets k-factor and its pT dependence.

### Top-jet (tt̄ CR, μ+e)

| pT [GeV] | α | N(Data) | N(Data, tag) | eff_MC | eff_Data raw | SF_Top raw | eff_Data SF_QCD-corr. | SF_Top SF_QCD-corr. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 250–325 | 0.941 | 578 | 119 | 0.689 | 0.741 ± 0.193 | **1.08 ± 0.28** | 0.539 ± 0.198 | **0.78 ± 0.29** |
| 325–400 | 0.822 | 1140 | 244 | 0.749 | 0.670 ± 0.120 | **0.89 ± 0.16** | 0.497 ± 0.124 | **0.66 ± 0.17** |
| 400–600 | 0.862 | 4406 | 1178 | 0.794 | 0.661 ± 0.071 | **0.83 ± 0.09** | 0.573 ± 0.073 | **0.72 ± 0.09** |
| 600–2000 | 0.880 | 1818 | 429 | 0.811 | 0.752 ± 0.114 | **0.93 ± 0.14** | 0.598 ± 0.117 | **0.74 ± 0.14** |

eff_Data (raw) uncertainty split: stat 0.147 / 0.079 / 0.032 / 0.062, syst 0.125 / 0.091 / 0.064 / 0.095. The syst part is dominated by SF_B = 1 ± 0.5 on the B-jet tagged term.

Cross-check, lepton-flavour split (SF_Top raw / SF_QCD-corrected):

| pT [GeV] | tmcr raw | tmcr corr. | tecr raw | tecr corr. |
|---|---:|---:|---:|---:|
| 250–325 | 0.98 ± 0.33 | 0.64 ± 0.33 | 1.24 ± 0.40 | 1.01 ± 0.39 |
| 325–400 | 0.92 ± 0.19 | 0.68 ± 0.19 | 0.86 ± 0.20 | 0.64 ± 0.20 |
| 400–600 | 0.86 ± 0.10 | 0.74 ± 0.10 | 0.80 ± 0.10 | 0.69 ± 0.10 |
| 600–2000 | 1.00 ± 0.16 | 0.80 ± 0.16 | 0.82 ± 0.17 | 0.65 ± 0.17 |

- μ and e agree within errors in every bin. The α difference between them (Section 3) cancels in the SF.
- The raw subtraction gives SF_Top ≈ 0.85–1.0 above 325 GeV. The SF_QCD-corrected subtraction multiplies the QCD-jet tagged term (6–14% of the tagged MC) by ~2.5, lowering SF_Top by 0.11–0.30 to ≈ 0.7.
- The data tag rate is below MC in the three upper bins: tagged data / tagged MC = 0.76, 0.75, 0.83 vs α = 0.82, 0.86, 0.88. This deficit is what makes SF_Top < 1.

Figures: `plots_effSF_TEST3/effdata/2023preBPix/effdata_qcdjet_gcr.pdf`, `effdata_topjet_ttcr_combined_zoom.pdf`, `effdata_topjet_tmcr_zoom.pdf`, `effdata_topjet_tecr_zoom.pdf`. The Top-jet figures show the SF_QCD-corrected Data only; the raw-subtraction values are in the tables above.

## 6. Conclusions and next steps

1. TEST3 fixes the MC normalization. α, eff_MC, eff_Data and SF are all physical for 2023preBPix, and the pipeline now follows the AN definitions (all-background α and Top-jet class yields).
2. SF_Top ≈ 0.8–0.9 (raw) or ≈ 0.7 (SF_QCD-corrected), with 10–30% per-bin uncertainty; SF_QCD ≈ 2.3–2.8.
3. Decide between the raw and SF_QCD-corrected subtraction for the nominal SF_Top. The corrected one is more consistent if SF_QCD is trusted, but it inherits the γ CR normalization question.
4. Investigate the 2× G + Jets excess in the γ CR before trusting SF_QCD ≈ 2.5.
5. Check the ~10–20% lower α in the electron tt̄ CR (trigger/ID SFs).
6. Missing: systematics (JES/JER, parton shower, top-pT, scale/PDF, WP), the other eras (2022preEE TEST3 data statistics are too low), and the correctionlib SF file (Stage 3).
