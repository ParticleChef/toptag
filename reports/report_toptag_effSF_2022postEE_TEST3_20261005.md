# Top-tag α, MC efficiency, Data efficiency, SF and uncertainties — 2022postEE (TEST3)

**Input:** `hadmonotop2022EE_toptagTEST3.scaled` (updated 2026-10-05 20:43)
**Report date:** 2026-10-05
**Tagger / WP:** particleTransformer `TvsQCD` > 0.33, leading AK15 jet
**Luminosity:** 2022postEE 26.68 fb⁻¹ (13.6 TeV)
**PDF version (with figures):** `reports/report_toptag_effSF_2022postEE_TEST3_20261005.pdf`
**Companions:** `report_toptag_effSF_2022preEE_TEST3_20260930.md`, `report_toptag_effSF_2023preBPix_TEST3_20260930.md` (same method)

Definitions follow the AN (AN-2020/061 Sec. 6.6): α and the Top-jet class yields use **all background MC**.

## Summary

- **MC is above data in every region.** α = 0.64–0.75 in the γ CR and 0.78–0.91 in the combined tt̄ CR (2022preEE: 0.90–1.05). Compared per fb⁻¹ with 2022preEE, the MC agrees to within 4% above 400 GeV (≤ 12% higher below, excluding the outlier), while the data is 8–31% lower. The deficit is therefore on the **data side** (lumi, golden JSON / missing files, or a postEE-specific veto applied differently to data and MC) and needs checking.
- **One QCD Multijet MC outlier spoils the lowest tt̄ bin.** In `tmcr` at 250–325 GeV, QCD Multijet gives 157 ± 106 untagged QCD-jet events (≈ 2 events of weight ≈ 75), 13.5% of the tt̄ CR MC in that bin. This makes σ(eff_Data) = 0.70 (raw 0.83) in that bin, and the error is strongly non-linear there. **SF_Top at 250–325 GeV is not usable as it stands.**
- **SF_QCD (γ CR):** 1.62–1.93 (±12–16%). As in 2023 and unlike 2022preEE, it comes mostly from the pre-tag normalization: tagged data / tagged MC = 1.05–1.13 at α ≈ 0.7.
- **SF_Top (tt̄ CR, μ+e), 325–2000 GeV:** 0.81–0.87 with raw subtraction, 0.68–0.76 with SF_QCD-corrected subtraction (±10–22%). The 400–600 GeV bin is the most precise (0.73 ± 0.09 corrected). The raw vs corrected choice moves SF_Top by 0.08–0.13, about 1σ.
- **Uncertainties:** above 325 GeV, **SF_B = 1 ± 0.5** is the largest source (0.065–0.110 on eff_Data), 1.8–2.9× the data statistics. Data statistics are no longer limiting in this era.
- **eff_MC is the same as in 2022preEE** (Top-jet 0.72–0.83, QCD-jet 2.6–3.7%), so the SF differences between the two 2022 eras come from data.

## 1. Steps

All steps run with the `plot_env` interpreter (Python 3.10) and take the `.scaled` file with `-f`. Outputs go to `plots_effSF_<version>/<kind>/<era>/`.

`studies/effsf_cli.py` now accepts an `EE` / `BPix` suffix after the year (`hadmonotop2022EE_…` → 2022postEE, `hadmonotop2023BPix_…` → 2023postBPix). Before this change the `2022EE` name failed era inference.

```bash
cd /Users/hongjieun/working/toptag
PY=/Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python
F=hadmonotop2022EE_toptagTEST3
$PY studies/plot_alpha_toptagTEST.py           -f $F   # Step 1: α
$PY studies/plot_toptag_eff_toptagTEST.py      -f $F   # Step 2: eff_MC
$PY studies/plot_toptag_eff_data_toptagTEST.py -f $F   # Step 3: eff_Data + SF (+ sf_json)
```

Full console output: `plots_effSF_toptagTEST3/run_log_2022postEE.txt`. SF file: `plots_effSF_toptagTEST3/sf_json/toptag_sf_2022postEE.json`. The per-source uncertainty breakdown (Section 6) and the QCD Multijet diagnostic (Section 5.3) were computed with a one-off script that re-uses `measure()` and the yield functions from Step 3; its totals reproduce `measure()` exactly.

| Step | Script | Reads | Computes | Output |
|---|---|---|---|---|
| 1 | `plot_alpha_toptagTEST.py` | `fjpt_bin` (bkg + data) | α per region and pT bin: γ CR, tt̄ CR μ+e, μ/e split | `alpha/2022postEE/` |
| 2 | `plot_toptag_eff_toptagTEST.py` | `fjpt_cat` (bkg) | eff_MC: Top-jet (tmcr, tecr, μ+e, all bkg), QCD-jet (gcr, G + Jets) | `eff/2022postEE/` |
| 3 | `plot_toptag_eff_data_toptagTEST.py` | both | eff_Data and SF: QCD-jet first, then Top-jet (raw and SF_QCD-corrected) | `effdata/2022postEE/` |

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
| TT | 32.6 / 37.0 / 47.3 / 41.5 | 0 |
| ST | 8.7 / 8.7 / 8.5 / 11.2 | 0 |
| W (ℓν) + Jets | 43.7 / 50.6 / 40.8 / 42.6 | 0.3 / 0.2 / 0.3 / 0.6 |
| QCD Multijet | **13.5** / 2.1 / 1.9 / 3.0 | 2.7 / 1.2 / 1.8 / 1.3 |
| G + Jets | 0 | 96.9 / 98.5 / 97.8 / 97.9 |
| other (VV, Z+jets) | 1.5 / 1.6 / 1.6 / 1.6 | < 0.1 |

The 13.5% QCD Multijet in the lowest tt̄ bin is the outlier in `tmcr` (Section 5.3). Without it the composition matches 2022preEE (QCD Multijet 1.7%).

Jet-class yields in the tt̄ CR (μ+e), all background MC:

| pT [GeV] | Top pre-tag | Top tag | B pre-tag | B tag | QCD pre-tag | QCD tag | MC tag | Data tag |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 250–325 | 195.6 | 140.0 | 132.6 | 51.0 | 853.7 | 35.2 | 226.2 | 170 |
| 325–400 | 557.0 | 429.0 | 259.7 | 99.5 | 1541.6 | 74.5 | 602.9 | 497 |
| 400–600 | 2703.6 | 2201.5 | 847.3 | 349.7 | 4397.2 | 201.9 | 2753.1 | 1923 |
| 600–2000 | 699.0 | 582.7 | 429.3 | 154.0 | 1859.7 | 79.2 | 815.9 | 548 |

- Top-jets are 17–34% of the pre-tag MC but 62–80% of the tagged MC.
- The subtracted tagged terms are 20–38% of the tagged MC (B-jet 13–23%, QCD-jet 7–16%); they are largest in the lowest pT bin.
- In the lowest bin the subtracted pre-tag term α·(N(B) + N(QCD)) = 770 is 83% of N(Data) = 923. This small denominator (153 events) is what makes that bin fragile.

## 3. Step 1 — α = N(Data) / N(all bkg MC)

| Region | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| γ CR (gcr) | 0.710 ± 0.020 | 0.753 ± 0.007 | 0.699 ± 0.004 | 0.642 ± 0.008 |
| tt̄(μ) CR (tmcr) | 0.729 ± 0.106 | 0.912 ± 0.028 | 0.847 ± 0.014 | 0.802 ± 0.022 |
| tt̄(e) CR (tecr) | 0.876 ± 0.048 | 0.904 ± 0.032 | 0.788 ± 0.016 | 0.742 ± 0.025 |
| **tt̄ CR (μ+e)** | **0.781 ± 0.075** | **0.909 ± 0.021** | **0.823 ± 0.011** | **0.776 ± 0.017** |

N(Data) / N(MC), tt̄ CR μ+e: 923/1181.9, 2143/2358.3, 6542/7948.1, 2320/2988.0; γ CR: 2746/3865.4, 12936/17174, 34546/49436, 7185/11186.

- α is 0.64–0.91, well below 1 everywhere. It falls with pT above 400 GeV in both CRs.
- The large error at 250–325 GeV in `tmcr` (±15% for 556 data events) is the QCD Multijet outlier. Without it, α(tmcr) = 0.92 in that bin.
- e and μ differ by up to 0.15. The electron CR is higher at low pT and lower above 400 GeV (2.7σ at 400–600 GeV). This pattern is not the same as in 2022preEE or 2023preBPix.

**Where the deficit comes from.** Yield per fb⁻¹, 2022postEE / 2022preEE (lumi ratio 26.68 / 7.99 = 3.34):

| | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| γ CR, data | 0.85 | 0.79 | 0.76 | 0.69 |
| γ CR, MC | 1.12 | 1.09 | 1.03 | 1.03 |
| tt̄ CR, data | 0.92 | 0.86 | 0.82 | 0.83 |
| tt̄ CR, MC | 1.14* | 1.00 | 0.97 | 0.96 |

\* QCD Multijet outlier.

The MC scales with luminosity as expected, but the data does not. Things to check: the 2022postEE luminosity (26.68 fb⁻¹) against the lumi mask used, missing data files in the postEE EGamma / MET datasets, and whether the postEE ECAL-endcap veto region is applied identically to data and MC. The deficit grows with pT in the γ CR, which points to an effect that depends on the jet or photon (for example the veto) rather than a flat luminosity offset.

Figures: `plots_effSF_toptagTEST3/alpha/2022postEE/alpha_gcr_zoom.pdf`, `alpha_ttcr_split_zoom.pdf`.

## 4. Step 2 — MC efficiency (eff_MC)

### Top-jet tag efficiency (tt̄ CR, all bkg)

| Region | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| tmcr | 0.718 ± 0.006 | 0.772 ± 0.003 | 0.815 ± 0.001 | 0.834 ± 0.003 |
| tecr | 0.713 ± 0.007 | 0.768 ± 0.004 | 0.814 ± 0.002 | 0.833 ± 0.003 |
| **μ+e** | **0.716 ± 0.005** | **0.770 ± 0.003** | **0.814 ± 0.001** | **0.834 ± 0.002** |

### QCD-jet mis-tag rate (γ CR, G + Jets)

| 250–325 | 325–400 | 400–600 | 600–2000 |
|---:|---:|---:|---:|
| 0.0366 ± 0.0009 | 0.0297 ± 0.0004 | 0.0276 ± 0.0002 | 0.0257 ± 0.0002 |

- The Top-jet efficiency rises from 0.72 to 0.83 with pT, within 0.003 of 2022preEE in every bin. μ and e agree.
- The QCD-jet mis-tag rate is 2.6–3.7%, the same as 2022preEE (2.6–3.9%).
- With 3.3× the luminosity, the eff_MC errors are about half those of 2022preEE.

Figures: `plots_effSF_toptagTEST3/eff/2022postEE/eff_toptag_combined.pdf`, `eff_qcdjet_gcr_zoom.pdf`.

## 5. Step 3 — Data efficiency and SF

### 5.1 QCD-jet (γ CR, B-jet subtracted)

| pT [GeV] | α | N(Data) | N(Data, tag) | eff_MC | eff_Data | SF_QCD |
|---|---:|---:|---:|---:|---:|---:|
| 250–325 | 0.710 | 2746 | 187 | 0.0366 | 0.0593 ± 0.0096 | **1.62 ± 0.26** |
| 325–400 | 0.753 | 12936 | 732 | 0.0297 | 0.0489 ± 0.0072 | **1.65 ± 0.24** |
| 400–600 | 0.699 | 34546 | 1966 | 0.0276 | 0.0505 ± 0.0061 | **1.83 ± 0.22** |
| 600–2000 | 0.642 | 7185 | 397 | 0.0257 | 0.0495 ± 0.0062 | **1.93 ± 0.24** |

Tagged data / tagged MC in the γ CR = 1.05, 1.11, 1.13, 1.08 at α = 0.64–0.75. Most of SF_QCD therefore comes from dividing by α, as in 2023preBPix (α ≈ 0.5, SF_QCD ≈ 2.5) and unlike 2022preEE (α ≈ 1, tagged excess 1.17–1.35). If the pre-tag data deficit (Section 3) is an artifact, SF_QCD will change substantially. The rise of SF_QCD with pT follows the fall of α.

### 5.2 Top-jet (tt̄ CR, μ+e)

| pT [GeV] | α | N(Data) | N(Data, tag) | eff_MC | eff_Data raw | SF_Top raw | eff_Data SF_QCD-corr. | SF_Top SF_QCD-corr. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 250–325 † | 0.781 | 923 | 170 | 0.716 | 0.672 ± 0.831 | 0.94 ± 1.16 | 0.561 ± 0.699 | 0.78 ± 0.98 |
| 325–400 | 0.909 | 2143 | 497 | 0.770 | 0.670 ± 0.106 | **0.87 ± 0.14** | 0.583 ± 0.110 | **0.76 ± 0.14** |
| 400–600 | 0.823 | 6542 | 1923 | 0.814 | 0.660 ± 0.069 | **0.81 ± 0.08** | 0.598 ± 0.071 | **0.73 ± 0.09** |
| 600–2000 | 0.776 | 2320 | 548 | 0.834 | 0.676 ± 0.121 | **0.81 ± 0.15** | 0.571 ± 0.124 | **0.68 ± 0.15** |

† Not usable: dominated by the QCD Multijet MC outlier (Section 5.3).

Lepton-flavour split (SF_Top raw / SF_QCD-corrected):

| pT [GeV] | tmcr raw | tmcr corr. | tecr raw | tecr corr. |
|---|---:|---:|---:|---:|
| 250–325 | 1.07 ± 17.98 † | 0.90 ± 15.10 † | 0.75 ± 0.26 | 0.62 ± 0.27 |
| 325–400 | 0.91 ± 0.15 | 0.79 ± 0.16 | 0.81 ± 0.16 | 0.71 ± 0.16 |
| 400–600 | 0.82 ± 0.09 | 0.74 ± 0.09 | 0.80 ± 0.09 | 0.73 ± 0.09 |
| 600–2000 | 0.79 ± 0.16 | 0.65 ± 0.16 | 0.85 ± 0.16 | 0.73 ± 0.16 |

- μ and e agree within errors in every usable bin.
- Because SF_QCD ≈ 1.6–1.9, the correction lowers SF_Top by 0.08–0.13. This is larger than in 2022preEE (0.04–0.19 on larger errors), smaller than in 2023preBPix, and about 1σ per bin here.
- Tagged data / tagged MC = 0.75, 0.82, 0.70, 0.67 vs α = 0.78, 0.91, 0.82, 0.78. SF_Top < 1 comes from tagged data falling even further below MC than the pre-tag data does, as in both other eras.
- SF_Top (325–2000 GeV) agrees with 2022preEE (raw 0.74–0.92) and 2023preBPix (corrected 0.66–0.78) within errors.

### 5.3 The QCD Multijet outlier at 250–325 GeV

QCD Multijet MC per jet class in `tmcr`, 250–325 GeV: QCD-jet pre-tag **156.8 ± 105.6**, QCD-jet tagged 2.5 ± 2.5, B-jet 2.9 ± 2.6, Top-jet 0. A variance this large for that yield implies about two events of weight ≈ 75. The same process gives 0–2.5 events per bin in `tecr` and 42–130 (±5–17) at higher pT in `tmcr`.

The events are untagged and land in the QCD-jet pre-tag term, which sits in both α and the eff_Data denominator. The denominator N(Data) − α·(N(B) + N(QCD)) is only 153 events in that bin, so a ±106-event MC uncertainty on N(QCD) dominates. The resulting eff_Data error of 0.67 from the QCD-jet pre-tag term alone (Section 6.3) is also highly asymmetric: a −1σ shift gives a total of 0.52 instead of 0.70.

Diagnostic: SF_Top with QCD Multijet removed from the tt̄ CR MC (α and all subtraction terms recomputed). Removing a whole process is biased, so this is a cross-check only:

| pT [GeV] | α (μ+e) | SF_Top raw | SF_Top corr. | tmcr only, raw | tmcr only, corr. |
|---|---:|---:|---:|---:|---:|
| 250–325 | 0.903 | 0.75 ± 0.22 | 0.61 ± 0.23 | 0.75 ± 0.25 | 0.60 ± 0.25 |
| 325–400 | 0.928 | 0.85 ± 0.14 | 0.74 ± 0.14 | 0.88 ± 0.15 | 0.76 ± 0.15 |
| 400–600 | 0.839 | 0.80 ± 0.08 | 0.72 ± 0.08 | 0.79 ± 0.08 | 0.71 ± 0.09 |
| 600–2000 | 0.800 | 0.79 ± 0.14 | 0.67 ± 0.14 | 0.75 ± 0.15 | 0.63 ± 0.15 |

- Above 325 GeV, SF_Top moves by ≤ 0.02, so the usable bins do not depend on the QCD Multijet MC.
- At 250–325 GeV the bin becomes measurable (±0.22). It then agrees with the e-only result (0.75 ± 0.26 raw) and with the higher bins.
- Possible fixes: drop or cap the outlier events at the processor level (preferred), use the e CR alone for this bin, or merge 250–325 into 325–400.

Figures: `plots_effSF_toptagTEST3/effdata/2022postEE/effdata_qcdjet_gcr.pdf`, `effdata_topjet_ttcr_combined_zoom.pdf`, `effdata_topjet_tmcr_zoom.pdf`, `effdata_topjet_tecr_zoom.pdf`.

## 6. Uncertainties

The horizontal bars in all plots are the bin half-widths, not uncertainties. Histogram variances are Σw² for MC and N (Poisson) for data.

### 6.1 How each error is computed

Same method as the 2022preEE / 2023preBPix TEST3 reports:

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
| N(Data, tag) | 0.0052 | 0.0022 | 0.0013 | 0.0029 |
| N(Data) | 0.0013 | 0.0005 | 0.0003 | 0.0007 |
| N(MC) (α) | 0.0002 | < 0.0001 | < 0.0001 | < 0.0001 |
| B-jet tagged: 50% norm. + MC stat | 0.0056 | 0.0048 | 0.0041 | 0.0038 |
|   of which MC stat | 0.0005 | 0.0002 | 0.0001 | 0.0001 |
| B-jet pre-tag | 0.0010 | 0.0008 | 0.0008 | 0.0009 |
| SF_B = 1 ± 0.5 | 0.0056 | 0.0048 | 0.0041 | 0.0038 |
| **Total** | **0.0096** | **0.0072** | **0.0061** | **0.0062** |
| stat / syst | 0.0053 / 0.0079 | 0.0022 / 0.0068 | 0.0014 / 0.0059 | 0.0029 / 0.0055 |

Relative σ(SF_QCD): from eff_Data 16 / 15 / 12 / 13%, from eff_MC 2.6 / 1.3 / 0.7 / 0.6%.

- The systematics dominate in every bin: the B-jet 50% normalization and SF_B ± 0.5 contribute equally. This includes the lowest bin, which was statistics-limited in 2022preEE.
- The total SF_QCD error does not cover the open normalization question of Section 3.

### 6.3 Breakdown — Top-jet eff_Data (tt̄ CR, μ+e, SF_QCD-corrected)

| Source | 250–325 † | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| N(Data, tag) | 0.085 | 0.044 | 0.020 | 0.043 |
| N(Data) | 0.036 | 0.021 | 0.011 | 0.021 |
| N(MC) (α) | 0.133 | 0.012 | 0.003 | 0.007 |
| B-jet tagged | 0.009 | 0.004 | 0.002 | 0.003 |
| B-jet pre-tag | 0.012 | 0.006 | 0.002 | 0.005 |
| QCD-jet tagged | 0.028 | 0.012 | 0.004 | 0.007 |
| QCD-jet pre-tag | **0.665** | 0.020 | 0.006 | 0.010 |
| SF_B = 1 ± 0.5 | 0.130 | 0.089 | 0.065 | 0.110 |
| SF_QCD (γ CR) | 0.048 | 0.032 | 0.017 | 0.027 |
| **Total** | **0.699** | **0.110** | **0.071** | **0.124** |

Grouped:

| Group | 250–325 † | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| Data stat | 0.092 | 0.049 | 0.022 | 0.048 |
| MC stat | 0.679 | 0.027 | 0.009 | 0.015 |
| SF_B | 0.130 | 0.089 | 0.065 | 0.110 |
| SF_QCD | 0.048 | 0.032 | 0.017 | 0.027 |

Relative σ(SF_Top): from eff_Data 125 / 19 / 12 / 22%, from eff_MC 0.6 / 0.3 / 0.1 / 0.2%.

- **SF_B ± 0.5 is the largest source in every usable bin**, 1.8–2.9× the data statistics. The data statistics improved by about √3.3 over 2022preEE, but SF_B did not, so a B-jet SF measurement or a tighter prior is now clearly the main improvement.
- The MC-stat terms (N(MC) and QCD-jet pre-tag) at 250–325 GeV are the QCD Multijet outlier. Above 325 GeV, MC statistics are at most 0.027.
- In the raw variant the SF_QCD term is absent; the totals are 0.831 / 0.106 / 0.069 / 0.121.
- eff_MC contributes ≤ 0.6% to σ(SF_Top), so the SF error is entirely the eff_Data error.

### 6.4 Validity checks

**Linearity.** Above 325 GeV, repeating the propagation with −1σ shifts changes the totals by < 0.1%. At 250–325 GeV the Top-jet totals drop from 0.699 to 0.518 (corrected) and from 0.831 to 0.647 (raw), so the linear error is not valid there. This is another reason to treat that bin as unmeasured until the outlier is fixed. The QCD-jet totals are unchanged.

**Correlation of tagged and total data.** N(Data, tag) is a subset of N(Data), but they are propagated as independent. Re-parametrizing with the independent counts N(Data, tag) and N(Data, fail) = N(Data) − N(Data, tag):

| | 250–325 † | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| Top-jet raw: data stat, current | 0.092 | 0.049 | 0.022 | 0.048 |
| Top-jet raw: data stat, (tag, fail) | 0.076 | 0.038 | 0.016 | 0.037 |
| Top-jet raw: total, current | 0.831 | 0.106 | 0.069 | 0.121 |
| Top-jet raw: total, (tag, fail) | 0.830 | 0.102 | 0.067 | 0.118 |
| QCD-jet: total, current | 0.0096 | 0.0072 | 0.0061 | 0.0062 |
| QCD-jet: total, (tag, fail) | 0.0094 | 0.0071 | 0.0061 | 0.0061 |

The data-stat part shrinks by 18–27%. Because SF_B dominates here, the Top-jet total shrinks by only 2–4% (6–10% in 2022preEE), and the QCD-jet total by ≤ 2%. The recommendation to switch to the (tag, fail) parametrization stands.

**eff_MC–eff_Data correlation.** N(MC) in α includes the target class, so eff_MC and eff_Data share MC events. With eff_MC contributing ≤ 2.6% (QCD) and ≤ 0.6% (Top) of σ(SF), neglecting this is safe.

**Not included:** JES/JER, parton shower, top-pT reweighting, scale/PDF, pileup, lepton/photon ID and trigger SFs, the WP choice, and the open data-normalization question of Section 3.

## 7. Conclusions and next steps

1. **Data normalization:** α ≈ 0.7–0.9 comes from a data deficit per fb⁻¹ relative to 2022preEE, not from the MC. Check the postEE lumi / lumi mask, data-file completeness, and the ECAL-endcap veto in data vs MC before using these SFs. SF_QCD (1.6–1.9) depends directly on this.
2. **QCD Multijet outlier:** about two high-weight events in `tmcr` 250–325 GeV make that SF_Top bin unusable. Fix them at the processor level (drop or cap), or use the e CR / merge bins.
3. **SF_Top, 325–2000 GeV:** 0.81–0.87 raw, 0.68–0.76 SF_QCD-corrected, ±10–22%, consistent with 2022preEE and 2023preBPix. The raw vs corrected choice (about 1σ) still needs a decision.
4. **SF_B = 1 ± 0.5** is now the dominant Top-jet uncertainty in every usable bin; constraining it is the most effective improvement.
5. Switch the data-stat propagation to (tagged, failed) counts (Top-jet bars −2–4% here).
6. Still missing: the systematics listed in 6.4, 2023postBPix, and the correctionlib SF file (Stage 3). `toptag_sf_2022postEE.json` was written but includes the unusable 250–325 GeV bin as is.
