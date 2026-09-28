# Top-tag α, MC efficiency, Data efficiency and SF — 2022preEE and 2023preBPix

**Inputs:** `hadmonotop2022_toptagTEST.scaled`, `hadmonotop2023_toptagTEST.scaled` (updated 2026-09-28)
**Report date:** 2026-09-28 (revision 2: α uses the calibrated process only)
**Tagger / WP:** particleTransformer `TvsQCD` > 0.33, leading AK15 jet
**Luminosity:** 2022preEE 7.99 fb⁻¹, 2023preBPix 17.96 fb⁻¹ (13.6 TeV)

**Revision 2 change:** α = N(Data) / N(MC) now uses the **total yield of the calibrated process** as N(MC), not the sum of all backgrounds: G + Jets in the γ CR, TT + ST in the tt̄ CRs. The Top-jet class yields (eff_MC and the subtracted B/QCD terms) now also come from TT + ST only, matching the QCD-jet calibration, which already used G + Jets only. Revision 1 (all-background α) is superseded.

## Summary

- Both eras are measured in **4 AK15 pT bins** (250–325, 325–400, 400–600, 600–2000 GeV) on the data side as well.
- **QCD-jet (mis-tag) SF from the γ CR:** 1.09–1.73 in 2022preEE and 2.16–2.70 in 2023preBPix, with 12–23% uncertainty per bin. It moved by 0–13% from revision 1.
- **Top-jet SF from the tt̄ CRs is still NOT usable.** With the TT+ST α, eff_Data (μ+e) is −0.9 to +0.25 (negative in 3 of 4 bins per era) with ±0.9 to ±2.7 uncertainty. Cause: TT + ST is only 3–16% of the tt̄-CR MC, so α ≈ 100–170 and α·N(TT+ST, B+QCD, tag) exceeds the tagged data in most bins.
- Taking |eff_Data| only flips the sign of the Top-jet SF; it does not fix the over-subtraction (section 5).

## 1. Inputs and method

### Histograms

| Histogram | Axes | Filled for | Used for |
|---|---|---|---|
| `fjpt_bin` | region × pT (250, 325, 400, 600, 2000) × TvsQCD | bkg + data | N(Data), N(Data, tag), N(MC) → α |
| `fjpt_cat` | region × pT (same edges) × TvsQCD × n_qinfj × n_binfj | bkg only | class-split MC yields, eff_MC |

TvsQCD axis edges [0, 0.33, 1]: bin 1 = tagged. For the selected processes, `fjpt_bin` MC equals `fjpt_cat` MC summed over all classes in every bin, region and era (checked on every run). The plain `fjpt` histogram is filled twice per event in both files and is not used.

Data dataset per region: `tmcr` → MET, `tecr` and `gcr` → EGamma. Both files contain all 10 background processes, so one file serves both calibrations in each era.

### Jet classes

Top-jet = classes 5–8 (≥2 top-decay quarks in the jet); B-jet = classes 2, 4; QCD-jet = classes 1, 3.

### Formulas

- **α** = N(Data) / N(MC), pre-tag, per region and pT bin, where N(MC) is the total (all jet classes) of the calibrated process: **G + Jets** in the γ CR, **TT + ST** in the tt̄ CRs.
- **QCD-jet (γ CR, G + Jets):**
  eff_Data = [N(Data, tag) − α·N(GJ, B, tag)] / [N(Data) − α·N(GJ, B)]
- **Top-jet (tt̄ CR, TT + ST):**
  eff_Data = [N(Data, tag) − α·N(TT+ST, B+QCD, tag)] / [N(Data) − α·N(TT+ST, B+QCD)]
- **SF** = eff_Data / eff_MC, with eff_MC = N(MC, class, tag) / N(MC, class).

Because α normalizes the calibrated process's total to data, the denominator equals α times that process's target-class yield: N(Data) − α·N(GJ, B) = α·N(GJ, QCD) (G + Jets has no Top-jets) and N(Data) − α·N(TT+ST, B+QCD) = α·N(TT+ST, Top). So SF = [N(Data, tag) − α·N(other classes, tag)] / [α·N(target, tag)], the README's SF formula.

SF_B = 1 (central). Uncertainties: statistics, SF_B = 1 ± 0.5 and (γ CR) a 50% normalization uncertainty on the G + Jets B-jet MC; inputs treated as uncorrelated. No JES, parton-shower, top-pT or scale/PDF systematics.

## 2. α = N(Data) / N(calibrated process)

| Region | N(MC) | Era | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---|---|---:|---:|---:|---:|
| γ CR (gcr) | G + Jets | 2022preEE | 69.2 ± 2.2 | 99.0 ± 1.4 | 13.48 ± 0.12 | 1.73 ± 0.03 |
| γ CR (gcr) | G + Jets | 2023preBPix | 37.8 ± 1.2 | 44.1 ± 0.7 | 8.15 ± 0.06 | 1.41 ± 0.02 |
| tt̄(μ) CR (tmcr) | TT + ST | 2022preEE | 124 ± 11 | 162 ± 9 | 138 ± 4 | 96 ± 5 |
| tt̄(e) CR (tecr) | TT + ST | 2022preEE | 154 ± 14 | 175 ± 11 | 140 ± 5 | 96 ± 5 |
| tt̄ CR (μ+e) | TT + ST | 2022preEE | 137 ± 8 | 167 ± 7 | 139 ± 3 | 96 ± 4 |
| tt̄(μ) CR (tmcr) | TT + ST | 2023preBPix | 149 ± 9 | 139 ± 6 | 137 ± 3 | 111 ± 4 |
| tt̄(e) CR (tecr) | TT + ST | 2023preBPix | 115 ± 9 | 108 ± 6 | 108 ± 3 | 80 ± 3 |
| tt̄ CR (μ+e) | TT + ST | 2023preBPix | 135 ± 6 | 126 ± 4 | 125 ± 2 | 97 ± 3 |

- In the γ CR, α falls steeply with pT. G + Jets is 47–97% of the γ-CR background MC in 2022 and 18–95% in 2023 (lowest at low pT); the rest, mostly QCD Multijet, is now absorbed into α.
- In the tt̄ CRs, α ≈ 80–175. TT + ST is only 3–16% of the tt̄-CR background MC (W+jets ~60% and QCD Multijet ~20% dominate), so α scales TT + ST up to all of the data.
- α errors are now data-statistics dominated; the G + Jets and TT + ST MC statistics are small.

## 3. MC efficiency (eff_MC)

Unchanged from revision 1: all Top-jet MC in the tt̄ CRs comes from TT + ST, and the QCD-jet eff_MC already used G + Jets only.

### Top-jet tag efficiency (tt̄ CR, TT + ST)

| Era | Region | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---|---:|---:|---:|---:|
| 2022preEE | tmcr | 0.693 ± 0.019 | 0.774 ± 0.010 | 0.826 ± 0.004 | 0.863 ± 0.006 |
| 2022preEE | tecr | 0.671 ± 0.023 | 0.765 ± 0.012 | 0.821 ± 0.005 | 0.858 ± 0.007 |
| 2022preEE | μ+e | 0.683 ± 0.015 | 0.770 ± 0.008 | 0.824 ± 0.003 | 0.861 ± 0.005 |
| 2023preBPix | tmcr | 0.683 ± 0.020 | 0.744 ± 0.011 | 0.803 ± 0.004 | 0.820 ± 0.007 |
| 2023preBPix | tecr | 0.632 ± 0.024 | 0.765 ± 0.012 | 0.803 ± 0.005 | 0.840 ± 0.008 |
| 2023preBPix | μ+e | 0.661 ± 0.016 | 0.753 ± 0.008 | 0.803 ± 0.003 | 0.828 ± 0.005 |

### QCD-jet mis-tag rate (γ CR, G + Jets)

| Era | 250–325 | 325–400 | 400–600 | 600–2000 |
|---|---:|---:|---:|---:|
| 2022preEE | 0.0378 ± 0.0009 | 0.0331 ± 0.0005 | 0.0277 ± 0.0001 | 0.0262 ± 0.0001 |
| 2023preBPix | 0.0285 ± 0.0007 | 0.0238 ± 0.0003 | 0.0198 ± 0.0001 | 0.0190 ± 0.0001 |

The Top-jet efficiency rises with pT (0.63–0.86) in both eras and channels. The QCD-jet mis-tag rate falls with pT and is 25–30% lower in 2023 than in 2022. Errors are MC statistics only.

## 4. QCD-jet: Data efficiency and SF (γ CR)

| Era | pT [GeV] | α | eff_MC | eff_Data | SF | SF (rev. 1) |
|---|---|---:|---:|---:|---:|---:|
| 2022preEE | 250–325 | 69.2 | 0.0378 | 0.0654 ± 0.0118 | 1.73 ± 0.32 | 1.81 |
| 2022preEE | 325–400 | 99.0 | 0.0331 | 0.0361 ± 0.0083 | 1.09 ± 0.25 | 1.25 |
| 2022preEE | 400–600 | 13.5 | 0.0277 | 0.0389 ± 0.0060 | 1.40 ± 0.22 | 1.45 |
| 2022preEE | 600–2000 | 1.73 | 0.0262 | 0.0346 ± 0.0064 | 1.32 ± 0.24 | 1.33 |
| 2023preBPix | 250–325 | 37.8 | 0.0285 | 0.0615 ± 0.0111 | 2.16 ± 0.40 | 2.40 |
| 2023preBPix | 325–400 | 44.1 | 0.0238 | 0.0571 ± 0.0079 | 2.40 ± 0.33 | 2.61 |
| 2023preBPix | 400–600 | 8.15 | 0.0198 | 0.0486 ± 0.0056 | 2.45 ± 0.28 | 2.52 |
| 2023preBPix | 600–2000 | 1.41 | 0.0190 | 0.0515 ± 0.0062 | 2.70 ± 0.33 | 2.72 |

γ CR yields per bin:

| Era | pT [GeV] | N(Data) | N(Data, tag) | N(GJ) | N(GJ, B) | N(GJ, B, tag) | α·N(GJ, B) | α·N(GJ, B, tag) | Numerator | Denominator = α·N(GJ, QCD) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2022preEE | 250–325 | 970 | 72 | 14.02 | 0.39 | 0.150 | 26.9 | 10.4 | 61.6 | 943.1 |
| 2022preEE | 325–400 | 4907 | 225 | 49.59 | 1.50 | 0.539 | 148.7 | 53.3 | 171.7 | 4758.3 |
| 2022preEE | 400–600 | 13672 | 630 | 1013.90 | 26.10 | 8.305 | 352.0 | 112.0 | 518.0 | 13320.0 |
| 2022preEE | 600–2000 | 3141 | 129 | 1814.45 | 47.91 | 13.432 | 82.9 | 23.3 | 105.7 | 3058.1 |
| 2023preBPix | 250–325 | 1003 | 70 | 26.52 | 0.71 | 0.264 | 26.7 | 10.0 | 60.0 | 976.3 |
| 2023preBPix | 325–400 | 4105 | 268 | 93.09 | 2.75 | 0.919 | 121.2 | 40.5 | 227.5 | 3983.8 |
| 2023preBPix | 400–600 | 16068 | 885 | 1971.52 | 50.32 | 15.236 | 410.1 | 124.2 | 760.8 | 15657.9 |
| 2023preBPix | 600–2000 | 3511 | 200 | 2498.06 | 63.83 | 16.954 | 89.7 | 23.8 | 176.2 | 3421.3 |

- The B-jet subtraction is small: α·N(GJ, B) is ≤3% of N(Data), and α·N(GJ, B, tag) is 12–24% of N(Data, tag).
- The SF dropped by 0–13% relative to revision 1, most at low pT, where α changed most. The larger α scales up the subtracted G + Jets B-jets.
- 2023 SFs (2.2–2.7) are near the upper end of the Run 2 range quoted in the AN (~2–2.5); 2022 SFs (1.1–1.7) are lower, mostly because eff_MC is lower in 2023 while eff_Data is similar.
- Uncertainties are 12–23% per bin. Statistics dominate at 250–325 GeV, the SF_B / B-jet normalization systematics dominate at 325–600 GeV, and the two are comparable above 600 GeV.
- Non-G+Jets MC (mostly QCD Multijet) no longer enters α or the subtraction; any data from those processes is attributed to G + Jets. It is 3–53% of the γ-CR background MC per bin in 2022 and 5–82% in 2023.

## 5. Top-jet: Data efficiency and SF (tt̄ CR)

### α-corrected subtraction (README formula), tt̄ CR μ+e

| Era | pT [GeV] | α | eff_MC | eff_Data | SF | SF with \|eff_Data\| |
|---|---|---:|---:|---:|---:|---:|
| 2022preEE | 250–325 | 137 | 0.683 | −0.91 ± 2.51 | −1.33 ± 3.68 | 1.33 ± 3.68 |
| 2022preEE | 325–400 | 167 | 0.770 | −0.37 ± 2.72 | −0.47 ± 3.54 | 0.47 ± 3.54 |
| 2022preEE | 400–600 | 139 | 0.824 | 0.25 ± 0.85 | 0.30 ± 1.03 | 0.30 ± 1.03 |
| 2022preEE | 600–2000 | 96 | 0.861 | −0.39 ± 1.47 | −0.45 ± 1.71 | 0.45 ± 1.71 |
| 2023preBPix | 250–325 | 135 | 0.661 | −0.77 ± 2.54 | −1.17 ± 3.84 | 1.17 ± 3.84 |
| 2023preBPix | 325–400 | 126 | 0.753 | −0.38 ± 2.28 | −0.51 ± 3.03 | 0.51 ± 3.03 |
| 2023preBPix | 400–600 | 125 | 0.803 | 0.11 ± 1.13 | 0.14 ± 1.41 | 0.14 ± 1.41 |
| 2023preBPix | 600–2000 | 97 | 0.828 | −0.07 ± 0.86 | −0.09 ± 1.04 | 0.09 ± 1.04 |

The μ-only and e-only channels give the same picture (eff_Data −1.15 to +0.26). **Do not use these numbers.**

- The uncertainties are dominated by SF_B = 1 ± 0.5 (syst 0.65–2.9 vs. stat 0.06–0.6 on eff_Data, all channels), because the subtracted B-jet term is the largest input.
- **Absolute values:** taking |eff_Data| flips the sign in the negative bins and leaves the magnitude and uncertainty unchanged. Every SF stays consistent with both 0 and 1. A negative eff_Data means the subtracted tagged yield exceeds the tagged data; the absolute value hides this rather than fixing it.
- **SF_QCD-corrected variant** (QCD-jet subtraction term scaled by SF_QCD, shown in `plots_effSF/effdata/`), μ+e: 2022 eff_Data −1.5 to +1.0 with ±0.6–4.2, 2023 eff_Data 0.03–0.49 with ±0.12–0.21. The 2023 values are positive but far below eff_MC (SF 0.03–0.73) and rest on the same over-subtraction, so they are not usable either.

### Why it fails

tt̄ CR μ+e, yields per bin:

| Era | pT [GeV] | N(Data) | N(Data, tag) | N(TT+ST) | N(all bkg) | N(TT+ST, B+QCD) | N(TT+ST, B+QCD, tag) | α·N(TT+ST, B+QCD) | α·N(TT+ST, B+QCD, tag) | Numerator: N(Data, tag) − α·N(B+QCD, tag) | Denominator = α·N(TT+ST, Top) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2022preEE | 250–325 | 302 | 67 | 2.21 | 15.30 | 1.96 | 0.72 | 267.4 | 98.4 | −31.4 | 34.6 |
| 2022preEE | 325–400 | 744 | 176 | 4.45 | 34.95 | 3.70 | 1.33 | 618.2 | 222.0 | −46.0 | 125.8 |
| 2022preEE | 400–600 | 2392 | 726 | 17.26 | 212.76 | 14.16 | 4.48 | 1963.4 | 620.4 | 105.6 | 428.6 |
| 2022preEE | 600–2000 | 839 | 186 | 8.73 | 259.38 | 7.61 | 2.37 | 731.2 | 228.0 | −42.0 | 107.8 |
| 2023preBPix | 250–325 | 578 | 119 | 4.30 | 27.67 | 3.74 | 1.31 | 503.6 | 176.4 | −57.4 | 74.4 |
| 2023preBPix | 325–400 | 1140 | 244 | 9.03 | 64.03 | 7.57 | 2.49 | 956.3 | 314.0 | −70.0 | 183.7 |
| 2023preBPix | 400–600 | 4401 | 1176 | 35.26 | 481.08 | 28.99 | 8.74 | 3618.5 | 1090.4 | 85.6 | 782.5 |
| 2023preBPix | 600–2000 | 1818 | 429 | 18.76 | 592.22 | 16.30 | 4.60 | 1580.4 | 446.1 | −17.1 | 237.6 |

- TT + ST is 3–16% of the tt̄-CR background MC. α therefore treats every data event as TT+ST-like and scales TT + ST by ~100–170.
- TT + ST events tag far more often (35–43% overall; B-jets alone 40–54%) than the data (21–30%), which are mostly W+jets and QCD-like. So α·N(TT+ST, B+QCD, tag) is as large as, or larger than, the tagged data, and the numerator is negative in 6 of 8 bins.
- Process makeup of the tt̄-CR background MC (μ+e, both eras): W(ℓν)+Jets ~60%, QCD Multijet ~20%, Z(ℓℓ)+Jets 6–7%, WW 5–6%, ST 4–5%, TT 1.5–1.6%. The CRs are W+jets dominated, not tt̄ dominated.

## 6. Caveats

- **Do not use the Top-jet eff_Data or SF** from these files, with or without absolute values.
- QCD-jet SF uncertainties are stat + SF_B + B-jet normalization only; no JES, PS, top-pT or scale/PDF systematics. Inputs are treated as uncorrelated.
- With α from the calibrated process only, contamination from other processes in the CR is absorbed into α instead of being subtracted. This matters in the γ CR at low pT (up to 53% non-G+Jets MC in 2022, 82% in 2023) and dominates the tt̄ CR.
- The 2022preEE numbers come from a single file and are not comparable with the two-file 2022 results in earlier reports.
- `fjpt` is double-filled in both files (not used here, but should be fixed in the processor).

## 7. Next steps

1. Check TT (and W+jets) cross sections and sum-of-weights normalization, and the tt̄-CR selection. TT should dominate a tt̄ CR; here TT + ST is 3–16%.
2. Consider subtracting the non-calibrated processes explicitly (MC-normalized, or with their own floated normalization) before computing α, or tighten the tt̄ CR (e.g. b-tag requirements) so that TT + ST dominates.
3. Check the QCD Multijet weights in the γ CR at low pT.
4. Revisit the Top-jet SF only after the tt̄ CR composition is fixed.

## 8. Files

All paths relative to `/Users/hongjieun/working/toptag`. Run with `/Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python`.

**Scripts** (`tmp_effSF_20260928/`; copies of `studies/` scripts, modified for the single 2022 file and for α from the calibrated process only):

- `plot_alpha_toptagTEST.py <era>` → `plots_effSF/alpha/<era>/`
- `plot_toptag_eff_toptagTEST.py <era>` → `plots_effSF/eff/<era>/`
- `plot_toptag_eff_data_toptagTEST.py <era>` → `plots_effSF/effdata/<era>/` (Top-jet: raw and SF_QCD-corrected subtraction)
- `plot_toptag_effSF_alpha.py <era>` → `plots_effSF/effSF_alpha/<era>/` (MC eff, α-corrected Data eff and SF only)
- `table_ttcr_yields.py`, `table_gcr_yields.py` → `tables_{ttcr,gcr}_yields_{<era>,all}.{md,csv}`

Logs: `tmp_effSF_20260928/run_log_v2.txt`. Plots and scripts from before this work: `tmp_backup_20260928/`.

**Google Sheets (revision 2):**

- [toptag ttbar CR yields per pT bin (2022preEE, 2023preBPix) — α from TT+ST](https://docs.google.com/spreadsheets/d/1NOY415x0aumCRA8r3uatvWbh_0tcD0kX7fVszxxYXUo/edit)
- [toptag gamma CR yields per pT bin (2022preEE, 2023preBPix) — α from G+Jets](https://docs.google.com/spreadsheets/d/1-oLNKkLn5TiBqc6mzztrrIxdhXNS-3M8yGUCO99fEBE/edit)
