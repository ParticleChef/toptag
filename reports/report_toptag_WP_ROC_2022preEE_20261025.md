# Top-tag working point from the ParT TvsQCD ROC — 2022preEE

**Input:** `hists/roc_arrays_2022pre{,_sr,_test}.coffea` (score arrays from `roc_curve.py`, produced 2026-06-30 / 07-01 on the hep server)
**Report date:** 2026-10-25
**Tagger:** particleTransformer `TvsQCD` = probT / (probT + probQCD), leading AK15 PFPuppi jet
**Samples:** signal mono-top `Mphi-1000_Mchi-150`, background `Zto2Nu` (2022preEE private NanoAOD, `2022_private_v2.json.gz`)
**PDF version (with figures):** `reports/report_toptag_WP_ROC_2022preEE_20261025.pdf`
**Companions:** `report_toptag_effSF_<era>_TEST3_<date>.md` (SF measurements at `TvsQCD > 0.33`)

## Summary

- **The working-point table printed by `roc_curve.py` has the wrong thresholds.** `compute_roc()` pairs each efficiency with a reversed threshold array, so every printed threshold t corresponds to the true cut 1 − 0.001 − t, and the printed AUC is 1 − AUC. The efficiency pairs (ε_sig, ε_bkg) and the ROC curve itself are correct. Recomputing from the saved arrays reproduces the logged numbers exactly with the bug and gives the corrected values below.
- **TvsQCD > 0.33 was taken from the 1% background-efficiency row of this table.** The true cut at ε_bkg = 1% (signal-region selection) is **TvsQCD > 0.666**, with ε_sig = 0.669 (rejection 100×).
- **The cut used everywhere now, TvsQCD > 0.33,** gives ε_sig = 0.813 and **ε_bkg = 3.26%** (rejection 31×) in the same selection, not the intended 1%.
- **AUC = 0.963** (signal-region selection) and 0.910 (no signal-region cuts), not 0.049 / 0.096. The AUC < 0.5 and the two diagonal lines in the old ROC plots, explained in the development log as physics, are artifacts of the same bug.
- **The TEST3 SF measurements are valid for the cut they used (0.33)**, but that cut is a ~3% background-efficiency WP, not the 1% one. Whether to keep 0.33 or move to ~0.67 is a decision for the analysis; moving needs new `.scaled` files (the cut is a bin edge of the `TvsQCD` axis).

## 1. Steps

ROC inputs are produced on the hep server (`~/envs/p38`; NanoAOD access). The corrected table and figures are recomputed on the laptop from the saved arrays.

```bash
# hep server: fill and save the score arrays (signal-region selection)
python roc_curve.py --metadata decaf/analysis/metadata/2022_private_v2.json.gz \
    --year 2022pre --workers 60 --save-arrays hists/roc_arrays_2022pre_sr.coffea --output plots/

# laptop: corrected WP table, ROC and score plots from the saved arrays
PY=/Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python
$PY studies/plot_roc_from_arrays.py hists/roc_arrays_2022pre_sr.coffea --tag sr
$PY studies/plot_roc_from_arrays.py hists/roc_arrays_2022pre.coffea    --tag incl
```

| Run | Arrays | Selection | Console log | Date |
|---|---|---|---|---|
| test | `roc_arrays_2022pre_test.coffea` | AK15 only, `--max-files 1` | `log/roc_curve_dev.log` | 2026-06-30 |
| inclusive | `roc_arrays_2022pre.coffea` | AK15 only | `log/roc_2022pre.log` | 2026-07-01 00:19 |
| signal region | `roc_arrays_2022pre_sr.coffea` | AK15 + SR cuts (Section 1.2) | `log/roc_2022pre_sr.log` | 2026-07-01 01:57 |

Old plots: `plots_roc/roc_TvsQCD_2022pre_{eff,rej}.{pdf,png}` (from the signal-region run; both full runs wrote to the same names, and the legend values ε_sig = 0.669 at 1% match the SR log). Corrected plots: `plots_roc/recomputed/`.

### 1.1 Inputs

| | Datasets | Files | Events selected (SR) | Σ sign(genWeight) (SR) | Events selected (incl.) |
|---|---|---|---|---|---|
| Signal `Mphi-1000_Mchi-150` | 4 chunks | 40 | 5,422 | 3,236 | 16,364 |
| Background `Zto2Nu` | 100 chunks | 962 | 2,070,490 | 1,171,470 | 5,933,357 |

3,079 other datasets in the metadata are skipped. All 1,002 files passed the open/`Events` check.

### 1.2 Selection

Event-level, applied in `ROCProcessor.process()`:

- ≥ 1 good AK15 jet (decaf `isGoodAK15`: pT > 160 GeV, |η| < 2.4, tight + tightLepVeto jet ID); the leading one has pT > 250 GeV.
- *Signal-region run only:*
  - MET > 350 GeV;
  - no loose electron (pT > 10 GeV, cutBased ≥ 1), muon (pT > 10 GeV, loose ID, pfIsoId ≥ 2) or photon (pT > 20 GeV, cutBased ≥ 1, electron veto);
  - no PNet-loose b-tagged AK4 jet (pT > 30 GeV, |η| < 2.4; 2022pre WP 0.047) at ΔR > 1.5 from the leading AK15 jet.

### 1.3 Discriminant and weights

- TvsQCD = (probTopbWqq + probTopbWcs) / (probTopbWqq + probTopbWcs + Σ probQCD{bb,cc,b,c,others}); jets with probT + probQCD = 0 get −1 and are dropped.
- Each event is weighted by **sign(genWeight)** (±1), not by cross section × genWeight. Efficiencies use the net signed sum in the denominator.
- Scores are histogrammed in 1000 bins on [0, 1]; ε(t) is the fraction of the net weight with TvsQCD ≥ t.

## 2. Development history (`log/roc_curve_dev.log`)

Fixed during development (2026-06-30):

1. `hist.Hist` in a `dict_accumulator` failed to merge → scores collected with `column_accumulator`, histograms built afterwards.
2. `mplhep.cms.label(energy=…)` → `com=13.6`.
3. Efficiency denominator used Σ|w| while the numerator used signed weights (ε at t = 0 was 0.58) → net signed sum.
4. Raw genWeights differ by orders of magnitude between the Zto2Nu slices → sign(genWeight).

The log then explains the AUC = 0.09 of the test run as physics (non-boosted tops in the inclusive signal sample) and adds the signal-region cuts to "restore AUC > 0.5". The signal-region run gave AUC = 0.049, lower still. Section 4 shows that both small AUCs and the diagonal lines come from a fifth bug that was not found.

## 3. Results

### 3.1 Signal-region selection (used for the WP)

| Target ε_bkg | Logged threshold | **True threshold** | ε_sig | ε_bkg | Rejection |
|---|---|---|---|---|---|
| 5.0% | 0.784 | **0.215** | 0.855 | 5.00% | 20× |
| 2.5% | 0.591 | **0.408** | 0.786 | 2.50% | 40× |
| 1.0% | 0.333 | **0.666** | 0.669 | 1.00% | 100× |
| 0.5% | 0.202 | **0.797** | 0.575 | 0.50% | 200× |
| 0.1% | 0.057 | **0.942** | 0.346 | 0.10% | 1000× |

AUC: logged 0.0485, **true 0.963**.

Fixed cuts (true):

| Cut | ε_sig | ε_bkg | Rejection |
|---|---|---|---|
| TvsQCD > 0.26 (Run 2 AN WP) | 0.843 | 4.20% | 24× |
| **TvsQCD > 0.33 (used now)** | **0.813** | **3.26%** | **31×** |

### 3.2 Without signal-region cuts

| Target ε_bkg | Logged threshold | **True threshold** | ε_sig | Rejection |
|---|---|---|---|---|
| 5.0% | 0.715 | **0.284** | 0.721 | 20× |
| 2.5% | 0.508 | **0.491** | 0.631 | 40× |
| 1.0% | 0.270 | **0.729** | 0.493 | 100× |
| 0.5% | 0.160 | **0.839** | 0.396 | 200× |
| 0.1% | 0.047 | **0.952** | 0.222 | 1004× |

AUC: logged 0.0957, **true 0.910**. TvsQCD > 0.33 gives ε_sig = 0.704, ε_bkg = 4.26%.

- The signal-region cuts raise ε_sig at fixed ε_bkg by 0.12–0.18: they remove signal events without a boosted hadronic top, which is the physics the development log describes. They do not change the AUC from < 0.5 to > 0.5; both AUCs are well above 0.5.
- The test run (1 file per chunk) agrees with the inclusive run within statistics (true AUC 0.913; 1% WP at 0.729). Its logged "WP = 0.26: ε_sig = 0.472, ε_bkg = 0.94%" is also mirrored: that point is the cut at ~0.74.

Figures: `plots_roc/recomputed/roc_TvsQCD_2022pre_sr_eff.pdf`, `roc_TvsQCD_2022pre_sr_rej.pdf`, `score_TvsQCD_2022pre_sr.pdf` (and the `_incl` versions).

## 4. The threshold bug in `roc_curve.py`

In `compute_roc()`:

```python
sig_pass = np.cumsum(sig_hist[::-1])[::-1]       # sig_pass[i] = weight with score >= edges[i]  (ascending t)
sig_eff  = np.concatenate([[0.], sig_pass / sig_total, [1.]])
thrs     = np.concatenate([[1.], edges[:-1][::-1], [0.]])   # descending t
```

`sig_pass` (and `bkg_pass`) are ordered by **ascending** threshold, but `thrs` is **descending**. Entry k + 1 therefore pairs the efficiency of the cut at k/1000 with the threshold (999 − k)/1000:

- printed threshold = 0.999 − true threshold (0.333 ↔ 0.666; 0.270 ↔ 0.729);
- the efficiency pairs are internally consistent, so the ROC curve and the ε_sig values at each target ε_bkg are right;
- with the padding points (0, 0) at the start and (1, 1) at the end, the curve runs (0,0) → (1,1) → along the ROC to (~0, ~0) → (1,1). The two straight segments are the "two diagonal lines" in the old plots, and `np.trapz` integrates the ROC backwards, giving 1 − AUC.

**Evidence.** Applying the same indexing to the saved arrays reproduces the logged values exactly (SR: threshold 0.3330, ε_sig 0.669, AUC 0.0485; inclusive: 0.2700, 0.493, 0.0957). The logged thresholds also *decrease* as the target ε_bkg decreases, which is impossible for a cut "TvsQCD > t". The score distributions (`score_TvsQCD_2022pre_sr.pdf`) show directly that the background tail above 0.33 is ~3%: the signal peaks at 1 and the background at 0, as expected.

**Fix** (applied to `roc_curve.py` on 2026-10-25): the cumulative sums are kept in descending-threshold order, so each efficiency is paired with its own threshold. Tested on the saved arrays: the fixed `compute_roc()` gives the corrected tables in Section 3 exactly (1% WP at 0.666, AUC 0.963). Not yet rerun on the server.

## 5. Consequences for the analysis

- **What the docs said:** "WP 0.33: Top-jet efficiency ~67%, QCD rejection ~99%" (README / `docs/method.md`). Those are the ε_sig and 1 − ε_bkg at the true cut 0.666. At 0.33, the numbers are ε_sig ≈ 81% and ε_bkg ≈ 3.3% for this signal and Z→νν.
- **TEST3 SFs** (`report_toptag_effSF_*_TEST3_*`) measure data/MC for TvsQCD > 0.33 and are self-consistent for that cut. The measured eff_MC (Top-jet 0.72–0.83, QCD-jet mis-tag 2–4%) fits a ~3% background WP better than a 1% one.
- **To move to the 1% WP** (TvsQCD > ~0.67): change `TVSQCD_WP` in `histmaker.py`, `WP` in `studies/plot_toptag_eff_data_toptagTEST.py`, and the TvsQCD axis edges of `fjpt_bin` / `fjpt_cat` in the histogram production, then regenerate the `.scaled` files and rerun the SF chain.

## 6. Caveats

- **Weights.** sign(genWeight) treats every event equally, so the Zto2Nu pT/HT slices enter in proportion to their generated statistics, not their cross sections. ε_bkg is therefore not the physical Z→νν mistag rate. A cross-section-weighted ROC (xs × genWeight / Σ genWeight per dataset) is needed before quoting background efficiencies.
- **One signal point** (m_φ = 1000, m_χ = 150 GeV). Other masses have different top pT spectra.
- **2022preEE only.** No ROC for 2022postEE or 2023.
- **Background = Z→νν only.** The SF measurement uses γ+jets (QCD-jet) and tt̄ (Top-jet), and the SR background also includes W+jets and tt̄.
- The luminosity label on the ROC plots (7.98 fb⁻¹) is cosmetic; the samples are not normalized to it.
- `plot_TvsQCD.py` (background shapes) never completed: both logs end in `ModuleNotFoundError: No module named 'hist'` on the server.

## 7. Next steps

1. ~~Fix `compute_roc()` in `roc_curve.py`~~ (done 2026-10-25). Rerun on the server to refresh `plots_roc/` and the logs.
2. Decide the WP: keep 0.33 (~3% ε_bkg, higher signal efficiency) or move to ~0.67 (1% ε_bkg). Update `docs/method.md` with the decision and the correct efficiencies either way.
3. Redo the ROC with cross-section weights, more signal mass points, and the other eras.
4. Install `hist` in `p38` and run `plot_TvsQCD.py`.
