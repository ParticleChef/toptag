# Top-Tagging Scale Factors for the CMS Run 3 Mono-top Analysis

## What this does

Measures data/MC **scale factors (SF) for the ParticleTransformer top tagger**
(`parT_TvsQCD` above a working point) on the leading AK15 jet, for the Run 3
mono-top analysis. It follows the Run 2 method of AN-2020/061, Section 6.6
(`docs/toptag_AN2020_061.pdf`):

- the working point is chosen from the `TvsQCD` ROC curve;
- each AK15 jet is classified from generator truth as **Top-jet**, **B-jet** or **QCD-jet**;
- the **QCD-jet** mistag SF is measured in the γ+jets CR (`gcr`);
- the **Top-jet** tag SF is measured in the tt̄ μ + e CRs (`tmcr` + `tecr`), after
  subtracting the B- and QCD-jet contributions;
- B-jets get `SF_B = 1 ± 0.5`;
- all in four pT bins: 250–325, 325–400, 400–600, 600+ GeV.

The SFs go into a correctionlib JSON, from which a per-event weight is built.
Formulas and definitions: [`docs/method.md`](docs/method.md). Current numbers:
[`docs/results.md`](docs/results.md).

## Workflow

| Stage | What | Runs on | Scripts | Status |
|-------|------|---------|---------|--------|
| 1. Working point | ROC of `TvsQCD` (mono-top vs Z→νν) → WP | hep server | `roc_curve.py`, `plot_TvsQCD.py` | done (2022pre) |
| 2. α, eff, SF | α, eff_MC, eff_Data, SF and SF JSON from the `.scaled` histogram files | laptop | `studies/*.py` | **current** |
| 3. Server pipeline | NanoAOD → histograms → SF JSON → validation → event weight | hep server | `histmaker.py` → `scalefactors.py` → `validate.py` → `weights.py` | implemented, **not used now** |

Environments (details in [`docs/syst_design_computing.md`](docs/syst_design_computing.md)):

```bash
# hep server (Stages 1 and 3)
source ~/envs/p38/bin/activate
# laptop (Stage 2); the system Python 3.7 cannot read .scaled files
PY=/Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python
```

---

## How to run

### Stage 2 — α, efficiencies and SF (current)

Input: a `.scaled` file (LZ4 cloudpickle) in `toptag/`, e.g.
`hadmonotop2023_toptagTEST3.scaled`. It must contain the `fjpt_bin`
(region × pT × TvsQCD) and `fjpt_cat` (+ `n_qinfj` × `n_binfj`, bkg only) histograms.

```bash
F=hadmonotop2023_toptagTEST3        # ".scaled" optional
$PY studies/plot_alpha_toptagTEST.py           -f $F   # 1. α = N(Data)/N(all bkg MC)
$PY studies/plot_toptag_eff_toptagTEST.py      -f $F   # 2. eff_MC
$PY studies/plot_toptag_eff_data_toptagTEST.py -f $F   # 3. eff_Data, SF, SF JSON
$PY studies/plot_class_counts_toptagTEST.py    -f $F   # (optional) jet-class composition
```

Scripts can be run from any directory. Options (shared code in `studies/effsf_cli.py`):

| Option | Meaning |
|--------|---------|
| `-f NAME` | input file; looked up in the current directory, then in `toptag/` |
| `-e ERA` | era; default from the file name: `2022` → 2022preEE, `2022EE` → 2022postEE, `2023` → 2023preBPix, `2023BPix` → 2023postBPix |
| `--alpha-var fjpt` | for old files without `fjpt_bin` |
| `--mc-eff`, `--data-eff` | (step 3) also write single-efficiency plots without the SF panel |
| `--json PATH` | (step 3) SF JSON path |

Step 3 measures SF_QCD first, then SF_Top (raw and SF_QCD-corrected subtraction).

### Stage 1 — Working point from the ROC (hep server)

The tagger working point is a cut on `parT_TvsQCD` of the leading AK15 jet. It is
chosen from a ROC curve that compares top jets from mono-top signal
(`Mphi-1000_Mchi-150`) with non-top jets from Z→νν background (`Zto2Nu`), reading
NanoAOD directly. Events are weighted by sign(genWeight) and must pass a
signal-region-like selection (constants at the top of `roc_curve.py`):

- AK15 jets pass the decaf `isGoodAK15` preselection (pT > 160 GeV, |η| < 2.4,
  tight + tightLepVeto jet ID), and the leading one must have pT > 250 GeV;
- MET > 350 GeV (`MET_PT_MIN`);
- no loose electron, muon or photon;
- no PNet-loose b-tagged AK4 jet outside ΔR < 1.5 of the leading AK15 jet.

Without these cuts the inclusive signal sample has many non-boosted tops and the
ROC falls onto the diagonal (see `log/roc_curve_dev.log`).

- `roc_curve.py` draws signal efficiency against background rejection and against
  background efficiency. For background efficiencies of 5%, 2.5%, 1%, 0.5% and 0.1%
  it prints the `TvsQCD` threshold, the signal efficiency and the rejection. Pick
  the WP from this table.
- `plot_TvsQCD.py` draws the normalized `TvsQCD` shape of each background MC process
  (AN Fig. 21 style), to see where each process sits relative to the cut.

```bash
Y=2022pre
python roc_curve.py --metadata decaf/analysis/metadata/2022_private_v2.json.gz \
    --year $Y --max-files 1 --workers 10 --output plots/                            # quick test
python roc_curve.py --metadata decaf/analysis/metadata/2022_private_v2.json.gz \
    --year $Y --workers 40 --save-arrays hists/roc_arrays_$Y.coffea --output plots/ # full run
python roc_curve.py --load-arrays hists/roc_arrays_$Y.coffea --output plots/        # replot only
python plot_TvsQCD.py --year $Y --workers 40 --save-hists hists/TvsQCD_$Y.coffea --output plots/
python plot_TvsQCD.py --year $Y --load-hists hists/TvsQCD_$Y.coffea --output plots/ # replot only
```

**Fixed 2026-10-25:** before this date `roc_curve.py` printed mirrored thresholds
(t → 0.999 − t) and 1 − AUC. The efficiencies were right. Logs and plots from earlier
runs (`log/roc_2022pre*.log`, `plots_roc/*.png`) still show the old values. The same
table and plots can also be made on the laptop from saved arrays:

```bash
$PY studies/plot_roc_from_arrays.py hists/roc_arrays_$Y.coffea --tag sr   # → plots_roc/recomputed/
```

`plot_TvsQCD.py` supports `--year 2022pre|2022post` only. Job splitting and resource
notes: [`docs/roc_curve_job_notes.txt`](docs/roc_curve_job_notes.txt).

The chosen WP has to be set in three places:
- `histmaker.py`, `TVSQCD_WP` (per year);
- `studies/plot_toptag_eff_data_toptagTEST.py`, `WP` (written into the SF JSON);
- the pass/fail edge of the `TvsQCD` axis of `fjpt_bin` / `fjpt_cat` in the `.scaled`
  files, which is set when the files are produced.

The current value and its efficiencies are in [`docs/method.md`](docs/method.md).

### Stage 3 — Server pipeline (not used now)

Implemented but never run on real NanoAOD; check the per-file docs before trusting it.

```bash
cd /home/jhong/run3Monotop/toptag && mkdir -p hists data plots log
Y=2022pre    # 2022pre | 2022post | 2023pre | 2023post
python classifiers.py                                               # self-check of the truth classification
nohup python histmaker.py --year $Y --workers 40 --output hists/toptag_hists_$Y.coffea \
    > log/histmaker_$Y.log 2>&1 &                                   # 1. fill gcr / tt̄ CR histograms
python scalefactors.py --year $Y --hists hists/toptag_hists_$Y.coffea \
    --output data/toptag_sf_$Y.json.gz                              # 2. SF (correctionlib)
python validate.py --year $Y --hists hists/toptag_hists_$Y.coffea \
    --sf data/toptag_sf_$Y.json.gz --outdir plots/validate_$Y/      # 3. Data/MC before/after
```

Choose `--workers` from the server load (20 / 40 / 60), see
[`docs/syst_design_computing.md`](docs/syst_design_computing.md). Verify
`histmaker.py`'s `METADATA_FILES` mapping against `analysis/run.py` first
([`docs/histmaker.md`](docs/histmaker.md)). Expected ballpark: SF_QCD ~2–2.5, SF_Top ~0.8–1.0.

### Using the SFs in an analysis

```python
import correctionlib
cset = correctionlib.CorrectionSet.from_file("plots_effSF_toptagTEST3/sf_json/toptag_sf_2023preBPix.json")
sf  = cset["toptag_sf"].evaluate(450.0, "top", "nominal")     # pt as float; jetclass top|qcd|b; nominal|up|down
eff = cset["toptag_efficiency"].evaluate(450.0, "top")

from weights import get_toptag_weight                         # per-event weights (server)
nominal, up, down = get_toptag_weight(year="2022pre", jet_pt=fatjet.pt[:, 0],
                                      jet_class=jet_class, is_tagged=fatjet_tagged[:, 0])
```

The weight is `SF` for a tagged jet and `(1 − ε·SF)/(1 − ε)` for an untagged one.

---

## Outputs

### Stage 2: `plots_effSF_<version>/` (e.g. `plots_effSF_toptagTEST3/`)

| Path | Content |
|------|---------|
| `alpha/<era>/` | α vs pT: γ CR, tt̄ CR combined, μ/e split (`alpha_*[_zoom].{pdf,png}`) |
| `eff/<era>/` | eff_MC: Top-jet (tmcr, tecr, combined), QCD-jet mis-tag rate |
| `effdata/<era>/` | eff_MC, eff_Data, SF; Top-jet plots show MC and SF_QCD-corrected Data |
| `sf_json/toptag_sf_<era>.json` | correctionlib SF file (below) |
| `run_log_<era>.txt` | console log, incl. raw-subtraction numbers |

`toptag_sf_<era>.json` (same layout as `scalefactors.py`, readable by `weights.py`):

| Correction | Inputs | Content |
|------------|--------|---------|
| `toptag_sf` | `pt`, `jetclass` (`top`/`qcd`/`b`), `systematic` (`nominal`/`up`/`down`) | SF_Top (tt̄ e+μ, SF_QCD-corrected), SF_QCD (γ CR), SF_B = 1 ± 0.5; up/down = ±1σ total, down clipped at 0 |
| `toptag_efficiency` | `pt`, `jetclass` | eff_MC: Top (tt̄ CR), QCD (γ CR, G + Jets), B (tt̄ CR) |

pT bins are the `fjpt_cat` edges with `flow="clamp"`. Bins with a non-finite SF are
set to 1 ± 1 (a warning is printed).

`plots_classcount_<version>/<era>/{inclusive,toptag_pass,toptag_fail}/`: MC yield per
jet class stacked by process, `classcount_<region>` (class 1–8) and `aggcount_<region>`
(QCD/B/Top), for `wmcr tmcr zmcr wecr tecr gcr zecr` and `ttcr` = `tmcr` + `tecr`.

Plot style: MC `#1f77b4`, Data `#ff7f0e`, SF black.

### Stage 1 and 3

| Path | Content |
|------|---------|
| `<output>/roc_TvsQCD_<year>{_eff,_rej}.{pdf,png}` | ROC vs bkg efficiency / rejection (`plots_roc/` predates the fix: wrong AUC and diagonals) |
| `plots_roc/recomputed/{roc,score}_TvsQCD_<year>_<tag>*` | corrected ROC and TvsQCD distributions (`plot_roc_from_arrays.py`) |
| `<output>/` (from `plot_TvsQCD.py`) | normalized `TvsQCD` shapes of the background MC |
| `hists/roc_arrays_<year>*.coffea` | saved TvsQCD score arrays (`--save-arrays`) |
| `hists/TvsQCD_<year>.coffea` | saved TvsQCD shape histograms (`--save-hists`) |
| `log/` | ROC / TvsQCD run logs |
| `hists/toptag_hists_<year>.coffea` | Stage 3 histograms (not produced yet) |
| `data/toptag_sf_<year>.json.gz` | Stage 3 SF file (not produced yet) |
| `plots/validate_<year>/` | Stage 3 validation plots |

### Reports

`reports/report_toptag_WP_ROC_<era>_<date>.{md,tex,pdf}`: working-point choice from the
ROC (selection, WP table, AUC, caveats).

`reports/report_toptag_effSF_<era>_TEST3_<date>.{md,tex,pdf}`, one per era: steps,
yields, α / eff / SF tables, uncertainty breakdown and caveats. Start here to see
what was measured. Summary table and history: [`docs/results.md`](docs/results.md).

---

## Repository layout

```
toptag/
├── README.md
├── roc_curve.py, plot_TvsQCD.py      # Stage 1
├── studies/                          # Stage 2 (laptop)
│   ├── effsf_cli.py                  #   shared -f / -e handling, era → lumi
│   ├── plot_alpha_toptagTEST.py      #   α
│   ├── plot_toptag_eff_toptagTEST.py #   eff_MC
│   ├── plot_toptag_eff_data_toptagTEST.py  # eff_Data + SF + JSON
│   ├── sf_json.py                    #   correctionlib writer
│   ├── plot_class_counts_toptagTEST.py     # jet-class composition
│   ├── plot_roc_from_arrays.py       #   corrected ROC / WP table from roc_curve.py arrays
│   ├── plot_fjpt_{class,agg}_toptagTEST.py # fjpt per class (no args; read/write two levels up)
│   └── plot_fjpt_toptagTEST.ipynb    #   exploration
├── classifiers.py                    # Stage 3: AK15 → Top/B/QCD (GenPart matching)
├── histmaker.py                      # Stage 3 step 1
├── scalefactors.py                   # Stage 3 step 2
├── validate.py                       # Stage 3 step 3
├── weights.py                        # get_toptag_weight()
├── hadmonotop*_toptagTEST3.scaled    # Stage 2 inputs (2022, 2022EE, 2023, 2023BPix)
├── plots_effSF_toptagTEST3/, plots_classcount_toptagTEST3/, plots_roc/   # outputs
├── hists/, log/                      # Stage 1 arrays and logs
├── reports/                          # current write-ups, one per era
├── samples/                          # local smoke-test ROOT files (2022_private_vToptag.json.gz)
├── analysis/                         # decaf framework, READ-ONLY reference (not modified)
├── old_results/                      # superseded inputs, scripts, plots, reports
└── docs/
```

The fjpt stack scripts and `plot_stack_toptag.py` were merged into `../scripts/plot_stack.py`
(e.g. `plot_stack.py -f <file> --tags fail,pass --unblind --signals none`); originals in
`../scripts/archive/`.

## Documentation

| File | Content |
|------|---------|
| [`docs/method.md`](docs/method.md) | tagger, jet classes, pT bins, α / eff / SF formulas, event weight |
| [`docs/results.md`](docs/results.md) | current SFs per era, open points, history of earlier attempts |
| [`docs/syst_design_computing.md`](docs/syst_design_computing.md) | Run 2 systematics reference, design rules, environments, hep-server resources |
| [`docs/classifiers.md`](docs/classifiers.md) | truth classification details and checks |
| [`docs/histmaker.md`](docs/histmaker.md) | Stage 3 step 1, histogram schema, caveats |
| [`docs/scalefactors.md`](docs/scalefactors.md) | Stage 3 step 2, formulas as implemented, correctionlib schema |
| [`docs/weights_and_validate.md`](docs/weights_and_validate.md) | Stage 3 steps 3–4 |
| [`docs/roc_curve_job_notes.txt`](docs/roc_curve_job_notes.txt) | ROC job splitting |
| `docs/toptag_AN2020_061.pdf` | reference AN (Run 2) |

## Reference

- CMS AN-2020/061, "Search for new physics in final states with large missing
  transverse momentum and a top quark in proton-proton collisions at 13 TeV",
  Section 6.6 (fatjet top-tagging).
- Jet collection: AK15 PFPuppi.
