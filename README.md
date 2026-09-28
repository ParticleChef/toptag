# Top-Tagging Weight Calculation for CMS Run 3 Mono-top Analysis

> **Note (2026-09-28):** `scalefactors.py`, `validate.py` and `weights.py` are not currently used; the α / efficiency / SF results come from the `studies/`-based scripts in `tmp_effSF_20260928/` (see [Current α / Efficiency / SF Workflow](#current-α--efficiency--sf-workflow-local-scaled-files) and `reports/report_toptag_effSF_alpha_20260928.md`).

## Overview

This project implements the **top-tagging scale factor (SF) and event weight** calculation
for the CMS mono-top analysis, adapting the Run 2 methodology documented in
`docs/toptag_AN2020_061.pdf` (Section 6.6) to the Run 3 framework used in `decaf/`.

The goal is to derive per-event reweighting factors that correct MC simulation to match
data for the `particleTransformer` top-tagger applied to large-radius (AK15) jets.

Top-tagging working point is measured: 0.33

### Workflow at a glance

| Stage | What | Where it runs | Scripts | Status |
|-------|------|---------------|---------|--------|
| 1. Working point | ROC of ParT `TvsQCD` (mono-top signal vs Z→νν) and TvsQCD shapes → WP 0.33 | hep server (`~/envs/p38`) | `roc_curve.py`, `plot_TvsQCD.py` | done for 2022pre — see [Working Point (ROC)](#working-point-roc) |
| 2. α, eff_MC, eff_Data, SF | α = N(Data)/N(MC), MC tag efficiency, data efficiency and SF per pT bin, from the `.scaled` files | laptop (`plot_env`) | `tmp_effSF_20260928/*.py` | **current** — see [Current α / Efficiency / SF Workflow](#current-α--efficiency--sf-workflow-local-scaled-files) |
| 3. SF file → event weight | correctionlib SF file, validation plots, per-event weight | hep server | `histmaker.py` → `scalefactors.py` → `validate.py` → `weights.py` | implemented, **not used now** — see [Running Commands](#running-commands-step-by-step) |

---

## Physics Background

### Tagger

- **Algorithm**: `parT` (particleTransformer)
- **Applied to**: leading AK15 fatjet in the event
- **Discriminant branch**: `parT_TvsQCD` (top vs QCD score)
- **Working point**: 0.33 measured for run3
  - Signal efficiency (Top-jet): ~67%
  - Background rejection (QCD-jet): ~99%

### AK15 Jet Classification

Each AK15 jet is categorized based on how many quarks from a hadronic top decay
are clustered into it (using ΔR ≤ 1.5 matching) and whether b-flavor is found
(via ghost-hadron matching):

| Class | Description |
|-------|-------------|
| 1 | No top-decay quark, no b-flavor |
| 2 | No top-decay quark, b-flavor found |
| 3 | One top-decay quark, no b-flavor |
| 4 | One top-decay quark, b-flavor found |
| 5 | Two top-decay quarks, no b-flavor |
| 6 | Two top-decay quarks, b-flavor found |
| 7 | Three top-decay quarks |
| 8 | More than three top-decay quarks |

These are merged into **three aggregated classes** for calibration:

| Aggregated Class | Original Classes | Description |
|-----------------|-----------------|-------------|
| **Top-jet** | 5, 6, 7, 8 | Fully/mostly merged hadronic top |
| **B-jet**   | 2, 4         | b-flavor present, no full top merge |
| **QCD-jet** | 1, 3         | Pure QCD, no b-flavor |

### pT Bins

Calibration is performed in four bins of the leading AK15 jet pT:

| Bin | pT range |
|-----|----------|
| 1 | 250 ≤ pT < 325 GeV |
| 2 | 325 ≤ pT < 400 GeV |
| 3 | 400 ≤ pT < 600 GeV |
| 4 | 600 GeV ≤ pT |

---

## Scale Factor Derivation

### Common: Data/MC Normalization Factor α

In each pT bin, before applying any top-tagging requirement:

```
α = N(Data) / N(MC)
```

This ratio is taken from the relevant control region and corrects for pre-existing
data/MC normalization differences independent of the tagger.

**Current definition (2026-09-28, report revision 2):** N(MC) is the total pre-tag
yield (all jet classes) of the **process being calibrated** in that region, not the
sum of all backgrounds:

| Calibration | Region | N(MC) in α | Class yields (eff_MC, subtracted terms) |
|-------------|--------|------------|-----------------------------------------|
| QCD-jet | `gcr` | G + Jets | G + Jets |
| Top-jet | `tmcr` + `tecr` | TT + ST | TT + ST |

N(Data) and N(MC) come from the `fjpt_bin` histogram, and the class-split yields
from `fjpt_cat` (see [Current α / Efficiency / SF Workflow](#current-α--efficiency--sf-workflow-local-scaled-files)).
Because α normalizes the calibrated process's total to data,
N(Data) − α·N(other classes) = α·N(target class), so ε_Data / ε_MC reduces to the
SF formulas below.

---

### QCD-jet Calibration (from γ+jets CR)

The γ+jets control region is dominated by QCD-jets and provides the QCD mistag
scale factor.

**MC efficiency**:
```
ε_QCD-jet,MC = N(γ+jets, QCD-jet, tagged) / N(γ+jets, QCD-jet)
```

**Data efficiency** (B-jet contribution subtracted):
```
ε_QCD-jet,Data = [N(Data, tagged) - α·N(γ+jets, B-jet, tagged)]
               / [N(Data) - α·N(γ+jets, B-jet)]
```

**Scale factor**:
```
SF_QCD-jet = ε_QCD-jet,Data / ε_QCD-jet,MC
           = [N(Data, tagged) - α·N(γ+jets, B-jet, tagged)]
           / [α·N(γ+jets, QCD-jet, tagged)]
```

---

### Top-jet Calibration (from tt̄ e + μ CRs combined)

The tt̄ electron and muon control regions are used together to derive the
Top-jet tag efficiency scale factor. In the current workflow, "MC" below means
**TT + ST** (see the α definition above).

**MC efficiency**:
```
ε_Top-jet,MC = N(MC, Top-jet, tagged) / N(MC, Top-jet)
```

**Data efficiency** (B/QCD-jet contribution subtracted using previously derived SFs):
```
ε_Top-jet,Data = [N(Data, tagged) - α·N(MC, B/QCD-jet, tagged)]
               / [N(Data) - α·N(MC, B/QCD-jet)]
```

**Scale factor**:
```
SF_Top-jet = ε_Top-jet,Data / ε_Top-jet,MC
           = [N(Data, tagged) - α·N(MC, B/QCD-jet, tagged)]
           / [α·N(MC, Top-jet, tagged)]
```

---

### B-jet Scale Factor

No dedicated control region is available for the B-jet class. A conservative
per-bin scale factor is assumed:

```
SF_B-jet = 1 ± 0.5
```

---

### Event Weight

The top-tagging weight is applied analogously to the AK4 b-tagging weight
(promotion/demotion method). For the leading AK15 jet:

```
w_top-tag = ∏_{j∈tagged}   ( ε_j,Sim · SF_j / ε_j,Sim )
          × ∏_{j∈untagged} ( (1 - ε_j,Sim · SF_j) / (1 - ε_j,Sim) )
```

Since only the leading AK15 jet enters the product, this simplifies to:

- **Tagged jet**:   `w = SF_j`
- **Untagged jet**: `w = (1 - ε_j,Sim · SF_j) / (1 - ε_j,Sim)`

where `j` is the aggregated class (Top/B/QCD) of the leading AK15 jet and
the SF and efficiency are taken from the appropriate pT bin.

---

## Run 2 Systematics, Design Rules and Computing Environment

Moved to [`docs/syst_design_computing.md`](docs/syst_design_computing.md):
the Run 2 (AN-2020/061) systematic-uncertainty tables for the QCD-jet and Top-jet
SFs, the project design rules, and the hep-server resources and worker-count
guidance for the coffea jobs.

---

## Working Point (ROC)

Stage 1: choose the `TvsQCD` working point. Runs on the hep server with `~/envs/p38`
(needs `uproot`/coffea NanoAOD processing, which the laptop `plot_env` does not have,
even for replotting).

| Script | What it does | Outputs |
|--------|--------------|---------|
| `roc_curve.py` | ROC of ParT `TvsQCD` on the leading AK15 jet (pT > 250 GeV; isGoodAK15-style preselection). Signal: mono-top `Mphi-1000_Mchi-150`; background: `Zto2Nu`. Prints the threshold, signal efficiency and rejection at background efficiencies 5%, 2.5%, 1%, 0.5%, 0.1%. | `roc_TvsQCD_<year>_rej.{pdf,png}` (vs rejection), `roc_TvsQCD_<year>_eff.{pdf,png}` (vs bkg efficiency); optional saved arrays |
| `plot_TvsQCD.py` | Normalized `TvsQCD` shape comparison of background MC (AN Fig. 21 style), `--year 2022pre` or `2022post`. | shape plots; optional saved histograms |

```bash
source ~/envs/p38/bin/activate
# full run, saving the score arrays for fast replotting
python roc_curve.py --metadata decaf/analysis/metadata/2022_private_v2.json.gz \
    --year 2022pre --workers 40 --save-arrays hists/roc_arrays_2022pre.coffea --output plots/
# replot only
python roc_curve.py --load-arrays hists/roc_arrays_2022pre.coffea --output plots/
python plot_TvsQCD.py --year 2022pre --workers 40 --output plots/
```

Existing results: `hists/roc_arrays_2022pre{,_sr,_test}.coffea` and `plots_roc/`
(`roc_TvsQCD_2022pre{,_eff,_rej}.{pdf,png}`). Job-splitting and resource notes:
`docs/roc_curve_job_notes.txt` (its example commands say `--save-hists`; the
`roc_curve.py` option is `--save-arrays`). The chosen WP, `TvsQCD > 0.33`, is the one
used everywhere below.

---

## Current α / Efficiency / SF Workflow (local, `.scaled` files)

Stage 2, the one used for the current results. Runs on the laptop with the
`plot_env` interpreter (the system Python 3.7 cannot read the `.scaled` pickles):

```bash
cd /Users/hongjieun/working/toptag
PY=/Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python
for era in 2022preEE 2023preBPix; do
  $PY tmp_effSF_20260928/plot_alpha_toptagTEST.py        $era   # α
  $PY tmp_effSF_20260928/plot_toptag_eff_toptagTEST.py   $era   # eff_MC
  $PY tmp_effSF_20260928/plot_toptag_eff_data_toptagTEST.py $era  # eff_Data + SF (raw and SF_QCD-corrected)
  $PY tmp_effSF_20260928/plot_toptag_effSF_alpha.py      $era   # eff_MC, α-corrected eff_Data, SF only
done
$PY tmp_effSF_20260928/table_ttcr_yields.py   # tt̄ CR yield tables, both eras
$PY tmp_effSF_20260928/table_gcr_yields.py    # γ CR yield tables, both eras
```

**Inputs** (in `toptag/`): `hadmonotop2022_toptagTEST.scaled` (2022preEE, 7.99 fb⁻¹)
and `hadmonotop2023_toptagTEST.scaled` (2023preBPix, 17.96 fb⁻¹). Each file holds
all 10 background processes, data (MET for `tmcr`, EGamma for `tecr`/`gcr`) and
these histograms:

| Histogram | Axes | Used for |
|-----------|------|----------|
| `fjpt_bin` | region × pT (250, 325, 400, 600, 2000) × TvsQCD | N(Data), N(Data, tag), N(MC) → α |
| `fjpt_cat` | region × pT (same edges) × TvsQCD × `n_qinfj` × `n_binfj` (bkg only) | jet-class yields, eff_MC |

`fjpt` is filled twice per event in both files and is not used.

**Scripts** (`tmp_effSF_20260928/`): copies of the three `studies/` α/eff scripts,
changed to read one file per era and to take N(MC) from the calibrated process
(G + Jets / TT + ST). Each takes the era as its only argument (default
`2022preEE`), reads the `.scaled` file from `toptag/` and writes plots into
`toptag/plots_effSF/`. The `studies/` originals still point at the old
two-file 2022 setup (`hadmonotop2022_toptagTEST2.scaled`, no longer present).

**Outputs**:

| Output | Content |
|--------|---------|
| `plots_effSF/alpha/<era>/` | α vs pT: γ CR, tt̄ CR combined, μ/e split |
| `plots_effSF/eff/<era>/` | eff_MC: Top-jet (tmcr, tecr, combined), QCD-jet mis-tag rate |
| `plots_effSF/effdata/<era>/` | eff_MC, eff_Data and SF; Top-jet shows raw and SF_QCD-corrected subtraction |
| `plots_effSF/effSF_alpha/<era>/` | eff_MC, α-corrected eff_Data and SF only |
| `tmp_effSF_20260928/tables_{ttcr,gcr}_yields_{<era>,all}.{md,csv}` | per-bin N(Data), N(Data, tag), α, class yields and their α-scaled values |
| `reports/report_toptag_effSF_alpha_20260928.md` | write-up of the results (revision 2) |

Google Drive copies (revision 2): [report](https://docs.google.com/document/d/1ch1BqwKtkDHoGtmisBU9RkQPfg9hQ_8hMm0CylD8FOI/edit),
[tt̄ CR sheet](https://docs.google.com/spreadsheets/d/1NOY415x0aumCRA8r3uatvWbh_0tcD0kX7fVszxxYXUo/edit),
[γ CR sheet](https://docs.google.com/spreadsheets/d/1-oLNKkLn5TiBqc6mzztrrIxdhXNS-3M8yGUCO99fEBE/edit).

**Status**: the QCD-jet SF (γ CR) is usable (2022preEE 1.1–1.7, 2023preBPix
2.2–2.7). The Top-jet SF is **not**: TT + ST is only 3–16% of the tt̄ CR MC, so
α ≈ 100–170 and the subtraction removes more tagged events than the data has. See
the report for details. Plots and scripts from before 2026-09-28 are kept in
`tmp_backup_20260928/`.

---

## Running Commands (Step by Step)

> Stage 3 (server pipeline). Implemented but **not used for the current results**;
> `scalefactors.py` still uses the all-background α and the 250–3000 GeV last bin.

### 0. Setup environment

```bash
source ~/envs/p38/bin/activate
cd /home/jhong/run3Monotop/toptag
mkdir -p hists data plots log
```

### Step 1 — Fill histograms (γ+jets and tt̄ CRs)

Runs a coffea processor over NanoAOD files. Automatically filters for only
the γ+jets and tt̄ datasets needed. Uses 70 workers by default.

```bash
# Default: 40 workers. Increase to 60 if server is quiet, decrease to 20 if busy.
python histmaker.py --year 2022pre  --workers 40 --output hists/toptag_hists_2022pre.coffea
python histmaker.py --year 2022post --workers 40 --output hists/toptag_hists_2022post.coffea
python histmaker.py --year 2023pre  --workers 40 --output hists/toptag_hists_2023pre.coffea
python histmaker.py --year 2023post --workers 40 --output hists/toptag_hists_2023post.coffea
```

To run in the background and keep a log (useful for long runs):
```bash
nohup python histmaker.py --year 2022pre --workers 40 --output hists/toptag_hists_2022pre.coffea \
    > log/histmaker_2022pre.log 2>&1 &
echo "PID: $!"
```

Check progress:
```bash
tail -f log/histmaker_2022pre.log
```

Input: `analysis/metadata/<file>.json.gz` per `histmaker.py`'s `METADATA_FILES`
CONFIG dict (files at `/data/mc/privatemc/` per each metadata entry) —
**verify this mapping against whatever `analysis/run.py` normally uses for
each year/era before trusting it**; see `docs/histmaker.md`.
Output: `hists/toptag_hists_<year>.coffea`

### Step 2 — Derive scale factors

Reads the histograms from Step 1. Runs in seconds (pure Python, no ROOT files).

```bash
python scalefactors.py --year 2022pre  --hists hists/toptag_hists_2022pre.coffea  --output data/toptag_sf_2022pre.json.gz
python scalefactors.py --year 2022post --hists hists/toptag_hists_2022post.coffea --output data/toptag_sf_2022post.json.gz
python scalefactors.py --year 2023pre  --hists hists/toptag_hists_2023pre.coffea  --output data/toptag_sf_2023pre.json.gz
python scalefactors.py --year 2023post --hists hists/toptag_hists_2023post.coffea --output data/toptag_sf_2023post.json.gz
```

Input: `hists/toptag_hists_<year>.coffea`  
Output: `data/toptag_sf_<year>.json.gz` (correctionlib format)

Each run prints a table like:
```
Year: 2022pre
pT bin        SF_QCD          SF_Top
250-325 GeV   2.31 ± 0.42     0.91 ± 0.15
325-400 GeV   2.18 ± 0.27     0.95 ± 0.11
400-600 GeV   2.05 ± 0.25     0.98 ± 0.09
600+ GeV      1.98 ± 0.33     1.02 ± 0.16
```

### Step 3 — Validate (optional)

```bash
python validate.py --year 2022pre \
    --hists hists/toptag_hists_2022pre.coffea \
    --sf    data/toptag_sf_2022pre.json.gz \
    --outdir plots/validate_2022pre/
```

Produces Data/MC comparison plots of `parT_TvsQCD` before and after
calibration (reproducing AN Fig. 29/30 style), saved to `plots/validate_<year>/`.

### Step 4 — Use the weight in your analysis

```python
from weights import get_toptag_weight

# nominal, up, down are per-event arrays
nominal, up, down = get_toptag_weight(
    year      = '2022pre',
    jet_pt    = fatjet.pt[:, 0],        # leading AK15 jet pT
    jet_class = jet_class,              # 'top', 'b', or 'qcd'
    is_tagged = fatjet_tagged[:, 0],    # bool: parT_TvsQCD > 0.33 (Run 3 WP, see below)
)
```

---

## Code Structure

```
toptag/
├── README.md               # This file — the complete reference
├── docs/
│   ├── toptag_AN2020_061.pdf           # Reference AN (Run 2, 2018)
│   ├── classifiers.md                  # classifiers.py details + verification steps
│   ├── histmaker.md                    # histmaker.py details + verification steps
│   ├── scalefactors.md                 # scalefactors.py formulas, correctionlib schema, verification
│   ├── syst_design_computing.md   # Run 2 systematics tables, design rules, server computing notes
│   └── weights_and_validate.md         # weights.py / validate.py details
├── analysis/                           # Main analysis framework (READ ONLY reference)
│   ├── libs/mycoffeav2.py              # CustomNanoAODSchema (AK15 branch aliasing)
│   ├── utils/ids.py                    # isGoodAK15 and other object ID helpers
│   ├── utils/corrections.py            # Pattern for SF/weight functions
│   ├── utils/common.py                 # BTag WPs, helpers
│   └── data/                           # Existing SF data files
│
├── roc_curve.py            # Stage 1: ROC of TvsQCD (mono-top vs Z→νν), WP table
├── plot_TvsQCD.py          # Stage 1: normalized TvsQCD shapes of background MC
├── plots_roc/              # Stage 1 outputs: roc_TvsQCD_2022pre{,_eff,_rej}.{pdf,png}
│
├── hadmonotop2022_toptagTEST.scaled  # Stage 2 input, 2022preEE (LZ4 cloudpickle)
├── hadmonotop2023_toptagTEST.scaled  # Stage 2 input, 2023preBPix
├── tmp_effSF_20260928/     # Stage 2 scripts (CURRENT): α, eff_MC, eff_Data/SF plots,
│                           #   yield-table scripts, tables_*.{md,csv}, run logs
├── plots_effSF/            # Stage 2 outputs: alpha/, eff/, effdata/, effSF_alpha/ per era
├── tmp_backup_20260928/    # plots and scripts from before 2026-09-28
│
├── classifiers.py          # [done, unverified] AK15 jet → Top/B/QCD class (GenPart matching)
├── histmaker.py            # [done, unverified] Stage 3 step 1: fill CR histograms → .coffea file
├── scalefactors.py         # [done, verified with synthetic histograms; NOT USED NOW] Step 2: histograms → SF JSON (correctionlib)
├── weights.py              # [done, unverified; NOT USED NOW] Step 4: get_toptag_weight() for use in analysis
├── validate.py             # [done, unverified; NOT USED NOW] Step 3: Data/MC comparison plots
├── samples/                 # local smoke-test ROOT files (one per dataset in
│                             #   2022_private_vToptag.json.gz)
│
├── hists/                  # Stage 1 ROC arrays (roc_arrays_2022pre*.coffea);
│   └── toptag_hists_<year>.coffea   #   Stage 3 step 1 outputs (not produced yet)
├── data/                   # Stage 3 step 2 outputs — SF files (not produced yet)
│   └── toptag_sf_<year>.json.gz
├── plots/                  # Stage 3 step 3 outputs (git-ignored)
│   └── validate_<year>/
│
├── studies/                # LOCAL plotting/validation scripts (see below)
│   ├── plot_alpha_toptagTEST.py
│   ├── plot_toptag_eff_toptagTEST.py
│   ├── plot_toptag_eff_data_toptagTEST.py
│   ├── plot_fjpt_inclusive_toptagTEST.py
│   ├── plot_fjpt_toptagTEST.py / .ipynb
│   ├── plot_fjpt_class_toptagTEST.py
│   ├── plot_fjpt_agg_toptagTEST.py
│   ├── plot_class_counts_toptagTEST.py
│   ├── plot_stack_toptag.py            # earlier 2023_0702 work
│   └── temp_plotting_2022EE_toptag*.py # earlier 2022EE binning work
└── reports/                # write-ups, one per study
    ├── report_toptag_effSF_alpha_20260928.md   # CURRENT α/eff/SF results, both eras
    ├── report_toptag_effSF_2023preBPix.txt     # 2026-09-23, 2023 (inclusive data bin)
    ├── report_alpha_toptagTEST.txt             # 2026-09-21/23, 2022 two-file setup
    ├── report_toptag_eff_toptagTEST.txt
    ├── report_toptag_eff_data_toptagTEST.txt
    ├── report_fjpt_inclusive_toptagTEST.txt
    ├── report_fjpt_toptagTEST.txt
    └── log_toptag_binning_2022EE.md
```

### `studies/` and `reports/` — a different environment

**These do NOT run on the hep server and do NOT use `~/envs/p38`.** They are
local analysis scripts that read the already-produced `.scaled` histogram
files (LZ4 cloudpickle) on the laptop, using the `plot_env` conda environment:

```bash
/Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python toptag/studies/<script>.py
```

The system default Python 3.7 **cannot** read those files (pickle protocol 5).

Each script is standalone, with a CONFIG block at the top and a
`JOBS`/`INPUT_FILES` dict driving the loop, and resolves its own paths
(`WORK_DIR` near the top of the file), so it can be run from any directory.
There are two groups:

- **α / efficiency scripts** (`plot_alpha_toptagTEST.py`,
  `plot_toptag_eff_toptagTEST.py`, `plot_toptag_eff_data_toptagTEST.py`): take
  the era as an optional argument (`2022preEE` default, or `2023preBPix`),
  read the `.scaled` inputs from `toptag/` (one level up) and write
  `toptag/plots_effSF/`. The current versions are the copies in
  `tmp_effSF_20260928/` (see [Current α / Efficiency / SF Workflow](#current-α--efficiency--sf-workflow-local-scaled-files));
  the `studies/` originals still expect the old two-file 2022 inputs.
- **All other scripts** (fjpt, class counts, stack, 2022EE binning): no
  arguments; read the `.scaled` inputs from, and write their `plots_*/`
  directories to, the analysis working directory two levels up, NOT `toptag/`.

Every report carries a `HOW TO RUN` section with the exact command, what it
produces, and any ordering constraint. Start with `reports/` to find out what
was measured and what the numbers mean — several carry CAVEATS sections
recording results that are known to be wrong and why.

**Status**: all five files (`classifiers.py`, `histmaker.py`, `scalefactors.py`,
`weights.py`, `validate.py`) are implemented. None have run against real
NanoAOD yet — this dev machine can't run `coffea.nanoevents`/`dask_awkward`/
`uproot` at all (see each file's doc below for exactly what *was* checked
locally: `scalefactors.py`'s math and correctionlib output were verified
against synthetic histograms with injected truth values). Before trusting
real numbers, run on the hep server in order: `python classifiers.py` →
`python histmaker.py --year <year> ...` → `python scalefactors.py --year
<year> ...`, and sanity-check the printed SF table against the AN's ballpark
(QCD-jet SF ~2-2.5, Top-jet SF ~0.8-1.0). Per-file docs, with exactly what's
verified vs. not and what to check first:

- [`docs/classifiers.md`](docs/classifiers.md) — AK15 → Top/B/QCD truth classification
- [`docs/histmaker.md`](docs/histmaker.md) — Step 1, histogram schema, known caveats (metadata mapping, cwd handling)
- [`docs/scalefactors.md`](docs/scalefactors.md) — Step 2, formulas as implemented, correctionlib schema reference, uncertainty scope
- [`docs/weights_and_validate.md`](docs/weights_and_validate.md) — Steps 3-4

---

## Environment

Two environments are used:

- **hep server** (Stages 1 and 3: `roc_curve.py`, `plot_TvsQCD.py`,
  `histmaker.py` → `weights.py`):

```bash
source ~/envs/p38/bin/activate   # Python 3.8
```

- **laptop** (Stage 2 and `studies/`: everything that reads `.scaled` files):
  `/Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python` (Python 3.10;
  `coffea`, `hist`, `mplhep`, `matplotlib`, `numpy`, `lz4`, `cloudpickle`; no
  `uproot`). The system Python 3.7 cannot read the `.scaled` files.

Required server packages (all available in the p38 environment, matching `decaf/`):

| Package | Purpose |
|---------|---------|
| `coffea` | Columnar HEP analysis, histogram objects |
| `awkward` | Jagged array operations |
| `correctionlib` | SF storage and evaluation |
| `uproot` | ROOT file I/O |
| `hist` | Histogram manipulation |
| `numpy` | Numerical operations |
| `matplotlib` | Validation plots |

---

## Reference

- **AN**: CMS AN2020/061 — "Search for new physics in final states with large
  missing transverse momentum and a top quark in proton-proton collisions at
  13 TeV" (Section 6.6: Fatjet top-tagging)
- **Working point**: the AN's 2018 (Run 2) calibration used `TvsQCD > 0.26`;
  this Run 3 codebase re-measured and uses `TvsQCD > 0.33`
  (`histmaker.py`'s `TVSQCD_WP`, matching `dev_run3.py:74-79`) — everywhere
  else in this README/code, 0.33 is the number that's actually applied.
- **Jet collection**: AK15 PFPuppi fatjets
