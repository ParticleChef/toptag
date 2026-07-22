# Top-Tagging Weight Calculation for CMS Run 3 Mono-top Analysis

## Overview

This project implements the **top-tagging scale factor (SF) and event weight** calculation
for the CMS mono-top analysis, adapting the Run 2 methodology documented in
`docs/toptag_AN2020_061.pdf` (Section 6.6) to the Run 3 framework used in `decaf/`.

The goal is to derive per-event reweighting factors that correct MC simulation to match
data for the `particleTransformer` top-tagger applied to large-radius (AK15) jets.

Top-tagging working point is measured: 0.33

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
Top-jet tag efficiency scale factor.

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

## Systematic Uncertainties

### QCD-jet mistagging SF (dominant sources)

| Source | 250–325 GeV | 325–400 GeV | 400–600 GeV | 600+ GeV |
|--------|-------------|-------------|-------------|----------|
| statistics | 15.1% | 6.2% | 4.1% | 6.4% |
| final-state PS | 0.7% | 5.0% | 7.8% | 13.3% |
| BJetMistag | 4.6% | 5.8% | 5.7% | 6.1% |
| BJetNorm | 2.7% | 3.9% | 3.8% | 3.9% |
| **total** | **18.1%** | **12.2%** | **12.1%** | **16.7%** |

### Top-jet tagging SF (dominant sources)

| Source | 250–325 GeV | 325–400 GeV | 400–600 GeV | 600+ GeV |
|--------|-------------|-------------|-------------|----------|
| statistics | 9.6% | 5.4% | 2.4% | 6.7% |
| TopPt reweighting | 5.6% | 6.1% | 6.7% | 9.4% |
| QCDMistag | 5.1% | 2.3% | 1.4% | 5.0% |
| BJetMistag | 7.0% | 4.3% | 3.4% | 7.6% |
| **total** | **16.4%** | **11.5%** | **8.8%** | **15.5%** |

---

## Design Rules

These rules ensure every step can be reproduced by anyone on the team,
without relying on any prior conversation or external context.

1. **Every script is standalone and runnable from the command line.**
   No step requires opening a notebook or calling internal functions manually.
   Each script uses `argparse` and prints its inputs, outputs, and progress.

2. **All configurable parameters live at the top of each script** in a clearly
   labeled `CONFIG` block. No magic numbers buried in logic.

3. **Each step reads from files and writes to files.** Intermediate results
   (histograms, SFs) are saved to disk so any step can be re-run independently.

4. **No implicit ordering.** The README lists the exact commands in order.
   Running them top-to-bottom from a clean checkout must always work.

5. **No `decaf/` code is modified.** This project only reads from `decaf/`
   as a reference. All new code lives under `toptag/`.

6. **Output files are versioned by year** (e.g. `toptag_sf_2022pre.json.gz`).
   Never overwrite without an explicit `--year` argument.

---

## Computing Environment

This is only for hep server.

This server has no batch/condor nodes. All jobs run locally.

| Resource | Value |
|----------|-------|
| CPU | AMD EPYC 9354 × 2 = **128 cores** |
| RAM | **251 GB** total (~239 GB available) |
| Local disk (`/`) | 1.8 TB NVMe, ~277 GB free |
| Home (`/home`) | NFS 2.3 TB, ~1.2 TB free |
| Input files | **Local** at `/data/mc/privatemc/<year>/...` |

The `2022_private_v1.json.gz` metadata contains ~2,994 datasets and ~29,643 files
in total. For the toptag calibration we only need:
- γ+jets CR: ~287 dataset chunks
- tt̄ CR: ~435 dataset chunks

### Recommended: coffea `futures_executor` with 40 workers (default)

Because all input files are on local NVMe (not XROOTD), the bottleneck is CPU,
not I/O. Running coffea's `futures_executor` with 40 workers in a single process
is the recommended approach for this subset of datasets:

- Uses ~40 cores out of 128 — leaves ~88 cores free for other users on the server
- Each worker uses ~2–3 GB RAM → ~80–120 GB total, well within the 239 GB available
- No monitoring loop or per-dataset log management needed
- If it crashes, coffea `.futures` checkpoint files allow resuming

**Always check current load before choosing a worker count:**
```bash
htop        # interactive view — press q to quit
# or a quick snapshot:
uptime      # load average over 1/5/15 min; keep load < ~100 on this machine
free -h     # check available RAM
```

Rule of thumb for `--workers N`:
- Server is quiet (load < 20):  `--workers 60`
- Normal use (load 20–60):      `--workers 40`  ← default
- Server is busy (load > 60):   `--workers 20`

Contrast with `nohup_job_new.py` (used for the full analysis): that runs 50 concurrent
single-worker processes, which is fine for 2,994 datasets but unnecessarily complex
for the ~700 datasets needed here.

---

## Running Commands (Step by Step)

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
│   └── weights_and_validate.md         # weights.py / validate.py details
├── analysis/                           # Main analysis framework (READ ONLY reference)
│   ├── libs/mycoffeav2.py              # CustomNanoAODSchema (AK15 branch aliasing)
│   ├── utils/ids.py                    # isGoodAK15 and other object ID helpers
│   ├── utils/corrections.py            # Pattern for SF/weight functions
│   ├── utils/common.py                 # BTag WPs, helpers
│   └── data/                           # Existing SF data files
│
├── classifiers.py          # [done, unverified] AK15 jet → Top/B/QCD class (GenPart matching)
├── histmaker.py            # [done, unverified] Step 1: fill CR histograms → .coffea file
├── scalefactors.py         # [done, verified with synthetic histograms] Step 2: histograms → SF JSON (correctionlib)
├── weights.py              # [done, unverified] Step 4: get_toptag_weight() for use in analysis
├── validate.py             # [done, unverified] Step 3: Data/MC comparison plots
├── samples/                 # local smoke-test ROOT files (one per dataset in
│                             #   2022_private_vToptag.json.gz)
│
├── hists/                  # Step 1 outputs (git-ignored, large)
│   └── toptag_hists_<year>.coffea
├── data/                   # Step 2 outputs — SF files (committed)
│   └── toptag_sf_<year>.json.gz
└── plots/                  # Step 3 outputs (git-ignored)
    └── validate_<year>/
```

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

```bash
source ~/envs/p38/bin/activate   # Python 3.8
```

Required packages (all available in the p38 environment, matching `decaf/`):

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
