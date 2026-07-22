# histmaker.py

Step 1 of the top-tagging SF pipeline (see root `README.md`). Runs a coffea
processor over NanoAOD in the `gcr` (γ+jets), `tmcr` (ttbar μ), and `tecr`
(ttbar e) control regions and fills one histogram per region:
`fjpt` (4 AN pT bins) × `TvsQCD` (untagged/tagged at WP 0.33) × `jetclass`
(`"top"`/`"b"`/`"qcd"` for MC via `classifiers.py`, `"data"` for data).

Reference: `docs/toptag_AN2020_061.pdf` Sec 6.6.1.3 (control regions, pT
bins); object IDs/selection/weights ported from
`analysis/processors/dev_run3.py` (read-only reference, not modified).

## What it does

1. **Object selection & weights** — muons/electrons/photons/AK15/AK4 jet
   definitions, `gcr`/`tmcr`/`tecr` region cuts, and the standard per-event
   MC weights (`genw`, `pileup`, `nlo_ewk`, `topptreweighting`, `btagSF`,
   `ids`, `reco`, `iso`, `trig`) are ported verbatim from `dev_run3.py`,
   trimmed to drop everything specific to `sr`/`wmcr`/`wecr`/`zmcr`/`zecr`.
   Systematic weight *variations* (btag up/down, etc.) are **not** carried —
   only the nominal per-event weight, matching the "stat + B-jet assumption
   only" scope for this calibration.
2. **Truth jet class** — for MC events, calls
   `classifiers.classify_ak15_jets(leading_fj_jagged, events.GenPart)` on the
   leading clean/good AK15 jet and fills the 4th histogram axis with the
   result (`"none"` if there's no leading jet in that event — harmless,
   never selected downstream). Data rows get `"data"`.
3. **In-processor xs·lumi scaling**, exactly like `dev_run3.py:1166-1173`.
4. **Driver loop** (`main()`) — filters `analysis/metadata/<file>.json.gz`
   down to datasets used by `gcr`/`tmcr`/`tecr` (`select_datasets`, matching
   `dev_run3.py`'s `self._samples` substrings), runs one
   `coffea.processor.run_uproot_job` per dataset (same pattern as
   `analysis/run.py`), then `group_and_normalize()` applies the `1/sumw` MC
   normalization and the `bkg_map`/`data_map` process grouping ported
   verbatim from `analysis/macros/scale.py:84-160`. Output:
   `{'bkg': {process: hist}, 'data': {process: hist}}` per region — the same
   shape convention as a `.scaled` file (see root `CLAUDE.md`), just
   restricted to these 3 regions and with the new axes.

## Known caveats / things to double-check before trusting a run

- **`METADATA_FILES`** (top of the file) is a *best guess* at which
  `analysis/metadata/*.json.gz` file covers each year/era, picked from what's
  present in this checkout — it was never confirmed against what
  `analysis/run.py` normally uses in production. Check before running.
- **cwd handling**: `ids.coffea`/`corrections.coffea`/`common.coffea` (and
  the correctionlib files they lazily open, e.g. `isGoodAK4`'s jet-ID file)
  resolve paths relative to `analysis/`, because `analysis/run.py` is always
  run from inside that directory. `histmaker.py` lives at the toptag root, so
  `main()` resolves the `--output` path to absolute *first*, then
  `os.chdir()`s into `analysis/` before loading anything. If you see a
  `FileNotFoundError` for something under `data/...`, this is the first
  place to look.
- **`PT_EDGES`'s last bin** uses `3000` as a stand-in for "600-∞ GeV" — fine
  for filling/deriving SFs, but note it if you ever read the histogram's raw
  bin edges directly.

## Status: implemented, not yet run against real data

This dev machine has no `coffea.nanoevents` stack at all (confirmed:
`plot_env` is missing `awkward`; no local `p38`/CMSSW env), so this file has
only been syntax-checked (`python -m py_compile histmaker.py`), not executed.

## How to run / verify

```bash
source ~/envs/p38/bin/activate
cd /home/jhong/run3Monotop/toptag   # or wherever this repo lives on the server
python histmaker.py --year 2022pre --workers 40 --output hists/toptag_hists_2022pre.coffea
```

**What to check first**, before trusting the output:
- It runs without a `FileNotFoundError`/`KeyError` — confirms the cwd
  handling and `METADATA_FILES` mapping above are correct for this checkout.
- The printed `"N datasets selected for gcr/tmcr/tecr"` count is in the
  right ballpark (README's own estimate: γ+jets ~287 chunks, ttbar ~435).
- Load the output and spot-check yields make sense before running
  `scalefactors.py` on it:
  ```python
  from coffea.util import load
  h = load("hists/toptag_hists_2022pre.coffea")
  print(h["bkg"].keys())   # should include "G + Jets", "TT", ...
  print(h["data"].keys())  # should include "MET", "EGamma"
  print(h["bkg"]["TT"][{"region": "tmcr"}].project("fjpt", "TvsQCD").values())
  ```
  Sanity checks: `"G + Jets"` should dominate `gcr`'s `"qcd"` jetclass;
  `"TT"` should dominate `tmcr`/`tecr`'s `"top"` jetclass (see AN Fig. 22 for
  the expected process composition per region).
