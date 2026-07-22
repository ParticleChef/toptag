# scalefactors.py

Step 2 of the top-tagging SF pipeline (see root `README.md`). Reads the
histograms from `histmaker.py` and derives the QCD-jet and Top-jet data/MC
scale factors per `docs/toptag_AN2020_061.pdf` Sec 6.6.1.3, writing the
result as a `correctionlib` `CorrectionSet` (`.json.gz`).

## What it does

### Formulas implemented (per pT bin)

- **α** = N(Data)/N(MC), pre-tag, over *all* processes/jetclasses in the
  relevant region (`gcr` for QCD-jet, `tmcr`+`tecr` combined for Top-jet).
- **QCD-jet** (from `gcr`, using the `"G + Jets"` bkg process specifically —
  the AN's γ+jets-dominance approximation):
  `SF_QCD = [N(Data,tagged) − α·N(GJ,B,tagged)] / [α·N(GJ,QCD,tagged)]`
- **Top-jet** (from `tmcr`+`tecr`, `MC` = all bkg processes summed):
  `SF_Top = [N(Data,tagged) − α·(SF_B·N(MC,B,tagged) + SF_QCD·N(MC,QCD,tagged))] / [α·N(MC,Top,tagged)]`

  **Note the `SF_QCD·N(MC,QCD,tagged)` term** — the AN's boxed formula
  (and the root README's literal transcription of it) shows this
  subtraction using the raw MC yield, but the surrounding AN prose says the
  previously-derived QCD-jet SF is used to *correct* this term, and AN
  Table 32 lists a nonzero "QCDMistag" uncertainty that would be
  meaningless if the central formula didn't depend on `SF_QCD` at all. This
  implementation applies the correction (`sf_qcd_` argument in
  `derive_top_sf`'s `formula()`); if you disagree with this reading of the
  AN, that's the one line to change.
- **B-jet**: fixed `SF_B = 1 ± 0.5` (`SF_B_CENTRAL`/`SF_B_UNCERTAINTY`), no CR.

### Uncertainty propagation (scope: statistics + B-jet assumption only)

`propagate()` is a generic finite-difference error propagator: given a
formula function and a list of `(nominal, sigma)` input pairs, it shifts
each input by its own σ (holding the rest fixed) and adds the resulting
shifts in quadrature — the standard linearized-uncertainty estimate. Inputs
are treated as uncorrelated, which is an approximation for terms drawn from
overlapping samples (e.g. a tagged yield vs. its own total); called out in
the docstring rather than hidden. Sources actually propagated, matching AN
Tables 31/32's row names:
- **statistics** — each histogram term's own `sqrt(variance)`.
- **BJetNorm** — 50% relative uncertainty on the γ+jets B-jet MC yield,
  folded into that term's effective σ (`b_tag_sigma` in `derive_qcd_sf`).
- **BJetMistag** — the assumed `SF_B = 1 ± 0.5`, passed as its own
  `(nominal, sigma)` pair into both `derive_qcd_sf` and `derive_top_sf`.
- **QCDMistag** (Top-jet SF only) — the previously-derived `(sf_qcd,
  sf_qcd_unc)` passed straight into `derive_top_sf`'s `propagate()` call.

**Not propagated** (out of scope for this version, per AN Tables 31/32):
JES (`jesRelativeBal`, `jesRelativeSample`, `jesFlavorQCD`), QCD scale/PDF
variations, parton-shower variations, `TopPt_reweighting` uncertainty. These
would need dedicated systematic-varied histograms from `histmaker.py` that
don't currently exist.

### correctionlib schema

Two `Correction`s in one `CorrectionSet` (`schema_version=2`):

- **`toptag_sf`** — inputs `pt` (real), `jetclass` (string: `"top"`/`"b"`/
  `"qcd"`), `systematic` (string: `"nominal"`/`"up"`/`"down"`); output `sf`
  (real). Binned in `pt` per `PT_EDGES = [250,325,400,600,3000]` (last edge
  is a stand-in for "∞"), `flow="clamp"` — pT outside range clamps to the
  nearest bin rather than erroring.
- **`toptag_efficiency`** — inputs `pt`, `jetclass`; output `efficiency`
  (real), no `systematic` input (only the nominal MC efficiency is needed —
  see `weights.py`'s promotion/demotion formula, which only takes the
  uncertainty through `sf`, not `efficiency`). `"b"`-class efficiency is
  read directly from MC truth in the combined ttbar CR
  (`derive_b_efficiency`) — there's no CR-measured *scale factor* for B-jet,
  but the MC efficiency itself is directly computable from simulation.

Evaluate like any correctionlib file:
```python
import correctionlib
cset = correctionlib.CorrectionSet.from_file("data/toptag_sf_2022pre.json.gz")
sf = cset["toptag_sf"].evaluate(pt, "top", "nominal")
eff = cset["toptag_efficiency"].evaluate(pt, "top")
```

## Status: implemented, verified against synthetic histograms (not real data)

This dev machine can't run `coffea.nanoevents`, so there's no real
`hists/toptag_hists_<year>.coffea` to test against. What *was* done: built
synthetic `hist.Hist` objects matching `histmaker.py`'s exact schema, with
injected truth `SF_QCD`, `SF_Top`, and efficiency values (including
realistic γ+jets B-jet contamination and TT QCD/B-jet contamination), ran
`derive_qcd_sf`/`derive_top_sf`/`derive_b_efficiency`/`build_correction_set`
directly, and confirmed:
- the recovered SF/efficiency arrays match the injected truth to `rtol` 2-5%
  (only imprecise because of a coarse synthetic-fill discretization, not a
  formula bug),
- the written `.json.gz` loads back via
  `correctionlib.CorrectionSet.from_file(...)` and `.evaluate(...)` returns
  the same numbers, with `up > nominal > down` as expected,
  and out-of-range `pt` correctly clamps to the edge bin's value.

This confirms the **math and the correctionlib construction/round-trip are
correct** — it does not confirm anything about `histmaker.py`'s physics
selection or real yields, which can only be checked on the hep server.

## How to run / verify

```bash
source ~/envs/p38/bin/activate
cd /home/jhong/run3Monotop/toptag
python scalefactors.py --year 2022pre \
    --hists hists/toptag_hists_2022pre.coffea \
    --output data/toptag_sf_2022pre.json.gz
```

**What to check**: the printed SF table's ballpark against AN Sec 6.6.1.3
"Results" — QCD-jet SF roughly 2-2.5, Top-jet SF roughly 0.8-1.0, both
decreasing/flattening slightly with pT. A QCD-jet SF near 1 would suggest
something's wrong (the whole point of this calibration is that the γ+jets
mistag rate is *not* well modeled in simulation, per the AN).
