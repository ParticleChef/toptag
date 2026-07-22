# classifiers.py

Truth-level AK15 fatjet classification (Top-jet / B-jet / QCD-jet) used by
`histmaker.py` to fill the per-class histograms that `scalefactors.py` turns
into efficiencies and scale factors. MC only — never call on data.

Reference: `docs/toptag_AN2020_061.pdf` Sec 6.6; scheme summarized in the
root `README.md` under "AK15 Jet Classification".

## What it does

For each AK15 jet, two truth-level signals are combined:

1. **How many quarks from a hadronically-decaying top land inside the jet.**
   A hadronic top decays to `b + q q'` (3 quarks). We walk
   `GenPart -> distinctChildren` from each `fromHardProcess, isLastCopy` top
   to find its `b` and `W`, then `W -> distinctChildren` filtered to
   `|pdgId| <= 4` to get the two light/charm quarks — a `W` with fewer than 2
   such children decayed leptonically and is dropped. The surviving b/W-quark
   lists are flat per event (summed over both tops in an all-hadronic ttbar
   event, since a jet can occasionally straddle both). A jet counts a quark as
   "matched" if `fj.metric_table(quark) <= 1.5` (AK15 cone size).

2. **Whether a B hadron is ghost-matched into the jet** — `fj.nBHadrons > 0`,
   read directly from the ntuple (already computed at NanoAOD level).

These combine into the AN's 8-class scheme, then collapse to 3 classes:

| n matched quarks | no b-flavor | b-flavor found |
|---|---|---|
| 0 | 1 → **qcd** | 2 → **b** |
| 1 | 3 → **qcd** | 4 → **b** |
| 2 | 5 → **top** | 6 → **top** |
| 3 | 7 → **top** | — |
| >3 | 8 → **top** | — |

(≥2 matched quarks is always "top" regardless of b-flavor, per the AN table.)

## API

```python
from classifiers import classify_ak15_jets, aggregate_class

class8, class3 = classify_ak15_jets(fj, genpart)
# class8: int array [1..8], shaped like fj
# class3: string array "top"/"b"/"qcd", shaped like fj
```

`fj` — any AK15PuppiJet-shaped collection, e.g. the leading clean fatjet kept
as `[event, 1]` (jagged, not flattened — `metric_table` needs the jet axis).
`genpart` — `events.GenPart`, loaded with `CustomNanoAODSchema`
(`analysis/libs/mycoffeav2.py`) so `genPartIdxMother`-based navigation works.

Also exported:
- `count_matched_top_quarks(fj, genpart)` — just the quark count, if you need
  it without the class labels.
- `get_hadronic_top_decay_quarks(genpart)` — the raw `(b_quarks, w_quarks)`
  gen-level lists, if you want to reuse the matching for something else.
- `aggregate_class(class8)` — remaps a saved int `class8` array back to
  `top`/`b`/`qcd`, for reading back histograms that store the 8-class axis.

## Status: implemented, not yet verified against real data

Branch/collection names (`AK15PuppiJet`, `ParT_probTopbWqq` /
`ParT_probQCD*`, `nBHadrons`) are confirmed against existing, working code in
`analysis/processors/hadmonotop_run3.py` and `analysis/libs/mycoffeav2.py` —
not guessed. The one part with no precedent elsewhere in this repo is the
`GenPart.distinctChildren` / `.distinctParent` walk (standard coffea
NanoAODSchema behavior, but untested against these specific private
NanoTuples).

This machine's conda envs only have coffea 0.6.37 / awkward 0.12 (pre-
`coffea.nanoevents`, pre-awkward2), so the self-test below cannot be run
locally — it needs the `~/envs/p38` stack on the hep server.

## How to verify

```bash
source ~/envs/p38/bin/activate
cd /home/jhong/run3Monotop/toptag   # or wherever this repo lives on the server
python classifiers.py
```

This opens the local smoke-test file
`samples/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/nano_109.root` (all-hadronic
ttbar), classifies the leading good AK15 jet per event, and prints the
class8/class3 breakdown.

**What to check:**
- It runs without a `KeyError` on `AK15PuppiJet`, `GenPart`, `nBHadrons`, or
  `metric_table` — confirms the branch names and the schema-based GenPart
  navigation both work on this ntuple production.
- The `top` count in the class3 breakdown is clearly nonzero — TTto4Q is
  all-hadronic, so a real fraction of leading fatjets should have ≥2 matched
  quarks. All-zero `top` (with GJ/JetMET plausibly all-`qcd`) would point to
  a broken match (wrong cone size, broken parent/child walk, or `nBHadrons`
  not populated) rather than physics.
- Spot-check that `class8` classes 7 and 8 are rare relative to 5/6 — a fully
  merged 3-quark match should dominate over the >3-quark case, which only
  fires when a jet catches quarks from both hadronic tops in the event.

Once this looks sane, `histmaker.py` can call `classify_ak15_jets()` on the
leading fatjet in the tt̄ e/μ CRs and the γ+jets CR to fill the per-class
histograms the README's SF derivation needs.
