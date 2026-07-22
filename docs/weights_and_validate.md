# weights.py / validate.py

Steps 3 (optional) and 4 of the top-tagging SF pipeline (see root
`README.md`). Both are downstream consumers of `scalefactors.py`'s
`data/toptag_sf_<year>.json.gz`.

## weights.py

`get_toptag_weight(year, jet_pt, jet_class, is_tagged)` — the promotion/
demotion event weight from the root README's "Event Weight" section:

```
tagged jet:   w = SF
untagged jet: w = (1 - eff*SF) / (1 - eff)
```

evaluated per-event from the `toptag_sf`/`toptag_efficiency` corrections in
`data/toptag_sf_<year>.json.gz` (`correctionlib.CorrectionSet.from_file`,
cached per year via `functools.lru_cache` so repeated calls in the same
process don't re-read the file). Returns `(nominal, up, down)` arrays.

**Note**: `sf_corr.evaluate(...)`/`eff_corr.evaluate(...)` are called in a
plain Python loop over events, not vectorized — correctionlib's nested
`Category`/`Binning` structure doesn't have a batch-evaluate path in the
version used here. Fine for the sizes this calibration deals with; if this
becomes a bottleneck in the main analysis, that loop is the place to
optimize (e.g. group by `jet_class` first, since correctionlib's binned
`evaluate` degenerates to array ops once the categorical branch is fixed).

`jet_pt`/`jet_class`/`is_tagged` are expected already flattened to one value
per event (e.g. `fatjet.pt[:, 0]`), matching the README's own usage example
— not jagged arrays.

## validate.py

Reproduces the AN's Fig. 29/30 style: two-bin (`[0,wp]`/`[wp,1]`) Data/MC
comparison, before vs. after applying the calibration, per region
(`gcr`/`tmcr`/`tecr`). `calibrated_mc()` reweights each jetclass's MC yield
using the same nominal-only formula as `weights.py`, applied directly to
histogram bin contents (no per-event loop needed since we're just scaling
bin values, evaluating `sf`/`eff` once per pT-bin center per jetclass).
Output: `plots/validate_<year>/toptag_validate_<region>_<year>.{pdf,png}`.

This is the least load-bearing of the four new files (README marks Step 3
"optional") — it's a visual sanity check, not part of the SF derivation
itself. Kept simpler than the full `docs/cms_plot_style.md` stack-plot
convention (no full process-by-process stack, just tagged/untagged data-vs-
MC); reuse `plot_stack.py`'s style if you want the full CMS look.

## Status: implemented, not yet run

Neither has been run against real data — both need `data/toptag_sf_<year>.
json.gz` from a real `scalefactors.py` run, which needs a real `histmaker.py`
run first (see `docs/histmaker.md`, `docs/scalefactors.md`). Both were only
syntax-checked (`python -m py_compile weights.py validate.py`) on this dev
machine.

## How to run / verify

```bash
source ~/envs/p38/bin/activate
cd /home/jhong/run3Monotop/toptag
python validate.py --year 2022pre \
    --hists hists/toptag_hists_2022pre.coffea \
    --sf data/toptag_sf_2022pre.json.gz \
    --outdir plots/validate_2022pre/
```

**What to check**: per AN Fig. 29/30, the "after calibration" panel's
Data/MC ratio in the tagged bin should move noticeably closer to 1.0 than
the "before calibration" panel — that's the whole point of the calibration,
and the most direct visual check that the derived SFs are doing their job.

For `weights.py`, a quick smoke test once a real SF file exists:
```python
from weights import get_toptag_weight
nominal, up, down = get_toptag_weight(
    year="2022pre", jet_pt=[300.0], jet_class=["top"], is_tagged=[True],
)
print(nominal, up, down)  # should equal SF_Top's [250,325) bin, +-unc
```
