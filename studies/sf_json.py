"""correctionlib output for the top-tag SFs from plot_toptag_eff_data_toptagTEST.py.

Writes a CorrectionSet with the same layout as toptag/scalefactors.py, so
toptag/weights.py can read it unchanged:

  toptag_sf          (pt, jetclass, systematic) -> SF
                     jetclass   "top" | "qcd" | "b"
                     systematic "nominal" | "up" | "down"  (up/down = +-1 sigma)
  toptag_efficiency  (pt, jetclass) -> MC tag efficiency

Both are binned in the leading AK15 jet pT, flow="clamp" (pT outside the
edges takes the first/last bin).
"""

import json
import os

import numpy as np

SCHEMA_VERSION = 2


def sanitize(values, errors, label, fallback=(1.0, 1.0)):
    """Replace non-finite (value, error) bins by `fallback` and warn.

    correctionlib needs finite numbers; empty bins (e.g. low-statistics
    eras) would otherwise give NaN/inf.
    """
    values, errors = np.array(values, float), np.array(errors, float)
    bad = ~(np.isfinite(values) & np.isfinite(errors))
    if bad.any():
        print(f"  WARNING: {label}: non-finite bins {np.flatnonzero(bad).tolist()}"
              f" set to {fallback[0]} +- {fallback[1]}")
        values[bad], errors[bad] = fallback
    return values, errors


def binning(edges, content):
    """correctionlib Binning node in pt."""
    return {"nodetype": "binning", "input": "pt",
            "edges": [float(e) for e in edges],
            "content": [float(c) for c in content], "flow": "clamp"}


def category(var, items):
    """correctionlib Category node from a {key: node} dict."""
    return {"nodetype": "category", "input": var,
            "content": [{"key": k, "value": v} for k, v in items.items()]}


def sf_with_syst(edges, sf, err):
    """nominal / up / down Binning nodes; down is clipped at 0."""
    return category("systematic", {
        "nominal": binning(edges, sf),
        "up": binning(edges, sf + err),
        "down": binning(edges, np.clip(sf - err, 0, None)),
    })


def build_correction_set(era, edges, sfs, effs, wp):
    """CorrectionSet dict.

    sfs:  {jetclass: (sf, sf_err)} per pT bin
    effs: {jetclass: eff_mc}       per pT bin
    """
    pt = {"name": "pt", "type": "real", "description": "Leading AK15 jet pT [GeV]"}
    jc = {"name": "jetclass", "type": "string", "description": '"top", "qcd" or "b"'}
    sf_corr = {
        "name": "toptag_sf", "version": 1,
        "description": (f"ParT TvsQCD > {wp} data/MC scale factor for the leading"
                        f" AK15 jet, {era} (AN-2020/061 Sec. 6.6.1, Run 3)."
                        " up/down = total (stat + SF_B + B-jet norm.) uncertainty."),
        "inputs": [pt, jc, {"name": "systematic", "type": "string",
                            "description": '"nominal", "up" or "down"'}],
        "output": {"name": "sf", "type": "real"},
        "data": category("jetclass", {k: sf_with_syst(edges, *v) for k, v in sfs.items()}),
    }
    eff_corr = {
        "name": "toptag_efficiency", "version": 1,
        "description": f"MC top-tag (TvsQCD > {wp}) efficiency of the leading AK15 jet, {era}.",
        "inputs": [pt, jc],
        "output": {"name": "efficiency", "type": "real"},
        "data": category("jetclass", {k: binning(edges, v) for k, v in effs.items()}),
    }
    return {"schema_version": SCHEMA_VERSION, "corrections": [sf_corr, eff_corr]}


def validate(cset):
    """Check the dict against the correctionlib schema; return the evaluator."""
    import correctionlib
    import correctionlib.schemav2 as cs
    model = cs.CorrectionSet.model_validate(cset)
    return correctionlib.CorrectionSet.from_string(model.model_dump_json())


def write(cset, path):
    """Validate and write the CorrectionSet as (plain) JSON."""
    validate(cset)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as f:
        json.dump(cset, f, indent=1)
    print(f"Saved {path}")
