#!/usr/bin/env python
"""weights.py - Step 4 of the top-tagging SF pipeline (see README.md).

Applies the promotion/demotion top-tagging event weight, using the
correctionlib file produced by scalefactors.py:

    tagged jet:   w = SF
    untagged jet: w = (1 - eff*SF) / (1 - eff)

per README.md "Event Weight" (analogous to the standard AK4 b-tagging
promotion/demotion method).
"""

import functools
import os

import awkward as ak
import numpy as np

TOPTAG_ROOT = os.path.dirname(os.path.abspath(__file__))


@functools.lru_cache(maxsize=None)
def _correction_set(year):
    import correctionlib

    path = os.path.join(TOPTAG_ROOT, "data", f"toptag_sf_{year}.json.gz")
    return correctionlib.CorrectionSet.from_file(path)


def get_toptag_weight(year, jet_pt, jet_class, is_tagged):
    """Per-event top-tagging weight.

    year      - e.g. "2022pre"
    jet_pt    - array-like, leading AK15 jet pT [GeV]
    jet_class - array-like of "top"/"b"/"qcd" (from classifiers.py, MC truth)
    is_tagged - array-like of bool, TvsQCD > working point

    Returns (nominal, up, down) arrays, same shape as the inputs.
    """
    cset = _correction_set(year)
    sf_corr = cset["toptag_sf"]
    eff_corr = cset["toptag_efficiency"]

    def to_numpy(arr):
        return ak.to_numpy(arr) if isinstance(arr, ak.Array) else np.asarray(arr)

    pt = to_numpy(jet_pt)
    jclass = to_numpy(jet_class)
    tagged = to_numpy(is_tagged)

    def weight_for(systematic):
        sf = np.array([sf_corr.evaluate(p, c, systematic) for p, c in zip(pt, jclass)])
        eff = np.array([eff_corr.evaluate(p, c) for p, c in zip(pt, jclass)])
        return np.where(tagged, sf, (1 - eff * sf) / (1 - eff))

    nominal = weight_for("nominal")
    up = weight_for("up")
    down = weight_for("down")
    return nominal, up, down
