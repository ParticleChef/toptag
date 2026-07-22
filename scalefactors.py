#!/usr/bin/env python
"""scalefactors.py - Step 2 of the top-tagging SF pipeline (see README.md).

Reads the region x fjpt x TvsQCD x jetclass histograms from histmaker.py and
derives the QCD-jet and Top-jet data/MC scale factors per AN-2020/061 Sec
6.6.1.3, propagating statistics plus the SF_B-jet=1+-0.5 assumption (and its
associated gamma+jets B-jet MC normalization uncertainty). JES/scale/parton-
shower systematics (AN Tables 31-32) are NOT propagated in this version.

Writes the result as a correctionlib CorrectionSet ("toptag_sf",
"toptag_efficiency") to a gzip-compressed JSON file, ready for
correctionlib.CorrectionSet.from_file(...).

Usage:
    python scalefactors.py --year 2022pre --hists hists/toptag_hists_2022pre.coffea \
        --output data/toptag_sf_2022pre.json.gz
"""

import argparse
import gzip
import os

import numpy as np
from coffea.util import load

# ============================================================
# CONFIG
# ============================================================

PT_EDGES = [250, 325, 400, 600, 3000]  # must match histmaker.py's PT_EDGES
PT_LABELS = ["250-325 GeV", "325-400 GeV", "400-600 GeV", "600+ GeV"]

GJ_PROCESS = "G + Jets"          # bkg process name used for the QCD-jet calibration
TTCR_DATA_PROCESS = "EGamma"     # data process common to tecr; tmcr uses "MET"

SF_B_CENTRAL = 1.0
SF_B_UNCERTAINTY = 0.5           # AN: no dedicated CR for the B-jet class
GJ_BJET_NORM_UNCERTAINTY = 0.5   # AN Table 31 "BJetNorm": 50% norm unc. on GJ B-jet MC


# ============================================================
# Histogram helpers
# ============================================================

def sum_processes(hist_dict, processes=None):
    """Sum hist.Hist objects across the given process keys (default: all)."""
    keys = processes if processes is not None else list(hist_dict)
    total = None
    for key in keys:
        if key not in hist_dict:
            continue
        total = hist_dict[key] if total is None else total + hist_dict[key]
    if total is None:
        raise KeyError(f"none of {keys} found in histogram (have {list(hist_dict)})")
    return total


def counts(h, region, jetclass=None):
    """(values, variances), each shape (n_pt_bins, 2) with [:,0]=untagged, [:,1]=tagged."""
    sliced = h[{"region": region}]
    if jetclass is not None:
        sliced = sliced[{"jetclass": jetclass}]
    sliced = sliced.project("fjpt", "TvsQCD")
    return sliced.values(), sliced.variances()


def tag_split(h, region, jetclass=None):
    """(tagged, tagged_var, total, total_var), each shape (n_pt_bins,)."""
    values, variances = counts(h, region, jetclass=jetclass)
    tagged, tagged_var = values[:, 1], variances[:, 1]
    total, total_var = values.sum(axis=1), variances.sum(axis=1)
    return tagged, tagged_var, total, total_var


# ============================================================
# Error propagation
# ============================================================

def propagate(formula, *args_with_unc):
    """Standard linearized error propagation via one-at-a-time finite differences:
    shift each input by its own 1-sigma uncertainty (holding the others at their
    central value), add the resulting shifts in quadrature. This treats inputs
    as uncorrelated, which is an approximation for terms drawn from overlapping
    event samples (e.g. a tagged yield vs. its own total) - acceptable for a
    first "stat + B-jet assumption" version, called out here rather than hidden.

    args_with_unc: each a (nominal, sigma) pair of same-shape numpy arrays.
    """
    nominals = [a[0] for a in args_with_unc]
    with np.errstate(divide="ignore", invalid="ignore"):
        central = formula(*nominals)
        var = np.zeros_like(np.asarray(central, dtype=float))
        for i in range(len(args_with_unc)):
            shifted = list(nominals)
            shifted[i] = nominals[i] + args_with_unc[i][1]
            var = var + (formula(*shifted) - central) ** 2
    return central, np.sqrt(var)


# ============================================================
# QCD-jet calibration (gamma+jets CR)
# ============================================================

def derive_qcd_sf(bkg, data):
    mc_total = sum_processes(bkg)
    data_total = sum_processes(data, processes=[TTCR_DATA_PROCESS] if TTCR_DATA_PROCESS in data else None)

    _, _, n_mc_tot, v_mc_tot = tag_split(mc_total, "gcr")
    d_tag, vd_tag, n_data_tot, v_data_tot = tag_split(data_total, "gcr")

    gj = bkg[GJ_PROCESS]
    q_tag, vq_tag, q_tot, _ = tag_split(gj, "gcr", jetclass="qcd")
    b_tag, vb_tag, _, _ = tag_split(gj, "gcr", jetclass="b")
    with np.errstate(divide="ignore", invalid="ignore"):
        eps_qcd_mc = q_tag / q_tot

    # fold the assumed 50% gamma+jets-B-jet-MC normalization uncertainty into
    # b_tag's effective sigma (AN Table 31 "BJetNorm")
    b_tag_sigma = np.sqrt(vb_tag + (GJ_BJET_NORM_UNCERTAINTY * b_tag) ** 2)

    def formula(d, b, q, n_data, n_mc, sf_b):
        alpha = n_data / n_mc
        return (d - alpha * sf_b * b) / (alpha * q)

    sf_qcd, sf_qcd_unc = propagate(
        formula,
        (d_tag, np.sqrt(vd_tag)),
        (b_tag, b_tag_sigma),
        (q_tag, np.sqrt(vq_tag)),
        (n_data_tot, np.sqrt(v_data_tot)),
        (n_mc_tot, np.sqrt(v_mc_tot)),
        (np.full_like(q_tag, SF_B_CENTRAL), np.full_like(q_tag, SF_B_UNCERTAINTY)),
    )
    return sf_qcd, sf_qcd_unc, eps_qcd_mc


# ============================================================
# Top-jet calibration (ttbar e+mu CR, combined)
# ============================================================

def derive_top_sf(bkg, data, sf_qcd, sf_qcd_unc):
    mc_all = sum_processes(bkg)
    data_tmcr = sum_processes(data, processes=["MET"] if "MET" in data else None)
    data_tecr = sum_processes(data, processes=[TTCR_DATA_PROCESS] if TTCR_DATA_PROCESS in data else None)

    _, _, n_mc_tot_m, v_mc_tot_m = tag_split(mc_all, "tmcr")
    _, _, n_mc_tot_e, v_mc_tot_e = tag_split(mc_all, "tecr")
    n_mc_tot, v_mc_tot = n_mc_tot_m + n_mc_tot_e, v_mc_tot_m + v_mc_tot_e

    d_tag_m, vd_tag_m, n_data_tot_m, v_data_tot_m = tag_split(data_tmcr, "tmcr")
    d_tag_e, vd_tag_e, n_data_tot_e, v_data_tot_e = tag_split(data_tecr, "tecr")
    d_tag, vd_tag = d_tag_m + d_tag_e, vd_tag_m + vd_tag_e
    n_data_tot, v_data_tot = n_data_tot_m + n_data_tot_e, v_data_tot_m + v_data_tot_e

    top_tag_m, vtop_m, top_tot_m, vtop_tot_m = tag_split(mc_all, "tmcr", jetclass="top")
    top_tag_e, vtop_e, top_tot_e, vtop_tot_e = tag_split(mc_all, "tecr", jetclass="top")
    top_tag, vtop_tag = top_tag_m + top_tag_e, vtop_m + vtop_e
    top_tot, vtop_tot = top_tot_m + top_tot_e, vtop_tot_m + vtop_tot_e
    eps_top_mc = top_tag / top_tot

    b_tag_m, vb_m, _, _ = tag_split(mc_all, "tmcr", jetclass="b")
    b_tag_e, vb_e, _, _ = tag_split(mc_all, "tecr", jetclass="b")
    b_tag, vb_tag = b_tag_m + b_tag_e, vb_m + vb_e

    qcd_tag_m, vq_m, _, _ = tag_split(mc_all, "tmcr", jetclass="qcd")
    qcd_tag_e, vq_e, _, _ = tag_split(mc_all, "tecr", jetclass="qcd")
    qcd_tag, vqcd_tag = qcd_tag_m + qcd_tag_e, vq_m + vq_e

    def formula(d, top, b, qcd, n_data, n_mc, sf_b, sf_qcd_):
        alpha = n_data / n_mc
        bqcd_corrected = sf_b * b + sf_qcd_ * qcd
        return (d - alpha * bqcd_corrected) / (alpha * top)

    sf_top, sf_top_unc = propagate(
        formula,
        (d_tag, np.sqrt(vd_tag)),
        (top_tag, np.sqrt(vtop_tag)),
        (b_tag, np.sqrt(vb_tag)),
        (qcd_tag, np.sqrt(vqcd_tag)),
        (n_data_tot, np.sqrt(v_data_tot)),
        (n_mc_tot, np.sqrt(v_mc_tot)),
        (np.full_like(qcd_tag, SF_B_CENTRAL), np.full_like(qcd_tag, SF_B_UNCERTAINTY)),  # BJetMistag
        (sf_qcd, sf_qcd_unc),  # QCDMistag: propagate the previously-derived SF_QCD
    )
    return sf_top, sf_top_unc, eps_top_mc


def derive_b_efficiency(bkg):
    """MC-truth B-jet-class tag efficiency, read directly from simulation
    (combined ttbar CR - same region used for the Top-jet calibration).
    No data-driven scale factor exists for this class (see SF_B_CENTRAL)."""
    mc = sum_processes(bkg)
    b_tag_m, _, b_tot_m, _ = tag_split(mc, "tmcr", jetclass="b")
    b_tag_e, _, b_tot_e, _ = tag_split(mc, "tecr", jetclass="b")
    b_tag, b_tot = b_tag_m + b_tag_e, b_tot_m + b_tot_e
    with np.errstate(divide="ignore", invalid="ignore"):
        return b_tag / b_tot


# ============================================================
# correctionlib output
# ============================================================

def build_correction_set(year, sf_qcd, sf_qcd_unc, sf_top, sf_top_unc, eps_qcd, eps_top, eps_b):
    import correctionlib.schemav2 as cs

    def binning(values):
        return cs.Binning(
            nodetype="binning",
            input="pt",
            edges=[float(e) for e in PT_EDGES],
            content=[float(v) for v in values],
            flow="clamp",
        )

    def sf_binning_with_syst(central, unc):
        return {
            "nominal": binning(central),
            "up": binning(central + unc),
            "down": binning(np.clip(central - unc, a_min=0, a_max=None)),
        }

    qcd_syst = sf_binning_with_syst(sf_qcd, sf_qcd_unc)
    top_syst = sf_binning_with_syst(sf_top, sf_top_unc)
    b_central = np.full(len(PT_EDGES) - 1, SF_B_CENTRAL)
    b_unc = np.full(len(PT_EDGES) - 1, SF_B_UNCERTAINTY)
    b_syst = sf_binning_with_syst(b_central, b_unc)

    sf_correction = cs.Correction(
        name="toptag_sf",
        description=(
            f"Top-tagging (TvsQCD>0.33) data/MC scale factor for the leading AK15 jet, "
            f"year {year}, per AN-2020/061 Sec 6.6.1.3 (adapted for Run 3)."
        ),
        version=1,
        inputs=[
            cs.Variable(name="pt", type="real", description="Leading AK15 jet pT [GeV]"),
            cs.Variable(name="jetclass", type="string", description='"top", "b", or "qcd"'),
            cs.Variable(name="systematic", type="string", description='"nominal", "up", or "down"'),
        ],
        output=cs.Variable(name="sf", type="real"),
        data=cs.Category(
            nodetype="category",
            input="jetclass",
            content=[
                cs.CategoryItem(key="top", value=cs.Category(
                    nodetype="category", input="systematic",
                    content=[cs.CategoryItem(key=k, value=v) for k, v in top_syst.items()],
                )),
                cs.CategoryItem(key="qcd", value=cs.Category(
                    nodetype="category", input="systematic",
                    content=[cs.CategoryItem(key=k, value=v) for k, v in qcd_syst.items()],
                )),
                cs.CategoryItem(key="b", value=cs.Category(
                    nodetype="category", input="systematic",
                    content=[cs.CategoryItem(key=k, value=v) for k, v in b_syst.items()],
                )),
            ],
        ),
    )

    eff_correction = cs.Correction(
        name="toptag_efficiency",
        description=f"MC top-tagging (TvsQCD>0.33) efficiency for the leading AK15 jet, year {year}.",
        version=1,
        inputs=[
            cs.Variable(name="pt", type="real", description="Leading AK15 jet pT [GeV]"),
            cs.Variable(name="jetclass", type="string", description='"top", "b", or "qcd"'),
        ],
        output=cs.Variable(name="efficiency", type="real"),
        data=cs.Category(
            nodetype="category",
            input="jetclass",
            content=[
                cs.CategoryItem(key="top", value=binning(eps_top)),
                cs.CategoryItem(key="qcd", value=binning(eps_qcd)),
                cs.CategoryItem(key="b", value=binning(eps_b)),
            ],
        ),
    )

    return cs.CorrectionSet(schema_version=2, corrections=[sf_correction, eff_correction])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", required=True)
    parser.add_argument("--hists", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    hists = load(args.hists)
    bkg, data = hists["bkg"], hists["data"]

    sf_qcd, sf_qcd_unc, eps_qcd = derive_qcd_sf(bkg, data)
    sf_top, sf_top_unc, eps_top = derive_top_sf(bkg, data, sf_qcd, sf_qcd_unc)
    eps_b = derive_b_efficiency(bkg)

    print(f"Year: {args.year}")
    print(f"{'pT bin':<14}{'SF_QCD':<18}{'SF_Top':<18}")
    for i, label in enumerate(PT_LABELS):
        print(f"{label:<14}{sf_qcd[i]:.2f} +- {sf_qcd_unc[i]:.2f}      {sf_top[i]:.2f} +- {sf_top_unc[i]:.2f}")

    cset = build_correction_set(args.year, sf_qcd, sf_qcd_unc, sf_top, sf_top_unc, eps_qcd, eps_top, eps_b)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with gzip.open(args.output, "wt") as fout:
        fout.write(cset.model_dump_json(exclude_unset=True))
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
