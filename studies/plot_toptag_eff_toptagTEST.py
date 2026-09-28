"""Top-tag efficiency vs AK15 jet pT, per aggregated jet class.

    eff(pT) = N(top-tagged, class) / N(top-tagged + untagged, class)

For the Top-jet class this is the tagging efficiency; for the QCD-jet class it
is the mis-tag rate. Each job declares its own input file, jet class, control
region(s) and process selection -- see JOBS below.
"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import mplhep as hep
import lz4.frame
import cloudpickle

# ─── Path resolution ────────────────────────────────────────────────────────
# This script lives in toptag/studies/, but it reads the .scaled files from,
# and writes its plot directory into, the toptag/ directory one level up.
# Everything below therefore uses paths relative to that root, so the script
# behaves identically no matter which directory it is run from.
WORK_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(WORK_DIR)


hep.style.use("CMS")

# Data-taking era: optional first command-line argument, default 2022preEE.
# Selects the input files, the luminosity label and the output sub-directory.
YEAR = sys.argv[1] if len(sys.argv) > 1 else "2022preEE"

# fb^-1, same values as toptag/histmaker.py LUMIS
LUMIS = {"2022preEE": 7.99, "2022postEE": 26.68,
         "2023preBPix": 17.96, "2023postBPix": 9.68}
LUMI = LUMIS[YEAR]

# Input files per era. "top" feeds the ttbar-CR (Top-jet) jobs, "qcd" the
# gamma-CR (QCD-jet) jobs -- the 2022 "top" file has no G + Jets sample, so
# 2022 needs a second file. "alpha_var" is the histogram the data-side yields
# (alpha, N(Data)) are read from: `fjpt_bin`, binned exactly like fjpt_cat.
# The 2022 test files predate fjpt_bin and use `fjpt`, which has the same
# edges there. Efficiencies always come from fjpt_cat.
ERAS = {
    "2022preEE":   dict(top="hadmonotop2022_toptagTEST.scaled",
                        qcd="hadmonotop2022_toptagTEST2.scaled", alpha_var="fjpt"),
    "2023preBPix": dict(top="hadmonotop2023_toptagTEST.scaled",
                        qcd="hadmonotop2023_toptagTEST.scaled", alpha_var="fjpt_bin"),
}
TEST, TEST2 = ERAS[YEAR]["top"], ERAS[YEAR]["qcd"]

VAR = "fjpt_cat"
VAR_LABEL = r"AK15 jet $p_{\mathrm{T}}$ [GeV]"
OUTDIR = os.path.join("plots_effSF", "eff", YEAR)  # plots_effSF/<kind>/<era>/

# The fjpt axis ends in a very wide 600-2000 GeV bin. It is drawn as an
# overflow bin closing at X_MAX: bin CONTENTS are untouched, only the drawn
# right edge changes.
X_MAX = 1000.0

# 8-class AK15 jet classification from toptag/README.md. The 8 classes are a
# strict partition of the (n_qinfj, n_binfj) plane, so an aggregated class is a
# plain sum over its constituent classes.
CLASS_DEF = {
    1: dict(nq=0,     nb=0),
    2: dict(nq=0,     nb="ge1"),
    3: dict(nq=1,     nb=0),
    4: dict(nq=1,     nb="ge1"),
    5: dict(nq=2,     nb=0),
    6: dict(nq=2,     nb="ge1"),
    7: dict(nq=3,     nb="any"),
    8: dict(nq="gt3", nb="any"),
}

# Aggregated calibration classes (README)
AGG_GROUPS = {
    "QCD-jet": [1, 3],
    "B-jet":   [2, 4],
    "Top-jet": [5, 6, 7, 8],
}

# TvsQCD axis edges are [0, 0.33, 1] -> bin 0 = fail, bin 1 = pass (WP = 0.33)
TAG_FAIL, TAG_PASS = 0, 1

# ─── Jobs ────────────────────────────────────────────────────────────────────
# procs=None means "sum over every background process in the file".
#
# NOTE ON INPUT FILES: Top-jet jobs read TEST, QCD-jet jobs read TEST2 (see
# ERAS). For 2022 these are different MC test samples -- do not compare their
# absolute yields.

JOBS = {
    "eff_toptag_tmcr": dict(
        file=TEST, group="Top-jet", regions=["tmcr"], procs=None,
        label=r"$t\bar{t}$($\mu$) CR", sublabel="Top-jet (class 5,6,7,8), all bkg",
        ylabel="Top-tag efficiency", ylim=(0, 1.35), ygrid=0.2),
    "eff_toptag_tecr": dict(
        file=TEST, group="Top-jet", regions=["tecr"], procs=None,
        label=r"$t\bar{t}$(e) CR", sublabel="Top-jet (class 5,6,7,8), all bkg",
        ylabel="Top-tag efficiency", ylim=(0, 1.35), ygrid=0.2),
    "eff_toptag_combined": dict(
        file=TEST, group="Top-jet", regions=["tmcr", "tecr"], procs=None,
        label=r"$t\bar{t}$ CR ($\mu$+e)", sublabel="Top-jet (class 5,6,7,8), all bkg",
        ylabel="Top-tag efficiency", ylim=(0, 1.35), ygrid=0.2),
    # QCD-jet mis-tag rate from the G+Jets sample in the photon CR.
    "eff_qcdjet_gcr": dict(
        file=TEST2, group="QCD-jet", regions=["gcr"], procs=[r"G + Jets"],
        label=r"$\gamma$ CR", sublabel=r"QCD-jet (class 1,3), $\gamma$+jets",
        ylabel="Top-tag mis-tag rate", ylim=(0, 1.35), ygrid=0.2),
    # Same numbers, zoomed: the mis-tag rate is ~0.03 and is unreadable on the
    # 0-1.35 axis used for direct comparison with the Top-jet plots.
    "eff_qcdjet_gcr_zoom": dict(
        file=TEST2, group="QCD-jet", regions=["gcr"], procs=[r"G + Jets"],
        label=r"$\gamma$ CR", sublabel=r"QCD-jet (class 1,3), $\gamma$+jets",
        ylabel="Top-tag mis-tag rate", ylim=(0, 0.06), ygrid=0.01),
}

_cache = {}


def load_bkg(filename):
    """Load and cache the bkg[fjpt_cat] dict from one .scaled file."""
    if filename not in _cache:
        with open(filename, "rb") as f:
            _cache[filename] = cloudpickle.loads(lz4.frame.decompress(f.read()))
    return _cache[filename]["bkg"][VAR]


def class_selector(class_id):
    """Build the fjpt_cat UHI selection dict for one of the 8 jet classes."""
    nq, nb = CLASS_DEF[class_id]["nq"], CLASS_DEF[class_id]["nb"]
    sel = {"n_qinfj": slice(4, None, sum) if nq == "gt3" else nq}
    if nb == "any":
        sel["n_binfj"] = sum
    elif nb == "ge1":
        sel["n_binfj"] = slice(1, None, sum)
    else:
        sel["n_binfj"] = nb
    return sel


def has_region(h, region):
    """True if this histogram's region axis actually contains `region`."""
    return region in list(h.axes["region"])


def class_spectrum(bkg, group, regions, procs, cat_idx):
    """Return (values, variances, edges) over fjpt for one aggregated class.

    Summed over the listed regions and over `procs` (all processes if None),
    restricted to the given TvsQCD bin.
    """
    names = list(bkg) if procs is None else [p for p in procs if p in bkg]
    edges = vals = vars_ = None
    for proc in names:
        h = bkg[proc]
        if edges is None:
            edges = h.axes["fjpt"].edges
            vals = np.zeros(len(edges) - 1)
            vars_ = np.zeros(len(edges) - 1)
        for region in regions:
            if not has_region(h, region):
                continue
            for class_id in AGG_GROUPS[group]:
                sel = {"region": region, "TvsQCD": cat_idx}
                sel.update(class_selector(class_id))
                hs = h[sel]
                v, e = hs.values(), hs.variances()
                vals += np.where(np.isfinite(v), v, 0)
                vars_ += np.where(np.isfinite(e), e, 0)
    return vals, vars_, edges


def efficiency(cfg):
    """Return (eff, eff_err, n_pass, n_tot, edges) vs fjpt for one job.

    Pass and fail are disjoint samples, so the error on eff = P/(P+F) is
    propagated as sigma^2 = [(1-eff)^2 var_P + eff^2 var_F] / N^2, the
    weighted-histogram generalisation of the binomial error.
    """
    bkg = load_bkg(cfg["file"])
    args = (bkg, cfg["group"], cfg["regions"], cfg["procs"])
    n_pass, var_pass, edges = class_spectrum(*args, TAG_PASS)
    n_fail, var_fail, _ = class_spectrum(*args, TAG_FAIL)

    n_tot = n_pass + n_fail
    with np.errstate(divide="ignore", invalid="ignore"):
        eff = np.where(n_tot > 0, n_pass / n_tot, np.nan)
        err = np.where(
            n_tot > 0,
            np.sqrt(((1 - eff) ** 2 * var_pass + eff ** 2 * var_fail)) / n_tot,
            np.nan,
        )
    return eff, err, n_pass, n_tot, edges


def make_eff_plot(key):
    """Efficiency vs fjpt for one job."""
    cfg = JOBS[key]
    eff, err, n_pass, n_tot, edges = efficiency(cfg)

    if not np.any(n_tot > 0):
        print(f"  Skipping {key}: no events in this class/region")
        return

    edges = edges.copy()
    edges[-1] = X_MAX  # draw the last bin as an overflow bin
    centers = 0.5 * (edges[:-1] + edges[1:])
    xerr = 0.5 * (edges[1:] - edges[:-1])

    fig, ax = plt.subplots(figsize=(10, 9))

    ax.errorbar(centers, eff, yerr=err, xerr=xerr,
                fmt="o", color="black", markersize=7,
                linewidth=1.5, capsize=0, label="MC efficiency", zorder=5)

    ax.set_ylim(*cfg["ylim"])
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel(VAR_LABEL)
    ax.set_ylabel(cfg["ylabel"])

    # Horizontal guide lines, behind the points
    ax.yaxis.set_major_locator(mticker.MultipleLocator(cfg["ygrid"]))
    ax.grid(axis="y", color="gray", linestyle="-", linewidth=0.6,
            alpha=0.4, zorder=0)
    ax.set_axisbelow(True)

    leg = ax.legend(loc="upper right", fontsize=14,
                    title=f"{cfg['label']}\n{cfg['sublabel']}")
    leg.get_title().set_ha("center")
    leg.get_title().set_fontweight("bold")
    leg.get_title().set_fontsize(10)

    hep.cms.label("Private Work", data=False, lumi=LUMI, com=13.6,
                  ax=ax, fontsize=20)

    os.makedirs(OUTDIR, exist_ok=True)
    base = os.path.join(OUTDIR, key)
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"  {cfg['label']} / {cfg['group']}  [{os.path.basename(cfg['file'])}]")
    for i in range(len(eff)):
        rng = f"{edges[i]:.0f}-{edges[i+1]:.0f}" + (" (of)" if i == len(eff) - 1 else "")
        if n_tot[i] > 0:
            print(f"    {rng:>14s}  pass={n_pass[i]:10.5g}  tot={n_tot[i]:10.5g}"
                  f"  eff={eff[i]:.4f} +/- {err[i]:.4f}")
        else:
            print(f"    {rng:>14s}  empty")


# ─── Run ─────────────────────────────────────────────────────────────────────
for key in JOBS:
    print(f"{key} ...")
    make_eff_plot(key)
print(f"Done. Plots saved in {OUTDIR}/")
