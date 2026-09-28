"""Data/MC normalization factor alpha per AK15 jet pT bin.

    alpha(pT) = N(Data) / N(MC)     [pre-tag: TvsQCD summed over fail+pass]

alpha is the pre-tag data/MC normalization correction of toptag/README.md
("Common: Data/MC Normalization Factor alpha"). It is defined PER CONTROL
REGION -- `gcr` for the QCD-jet calibration, `tmcr`+`tecr` summed for the
Top-jet calibration -- and per pT bin, summed over every background process
and every jet class. It is NOT an inclusive all-event quantity.

Reads the `fjpt_bin` histogram (region x fjpt_bin x TvsQCD), which has data
filled and the same pT edges as `fjpt_cat`; fjpt_cat itself has no data, since
truth matching is undefined for data. (2022 test files: `fjpt`, see ERAS.)
"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
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
ALPHA_VAR = ERAS[YEAR]["alpha_var"]

VAR_LABEL = r"AK15 jet $p_{\mathrm{T}}$ [GeV]"
OUTDIR = os.path.join("plots_effSF", "alpha", YEAR)  # plots_effSF/<kind>/<era>/

# The fjpt axis ends in a very wide 600-2000 GeV bin. It is drawn as an
# overflow bin closing at X_MAX: bin CONTENTS are untouched, only the drawn
# right edge changes.
X_MAX = 1000.0

# Which data dataset is the one recorded for each control region.
REGION_DATASET = {
    "sr": "MET", "wmcr": "MET", "tmcr": "MET", "zmcr": "MET",
    "wecr": "EGamma", "tecr": "EGamma", "gcr": "EGamma", "zecr": "EGamma",
}

# NOTE ON INPUT FILES: each alpha is computed from the SAME file as the MC it
# is meant to correct -- TEST2 for the gcr / QCD-jet calibration, TEST for the
# ttbar CRs / Top-jet calibration, matching plot_toptag_eff_toptagTEST.py. For
# 2022 these are different MC samples; do not compare yields across them.

JOBS = {
    # QCD-jet calibration: alpha from the photon CR alone.
    "alpha_gcr": dict(
        file=TEST2,
        series=[dict(regions=["gcr"], label=r"$\gamma$ CR",
                     color="black", marker="o")],
        title=r"$\gamma$ CR", sublabel="pre-tag, all bkg processes"),
    # Top-jet calibration: ONE alpha from the summed tmcr+tecr yields.
    "alpha_ttcr_combined": dict(
        file=TEST,
        series=[dict(regions=["tmcr", "tecr"], label=r"$t\bar{t}$ CR ($\mu$+e)",
                     color="black", marker="o")],
        title=r"$t\bar{t}$ CR ($\mu$+e)", sublabel="pre-tag, all bkg processes"),
    # Cross-check: is the combined alpha hiding a mu/e difference?
    "alpha_ttcr_split": dict(
        file=TEST,
        series=[dict(regions=["tmcr"], label=r"$t\bar{t}$($\mu$) CR",
                     color="#1f77b4", marker="s"),
                dict(regions=["tecr"], label=r"$t\bar{t}$(e) CR",
                     color="#ff7f0e", marker="^"),
                dict(regions=["tmcr", "tecr"], label=r"combined ($\mu$+e)",
                     color="black", marker="o")],
        title=r"$t\bar{t}$ CR", sublabel="pre-tag, all bkg processes"),
}

_cache = {}


def load(filename):
    """Load and cache one .scaled file."""
    if filename not in _cache:
        with open(filename, "rb") as f:
            _cache[filename] = cloudpickle.loads(lz4.frame.decompress(f.read()))
    return _cache[filename]


def data_side_hists(raw):
    """Return the (bkg, data) dicts of ALPHA_VAR, after checking they are usable.

    Its yields are combined bin-by-bin with the class-split fjpt_cat MC, so it
    must be filled and must have exactly the fjpt_cat pT edges.
    """
    if ALPHA_VAR not in raw["bkg"] or not any(
            len(h.axes["region"]) for h in raw["bkg"][ALPHA_VAR].values()):
        raise SystemExit(f"'{ALPHA_VAR}' is missing or empty in this {YEAR} file"
                         " -- it must be filled by the processor")
    cat_edges = next(iter(raw["bkg"]["fjpt_cat"].values())).axes["fjpt"].edges
    for kind in ("bkg", "data"):
        for name, h in raw[kind][ALPHA_VAR].items():
            edges = h.axes[ALPHA_VAR].edges
            if len(edges) != len(cat_edges) or not np.allclose(edges, cat_edges):
                raise SystemExit(f"{kind}/{ALPHA_VAR}/{name}: pT edges {edges}"
                                 f" differ from fjpt_cat {cat_edges}")
    return raw["bkg"][ALPHA_VAR], raw["data"][ALPHA_VAR]


def pretag_yield(hist_dict, region, procs=None):
    """(values, variances, edges) vs pT, pre-tag, summed over processes.

    Pre-tag means the TvsQCD axis is summed over both the fail and the pass
    bin: alpha is defined before any top-tagging requirement.
    """
    names = list(hist_dict) if procs is None else [p for p in procs if p in hist_dict]
    edges = vals = vars_ = None
    for proc in names:
        h = hist_dict[proc]
        if region not in list(h.axes["region"]):
            continue
        if edges is None:
            edges = h.axes[ALPHA_VAR].edges
            vals = np.zeros(len(edges) - 1)
            vars_ = np.zeros(len(edges) - 1)
        hs = h[{"region": region, "TvsQCD": sum}]
        v, e = hs.values(), hs.variances()
        vals += np.where(np.isfinite(v), v, 0)
        vars_ += np.where(np.isfinite(e), e, 0)
    return vals, vars_, edges


def alpha(filename, regions):
    """Return (alpha, alpha_err, n_data, n_mc, edges) vs fjpt.

    Regions are SUMMED before the ratio is taken -- one alpha for the region
    set, not an average of per-region ratios. Data and MC are independent
    samples, so sigma_alpha/alpha = sqrt(var_D/D^2 + var_M/M^2).
    """
    bkg, data = data_side_hists(load(filename))

    edges = n_data = n_mc = v_data = v_mc = None
    for region in regions:
        d, vd, edges = pretag_yield(data, region, [REGION_DATASET[region]])
        m, vm, _ = pretag_yield(bkg, region)
        n_data = d if n_data is None else n_data + d
        v_data = vd if v_data is None else v_data + vd
        n_mc = m if n_mc is None else n_mc + m
        v_mc = vm if v_mc is None else v_mc + vm

    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.where(n_mc > 0, n_data / n_mc, np.nan)
        err = np.where(
            (n_mc > 0) & (n_data > 0),
            a * np.sqrt(v_data / n_data ** 2 + v_mc / n_mc ** 2),
            np.nan,
        )
    return a, err, n_data, n_mc, edges


def make_alpha_plot(key, zoom=False):
    """alpha vs fjpt for one job, one or more region series overlaid.

    zoom=False draws the alpha=1 reference line, which sets the y-range and
    shows how far the normalization sits from unity. zoom=True drops the line
    and lets the y-axis auto-range on the points, so the pT dependence and the
    error bars are actually readable.
    """
    cfg = JOBS[key]
    results = [(s, alpha(cfg["file"], s["regions"])) for s in cfg["series"]]

    edges = results[0][1][4].copy()
    edges[-1] = X_MAX  # draw the last bin as an overflow bin
    centers = 0.5 * (edges[:-1] + edges[1:])
    xerr = 0.5 * (edges[1:] - edges[:-1])

    fig, ax = plt.subplots(figsize=(10, 9))

    # Offset overlaid series slightly so their error bars stay readable.
    n_series = len(results)
    for i, (s, (a, err, _, _, _)) in enumerate(results):
        dx = 0.0 if n_series == 1 else (i - (n_series - 1) / 2) * 0.10 * 2 * xerr
        ax.errorbar(centers + dx, a, yerr=err, xerr=xerr,
                    fmt=s["marker"], color=s["color"], markersize=7,
                    linewidth=1.5, capsize=0, label=s["label"], zorder=5)

    if not zoom:
        ax.axhline(1.0, color="red", linestyle="--", linewidth=1.5,
                   zorder=1, label=r"$\alpha = 1$")

    ax.set_yscale("log")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel(VAR_LABEL)
    ax.set_ylabel(r"$\alpha = N(\mathrm{Data})\,/\,N(\mathrm{MC})$")

    ax.grid(axis="y", color="gray", linestyle="-", linewidth=0.6,
            alpha=0.4, zorder=0)
    ax.set_axisbelow(True)

    leg = ax.legend(loc="upper right", fontsize=14,
                    title=f"{cfg['title']}\n{cfg['sublabel']}")
    leg.get_title().set_ha("center")
    leg.get_title().set_fontweight("bold")
    leg.get_title().set_fontsize(10)

    hep.cms.label("Private Work", data=True, lumi=LUMI, com=13.6,
                  ax=ax, fontsize=20)

    os.makedirs(OUTDIR, exist_ok=True)
    base = os.path.join(OUTDIR, key + ("_zoom" if zoom else ""))
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    if zoom:  # the table is identical for both variants
        return
    print(f"  [{os.path.basename(cfg['file'])}]")
    for s, (a, err, n_data, n_mc, _) in results:
        print(f"    {'+'.join(s['regions'])}")
        for i in range(len(a)):
            rng = f"{edges[i]:.0f}-{edges[i+1]:.0f}" + (" (of)" if i == len(a) - 1 else "")
            print(f"      {rng:>14s}  N_data={n_data[i]:10.5g}  N_mc={n_mc[i]:10.5g}"
                  f"  alpha={a[i]:12.5g} +/- {err[i]:.5g}")


# ─── Run ─────────────────────────────────────────────────────────────────────
for key in JOBS:
    print(f"{key} ...")
    make_alpha_plot(key, zoom=False)
    make_alpha_plot(key, zoom=True)
print(f"Done. Plots saved in {OUTDIR}/")
