"""AK15 jet pT (fjpt) stacks per aggregated jet class, from fjpt_cat.

One plot per region x top-tag mode x aggregated calibration class
(QCD-jet / B-jet / Top-jet). Backgrounds only -- fjpt_cat carries truth
quark/b-flavor matching (n_qinfj / n_binfj), so data and signal are empty for
this variable by construction and no ratio panel is drawn.
"""

import warnings
warnings.filterwarnings("ignore")

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import mplhep as hep
import lz4.frame
import cloudpickle

# ─── Path resolution ────────────────────────────────────────────────────────
# This script lives in toptag/studies/, but it reads the .scaled files from,
# and writes its plot directory into, the analysis working directory two
# levels up. Everything below therefore uses paths relative to that root, so
# the script behaves identically no matter which directory it is run from.
WORK_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(WORK_DIR)


hep.style.use("CMS")

# ─── Load ───────────────────────────────────────────────────────────────────
INPUT_FILE = "hadmonotop2022_toptagTEST.scaled"
LUMI = 7.99  # 2022preEE, fb^-1

with open(INPUT_FILE, "rb") as f:
    raw = cloudpickle.loads(lz4.frame.decompress(f.read()))

bkg_data = raw["bkg"]

VAR = "fjpt_cat"
VAR_LABEL = r"AK15 jet $p_{\mathrm{T}}$ [GeV]"

# The fjpt axis ends in a very wide 600-2000 GeV bin that dominates the x-range
# while holding few events. It is redrawn as an overflow bin closing at 1000 GeV:
# bin CONTENTS are untouched, only the drawn right edge changes.
X_MAX = 1000.0
OUTDIR_BASE = "plots_2022_toptagTEST_0921_fjptcat_agg"

# ─── Config ──────────────────────────────────────────────────────────────────
REGIONS = ["sr", "wmcr", "tmcr", "zmcr", "wecr", "tecr", "gcr", "zecr"]

REGION_LABEL = {
    "sr":   "SR",
    "wmcr": "W(mu) CR",
    "tmcr": r"$t\bar{t}$(mu) CR",
    "zmcr": "Z(mumu) CR",
    "wecr": "W(e) CR",
    "tecr": r"$t\bar{t}$(e) CR",
    "gcr":  "gamma CR",
    "zecr": "Z(ee) CR",
}

# Per-region stack order (bottom → top); VV = WW+WZ+ZZ combined.
# Only processes actually present in the input file are drawn (see
# get_bkg_arrays), so entries missing from this file are simply skipped.
BKG_ORDER = {
    "sr":   [r"Z ($\nu\nu$) + Jets", r"W ($\ell\nu$) + Jets", "VV", "TT", r"G + Jets", "ST", r"Z ($\ell\ell$) + Jets", "QCD Multijet"][::-1],
    "wmcr": [r"W ($\ell\nu$) + Jets", "TT", "QCD Multijet", "VV", "ST", r"Z ($\ell\ell$) + Jets", r"Z ($\nu\nu$) + Jets", r"G + Jets"][::-1],
    "tmcr": [r"W ($\ell\nu$) + Jets", "TT", "QCD Multijet", "VV", "ST", r"Z ($\ell\ell$) + Jets", r"Z ($\nu\nu$) + Jets", r"G + Jets"][::-1],
    "wecr": [r"W ($\ell\nu$) + Jets", "TT", "QCD Multijet", "VV", "ST", r"Z ($\ell\ell$) + Jets", r"Z ($\nu\nu$) + Jets", r"G + Jets"][::-1],
    "tecr": [r"W ($\ell\nu$) + Jets", "TT", "QCD Multijet", "VV", "ST", r"Z ($\ell\ell$) + Jets", r"Z ($\nu\nu$) + Jets", r"G + Jets"][::-1],
    "zecr": [r"Z ($\ell\ell$) + Jets", r"W ($\ell\nu$) + Jets", "TT", "QCD Multijet", "VV", "ST", r"Z ($\nu\nu$) + Jets", r"G + Jets"][::-1],
    "zmcr": [r"Z ($\ell\ell$) + Jets", r"W ($\ell\nu$) + Jets", "TT", "QCD Multijet", "VV", "ST", r"Z ($\nu\nu$) + Jets", r"G + Jets"][::-1],
    "gcr":  [r"G + Jets", r"W ($\ell\nu$) + Jets", "TT", "QCD Multijet", "VV", "ST", r"Z ($\ell\ell$) + Jets", r"Z ($\nu\nu$) + Jets"][::-1],
}

# WW, WZ, ZZ are merged into VV
VV_PROCS = ["WW", "WZ", "ZZ"]

BKG_COLOR = {
    "QCD Multijet":            "#2ca02c",
    r"Z ($\ell\ell$) + Jets":  "#aec7e8",
    r"Z ($\nu\nu$) + Jets":    "#ffbb78",
    r"G + Jets":               "#98df8a",
    r"W ($\ell\nu$) + Jets":   "#1f77b4",
    "ST":                      "#e377c2",
    "TT":                      "#ff7f0e",
    "VV":                      "#9467bd",
}

BKG_LEGEND = {
    "QCD Multijet":            "QCD Multijet",
    r"Z ($\ell\ell$) + Jets":  r"Z($\ell\ell$)+jets",
    r"Z ($\nu\nu$) + Jets":    r"Z($\nu\nu$)+jets",
    r"G + Jets":               r"$\gamma$+jets",
    r"W ($\ell\nu$) + Jets":   r"W($\ell\nu$)+jets",
    "ST":                      "Single top",
    "TT":                      r"$t\bar{t}$",
    "VV":                      "Diboson",
}

# 8-class AK15 jet classification from toptag/README.md, built from the
# n_qinfj / n_binfj axes on fjpt_cat. "ge1" = n_binfj >= 1, "gt3" = n_qinfj > 3,
# "any" = summed over that axis. The 8 classes are a strict partition of the
# (n_qinfj, n_binfj) plane, so the aggregated groups below are plain sums.
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

# Aggregated calibration classes (README: QCD-jet=1,3; B-jet=2,4; Top-jet=5,6,7,8)
AGG_GROUPS = {
    "QCD-jet": [1, 3],
    "B-jet":   [2, 4],
    "Top-jet": [5, 6, 7, 8],
}

AGG_LABEL = {
    "QCD-jet": "QCD-jet (class 1,3)",
    "B-jet":   "B-jet (class 2,4)",
    "Top-jet": "Top-jet (class 5,6,7,8)",
}

# TvsQCD axis edges are [0, 0.33, 1] -> bin 0 = fail, bin 1 = pass (WP = 0.33)
# None means "inclusive": sum over the TvsQCD axis.
TAG_MODES = {
    "inclusive":   None,
    "toptag_fail": 0,
    "toptag_pass": 1,
}

TAG_TITLE = {
    "inclusive":   "",
    "toptag_fail": ", top-tag fail",
    "toptag_pass": ", top-tag pass",
}


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
    """True if this histogram's region axis actually contains `region`.

    Not every process is filled in every region in these test files.
    """
    return region in list(h.axes["region"])


def group_spectrum(h, region, cat_idx, group):
    """Return (values, variances, edges) over fjpt for one aggregated class.

    cat_idx selects the TvsQCD bin, or None to sum over the whole axis. The
    classes in a group are disjoint, so values and variances simply add.
    """
    edges = h.axes["fjpt"].edges
    vals = np.zeros(len(edges) - 1)
    vars_ = np.zeros(len(edges) - 1)
    if not has_region(h, region):
        return vals, vars_, edges

    for class_id in AGG_GROUPS[group]:
        sel = {"region": region, "TvsQCD": sum if cat_idx is None else cat_idx}
        sel.update(class_selector(class_id))
        hs = h[sel]
        v, e = hs.values(), hs.variances()
        vals += np.where(np.isfinite(v), v, 0)
        vars_ += np.where(np.isfinite(e), e, 0)
    return vals, vars_, edges


def get_bkg_arrays(region, cat_idx, group):
    """Return list of (label, values, variances, edges), bottom → top.

    VV is WW+WZ+ZZ summed. Only processes actually present in bkg_data[VAR]
    are included, so this works unchanged as more backgrounds are added.
    """
    out = []
    for label in BKG_ORDER[region]:
        if label == "VV":
            present = [p for p in VV_PROCS
                       if p in bkg_data[VAR] and has_region(bkg_data[VAR][p], region)]
            if not present:
                continue
            vals = vars_ = edges = None
            for p in present:
                v, e, ed = group_spectrum(bkg_data[VAR][p], region, cat_idx, group)
                vals = v if vals is None else vals + v
                vars_ = e if vars_ is None else vars_ + e
                edges = ed
        else:
            if label not in bkg_data[VAR]:
                continue
            vals, vars_, edges = group_spectrum(
                bkg_data[VAR][label], region, cat_idx, group)
        if vals.sum() <= 0:
            continue
        out.append((label, vals, vars_, edges))
    return out


def make_stack_plot(region, tag_mode, group, scale, outdir, outname):
    """Background-only fjpt stack for one region / top-tag mode / jet class.

    `scale` is "lin" or "log" and sets the y-axis scaling.
    """
    bkg_list = get_bkg_arrays(region, TAG_MODES[tag_mode], group)
    if not bkg_list:
        print(f"  Skipping {region}/{tag_mode}/{group}: no bkg entries")
        return

    edges = bkg_list[0][3].copy()
    edges[-1] = X_MAX  # draw the last bin as an overflow bin ending at X_MAX
    cumsum = np.cumsum([v for _, v, _, _ in bkg_list], axis=0)
    total_bkg = cumsum[-1]
    tot_err = np.sqrt(np.sum([e for _, _, e, _ in bkg_list], axis=0))

    fig, ax = plt.subplots(figsize=(10, 9))

    prev = np.zeros(len(total_bkg))
    for (label, _, _, _), top in zip(bkg_list, cumsum):
        ax.stairs(top, edges, baseline=prev, fill=True,
                  color=BKG_COLOR[label], edgecolor="black", linewidth=1.5,
                  label=BKG_LEGEND[label])
        prev = top.copy()

    # Total background uncertainty band from the summed per-process variances
    ax.stairs(total_bkg + tot_err, edges, baseline=total_bkg - tot_err,
              fill=True, color="none", edgecolor="gray", hatch="///",
              linewidth=0, label="Bkg. unc.")

    top_env = total_bkg + tot_err
    if scale == "log":
        ax.set_yscale("log")
        # Decade-rounded ylim + decade-only major ticks: avoids matplotlib's
        # "6 x 10^0"-style multiplier labels on narrow-range log axes, which
        # widen the left margin and push the CMS header into the lumi label.
        # Floor from the smallest nonzero STACK BOUNDARY, not from the total:
        # flooring on the total pushes the smaller processes below the axis and
        # they vanish from the plot.
        nonzero = cumsum[cumsum > 0]
        ylo = 10 ** np.floor(np.log10(nonzero.min())) if len(nonzero) else 1e-2
        env_nz = top_env[top_env > 0]
        yhi = 10 ** np.ceil(np.log10(env_nz.max()) + 1) if len(env_nz) else 1e2
        ax.set_ylim(ylo, yhi)
        ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0))
        ax.yaxis.set_minor_locator(
            mticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    else:
        ymax = top_env.max()
        ax.set_ylim(0, ymax * 1.6 if ymax > 0 else 1)

    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel(VAR_LABEL)
    ax.set_ylabel("Events / bin")

    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles[::-1], labels[::-1], ncol=2, fontsize=14,
                    loc="upper right",
                    title=f"{REGION_LABEL[region]}{TAG_TITLE[tag_mode]}\n{AGG_LABEL[group]}")
    leg.get_title().set_ha("center")
    leg.get_title().set_fontweight("bold")
    leg.get_title().set_fontsize(10)

    hep.cms.label("Private Work", data=False, lumi=LUMI, com=13.6,
                  ax=ax, fontsize=20)

    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, outname + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(outdir, outname + ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


# ─── Run ─────────────────────────────────────────────────────────────────────
# One subdirectory per top-tag mode and aggregated class, one fjpt stack per
# region inside it.
DIRNAME = {"QCD-jet": "qcdjet", "B-jet": "bjet", "Top-jet": "topjet"}
SUFFIX = {"lin": "", "log": "_log"}

for tag_mode in TAG_MODES:
    for group in AGG_GROUPS:
        outdir = os.path.join(OUTDIR_BASE, tag_mode, DIRNAME[group])
        for scale in SUFFIX:
            for region in REGIONS:
                print(f"Plotting fjpt in {region} [{tag_mode}, {group}, {scale}] ...")
                make_stack_plot(region, tag_mode, group, scale, outdir,
                                f"stack_{region}_fjpt{SUFFIX[scale]}")
        print(f"Done: {outdir}/")
