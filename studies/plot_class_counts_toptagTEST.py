"""Stacked background event counts per AK15 jet class (fjpt_cat).

One plot per region: x-axis = the 8 truth-level jet classes (class1..class8)
or the 3 aggregated classes, y-axis = number of events, stacked by background
process. Backgrounds only -- fjpt_cat carries truth quark/b-flavor matching
(n_qinfj / n_binfj), so data and signal are empty for this variable by
construction and no ratio panel is drawn.

Regions: every CR plus the combined tt(e+mu) CR "ttcr" (tmcr + tecr summed).

    python studies/plot_class_counts_toptagTEST.py -f hadmonotop2023_toptagTEST3

Outputs: plots_classcount_<version>/<era>/<tag mode>/{classcount,aggcount}_<region>.{pdf,png}
"""

import warnings
warnings.filterwarnings("ignore")

import os
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import lz4.frame
import cloudpickle

from effsf_cli import base_parser, setup

hep.style.use("CMS")

# ─── Load ───────────────────────────────────────────────────────────────────
_args = base_parser(__doc__.splitlines()[0], alpha_var=False).parse_args()
_cfg = setup(_args)  # also changes into toptag/
LUMI = _cfg.lumi

with open(_cfg.file, "rb") as f:
    raw = cloudpickle.loads(lz4.frame.decompress(f.read()))

bkg_data = raw["bkg"]

VAR = "fjpt_cat"
OUTDIR_BASE = os.path.join(_cfg.outbase.replace("plots_effSF_", "plots_classcount_"), _cfg.era)

# ─── Config ──────────────────────────────────────────────────────────────────
REGIONS = ["wmcr", "tmcr", "zmcr", "wecr", "tecr", "gcr", "zecr", "ttcr"]

# Region -> regions in the file that are summed for it
REGION_SOURCES = {r: [r] for r in REGIONS}
REGION_SOURCES["ttcr"] = ["tmcr", "tecr"]

REGION_LABEL = {
    "sr":   "SR",
    "wmcr": "W(mu) CR",
    "tmcr": r"$t\bar{t}$(mu) CR",
    "zmcr": "Z(mumu) CR",
    "wecr": "W(e) CR",
    "tecr": r"$t\bar{t}$(e) CR",
    "gcr":  "gamma CR",
    "zecr": "Z(ee) CR",
    "ttcr": r"$t\bar{t}$(e+mu) CR",
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
    "ttcr": [r"W ($\ell\nu$) + Jets", "TT", "QCD Multijet", "VV", "ST", r"Z ($\ell\ell$) + Jets", r"Z ($\nu\nu$) + Jets", r"G + Jets"][::-1],
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
# "any" = summed over that axis.
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

CLASS_LABEL = {
    1: "no q, no b",
    2: "no q, b",
    3: "1q, no b",
    4: "1q, b",
    5: "2q, no b",
    6: "2q, b",
    7: "3q",
    8: ">3q",
}

# Aggregated calibration group per class (README: QCD-jet=1,3; B-jet=2,4; Top-jet=5,6,7,8)
CLASS_GROUP = {1: "QCD-jet", 2: "B-jet", 3: "QCD-jet", 4: "B-jet",
               5: "Top-jet", 6: "Top-jet", 7: "Top-jet", 8: "Top-jet"}

CLASS_IDS = list(CLASS_DEF)

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


def axis_mask(n_bins_flow, spec):
    """Mask over an Integer(0, 10) axis incl. flow bins (index = value + 1)."""
    m = np.zeros(n_bins_flow, dtype=bool)
    if spec == "any":
        m[:] = True
    elif spec == "ge1":
        m[2:] = True           # values >= 1, overflow included
    elif spec == "gt3":
        m[5:] = True           # values >= 4, overflow included
    else:
        m[spec + 1] = True
    return m


def class_mask(class_id, shape):
    """2D (n_qinfj, n_binfj) mask selecting one jet class, flow bins included."""
    d = CLASS_DEF[class_id]
    return np.outer(axis_mask(shape[0], d["nq"]), axis_mask(shape[1], d["nb"]))


def clean(a):
    """Replace non-finite entries (bad MC weights) by 0."""
    return np.where(np.isfinite(a), a, 0)


def has_region(h, region):
    """True if this histogram's region axis actually contains `region`.

    Not every process is filled in every region in these test files.
    """
    return region in list(h.axes["region"])


def class_yields_one(h, region, cat_idx):
    """(values, variances) over the 8 classes for one file region, fjpt summed.

    cat_idx selects the TvsQCD bin, or None to sum over the whole axis.
    """
    if not has_region(h, region):
        return np.zeros(len(CLASS_IDS)), np.zeros(len(CLASS_IDS))
    hr = h[{"region": region}]
    # axes: fjpt, TvsQCD (no flow), n_qinfj, n_binfj; fjpt flow excluded
    v = clean(hr.values(flow=True))[1:-1].sum(axis=0)
    w = clean(hr.variances(flow=True))[1:-1].sum(axis=0)
    v, w = (v.sum(axis=0), w.sum(axis=0)) if cat_idx is None else (v[cat_idx], w[cat_idx])
    masks = [class_mask(cid, v.shape) for cid in CLASS_IDS]
    return np.array([v[m].sum() for m in masks]), np.array([w[m].sum() for m in masks])


def class_yields(h, region, cat_idx):
    """(values, variances) over the 8 classes, summed over REGION_SOURCES[region]."""
    vals, vars_ = np.zeros(len(CLASS_IDS)), np.zeros(len(CLASS_IDS))
    for src in REGION_SOURCES[region]:
        v, w = class_yields_one(h, src, cat_idx)
        vals, vars_ = vals + v, vars_ + w
    return vals, vars_


def get_bkg_arrays(region, cat_idx):
    """Return list of (label, values, variances) per stack entry, bottom → top.

    VV is WW+WZ+ZZ summed. Only processes actually present in bkg_data[VAR]
    are included, so this works unchanged as more backgrounds are added.
    """
    out = []
    for label in BKG_ORDER[region]:
        if label == "VV":
            present = [p for p in VV_PROCS
                       if p in bkg_data[VAR]]
            if not present:
                continue
            vals = np.zeros(len(CLASS_IDS))
            vars_ = np.zeros(len(CLASS_IDS))
            for p in present:
                v, e = class_yields(bkg_data[VAR][p], region, cat_idx)
                vals += v
                vars_ += e
        else:
            if label not in bkg_data[VAR]:
                continue
            vals, vars_ = class_yields(bkg_data[VAR][label], region, cat_idx)
        if vals.sum() <= 0:
            continue
        out.append((label, vals, vars_))
    return out


# Aggregated calibration classes (README): the 8 classes are a strict partition,
# so group yields and variances are just sums over disjoint class entries.
AGG_GROUPS = {
    "QCD-jet": [1, 3],
    "B-jet":   [2, 4],
    "Top-jet": [5, 6, 7, 8],
}
AGG_ORDER = list(AGG_GROUPS)


def aggregate(bkg_list):
    """Collapse per-class (values, variances) arrays into the 3 calibration groups."""
    idx = {g: [CLASS_IDS.index(c) for c in cs] for g, cs in AGG_GROUPS.items()}
    out = []
    for label, vals, vars_ in bkg_list:
        out.append((
            label,
            np.array([vals[idx[g]].sum() for g in AGG_ORDER]),
            np.array([vars_[idx[g]].sum() for g in AGG_ORDER]),
        ))
    return out


# Binning schemes: x-axis tick labels and the calibration group of each bin.
BINNINGS = {
    "class": dict(
        ticks=[f"class{c}\n{CLASS_LABEL[c]}" for c in CLASS_IDS],
        groups=[CLASS_GROUP[c] for c in CLASS_IDS],
        xlabel="AK15 jet class",
        fontsize=13,
    ),
    "group": dict(
        ticks=[f"{g}\nclass " + ", ".join(str(c) for c in AGG_GROUPS[g])
               for g in AGG_ORDER],
        groups=AGG_ORDER,
        xlabel="AK15 jet class (aggregated)",
        fontsize=16,
    ),
}


def make_plot(region, tag_mode, binning, outdir, outname):
    """Stacked background event counts per jet category for one region/tag mode.

    `binning` is "class" (the 8 jet classes) or "group" (the 3 aggregated
    calibration classes QCD-jet / B-jet / Top-jet).
    """
    bkg_list = get_bkg_arrays(region, TAG_MODES[tag_mode])
    if not bkg_list:
        print(f"  Skipping {region}/{tag_mode}/{binning}: no bkg entries")
        return
    if binning == "group":
        bkg_list = aggregate(bkg_list)

    cfg = BINNINGS[binning]
    x = np.arange(len(cfg["ticks"]))
    fig, ax = plt.subplots(figsize=(12, 9))

    bottom = np.zeros(len(x))
    for label, vals, _ in bkg_list:
        ax.bar(x, vals, bottom=bottom, width=1.0,
               color=BKG_COLOR[label], edgecolor="black", linewidth=1.5,
               label=BKG_LEGEND[label], zorder=2)
        bottom += vals

    # Total background uncertainty band from the summed per-process variances
    tot_var = np.sum([v for _, _, v in bkg_list], axis=0)
    tot_err = np.sqrt(tot_var)
    ax.bar(x, 2 * tot_err, bottom=bottom - tot_err, width=1.0,
           color="none", edgecolor="gray", hatch="///", linewidth=0.0,
           label="Bkg. unc.", zorder=3)

    ymax = bottom.max() + tot_err.max()
    ax.set_ylim(0, ymax * 1.6 if ymax > 0 else 1)
    ax.set_xlim(-0.5, len(x) - 0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(cfg["ticks"], fontsize=cfg["fontsize"])
    ax.set_xlabel(cfg["xlabel"])
    ax.set_ylabel("Events")

    # Shade the QCD-jet / B-jet / Top-jet calibration groups behind the bars
    for i, grp in enumerate(cfg["groups"]):
        if grp == "B-jet":
            ax.axvspan(i - 0.5, i + 0.5, color="black", alpha=0.04, zorder=0)
        elif grp == "Top-jet":
            ax.axvspan(i - 0.5, i + 0.5, color="black", alpha=0.08, zorder=0)

    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles[::-1], labels[::-1], ncol=2, fontsize=14,
                    loc="upper right",
                    title=f"{REGION_LABEL[region]}{TAG_TITLE[tag_mode]}")
    leg.get_title().set_ha("center")
    leg.get_title().set_fontweight("bold")
    leg.get_title().set_fontsize(10)

    hep.cms.label("Private Work", data=False, lumi=LUMI, com=13.6,
                  ax=ax, fontsize=20)

    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, outname + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(outdir, outname + ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


# ─── Run ─────────────────────────────────────────────────────
# One subdirectory per top-tag mode; per region, one 8-class plot and one
# aggregated 3-class plot.
PREFIX = {"class": "classcount", "group": "aggcount"}

for tag_mode in TAG_MODES:
    outdir = os.path.join(OUTDIR_BASE, tag_mode)
    for binning in BINNINGS:
        for region in REGIONS:
            print(f"Plotting {binning} counts in {region} [{tag_mode}] ...")
            make_plot(region, tag_mode, binning, outdir,
                      f"{PREFIX[binning]}_{region}")
    print(f"Done with {tag_mode}. Plots saved in {outdir}/")
