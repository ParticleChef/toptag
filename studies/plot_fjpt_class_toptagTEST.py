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
dat_data = raw["data"]

# fjpt_cat has the same region/fjpt/TvsQCD axes as fjpt, plus n_qinfj and
# n_binfj — the classification is only available there, not on plain fjpt.
VAR = "fjpt"
VAR_CAT = "fjpt_cat"
VAR_LABEL = r"AK15 jet $p_{\mathrm{T}}$ [GeV]"

# ─── Output directory ────────────────────────────────────────────────────────
OUTDIR_BASE = "plots_2022_toptagTEST_0917"

# ─── Config ──────────────────────────────────────────────────────────────────
REGION_DATASET = {
    "sr":   "MET",
    "wmcr": "MET",
    "tmcr": "MET",
    "zmcr": "MET",
    "wecr": "EGamma",
    "tecr": "EGamma",
    "gcr":  "EGamma",
    "zecr": "EGamma",
}

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
# get_bkg_arrays), so this list is safe to keep even when most of these
# aren't in this small TT-only test file yet.
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


def slice_region_class(h, region, class_id):
    """Select region, sum over TvsQCD, and select the given jet class."""
    sel = {"region": region, "TvsQCD": sum}
    sel.update(class_selector(class_id))
    return h[sel]


def get_bkg_arrays(region, class_id):
    """Return list of (label, values, edges) per stack entry; VV is WW+WZ+ZZ summed.

    Only processes actually present in bkg_data[VAR_CAT] are included, so this
    works unchanged as more backgrounds are added to the input file.
    """
    result = []
    for proc in BKG_ORDER[region]:
        if proc == "VV":
            vv_vals, vv_edges = None, None
            for vv in VV_PROCS:
                if vv not in bkg_data.get(VAR_CAT, {}):
                    continue
                h = bkg_data[VAR_CAT][vv]
                if region not in list(h.axes["region"]):
                    continue
                hs = slice_region_class(h, region, class_id)
                v = np.where(np.isfinite(hs.values()), hs.values(), 0)
                e = np.array(list(hs.axes[0].edges))
                vv_vals = v if vv_vals is None else vv_vals + v
                vv_edges = e
            if vv_vals is not None:
                result.append(("VV", vv_vals, vv_edges))
        else:
            if proc not in bkg_data.get(VAR_CAT, {}):
                continue
            h = bkg_data[VAR_CAT][proc]
            if region not in list(h.axes["region"]):
                continue
            hs = slice_region_class(h, region, class_id)
            vals = np.where(np.isfinite(hs.values()), hs.values(), 0)
            edges = np.array(list(hs.axes[0].edges))
            result.append((proc, vals, edges))
    return result


def get_data_arrays(region, class_id):
    """Return (values, bin centers, errors) for the data histogram, or (None, None, None)."""
    dataset = REGION_DATASET[region]
    if dataset not in dat_data.get(VAR_CAT, {}):
        return None, None, None
    h = dat_data[VAR_CAT][dataset]
    if region not in list(h.axes["region"]):
        return None, None, None
    hs = slice_region_class(h, region, class_id)
    vals = hs.values()
    edges = np.array(list(hs.axes[0].edges))
    centers = 0.5 * (edges[:-1] + edges[1:])
    return vals, centers, np.sqrt(vals)


def make_stack_plot(region, class_id, outdir, outname):
    """Data/MC stack + ratio plot of fjpt for one region and one jet class."""
    bkg_list = get_bkg_arrays(region, class_id)
    dat_vals, dat_centers, dat_errs = get_data_arrays(region, class_id)

    if not bkg_list:
        print(f"  Skipping class {class_id}/{region}: no bkg data")
        return

    total_sum = sum(np.sum(v) for _, v, _ in bkg_list)
    data_sum = np.sum(dat_vals) if dat_vals is not None else 0
    if total_sum == 0 and data_sum == 0:
        print(f"  Skipping class {class_id}/{region}: empty")
        return

    edges = bkg_list[0][2]
    bkg_vals_list = [v for _, v, _ in bkg_list]
    bkg_cumsum = np.cumsum(bkg_vals_list, axis=0)
    total_bkg = bkg_cumsum[-1]

    # ─── Figure layout ───────────────────────────────────────────────────────
    fig, (ax_main, ax_ratio) = plt.subplots(
        2, 1, figsize=(8, 8),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        sharex=True,
    )

    # ─── Main panel: stepfilled stack ────────────────────────────────────────
    prev = np.zeros(len(total_bkg))
    for (proc, vals, _), top in zip(bkg_list, bkg_cumsum):
        ax_main.stairs(
            top, edges,
            baseline=prev,
            fill=True,
            color=BKG_COLOR.get(proc, "gray"),
            edgecolor="black",
            linewidth=1.5,
            label=BKG_LEGEND.get(proc, proc),
        )
        prev = top.copy()

    # total bkg uncertainty band
    bkg_err = np.sqrt(total_bkg)
    ax_main.stairs(
        total_bkg + bkg_err, edges,
        baseline=total_bkg - bkg_err,
        fill=True,
        color="none", edgecolor="gray", hatch="///",
        label="Bkg unc.", linewidth=0,
    )

    # data points
    if dat_vals is not None:
        ax_main.errorbar(
            dat_centers, dat_vals, yerr=dat_errs,
            fmt="o", color="black", markersize=4,
            linewidth=1, capsize=2, label="Data",
            zorder=5,
        )

    ax_main.set_yscale("log")

    # Decade-rounded ylim + decade-only major ticks: avoids matplotlib's
    # "6 x 10^0"-style multiplier tick labels on narrow-range log axes, which
    # widen the left margin and push the CMS header text into the lumi label.
    combined = total_bkg + bkg_err
    if dat_vals is not None:
        combined = np.concatenate([combined, dat_vals])
    nonzero = combined[combined > 0]
    ylo = 10 ** np.floor(np.log10(nonzero.min())) if len(nonzero) else 1e-2
    yhi = 10 ** np.ceil(np.log10(nonzero.max()) + 1) if len(nonzero) else 1e2
    ax_main.set_ylim(ylo, yhi)
    ax_main.yaxis.set_major_locator(mticker.LogLocator(base=10.0))
    ax_main.yaxis.set_minor_locator(mticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax_main.yaxis.set_minor_formatter(mticker.NullFormatter())

    ax_main.set_ylabel("Events / bin")
    ax_main.set_xlim(edges[0], edges[-1])

    # Legend: Data first, bkg reversed (top of stack → bottom), Bkg unc. last
    handles, labels = ax_main.get_legend_handles_labels()
    order = []
    for i, lbl in enumerate(labels):
        if lbl == "Data":
            order.insert(0, i)
        elif lbl == "Bkg unc.":
            order.append(i)
        else:
            order.insert(1 if order and labels[order[0]] == "Data" else 0, i)

    leg = ax_main.legend(
        [handles[i] for i in order],
        [labels[i] for i in order],
        title=f"{REGION_LABEL[region]}, class {class_id}: {CLASS_LABEL[class_id]} ({CLASS_GROUP[class_id]})",
        title_fontsize=9,
        ncol=2, fontsize=9, framealpha=0.8,
        loc="upper right",
    )
    leg.get_title().set_ha("center")
    leg.get_title().set_fontweight("bold")

    hep.cms.label(
        "Private Work", data=True,
        lumi=LUMI, com=13.6,
        ax=ax_main,
    )

    # ─── Ratio panel ─────────────────────────────────────────────────────────
    if dat_vals is not None:
        ratio     = np.where(total_bkg > 0, dat_vals / total_bkg, np.nan)
        ratio_err = np.where(total_bkg > 0, dat_errs / total_bkg, np.nan)
        bkg_rel   = np.where(total_bkg > 0, bkg_err / total_bkg, 0)
        bin_w     = edges[1:] - edges[:-1]

        ax_ratio.bar(
            edges[:-1], 2 * bkg_rel, width=bin_w,
            bottom=1 - bkg_rel, align="edge",
            color="gray", alpha=0.4,
        )
        ax_ratio.errorbar(
            dat_centers, ratio, yerr=ratio_err,
            fmt="o", color="black", markersize=4,
            linewidth=1, capsize=2, zorder=5,
        )
        ax_ratio.axhline(1, color="gray", linestyle="--", linewidth=1)

    ax_ratio.set_xlabel(VAR_LABEL)
    ax_ratio.set_ylabel("Data / MC")
    ax_ratio.set_ylim(0.5, 1.5)
    ax_ratio.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])

    base = os.path.join(outdir, outname)
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {base}.pdf / .png")


# ─── Produce plots ───────────────────────────────────────────────────────────
# One subdirectory per jet class under OUTDIR_BASE, one stack+ratio plot per
# region inside it — same layout convention as plot_stack_toptag.py's
# toptag_fail/toptag_pass split, just nested under a single parent directory.
for class_id in CLASS_DEF:
    outdir = os.path.join(OUTDIR_BASE, f"class{class_id}")
    os.makedirs(outdir, exist_ok=True)
    for region in REGION_DATASET:
        print(f"Plotting {VAR} in {region} [class {class_id}] ...")
        make_stack_plot(region, class_id, outdir, f"stack_{region}_{VAR}")
    print(f"Done with class {class_id}. Plots saved in {outdir}/")
