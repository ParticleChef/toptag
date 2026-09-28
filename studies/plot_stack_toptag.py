import warnings
warnings.filterwarnings("ignore")

import os
import re
import numpy as np
import matplotlib.pyplot as plt
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
INPUT_FILE = "hadmonotop2023_0702.scaled"
LUMI = 17.96  # 2023preBPix, fb^-1

with open(INPUT_FILE, "rb") as f:
    raw = cloudpickle.loads(lz4.frame.decompress(f.read()))

bkg_data = raw["bkg"]
dat_data = raw["data"]

# ─── Output directory ────────────────────────────────────────────────────────
_stem  = os.path.splitext(os.path.basename(INPUT_FILE))[0]
_year  = re.search(r'\d{4}[A-Za-z]*', _stem).group()
_last4 = _stem.rsplit("_", 1)[-1]
OUTDIR_BASE = f"plots_{_year}_{_last4}"

# Top-tag categories: TvsQCD axis has edges [0, 0.33, 1] -> bin 0 = fail, bin 1 = pass
TAG_CATEGORIES = {
    "toptag_fail": 0,
    "toptag_pass": 1,
}

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

# Per-region stack order (bottom → top); VV = WW+WZ+ZZ combined
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

# TvsQCD is excluded here: it IS the tagger score being cut on, so it has no
# separate TvsQCD-cut axis to slice and cannot be split into pass/fail itself.
PLOT_VARS = [
    "ut", "ut_bin2",
    "fjpt", "fjeta", "nfjets",
    "jpt", "jeta", "njets", "njbtagL",
    "mupt", "mueta",
    "phopt", "phoeta",
    "elept", "eleeta",
    "mindphirecoil", "minDphirecoil",
]

VAR_LABEL = {
    "ut":            r"$p_{\mathrm{T}}^{\mathrm{miss}}$ [GeV]",
    "ut_bin2":       r"$p_{\mathrm{T}}^{\mathrm{miss}}$ [GeV]",
    "fjpt":          r"AK15 jet $p_{\mathrm{T}}$ [GeV]",
    "fjeta":         r"AK15 jet $\eta$",
    "nfjets":        "Number of AK15 jets",
    "jpt":           r"Leading AK4 jet $p_{\mathrm{T}}$ [GeV]",
    "jeta":          r"Leading AK4 jet $\eta$",
    "njets":         "Number of AK4 jets",
    "njbtagL":       "Number of b-tagged jets (loose)",
    "mupt":          r"Muon $p_{\mathrm{T}}$ [GeV]",
    "mueta":         r"Muon $\eta$",
    "phopt":         r"Photon $p_{\mathrm{T}}$ [GeV]",
    "phoeta":        r"Photon $\eta$",
    "elept":         r"Electron $p_{\mathrm{T}}$ [GeV]",
    "eleeta":        r"Electron $\eta$",
    "mindphirecoil": r"min $\Delta\phi$(jets, recoil)",
    "minDphirecoil": r"min $\Delta\phi$(jets, recoil)",
}

TAG_LABEL = {
    "toptag_fail": "fail",
    "toptag_pass": "pass",
}


def slice_region_cat(h, region, cat_idx):
    """Select region, and the TvsQCD-cut bin if this histogram has that axis."""
    sel = {"region": region}
    if "TvsQCD" in h.axes.name and len(h.axes) >= 3:
        sel["TvsQCD"] = cat_idx
    return h[sel]


def get_bkg_arrays(var, region, cat_idx):
    """Return list of (label, values, edges) per stack entry; VV is WW+WZ+ZZ summed."""
    result = []
    for proc in BKG_ORDER[region]:
        if proc == "VV":
            vv_vals, vv_edges = None, None
            for vv in VV_PROCS:
                if vv not in bkg_data.get(var, {}):
                    continue
                h = bkg_data[var][vv]
                if region not in list(h.axes["region"]):
                    continue
                hs = slice_region_cat(h, region, cat_idx)
                v = np.where(np.isfinite(hs.values()), hs.values(), 0)
                e = np.array(list(hs.axes[0].edges))
                vv_vals = v if vv_vals is None else vv_vals + v
                vv_edges = e
            if vv_vals is not None:
                result.append(("VV", vv_vals, vv_edges))
        else:
            if proc not in bkg_data.get(var, {}):
                continue
            h = bkg_data[var][proc]
            if region not in list(h.axes["region"]):
                continue
            hs = slice_region_cat(h, region, cat_idx)
            vals = np.where(np.isfinite(hs.values()), hs.values(), 0)
            edges = np.array(list(hs.axes[0].edges))
            result.append((proc, vals, edges))
    return result


def get_data_arrays(var, region, cat_idx):
    """Return (values, bin centers, errors) for the data histogram."""
    dataset = REGION_DATASET[region]
    if dataset not in dat_data.get(var, {}):
        return None, None, None
    h = dat_data[var][dataset]
    if region not in list(h.axes["region"]):
        return None, None, None
    hs = slice_region_cat(h, region, cat_idx)
    vals = hs.values()
    edges = np.array(list(hs.axes[0].edges))
    centers = 0.5 * (edges[:-1] + edges[1:])
    return vals, centers, np.sqrt(vals)


def make_stack_plot(var, region, cat_name, cat_idx, outdir, outname):
    bkg_list = get_bkg_arrays(var, region, cat_idx)
    dat_vals, dat_centers, dat_errs = get_data_arrays(var, region, cat_idx)

    if not bkg_list:
        print(f"  Skipping {var}/{region}: no bkg data")
        return

    # Skip if histogram is empty (handles e/mu/pho in irrelevant regions)
    total_sum = sum(np.sum(v) for _, v, _ in bkg_list)
    data_sum  = np.sum(dat_vals) if dat_vals is not None else 0
    if total_sum == 0 and data_sum == 0:
        print(f"  Skipping {var}/{region}: empty")
        return

    edges = bkg_list[0][2]
    bkg_vals_list = [v for _, v, _ in bkg_list]
    bkg_cumsum    = np.cumsum(bkg_vals_list, axis=0)
    total_bkg     = bkg_cumsum[-1]

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
    ax_main.set_ylim(1e-2, 1e6)
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
        title=f"{REGION_LABEL[region]}, {TAG_LABEL[cat_name]}",
        title_fontsize=10,
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

    ax_ratio.set_xlabel(VAR_LABEL.get(var, var))
    ax_ratio.set_ylabel("Data / MC")
    ax_ratio.set_ylim(0.5, 1.5)
    ax_ratio.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])

    base = os.path.join(outdir, outname)
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {base}.pdf / .png")


# ─── Produce plots ───────────────────────────────────────────────────────────
for cat_name, cat_idx in TAG_CATEGORIES.items():
    outdir = f"{OUTDIR_BASE}_{cat_name}"
    os.makedirs(outdir, exist_ok=True)
    for region in REGION_DATASET:
        for var in PLOT_VARS:
            if var not in bkg_data:
                print(f"  Skipping {var}: not in data file")
                continue
            print(f"Plotting {var} in {region} [{cat_name}] ...")
            make_stack_plot(var, region, cat_name, cat_idx, outdir, f"stack_{region}_{var}")

    print(f"\nDone with {cat_name}. Plots saved in {outdir}/")
