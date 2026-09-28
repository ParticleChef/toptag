"""fjpt Data/MC stack plots, TvsQCD-inclusive, one per control region.

The TvsQCD axis is SUMMED over both bins (fail + pass), i.e. no top-tagging
requirement is applied -- this is the pre-tag AK15 leading-jet pT spectrum.
No jet-class split: `fjpt` has no class axes, and the class-split histogram
`fjpt_cat` is background-only, so a data/MC comparison is only possible here.

The fjpt axis ends in a very wide 600-2000 GeV bin. It is drawn as an OVERFLOW
BIN closing at X_MAX: bin CONTENTS are untouched, only the drawn right edge
moves. The last bin therefore represents all jets with fjpt >= 600 GeV. Same
convention as the alpha and efficiency plots in this series.

Same style and layout as plot_fjpt_toptagTEST.py, which makes the same plots
split by top-tag category (toptag_fail / toptag_pass). This script is the
TvsQCD-inclusive counterpart and is otherwise deliberately identical, so the
three sets can be read side by side.

Run with the plot_env interpreter (the system Python 3.7 cannot read these
files -- pickle protocol 5):
  /Users/hongjieun/opt/anaconda3/envs/plot_env/bin/python plot_fjpt_inclusive_toptagTEST.py
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

LUMI = 7.99  # 2022preEE, fb^-1 (both input files are 2022preEE)

# fjpt (region, fjpt, TvsQCD) is the data/MC comparison histogram -- it has
# both MC and real MET/EGamma data. fjpt_cat additionally carries n_qinfj /
# n_binfj (truth quark / b-flavor matching) and is bkg-only; it is not used
# here and is left untouched.
VAR = "fjpt"
VAR_LABEL = r"AK15 jet $p_{\mathrm{T}}$ [GeV]"

# Right edge the final 600-2000 GeV bin is DRAWN at; contents are unchanged.
X_MAX = 1000.0

# Both test files are plotted, into separate directories. They are DIFFERENT
# MC test samples with different process lists -- (b) has QCD Multijet and
# G + Jets, (a) does not. Compare each against its own data, never their
# absolute yields against each other.
INPUT_FILES = {
    "hadmonotop2022_toptagTEST.scaled":  "plots_2022_toptagTEST_0922_fjpt_incl",
    "hadmonotop2022_toptagTEST2.scaled": "plots_2022_toptagTEST2_0922_fjpt_incl",
}

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

# Per-region stack order (bottom -> top); VV = WW+WZ+ZZ combined.
# Only processes actually present in the file are drawn (see get_bkg_arrays),
# so this list is safe to keep for a file that lacks some of them.
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

VV_PROCS = ["WW", "WZ", "ZZ"]  # merged into one "VV" stack entry

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


def slice_region_inclusive(h, region):
    """Select one region and SUM the TvsQCD axis (pre-tag, no tagging cut)."""
    sel = {"region": region}
    if "TvsQCD" in h.axes.name:
        sel["TvsQCD"] = sum
    return h[sel]


def get_bkg_arrays(bkg_data, region):
    """Return [(label, values, variances, edges)] per stack entry, VV summed.

    Only processes actually present in the file are included, so this works
    unchanged on a file that lacks QCD Multijet / G + Jets.
    """
    result = []
    for proc in BKG_ORDER[region]:
        procs = VV_PROCS if proc == "VV" else [proc]
        vals = vars_ = edges = None
        for name in procs:
            if name not in bkg_data.get(VAR, {}):
                continue
            h = bkg_data[VAR][name]
            if region not in list(h.axes["region"]):
                continue
            hs = slice_region_inclusive(h, region)
            v = np.where(np.isfinite(hs.values()), hs.values(), 0)
            e = np.where(np.isfinite(hs.variances()), hs.variances(), 0)
            vals = v if vals is None else vals + v
            vars_ = e if vars_ is None else vars_ + e
            edges = np.array(list(hs.axes[0].edges))
        if vals is not None:
            result.append((proc, vals, vars_, edges))
    return result


def get_data_arrays(dat_data, region):
    """Return (values, errors) for the data histogram, or two Nones.

    Bin centers are NOT returned: they are recomputed in make_stack_plot from
    the drawn edges, so the last data point sits inside the overflow bin as it
    is drawn rather than at the true 1300 GeV center.
    """
    dataset = REGION_DATASET[region]
    if dataset not in dat_data.get(VAR, {}):
        return None, None
    h = dat_data[VAR][dataset]
    if region not in list(h.axes["region"]):
        return None, None
    hs = slice_region_inclusive(h, region)
    vals = hs.values()
    return vals, np.sqrt(vals)


def make_stack_plot(bkg_data, dat_data, region, outdir, outname):
    """TvsQCD-inclusive data/MC stack + ratio plot of fjpt for one region."""
    bkg_list = get_bkg_arrays(bkg_data, region)
    dat_vals, dat_errs = get_data_arrays(dat_data, region)

    if not bkg_list:
        print(f"  Skipping {region}: no bkg in this file")
        return

    total_sum = sum(np.sum(v) for _, v, _, _ in bkg_list)
    data_sum = np.sum(dat_vals) if dat_vals is not None else 0
    if total_sum == 0 and data_sum == 0:
        print(f"  Skipping {region}: empty")
        return

    edges = bkg_list[0][3].copy()
    edges[-1] = X_MAX  # draw the last bin as an overflow bin
    dat_centers = 0.5 * (edges[:-1] + edges[1:])
    bkg_cumsum = np.cumsum([v for _, v, _, _ in bkg_list], axis=0)
    total_bkg = bkg_cumsum[-1]
    # MC entries are weighted, so the band comes from the stored sum of
    # squared weights, not sqrt(N).
    bkg_err = np.sqrt(np.sum([e for _, _, e, _ in bkg_list], axis=0))

    fig, (ax_main, ax_ratio) = plt.subplots(
        2, 1, figsize=(8, 8),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
        sharex=True,
    )

    prev = np.zeros(len(total_bkg))
    for (proc, vals, _, _), top in zip(bkg_list, bkg_cumsum):
        ax_main.stairs(
            top, edges, baseline=prev, fill=True,
            color=BKG_COLOR.get(proc, "gray"),
            edgecolor="black", linewidth=1.5,
            label=BKG_LEGEND.get(proc, proc),
        )
        prev = top.copy()

    ax_main.stairs(
        total_bkg + bkg_err, edges, baseline=total_bkg - bkg_err,
        fill=True, color="none", edgecolor="gray", hatch="///",
        label="Bkg unc.", linewidth=0,
    )

    if dat_vals is not None:
        ax_main.errorbar(
            dat_centers, dat_vals, yerr=dat_errs,
            fmt="o", color="black", markersize=4,
            linewidth=1, capsize=2, label="Data", zorder=5,
        )

    ax_main.set_yscale("log")

    # Decade-rounded ylim + decade-only major ticks: avoids matplotlib's
    # "6 x 10^0"-style multiplier tick labels on narrow-range log axes, which
    # widen the left margin and push the CMS header into the lumi label.
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

    # Legend: Data first, bkg reversed (top of stack -> bottom), Bkg unc. last
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
        [handles[i] for i in order], [labels[i] for i in order],
        title=f"{REGION_LABEL[region]}, TvsQCD inclusive\n"
              r"last bin: $p_{\mathrm{T}} \geq 600$ GeV",
        title_fontsize=10, ncol=2, fontsize=9, framealpha=0.8,
        loc="upper right",
    )
    leg.get_title().set_ha("center")
    leg.get_title().set_fontweight("bold")

    hep.cms.label("Private Work", data=True, lumi=LUMI, com=13.6, ax=ax_main)

    if dat_vals is not None:
        ratio = np.where(total_bkg > 0, dat_vals / total_bkg, np.nan)
        ratio_err = np.where(total_bkg > 0, dat_errs / total_bkg, np.nan)
        bkg_rel = np.where(total_bkg > 0, bkg_err / total_bkg, 0)
        bin_w = edges[1:] - edges[:-1]

        ax_ratio.bar(
            edges[:-1], 2 * bkg_rel, width=bin_w,
            bottom=1 - bkg_rel, align="edge", color="gray", alpha=0.4,
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

    ratio_txt = "no data"
    if dat_vals is not None and total_sum > 0:
        ratio_txt = f"Data/MC = {data_sum / total_sum:.4g}"
    print(f"  {region:5s}  N(Data)={data_sum:9.0f}  N(MC)={total_sum:11.5g}"
          f"  {ratio_txt}   -> {base}.pdf/.png")


# ─── Produce plots ───────────────────────────────────────────────────────────
for input_file, outdir in INPUT_FILES.items():
    print(f"{input_file} -> {outdir}/")
    with open(input_file, "rb") as f:
        raw = cloudpickle.loads(lz4.frame.decompress(f.read()))
    os.makedirs(outdir, exist_ok=True)
    for region in REGION_DATASET:
        make_stack_plot(raw["bkg"], raw["data"], region, outdir,
                        f"stack_{region}_{VAR}_incl")
    print()

print("Done.")
