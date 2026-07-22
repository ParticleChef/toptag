#!/usr/bin/env python
"""validate.py - Step 3 (optional) of the top-tagging SF pipeline (see README.md).

Data/MC comparison of TvsQCD (collapsed to the 2-bin below/above-WP view, per
AN-2020/061 Fig. 29/30) in gcr/tmcr/tecr, before and after applying the
top-tagging calibration from scalefactors.py - a quick visual check that the
calibration corrects the pre-tag/post-tag data/MC discrepancy it was derived
to fix.

Usage:
    python validate.py --year 2022pre --hists hists/toptag_hists_2022pre.coffea \
        --sf data/toptag_sf_2022pre.json.gz --outdir plots/validate_2022pre/
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from coffea.util import load

from scalefactors import PT_EDGES, sum_processes, tag_split, GJ_PROCESS  # noqa: F401

# ============================================================
# CONFIG
# ============================================================

LUMIS = {"2022pre": 7.99, "2022post": 26.68, "2023pre": 17.96, "2023post": 9.68}
REGION_DATA = {"gcr": "EGamma", "tmcr": "MET", "tecr": "EGamma"}
REGION_LABEL = {"gcr": r"$\gamma$-enriched CR", "tmcr": r"t$\bar{t}$-enriched CR ($\mu$)", "tecr": r"t$\bar{t}$-enriched CR (e)"}
PT_CENTERS = [0.5 * (PT_EDGES[i] + min(PT_EDGES[i + 1], 1000)) for i in range(len(PT_EDGES) - 1)]


def mc_yields_by_class(bkg, region):
    """dict[jetclass] -> (tagged, untagged) arrays over pT bins."""
    mc = sum_processes(bkg)
    out = {}
    for jetclass in ("top", "b", "qcd"):
        tagged, _, total, _ = tag_split(mc, region, jetclass=jetclass)
        out[jetclass] = (tagged, total - tagged)
    return out


def calibrated_mc(bkg, region, cset):
    sf_corr = cset["toptag_sf"]
    eff_corr = cset["toptag_efficiency"]
    yields = mc_yields_by_class(bkg, region)
    tagged_total = np.zeros(len(PT_CENTERS))
    untagged_total = np.zeros(len(PT_CENTERS))
    for jetclass, (tagged, untagged) in yields.items():
        for i, pt in enumerate(PT_CENTERS):
            sf = sf_corr.evaluate(pt, jetclass, "nominal")
            eff = eff_corr.evaluate(pt, jetclass)
            tagged_total[i] += tagged[i] * sf
            denom = 1 - eff
            untagged_total[i] += untagged[i] * ((1 - eff * sf) / denom) if denom > 0 else untagged[i]
    return tagged_total, untagged_total


def raw_mc(bkg, region):
    yields = mc_yields_by_class(bkg, region)
    tagged_total = sum(t for t, _ in yields.values())
    untagged_total = sum(u for _, u in yields.values())
    return tagged_total, untagged_total


def plot_region(region, data_h, mc_before, mc_after, year, outdir):
    data_tagged, _, data_total, data_var = tag_split(data_h, region)
    data_untagged = data_total - data_tagged

    hep.style.use("CMS")
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), gridspec_kw={"height_ratios": [3, 1]}, sharex="col")
    for col, (mc_tagged, mc_untagged, title) in enumerate(
        [(mc_before[0], mc_before[1], "before calibration"), (mc_after[0], mc_after[1], "after calibration")]
    ):
        ax, rax = axes[0, col], axes[1, col]
        mc_stack = np.array([mc_untagged.sum(), mc_tagged.sum()])
        data_stack = np.array([data_untagged.sum(), data_tagged.sum()])
        x = np.array([0, 1])
        ax.stairs(mc_stack, edges=[-0.5, 0.5, 1.5], fill=True, color="#1f77b4", edgecolor="black", linewidth=1.5, label="MC")
        ax.errorbar(x, data_stack, yerr=np.sqrt(data_stack), fmt="ko", label="Data")
        ax.set_yscale("log")
        ax.set_title(title, fontsize=11)
        ax.legend()
        ratio = np.divide(data_stack, mc_stack, out=np.full_like(data_stack, np.nan, dtype=float), where=mc_stack > 0)
        rax.errorbar(x, ratio, fmt="ko")
        rax.axhline(1.0, color="gray", linestyle="--")
        rax.set_ylim(0.5, 1.5)
        rax.set_xticks(x)
        rax.set_xticklabels(["untagged", "tagged"])
        rax.set_ylabel("Data/MC")
    hep.cms.label(ax=axes[0, 0], label="Preliminary", data=True, lumi=LUMIS[year], year=year, com=13.6)
    fig.suptitle(REGION_LABEL[region])
    fig.tight_layout()

    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"toptag_validate_{region}_{year}.{ext}"), dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", required=True)
    parser.add_argument("--hists", required=True)
    parser.add_argument("--sf", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    import correctionlib

    hists = load(args.hists)
    bkg, data = hists["bkg"], hists["data"]
    cset = correctionlib.CorrectionSet.from_file(args.sf)

    for region, data_process in REGION_DATA.items():
        data_h = data[data_process] if data_process in data else sum_processes(data)
        mc_before = raw_mc(bkg, region)
        mc_after = calibrated_mc(bkg, region, cset)
        plot_region(region, data_h, mc_before, mc_after, args.year, args.outdir)
        print(f"Saved plots for {region}")


if __name__ == "__main__":
    main()
