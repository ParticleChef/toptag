#!/usr/bin/env python3
"""
plot_roc_from_arrays.py

Recompute the ParT TvsQCD ROC and working-point table from the score arrays
saved by roc_curve.py (--save-arrays), on the laptop (plot_env).

roc_curve.py's compute_roc() pairs the efficiencies with a reversed threshold
array, so its printed thresholds are mirrored (t -> 1 - 0.001 - t) and its AUC
is 1 - AUC. This script pairs each efficiency with its own threshold.

Usage:
    $PY studies/plot_roc_from_arrays.py hists/roc_arrays_2022pre_sr.coffea --tag sr
    $PY studies/plot_roc_from_arrays.py hists/roc_arrays_2022pre.coffea --tag incl

Outputs (plots_roc/recomputed/):
    roc_TvsQCD_<year>_<tag>_{eff,rej}.{pdf,png}   ROC vs bkg efficiency / rejection
    score_TvsQCD_<year>_<tag>.{pdf,png}           normalized TvsQCD, signal vs bkg
"""

import argparse
import os
import pickle

import lz4.frame
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep

# ============================================================
# CONFIG
# ============================================================

WORK_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = "plots_roc/recomputed"
NBINS = 1000                                    # same as roc_curve.py
BKG_EFF_WPS = [0.05, 0.025, 0.01, 0.005, 0.001]
WP_COLORS = ["#e41a1c", "#ff7f00", "#4daf4a", "#984ea3", "#377eb8"]
FIXED_CUTS = [0.26, 0.33]                       # Run 2 WP, cut used now
LUMI = {"2022pre": 7.98, "2022post": 26.67, "2023pre": 17.96, "2023post": 9.68}


class _Acc:
    """Stand-in for coffea accumulators, so coffea.processor is not imported."""

    def __setstate__(self, state):
        self.__dict__.update(state if isinstance(state, dict) else {"_value": state})

    @property
    def value(self):
        return self.__dict__.get("_value")


class _Unpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module.startswith("coffea"):
            return dict if name == "dict_accumulator" else _Acc
        return super().find_class(module, name)


def load_arrays(path):
    """Read a roc_curve.py --save-arrays file; return (sig, sigw, bkg, bkgw, year)."""
    with lz4.frame.open(path) as f:
        bundle = _Unpickler(f).load()
    out = bundle["out"]
    arrs = [out[k].value for k in ("sig_scores", "sig_weights", "bkg_scores", "bkg_weights")]
    return (*arrs, bundle.get("year", "2022pre"))


def pass_fraction(scores, weights, edges):
    """Weighted fraction of events with score >= edges[i], for each lower edge i."""
    ok = scores >= 0.                           # -1 = invalid jet sentinel
    h, _ = np.histogram(scores[ok], bins=edges, weights=weights[ok])
    return np.cumsum(h[::-1])[::-1] / weights[ok].sum()


def compute_roc(sig, sigw, bkg, bkgw):
    """Return (thresholds, sig_eff, bkg_eff, auc); index i means cut TvsQCD >= thr[i]."""
    edges = np.linspace(0., 1., NBINS + 1)
    se = pass_fraction(sig, sigw, edges)
    be = pass_fraction(bkg, bkgw, edges)
    auc = float(np.trapz(np.r_[se, 0.][::-1], np.r_[be, 0.][::-1]))
    return edges[:-1], se, be, auc


def wp_table(thr, se, be):
    """Rows (target, thr, sig_eff, bkg_eff) at the target bkg efficiencies and fixed cuts."""
    rows = []
    for t in BKG_EFF_WPS:
        i = int(np.argmin(np.abs(be - t)))
        rows.append((f"{t * 100:.1f}%", thr[i], se[i], be[i]))
    for c in FIXED_CUTS:
        i = int(round(c * NBINS))
        rows.append((f"cut {c}", thr[i], se[i], be[i]))
    return rows


def print_table(rows, auc, nsig, nbkg):
    """Print the working-point table in the roc_curve.py format."""
    print(f"  Signal events: {nsig:,}   Background events: {nbkg:,}")
    print(f"  AUC = {auc:.4f}\n")
    print(f"  {'Target':>10}  {'Threshold':>9}  {'eff_sig':>7}  {'eff_bkg':>8}  {'Rejection':>9}")
    for name, t, s, b in rows:
        print(f"  {name:>10}  {t:>9.3f}  {s:>7.3f}  {b * 100:>7.2f}%  {1 / b:>8.1f}x")


def _style(ax, year, xlabel):
    """Common axis labels and CMS header."""
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"Signal efficiency ($m_\phi=1000,\ m_\chi=150$)")
    ax.set_ylim(0., 1.)
    hep.cms.label(ax=ax, text="Private Work", data=False, lumi=LUMI.get(year), com=13.6, fontsize=18)


def _save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUTDIR, f"{name}.{ext}"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {OUTDIR}/{name}.{{pdf,png}}")


def plot_roc(thr, se, be, auc, rows, year, tag):
    """ROC vs bkg efficiency and vs bkg rejection, with the WP markers."""
    wps = rows[:len(BKG_EFF_WPS)]
    for kind in ("eff", "rej"):
        fig, ax = plt.subplots(figsize=(10, 9))
        x = be if kind == "eff" else np.where(be > 0, 1. / np.maximum(be, 1e-12), np.nan)
        ax.plot(x, se, lw=2, color="steelblue", label=f"ParT TvsQCD (AUC = {auc:.3f})")
        for (name, t, s, b), col in zip(wps, WP_COLORS):
            xb = b if kind == "eff" else 1. / b
            ax.scatter([xb], [s], color=col, s=60, zorder=5,
                       label=rf"$\varepsilon_{{bkg}}$={name}: TvsQCD > {t:.3f}, $\varepsilon_{{sig}}$={s:.3f}")
        if kind == "eff":
            ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5, label="Random")
            ax.set_xlim(0., 1.)
            _style(ax, year, r"Background efficiency ($Z\to\nu\nu$)")
            ax.legend(fontsize=13, loc="lower right")
        else:
            ax.set_xscale("log")
            _style(ax, year, r"Background rejection ($Z\to\nu\nu$)")
            ax.legend(fontsize=13, loc="upper right")
        _save(fig, f"roc_TvsQCD_{year}_{tag}_{kind}")


def plot_scores(sig, sigw, bkg, bkgw, rows, year, tag):
    """Normalized TvsQCD distributions of signal and background, with cut lines."""
    edges = np.linspace(0., 1., 51)
    fig, ax = plt.subplots(figsize=(10, 9))
    for s, w, lab, col in ((sig, sigw, "Signal (mono-top)", "#ff7f0e"),
                           (bkg, bkgw, r"Background ($Z\to\nu\nu$)", "#1f77b4")):
        ok = s >= 0.
        h, _ = np.histogram(s[ok], bins=edges, weights=w[ok])
        hep.histplot(h / h.sum(), edges, ax=ax, lw=2, color=col, label=lab)
    lines = [(rows[2][1], "-", r"$\varepsilon_{bkg}$ = 1%"), (0.33, ":", "TvsQCD = 0.33")]
    for x, ls, lab in lines:
        ax.axvline(x, color="k", ls=ls, lw=1.2, label=f"{lab} ({x:.3f})")
    ax.set_yscale("log")
    ax.set_xlabel("ParT TvsQCD (leading AK15 jet)")
    ax.set_ylabel("Fraction of events")
    hep.cms.label(ax=ax, text="Private Work", data=False, lumi=LUMI.get(year), com=13.6, fontsize=18)
    ax.legend(fontsize=14, loc="upper center")
    _save(fig, f"score_TvsQCD_{year}_{tag}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("arrays", help="roc_curve.py --save-arrays file")
    p.add_argument("--tag", default="sr", help="suffix for the output names (default: sr)")
    args = p.parse_args()

    path = os.path.abspath(args.arrays)
    os.chdir(WORK_DIR)
    os.makedirs(OUTDIR, exist_ok=True)
    hep.style.use("CMS")

    sig, sigw, bkg, bkgw, year = load_arrays(path)
    thr, se, be, auc = compute_roc(sig, sigw, bkg, bkgw)
    rows = wp_table(thr, se, be)
    print(f"== {path}")
    print_table(rows, auc, len(sig), len(bkg))
    plot_roc(thr, se, be, auc, rows, year, args.tag)
    plot_scores(sig, sigw, bkg, bkgw, rows, year, args.tag)


if __name__ == "__main__":
    main()
