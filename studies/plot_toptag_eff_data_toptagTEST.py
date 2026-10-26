"""Data-side top-tag efficiency, and the resulting scale factor.

Per toptag/README.md ("Scale Factor Derivation"), the data efficiency for an
aggregated jet class is the data yield with the OTHER classes subtracted, each
scaled to data by the pre-tag normalization factor alpha = N(Data)/N(MC):

  QCD-jet (gamma CR, gamma+jets sample, B-jet subtracted):
    eff_Data = [N(Data,tag) - alpha*SF_B*N(GJ,B,tag)]
             / [N(Data)     - alpha*N(GJ,B)]

  Top-jet (ttbar CRs, B-jet and QCD-jet subtracted):
    eff_Data = [N(Data,tag) - alpha*(SF_B*N(B,tag) + SF_QCD*N(QCD,tag))]
             / [N(Data)     - alpha*(N(B)          + N(QCD))]

  SF_B / SF_QCD correct the TAG RATE of the subtracted classes, so they
  multiply only the tagged terms; the pre-tag terms are fixed by alpha alone.

  SF = eff_Data / eff_MC

Data yields and the all-class MC total (alpha) come from `fjpt_bin`, which has
data; the class-split MC yields come from `fjpt_cat`, which does not -- truth
matching is undefined for data. Both have the same pT edges, and fjpt_bin MC
must equal fjpt_cat MC summed over its class axes (checked on every run), so
they combine bin-by-bin. (Old 2022 test files: `--alpha-var fjpt`.)

alpha is taken PER CONTROL REGION, per pT bin -- see plot_alpha_toptagTEST.py
and report_alpha_toptagTEST.txt.

Usage:
  python studies/plot_toptag_eff_data_toptagTEST.py -f hadmonotop2023_toptagTEST4 \\
      [-e ERA] [--mc-eff] [--data-eff]
Always writes the eff_MC + eff_Data + SF plots; --mc-eff / --data-eff add
single-efficiency plots without the SF panel.
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
import sf_json

hep.style.use("CMS")

# Input file, era and output base from the command line (see effsf_cli.py).
# setup() also changes into toptag/, so relative paths below resolve there.
_parser = base_parser(__doc__.splitlines()[0])
_parser.add_argument("--mc-eff", action="store_true",
                     help="also write MC-only efficiency plots (<name>_mc)")
_parser.add_argument("--data-eff", action="store_true",
                     help="also write Data-only efficiency plots (<name>_data)")
_parser.add_argument("--json", metavar="PATH",
                     help="correctionlib SF file (default:"
                          " plots_effSF_<version>/sf_json/toptag_sf_<era>.json)")
ARGS = _parser.parse_args()
CFG = setup(ARGS)
FILE, ERA, LUMI = CFG.file, CFG.era, CFG.lumi
# Histogram the data-side yields (alpha, N(Data)) are read from: `fjpt_bin`,
# binned exactly like fjpt_cat (old 2022 test files: `fjpt`).
ALPHA_VAR = ARGS.alpha_var

VAR_LABEL = r"AK15 jet $p_{\mathrm{T}}$ [GeV]"
OUTDIR = os.path.join(CFG.outbase, "effdata", ERA)  # plots_effSF_<version>/effdata/<era>/
X_MAX = 1000.0  # the 600-2000 GeV bin is drawn as an overflow bin closing here

TAG_PASS = 1  # TvsQCD axis edges [0, 0.33, 1] -> bin 1 = tagged (WP = 0.33)
WP = 0.33

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
AGG_GROUPS = {"QCD-jet": [1, 3], "B-jet": [2, 4], "Top-jet": [5, 6, 7, 8]}

REGION_DATASET = {
    "sr": "MET", "wmcr": "MET", "tmcr": "MET", "zmcr": "MET",
    "wecr": "EGamma", "tecr": "EGamma", "gcr": "EGamma", "zecr": "EGamma",
}

# AN Table 31/32 inputs that are not histogram statistics
SF_B_CENTRAL = 1.0
SF_B_UNCERTAINTY = 0.5          # no dedicated CR exists for the B-jet class
GJ_BJET_NORM_UNCERTAINTY = 0.5  # 50% norm. uncertainty on the GJ B-jet MC

# Plot colors, shared with plot_toptag_eff_toptagTEST.py (MC_COLOR)
MC_COLOR = "#1f77b4"
DATA_COLOR = "#ff7f0e"
SF_COLOR = "black"


_cache = {}


def load(filename):
    """Load and cache one .scaled file."""
    if filename not in _cache:
        with open(filename, "rb") as f:
            _cache[filename] = cloudpickle.loads(lz4.frame.decompress(f.read()))
    return _cache[filename]


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


def data_side_hists(raw):
    """Return the (bkg, data) dicts of ALPHA_VAR, after checking they are usable.

    Its yields are combined bin-by-bin with the class-split fjpt_cat MC, so it
    must be filled and must have exactly the fjpt_cat pT edges.
    """
    if ALPHA_VAR not in raw["bkg"] or not any(
            len(h.axes["region"]) for h in raw["bkg"][ALPHA_VAR].values()):
        raise SystemExit(f"'{ALPHA_VAR}' is missing or empty in {os.path.basename(FILE)}"
                         " -- it must be filled by the processor")
    cat_edges = next(iter(raw["bkg"]["fjpt_cat"].values())).axes["fjpt"].edges
    for kind in ("bkg", "data"):
        for name, h in raw[kind][ALPHA_VAR].items():
            edges = h.axes[ALPHA_VAR].edges
            if len(edges) != len(cat_edges) or not np.allclose(edges, cat_edges):
                raise SystemExit(f"{kind}/{ALPHA_VAR}/{name}: pT edges {edges}"
                                 f" differ from fjpt_cat {cat_edges}")
    return raw["bkg"][ALPHA_VAR], raw["data"][ALPHA_VAR]


def _accumulate(hs, vals, vars_):
    """Add one sliced histogram's values/variances into the running sums."""
    v, e = hs.values(), hs.variances()
    vals += np.where(np.isfinite(v), v, 0)
    vars_ += np.where(np.isfinite(e), e, 0)


def mc_class_yield(bkg_cat, group, regions, procs, tvsqcd):
    """(values, variances) vs fjpt for one aggregated class, from fjpt_cat.

    Summed over regions, over `procs` (all processes if None) and over the
    constituent classes -- the 8 classes strictly partition the
    (n_qinfj, n_binfj) plane, so this is a plain sum with no double counting.
    """
    names = list(bkg_cat) if procs is None else [p for p in procs if p in bkg_cat]
    vals = vars_ = None
    for proc in names:
        h = bkg_cat[proc]
        if vals is None:
            n = len(h.axes["fjpt"].edges) - 1
            vals, vars_ = np.zeros(n), np.zeros(n)
        for region in regions:
            if region not in list(h.axes["region"]):
                continue
            for class_id in AGG_GROUPS[group]:
                sel = {"region": region, "TvsQCD": tvsqcd}
                sel.update(class_selector(class_id))
                _accumulate(h[sel], vals, vars_)
    return vals, vars_


def mc_total_yield(bkg_bin, regions, tvsqcd):
    """(values, variances) vs pT over ALL processes and classes, from ALPHA_VAR."""
    vals = vars_ = None
    for proc, h in bkg_bin.items():
        if vals is None:
            n = len(h.axes[ALPHA_VAR].edges) - 1
            vals, vars_ = np.zeros(n), np.zeros(n)
        for region in regions:
            if region not in list(h.axes["region"]):
                continue
            _accumulate(h[{"region": region, "TvsQCD": tvsqcd}], vals, vars_)
    return vals, vars_


def data_yield(data_bin, regions, tvsqcd):
    """(values, variances, edges) vs pT from ALPHA_VAR, each region's own dataset."""
    vals = vars_ = edges = None
    for region in regions:
        h = data_bin[REGION_DATASET[region]]
        if vals is None:
            edges = h.axes[ALPHA_VAR].edges
            vals, vars_ = np.zeros(len(edges) - 1), np.zeros(len(edges) - 1)
        _accumulate(h[{"region": region, "TvsQCD": tvsqcd}], vals, vars_)
    return vals, vars_, edges


def propagate(formula, *args_with_unc):
    """Linearized error propagation by one-at-a-time finite differences.

    Shifts each input by its own 1-sigma (holding the rest at their central
    value) and adds the resulting shifts in quadrature. Inputs are treated as
    uncorrelated, which is an approximation for terms drawn from overlapping
    event samples (a tagged yield vs. its own total). Same convention, and
    same caveat, as toptag/scalefactors.py's propagate().
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


# ─── Measurements ────────────────────────────────────────────────────────────

def measure(filename, group, regions, procs, sf_qcd=None):
    """Return a dict of everything measured for one (class, region set).

    `group` is the target class; every OTHER class listed in `subtract` is
    removed from data, scaled by alpha. `procs` restricts the MC subtraction
    terms to a process list (the QCD-jet calibration uses G + Jets alone,
    the AN's gamma+jets-dominance approximation); alpha itself is ALWAYS
    computed over every process in the region.

    `sf_qcd` is the (value, sigma) pair from a previous QCD-jet measurement.
    When given, the QCD-jet subtraction term is corrected by it -- see
    toptag/docs/scalefactors.md on why the implemented formula uses SF_QCD
    here even though the README's boxed formula does not.
    """
    raw = load(filename)
    bkg_bin, data_bin = data_side_hists(raw)
    bkg_cat = raw["bkg"]["fjpt_cat"]

    d_tag, vd_tag, edges = data_yield(data_bin, regions, TAG_PASS)
    d_tot, vd_tot, _ = data_yield(data_bin, regions, sum)
    m_tot, vm_tot = mc_total_yield(bkg_bin, regions, sum)

    tgt_tag, vtgt_tag = mc_class_yield(bkg_cat, group, regions, procs, TAG_PASS)
    tgt_tot, vtgt_tot = mc_class_yield(bkg_cat, group, regions, procs, sum)

    # ALPHA_VAR and fjpt_cat must describe the same MC events (the 2023 test
    # file double-filled fjpt, which this catches).
    cat_all = sum(mc_class_yield(bkg_cat, g, regions, None, sum)[0] for g in AGG_GROUPS)
    if not np.allclose(cat_all, m_tot, rtol=1e-6):
        print(f"  WARNING: {ALPHA_VAR} / fjpt_cat MC mismatch, ratio {m_tot / cat_all}")

    subtract = [g for g in AGG_GROUPS if g != group] if group == "Top-jet" else ["B-jet"]
    terms = {g: dict(zip(("tag", "vtag", "tot", "vtot"),
                         mc_class_yield(bkg_cat, g, regions, procs, TAG_PASS)
                         + mc_class_yield(bkg_cat, g, regions, procs, sum)))
             for g in subtract}

    with np.errstate(divide="ignore", invalid="ignore"):
        alpha = np.where(m_tot > 0, d_tot / m_tot, np.nan)
        eff_mc = np.where(tgt_tot > 0, tgt_tag / tgt_tot, np.nan)

    b = terms["B-jet"]
    # AN Table 31 "BJetNorm": fold the 50% normalization uncertainty on the
    # gamma+jets B-jet MC into that term's effective sigma.
    b_norm = GJ_BJET_NORM_UNCERTAINTY if group == "QCD-jet" else 0.0
    b_tag_sig = np.sqrt(b["vtag"] + (b_norm * b["tag"]) ** 2)
    b_tot_sig = np.sqrt(b["vtot"] + (b_norm * b["tot"]) ** 2)

    if group == "QCD-jet":
        def formula(dt, dT, mT, btag, btot, sfb):
            a = dT / mT
            # SF_B corrects the tag rate only: pre-tag B yield is alpha*btot
            return (dt - a * sfb * btag) / (dT - a * btot)

        inputs = [(d_tag, np.sqrt(vd_tag)), (d_tot, np.sqrt(vd_tot)),
                  (m_tot, np.sqrt(vm_tot)),
                  (b["tag"], b_tag_sig), (b["tot"], b_tot_sig),
                  (np.full_like(d_tag, SF_B_CENTRAL),
                   np.full_like(d_tag, SF_B_UNCERTAINTY))]
    else:
        q = terms["QCD-jet"]
        sf_q, sf_q_sig = sf_qcd if sf_qcd is not None else (
            np.ones_like(d_tag), np.zeros_like(d_tag))

        def formula(dt, dT, mT, btag, btot, qtag, qtot, sfb, sfq):
            a = dT / mT
            # SFs correct the tag rate only: pre-tag yields are alpha*(btot+qtot)
            return ((dt - a * (sfb * btag + sfq * qtag))
                    / (dT - a * (btot + qtot)))

        inputs = [(d_tag, np.sqrt(vd_tag)), (d_tot, np.sqrt(vd_tot)),
                  (m_tot, np.sqrt(vm_tot)),
                  (b["tag"], b_tag_sig), (b["tot"], b_tot_sig),
                  (q["tag"], np.sqrt(q["vtag"])), (q["tot"], np.sqrt(q["vtot"])),
                  (np.full_like(d_tag, SF_B_CENTRAL),
                   np.full_like(d_tag, SF_B_UNCERTAINTY)),
                  (sf_q, sf_q_sig)]

    eff_data, eff_data_err = propagate(formula, *inputs)

    # Statistics-only variant: same formula, only the histogram sigmas kept.
    # The B-jet terms must fall back to their pure sqrt(variance) here, since
    # b_tag_sig / b_tot_sig also carry the BJetNorm systematic.
    stat_inputs = list(inputs)
    stat_inputs[3] = (b["tag"], np.sqrt(b["vtag"]))
    stat_inputs[4] = (b["tot"], np.sqrt(b["vtot"]))
    stat_inputs[-2:] = [(p[0], np.zeros_like(p[0])) for p in inputs[-2:]]
    _, eff_data_stat = propagate(formula, *stat_inputs)

    # eff_MC error: pass/fail are disjoint, same formula as the eff_MC script.
    vtgt_fail = vtgt_tot - vtgt_tag
    with np.errstate(divide="ignore", invalid="ignore"):
        eff_mc_err = np.sqrt((1 - eff_mc) ** 2 * vtgt_tag
                             + eff_mc ** 2 * vtgt_fail) / tgt_tot
        sf = eff_data / eff_mc
        sf_err = np.abs(sf) * np.sqrt((eff_data_err / eff_data) ** 2
                                      + (eff_mc_err / eff_mc) ** 2)

    return dict(edges=edges, alpha=alpha, eff_mc=eff_mc, eff_mc_err=eff_mc_err,
                eff_data=eff_data, eff_data_err=eff_data_err,
                eff_data_stat=eff_data_stat, sf=sf, sf_err=sf_err,
                d_tag=d_tag, d_tot=d_tot, m_tot=m_tot,
                tgt_tag=tgt_tag, tgt_tot=tgt_tot, terms=terms)


# ─── Plot ────────────────────────────────────────────────────────────────────

def make_plot(key, cfg, series, zoom=False, mc_color=MC_COLOR,
              mc_label="MC", sf_color=SF_COLOR):
    """Efficiency (MC + one or more Data variants) with an SF ratio panel.

    zoom=True re-draws the same points on the physical 0-1 efficiency range;
    points outside it simply fall off the axis. sf_color=None draws each SF
    series in its Data series' color (needed when several Data series share
    the ratio panel).
    """
    res0 = series[0]["res"]
    edges = res0["edges"].copy()
    edges[-1] = X_MAX  # draw the last bin as an overflow bin
    centers = 0.5 * (edges[:-1] + edges[1:])
    xerr = 0.5 * (edges[1:] - edges[:-1])

    fig, (ax, rax) = plt.subplots(
        2, 1, figsize=(10, 11), sharex=True,
        gridspec_kw=dict(height_ratios=[3, 1], hspace=0.08))

    n = len(series) + 1
    offs = [(i - (n - 1) / 2) * 0.12 * 2 * xerr for i in range(n)]

    ax.errorbar(centers + offs[0], res0["eff_mc"], yerr=res0["eff_mc_err"], xerr=xerr,
                fmt="o", color=mc_color, markersize=7, linewidth=1.5,
                capsize=0, label=mc_label, zorder=5)
    for i, s in enumerate(series):
        ax.errorbar(centers + offs[i + 1], s["res"]["eff_data"],
                    yerr=s["res"]["eff_data_err"], xerr=xerr,
                    fmt=s["marker"], color=s["color"], markersize=7,
                    linewidth=1.5, capsize=0, label=s["label"], zorder=6 + i)
        rax.errorbar(centers + offs[i + 1], s["res"]["sf"],
                     yerr=s["res"]["sf_err"], xerr=xerr,
                     fmt=s["marker"], color=sf_color or s["color"], markersize=6,
                     linewidth=1.5, capsize=0, zorder=5 + i)

    # An efficiency above this line is unphysical.
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=2.0, zorder=1)

    ax.set_ylim(*(cfg["zoomlim"] if zoom else cfg["ylim"]))
    ax.set_ylabel(cfg["ylabel"])
    ax.grid(axis="y", color="gray", linestyle="-", linewidth=0.6,
            alpha=0.4, zorder=0)
    ax.set_axisbelow(True)

    leg = ax.legend(loc=cfg["legloc"], fontsize=13,
                    title=f"{cfg['label']}\n{cfg['sublabel']}")
    leg.get_title().set_ha("center")
    leg.get_title().set_fontweight("bold")
    leg.get_title().set_fontsize(10)

    rax.axhline(1.0, color="red", linestyle="--", linewidth=1.5, zorder=1)
    rax.set_ylim(*(cfg["zoomsflim"] if zoom else cfg["sflim"]))
    rax.set_ylabel("SF")
    rax.set_xlabel(VAR_LABEL)
    rax.set_xlim(edges[0], edges[-1])
    rax.grid(axis="y", color="gray", linestyle="-", linewidth=0.6,
             alpha=0.4, zorder=0)
    rax.set_axisbelow(True)

    hep.cms.label("Private Work", data=True, lumi=LUMI, com=13.6,
                  ax=ax, fontsize=20)

    os.makedirs(OUTDIR, exist_ok=True)
    base = os.path.join(OUTDIR, key + ("_zoom" if zoom else ""))
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_eff_only_plot(key, cfg, res, which, data_label="Data",
                       data_marker="s", zoom=False):
    """Single efficiency (which = "mc" or "data") vs pT, no SF panel.

    Saved as <key>_<which>[_zoom]. The MC efficiency is always physical, so
    its plot uses the zoomed 0-1 range; the Data one follows `zoom` like
    make_plot.
    """
    if which == "mc":
        eff, err, color, label = res["eff_mc"], res["eff_mc_err"], MC_COLOR, "MC"
        ylim, marker = cfg["zoomlim"], "o"
    else:
        eff, err, color, label = res["eff_data"], res["eff_data_err"], DATA_COLOR, data_label
        ylim, marker = (cfg["zoomlim"] if zoom else cfg["ylim"]), data_marker

    edges = res["edges"].copy()
    edges[-1] = X_MAX  # draw the last bin as an overflow bin
    centers = 0.5 * (edges[:-1] + edges[1:])
    xerr = 0.5 * (edges[1:] - edges[:-1])

    fig, ax = plt.subplots(figsize=(10, 9))
    ax.errorbar(centers, eff, yerr=err, xerr=xerr, fmt=marker, color=color,
                markersize=7, linewidth=1.5, capsize=0, label=label, zorder=5)
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=2.0, zorder=1)

    ax.set_ylim(*ylim)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel(VAR_LABEL)
    ax.set_ylabel(cfg["ylabel"])
    ax.grid(axis="y", color="gray", linestyle="-", linewidth=0.6,
            alpha=0.4, zorder=0)
    ax.set_axisbelow(True)

    leg = ax.legend(loc=cfg["legloc"], fontsize=13,
                    title=f"{cfg['label']}\n{cfg['sublabel']}")
    leg.get_title().set_ha("center")
    leg.get_title().set_fontweight("bold")
    leg.get_title().set_fontsize(10)

    hep.cms.label("Private Work", data=(which == "data"), lumi=LUMI, com=13.6,
                  ax=ax, fontsize=20)

    os.makedirs(OUTDIR, exist_ok=True)
    base = os.path.join(OUTDIR, f"{key}_{which}" + ("_zoom" if zoom else ""))
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def report(cfg, label, res):
    """Print the full per-bin table for one measurement variant."""
    edges = res["edges"].copy()
    edges[-1] = X_MAX
    print(f"  {cfg['label']} / {cfg['group']} / {label}"
          f"  [{os.path.basename(cfg['file'])}]")
    print(f"    {'pT [GeV]':>14s} {'alpha':>10s} {'N(Data)':>9s} {'N(D,tag)':>9s}"
          f" {'eff_MC':>8s} {'eff_Data':>20s} {'SF':>19s}")
    for i in range(len(res["eff_mc"])):
        rng = f"{edges[i]:.0f}-{edges[i+1]:.0f}" + (" (of)" if i == len(res["eff_mc"]) - 1 else "")
        print(f"    {rng:>14s} {res['alpha'][i]:10.4g} {res['d_tot'][i]:9.0f}"
              f" {res['d_tag'][i]:9.0f} {res['eff_mc'][i]:8.4f}"
              f" {res['eff_data'][i]:9.4f} +/- {res['eff_data_err'][i]:6.4f}"
              f" {res['sf'][i]:9.4f} +/- {res['sf_err'][i]:6.4f}")
    stat = res["eff_data_stat"]
    syst = np.sqrt(np.maximum(res["eff_data_err"] ** 2 - stat ** 2, 0))
    print(f"      eff_Data unc. breakdown  stat={np.array2string(stat, precision=4)}"
          f"  syst={np.array2string(syst, precision=4)}")


# ─── Run ─────────────────────────────────────────────────────────────────────
# ORDERING CONSTRAINT: SF_QCD must be derived before the Top-jet measurement,
# which uses it to correct its QCD-jet subtraction term.

JOBS_QCD = {
    "effdata_qcdjet_gcr": dict(
        file=FILE, group="QCD-jet", regions=["gcr"], procs=[r"G + Jets"],
        label=r"$\gamma$ CR", sublabel=r"QCD-jet (class 1,3), $\gamma$+jets",
        ylabel="Top-tag mis-tag rate", ylim=(0, 0.08), sflim=(0, 4),
        zoomlim=(0, 0.08), zoomsflim=(0, 4), legloc="upper right"),
}

JOBS_TOP = {
    "effdata_topjet_ttcr_combined": dict(
        file=FILE, group="Top-jet", regions=["tmcr", "tecr"], procs=None,
        label=r"$t\bar{t}$ CR ($\mu$+e)", sublabel="Top-jet (class 5,6,7,8), all bkg",
        ylabel="Top-tag efficiency", ylim=(-5, 25), sflim=(-6, 27),
        zoomlim=(0, 1.5), zoomsflim=(0, 3), legloc="upper left"),
    "effdata_topjet_tmcr": dict(
        file=FILE, group="Top-jet", regions=["tmcr"], procs=None,
        label=r"$t\bar{t}$($\mu$) CR", sublabel="Top-jet (class 5,6,7,8), all bkg",
        ylabel="Top-tag efficiency", ylim=(-5, 25), sflim=(-6, 27),
        zoomlim=(0, 1.5), zoomsflim=(0, 3), legloc="upper left"),
    "effdata_topjet_tecr": dict(
        file=FILE, group="Top-jet", regions=["tecr"], procs=None,
        label=r"$t\bar{t}$(e) CR", sublabel="Top-jet (class 5,6,7,8), all bkg",
        ylabel="Top-tag efficiency", ylim=(-5, 25), sflim=(-6, 27),
        zoomlim=(0, 1.5), zoomsflim=(0, 3), legloc="upper left"),
}

sf_qcd = None
results = {}  # key -> SF-corrected measurement, for the correctionlib file
for key, cfg in JOBS_QCD.items():
    print(f"{key} ...")
    res = measure(cfg["file"], cfg["group"], cfg["regions"], cfg["procs"])
    report(cfg, "B-jet subtracted", res)
    series = [dict(res=res, label="Data (B-jet subtracted)",
                   color=DATA_COLOR, marker="s")]
    make_plot(key, cfg, series)
    if ARGS.mc_eff:
        make_eff_only_plot(key, cfg, res, "mc")
    if ARGS.data_eff:
        make_eff_only_plot(key, cfg, res, "data", series[0]["label"], "s")
    sf_qcd = (res["sf"], res["sf_err"])  # feeds the Top-jet subtraction
    results[key] = res

# Two subtraction variants for the Top-jet class:
#   raw       - README's boxed eff_Data formula, uncorrected MC B/QCD yields
#   corrected - the same with SF_B and the SF_QCD just derived applied to the
#               subtraction terms, per toptag/docs/scalefactors.md
# Both are printed; only the SF-corrected one is plotted.
for key, cfg in JOBS_TOP.items():
    print(f"{key} ...")
    res_raw = measure(cfg["file"], cfg["group"], cfg["regions"], cfg["procs"])
    res_cor = measure(cfg["file"], cfg["group"], cfg["regions"], cfg["procs"],
                      sf_qcd=sf_qcd)
    report(cfg, "raw subtraction", res_raw)
    report(cfg, "SF-corrected subtraction", res_cor)
    results[key] = res_cor
    series = [dict(res=res_cor, label=r"Data ($SF$-corrected subtraction)",
                   color=DATA_COLOR, marker="^")]
    make_plot(key, cfg, series, zoom=False)
    make_plot(key, cfg, series, zoom=True)
    if ARGS.mc_eff:
        make_eff_only_plot(key, cfg, res_cor, "mc")
    if ARGS.data_eff:
        for zoom in (False, True):
            make_eff_only_plot(key, cfg, res_cor, "data", series[0]["label"],
                               "^", zoom=zoom)



# ─── correctionlib output ────────────────────────────────────────────────────

def b_jet_eff_mc(res):
    """MC B-jet tag efficiency from the B-jet subtraction terms of a measurement."""
    b = res["terms"]["B-jet"]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(b["tot"] > 0, b["tag"] / b["tot"], np.nan)


def write_sf_json(qcd, top, path):
    """correctionlib file: SF_QCD (gamma CR), SF_Top (ttbar e+mu, SF-corrected),
    SF_B = 1 +- 0.5, and the MC efficiencies of the same measurements."""
    nbin = len(top["sf"])
    sfs = {"top": sf_json.sanitize(top["sf"], top["sf_err"], "SF_Top"),
           "qcd": sf_json.sanitize(qcd["sf"], qcd["sf_err"], "SF_QCD"),
           "b": (np.full(nbin, SF_B_CENTRAL), np.full(nbin, SF_B_UNCERTAINTY))}
    effs = {"top": sf_json.sanitize(top["eff_mc"], top["eff_mc_err"], "eff_MC Top", (0.0, 0.0))[0],
            "qcd": sf_json.sanitize(qcd["eff_mc"], qcd["eff_mc_err"], "eff_MC QCD", (0.0, 0.0))[0],
            "b": sf_json.sanitize(b_jet_eff_mc(top), np.zeros(nbin), "eff_MC B", (0.0, 0.0))[0]}
    cset = sf_json.build_correction_set(ERA, top["edges"], sfs, effs, WP)
    sf_json.write(cset, path)


write_sf_json(results["effdata_qcdjet_gcr"], results["effdata_topjet_ttcr_combined"],
              ARGS.json or os.path.join(CFG.outbase, "sf_json", f"toptag_sf_{ERA}.json"))
print(f"Done. Plots saved in {OUTDIR}/")
