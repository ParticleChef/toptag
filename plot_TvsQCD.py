#!/usr/bin/env python3
"""
plot_TvsQCD.py

Normalized TvsQCD shape comparison for background MC (Figure 21 style, AN2020_061).
Uses AK15 ParticleTransformer (ParT) scores from Run 3 private NanoAOD.

Run from toptag/ directory after: source ~/envs/p38/bin/activate

Quick test (1 file per dataset chunk, ~10 min with 40 workers):
    python plot_TvsQCD.py --year 2022pre  --workers 40 --max-files 1 --output plots/
    python plot_TvsQCD.py --year 2022post --workers 40 --max-files 1 --output plots/

Full run (all files, save histograms to reuse later):
    python plot_TvsQCD.py --year 2022pre  --workers 40 --save-hists hists/TvsQCD_2022pre.coffea  --output plots/
    python plot_TvsQCD.py --year 2022post --workers 40 --save-hists hists/TvsQCD_2022post.coffea --output plots/

Replot from saved histograms (instant):
    python plot_TvsQCD.py --year 2022pre  --load-hists hists/TvsQCD_2022pre.coffea  --output plots/
    python plot_TvsQCD.py --year 2022post --load-hists hists/TvsQCD_2022post.coffea --output plots/
"""

import argparse
import sys
import os
import gzip
import json
import re
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import awkward as ak
import hist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from concurrent.futures import ThreadPoolExecutor
import uproot

from coffea import processor
from coffea.processor import dict_accumulator, value_accumulator
from coffea.util import save, load

# Custom AK15 NanoAOD schema from the main analysis framework
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'decaf/analysis'))
from libs.mycoffeav2 import CustomNanoAODSchema


# ============================================================
# CONFIG
# ============================================================

YEAR_CONFIG = {
    '2022pre':  {'meta': '2022_private_v1',   'lumi': 7.980},   # fb^-1
    '2022post': {'meta': '2022EE_private_v1', 'lumi': 26.67},   # fb^-1
}

# Process groups: (group_name, [substrings that match the dataset name]).
# First match wins, so put more specific patterns before broader ones.
PROCESS_KEYS = [
    ('tt',        ['TTto']),
    ('wjets',     ['WtoLNu-2Jets']),
    ('zjets',     ['Zto2Nu-2Jets']),
    ('dyjets',    ['DYto2L-2Jets']),
    ('gjets',     ['GJ_PTG']),
    ('singletop', ['ST-T']),
    ('diboson',   ['WW_TuneCP5', 'WZ_TuneCP5', 'ZZ_TuneCP5']),
    ('qcd',       ['QCD_PT']),
]

# Plot colors and labels for the legend (Fig 21 style)
PLOT_STYLE = {
    'tt':        {'color': '#d62728', 'label': r'$t\bar{t}$',           'zorder': 8},
    'wjets':     {'color': '#1f77b4', 'label': 'W+Jets',                'zorder': 5},
    'zjets':     {'color': '#2ca02c', 'label': r'$Z(\nu\nu)$+Jets',     'zorder': 6},
    'dyjets':    {'color': '#ff7f0e', 'label': 'DY+Jets',               'zorder': 4},
    'gjets':     {'color': '#e377c2', 'label': r'$\gamma$+Jets',        'zorder': 7},
    'singletop': {'color': '#8c564b', 'label': 'Single Top',            'zorder': 9},
    'diboson':   {'color': '#17becf', 'label': 'Diboson',               'zorder': 3},
    'qcd':       {'color': '#7f7f7f', 'label': 'QCD',                   'zorder': 2},
}

# AK15 jet and event preselection (signal region)
AK15_PT_GOOD   = 160.   # GeV — isGoodAK15 threshold (ids.py:309)
AK15_ETA_MAX   = 2.4
AK15_PT_LEAD   = 250.   # GeV — leading jet required
PUPPIMET_MIN   = 200.   # GeV — simplified MET cut (full analysis uses recoil > 250)
MU_PT_MIN      = 20.    # GeV — tight muon veto threshold
ELE_PT_MIN     = 40.    # GeV — tight electron veto threshold

# TvsQCD histogram binning and working point
TVSQCD_NBINS = 50
TVSQCD_WP    = 0.26


# ============================================================
# Helper: dataset name → process group
# ============================================================

def dataset_to_process(name):
    """Strip the ____N_ chunk suffix and match to a process group."""
    base = re.sub(r'____\d+_$', '', name)
    for proc, keys in PROCESS_KEYS:
        if any(k in base for k in keys):
            return proc
    return None


# ============================================================
# File validation — filters corrupt/empty files before coffea preprocessing
# ============================================================

def _check_one_file(fpath):
    """Return fpath if the file is readable and has Events entries, else None."""
    try:
        with uproot.open(fpath) as f:
            if 'Events' not in f:
                return None
            if f['Events'].num_entries == 0:
                return None
        return fpath
    except Exception:
        return None


def validate_fileset(fileset, workers=40):
    """
    Open every file with uproot to check it is not corrupt or empty.
    Uses a thread pool for speed — typically adds <2 min for ~3000 files.
    Returns (good_fileset, n_bad).
    """
    all_pairs = [(ds, f) for ds, files in fileset.items() for f in files]
    n_total   = len(all_pairs)
    print(f'  Validating {n_total} files with {workers} threads...')

    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(lambda p: (p[0], _check_one_file(p[1])), all_pairs))

    good_fileset = {}
    n_bad = 0
    for ds, fpath in results:
        if fpath is not None:
            good_fileset.setdefault(ds, []).append(fpath)
        else:
            n_bad += 1

    # Remove dataset entries that lost all files
    good_fileset = {ds: files for ds, files in good_fileset.items() if files}
    print(f'  {n_total - n_bad} good, {n_bad} bad/empty files removed')
    return good_fileset, n_bad


# ============================================================
# Coffea processor
# ============================================================

class TvsQCDProcessor(processor.ProcessorABC):
    """
    Fills one TvsQCD histogram and sumw per dataset chunk.
    Histograms are weighted by genWeight only; xs scaling is
    applied in post-processing.
    """

    def process(self, events):
        dataset = events.metadata['dataset']

        # Skip data (no genWeight) and any samples not in our process groups
        if not hasattr(events, 'genWeight') or dataset_to_process(dataset) is None:
            return {
                'sumw':  dict_accumulator({}),
                'hists': dict_accumulator({}),
            }

        gw = ak.to_numpy(events.genWeight)

        # --- AK15 jet selection (isGoodAK15 from decaf/analysis/utils/ids.py) ---
        fj = events.AK15PuppiJet
        fj_good = fj[(fj.pt > AK15_PT_GOOD) & (np.abs(fj.eta) < AK15_ETA_MAX) & ((fj.jetId & 6) == 6)]
        leading_fj = ak.firsts(fj_good)

        has_leading_fj = (
            (ak.num(fj_good) > 0) &
            (ak.fill_none(leading_fj.pt, 0.) > AK15_PT_LEAD)
        )

        # --- MET cut ---
        met_pass = events.PuppiMET.pt > PUPPIMET_MIN

        # --- Zero lepton veto ---
        mu  = events.Muon
        ele = events.Electron
        n_tight_mu  = ak.sum((mu.pt > MU_PT_MIN)  & (np.abs(mu.eta)  < 2.4) & mu.tightId,        axis=1)
        n_tight_ele = ak.sum((ele.pt > ELE_PT_MIN) & (np.abs(ele.eta) < 2.5) & (ele.cutBased >= 4), axis=1)
        zero_lepton = (n_tight_mu == 0) & (n_tight_ele == 0)

        sel = has_leading_fj & met_pass & zero_lepton

        # --- TvsQCD from ParT scores (hadmonotop_run3.py lines 559-561) ---
        lf      = leading_fj
        probQCD = (lf.ParT_probQCDbb + lf.ParT_probQCDcc +
                   lf.ParT_probQCDb  + lf.ParT_probQCDc  + lf.ParT_probQCDothers)
        probT   = lf.ParT_probTopbWqq + lf.ParT_probTopbWcs
        TvsQCD  = ak.fill_none(probT / (probT + probQCD), -1.)

        tvsqcd_vals = ak.to_numpy(TvsQCD[sel])
        gw_sel      = gw[ak.to_numpy(sel)]

        h = hist.Hist(
            hist.axis.Regular(TVSQCD_NBINS, 0, 1, name='TvsQCD',
                              label='ParT TvsQCD score (leading AK15 jet)'),
            storage=hist.storage.Weight(),
        )
        h.fill(TvsQCD=tvsqcd_vals, weight=gw_sel)

        return {
            'sumw':  dict_accumulator({dataset: value_accumulator(float, float(np.sum(np.abs(gw))))}),
            'hists': dict_accumulator({dataset: h}),
        }

    def postprocess(self, accumulator):
        return accumulator


# ============================================================
# File loading
# ============================================================

def load_fileset(meta_name, max_files):
    """Build coffea fileset and xsec map from metadata, MC only."""
    meta_path = os.path.join('decaf/analysis/metadata', meta_name + '.json.gz')
    with gzip.open(meta_path) as f:
        meta = json.load(f)

    fileset = {}
    xsec_map = {}
    skipped  = 0
    for ds, info in meta.items():
        if dataset_to_process(ds) is None:
            skipped += 1
            continue
        files = info['files'][:max_files] if max_files > 0 else info['files']
        fileset[ds]  = files
        xsec_map[ds] = info['xs']

    n_files = sum(len(v) for v in fileset.values())
    print(f'  {len(fileset)} dataset chunks, {n_files} total files  ({skipped} skipped)')
    return fileset, xsec_map


# ============================================================
# Post-processing: xs scaling and process grouping
# ============================================================

def build_process_hists(out, xsec_map, lumi_pb):
    """
    Scale each dataset histogram by xs * lumi / sumw, then sum by process group.
    Returns dict: process_name → scaled hist.Hist
    """
    proc_hists = {}

    for ds, h in out['hists'].items():
        proc = dataset_to_process(ds)
        if proc is None:
            continue
        sumw = out['sumw'][ds].value
        if sumw == 0.:
            continue
        xs    = xsec_map.get(ds, 1.)
        scale = xs * lumi_pb / sumw

        scaled = h * scale
        if proc in proc_hists:
            proc_hists[proc] = proc_hists[proc] + scaled
        else:
            proc_hists[proc] = scaled

    return proc_hists


# ============================================================
# Plotting
# ============================================================

def make_plot(proc_hists, year, outdir):
    os.makedirs(outdir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 6))

    lumi_str = {'2022pre': r'7.98 fb$^{-1}$', '2022post': r'26.67 fb$^{-1}$'}
    era_str  = {'2022pre': '2022 (pre-EE, Run CD)',  '2022post': '2022 (post-EE, Run EFG)'}

    plotted = 0
    for proc, style in PLOT_STYLE.items():
        if proc not in proc_hists:
            continue
        h      = proc_hists[proc]
        vals   = h.values()
        edges  = h.axes[0].edges
        widths = edges[1:] - edges[:-1]
        total  = np.sum(vals * widths)
        if total <= 0.:
            continue

        norm_vals = vals / total
        # step plot: repeat last value for the final bin edge
        ax.step(edges, np.append(norm_vals, norm_vals[-1]),
                where='post',
                color=style['color'],
                label=style['label'],
                linewidth=1.6,
                zorder=style.get('zorder', 1))
        plotted += 1

    ax.axvline(TVSQCD_WP, color='black', linestyle='--', linewidth=1.2,
               alpha=0.8, label=f'WP = {TVSQCD_WP}', zorder=10)

    ax.set_yscale('log')
    ax.set_xlim(0, 1)
    ax.set_xlabel('particleNet TvsQCD (leading AK15 jet)', fontsize=13)
    ax.set_ylabel('Normalized Events / bin width', fontsize=13)

    # CMS style header
    ax.text(0.01, 1.025, 'CMS', transform=ax.transAxes,
            fontsize=14, fontweight='bold', va='bottom')
    ax.text(0.10, 1.025, 'private work', transform=ax.transAxes,
            fontsize=12, style='italic', va='bottom')
    ax.text(0.99, 1.025, f'{lumi_str.get(year, "")} (13.6 TeV)',
            transform=ax.transAxes, fontsize=11, ha='right', va='bottom')

    ax.set_title(f'Signal region preselection — {era_str.get(year, year)}',
                 fontsize=10, pad=22)

    ax.legend(fontsize=9, ncol=2, loc='upper left',
              framealpha=0.85, handlelength=1.5)

    for ext in ('pdf', 'png'):
        outfile = os.path.join(outdir, f'TvsQCD_shape_{year}.{ext}')
        fig.savefig(outfile, bbox_inches='tight', dpi=150)
        print(f'Saved: {outfile}')
    plt.close(fig)


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Make normalized TvsQCD shape plot (Fig 21 style, AN2020_061)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--year',       required=True, choices=['2022pre', '2022post'],
                        help='Era: 2022pre (2022_private_v1) or 2022post (2022EE_private_v1)')
    parser.add_argument('--workers',    type=int, default=40,
                        help='Parallel coffea workers (default: 40). Reduce if server is busy.')
    parser.add_argument('--max-files',  type=int, default=0,
                        help='Max ROOT files per dataset chunk (0 = all; use 1 for quick test)')
    parser.add_argument('--output',     default='plots/',
                        help='Directory for output plots (default: plots/)')
    parser.add_argument('--save-hists', default='',
                        help='Save raw histograms to this .coffea file for later reuse')
    parser.add_argument('--load-hists', default='',
                        help='Skip running and load histograms from this .coffea file')
    args = parser.parse_args()

    cfg      = YEAR_CONFIG[args.year]
    lumi     = cfg['lumi']          # fb^-1
    lumi_pb  = lumi * 1000.         # pb^-1

    if args.load_hists:
        print(f'Loading histograms from {args.load_hists}')
        bundle = load(args.load_hists)
        out      = bundle['out']
        xsec_map = bundle['xsec_map']
    else:
        print(f'Year: {args.year}  |  Lumi: {lumi} fb^-1  |  Workers: {args.workers}')
        if args.max_files > 0:
            print(f'  [quick mode] max {args.max_files} file(s) per dataset chunk')
        print(f'Loading metadata: {cfg["meta"]}')
        fileset, xsec_map = load_fileset(cfg['meta'], args.max_files)

        fileset, n_bad = validate_fileset(fileset, workers=args.workers)
        if not fileset:
            print('ERROR: no valid files found. Exiting.')
            sys.exit(1)

        print('Running coffea processor...')
        out = processor.run_uproot_job(
            fileset,
            treename           = 'Events',
            processor_instance = TvsQCDProcessor(),
            executor           = processor.futures_executor,
            executor_args      = {
                'schema':       CustomNanoAODSchema,
                'workers':      args.workers,
                'skipbadfiles': True,
            },
        )

        if args.save_hists:
            os.makedirs(os.path.dirname(args.save_hists) or '.', exist_ok=True)
            save({'out': out, 'xsec_map': xsec_map}, args.save_hists)
            print(f'Histograms saved to {args.save_hists}')

    print('Building process histograms...')
    proc_hists = build_process_hists(out, xsec_map, lumi_pb)
    for proc, h in proc_hists.items():
        n = np.sum(h.values())
        print(f'  {proc:12s}: {n:.2e} weighted events after preselection')

    print(f'Making plot → {args.output}')
    make_plot(proc_hists, args.year, args.output)
    print('Done.')


if __name__ == '__main__':
    main()
