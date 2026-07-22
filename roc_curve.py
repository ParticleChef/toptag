#!/usr/bin/env python3
"""
roc_curve.py

ROC curve for the ParT TvsQCD discriminator on AK15 jets.
  Signal     : Mphi-1000_Mchi-150 (mono-top)
  Background : Zto2Nu (Z → νν)

Jet preselection (matches decaf ids.py isGoodAK15):
  - pT > 160 GeV, |η| < 2.4, tight+tightLepVeto jetId
  - Leading jet required with pT > 250 GeV

Usage
-----
Quick test (1 file per chunk, ~few minutes):
    python roc_curve.py \\
        --metadata decaf/analysis/metadata/2022_private_v2.json.gz \\
        --year 2022pre --max-files 1 --workers 10 --output plots/

Full run + save arrays for instant replotting:
    python roc_curve.py \\
        --metadata decaf/analysis/metadata/2022_private_v2.json.gz \\
        --year 2022pre --workers 40 \\
        --save-arrays hists/roc_arrays_2022pre.coffea --output plots/

Replot from saved arrays (no reprocessing needed):
    python roc_curve.py --load-arrays hists/roc_arrays_2022pre.coffea --output plots/
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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mplhep

from concurrent.futures import ThreadPoolExecutor
import uproot

from coffea import processor
from coffea.processor import column_accumulator
from coffea.util import save, load

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'decaf/analysis'))
from libs.mycoffeav2 import CustomNanoAODSchema


# ============================================================
# CONFIG
# ============================================================

SIGNAL_PATTERN = 'Mphi-1000_Mchi-150'
BKG_PATTERN    = 'Zto2Nu'

# AK15 jet preselection (matches decaf/analysis/utils/ids.py isGoodAK15)
AK15_PT_GOOD = 160.
AK15_ETA_MAX = 2.4
AK15_PT_LEAD = 250.

# ROC histogram resolution
TVSQCD_NBINS = 1000

# Working points defined by target background efficiency
BKG_EFF_WPS = [0.05, 0.025, 0.01, 0.005, 0.001]
WP_COLORS   = ['#e41a1c', '#ff7f00', '#4daf4a', '#984ea3', '#377eb8']

# Signal-region event-level cuts
MET_PT_MIN   = 350.
BJET_DR_CONE = 1.5

# PNet b-tag loose WPs (from decaf/analysis/utils/common.py PNetUParTWPs)
PNET_BTAG_LOOSE = {
    '2022pre':  0.047,
    '2022post': 0.0499,
    '2023pre':  0.0358,
    '2023post': 0.0359,
}

LUMI_MAP = {
    '2022pre':  (7.98,  '2022 (pre-EE)'),
    '2022post': (26.67, '2022 (post-EE)'),
    '2023pre':  (17.96, '2023 (pre-BPix)'),
    '2023post': (9.68,  '2023 (post-BPix)'),
}


# ============================================================
# Helpers
# ============================================================

def classify_dataset(name):
    """Return 'signal', 'background', or None."""
    base = re.sub(r'____\d+_$', '', name)
    if SIGNAL_PATTERN in base:
        return 'signal'
    if BKG_PATTERN in base:
        return 'background'
    return None


def _check_file(fpath):
    try:
        with uproot.open(fpath) as f:
            if 'Events' not in f or f['Events'].num_entries == 0:
                return None
        return fpath
    except Exception:
        return None


def validate_fileset(fileset, workers):
    pairs   = [(ds, f) for ds, files in fileset.items() for f in files]
    n_total = len(pairs)
    print(f'  Validating {n_total} files with {workers} threads ...')
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(lambda p: (p[0], _check_file(p[1])), pairs))
    good  = {}
    n_bad = 0
    for ds, fpath in results:
        if fpath is not None:
            good.setdefault(ds, []).append(fpath)
        else:
            n_bad += 1
    good = {ds: fs for ds, fs in good.items() if fs}
    print(f'  {n_total - n_bad} good, {n_bad} removed')
    return good, n_bad


# ============================================================
# Coffea processor
# ============================================================

_EMPTY = np.array([], dtype=np.float32)


class ROCProcessor(processor.ProcessorABC):
    """
    Collects TvsQCD scores and genWeights for signal and background.
    Uses column_accumulator (native coffea type) so arrays are
    concatenated correctly across workers.
    """

    def __init__(self, year):
        self.year = year

    def process(self, events):
        dataset = events.metadata['dataset']
        label   = classify_dataset(dataset)

        out = {
            'sig_scores':  column_accumulator(_EMPTY.copy()),
            'sig_weights': column_accumulator(_EMPTY.copy()),
            'bkg_scores':  column_accumulator(_EMPTY.copy()),
            'bkg_weights': column_accumulator(_EMPTY.copy()),
        }

        if label is None or not hasattr(events, 'genWeight'):
            return out

        gw = ak.to_numpy(events.genWeight)

        # Good AK15 jets: tight + tightLepVeto (jetId bits 1,2 → mask = 6)
        fj      = events.AK15PuppiJet
        fj_good = fj[
            (fj.pt > AK15_PT_GOOD) &
            (np.abs(fj.eta) < AK15_ETA_MAX) &
            ((fj.jetId & 6) == 6)
        ]
        leading = ak.firsts(fj_good)

        # MET cut
        met_cut = events.MET.pt > MET_PT_MIN

        # Loose electron veto (ids.py isLooseElectron, 2022/2023)
        e      = events.Electron
        e_aeta = np.abs(e.eta + e.deltaEtaSC)
        loose_e = (
            ((e.pt > 10) & (e_aeta < 1.4442) & (e.cutBased >= 1)) |
            ((e.pt > 10) & (e_aeta > 1.5660) & (e_aeta < 2.5) & (e.cutBased >= 1))
        )

        # Loose muon veto (ids.py isLooseMuon)
        mu = events.Muon
        loose_mu = (
            (mu.pt > 10) & (np.abs(mu.eta) < 2.4) &
            mu.looseId & mu.isTracker & mu.isPFcand & mu.isGlobal & (mu.pfIsoId >= 2)
        )

        # Loose photon veto (ids.py isLoosePhoton, 2022/2023)
        pho      = events.Photon
        pho_aeta = np.abs(pho.eta)
        loose_pho = (
            (pho.pt > 20) &
            (~(pho_aeta > 1.4442) | (pho_aeta > 1.5660)) &
            (pho_aeta < 2.5) &
            (pho.cutBased >= 1) &
            pho.electronVeto
        )

        lepton_veto = (
            (ak.sum(loose_e,   axis=1) == 0) &
            (ak.sum(loose_mu,  axis=1) == 0) &
            (ak.sum(loose_pho, axis=1) == 0)
        )

        # b-jet veto: no PNet-loose b-tagged AK4 jet outside dR < 1.5 of leading AK15
        ak4      = events.Jet
        good_ak4 = ak4[(ak4.pt > 30) & (np.abs(ak4.eta) < 2.4) & ((ak4.jetId & 6) == 6)]
        bjets    = good_ak4[good_ak4.btagPNetB > PNET_BTAG_LOOSE[self.year]]
        lead_eta = ak.fill_none(leading.eta, 999.)
        lead_phi = ak.fill_none(leading.phi, 999.)
        deta     = bjets.eta - lead_eta
        dphi     = bjets.phi - lead_phi
        dphi     = ak.where(dphi >  np.pi, dphi - 2*np.pi, dphi)
        dphi     = ak.where(dphi < -np.pi, dphi + 2*np.pi, dphi)
        dr       = np.sqrt(deta**2 + dphi**2)
        bjet_veto = ak.sum(dr > BJET_DR_CONE, axis=1) == 0

        sel = (
            (ak.num(fj_good) > 0) &
            (ak.fill_none(leading.pt, 0.) > AK15_PT_LEAD) &
            met_cut &
            lepton_veto &
            bjet_veto
        )

        # TvsQCD discriminant (toptagging_guide.md)
        lf      = leading
        probQCD = (lf.ParT_probQCDbb + lf.ParT_probQCDcc +
                   lf.ParT_probQCDb  + lf.ParT_probQCDc  + lf.ParT_probQCDothers)
        probT   = lf.ParT_probTopbWqq + lf.ParT_probTopbWcs
        TvsQCD  = ak.fill_none(probT / (probT + probQCD), -1.)

        scores  = ak.to_numpy(TvsQCD[sel]).astype(np.float32)
        # Use sign(genWeight) so each event counts as +1 or -1.
        # Raw genWeights vary by orders of magnitude across Zto2Nu pT-slices;
        # sign-weighting treats all events equally while still cancelling NLO
        # negative-weight events within each sample.
        weights = np.sign(gw[ak.to_numpy(sel)]).astype(np.float32)

        if label == 'signal':
            out['sig_scores']  = column_accumulator(scores)
            out['sig_weights'] = column_accumulator(weights)
        else:
            out['bkg_scores']  = column_accumulator(scores)
            out['bkg_weights'] = column_accumulator(weights)

        return out

    def postprocess(self, accumulator):
        return accumulator


# ============================================================
# File loading
# ============================================================

def load_fileset(meta_path, max_files):
    with gzip.open(meta_path, 'rt') as f:
        meta = json.load(f)

    fileset = {}
    skipped = 0
    for ds, info in meta.items():
        if classify_dataset(ds) is None:
            skipped += 1
            continue
        files = info['files'][:max_files] if max_files > 0 else info['files']
        fileset[ds] = files

    sig_ds = [k for k in fileset if classify_dataset(k) == 'signal']
    bkg_ds = [k for k in fileset if classify_dataset(k) == 'background']
    print(f'  Signal:     {len(sig_ds)} chunks, '
          f'{sum(len(fileset[k]) for k in sig_ds)} files')
    print(f'  Background: {len(bkg_ds)} chunks, '
          f'{sum(len(fileset[k]) for k in bkg_ds)} files')
    print(f'  Skipped:    {skipped} unrelated datasets')
    return fileset


# ============================================================
# ROC computation
# ============================================================

def compute_roc(sig_scores, sig_weights, bkg_scores, bkg_weights):
    """
    Build weighted histograms then sweep threshold right→left to
    get signal efficiency and background efficiency curves.
    Returns (sig_eff, bkg_eff, thresholds, auc).

    NLO samples have ~20% negative genWeights.  The denominator must
    use the net sum (not abs), otherwise efficiency at threshold=0
    does not reach 1.0 and the AUC integral is wrong.
    """
    # Drop events with score < 0 (fill_none sentinel for invalid jets)
    sig_mask = sig_scores >= 0.
    bkg_mask = bkg_scores >= 0.
    sig_scores, sig_weights = sig_scores[sig_mask], sig_weights[sig_mask]
    bkg_scores, bkg_weights = bkg_scores[bkg_mask], bkg_weights[bkg_mask]

    edges = np.linspace(0., 1., TVSQCD_NBINS + 1)
    sig_hist, _ = np.histogram(sig_scores, bins=edges, weights=sig_weights)
    bkg_hist, _ = np.histogram(bkg_scores, bins=edges, weights=bkg_weights)

    # Net totals: positive for physical NLO samples even with negative weights
    sig_total = np.sum(sig_weights)
    bkg_total = np.sum(bkg_weights)

    # Cumulative sum from high score to low → events passing threshold t
    sig_pass = np.cumsum(sig_hist[::-1])[::-1]
    bkg_pass = np.cumsum(bkg_hist[::-1])[::-1]

    sig_eff = np.concatenate([[0.], sig_pass / sig_total, [1.]])
    bkg_eff = np.concatenate([[0.], bkg_pass / bkg_total, [1.]])
    thrs    = np.concatenate([[1.], edges[:-1][::-1], [0.]])

    auc = float(np.trapz(sig_eff, bkg_eff))
    return sig_eff, bkg_eff, thrs, auc


def find_wp_at_bkg_eff(bkg_eff, sig_eff, thrs, target):
    """Return (threshold, sig_eff, actual_bkg_eff) at the point closest to target bkg_eff."""
    idx = np.argmin(np.abs(bkg_eff - target))
    return thrs[idx], sig_eff[idx], bkg_eff[idx]


# ============================================================
# Plotting
# ============================================================

def make_roc_plot(out, outdir, year):
    sig_scores  = out['sig_scores'].value
    sig_weights = out['sig_weights'].value
    bkg_scores  = out['bkg_scores'].value
    bkg_weights = out['bkg_weights'].value

    print(f'  Signal events passing presel:     {len(sig_scores):,}')
    print(f'  Background events passing presel: {len(bkg_scores):,}')

    if len(sig_scores) == 0 or len(bkg_scores) == 0:
        print('ERROR: empty signal or background — cannot plot ROC')
        return

    sig_eff, bkg_eff, thrs, auc = compute_roc(
        sig_scores, sig_weights, bkg_scores, bkg_weights)

    # Find and print working points
    wps = []
    print(f'\n  AUC = {auc:.4f}')
    print(f'\n  {"Target ε_bkg":>12}  {"Actual ε_bkg":>12}  {"Threshold":>9}  {"ε_sig":>7}  {"Rejection":>10}')
    print(f'  {"-"*12}  {"-"*12}  {"-"*9}  {"-"*7}  {"-"*10}')
    for target in BKG_EFF_WPS:
        thr, se, be = find_wp_at_bkg_eff(bkg_eff, sig_eff, thrs, target)
        rej_str = f'{1/be:.1f}x' if be > 0 else 'inf'
        print(f'  {target*100:>11.1f}%  {be*100:>11.2f}%  {thr:>9.4f}  {se:>7.3f}  {rej_str:>10}')
        wps.append((target, thr, se, be))

    lumi_val, era_label = LUMI_MAP.get(year, (0., year))
    os.makedirs(outdir, exist_ok=True)

    # ── Plot 1: x = background rejection (log scale) ──────────────────────
    rej_curve = np.full_like(bkg_eff, np.nan)
    valid = bkg_eff > 0
    rej_curve[valid] = 1. / bkg_eff[valid]

    mplhep.style.use('CMS')
    fig1, ax1 = plt.subplots(figsize=(8, 7))
    ax1.plot(rej_curve[valid], sig_eff[valid], color='steelblue', linewidth=2,
             label=f'ParT TvsQCD  (AUC = {auc:.3f})')
    for (target, thr, se, be), col in zip(wps, WP_COLORS):
        if be > 0:
            ax1.scatter([1./be], [se], color=col, zorder=5, s=60,
                        label=f'ε_bkg={target*100:.1f}%: ε_sig={se:.3f}')
    ax1.set_xlabel(r'Background rejection  ($Z\to\nu\nu$)', fontsize=14)
    ax1.set_ylabel(r'Signal efficiency  ($m_\phi\!=\!1000,\;m_\chi\!=\!150$)', fontsize=14)
    ax1.set_xscale('log')
    ax1.set_ylim(0., 1.)
    ax1.legend(fontsize=9, loc='upper right', framealpha=0.9)
    mplhep.cms.label(ax=ax1, label='Private Work',
                     data=False, lumi=lumi_val, com=13.6, fontsize=13)
    ax1.set_title(f'ParT AK15 TvsQCD ROC — {era_label}', fontsize=11, pad=28)
    for ext in ('pdf', 'png'):
        path = os.path.join(outdir, f'roc_TvsQCD_{year}_rej.{ext}')
        fig1.savefig(path, bbox_inches='tight', dpi=300)
        print(f'Saved: {path}')
    plt.close(fig1)

    # ── Plot 2: x = background efficiency ─────────────────────────────────
    mplhep.style.use('CMS')
    fig2, ax2 = plt.subplots(figsize=(8, 7))
    ax2.plot(bkg_eff, sig_eff, color='steelblue', linewidth=2,
             label=f'ParT TvsQCD  (AUC = {auc:.3f})')
    ax2.plot([0, 1], [0, 1], 'k--', linewidth=0.8, alpha=0.5, label='Random')
    for (target, thr, se, be), col in zip(wps, WP_COLORS):
        ax2.scatter([be], [se], color=col, zorder=5, s=60,
                    label=f'ε_bkg={target*100:.1f}%: ε_sig={se:.3f}')
    ax2.set_xlabel(r'Background efficiency  ($Z\to\nu\nu$)', fontsize=14)
    ax2.set_ylabel(r'Signal efficiency  ($m_\phi\!=\!1000,\;m_\chi\!=\!150$)', fontsize=14)
    ax2.set_xlim(0., 1.)
    ax2.set_ylim(0., 1.)
    ax2.legend(fontsize=9, loc='lower right', framealpha=0.9)
    mplhep.cms.label(ax=ax2, label='Private Work',
                     data=False, lumi=lumi_val, com=13.6, fontsize=13)
    ax2.set_title(f'ParT AK15 TvsQCD ROC — {era_label}', fontsize=11, pad=28)
    for ext in ('pdf', 'png'):
        path = os.path.join(outdir, f'roc_TvsQCD_{year}_eff.{ext}')
        fig2.savefig(path, bbox_inches='tight', dpi=300)
        print(f'Saved: {path}')
    plt.close(fig2)


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='ParT TvsQCD ROC curve: Mphi-1000_Mchi-150 vs Zto2Nu',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--metadata',     default='',
                        help='Path to metadata .json.gz (required unless --load-arrays)')
    parser.add_argument('--year',         default='2022pre',
                        choices=list(LUMI_MAP.keys()),
                        help='Era label for the plot header (default: 2022pre)')
    parser.add_argument('--workers',      type=int, default=40,
                        help='Coffea futures workers (default: 40)')
    parser.add_argument('--max-files',    type=int, default=0,
                        help='Max ROOT files per dataset chunk (0=all; 1=quick test)')
    parser.add_argument('--output',       default='plots/',
                        help='Output directory for plots (default: plots/)')
    parser.add_argument('--save-arrays',  default='',
                        help='Save score arrays to this .coffea file for later reuse')
    parser.add_argument('--load-arrays',  default='',
                        help='Skip processing; load arrays from this .coffea file')
    args = parser.parse_args()

    if args.load_arrays:
        print(f'Loading arrays from {args.load_arrays}')
        bundle = load(args.load_arrays)
        out  = bundle['out']
        year = bundle.get('year', args.year)
    else:
        if not args.metadata:
            parser.error('--metadata is required unless --load-arrays is given')

        print(f'Year: {args.year}  |  Workers: {args.workers}'
              + (f'  |  max-files: {args.max_files}' if args.max_files > 0 else ''))
        print(f'Metadata: {args.metadata}')
        fileset = load_fileset(args.metadata, args.max_files)

        fileset, _ = validate_fileset(fileset, args.workers)
        if not fileset:
            print('ERROR: no valid files found')
            sys.exit(1)

        print('Running coffea processor ...')
        out = processor.run_uproot_job(
            fileset,
            treename           = 'Events',
            processor_instance = ROCProcessor(year=args.year),
            executor           = processor.futures_executor,
            executor_args      = {
                'schema':       CustomNanoAODSchema,
                'workers':      args.workers,
                'skipbadfiles': True,
            },
        )
        year = args.year

        if args.save_arrays:
            os.makedirs(os.path.dirname(args.save_arrays) or '.', exist_ok=True)
            save({'out': out, 'year': year}, args.save_arrays)
            print(f'Arrays saved to {args.save_arrays}')

    print(f'Making ROC plot → {args.output}')
    make_roc_plot(out, args.output, year)
    print('Done.')


if __name__ == '__main__':
    main()
