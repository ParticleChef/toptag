######
## Run Recoil binning playground - hadmonotop2022EE_0702.scaled, TvsQCD pass/fail split
######
import os
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import hist
import lz4.frame
import cloudpickle
from cycler import cycler

# ─── Path resolution ────────────────────────────────────────────────────────
# This script lives in toptag/studies/, but it reads the .scaled files from,
# and writes its plot directory into, the analysis working directory two
# levels up. Everything below therefore uses paths relative to that root, so
# the script behaves identically no matter which directory it is run from.
WORK_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(WORK_DIR)

plt.style.use(hep.style.CMS)

# ─── Load ───────────────────────────────────────────────────────────────────
INPUT_FILE = "hadmonotop2022EE_0702.scaled"
with open(INPUT_FILE, "rb") as f:
    hists = cloudpickle.loads(lz4.frame.decompress(f.read()))

year = "2022post"
whathist_tag = "0702"
lumi = {"2022post": 26.68}

## Merge VV
for dataset in ['bkg']:
    print(f'Now {dataset} processing')
    for key in hists[dataset].keys():
        hists[dataset][key]['VV'] = hists[dataset][key]['WW'] + hists[dataset][key]['WZ'] + hists[dataset][key]['ZZ']
        for p in ['WW', 'WZ', 'ZZ']:
            del hists[dataset][key][p]

stack = {
    'sr':[r'Z ($\nu\nu$) + Jets',r'W ($\ell\nu$) + Jets', 'VV','TT','G + Jets','ST', 'Z ($\ell\ell$) + Jets', 'QCD Multijet' ],
    'tw':[r'W ($\ell\nu$) + Jets', 'TT','QCD Multijet', 'VV','ST', 'Z ($\ell\ell$) + Jets',r'Z ($\nu\nu$) + Jets', 'G + Jets'],
    'zr':['Z ($\ell\ell$) + Jets',r'W ($\ell\nu$) + Jets', 'TT','QCD Multijet', 'VV','ST',  r'Z ($\nu\nu$) + Jets', 'G + Jets'],
    'gr':['G + Jets',r'W ($\ell\nu$) + Jets', 'TT','QCD Multijet', 'VV', 'ST',  'Z ($\ell\ell$) + Jets',r'Z ($\nu\nu$) + Jets'],
}

cat_colors = {
    r'W ($\ell\nu$) + Jets': "#1f77b4",
    'TT': "#ff7f0e",
    'QCD Multijet': "#2ca02c",
    'VV': "#9467bd",
    'ST': "#aec7e8",
    r'Z ($\ell\ell$) + Jets': "#ffbb78",
    r'Z ($\nu\nu$) + Jets': "#98df8a",
    'G + Jets': "#ff9896"
}

leg_reg = {
    'sr': 'SR', 'wmcr': r'W($\mu$) CR', 'wecr': r'W(e) CR', 'tmcr': r'Top($\mu$) CR', 'gcr': 'G CR',
    'tecr': r'Top(e) CR', 'zecr': r'Z(ee) CR', 'zmcr': r'Z($\mu\mu$) CR'
}

# region -> stack_map, per the notebook's per-region cells (sr already run/blinded;
# wmcr/tmcr/wecr/tecr share the 'tw' order, gcr uses 'gr', zmcr/zecr use 'zr')
REGION_STACK_MAP = {
    'sr':   'sr',
    'wmcr': 'tw',
    'tmcr': 'tw',
    'wecr': 'tw',
    'tecr': 'tw',
    'gcr':  'gr',
    'zmcr': 'zr',
    'zecr': 'zr',
}

# TvsQCD axis has edges [0, 0.33, 1] -> bin index 0 = fail, bin index 1 = pass
TAG_CATEGORIES = {
    'toptag_fail': 0,
    'toptag_pass': 1,
}
tag_label = {'toptag_fail': 'fail', 'toptag_pass': 'pass'}

new_bin = [250, 400, 600, 1000]
bin_tag = "-".join(str(b) for b in new_bin)  # encodes new_bin into output filenames
keys = ['ut'
       ]
# skipped automatically if not found below (mirrors temp_plotting_2023.py in case
# this file also ships without signal samples).
sig_key1 = 'sig_Mphi-2000_Mchi-150'
sig_key2 = 'sig_Mphi-200_Mchi-150'

for tag_name, tag_idx in TAG_CATEGORIES.items():
    outdir = f"plots{year}_{whathist_tag}_{tag_name}"
    os.makedirs(outdir, exist_ok=True)

    for region, stack_map in REGION_STACK_MAP.items():
        for key in keys:
            fig, (ax, rax) = plt.subplots(
                nrows=2,
                ncols=1,
                figsize=(10,10),
                gridspec_kw={"height_ratios": (3, 1)},
                sharex=True
            )

            fig.subplots_adjust(hspace=.07)
            hep.cms.label(ax=ax, llabel='Private Work', rlabel=str(lumi[year])+' fb$^{-1}$ (13.6 TeV)')
            #ax.set_prop_cycle(cycler(color=colors[stack_map]))

            bkg  = hists['bkg']
            data = hists['data']
            sig  = hists['sig']

            mcstack = None

            for sample in stack[stack_map]:
                try:
                    target = bkg[key][sample][{'region': region,'TvsQCD': tag_idx,'ut': hist.rebin(edges=new_bin)}]
                    bins = target.axes.edges[0]

                except:
                    print(f'In {region} sample of {sample} is Not exist')
                    continue

                try:
                    mcstack += bkg[key][sample][{'region': region,'TvsQCD':tag_idx,'ut': hist.rebin(edges=new_bin)}]
                    drawstack += bkg[key][sample][{'region': region,'TvsQCD':tag_idx,'ut': hist.rebin(edges=new_bin)}].values()
                except:
                    print("Start adding the stack mc & data")
                    mcstack = bkg[key][sample][{'region': region,'TvsQCD':tag_idx,'ut': hist.rebin(edges=new_bin)}]
                    drawstack = bkg[key][sample][{'region': region,'TvsQCD':tag_idx,'ut': hist.rebin(edges=new_bin)}].values()
            print(target.axes.edges[0])
            for sample in stack[stack_map]:
                if region not in bkg[key][sample].axes["region"]:
                    print(f"{region} is not in {sample}'s hist axis")
                    continue
                ax.hist(bins[:-1], bins, weights=drawstack, histtype='stepfilled', label=sample, linewidth=1.5, color=cat_colors[sample], edgecolor=(0,0,0,0.3), alpha=1.0)
                drawstack -= bkg[key][sample][{'region': region,'TvsQCD':tag_idx,'ut': hist.rebin(edges=new_bin)}].values()

            error_opts = {
                'step': 'post',
                'label': 'Stat. Unc.',
                'hatch': '///',
                'facecolor': 'None',
                'edgecolor': (0, 0, 0, 0.3),
                'linewidth': 0
            }
            if region != 'sr':
                if region == 'sr' or 'mcr' in region:
                    data[key]['MET'][{'region': region,'TvsQCD':tag_idx,'ut': hist.rebin(edges=new_bin)}].plot1d(ax=ax, histtype='errorbar', color='k', stack=False, label='Data')
                    datastack = data[key]['MET'][{'region': region,'TvsQCD':tag_idx,'ut': hist.rebin(edges=new_bin)}]
                else:
                    data[key]['EGamma'][{'region': region,'TvsQCD':tag_idx,'ut': hist.rebin(edges=new_bin)}].plot1d(ax=ax, histtype='errorbar', color='k', stack=False, label='Data')
                    datastack = data[key]['EGamma'][{'region': region,'TvsQCD':tag_idx,'ut': hist.rebin(edges=new_bin)}]

            # Signal overlay: skip silently if a sample (or region within it) isn't present.
            for sig_key, color in [(sig_key1, 'c'), (sig_key2, 'r')]:
                if sig_key not in sig.get(key, {}):
                    print(f'Signal sample {sig_key} not found for {key}, skipping')
                    continue
                if region not in sig[key][sig_key].axes['region']:
                    print(f'{region} is not in signal sample {sig_key}\'s hist axis, skipping')
                    continue
                sig[key][sig_key][{'region': region,'TvsQCD':tag_idx,'ut': hist.rebin(edges=new_bin)}].plot1d(ax=ax, histtype='step', color=color, stack=False, label=sig_key)
            leg = ax.legend(ncol=3, loc='upper left', fontsize=12, title=f"{leg_reg[region]}, {tag_label[tag_name]}")
            ax.set_yscale('log')
            ax.set_ylim(0.01, 1000000)
            ax.set_xlabel(None)
            rax.set_xlabel(key)

            if region != 'sr':
                ds = datastack.view().value
                ms = mcstack.view().value
                ratio = ds / ms

                rax.errorbar(
                    x=target.axes.edges[0][:-1] + np.diff(target.axes.edges[0]) / 2,
                    y=ratio,
                    yerr=np.sqrt(datastack.view().variance) / mcstack.view().value,
                    fmt="o",
                    color="k",
                )

                rax.fill_between(
                    bins,
                    # append the last edge to the bins
                    (1.0 - np.sqrt(mcstack.view().variance) / mcstack.view().value).tolist() + [1.0],
                    (1.0 + np.sqrt(mcstack.view().variance) / mcstack.view().value).tolist() + [1.0],
                    step="post",
                    color="k",
                    alpha=0.5,
                )

            if not key == 'template':
                lab = bkg[key][sample].axes[key].label
                rax.set_xlabel(lab)
            rax.set_ylabel('Data/MC')
            rax.set_ylim(0.5, 1.5)
            if 'ut' in key:
                ax.set_xlim(300, 1000)
                rax.set_xlabel('$U_{T}$ [GeV]')

            rax.axhline(1, color='k', linestyle='--',alpha=0.5)
            rax.grid(True)

            fig.savefig(f"{outdir}/{region}_{key}_{whathist_tag}_bins{bin_tag}_{tag_name}.png")
            plt.close(fig)

    print(f"\nDone with {tag_name}. Plots saved in {outdir}/")
