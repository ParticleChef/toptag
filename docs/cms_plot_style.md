# CMS Analysis Style Guide

## Plot Style
- Library: matplotlib + mplhep (CMS style: `mplhep.style.use("CMS")`)
- Always call `mplhep.cms.label()` with luminosity and energy on every plot

## Stack Plot Style
- Histogram style: `stepfilled` via `matplotlib.axes.Axes.stairs(fill=True)`
- Bin outline: `edgecolor="black"`, `linewidth=1.5`

## Stack Plot Defaults
- Collision energy of Run3 CMS: √s = 13.6 TeV
- Luminosity: 
  - 2022preEE: 7.99 fb⁻¹
  - 2022postEE: 26.68 fb⁻¹
  - 2023preBPix: 17.96 fb⁻¹
  - 2023postBPix: 9.68 fb⁻¹
- Process color scheme :
  - QCD multijet: `#2ca02c`
  - W+jets: `#1f77b4`
  - Z(ll)+jets: `#aec7e8`
  - Z(nunu)+jets: `#ffbb78`
  - tt̄: `#ff7f0e`
  - Single top: `#e377c2`
  - Diboson: `#9467bd`
  - Gamma+jets: ``#98df8a"
  - Signal: dashed black line, not stacked
- y-axis range: 10^-2 to 10^6
- Order of stacks(First is top, Last is bottom):
  - 'sr':['Z ($\\nu\\nu$) + Jets','W ($\\ell\\nu$) + Jets', 'VV','TT','G + Jets','ST', 'Z ($\\ell\\ell$) + Jets', 'QCD Multijet' ]
  - 'tecr', 'tmcr', 'wecr', 'wmcr' :['W ($\\ell\\nu$) + Jets', 'TT','QCD Multijet', 'VV','ST', 'Z ($\\ell\\ell$) + Jets','Z ($\\nu\\nu$) + Jets', 'G + Jets']
  - 'zecr', 'zmcr':['Z ($\\ell\\ell$) + Jets','W ($\\ell\\nu$) + Jets', 'TT','QCD Multijet', 'VV','ST',  'Z ($\\nu\\nu$) + Jets', 'G + Jets']
  - 'gcr':['G + Jets','W ($\\ell\\nu$) + Jets', 'TT','QCD Multijet', 'VV', 'ST',  'Z ($\\ell\\ell$) + Jets','Z ($\\nu\\nu$) + Jets']

## Legend
- Title: CR/SR name (use Roman alphabet for mu and gamma): SR, tt̄(e) CR, tt̄(mu) CR, W(e) CR, W(mu) CR, Z(ee) CR, Z(mumu) CR, gamma CR
- Title position: top-center inside the legend box (`legend.get_title().set_ha("center")`)
- Title style: bold (`fontweight="bold"`), fontsize 10

## Data/MC Ratio Panel
- Always include ratio panel below stack plot
- y-axis range: 0.5–1.5
- Horizontal line at 1.0

## General
- Save plots as PDF + PNG (300 dpi)
- Output directory: `plots_{year}_{last 4 number of input file}/`
