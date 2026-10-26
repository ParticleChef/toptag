# Top-Tagging SF Method

Physics background and formulas for the top-tagging scale factors, adapted from
the Run 2 method in `toptag_AN2020_061.pdf` (Section 6.6). Moved out of the main
[`README.md`](../README.md) on 2026-10-25. For how the formulas are implemented in
the server pipeline, see [`scalefactors.md`](scalefactors.md); for the measured
numbers, see [`results.md`](results.md).

---

## Tagger

- **Algorithm**: `parT` (particleTransformer)
- **Applied to**: leading AK15 PFPuppi fatjet in the event
- **Discriminant branch**: `parT_TvsQCD` (top vs QCD score)
- **Working point**: `TvsQCD > 0.33`, taken from the Run 3 ROC (see the ROC step in the README)
  - It was meant to be the 1% background-efficiency WP (~67% signal efficiency),
    but `roc_curve.py` printed mirrored thresholds. That WP is really `TvsQCD > 0.666`.
  - At 0.33 (2022preEE, signal-region selection, mono-top vs Z→νν): signal efficiency
    ~81%, background efficiency ~3.3%. See `reports/report_toptag_WP_ROC_2022preEE_20261025.md`.
  - Keeping 0.33 or moving to ~0.67 is still open.
- The AN's 2018 (Run 2) calibration used `TvsQCD > 0.26`. This codebase uses 0.33
  everywhere (`histmaker.py`'s `TVSQCD_WP`, matching `dev_run3.py:74-79`).

## AK15 Jet Classification

Each AK15 jet is categorized based on how many quarks from a hadronic top decay
are clustered into it (using ΔR ≤ 1.5 matching) and whether b-flavor is found
(via ghost-hadron matching):

| Class | Description |
|-------|-------------|
| 1 | No top-decay quark, no b-flavor |
| 2 | No top-decay quark, b-flavor found |
| 3 | One top-decay quark, no b-flavor |
| 4 | One top-decay quark, b-flavor found |
| 5 | Two top-decay quarks, no b-flavor |
| 6 | Two top-decay quarks, b-flavor found |
| 7 | Three top-decay quarks |
| 8 | More than three top-decay quarks |

These are merged into **three aggregated classes** for calibration:

| Aggregated Class | Original Classes | Description |
|-----------------|-----------------|-------------|
| **Top-jet** | 5, 6, 7, 8 | Fully/mostly merged hadronic top |
| **B-jet**   | 2, 4         | b-flavor present, no full top merge |
| **QCD-jet** | 1, 3         | Pure QCD, no b-flavor |

In the `.scaled` files the class is built from the `n_qinfj` / `n_binfj` axes of
`fjpt_cat`, including their overflow bins (class 8 = n_qinfj ≥ 4; "b" = n_binfj ≥ 1).

## pT Bins

Calibration is performed in four bins of the leading AK15 jet pT:

| Bin | pT range |
|-----|----------|
| 1 | 250 ≤ pT < 325 GeV |
| 2 | 325 ≤ pT < 400 GeV |
| 3 | 400 ≤ pT < 600 GeV |
| 4 | 600 GeV ≤ pT (2000 GeV upper edge in `fjpt_bin` / `fjpt_cat`; 3000 GeV in `scalefactors.py`) |

---

## Scale Factor Derivation

### Common: Data/MC Normalization Factor α

In each pT bin, before applying any top-tagging requirement:

```
α = N(Data) / N(MC)
```

This ratio is taken from the relevant control region and corrects for pre-existing
data/MC normalization differences independent of the tagger.

**Definition (AN2020-061 Sec. 6.6.1):** N(MC) is the
**total expected yield of all background processes** in that region, independent of
the top-tagging status and the jet class. Data in a CR contains every process, so the
pre-tag normalization must compare it with the full MC:

| Calibration | Region | N(MC) in α | Class yields (eff_MC, subtracted terms) |
|-------------|--------|------------|-----------------------------------------|
| QCD-jet | `gcr` | all bkg | G + Jets (AN's γ+jets-dominance approximation) |
| Top-jet | `tmcr` + `tecr` | all bkg | all bkg |

N(Data) and N(MC) come from the `fjpt_bin` histogram, and the class-split yields
from `fjpt_cat`. Because α normalizes the total MC to data,
N(Data) − α·N(other classes) = α·N(target class), so ε_Data / ε_MC reduces to the
SF formulas below. For the QCD-jet case this holds only approximately, to the extent
that γ+jets dominates the γ CR.

> **Superseded (2026-09-28, report revision 2):** N(MC) was temporarily taken as the
> calibrated process only (G + Jets / TT + ST). This does not follow the AN: it scales
> the calibrated process up to absorb all other processes in data. The
> `old_results/tmp_effSF_20260928/` scripts and the revision-2 report still use it.

### QCD-jet Calibration (from γ+jets CR)

The γ+jets control region is dominated by QCD-jets and provides the QCD mistag
scale factor.

**MC efficiency**:
```
ε_QCD-jet,MC = N(γ+jets, QCD-jet, tagged) / N(γ+jets, QCD-jet)
```

**Data efficiency** (B-jet contribution subtracted):
```
ε_QCD-jet,Data = [N(Data, tagged) - α·N(γ+jets, B-jet, tagged)]
               / [N(Data) - α·N(γ+jets, B-jet)]
```

**Scale factor**:
```
SF_QCD-jet = ε_QCD-jet,Data / ε_QCD-jet,MC
           = [N(Data, tagged) - α·N(γ+jets, B-jet, tagged)]
           / [α·N(γ+jets, QCD-jet, tagged)]
```

### Top-jet Calibration (from tt̄ e + μ CRs combined)

The tt̄ electron and muon control regions are used together to derive the
Top-jet tag efficiency scale factor. "MC" below means **all background processes**
in the tt̄ CRs (see the α definition above).

**MC efficiency**:
```
ε_Top-jet,MC = N(MC, Top-jet, tagged) / N(MC, Top-jet)
```

**Data efficiency** (B/QCD-jet contribution subtracted using previously derived SFs):
```
ε_Top-jet,Data = [N(Data, tagged) - α·N(MC, B/QCD-jet, tagged)]
               / [N(Data) - α·N(MC, B/QCD-jet)]
```
The SFs correct the **tagged** subtraction term only,
α·(SF_B·N(MC, B-jet, tagged) + SF_QCD·N(MC, QCD-jet, tagged)). They correct the tag
rate, not the pre-tag yield, so the denominator uses the uncorrected α·N(MC, B/QCD-jet).
The same holds for SF_B in the QCD-jet formula.

Two variants are computed: *raw* (SF_QCD = 1) and *SF_QCD-corrected* (SF_QCD from
the γ CR, per pT bin). The corrected one is written to the SF JSON.

**Scale factor**:
```
SF_Top-jet = ε_Top-jet,Data / ε_Top-jet,MC
           = [N(Data, tagged) - α·N(MC, B/QCD-jet, tagged)]
           / [α·N(MC, Top-jet, tagged)]
```

### B-jet Scale Factor

No dedicated control region is available for the B-jet class. A conservative
per-bin scale factor is assumed:

```
SF_B-jet = 1 ± 0.5
```

### Uncertainties

Propagated now: statistics (data and MC), SF_B = 1 ± 0.5, a 50% normalization
uncertainty on the G + Jets B-jet MC (γ CR) and, in the corrected variant, the
SF_QCD uncertainty. Inputs are treated as uncorrelated (one-at-a-time
finite-difference propagation). No JES, parton-shower, top-pT or scale/PDF
systematics yet. The Run 2 reference tables are in
[`syst_design_computing.md`](syst_design_computing.md).

---

## Event Weight

The top-tagging weight is applied analogously to the AK4 b-tagging weight
(promotion/demotion method). For the leading AK15 jet:

```
w_top-tag = ∏_{j∈tagged}   ( ε_j,Sim · SF_j / ε_j,Sim )
          × ∏_{j∈untagged} ( (1 - ε_j,Sim · SF_j) / (1 - ε_j,Sim) )
```

Since only the leading AK15 jet enters the product, this simplifies to:

- **Tagged jet**:   `w = SF_j`
- **Untagged jet**: `w = (1 - ε_j,Sim · SF_j) / (1 - ε_j,Sim)`

where `j` is the aggregated class (Top/B/QCD) of the leading AK15 jet and
the SF and efficiency are taken from the appropriate pT bin.
