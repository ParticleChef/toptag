# Run 2 Systematics, Design Rules and Computing Environment

Reference material moved out of the main [`README.md`](../README.md) on 2026-09-28.

- **Systematic Uncertainties**: Run 2 (2018) reference values, taken from
  AN-2020/061 (`toptag_AN2020_061.pdf`, Section 6.6, most likely Tables 31/32; not
  re-checked against the PDF). They are **not** Run 3 measurements from this
  project. Only the "dominant sources" are listed, so the rows do not add up in
  quadrature to the totals. The current Run 3 results propagate statistics,
  BJetMistag and BJetNorm only (see `reports/report_toptag_effSF_alpha_20260928.md`).
- **Design Rules** and **Computing Environment**: conventions for the server
  pipeline (`histmaker.py` → `scalefactors.py` → `validate.py` → `weights.py`)
  and for running on the hep server.

---

## Systematic Uncertainties (Run 2 reference, AN-2020/061)

### QCD-jet mistagging SF (dominant sources)

| Source | 250–325 GeV | 325–400 GeV | 400–600 GeV | 600+ GeV |
|--------|-------------|-------------|-------------|----------|
| statistics | 15.1% | 6.2% | 4.1% | 6.4% |
| final-state PS | 0.7% | 5.0% | 7.8% | 13.3% |
| BJetMistag | 4.6% | 5.8% | 5.7% | 6.1% |
| BJetNorm | 2.7% | 3.9% | 3.8% | 3.9% |
| **total** | **18.1%** | **12.2%** | **12.1%** | **16.7%** |

### Top-jet tagging SF (dominant sources)

| Source | 250–325 GeV | 325–400 GeV | 400–600 GeV | 600+ GeV |
|--------|-------------|-------------|-------------|----------|
| statistics | 9.6% | 5.4% | 2.4% | 6.7% |
| TopPt reweighting | 5.6% | 6.1% | 6.7% | 9.4% |
| QCDMistag | 5.1% | 2.3% | 1.4% | 5.0% |
| BJetMistag | 7.0% | 4.3% | 3.4% | 7.6% |
| **total** | **16.4%** | **11.5%** | **8.8%** | **15.5%** |

---

## Design Rules

These rules ensure every step can be reproduced by anyone on the team,
without relying on any prior conversation or external context.

1. **Every script is standalone and runnable from the command line.**
   No step requires opening a notebook or calling internal functions manually.
   Each script uses `argparse` and prints its inputs, outputs, and progress.

2. **All configurable parameters live at the top of each script** in a clearly
   labeled `CONFIG` block. No magic numbers buried in logic.

3. **Each step reads from files and writes to files.** Intermediate results
   (histograms, SFs) are saved to disk so any step can be re-run independently.

4. **No implicit ordering.** The README lists the exact commands in order.
   Running them top-to-bottom from a clean checkout must always work.

5. **No `decaf/` code is modified.** This project only reads from `decaf/`
   as a reference. All new code lives under `toptag/`.

6. **Output files are versioned by year** (e.g. `toptag_sf_2022pre.json.gz`).
   Never overwrite without an explicit `--year` argument.

---

## Computing Environment

This is only for hep server.

This server has no batch/condor nodes. All jobs run locally.

| Resource | Value |
|----------|-------|
| CPU | AMD EPYC 9354 × 2 = **128 cores** |
| RAM | **251 GB** total (~239 GB available) |
| Local disk (`/`) | 1.8 TB NVMe, ~277 GB free |
| Home (`/home`) | NFS 2.3 TB, ~1.2 TB free |
| Input files | **Local** at `/data/mc/privatemc/<year>/...` |

The `2022_private_v1.json.gz` metadata contains ~2,994 datasets and ~29,643 files
in total. For the toptag calibration we only need:
- γ+jets CR: ~287 dataset chunks
- tt̄ CR: ~435 dataset chunks

### Recommended: coffea `futures_executor` with 40 workers (default)

Because all input files are on local NVMe (not XROOTD), the bottleneck is CPU,
not I/O. Running coffea's `futures_executor` with 40 workers in a single process
is the recommended approach for this subset of datasets:

- Uses ~40 cores out of 128 — leaves ~88 cores free for other users on the server
- Each worker uses ~2–3 GB RAM → ~80–120 GB total, well within the 239 GB available
- No monitoring loop or per-dataset log management needed
- If it crashes, coffea `.futures` checkpoint files allow resuming

**Always check current load before choosing a worker count:**
```bash
htop        # interactive view — press q to quit
# or a quick snapshot:
uptime      # load average over 1/5/15 min; keep load < ~100 on this machine
free -h     # check available RAM
```

Rule of thumb for `--workers N`:
- Server is quiet (load < 20):  `--workers 60`
- Normal use (load 20–60):      `--workers 40`  ← default
- Server is busy (load > 60):   `--workers 20`

Contrast with `nohup_job_new.py` (used for the full analysis): that runs 50 concurrent
single-worker processes, which is fine for 2,994 datasets but unnecessarily complex
for the ~700 datasets needed here.
