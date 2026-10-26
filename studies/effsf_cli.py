"""Command-line options shared by the alpha / eff_MC / eff_Data+SF scripts.

Each script takes the .scaled file name with -f, e.g.

    -f hadmonotop2023_toptagTEST4        (".scaled" is optional)

The era (luminosity label, output sub-directory) is inferred from the year in
the name unless -e is given, and the plots go to plots_effSF_<version>/, where
<version> is the part of the name after "hadmonotop<year>_" (toptagTEST4).
"""

import argparse
import os
import re
from types import SimpleNamespace

# toptag/ -- the .scaled files and the plot directories live here
WORK_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# fb^-1, same values as toptag/histmaker.py LUMIS
LUMIS = {"2022preEE": 7.99, "2022postEE": 26.68,
         "2023preBPix": 17.96, "2023postBPix": 9.68}

# Era assumed when -e is not given, keyed by the year (+ EE/BPix suffix) in
# the file name, e.g. hadmonotop2022EE_toptagTEST3 -> 2022postEE
DEFAULT_ERA = {"2022": "2022preEE", "2022EE": "2022postEE",
               "2023": "2023preBPix", "2023BPix": "2023postBPix"}

NAME_RE = re.compile(r"hadmonotop(?P<year>\d{4}(?:EE|BPix)?)_(?P<version>.+)$")


def base_parser(description, alpha_var=True):
    """ArgumentParser with the -f / -e (and optionally --alpha-var) options."""
    p = argparse.ArgumentParser(description=description)
    p.add_argument("-f", "--file", required=True,
                   help="input .scaled file, e.g. hadmonotop2023_toptagTEST4"
                        " (extension optional; looked up in the current"
                        " directory, then in toptag/)")
    p.add_argument("-e", "--era", choices=sorted(LUMIS),
                   help="data-taking era (default: from the year in the file"
                        f" name, {DEFAULT_ERA})")
    if alpha_var:
        p.add_argument("--alpha-var", default="fjpt_bin",
                       help="histogram for the data-side yields and alpha"
                            " (default fjpt_bin; old 2022 test files: fjpt)")
    return p


def resolve_file(name):
    """Absolute path of the .scaled file: as given, else inside toptag/."""
    if not name.endswith(".scaled"):
        name += ".scaled"
    for cand in (os.path.abspath(name), os.path.join(WORK_DIR, name)):
        if os.path.isfile(cand):
            return cand
    raise SystemExit(f"{name}: not found in {os.getcwd()} or {WORK_DIR}")


def setup(args):
    """Resolve file, era, lumi and output base from parsed args.

    Also changes into toptag/, so output paths are the same no matter where
    the script is run from.
    """
    path = resolve_file(args.file)
    stem = os.path.basename(path)[:-len(".scaled")]
    m = NAME_RE.match(stem)
    era = args.era or (DEFAULT_ERA.get(m.group("year")) if m else None)
    if era is None:
        raise SystemExit(f"cannot infer the era from '{stem}' -- pass -e")
    os.chdir(WORK_DIR)
    return SimpleNamespace(file=path, era=era, lumi=LUMIS[era],
                           outbase=f"plots_effSF_{m.group('version') if m else stem}")
