"""AK15 fatjet -> truth-level Top/B/QCD classification for top-tagging calibration.

Implements the 8-class -> 3-class scheme of docs/toptag_AN2020_061.pdf Sec 6.6
(see README.md "AK15 Jet Classification"): each AK15 jet is categorized by how
many quarks from a hadronically-decaying generator-level top (b + the two
light/charm quarks from W->qq') fall within DeltaR <= 1.5, combined with
whether a ghost-matched B hadron is found in the jet (nBHadrons).

Inputs are coffea NanoEvents collections loaded with the repo's
CustomNanoAODSchema (analysis/libs/mycoffeav2.py):
  fj      - AK15PuppiJet collection (any shape [event, jet], e.g. leading fatjet)
  genpart - GenPart collection (needs genPartIdxMother, i.e. .distinctChildren /
            .distinctParent, standard NanoAOD content)

MC truth only: do not call on data.
"""

import awkward as ak

DELTA_R_MATCH = 1.5

# class8 -> aggregated 3-class label, per README.md "AK15 Jet Classification" table
AGGREGATE_MAP = {
    1: "qcd",  # no top-decay quark, no b-flavor
    2: "b",    # no top-decay quark, b-flavor found
    3: "qcd",  # one top-decay quark, no b-flavor
    4: "b",    # one top-decay quark, b-flavor found
    5: "top",  # two top-decay quarks, no b-flavor
    6: "top",  # two top-decay quarks, b-flavor found
    7: "top",  # three top-decay quarks
    8: "top",  # more than three top-decay quarks
}


def get_hadronic_top_decay_quarks(genpart):
    """Return (b_quarks, w_quarks) from hadronically-decaying top quarks.

    Both are jagged per-event GenPart collections (shape [event, n]), flat
    across however many hadronic tops are in the event (0, 1, or 2 for
    semileptonic/all-hadronic ttbar). w_quarks holds both daughters of each
    hadronic W, so len(w_quarks) per event == 2 * len(b_quarks) per event.
    """
    is_top = (abs(genpart.pdgId) == 6) & genpart.hasFlags(["fromHardProcess", "isLastCopy"])
    tops = genpart[is_top]  # [event, ntop]

    children = tops.distinctChildren  # [event, ntop, nchild]
    b_of_top = ak.firsts(children[abs(children.pdgId) == 5], axis=2)  # [event, ntop]
    w_of_top = ak.firsts(children[abs(children.pdgId) == 24], axis=2)  # [event, ntop]

    w_children = w_of_top.distinctChildren  # [event, ntop, nwchild]
    w_quarks = w_children[abs(w_children.pdgId) <= 4]  # light/charm only -> hadronic W
    is_hadronic_top = ak.num(w_quarks, axis=2) >= 2

    b_quarks = b_of_top[is_hadronic_top]  # [event, n_hadronic_top]
    w_quarks = ak.flatten(w_quarks[is_hadronic_top], axis=1)  # [event, 2 * n_hadronic_top]

    return b_quarks, w_quarks


def count_matched_top_quarks(fj, genpart):
    """Number of hadronic-top-decay quarks within DeltaR <= 1.5 of each jet.

    Returns an int array shaped like fj (e.g. [event, jet]), summed over all
    hadronic tops in the event (an all-hadronic ttbar event can have two).
    """
    b_quarks, w_quarks = get_hadronic_top_decay_quarks(genpart)
    n_b = ak.sum(fj.metric_table(b_quarks) <= DELTA_R_MATCH, axis=2)
    n_w = ak.sum(fj.metric_table(w_quarks) <= DELTA_R_MATCH, axis=2)
    return n_b + n_w


def classify_ak15_jets(fj, genpart):
    """Classify each AK15 jet into the 8-class scheme and the aggregated 3-class label.

    Returns (class8, class3):
      class8 - int array [1..8], shaped like fj, per the AN table.
      class3 - string array ("top"/"b"/"qcd"), shaped like fj.
    """
    n = count_matched_top_quarks(fj, genpart)
    has_b = fj.nBHadrons > 0

    class8 = ak.where(
        n == 0, ak.where(has_b, 2, 1),
        ak.where(
            n == 1, ak.where(has_b, 4, 3),
            ak.where(
                n == 2, ak.where(has_b, 6, 5),
                ak.where(n == 3, 7, 8),
            ),
        ),
    )

    class3 = ak.where(n >= 2, "top", ak.where(has_b, "b", "qcd"))

    return class8, class3


def aggregate_class(class8):
    """Map an int class8 array/value (1..8) to its 3-class label, via AGGREGATE_MAP."""
    class3 = ak.where(
        (class8 == 5) | (class8 == 6) | (class8 == 7) | (class8 == 8), "top",
        ak.where((class8 == 2) | (class8 == 4), "b", "qcd"),
    )
    return class3


if __name__ == "__main__":
    # Self-test: classify the leading good AK15 jet in the local TTto4Q sample
    # (all-hadronic ttbar -> should see a sizable "top" fraction). Run this
    # after any change to confirm the GenPart navigation / branch names still
    # match the ntuples before wiring classifiers.py into histmaker.py.
    import os
    import sys

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "analysis"))
    from coffea.nanoevents import NanoEventsFactory
    from libs.mycoffeav2 import CustomNanoAODSchema
    from utils.ids import isGoodAK15

    sample = os.path.join(
        os.path.dirname(__file__),
        "samples/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/nano_109.root",
    )
    events = NanoEventsFactory.from_root(
        sample, treepath="Events", schemaclass=CustomNanoAODSchema
    ).events()

    fj = events.AK15PuppiJet
    fj_good = fj[isGoodAK15(fj)]
    leading_fj = ak.firsts(fj_good)[:, None]  # keep jagged shape [event, 1] for metric_table

    class8, class3 = classify_ak15_jets(leading_fj, events.GenPart)
    class8 = ak.flatten(class8)
    class3 = ak.flatten(class3)

    print(f"n events        : {len(events)}")
    print(f"n with leading fj: {ak.sum(~ak.is_none(ak.firsts(leading_fj)))}")
    print("class3 breakdown :")
    for label in ("top", "b", "qcd"):
        print(f"  {label:4s}: {ak.sum(class3 == label)}")
    print("class8 breakdown :")
    for c in range(1, 9):
        print(f"  {c}: {ak.sum(class8 == c)}")
