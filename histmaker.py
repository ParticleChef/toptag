#!/usr/bin/env python
"""histmaker.py - Step 1 of the top-tagging SF pipeline (see README.md).

Runs a coffea processor over NanoAOD in the gamma+jets (gcr) and ttbar e/mu
(tecr/tmcr) control regions and fills a region x fjpt x TvsQCD x jetclass
histogram, which scalefactors.py turns into efficiencies and scale factors.

Object IDs, region cuts, and event weights are ported from
analysis/processors/dev_run3.py (READ ONLY reference; not modified here),
restricted to the gcr/tmcr/tecr regions. Jet truth classification
(Top-jet/B-jet/QCD-jet) is done by classifiers.py.

Usage:
    python histmaker.py --year 2022pre --workers 40 --output hists/toptag_hists_2022pre.coffea
"""

import argparse
import gzip
import json
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

import awkward as ak
import hist
import numpy as np
from coffea import processor
from coffea.analysis_tools import PackedSelection, Weights
from coffea.nanoevents.methods import vector
from coffea.util import load, save

TOPTAG_ROOT = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.join(TOPTAG_ROOT, "analysis")
sys.path.insert(0, ANALYSIS_DIR)

from classifiers import classify_ak15_jets  # noqa: E402
from libs.mycoffeav2 import CustomNanoAODSchema  # noqa: E402

# ============================================================
# CONFIG
# ============================================================

# fb^-1, https://twiki.cern.ch/twiki/bin/viewauth/CMS/LumiRecommendationsRun3
LUMIS = {
    "2022pre": 7.99,
    "2022post": 26.68,
    "2023pre": 17.96,
    "2023post": 9.68,
}

# Run 3 particleTransformer TvsQCD working point, dev_run3.py:74-79
TVSQCD_WP = {
    "2022pre": 0.33,
    "2022post": 0.33,
    "2023pre": 0.33,
    "2023post": 0.33,
}

# AK15 leading-jet pT bins, AN-2020/061 Sec 6.6.1.3 / README.md; 3000 stands in for "600-Inf"
PT_EDGES = [250, 325, 400, 600, 3000]

MET_TRIGGERS = ["PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", "PFMETNoMu120_PFMHTNoMu120_IDTight"]
ELECTRON_TRIGGERS = ["Ele30_WPTight_Gsf", "Photon200"]
PHOTON_TRIGGERS = ["Photon200"]

MET_FILTERS = [
    "goodVertices", "globalSuperTightHalo2016Filter", "EcalDeadCellTriggerPrimitiveFilter",
    "BadPFMuonFilter", "BadPFMuonDzFilter", "hfNoisyHitsFilter", "eeBadScFilter",
]

LUMI_MASK_FILES = {
    "2022pre": "data/lumiMask/Cert_Collisions2022_355100_362760_Golden.json",
    "2022post": "data/lumiMask/Cert_Collisions2022_355100_362760_Golden.json",
    "2023pre": "data/lumiMask/Cert_Collisions2023_366442_370790_Golden.json",
    "2023post": "data/lumiMask/Cert_Collisions2023_366442_370790_Golden.json",
}

# Datasets that populate each CR, matched by substring against the dataset key -
# same lists as AnalysisProcessor._samples in dev_run3.py:63-72, restricted to
# the 3 regions this calibration needs.
REGION_SAMPLES = {
    "gcr": ("TT", "WtoLNu-2Jets", "DYto2L-2Jets", "Zto2Nu-2Jets", "GJ", "ST", "WW", "ZZ", "WZ", "QCD", "EGamma"),
    "tmcr": ("TT", "WtoLNu-2Jets", "DYto2L-2Jets", "Zto2Nu-2Jets", "GJ", "ST", "WW", "ZZ", "WZ", "QCD", "JetMET"),
    "tecr": ("TT", "WtoLNu-2Jets", "DYto2L-2Jets", "Zto2Nu-2Jets", "GJ", "ST", "WW", "ZZ", "WZ", "QCD", "EGamma"),
}

# dataset -> metadata file under analysis/metadata/ (no decaf/ checkout in this repo).
# VERIFY against analysis/run.py's usual --metadata argument for a given year/era
# before trusting these on the hep server - not all eras have a single combined
# MC+data fileset (2023post ships separate vMC/vDATA files, hence the list).
METADATA_FILES = {
    "2022pre": ["2022_private_v2.json.gz"],
    "2022post": ["2022EE_private_v2.json.gz"],
    "2023pre": ["2023_privateLPC_v1.json.gz"],
    "2023post": ["2023BPix_private_vMC.json.gz", "2023BPix_private_vDATA.json.gz"],
}


class TopTagHistProcessor(processor.ProcessorABC):
    """gcr/tmcr/tecr selection + region x fjpt x TvsQCD x jetclass histogram.

    Trimmed from analysis/processors/dev_run3.py: same object IDs, region
    cuts, and event weights, restricted to the 3 regions this calibration
    needs (drops sr/wmcr/wecr/zmcr/zecr and their unique variables).
    """

    def __init__(self, year, xsec, corrections, ids, common):
        self._year = year
        self._lumi = 1000.0 * float(LUMIS[year])
        self._xsec = xsec
        self._wp = TVSQCD_WP[year]
        self._lumimask = load_lumimask(year)
        self._corrections = corrections
        self._ids = ids
        self._common = common

        self._samples = REGION_SAMPLES

        self.make_output = lambda: {
            "sumw": 0.0,
            "jetpt_class": hist.Hist(
                hist.axis.StrCategory([], name="region", growth=True),
                hist.axis.Variable(PT_EDGES, name="fjpt", label="AK15 Leading Jet $p_{T}$"),
                hist.axis.Variable([0, self._wp, 1], name="TvsQCD"),
                hist.axis.StrCategory([], name="jetclass", growth=True),
                storage=hist.storage.Weight(),
            ),
        }

    def process(self, events):
        return self.process_shift(events)

    def process_shift(self, events):
        dataset = events.metadata["dataset"]
        caller = dataset.split("____")[0]
        isData = not hasattr(events, "genWeight")

        selected_regions = [
            region for region, samples in self._samples.items()
            if any(sample in dataset for sample in samples)
        ]

        selection = PackedSelection(dtype="uint64")
        weights = Weights(len(events), storeIndividual=True)
        output = self.make_output()

        if not isData:
            output["sumw"] = ak.sum(events.genWeight)

        get_pu_weight = self._corrections["get_pu_weight"]
        get_mu_loose_id_sf = self._corrections["get_mu_loose_id_sf"]
        get_mu_tight_id_sf = self._corrections["get_mu_tight_id_sf"]
        get_mu_loose_iso_sf = self._corrections["get_mu_loose_iso_sf"]
        get_mu_tight_iso_sf = self._corrections["get_mu_tight_iso_sf"]
        get_met_xy_correction = self._corrections["get_met_xy_correction"]
        get_ele_loose_id_sf = self._corrections["get_ele_loose_id_sf"]
        get_ele_tight_id_sf = self._corrections["get_ele_tight_id_sf"]
        get_ele_reco_sf_below20 = self._corrections["get_ele_reco_sf_below20"]
        get_ele_reco_sf_20to75 = self._corrections["get_ele_reco_sf_20to75"]
        get_ele_reco_sf_Above75 = self._corrections["get_ele_reco_sf_Above75"]
        get_photon_id_sf = self._corrections["get_photon_id_sf"]
        get_met_trig_weight = self._corrections["get_met_trig_weight"]
        get_ele_trig_weight = self._corrections["get_ele_trig_weight"]
        get_pho_trig_weight = self._corrections["get_pho_trig_weight"]
        get_nlo_ewk_weight = self._corrections["get_nlo_ewk_weight"]
        get_btag_weight = self._corrections["get_btag_weight"]
        get_ttbar_weight = self._corrections["get_ttbar_weight"]
        get_jec_correction = self._corrections["get_jec_correction"]
        get_fjec_correction = self._corrections["get_fjec_correction"]

        isLooseMuon = self._ids["isLooseMuon"]
        isTightMuon = self._ids["isTightMuon"]
        isLooseElectron = self._ids["isLooseElectron"]
        isTightElectron = self._ids["isTightElectron"]
        isLoosePhoton = self._ids["isLoosePhoton"]
        isTightPhoton = self._ids["isTightPhoton"]
        isGoodAK4 = self._ids["isGoodAK4"]
        isGoodAK15 = self._ids["isGoodAK15"]
        isJetVeto = self._ids["isJetVeto"]

        ###
        # Global quantities
        ###
        npv = events.PV.npvsGood
        run = events.run
        met = events.MET
        met["pt"], met["phi"] = get_met_xy_correction(self._year, "MET", isData, met.pt, met.phi, npv)

        ###
        # Objects
        ###
        mu = events.Muon
        mu["isloose"] = isLooseMuon(mu, self._year)
        mu["id_sf"] = ak.where(mu.isloose, get_mu_loose_id_sf(self._year, abs(mu.eta), mu.pt), ak.ones_like(mu.pt))
        mu["iso_sf"] = ak.where(mu.isloose, get_mu_loose_iso_sf(self._year, abs(mu.eta), mu.pt), ak.ones_like(mu.pt))
        mu["istight"] = isTightMuon(mu, self._year)
        mu["id_sf"] = ak.where(mu.istight, get_mu_tight_id_sf(self._year, abs(mu.eta), mu.pt), mu.id_sf)
        mu["iso_sf"] = ak.where(mu.istight, get_mu_tight_iso_sf(self._year, abs(mu.eta), mu.pt), mu.iso_sf)
        mu["T"] = ak.zip({"r": mu.pt, "phi": mu.phi}, with_name="PolarTwoVector", behavior=vector.behavior)
        mu_loose = mu[mu.isloose]
        mu_tight = mu[mu.istight]
        mu_nloose = ak.num(mu_loose, axis=1)
        mu_ntight = ak.num(mu_tight, axis=1)
        leading_mu = ak.firsts(mu_tight)

        e = events.Electron
        sf_below20, _, _ = get_ele_reco_sf_below20(self._year, abs(e.eta), e.pt, e.phi)
        sf_20to75, _, _ = get_ele_reco_sf_20to75(self._year, abs(e.eta), e.pt, e.phi)
        sf_Above75, _, _ = get_ele_reco_sf_Above75(self._year, abs(e.eta), e.pt, e.phi)
        e["reco_sf"] = ak.where(e.pt < 20, sf_below20, ak.where(e.pt < 75, sf_20to75, sf_Above75))
        sf_loose, _, _ = get_ele_loose_id_sf(self._year, abs(e.eta), e.pt, e.phi)
        sf_tight, _, _ = get_ele_tight_id_sf(self._year, abs(e.eta), e.pt, e.phi)
        eleHlt_loose, _, _ = get_ele_trig_weight(self._year, abs(e.eta), e.pt, "Loose")
        eleHlt_tight, _, _ = get_ele_trig_weight(self._year, abs(e.eta), e.pt, "Tight")
        e["isloose"] = isLooseElectron(e, self._year)
        e["id_sf"] = ak.where(e.isloose, sf_loose, ak.ones_like(e.pt))
        e["hlt_sf"] = ak.where(e.isloose, eleHlt_loose, ak.ones_like(e.pt))
        e["istight"] = isTightElectron(e, self._year)
        e["id_sf"] = ak.where(e.istight, sf_tight, e.id_sf)
        e["hlt_sf"] = ak.where(e.istight, eleHlt_tight, e.hlt_sf)
        e["T"] = ak.zip({"r": e.pt, "phi": e.phi}, with_name="PolarTwoVector", behavior=vector.behavior)
        e_loose = e[e.isloose]
        e_tight = e[e.istight]
        e_nloose = ak.num(e_loose, axis=1)
        e_ntight = ak.num(e_tight, axis=1)
        leading_e = ak.firsts(e_tight)

        pho = events.Photon
        pho_sf_loose, _, _ = get_photon_id_sf(self._year, "Loose", abs(pho.eta), pho.pt, pho.phi)
        pho_sf_tight, _, _ = get_photon_id_sf(self._year, "Tight", abs(pho.eta), pho.pt, pho.phi)
        pho["isloose"] = isLoosePhoton(pho, self._year)
        pho["id_sf"] = ak.where(pho.isloose, pho_sf_loose, ak.ones_like(pho.pt))
        pho["istight"] = isTightPhoton(pho, self._year)
        pho["id_sf"] = ak.where(pho.istight, pho_sf_tight, pho.id_sf)
        pho["T"] = ak.zip({"r": pho.pt, "phi": pho.phi}, with_name="PolarTwoVector", behavior=vector.behavior)
        pho_loose = pho[pho.isloose]
        pho_tight = pho[pho.istight]
        pho_nloose = ak.num(pho_loose, axis=1)
        pho_ntight = ak.num(pho_tight, axis=1)
        leading_pho = ak.firsts(pho_tight)

        fj = events.AK15PuppiJet
        rho_density = events.Rho.fixedGridRhoFastjetAll
        fjec_corr = get_fjec_correction(self._year, fj.pt, fj.eta, fj.phi, rho_density, fj.area, run, isData)
        fj["pt"] = fj.pt * fjec_corr
        fj["mass"] = fj.mass * fjec_corr
        fj["isclean"] = (
            ak.all(fj.metric_table(mu_loose) > 1.5, axis=2)
            & ak.all(fj.metric_table(e_loose) > 1.5, axis=2)
            & ak.all(fj.metric_table(pho_loose) > 1.5, axis=2)
        )
        fj["isgood"] = isGoodAK15(fj)
        fj["T"] = ak.zip({"r": fj.pt, "phi": fj.phi}, with_name="PolarTwoVector", behavior=vector.behavior)
        probQCD = fj.ParT_probQCDbb + fj.ParT_probQCDcc + fj.ParT_probQCDb + fj.ParT_probQCDc + fj.ParT_probQCDothers
        probT = fj.ParT_probTopbWqq + fj.ParT_probTopbWcs
        fj["TvsQCD"] = probT / (probT + probQCD)
        fj_good = fj[fj.isgood]
        fj_clean = fj_good[fj_good.isclean]
        fj_nclean = ak.num(fj_clean, axis=1)
        leading_fj = ak.firsts(fj_clean)
        # jagged [event, 0-or-1] shape for classify_ak15_jets' metric_table calls
        leading_fj_jagged = fj_clean[:, :1]

        j = events.Jet
        jec_corr = get_jec_correction(self._year, j.pt, j.eta, j.phi, rho_density, j.area, run, isData)
        j["pt"] = j.pt * jec_corr
        j["mass"] = j.mass * jec_corr
        j["isveto"] = isJetVeto(j, self._year)
        n_j_veto = ak.num(j[j.isveto], axis=1)
        j["T"] = ak.zip({"r": j.pt, "phi": j.phi}, with_name="PolarTwoVector", behavior=vector.behavior)
        j["isgood"] = isGoodAK4(j, self._year)
        j["isclean"] = (
            ak.all(j.metric_table(mu_loose) > 0.4, axis=2)
            & ak.all(j.metric_table(e_loose) > 0.4, axis=2)
        )
        j["isiso"] = ak.all(j.metric_table(leading_fj) > 1.5, axis=2)
        PNetUParTWPs = self._common["btagWPs"]["PNetUParT"][self._year]
        j["isbtagvL"] = j.btagPNetB > PNetUParTWPs["loose"]
        j_good = j[j.isgood]
        j_clean = j_good[j_good.isclean]
        j_iso = j_clean[j_clean.isiso]
        j_btagvL = j_iso[j_iso.isbtagvL]
        j_nclean = ak.num(j_clean, axis=1)
        j_nbtagvL = ak.num(j_btagvL, axis=1)

        ###
        # Recoil (u['sr'] is also the common reference recoil used by the
        # mindphi_<region> cut below, matching dev_run3.py:672-680,1141)
        ###
        u = {
            "sr": met,
            "tmcr": met + leading_mu.T,
            "tecr": met + leading_e.T,
            "gcr": met + leading_pho.T,
        }

        ###
        # MC-only truth: top-pt reweighting, NLO EWK k-factor, PU weight,
        # b-tag SF, and the AK15 jet-class truth label (classifiers.py)
        ###
        class3 = None
        if not isData:
            gen = events.GenPart
            gen["isTop"] = (abs(gen.pdgId) == 6) & gen.hasFlags(["fromHardProcess", "isLastCopy"])
            genTops = gen[gen.isTop]
            nlo = np.ones(len(events), dtype="float")
            if "TTto" in dataset:
                nlo = np.sqrt(get_ttbar_weight(genTops[:, 0].pt) * get_ttbar_weight(genTops[:, 1].pt))
            genW = gen[gen.hasFlags(["fromHardProcess", "isFirstCopy", "isPrompt"]) & (abs(gen.pdgId) == 24)]
            genZ = gen[(abs(gen.pdgId) == 23) & gen.hasFlags(["fromHardProcess", "isLastCopy"])]
            genDY = genZ[genZ.mass > 30]
            genA = gen[
                (abs(gen.pdgId) == 22)
                & gen.hasFlags(["isPrompt", "fromHardProcess", "isLastCopy"])
                & (gen.status == 1)
                & (gen.pt > 100)
            ]
            GenIsoPho = events.GenIsolatedPhoton
            nGenIsoPho = ak.num(GenIsoPho, axis=1)

            nlo_ewk = np.ones(len(events), dtype="float")
            if "WtoLNu" in dataset:
                nlo_ewk = get_nlo_ewk_weight["w"](ak.max(genW.pt, axis=1))
            elif "DY" in dataset:
                nlo_ewk = get_nlo_ewk_weight["dy"](ak.max(genDY.pt, axis=1))
            elif "Zto2Nu" in dataset:
                nlo_ewk = get_nlo_ewk_weight["z"](ak.max(genZ.pt, axis=1))
            elif "GJ" in dataset:
                nlo_ewk = get_nlo_ewk_weight["a"](ak.max(genA.pt, axis=1))

            pu = get_pu_weight(self._year, events.Pileup.nTrueInt)

            btagSF = get_btag_weight("PNetUParT", self._year, "loose", caller).btag_weight(
                j_iso.pt, j_iso.eta, j_iso.hadronFlavour, j_iso.isbtagvL
            )[0]

            weights.add("genw", events.genWeight)
            weights.add("pileup", pu)
            weights.add("nlo_ewk", nlo_ewk)
            weights.add("topptreweighting", nlo)
            weights.add("btagSF", btagSF)

            mask_isoneE = ak.to_numpy((e_ntight == 1) & (e_nloose == 1) & (mu_nloose == 0) & (pho_nloose == 0))
            mask_isoneM = ak.to_numpy((e_nloose == 0) & (mu_ntight == 1) & (mu_nloose == 1) & (pho_nloose == 0))
            mask_isoneG = ak.to_numpy((e_nloose == 0) & (mu_nloose == 0) & (pho_nloose == 1) & (pho_ntight == 1))

            id_sf = np.ones(len(events), dtype="float")
            id_sf_isoneE = ak.to_numpy(ak.fill_none(ak.prod(e_tight.id_sf, axis=1), 1.0))
            id_sf_isoneM = ak.to_numpy(ak.fill_none(ak.prod(mu_tight.id_sf, axis=1), 1.0))
            id_sf_isoneG = ak.to_numpy(ak.fill_none(ak.prod(pho_tight.id_sf, axis=1), 1.0))
            id_sf[mask_isoneE] = id_sf_isoneE[mask_isoneE]
            id_sf[mask_isoneM] = id_sf_isoneM[mask_isoneM]
            id_sf[mask_isoneG] = id_sf_isoneG[mask_isoneG]
            weights.add("ids", id_sf)

            reco_sf = np.ones(len(events), dtype="float")
            reco_sf_isoneE = ak.to_numpy(ak.fill_none(ak.prod(e_tight.reco_sf, axis=1), 1.0))
            reco_sf[mask_isoneE] = reco_sf_isoneE[mask_isoneE]
            weights.add("reco", reco_sf)

            iso_sf = np.ones(len(events), dtype="float")
            iso_sf_isoneM = ak.to_numpy(ak.fill_none(ak.prod(mu_tight.iso_sf, axis=1), 1.0))
            iso_sf[mask_isoneM] = iso_sf_isoneM[mask_isoneM]
            weights.add("iso", iso_sf)

            met_sf_m, _, _ = get_met_trig_weight(self._year, u["tmcr"].r)
            pho_sf, _, _ = get_pho_trig_weight(self._year, ak.firsts(pho_tight).pt)
            trg_sf = np.ones(len(events), dtype="float")
            trg_sf_isoneE = ak.to_numpy(ak.fill_none(ak.prod(e_tight.hlt_sf, axis=1), 1.0))
            trg_sf_isoneM = ak.to_numpy(ak.fill_none(met_sf_m, 1.0))
            trg_sf_isoneG = ak.to_numpy(ak.fill_none(pho_sf, 1.0))
            trg_sf[mask_isoneE] = trg_sf_isoneE[mask_isoneE]
            trg_sf[mask_isoneM] = trg_sf_isoneM[mask_isoneM]
            trg_sf[mask_isoneG] = trg_sf_isoneG[mask_isoneG]
            weights.add("trig", trg_sf)

            _, class3 = classify_ak15_jets(leading_fj_jagged, gen)
            class3 = ak.fill_none(ak.firsts(class3), "none")
        else:
            nGenIsoPho = ak.zeros_like(met.pt, dtype="int64")

        ###
        # Selections
        ###
        lumimask = np.ones(len(events), dtype="bool")
        if isData:
            lumimask = self._lumimask(events.run, events.luminosityBlock)
        selection.add("lumimask", lumimask)

        met_filters = np.ones(len(events), dtype="bool")
        for flag in MET_FILTERS:
            met_filters = met_filters & events.Flag[flag]
        selection.add("met_filters", met_filters)

        for name, paths in (
            ("met_triggers", MET_TRIGGERS),
            ("single_electron_triggers", ELECTRON_TRIGGERS),
            ("single_photon_triggers", PHOTON_TRIGGERS),
        ):
            triggers = np.zeros(len(events), dtype="bool")
            for path in paths:
                if path not in events.HLT.fields:
                    continue
                triggers = triggers | events.HLT[path]
            selection.add(name, ak.to_numpy(triggers))

        if isData:
            selection.add("exclude_wjets_greater_120", np.full(len(events), True))
            selection.add("exclude_wjets_less_120", np.full(len(events), True))
        else:
            genW_all = events.GenPart[
                events.GenPart.hasFlags(["fromHardProcess", "isFirstCopy", "isPrompt"])
                & (abs(events.GenPart.pdgId) == 24)
            ]
            if ("WtoLNu-2Jets" in dataset) and ("PT" in dataset):
                selection.add("exclude_wjets_greater_120", ak.to_numpy(ak.all(genW_all.pt > 120, axis=1)))
            else:
                selection.add("exclude_wjets_greater_120", np.full(len(events), True))
            if ("WtoLNu-2Jets" in dataset) and ("PT" not in dataset):
                selection.add("exclude_wjets_less_120", ak.to_numpy(ak.all(genW_all.pt <= 120, axis=1)))
            else:
                selection.add("exclude_wjets_less_120", np.full(len(events), True))

        if "QCD_PT" in dataset:
            selection.add("QCD_NoGenIsoPho", ak.to_numpy(nGenIsoPho == 0))
        else:
            selection.add("QCD_NoGenIsoPho", np.full(len(events), True))
        if "GJ_PTG" in dataset:
            selection.add("GJet_GenIsoPho", ak.to_numpy(nGenIsoPho >= 1))
        else:
            selection.add("GJet_GenIsoPho", np.full(len(events), True))

        selection.add("isoneM", (e_nloose == 0) & (mu_ntight == 1) & (mu_nloose == 1) & (pho_nloose == 0))
        selection.add("isoneE", (e_ntight == 1) & (e_nloose == 1) & (mu_nloose == 0) & (pho_nloose == 0))
        selection.add("isoneG", (e_nloose == 0) & (mu_nloose == 0) & (pho_nloose == 1) & (pho_ntight == 1))
        selection.add("one_ak15", fj_nclean > 0)
        selection.add("leading_fj250", leading_fj.pt > 250)
        selection.add("noextrab", j_nbtagvL == 0)
        selection.add("extrab", j_nbtagvL > 0)
        selection.add("jetveto", n_j_veto == 0)
        selection.add("met150", met.pt > 150)

        regions = {
            "tmcr": [
                "lumimask", "met_filters", "met_triggers",
                "exclude_wjets_greater_120", "exclude_wjets_less_120",
                "QCD_NoGenIsoPho", "GJet_GenIsoPho",
                "jetveto", "one_ak15", "leading_fj250",
                "isoneM", "met150", "extrab",
            ],
            "tecr": [
                "lumimask", "met_filters", "single_electron_triggers",
                "exclude_wjets_greater_120", "exclude_wjets_less_120",
                "QCD_NoGenIsoPho", "GJet_GenIsoPho",
                "jetveto", "one_ak15", "leading_fj250",
                "isoneE", "met150", "extrab",
            ],
            "gcr": [
                "lumimask", "met_filters", "single_photon_triggers",
                "exclude_wjets_greater_120", "exclude_wjets_less_120",
                "QCD_NoGenIsoPho", "GJet_GenIsoPho",
                "jetveto", "one_ak15", "leading_fj250",
                "isoneG", "noextrab",
            ],
        }

        def normalize(val, cut):
            return ak.to_numpy(ak.fill_none(val[cut], np.nan))

        def normalize_str(val, cut):
            return ak.to_numpy(val[cut])

        for region in list(regions):
            if region not in selected_regions:
                del regions[region]
                continue
            selection.add("recoil_" + region, u[region].r > 350)
            selection.add(
                "mindphi_" + region,
                ak.min(abs(u["sr"].delta_phi(j_clean.T)), axis=1, mask_identity=False) > 0.5,
            )
            selection.add("minDphi_" + region, abs(u[region].delta_phi(leading_fj.T)) > 1.5)
            regions[region] += ["recoil_" + region, "mindphi_" + region, "minDphi_" + region]

        for region, cuts in regions.items():
            cut = selection.all(*cuts)
            weight = weights.weight()[cut]
            class_arr = np.full(int(ak.sum(cut)), "data") if isData else normalize_str(class3, cut)
            output["jetpt_class"].fill(
                region=region,
                fjpt=normalize(leading_fj.pt, cut),
                TvsQCD=normalize(leading_fj.TvsQCD, cut),
                jetclass=class_arr,
                weight=weight,
            )

        scale = 1
        if self._xsec[dataset] != -1:
            scale = self._lumi * self._xsec[dataset]
        output["jetpt_class"] *= scale

        return output

    def postprocess(self, accumulator):
        return accumulator


def load_lumimask(year):
    from coffea.lumi_tools import LumiMask

    return LumiMask(os.path.join(ANALYSIS_DIR, LUMI_MASK_FILES[year]))


# dataset -> physics process grouping, ported verbatim from
# analysis/macros/scale.py bkg_map/data_map (lines 84-114)
BKG_MAP = {
    "WW": ["WW"], "WZ": ["WZ"], "ZZ": ["ZZ"], "QCD Multijet": ["QCD"],
    "TT": ["TT"], "ST": ["ST"],
    r"Z ($\nu\nu$) + Jets": ["Zto2Nu"], r"Z ($\ell\ell$) + Jets": ["DYto2L"],
    r"W ($\ell\nu$) + Jets": ["WtoLNu"], "G + Jets": ["GJ_PTG"],
}
DATA_MAP = {"MET": ["JetMET"], "EGamma": ["EGamma"]}


def select_datasets(samplefiles):
    """Restrict a fileset dict to datasets used by gcr/tmcr/tecr (README Step 1)."""
    wanted = set()
    for samples in REGION_SAMPLES.values():
        wanted.update(samples)
    return {
        dataset: info for dataset, info in samplefiles.items()
        if any(sample in dataset for sample in wanted)
    }


def run_one_dataset(processor_instance, dataset, files, workers):
    filelist = {dataset: files}
    return processor.run_uproot_job(
        filelist,
        "Events",
        processor_instance=processor_instance,
        executor=processor.futures_executor,
        executor_args={"schema": CustomNanoAODSchema, "workers": workers, "skipbadfiles": False},
    )


def group_and_normalize(per_dataset_output, xsec):
    """1/sumw MC normalization + bkg_map/data_map grouping, ported from
    analysis/macros/scale.py:65-160, restricted to the jetpt_class histogram."""
    bkg_hists = {}
    data_hists = {}
    for dataset, output in per_dataset_output.items():
        sumw = output["sumw"]
        h = output["jetpt_class"]
        if xsec.get(dataset, -1) != -1:  # MC
            if sumw == 0:
                print(f"WARNING: sumw==0 for {dataset}, skipping (no events passed genWeight sum)")
                continue
            h = h * (1.0 / sumw)
            for process, patterns in BKG_MAP.items():
                if not any(p in dataset for p in patterns):
                    continue
                bkg_hists[process] = h if process not in bkg_hists else bkg_hists[process] + h
        else:  # data
            for process, patterns in DATA_MAP.items():
                if not any(p in dataset for p in patterns):
                    continue
                data_hists[process] = h if process not in data_hists else data_hists[process] + h
    return {"bkg": bkg_hists, "data": data_hists}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", required=True, choices=list(LUMIS))
    parser.add_argument("--workers", type=int, default=40)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    output_path = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # ids/corrections/common.coffea, and correctionlib files loaded lazily
    # inside them, resolve relative paths against analysis/ - chdir there.
    os.chdir(ANALYSIS_DIR)

    corrections = load("data/corrections.coffea")
    ids = load("data/ids.coffea")
    common = load("data/common.coffea")

    samplefiles = {}
    for metadata_file in METADATA_FILES[args.year]:
        with gzip.open(os.path.join("metadata", metadata_file)) as fin:
            samplefiles.update(json.load(fin))
    xsec = {k: v["xs"] for k, v in samplefiles.items()}
    samplefiles = select_datasets(samplefiles)
    print(f"{len(samplefiles)} datasets selected for gcr/tmcr/tecr")

    processor_instance = TopTagHistProcessor(
        year=args.year, xsec=xsec, corrections=corrections, ids=ids, common=common
    )

    per_dataset_output = {}
    for i, (dataset, info) in enumerate(samplefiles.items(), 1):
        print(f"[{i}/{len(samplefiles)}] Processing {dataset}")
        tstart = time.time()
        per_dataset_output[dataset] = run_one_dataset(processor_instance, dataset, info["files"], args.workers)
        print(f"  {time.time() - tstart:.1f}s")

    result = group_and_normalize(per_dataset_output, xsec)
    save(result, output_path)
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
