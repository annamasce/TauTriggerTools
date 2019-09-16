import FWCore.ParameterSet.Config as cms

print "Running on MC"


HLTLIST = cms.VPSet(
#    cms.PSet (
#        HLT = cms.string("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v"),
#        path1 = cms.vstring ("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09", "hltOverlapFilterIsoMu17LooseIsoPFTau20"),
#        path2 = cms.vstring ("hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterIsoMu17LooseIsoPFTau20"),
#        leg1 = cms.int32(13),
#        leg2 = cms.int32(15)
#    )

    #Mu-Tau20 (VBF monitoring)
    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_LooseChargedIsoPFTau20_Trk1_eta2p1_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu27LooseChargedIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackLooseChargedIsoAgainstMuon", "hltOverlapFilterIsoMu27LooseChargedIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_MediumChargedIsoPFTau20_Trk1_eta2p1_SingleL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu27MediumChargedIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackMediumChargedIsoAgainstMuon", "hltOverlapFilterIsoMu27MediumChargedIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu27_TightChargedIsoPFTau20_Trk1_eta2p1_SingleL_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu27TightChargedIsoPFTau20"),
        path2 = cms.vstring ("hltPFTau20TrackTightChargedIsoAgainstMuon", "hltOverlapFilterIsoMu27TightChargedIsoPFTau20"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    #Mu-Tau35 (di-tau monitoring)
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg"),
        path2 = cms.vstring ("hltSelectedPFTau35TrackPt1MediumChargedIsolationL1HLTMatchedReg", "hltOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24MediumChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
        path2 = cms.vstring ("hltSelectedPFTau35TrackPt1MediumChargedIsolationAndTightOOSCPhotonsL1HLTMatchedReg", "hltOverlapFilterIsoMu24MediumChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24TightChargedIsoPFTau35MonitoringReg"),
        path2 = cms.vstring ("hltSelectedPFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg", "hltOverlapFilterIsoMu24TightChargedIsoPFTau35MonitoringReg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f24QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu24TightChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
        path2 = cms.vstring ("hltSelectedPFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsL1HLTMatchedReg", "hltOverlapFilterIsoMu24TightChargedIsoAndTightOOSCPhotonsPFTau35MonitoringReg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    #Mu-Tau50 (Tau+MET monitoring)
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu22erIsoTau40erL1f0L2f10QL3f24QL3trkIsoFiltered0p07"),
        path2 = cms.vstring ("hltSelectedPFTau50MediumChargedIsolationL1HLTMatchedMu22IsoTau40"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),

    #Mu-Tau27 (signal path)
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched", "hltOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27MediumChargedIsolationAgainstMuonL1HLTMatched", "hltOverlapFilterIsoMu20MediumChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27TightChargedIsolationAgainstMuonL1HLTMatched", "hltOverlapFilterIsoMu20TightChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27LooseChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched", "hltOverlapFilterIsoMu20LooseChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27MediumChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched", "hltOverlapFilterIsoMu20MediumChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        path2 = cms.vstring ("hltSelectedPFTau27TightChargedIsolationTightOOSCPhotonsAgainstMuonL1HLTMatched", "hltOverlapFilterIsoMu20TightChargedIsoTightOOSCPhotonsPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),

    #Mu+Tau HPS
    cms.PSet (
        HLT = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltHpsOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        path2 = cms.vstring ("hltHpsSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched", "hltHpsOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),    
    cms.PSet (
        HLT = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v"),
        path1 = cms.vstring ("hltL3crIsoL1sBigOrMuXXerIsoTauYYerL1f0L2f10QL3f20QL3trkIsoFiltered0p07", "hltHpsOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg"),
        path2 = cms.vstring ("hltHpsSelectedPFTau35TrackPt1MediumChargedIsolationL1HLTMatchedReg", "hltHpsOverlapFilterIsoMu24MediumChargedIsoPFTau35MonitoringReg"),
        leg1 = cms.int32(13),
        leg2 = cms.int32(15)
    ),

    #SingleTau
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),

    #DoubleTau  
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau35TrackPt1TightChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau40TrackPt1MediumChargedIsolationAndTightOOSCPhotonsDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau40TrackPt1TightChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    ),
    cms.PSet (
        HLT = cms.string("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v"),
        path1 = cms.vstring ("hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        path2 = cms.vstring ("hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        leg1 = cms.int32(15),
        leg2 = cms.int32(15)
    )
)

# filter HLT paths for T&P
#import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
#hltFilter = hlt.hltHighLevel.clone(
#    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
#    #TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
#    #TriggerResultsTag = cms.InputTag("TriggerResults","","RECO"),
#    #TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
#    #TriggerResultsTag = cms.InputTag("TriggerResults","","HLT2"),
#    HLTPaths = ['HLT_IsoMu18_v5'],
#    #HLTPaths = ['HLT_Mu7p5_L2Mu2_Jpsi_v4'],
#    #HLTPaths = ['HLT_IsoMu18_v3'],
#    #HLTPaths = ['HLT_IsoMu18_v3'],
#    andOr = cms.bool(True), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
#    throw = cms.bool(True) #if True: throws exception if a trigger path is invalid
#)

## good taus - apply analysis selection
goodTaus = cms.EDFilter("PATTauRefSelector",
        src = cms.InputTag("slimmedTausNewID"),
        cut = cms.string(
                #'pt > 18 && abs(eta) < 2.5 ' #kinematics
                'pt > 20 && abs(eta) < 2.1 ' #kinematics
                '&& abs(charge) > 0 && abs(charge) < 2 ' #sometimes 2 prongs have charge != 1
                '&& tauID("decayModeFinding") > 0.5 ' # tau ID
                #'&& tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 2.5 ' # tau iso - NOTE: can as well use boolean discriminators with WP
                #'&& tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 2.5 ' # tau iso - NOTE: can as well use boolean discriminators with WP
                #'&& tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 1.0 ' # tau iso - NOTE: can as well use boolean discriminators with WP
                '&& tauID("againstMuonTight3") > 0.5 ' # anti Muon tight
                '&& tauID("againstElectronVLooseMVA6") > 0.5 ' # anti-Ele loose
        ),
        filter = cms.bool(True)
)

genMatchedTaus = cms.EDFilter("genMatchTauFilter",
        taus = cms.InputTag("goodTaus")
    )

# Ntuplizer.taus = cms.InputTag("genMatchedTaus")
Ntuplizer_noTagAndProbe_multipleTaus = cms.EDAnalyzer("Ntuplizer_noTagAndProbe_multipleTaus",
    treeName = cms.string("TagAndProbe"),
    genCollection = cms.InputTag("generator"),
    taus  = cms.InputTag("genMatchedTaus"),
    triggerSet = cms.InputTag("slimmedPatTrigger"),
    triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
    storeAllTriggers = cms.bool(True),
    #triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
    #triggerResultsLabel = cms.InputTag("TriggerResults", "", "RECO"),
    #triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT2"),
    #triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
    #L1Tau = cms.InputTag("caloStage2Digis", "Tau", "HLT"),
    #L1Tau = cms.InputTag("caloStage2Digis", "Tau", "HLT2"),
    L1Tau = cms.InputTag("caloStage2Digis", "Tau", "RECO"),
    L1EmuTau = cms.InputTag("simCaloStage2Digis", "MP"),
    jetCollection = cms.InputTag("slimmedJets"),
    l1tJetCollection = cms.InputTag("caloStage2Digis","Jet"),
    Vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
    triggerList = HLTLIST,
    L2CaloJet_ForIsoPix_Collection = cms.InputTag("hltL2TausForPixelIsolation", "", "TEST"),
    L2CaloJet_ForIsoPix_IsoCollection = cms.InputTag("hltL2TauPixelIsoTagProducer", "", "TEST")   
)

TAndPseq = cms.Sequence(
    #hltFilter      +
    #goodMuons      +
    goodTaus       +
    genMatchedTaus 
)

NtupleSeq = cms.Sequence(
    Ntuplizer_noTagAndProbe_multipleTaus
)
