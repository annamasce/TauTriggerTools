import FWCore.ParameterSet.Config as cms

from RecoTauTag.RecoTau.PFRecoTauDiscriminationByIsolation_cfi import *
from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts
from RecoTauTag.RecoTau.PFRecoTauDiscriminationByHPSSelection_cfi import hpsSelectionDiscriminator

from RecoTauTag.RecoTau.PFTauPrimaryVertexProducer_cfi      import *
from RecoTauTag.RecoTau.PFTauSecondaryVertexProducer_cfi    import *
from RecoTauTag.RecoTau.PFTauTransverseImpactParameters_cfi import *

from RecoTauTag.RecoTau.DeepTau_cfi import *

from RecoTauTag.RecoTau.PFRecoTauPFJetInputs_cfi import PFRecoTauPFJetInputs
## DeltaBeta correction factor
_ak4dBetaCorrection = 0.20

def add_reg_tag(name, useReg):
    if useReg:
        return name + "Reg"
    return name

def add_deepTau_sequence(process, working_points, useReg):


    # New L1 filter corresponding to BigOR filters in tau paths
    process.hltL1sTauDeepTauOR = process.hltL1sTauVeryBigOR.clone(
        L1SeedsLogicalExpression = cms.string("L1_DoubleIsoTau32er2p1 OR L1_DoubleIsoTau34er2p1 OR L1_DoubleIsoTau36er2p1 OR L1_DoubleTau70er2p1 OR L1_IsoTau40er2p1_ETMHF100 OR L1_IsoTau40er2p1_ETMHF110 OR L1_IsoTau40er2p1_ETMHF80 OR L1_IsoTau40er2p1_ETMHF90 OR L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3 OR L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3 OR L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3 OR L1_Mu18er2p1_Tau24er2p1 OR L1_Mu18er2p1_Tau26er2p1 OR L1_SingleTau120er2p1 OR L1_SingleTau130er2p1"),
    )

    process.hltHpsL1JetsHLTForDeepTauInput = process.hltHpsL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationMatchReg.clone(
        L1TauTrigger = cms.InputTag( "hltL1sTauDeepTauOR" ),
        JetSrc = cms.InputTag(add_reg_tag('hltHpsPFTauProducer', useReg)),
        ReduceTauContent = cms.bool(False),
        KeepOriginalVertex = cms.bool(True),
    )

    process.hltHpsPFTauDiscriminationByDecayModeFindingNewDMsForDeepTau = process.hltHpsPFTauDiscriminationByDecayModeFindingNewDMs.clone(
        PFTauProducer = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" )
    )

    process.hpsPFTauPrimaryVertexProducerForDeepTau = PFTauPrimaryVertexProducer.clone(
        # PFTauTag = add_reg_tag("hltHpsPFTauProducer", useReg),
        PFTauTag = "hltHpsL1JetsHLTForDeepTauInput",
        ElectronTag = "hltEgammaCandidates",
        MuonTag = add_reg_tag("hltMuons", useReg),
        PVTag = "hltPixelVertices",
        beamSpot = "hltOnlineBeamSpot",
        discriminators = [
            cms.PSet(
                discriminator = cms.InputTag('hltHpsPFTauDiscriminationByDecayModeFindingNewDMsForDeepTau'),
                selectionCut = cms.double(0.5)
            )
        ],
        cut = "pt > 18.0 & abs(eta) < 2.4",
        qualityCuts = PFTauQualityCuts
    )

    process.hpsPFTauSecondaryVertexProducerForDeepTau = PFTauSecondaryVertexProducer.clone(
        # PFTauTag = add_reg_tag("hltHpsPFTauProducer", useReg),
        PFTauTag = "hltHpsL1JetsHLTForDeepTauInput",
    )
    process.hpsPFTauTransverseImpactParametersForDeepTau = PFTauTransverseImpactParameters.clone(
        # PFTauTag = add_reg_tag("hltHpsPFTauProducer", useReg),
        PFTauTag = "hltHpsL1JetsHLTForDeepTauInput",
        PFTauPVATag = "hpsPFTauPrimaryVertexProducerForDeepTau",
        PFTauSVATag = "hpsPFTauSecondaryVertexProducerForDeepTau",
        useFullCalculation = True
    )

    process.hltFixedGridRhoFastjetAllTau = cms.EDProducer( "FixedGridRhoProducerFastjet",
        gridSpacing = cms.double( 0.55 ),
        maxRapidity = cms.double( 5.0 ),
        # pfCandidatesTag = cms.InputTag( "hltParticleFlowForTausReg" )
        pfCandidatesTag = cms.InputTag(add_reg_tag("hltParticleFlowForTaus", useReg))
    )

    chargedIsolationQualityCuts = PFTauQualityCuts.clone(
        isolationQualityCuts = cms.PSet( 
            maxDeltaZ = cms.double( 0.2 ),
            minTrackPt = cms.double( 0.5 ),
            minGammaEt = cms.double( 0.5 ),
            minTrackHits = cms.uint32( 3 ),
            minTrackPixelHits = cms.uint32( 0 ),
            maxTrackChi2 = cms.double( 100.0 ),
            maxTransverseImpactParameter = cms.double( 0.1 ),
            useTracksInsteadOfPFHadrons = cms.bool( False )
        ),
        primaryVertexSrc = cms.InputTag( "hltPixelVertices" ),
        signalQualityCuts = cms.PSet( 
            maxDeltaZ = cms.double( 0.2 ),
            minTrackPt = cms.double( 0.0 ),
            minGammaEt = cms.double( 0.5 ),
            minTrackHits = cms.uint32( 3 ),
            minTrackPixelHits = cms.uint32( 0 ),
            maxTrackChi2 = cms.double( 1000.0 ),
            maxTransverseImpactParameter = cms.double( 0.2 ),
            useTracksInsteadOfPFHadrons = cms.bool( False ),
            minNeutralHadronEt = cms.double( 1.0 )
        ),
        vxAssocQualityCuts = cms.PSet( 
            minTrackPt = cms.double( 0.0 ),
            minGammaEt = cms.double( 0.5 ),
            minTrackHits = cms.uint32( 3 ),
            minTrackPixelHits = cms.uint32( 0 ),
            maxTrackChi2 = cms.double( 1000.0 ),
            maxTransverseImpactParameter = cms.double( 0.2 ),
            useTracksInsteadOfPFHadrons = cms.bool( False )
        ),
    )

    PFTauQualityCuts.primaryVertexSrc = cms.InputTag("hltPixelVertices")

    ## Decay mode prediscriminant
    requireDecayMode = cms.PSet(
        BooleanOperator = cms.string("and"),
        decayMode = cms.PSet(
            Producer = cms.InputTag('hltHpsPFTauDiscriminationByDecayModeFindingNewDMsForDeepTau'),
            cut = cms.double(0.5)
        )
    )

    ## Cut based isolations dR=0.5
    process.hpsPFTauBasicDiscriminatorsForDeepTau = pfRecoTauDiscriminationByIsolation.clone(
        PFTauProducer = 'hltHpsL1JetsHLTForDeepTauInput',
        Prediscriminants = cms.PSet(  BooleanOperator = cms.string( "and" ) ),
        deltaBetaPUTrackPtCutOverride     = True, # Set the boolean = True to override.
        deltaBetaPUTrackPtCutOverride_val = 0.5,  # Set the value for new value.
        particleFlowSrc = add_reg_tag("hltParticleFlowForTaus", useReg),
        vertexSrc = PFTauQualityCuts.primaryVertexSrc,
        customOuterCone = PFRecoTauPFJetInputs.isolationConeSize,
        isoConeSizeForDeltaBeta = 0.8,
        deltaBetaFactor = "%0.4f"%(_ak4dBetaCorrection),
        qualityCuts = dict(isolationQualityCuts = dict(minTrackHits = 3, minGammaEt = 1.0, minTrackPt = 0.5)),
        IDdefinitions = [
            cms.PSet(
                IDname = cms.string("ChargedIsoPtSum"),
                ApplyDiscriminationByTrackerIsolation = cms.bool(True),
                storeRawSumPt = cms.bool(True)
            ),
            cms.PSet(
                IDname = cms.string("NeutralIsoPtSum"),
                ApplyDiscriminationByECALIsolation = cms.bool(True),
                storeRawSumPt = cms.bool(True)
            ),
            cms.PSet(
                IDname = cms.string("NeutralIsoPtSumWeight"),
                ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
                storeRawSumPt = cms.bool(True),
                UseAllPFCandsForWeights = cms.bool(True)
            ),
            cms.PSet(
                IDname = cms.string("TauFootprintCorrection"),
                storeRawFootprintCorrection = cms.bool(True)
            ),
            cms.PSet(
                IDname = cms.string("PhotonPtSumOutsideSignalCone"),
                storeRawPhotonSumPt_outsideSignalCone = cms.bool(True)
            ),
            cms.PSet(
                IDname = cms.string("PUcorrPtSum"),
                applyDeltaBetaCorrection = cms.bool(True),
                storeRawPUsumPt = cms.bool(True)
            ),
        ],
    )

    ## Cut based isolations dR=0.3
    process.hpsPFTauBasicDiscriminatorsdR03ForDeepTau = process.hpsPFTauBasicDiscriminatorsForDeepTau.clone(
        customOuterCone = 0.3
    )

    file_names = [
                    'core:RecoTauTag/TrainingFiles/data/DeepTauId/deepTau_2017v2p6_e6_core.pb',
                    'inner:RecoTauTag/TrainingFiles/data/DeepTauId/deepTau_2017v2p6_e6_inner.pb',
                    'outer:RecoTauTag/TrainingFiles/data/DeepTauId/deepTau_2017v2p6_e6_outer.pb',
                ]

    process.deepTauProducer = DeepTau.clone(
        taus = 'hltHpsL1JetsHLTForDeepTauInput',
        pfcands = add_reg_tag('hltParticleFlowForTaus', useReg),
        vertices = 'hltPixelVertices',
        rho = 'hltFixedGridRhoFastjetAllTau',
        graph_file = file_names,
        disable_dxy_pca = cms.bool(True),
        is_online = cms.bool(True),
        pfTauTransverseImpactParameters = 'hpsPFTauTransverseImpactParametersForDeepTau',
        basicTauDiscriminators = 'hpsPFTauBasicDiscriminatorsForDeepTau',
        basicTauDiscriminatorsdR03 = 'hpsPFTauBasicDiscriminatorsdR03ForDeepTau',
        Prediscriminants = cms.PSet(  BooleanOperator = cms.string( "and" ) ),  
        VSeWP = working_points,
        VSmuWP = working_points,
        VSjetWP = working_points     
    )	

    # Add DeepTauProducer
    process.HLTHPSDeepTauIsoPFTauSequence = cms.Sequence(process.hltL1sTauDeepTauOR + process.hltHpsL1JetsHLTForDeepTauInput + process.hltHpsPFTauDiscriminationByDecayModeFindingNewDMsForDeepTau + process.hpsPFTauPrimaryVertexProducerForDeepTau + process.hpsPFTauSecondaryVertexProducerForDeepTau + process.hpsPFTauTransverseImpactParametersForDeepTau + process.hltFixedGridRhoFastjetAllTau + process.hpsPFTauBasicDiscriminatorsForDeepTau + process.hpsPFTauBasicDiscriminatorsdR03ForDeepTau + process.deepTauProducer)
    return process

def customiseDiTauForDeepTau(process, useReg, working_points, addCounters, isMC):

    ## Gen counter
    process.genCounter = cms.EDFilter( "CounterFilter",
        isMC = isMC, #from outside
        store_hist = cms.bool(False), #from outside
        store_both = cms.bool(True), #from outside
        store_MET = cms.bool(False),
        position = cms.string("gen"),
        use_deepTau = cms.bool(False),
        deepTauVSe = cms.InputTag('try1'),
        deepTauVSmu = cms.InputTag('try2'),
        deepTauVSjet = cms.InputTag('try3'),
        original_taus = cms.InputTag('try4'),
        taus = cms.InputTag('try5'),
        track_taus = cms.InputTag('try6'),
        MET = cms.InputTag('try7'),
        puInfo = cms.InputTag('try8'),
        vertices = cms.InputTag('try9'),
        genParticles = cms.InputTag('genParticles')
    )

    if addCounters:
        process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.insert(1, process.genCounter)

    ## Final counter
    process.jetsFilterDiTau = process.jetsFilter.clone(
        position = cms.string("final_DiTau"),
        taus = cms.InputTag("hltHpsDoublePFTau35TrackPt1DeepTau35IsolationDz02"),
        track_taus = cms.InputTag("hltHpsL1JetsHLTForDeepTauInput")        
    )
    
    process.hltHpsSelectedPFTausTrackPt1DeepTau35Isolation = process.hltHpsSelectedPFTausTrackPt1MediumChargedIsolationReg.clone(
        src = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
        discriminators = [
            # cms.PSet(  
            #     discriminator = cms.InputTag( "hltHpsPFTauTrackPt1Discriminator" ),
            #     selectionCut = cms.double( 0.5 )
            # )
        ],
        discriminatorContainers = [
            cms.PSet(  
                discriminator = cms.InputTag( "deepTauProducer", "VSjet" ),
                rawValues = cms.vstring(),
                selectionCuts = cms.vdouble(),
                workingPoints = cms.vstring(working_points),
            )
        ]
    )

    process.hltHpsDoublePFTau35TrackPt1DeepTau35Isolation = process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationReg.clone(
        inputTag = "hltHpsSelectedPFTausTrackPt1DeepTau35Isolation",
    )

    process.hltHpsL1JetsHLTDoublePFTauTrackPt1DeepTauMatch = process.hltHpsL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationMatchReg.clone(
        JetSrc = "hltHpsSelectedPFTausTrackPt1DeepTau35Isolation",
    )

    process.hltHpsDoublePFTau35TrackPt1DeepTauL1HLTMatched = process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationL1HLTMatchedReg.clone(
        inputTag = "hltHpsL1JetsHLTDoublePFTauTrackPt1DeepTauMatch",
    )

    process.hltHpsDoublePFTau35TrackPt1DeepTau35IsolationDz02 = process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg.clone(
        JetSrc = "hltHpsL1JetsHLTDoublePFTauTrackPt1DeepTauMatch"
    )

    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.remove(process.HLTHPSMediumChargedIsoPFTauSequenceReg)
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.remove(process.hltHpsSelectedPFTausTrackPt1MediumChargedIsolationReg)
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.remove(process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationReg)
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.remove(process.hltHpsL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationMatchReg)
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.remove(process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationL1HLTMatchedReg)
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.remove(process.hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg)
    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.remove(process.HLTEndSequence)

    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4 += (process.HLTHPSDeepTauIsoPFTauSequence + process.hltHpsSelectedPFTausTrackPt1DeepTau35Isolation + process.hltHpsDoublePFTau35TrackPt1DeepTau35Isolation + process.hltHpsL1JetsHLTDoublePFTauTrackPt1DeepTauMatch + process.hltHpsDoublePFTau35TrackPt1DeepTauL1HLTMatched + process.hltHpsDoublePFTau35TrackPt1DeepTau35IsolationDz02 + process.HLTEndSequence)

    if addCounters:
        process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.insert(-1, process.jetsFilterDiTau)

    return process

def customiseEleTauForDeepTau(process, working_points, addCounters):

    ## Final counter
    process.jetsFilterEleTau = process.jetsFilter.clone(
        position = cms.string("final_EleTau"),
        taus = cms.InputTag("hltHpsOverlapFilterIsoEle24WPTightGsfDeepTauPFTau30"),
        track_taus = cms.InputTag("hltHpsSelectedPFTausTrackFindingForDeepTau")
    )

    process.hltHpsPFTau30ForDeepTau = process.hltHpsPFTau30.clone(
        inputTag = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
    )

    process.hltHpsSelectedPFTausTrackFindingDeepTau = process.hltHpsSelectedPFTausTrackFindingLooseChargedIsolation.clone(
        src = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
        discriminators = cms.VPSet( 
            # cms.PSet(  discriminator = cms.InputTag( "hltHpsPFTauTrackFindingDiscriminator" ),
            #             selectionCut = cms.double( 0.5 )
            #  ),
            ),
        discriminatorContainers = cms.VPSet( 
            cms.PSet(  
                discriminator = cms.InputTag( "deepTauProducer", "VSjet" ),
                rawValues = cms.vstring(),
                selectionCuts = cms.vdouble(),
                workingPoints = cms.vstring(working_points),
            )
        )
    )

    process.hltHpsPFTau30TrackDeepTau = process.hltHpsPFTau30TrackLooseChargedIso.clone(
        inputTag = cms.InputTag( "hltHpsSelectedPFTausTrackFindingDeepTau" ),
    )

    process.hltHpsL1JetsHLTPFTauTrackDeepTauMatch = process.hltHpsL1JetsHLTPFTauTrackLooseChargedIsolationMatch.clone(
        JetSrc = cms.InputTag( "hltHpsSelectedPFTausTrackFindingDeepTau" ),
    )

    process.hltHpsSelectedPFTau30DeepTauL1HLTMatched = process.hltHpsSelectedPFTau30LooseChargedIsolationL1HLTMatched.clone(
        inputTag = cms.InputTag( "hltHpsL1JetsHLTPFTauTrackDeepTauMatch" ),
    )

    process.hltHpsOverlapFilterIsoEle24WPTightGsfDeepTauPFTau30 = process.hltHpsOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30.clone(
        originTag2 = cms.VInputTag( "hltHpsSelectedPFTausTrackFindingDeepTau" ),
        inputTag2 = cms.InputTag( "hltHpsSelectedPFTau30DeepTauL1HLTMatched" ),
    )

    process.HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v1.remove(process.HLTHPSLooseChargedIsoPFTau30Sequence)
    process.HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v1.remove(process.hltHpsL1JetsHLTPFTauTrackLooseChargedIsolationMatch)
    process.HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v1.remove(process.hltHpsSelectedPFTau30LooseChargedIsolationL1HLTMatched)
    process.HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v1.remove(process.hltHpsOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30)
    process.HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v1.remove(process.HLTEndSequence)

    process.HLTHPSDeepTauPFTau30Sequence = cms.Sequence(process.HLTHPSDeepTauIsoPFTauSequence + process.hltHpsPFTau30ForDeepTau + process.hltHpsPFTauTrackFindingDiscriminatorForDeepTau + cms.ignore(process.hltHpsSelectedPFTausTrackFindingForDeepTau) + process.hltHpsSelectedPFTausTrackFindingDeepTau + process.hltHpsPFTau30TrackDeepTau)
    process.HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v1 += (process.HLTHPSDeepTauPFTau30Sequence + process.hltHpsL1JetsHLTPFTauTrackDeepTauMatch + process.hltHpsSelectedPFTau30DeepTauL1HLTMatched + process.hltHpsOverlapFilterIsoEle24WPTightGsfDeepTauPFTau30 + process.HLTEndSequence)

    if addCounters:
        process.HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v1.insert(-1, process.jetsFilterEleTau)

    return process

def customiseMuTauForDeepTau(process, working_points, addCounters):

    ## Final counter
    process.jetsFilterMuTau = process.jetsFilter.clone(
        position = cms.string("final_MuTau"),
        taus = cms.InputTag("hltHpsOverlapFilterIsoMu20DeepTauPFTau27L1Seeded"),
        track_taus = cms.InputTag("hltHpsSelectedPFTausTrackFindingForDeepTau")
    )

    process.hltHpsPFTau27ForDeepTau = process.hltHpsPFTau27.clone(
         inputTag = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
    )

    process.hltHpsSelectedPFTausTrackFindingDeepTau = process.hltHpsSelectedPFTausTrackFindingLooseChargedIsolation.clone(
        src = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
        discriminators = cms.VPSet( 
            # cms.PSet(  discriminator = cms.InputTag( "hltHpsPFTauTrackFindingDiscriminator" ),
            #             selectionCut = cms.double( 0.5 )
            #  ),
            ),
        discriminatorContainers = [ 
            cms.PSet(  
                discriminator = cms.InputTag( "deepTauProducer", "VSjet" ),
                rawValues = cms.vstring(),
                selectionCuts = cms.vdouble(),
                workingPoints = cms.vstring(working_points),
            )
        ]
    )

    process.hltHpsPFTau27TrackDeepTau = process.hltHpsPFTau27TrackLooseChargedIso.clone(
        inputTag = cms.InputTag( "hltHpsSelectedPFTausTrackFindingDeepTau" ),
    )

    process.hltHpsPFTauDeepTauAgainstMuonDiscriminator = process.hltHpsPFTauAgainstMuonDiscriminator.clone(
        PFTauProducer = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
    )

    process.hltHpsSelectedPFTausTrackFindingDeepTauAgainstMuon = process.hltHpsSelectedPFTausTrackFindingLooseChargedIsolationAgainstMuon.clone(
        src = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
        discriminators = cms.VPSet( 
            # cms.PSet(  discriminator = cms.InputTag( "hltHpsPFTauTrackFindingDiscriminator" ),
            # selectionCut = cms.double( 0.5 )
            # ),
            cms.PSet(  discriminator = cms.InputTag( "hltHpsPFTauDeepTauAgainstMuonDiscriminator" ),
                        selectionCut = cms.double( 0.5 )
            )
        ),
        discriminatorContainers = [ 
            cms.PSet(  
                discriminator = cms.InputTag( "deepTauProducer", "VSjet" ),
                rawValues = cms.vstring(),
                selectionCuts = cms.vdouble(),
                workingPoints = cms.vstring(working_points),
            )
        ]
    )

    process.hltHpsPFTau27TrackDeepTauAgainstMuon = process.hltHpsPFTau27TrackLooseChargedIsoAgainstMuon.clone(
        inputTag = cms.InputTag( "hltHpsSelectedPFTausTrackFindingDeepTauAgainstMuon" ),
    )

    process.hltHpsL1JetsHLTPFTauTrackDeepTauAgainstMuonMatch = process.hltHpsL1JetsHLTPFTauTrackLooseChargedIsolationAgainstMuonMatch.clone(
        JetSrc = cms.InputTag( "hltHpsSelectedPFTausTrackFindingDeepTauAgainstMuon" ),
    )

    process.hltHpsSelectedPFTau27DeepTauAgainstMuonL1HLTMatched = process.hltHpsSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched.clone(
        inputTag = cms.InputTag( "hltHpsL1JetsHLTPFTauTrackDeepTauAgainstMuonMatch" ),
    )

    process.hltHpsOverlapFilterIsoMu20DeepTauPFTau27L1Seeded = process.hltHpsOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded.clone(
        originTag2 = cms.VInputTag( 'hltHpsSelectedPFTausTrackFindingDeepTauAgainstMuon' ),
        inputTag2 = cms.InputTag( "hltHpsSelectedPFTau27DeepTauAgainstMuonL1HLTMatched" ),
    )

    process.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v4.remove(process.HLTHPSLooseChargedIsoAntiMuonPFTau27Sequence)
    process.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v4.remove(process.hltHpsL1JetsHLTPFTauTrackLooseChargedIsolationAgainstMuonMatch)
    process.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v4.remove(process.hltHpsSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched)
    process.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v4.remove(process.hltHpsOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded)
    process.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v4.remove(process.HLTEndSequence)

    process.HLTHPSDeepTauAntiMuonPFTau27Sequence = cms.Sequence(process.HLTHPSDeepTauIsoPFTauSequence + process.hltHpsPFTau27ForDeepTau + process.hltHpsPFTauTrackFindingDiscriminatorForDeepTau + cms.ignore(process.hltHpsSelectedPFTausTrackFindingForDeepTau) + process.hltHpsSelectedPFTausTrackFindingDeepTau + process.hltHpsPFTau27TrackDeepTau + process.hltHpsPFTauDeepTauAgainstMuonDiscriminator + process.hltHpsSelectedPFTausTrackFindingDeepTauAgainstMuon + process.hltHpsPFTau27TrackDeepTauAgainstMuon)
    process.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v4 += (process.HLTHPSDeepTauAntiMuonPFTau27Sequence + process.hltHpsL1JetsHLTPFTauTrackDeepTauAgainstMuonMatch + process.hltHpsSelectedPFTau27DeepTauAgainstMuonL1HLTMatched + process.hltHpsOverlapFilterIsoMu20DeepTauPFTau27L1Seeded + process.HLTEndSequence)

    if addCounters:
        process.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v4.insert(-1, process.jetsFilterMuTau)

    return process

def customiseTauMETForDeepTau(process, working_points, addCounters):

    ## Final counter
    process.jetsFilterTauMET = process.jetsFilter.clone(
        store_MET = cms.bool(True),
        position = cms.string("final_TauMET"),
        taus = cms.InputTag("hltSelectedPFTau50DeepTauL1HLTMatched"),
        track_taus = cms.InputTag("hltSelectedPFTausTrackPt30"),
        MET = cms.InputTag("hltMet")
    )

    process.hltMET100.MinPt = cms.double( 80.0 )

    process.hltPFTau50ForDeepTau = process.hltPFTau50.clone(
        inputTag = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
    )

    process.hltPFTauTrackPt30DiscriminatorForDeepTau = process.hltPFTauTrackPt30Discriminator.clone(
        PFTauProducer = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" )
    )

    process.hltSelectedPFTausTrackPt30 = process.hltSelectedPFTausTrackPt30MediumAbsOrRelIsolation1Prong.clone(
        src = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
        discriminators = cms.VPSet( 
            cms.PSet(  discriminator = cms.InputTag( "hltPFTauTrackPt30DiscriminatorForDeepTau" ),
                        selectionCut = cms.double( 0.5 )
            ),
        ),
    )

    process.hltSelectedPFTausDeepTau = process.hltSelectedPFTausTrackPt30MediumAbsOrRelIsolation1Prong.clone(
        src = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
        discriminators = cms.VPSet( 
        ),
        discriminatorContainers = [ 
            cms.PSet(  
                discriminator = cms.InputTag( "deepTauProducer", "VSjet" ),
                rawValues = cms.vstring(),
                selectionCuts = cms.vdouble(),
                workingPoints = cms.vstring(working_points),
            )
        ]
    )

    process.hltPFTau50DeepTau = process.hltPFTau50TrackPt30MediumAbsOrRelIso1Prong.clone(
        inputTag = cms.InputTag( "hltSelectedPFTausDeepTau" ),
    )

    process.hltL1JetsHLTPFTau50DeepTauMatch = process.hltL1JetsHLTPFTau50Track30MediumChargedIsolationMatch.clone(
        JetSrc = cms.InputTag( "hltSelectedPFTausDeepTau" ),
    )

    process.hltSelectedPFTau50DeepTauL1HLTMatched = process.hltSelectedPFTau50MediumChargedIsolationL1HLTMatched.clone(
        inputTag = cms.InputTag( "hltL1JetsHLTPFTau50DeepTauMatch" ),
    )

    process.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v12.remove(process.HLTMediumChargedIsoPFTau50Trk30OnePrSequence)
    process.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v12.remove(process.hltL1JetsHLTPFTau50Track30MediumChargedIsolationMatch)
    process.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v12.remove(process.hltSelectedPFTau50MediumChargedIsolationL1HLTMatched)
    process.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v12.remove(process.HLTEndSequence)

    process.HLTDeepTauPFTau50Sequence = cms.Sequence(process.HLTHPSDeepTauIsoPFTauSequence + process.hltPFTau50ForDeepTau + process.hltPFTauTrackPt30DiscriminatorForDeepTau + cms.ignore(process.hltSelectedPFTausTrackPt30) + process.hltSelectedPFTausDeepTau + process.hltPFTau50DeepTau)
    process.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v12 += (process.HLTDeepTauPFTau50Sequence + process.hltL1JetsHLTPFTau50DeepTauMatch + process.hltSelectedPFTau50DeepTauL1HLTMatched + process.HLTEndSequence)
    
    if addCounters:
        process.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v12.insert(-1, process.jetsFilterTauMET)

    return process

def customiseHighPtTauForDeepTau(process, working_points, addCounters):

    ## Final counter
    process.jetsFilterHighPtTau = process.jetsFilter.clone(
        position = cms.string("final_HighPtTau"),
        taus = cms.InputTag("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched"),
        track_taus = cms.InputTag("hltSelectedPFTausTrackPt50")
    )

    process.hltPFTauTrackPt50DiscriminatorForDeepTau = process.hltPFTauTrackPt50Discriminator.clone(
        PFTauProducer = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" )
    )

    process.hltSelectedPFTausTrackPt50 = process.hltSelectedPFTausTrackPt50AbsOrRelMediumHighPtRelaxedIsoIsolation.clone(
        src = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
        discriminators = cms.VPSet( 
            cms.PSet(  discriminator = cms.InputTag( "hltPFTauTrackPt50DiscriminatorForDeepTau" ),
                        selectionCut = cms.double( 0.5 )
            ),
        ),
    )

    process.hltSelectedPFTausDeepTauHighPtRelaxedIsoIsolation = process.hltSelectedPFTausTrackPt50AbsOrRelMediumHighPtRelaxedIsoIsolation.clone(
        src = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
        discriminators = cms.VPSet( 
        ),
        discriminatorContainers = [ 
            cms.PSet(  
                discriminator = cms.InputTag( "deepTauProducer", "VSjet" ),
                rawValues = cms.vstring(),
                selectionCuts = cms.vdouble(),
                workingPoints = cms.vstring(working_points),
            )
        ]
    )

    process.hltPFTau180DeepTauHighPtRelaxedIso = process.hltPFTau180TrackPt50LooseAbsOrRelMediumHighPtRelaxedIsoIso.clone(
        inputTag = cms.InputTag( "hltSelectedPFTausDeepTauHighPtRelaxedIsoIsolation" ),
    )

    process.hltL1JetsHLTSinglePFTauDeepTauMatch = process.hltL1JetsHLTSinglePFTauTrackMediumChargedIsolationMatch.clone(
        JetSrc = cms.InputTag( "hltSelectedPFTausDeepTauHighPtRelaxedIsoIsolation" ),
    )

    process.hltSelectedPFTau180MediumChargedIsolationL1HLTMatched = process.hltSelectedPFTau180MediumChargedIsolationL1HLTMatched.clone(
        inputTag = cms.InputTag( "hltL1JetsHLTSinglePFTauDeepTauMatch" ),
    )

    process.HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v12.remove(process.hltPFTau180Track)
    process.HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v12.remove(process.HLTMediumChargedIsoPFTauHighPtRelaxedIsoTrk50)
    process.HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v12.remove(process.hltPFTau180TrackPt50LooseAbsOrRelMediumHighPtRelaxedIsoIso)
    process.HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v12.remove(process.hltL1JetsHLTSinglePFTauTrackMediumChargedIsolationMatch)
    process.HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v12.remove(process.hltSelectedPFTau180MediumChargedIsolationL1HLTMatched)
    process.HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v12.remove(process.HLTEndSequence)

    process.HLTDeepTauPFTauHighPtRelaxedIso = cms.Sequence(process.HLTHPSDeepTauIsoPFTauSequence + process.hltPFTauTrackPt50DiscriminatorForDeepTau + cms.ignore(process.hltSelectedPFTausTrackPt50) + process.hltSelectedPFTausDeepTauHighPtRelaxedIsoIsolation)
    process.HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v12 += (process.HLTDeepTauPFTauHighPtRelaxedIso + process.hltPFTau180DeepTauHighPtRelaxedIso + process.hltL1JetsHLTSinglePFTauDeepTauMatch + process.hltSelectedPFTau180MediumChargedIsolationL1HLTMatched + process.HLTEndSequence)
    
    if addCounters:
        process.HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v12.insert(-1, process.jetsFilterHighPtTau)

    return process

def update(process, useReg=True, resetWP=False, addCounters=False, isData=False):
    process.options.wantSummary = cms.untracked.bool(True)

    isMC = cms.bool(True)
    if isData:
        isMC = cms.bool(False)

    ## Final counter -> taus InputTag to be customised for each path
    process.jetsFilter = cms.EDFilter( "CounterFilter",
        isMC = isMC, #from outside
        store_hist = cms.bool(False), #from outside
        store_both = cms.bool(False), #from outside
        position = cms.string("final"),
        use_deepTau = cms.bool(True),
        store_MET = cms.bool(False),
        deepTauVSe = cms.InputTag('deepTauProducer', 'VSe'),
        deepTauVSmu = cms.InputTag('deepTauProducer', 'VSmu'),
        deepTauVSjet = cms.InputTag('deepTauProducer', 'VSjet'),
        original_taus = cms.InputTag('hltHpsL1JetsHLTForDeepTauInput'),
        taus = cms.InputTag(''),
        track_taus = cms.InputTag(''),
        MET = cms.InputTag(''),
        puInfo = cms.InputTag('addPileupInfo','','HLT'),
        vertices = cms.InputTag('hltPixelVertices'),
        genParticles = cms.InputTag('genParticles')
    )

    ## Define working points for deepTau
    def getLinExpression(x1, x2, y1, y2):
        return "(((({3}-{2})/({1}-{0}))*(pt-{0}))+{2})".format(x1, x2, y1, y2)

    val1, val2 = ("0.49948551", "0.125")
    working_points = ["{0}*(pt < 35)+".format(val1)+getLinExpression("35", "300", val1, val2)+ "*(35 <= pt && pt < 300) + {0}*(pt >= 300)".format(val2)]
    if resetWP:
        working_points = ["-1."]
    print("WP:", working_points)

    ## Add deepTau sequence
    process = add_deepTau_sequence(process, working_points=working_points, useReg=useReg)

    ## Add track-finding discriminator 
    process.hltHpsPFTauTrackFindingDiscriminatorForDeepTau = process.hltHpsPFTauTrackFindingDiscriminator.clone(
        PFTauProducer = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" )
    )
    process.hltHpsSelectedPFTausTrackFindingForDeepTau = process.hltHpsSelectedPFTausTrackFinding.clone(
        src = cms.InputTag( "hltHpsL1JetsHLTForDeepTauInput" ),
        discriminators = cms.VPSet( 
            cms.PSet(  discriminator = cms.InputTag( "hltHpsPFTauTrackFindingDiscriminatorForDeepTau" ),
                        selectionCut = cms.double( 0.5 )
            ),
        ),
    )

    ## Customise di-tau path
    process = customiseDiTauForDeepTau(process, useReg=useReg, working_points=working_points, addCounters=addCounters, isMC=isMC)

    ## Customise ele+tau path
    process = customiseEleTauForDeepTau(process, working_points=working_points, addCounters=addCounters)

    ## Customise mu+tau path
    process = customiseMuTauForDeepTau(process, working_points=working_points, addCounters=addCounters)

    ## Customise tau+MET path
    process = customiseTauMETForDeepTau(process, working_points=working_points, addCounters=addCounters)

    ## Customise High-pt path
    process = customiseHighPtTauForDeepTau(process, working_points=working_points, addCounters=addCounters)

    if addCounters:
        process.TFileService = cms.Service("TFileService", fileName = cms.string("histo.root"))

    return process

def update_oldHLT(process, isData=False):
    process.options.wantSummary = cms.untracked.bool(True)

    isMC = cms.bool(True)
    if isData:
        isMC = cms.bool(False)

    ## Gen counter
    process.genCounter = cms.EDFilter( "CounterFilter",
        isMC = isMC, #from outside
        store_hist = cms.bool(False), #from outside
        store_both = cms.bool(True), #from outside
        store_MET = cms.bool(False),
        position = cms.string("gen"),
        use_deepTau = cms.bool(False),
        deepTauVSe = cms.InputTag('try1'),
        deepTauVSmu = cms.InputTag('try2'),
        deepTauVSjet = cms.InputTag('try3'),
        original_taus = cms.InputTag('try4'),
        taus = cms.InputTag('try5'),
        track_taus = cms.InputTag('try6'),
        MET = cms.InputTag('try7'),
        puInfo = cms.InputTag('try8'),
        vertices = cms.InputTag('try9'),
        genParticles = cms.InputTag('genParticles')
    )

    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.insert(1, process.genCounter)

    ## Final counter -> taus InputTag to be customised for each path
    process.jetsFilter = cms.EDFilter( "CounterFilter",
        isMC = isMC, #from outside
        store_hist = cms.bool(False), #from outside
        store_both = cms.bool(False), #from outside
        position = cms.string("final"),
        use_deepTau = cms.bool(False),
        store_MET = cms.bool(False),
        deepTauVSe = cms.InputTag(""),
        deepTauVSmu = cms.InputTag(""),
        deepTauVSjet = cms.InputTag(""),
        original_taus = cms.InputTag(""),
        taus = cms.InputTag(''),
        track_taus = cms.InputTag(''),
        MET = cms.InputTag(''),
        puInfo = cms.InputTag('addPileupInfo','','HLT'),
        vertices = cms.InputTag('hltPixelVertices'),
        genParticles = cms.InputTag('genParticles')
    )

    ## Final counter Di-Tau
    process.jetsFilterDiTau = process.jetsFilter.clone(
        position = cms.string("final_DiTau"),
        original_taus = cms.InputTag("hltHpsL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationMatchReg"),
        taus = cms.InputTag("hltHpsDoublePFTau35TrackPt1MediumChargedIsolationDz02Reg"),
        track_taus = cms.InputTag("hltHpsL1JetsHLTDoublePFTauTrackPt1MediumChargedIsolationMatchReg")        
    )

    process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4.insert(-1, process.jetsFilterDiTau)

    ## Final counter EleTau
    process.jetsFilterEleTau = process.jetsFilter.clone(
        position = cms.string("final_EleTau"),
        original_taus = cms.InputTag("hltHpsL1JetsHLTPFTauTrackLooseChargedIsolationMatch"),
        taus = cms.InputTag("hltHpsOverlapFilterIsoEle24WPTightGsfLooseIsoPFTau30"),
        track_taus = cms.InputTag("hltHpsL1JetsHLTPFTauTrackLooseChargedIsolationMatch")
    )

    process.HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v1.insert(-1, process.jetsFilterEleTau)

    ## Final counter MuTau
    process.jetsFilterMuTau = process.jetsFilter.clone(
        position = cms.string("final_MuTau"),
        original_taus = cms.InputTag("hltHpsL1JetsHLTPFTauTrackLooseChargedIsolationAgainstMuonMatch"),
        taus = cms.InputTag("hltHpsOverlapFilterIsoMu20LooseChargedIsoPFTau27L1Seeded"),
        track_taus = cms.InputTag("hltHpsL1JetsHLTPFTauTrackLooseChargedIsolationAgainstMuonMatch")
    )

    process.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v4.insert(-1, process.jetsFilterMuTau)

    ## Final counter TauMET
    process.jetsFilterTauMET = process.jetsFilter.clone(
        store_MET = cms.bool(True),
        position = cms.string("final_TauMET"),
        original_taus = cms.InputTag("hltL1JetsHLTPFTau50Track30MediumChargedIsolationMatch"),
        taus = cms.InputTag("hltSelectedPFTau50MediumChargedIsolationL1HLTMatched"),
        track_taus = cms.InputTag("hltL1JetsHLTPFTau50Track30MediumChargedIsolationMatch"),
        MET = cms.InputTag("hltMet")
    )

    process.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v12.insert(-1, process.jetsFilterTauMET)

    ## Final counter High-pt Tau
    process.jetsFilterHighPtTau = process.jetsFilter.clone(
        position = cms.string("final_HighPtTau"),
        original_taus = cms.InputTag("hltL1JetsHLTSinglePFTauTrackMediumChargedIsolationMatch"),
        taus = cms.InputTag("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched"),
        track_taus = cms.InputTag("hltL1JetsHLTSinglePFTauTrackMediumChargedIsolationMatch")
    )

    process.HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v12.insert(-1, process.jetsFilterHighPtTau)

    process.TFileService = cms.Service("TFileService", fileName = cms.string("histo.root"))

    return process

