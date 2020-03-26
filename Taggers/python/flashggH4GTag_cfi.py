import FWCore.ParameterSet.Config as cms

from flashgg.Taggers.globalVariables_cff import globalVariables
from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag # should include jet systematics
from flashgg.MicroAOD.flashggJets_cfi import  maxJetCollections, flashggDeepCSV

# # to not run jets systematics
# UnpackedJetCollectionVInputTag =  cms.VInputTag(
#     cms.InputTag("flashggUnpackedJets","0"), cms.InputTag("flashggUnpackedJets","1"), cms.InputTag("flashggUnpackedJets","2"), cms.InputTag("flashggUnpackedJets","3"), cms.InputTag("flashggUnpackedJets","4"),
#     cms.InputTag("flashggUnpackedJets","5"), cms.InputTag("flashggUnpackedJets","6"), cms.InputTag("flashggUnpackedJets","7"), cms.InputTag("flashggUnpackedJets","8"), cms.InputTag("flashggUnpackedJets","9"),
#     cms.InputTag("flashggUnpackedJets","10"), cms.InputTag("flashggUnpackedJets","11")
# )


# cfi = configuration fragment include
# Clone these params into _cfg
flashggH4GTag = cms.EDProducer("FlashggH4GTagProducer",
                                    globalVariables=globalVariables,
                                    PhotonTag = cms.InputTag('flashggRandomizedPhotons'),
                                    # DiPhotonTag = cms.InputTag('flashggDiPhotonSystematics'),
                                    DiPhotonTag = cms.InputTag('flashggPreselectedDiPhotons'),
                                    # DiPhotonTag = cms.InputTag('flashggDiPhotons'),
                                    SystLabel = cms.string(""),
                                    JetsName = cms.string("bRegProducer"), #
                                    JetsCollSize = cms.uint32(maxJetCollections), #
                                    JetsSuffixes = cms.vstring(''), #nominal and systematic variations
                                    DiPhotonName = cms.string('flashggPreselectedDiPhotons'), #
                                    VertexTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                    GenParticleTag         = cms.InputTag('flashggPrunedGenParticles'),
                                    ElectronTag            = cms.InputTag('flashggSelectedElectrons'),
                                    MuonTag                = cms.InputTag('flashggSelectedMuons'),
                                    METTag                 = cms.InputTag('flashggMets'),
                                    # METTag                 = cms.InputTag('flashggMetsCorr'), # RunIIFall17-3-2-0 contains these and NOT flashggMets
                                    JetTags                = UnpackedJetCollectionVInputTag,
                                    DiPhotonSuffixes = cms.vstring(''), #nominal and systematic variations
                                    useVertex0only=cms.bool(False),
                                    MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                    leptonPtThreshold = cms.double(10.),
                                    muonEtaThreshold = cms.double(2.4),
                                    leadPhoOverMassThreshold = cms.double(0.35), # was 0.375
                                    subleadPhoOverMassThreshold = cms.double(0.25),
                                    MVAThreshold = cms.double(0.0),
                                    deltaRMuonPhoThreshold = cms.double(0.4),
                                    jetsNumberThreshold = cms.double(99.), # originially 3.
                                    jetPtThreshold = cms.double(25.),
                                    jetEtaThreshold= cms.double(2.4),
                                    deltaRPhoLeadJet = cms.double(0.4),
                                    deltaRPhoSubLeadJet = cms.double(0.4),
                                    muPFIsoSumRelThreshold = cms.double(0.15),
                                    deltaRJetMuonThreshold = cms.double(0.4),
                                    PuIDCutoffThreshold = cms.double(0.8),
                                    PhoMVAThreshold = cms.double(-0.9),
                                    METThreshold = cms.double(0.),
                                    DeltaRTrkElec = cms.double(.4),
                                    TransverseImpactParam = cms.double(0.02),
                                    LongitudinalImpactParam = cms.double(0.2),
                                    deltaRPhoElectronThreshold = cms.double(0.4), # was 1
                                    deltaMassElectronZThreshold = cms.double(10.),
                                    electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                    nonTrigMVAThresholds = cms.vdouble(0.913286,0.805013,0.358969),
                                    nonTrigMVAEtaCuts = cms.vdouble(0.8,1.479,2.5),
                                    electronIsoThreshold = cms.double(0.15),
                                    electronNumOfHitsThreshold = cms.double(1),
                                    useElectronMVARecipe = cms.bool(False),
                                    useElectronLooseID = cms.bool(True),
                                    rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                    beamSpotTag            = cms.InputTag('offlineBeamSpot'),
                                    RECOfilters = cms.InputTag('TriggerResults::RECO'),
                                    PATfilters = cms.InputTag('TriggerResults::PAT'),
                                    FLASHfilters = cms.InputTag('TriggerResults::FLASHggMicroAOD'),
                                    # bTag = cms.string(flashggDeepCSV),
                                    # btagThresh = cms.double(100)     # no btag (Save all btags < 100)
                                    btagThresh = cms.double(0.45),
                                    doHHWWggTagCutFlowAnalysis = cms.bool(False), # save events for cut flow analysis,
                                    idSelection     = cms.PSet(),
                                    VertexSelectorName     = cms.string("FlashggLegacyVertexSelector"),
                                    VertexCandidateMapTag  = cms.InputTag("flashggVertexMapUnique"),
                                    useSingleLeg               = cms.bool(True),
                                    conversionTag          = cms.InputTag('reducedEgamma', 'reducedConversions'),
                                    conversionTagSingleLeg = cms.InputTag("reducedEgamma","reducedSingleLegConversions"),
                                    vertexIdMVAweightfile      = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_SL_2016.xml"),
                                    vertexProbMVAweightfile    = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxProb_SL_2016.xml"),
                                    nVtxSaveInfo            = cms.untracked.uint32(99),
                                    trackHighPurity         = cms.bool(False),
                                    pureGeomConvMatching    = cms.bool(True),
                                    dRexclude               = cms.double(0.05),
                                    #new reso:
                                    sigma1Pix               = cms.double(0.0125255),
                                    sigma1Tib               = cms.double(0.716301),
                                    sigma1Tob               = cms.double(3.17615),
                                    sigma1PixFwd            = cms.double(0.0581667),
                                    sigma1Tid               = cms.double(0.38521),
                                    sigma1Tec               = cms.double(1.67937),
                                    sigma2Pix               = cms.double(0.0298574),
                                    sigma2Tib               = cms.double(0.414393),
                                    sigma2Tob               = cms.double(1.06805),
                                    sigma2PixFwd            = cms.double(0.180419),
                                    sigma2Tid               = cms.double(0.494722),
                                    sigma2Tec               = cms.double(1.21941),
                                    singlelegsigma1Pix      = cms.double(0.0178107),
                                    singlelegsigma1Tib      = cms.double(1.3188),
                                    singlelegsigma1Tob      = cms.double(2.23662),
                                    singlelegsigma1PixFwd   = cms.double(0.152157),
                                    singlelegsigma1Tid      = cms.double(0.702755),
                                    singlelegsigma1Tec      = cms.double(2.46599),
                                    singlelegsigma2Pix      = cms.double(0.0935307),
                                    singlelegsigma2Tib      = cms.double(0.756568),
                                    singlelegsigma2Tob      = cms.double(0.62143),
                                    singlelegsigma2PixFwd   = cms.double(0.577081),
                                    singlelegsigma2Tid      = cms.double(0.892751),
                                    singlelegsigma2Tec      = cms.double(1.56638)
                                    )
# flashggHHWWggTagSequence = cms.Sequence( flashggHHWWggTag ) # not used
