import FWCore.ParameterSet.Config as cms
import flashgg.Taggers.dumperConfigTools as cfgTools

from flashgg.Taggers.globalVariables_cff import globalVariables
from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag # should include jet systematics
from flashgg.MicroAOD.flashggJets_cfi import  maxJetCollections, flashggDeepCSV

vertexIdMVAweightfileH4G = ""

# cfi = configuration fragment include
# Clone these params into _cfg
flashggH4GTag = cms.EDProducer("FlashggH4GTagProducer",
                                    globalVariables=globalVariables,
                                    PhotonTag = cms.InputTag('flashggRandomizedPhotons'),
                                    # DiPhotonTag = cms.InputTag('flashggPreselectedDiPhotons'),
                                    #DiPhotonName = cms.string('flashggDiPhotonSystematics'),
                                    DiPhotonName = cms.string('flashggPreselectedDiPhotons'),
                                    idSelection     = cms.PSet(),
                                    SystLabel = cms.string(""),
                                    VertexTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                    GenParticleTag         = cms.InputTag('flashggPrunedGenParticles'),
                                    DiPhotonSuffixes = cms.vstring(''), #nominal and systematic variations
                                    useVertex0only=cms.bool(False),
                                    MVAResultTag=cms.InputTag('flashggDiPhotonMVA'),
                                    rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                    beamSpotTag            = cms.InputTag('offlineBeamSpot'),
                                    conversionTag          = cms.InputTag('reducedEgamma', 'reducedConversions'),
                                    conversionTagSingleLeg = cms.InputTag("reducedEgamma","reducedSingleLegConversions"),
                                    VertexSelectorName     = cms.string("FlashggLegacyVertexSelector"),
                                    VertexCandidateMapTag  = cms.InputTag("flashggVertexMapUnique"),
                                     ##Parameters for Legacy Vertex Selector
                                     vertexIdMVAweightfile      = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_SL_2016.xml"),
                                     vertexProbMVAweightfile    = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxProb_SL_2016.xml"),
                                     vertexIdMVAweightfileH4G = cms.untracked.FileInPath("%s"%vertexIdMVAweightfileH4G),
                                     # vertexIdMVAweightfileH4G = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_H4G_Total_2016.xml"),
                                     vertexProbMVAweightfileH4G = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxProb_H4G_Total_2016.xml"),
                                     diphoPairMVAweightfileH4G = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_diphoPairTMVA_H4G_2016_Combined_GB.weights.xml"),
                                     useSingleLeg               = cms.bool(True),
                                     doH4GVertex                = cms.bool(True),
                                     # saveDiphoPairingTree       = cms.bool(''),
                                     saveDiphoPairingTree       = cms.bool(True),
                                     useZerothVertexFromMicro   = cms.bool(False),

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
