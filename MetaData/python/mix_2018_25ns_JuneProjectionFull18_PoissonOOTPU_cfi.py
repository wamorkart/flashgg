import FWCore.ParameterSet.Config as cms

# configuration to model pileup for initial physics phase
from SimGeneral.MixingModule.mixObjects_cfi import theMixObjects
from SimGeneral.MixingModule.mixPoolSource_cfi import *
from SimGeneral.MixingModule.digitizers_cfi import *

mix = cms.EDProducer("MixingModule",
    digitizers = cms.PSet(theDigitizers),
    LabelPlayback = cms.string(''),
    maxBunch = cms.int32(3),
    minBunch = cms.int32(-12), ## in terms of 25 nsec

    bunchspace = cms.int32(25), ##ns
    mixProdStep1 = cms.bool(False),
    mixProdStep2 = cms.bool(False),

    playback = cms.untracked.bool(False),
    useCurrentProcessOnly = cms.bool(False),

    input = cms.SecSource("EmbeddedRootSource",
        type = cms.string('probFunction'),
        nbPileupEvents = cms.PSet(
            probFunctionVariable = cms.vint32(
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99
            ),
            probValue = cms.vdouble(
                4.695341e-10, 1.206213e-06, 1.162593e-06, 6.118058e-06, 1.626767e-05,
                3.508135e-05, 7.12608e-05, 0.0001400641, 0.0002663403, 0.0004867473,
                0.0008469, 0.001394142, 0.002169081, 0.003198514, 0.004491138,
                0.006036423, 0.007806509, 0.00976048, 0.0118498, 0.01402411,
                0.01623639, 0.01844593, 0.02061956, 0.02273221, 0.02476554,
                0.02670494, 0.02853662, 0.03024538, 0.03181323, 0.03321895,
                0.03443884, 0.035448, 0.03622242, 0.03674106, 0.0369877,
                0.03695224, 0.03663157, 0.03602986, 0.03515857, 0.03403612,
                0.0326868, 0.03113936, 0.02942582, 0.02757999, 0.02563551,
                0.02362497, 0.02158003, 0.01953143, 0.01750863, 0.01553934,
                0.01364905, 0.01186035, 0.01019246, 0.008660705, 0.007275915,
                0.006043917, 0.004965276, 0.004035611, 0.003246373, 0.002585932,
                0.002040746, 0.001596402, 0.001238498, 0.0009533139, 0.0007282885,
                0.000552306, 0.0004158005, 0.0003107302, 0.0002304612, 0.0001696012,
                0.0001238161, 8.96531e-05, 6.438087e-05, 4.585302e-05, 3.23949e-05,
                2.271048e-05, 1.580622e-05, 1.09286e-05, 7.512748e-06, 5.140304e-06,
                3.505254e-06, 2.386437e-06, 1.625859e-06, 1.111865e-06, 7.663272e-07,
                5.350694e-07, 3.808318e-07, 2.781785e-07, 2.098661e-07, 1.642811e-07,
                1.312835e-07, 1.081326e-07, 9.141993e-08, 7.890983e-08, 6.91468e-08,
                6.119019e-08, 5.443693e-08, 4.85036e-08, 4.31486e-08, 3.822112e-08
            ),
            histoFileName = cms.untracked.string('histProbFunction.root'),
        ),
        sequential = cms.untracked.bool(False),
        manage_OOT = cms.untracked.bool(True),  ## manage out-of-time pileup
        ## setting this to True means that the out-of-time pileup
        ## will have a different distribution than in-time, given
        ## by what is described on the next line:
        OOT_type = cms.untracked.string('Poisson'),  ## generate OOT with a Poisson matching the number chosen for in-time
        #OOT_type = cms.untracked.string('fixed'),  ## generate OOT with a fixed distribution
        #intFixed_OOT = cms.untracked.int32(2),
        fileNames = FileNames
    ),
    mixObjects = cms.PSet(theMixObjects)
)
