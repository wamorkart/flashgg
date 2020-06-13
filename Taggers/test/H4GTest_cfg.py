import importlib
import subprocess
import FWCore.ParameterSet.Config as cms
import os
import FWCore.ParameterSet.VarParsing as VarParsing

from flashgg.Taggers.flashggH4GCandidate_cfi import FlashggH4GCandidate
from flashgg.Taggers.flashggPreselectedDiPhotons_cfi import flashggPreselectedDiPhotons
from flashgg.Taggers.flashggPreselectedDiPhotons_LowMass_cfi import flashggPreselectedDiPhotonsLowMass
from flashgg.Taggers.flashggPreselectedDiPhotons_cfi_LowMass2017 import flashggPreselectedDiPhotons_LowMass2017
import flashgg.Taggers.dumperConfigTools as cfgTools
from flashgg.MetaData.MetaConditionsReader import *
from flashgg.Taggers.flashggUpdatedIdMVADiPhotons_cfi import flashggUpdatedIdMVADiPhotons

from flashgg.Systematics.flashggDiPhotonSystematics_cfi import flashggDiPhotonSystematics



process = cms.Process("FLASHggH4GTest")


from flashgg.MetaData.JobConfig import customize

customize.options.register('stdDumper',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "stdDumper")

customize.options.register('vtxBDTDumper',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "vtxBDTDumper")

customize.options.register('vtxProbDumper',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "vtxProbDumper")

customize.options.register('isCondor',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "isCondor")

customize.options.register('year',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "year")

customize.options.register('isSignal',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "isSignal")

customize.options.register('mass',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "mass")

customize.options.register('doH4GVertex',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "doH4GVertex")

from flashgg.MetaData.JobConfig import customize
customize.parse()

###--Scales and smearing stuff
customize.metaConditions = MetaConditionsReader(customize.metaConditions)

###---H4G candidates production
process.FlashggH4GCandidate = FlashggH4GCandidate.clone()

if (customize.year == "2016"):
    print "2016 diphoton ID selection"
    process.FlashggH4GCandidate.idSelection = cms.PSet(
            rho = flashggPreselectedDiPhotonsLowMass.rho,
            cut = flashggPreselectedDiPhotonsLowMass.cut,
            variables = flashggPreselectedDiPhotonsLowMass.variables,
            categories = flashggPreselectedDiPhotonsLowMass.categories
            )
elif (customize.year == "2017"):
        print "2017 diphoton ID selection"
        process.FlashggH4GCandidate.idSelection = cms.PSet(
                rho = flashggPreselectedDiPhotons_LowMass2017.rho,
                cut = flashggPreselectedDiPhotons_LowMass2017.cut,
                variables = flashggPreselectedDiPhotons_LowMass2017.variables,
                categories = flashggPreselectedDiPhotons_LowMass2017.categories
                )


###--- get the variables
import flashgg.Taggers.H4GTagVariables as var
vtx_BDT_variables_sig = var.vtx_BDT_variables_sig
vtx_BDT_variables_bkg = var.vtx_BDT_variables_bkg
vtxProb_BDT_variables = var.vtx_variables + var.vtxProb_BDT_variables
all_variables = var.pho_variables + var.dipho_variables + var.vtx_variables + var.tp_variables + var.gen_variables

from flashgg.Taggers.h4gCandidateDumper_cfi import h4gCandidateDumper

process.h4gCandidateDumper_vtxBDT_sig = h4gCandidateDumper.clone()
process.h4gCandidateDumper_vtxBDT_sig.dumpTrees = True
process.h4gCandidateDumper_vtxBDT_sig.dumpWorkspace = True


cfgTools.addCategories(process.h4gCandidateDumper_vtxBDT_sig,
                        [
                            ("Reject", "", -1),
                            ("4photons_sig","phoP4Corrected_dp.size() > 3"),
                            ("3photons_sig","phoP4Corrected_dp.size() == 3", 0),
                            ("2photons_sig","phoP4Corrected_dp.size() == 2", 0)
                        ],
                        variables = vtx_BDT_variables_sig,
                        histograms=[]
                        )


process.h4gCandidateDumper_vtxBDT_bkg = h4gCandidateDumper.clone()
process.h4gCandidateDumper_vtxBDT_bkg.dumpTrees = True
process.h4gCandidateDumper_vtxBDT_bkg.dumpWorkspace = True

cfgTools.addCategories(process.h4gCandidateDumper_vtxBDT_bkg,
                        [
                            ("Reject", "", -1),
                            ("4photons_bkg","phoP4Corrected_dp.size() > 3"),
                            ("3photons_bkg","phoP4Corrected_dp.size() == 3", 0),
                            ("2photons_bkg","phoP4Corrected_dp.size() == 2", 0)
                        ],
                        variables = vtx_BDT_variables_bkg,
                        histograms=[]
                        )

process.h4gCandidateDumper_vtxProb = h4gCandidateDumper.clone()
process.h4gCandidateDumper_vtxProb.dumpTrees = True
process.h4gCandidateDumper_vtxProb.dumpWorkspace = True

cfgTools.addCategories(process.h4gCandidateDumper_vtxProb,
                        [
                            ("Reject", "", -1),
                            ("4photons","phoP4Corrected_dp.size() > 3"),
                            ("3photons","phoP4Corrected_dp.size() == 3", 0),
                            ("2photons","phoP4Corrected_dp.size() == 2", 0)
                        ],
                        variables = vtxProb_BDT_variables,
                        histograms=[]
                        )

process.h4gCandidateDumper = h4gCandidateDumper.clone()
process.h4gCandidateDumper.dumpTrees = True
process.h4gCandidateDumper.dumpWorkspace = True

cfgTools.addCategories(process.h4gCandidateDumper,
                       [
                            ("Reject", "", -1),
                            ("4photons","phoP4Corrected_dp.size() > 3 "),
                            ("3photons","phoP4Corrected_dp.size() == 3", 0),
                            ("2photons","phoP4Corrected_dp.size() == 2", 0)
                        ],
                        variables = all_variables,
                        histograms=[]
                        )


if (customize.year == "2016"):
    if (customize.doH4GVertex == 0):
        process.source = cms.Source ("PoolSource",
            fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/group/phys_higgs/HiggsExo/H4Gamma/H4G_2016samples_producedwithBkgCustomization/H4GandHH4G_2016_27Sep2019/RunIIFall18-4_0_0-119-g2d54185d/SUSYGluGluToHToAA_AToGG_M-60_TuneCUETP8M1_13TeV_pythia8/H4GandHH4G_2016_27Sep2019-RunIIFall18-4_0_0-119-g2d54185d-v0-RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/190927_180219/0000/myMicroAODOutputFile_211.root"))
            # fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/user/spigazzi/flashgg/Era2016_RR-17Jul2018_v2/legacyRun2FullV1/DoubleEG/Era2016_RR-17Jul2018_v2-legacyRun2FullV1-v0-Run2016B-17Jul2018_ver2-v1/190605_220256/0000/myMicroAODOutputFile_932.root"))
    else:
        process.source = cms.Source ("PoolSource",
            fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/group/phys_higgs/HiggsExo/H4Gamma/H4G_2016samples_producedwithBkgCustomization/H4GandHH4G_2016_27Sep2019/RunIIFall18-4_0_0-119-g2d54185d/SUSYGluGluToHToAA_AToGG_M-60_TuneCUETP8M1_13TeV_pythia8/H4GandHH4G_2016_27Sep2019-RunIIFall18-4_0_0-119-g2d54185d-v0-RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/190927_180219/0000/myMicroAODOutputFile_211.root"),
            secondaryFileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/SUSYGluGluToHToAA_AToGG_M-60_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/270000/D489FF02-343F-E911-8857-AC1F6BAC7C2A.root")
            # fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/user/spigazzi/flashgg/Era2016_RR-17Jul2018_v2/legacyRun2FullV1/DoubleEG/Era2016_RR-17Jul2018_v2-legacyRun2FullV1-v0-Run2016B-17Jul2018_ver2-v1/190605_220256/0000/myMicroAODOutputFile_932.root"),
            # secondaryFileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/data/Run2016B/DoubleEG/MINIAOD/17Jul2018_ver2-v1/20000/D03AED69-308D-E811-AFFC-008CFA197CD0.root","root://cms-xrd-global.cern.ch//store/data/Run2016B/DoubleEG/MINIAOD/17Jul2018_ver2-v1/20000/FCFDB07D-378D-E811-89A0-008CFAE45144.root")

            )


elif (customize.year == "2017"):
      process.source = cms.Source ("PoolSource",
              # print "reading 2017 files"
              fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/group/phys_higgs/HiggsExo/H4Gamma/H4G_2017samples_producedwithBkgCustomization/H4G_2017_27Sep2019/RunIIFall18-4_0_0-119-g2d54185d/SUSYGluGluToHToAA_AToGG_M-50_TuneCP5_13TeV_pythia8/H4G_2017_27Sep2019-RunIIFall18-4_0_0-119-g2d54185d-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/190927_153802/0000/myMicroAODOutputFile_2.root"),
              secondaryFileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/SUSYGluGluToHToAA_AToGG_M-50_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/E2F1936C-FE74-E811-96F8-0CC47A2AECFA.root")
              )
if (customize.year == "2016"):
    print "2016 vtx ID"
    process.FlashggH4GCandidate.vertexIdMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_H4G_Total_2016.xml")
    process.FlashggH4GCandidate.vertexProbMVAweightfileH4G = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxProb_H4G_Total_2016.xml")
    print "2016 diphoton Pairing"
    print "mass = ", customize.mass
    process.FlashggH4GCandidate.mass = cms.untracked.double(customize.mass)
    process.FlashggH4GCandidate.doH4GVertex = cms.untracked.bool(customize.doH4GVertex)
    if (customize.doH4GVertex == 0):
        print "Zeroth Vertex"
    else:
        print "H4G Vertex"
    if (customize.mass == 60):
        process.FlashggH4GCandidate.diphoPairMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_diphoPairTMVA_H4G_M60.weights.xml")
        # process.FlashggH4GCandidate.CatMVAweightfileH4G = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_CatMVA_H4G_60.weights.xml")
        print "BDT mass training ", process.FlashggH4GCandidate.diphoPairMVAweightfileH4G
    elif (customize.mass == 55):
        process.FlashggH4GCandidate.diphoPairMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_diphoPairTMVA_H4G_M55.weights.xml")
        print "BDT mass training ", process.FlashggH4GCandidate.diphoPairMVAweightfileH4G
    elif (customize.mass == 50):
        process.FlashggH4GCandidate.diphoPairMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_diphoPairTMVA_H4G_M50.weights.xml")
        print "BDT mass training ", process.FlashggH4GCandidate.diphoPairMVAweightfileH4G
        # process.FlashggH4GCandidate.CatMVAweightfileH4G = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_CatMVA_H4G_60.weights.xml")
    elif (customize.mass == 45):
        process.FlashggH4GCandidate.diphoPairMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_diphoPairTMVA_H4G_M45.weights.xml")
        print "BDT mass training ", process.FlashggH4GCandidate.diphoPairMVAweightfileH4G
    elif (customize.mass == 40):
        process.FlashggH4GCandidate.diphoPairMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_diphoPairTMVA_H4G_M40.weights.xml")
        print "BDT mass training ", process.FlashggH4GCandidate.diphoPairMVAweightfileH4G
    elif (customize.mass == 35):
        process.FlashggH4GCandidate.diphoPairMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_diphoPairTMVA_H4G_M35.weights.xml")
        print "BDT mass training ", process.FlashggH4GCandidate.diphoPairMVAweightfileH4G
    elif (customize.mass == 30):
        process.FlashggH4GCandidate.diphoPairMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_diphoPairTMVA_H4G_M30.weights.xml")
        print "BDT mass training ", process.FlashggH4GCandidate.diphoPairMVAweightfileH4G
    elif (customize.mass == 25):
        process.FlashggH4GCandidate.diphoPairMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_diphoPairTMVA_H4G_M25.weights.xml")
        print "BDT mass training ", process.FlashggH4GCandidate.diphoPairMVAweightfileH4G
    elif (customize.mass == 20):
        process.FlashggH4GCandidate.diphoPairMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_diphoPairTMVA_H4G_M20.weights.xml")
        print "BDT mass training ", process.FlashggH4GCandidate.diphoPairMVAweightfileH4G
    elif (customize.mass == 15):
        process.FlashggH4GCandidate.diphoPairMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_diphoPairTMVA_H4G_M15.weights.xml")
        print "BDT mass training ", process.FlashggH4GCandidate.diphoPairMVAweightfileH4G
    # process.FlashggH4GCandidate.diphoPairMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_diPairBDT_H4G_M-60_2016.xml")


elif (customize.year == "2017"):
     print "2017 vtx ID"
     process.FlashggH4GCandidate.vertexIdMVAweightfileH4G= cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxId_H4G_Total_2017.xml")
     process.FlashggH4GCandidate.vertexProbMVAweightfileH4G = cms.FileInPath("flashgg/MicroAOD/data/TMVAClassification_BDTVtxProb_H4G_Total_2017.xml")

if customize.isSignal:
   process.FlashggH4GCandidate.saveDiphoPairingTree = cms.bool(True)
elif not customize.isSignal:
   process.FlashggH4GCandidate.saveDiphoPairingTree = cms.bool(False)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root"),
                                   closeFileFast = cms.untracked.bool(True)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

if customize.inputFiles:
    inputFile = customize.inputFiles

# Require low mass diphoton triggers
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
if (customize.year == "2016"):
   print "applying 2016 triggers"
   process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(
                                                              # "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*",
                                                              # "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", ##--std Hgg diphoton trigger
                                                              "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v*", ##--low mass trigger
                                                              "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v*"   ##--low mass diphoton trigger
                                                               ))

elif (customize.year == "2017"):
     print "applying 2017 triggers"
     process.hltHighLevel = hltHighLevel.clone(HLTPaths = cms.vstring(
                                                                "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55*"
                                                                ))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#############   Geometry  ###############
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.load("Configuration.Geometry.GeometryECALHCAL_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)

process.dataRequirements = cms.Sequence()
# process.load("flashgg/Taggers/vtxH4GSequence")
process.dataRequirements += process.hltHighLevel
if customize.processId == "Data":
   process.dataRequirements += process.hltHighLevel
   process.dataRequirements += process.eeBadScFilter

process.load("flashgg/Taggers/vtxH4GSequence")

process.flashggUpdatedIdMVADiPhotons = flashggUpdatedIdMVADiPhotons
process.load("flashgg.Systematics."+customize.metaConditions["flashggDiPhotonSystematics"])

sysmodule = importlib.import_module(
    "flashgg.Systematics."+customize.metaConditions["flashggDiPhotonSystematics"])
systModules2D = cms.VPSet()
systModules = cms.VPSet()

if customize.processId == "Data":
    print'Data'
    systModules.append(sysmodule.MCScaleHighR9EB_EGM)
    systModules.append(sysmodule.MCScaleLowR9EB_EGM)
    systModules.append(sysmodule.MCScaleHighR9EE_EGM)
    systModules.append(sysmodule.MCScaleLowR9EE_EGM)
    # systModules.append(sysmodule.MCScaleGain6EB_EGM)
    # systModules.append(sysmodule.MCScaleGain1EB_EGM)

    for module in systModules:
        module.ApplyCentralValue = cms.bool(True)

else:
    print'Not Data'
    systModules.append(sysmodule.MCScaleHighR9EB_EGM)
    systModules.append(sysmodule.MCScaleLowR9EB_EGM)
    systModules.append(sysmodule.MCScaleHighR9EE_EGM)
    systModules.append(sysmodule.MCScaleLowR9EE_EGM)

    systModules2D.append(sysmodule.MCSmearHighR9EE_EGM)
    systModules2D.append(sysmodule.MCSmearLowR9EE_EGM)
    systModules2D.append(sysmodule.MCSmearHighR9EB_EGM)
    systModules2D.append(sysmodule.MCSmearLowR9EB_EGM)

    for module in systModules:
        module.ApplyCentralValue = cms.bool(False)

process.flashggDiPhotonSystematics = flashggDiPhotonSystematics
process.flashggDiPhotonSystematics.src = "flashggDiPhotons"
process.flashggDiPhotonSystematics.SystMethods = systModules
process.flashggDiPhotonSystematics.SystMethods2D = systModules2D

if customize.stdDumper:
   #standard dumper sequence
   if (customize.doH4GVertex):
       process.path = cms.Path(process.vtxH4GSequence*process.dataRequirements*process.flashggDiPhotonSystematics*process.FlashggH4GCandidate*process.h4gCandidateDumper)
   else:
       process.path = cms.Path(process.dataRequirements*process.flashggDiPhotonSystematics*process.FlashggH4GCandidate*process.h4gCandidateDumper)
   #process.path = cms.Path(process.vtxH4GSequence*process.dataRequirements*process.FlashggH4GCandidate*process.h4gCandidateDumper)

if customize.vtxBDTDumper:
   #vtxBDT dumper sequence
   process.path = cms.Path(process.vtxH4GSequence*process.dataRequirements*process.flashggDiPhotonSystematics*process.FlashggH4GCandidate*process.h4gCandidateDumper_vtxBDT_sig*process.h4gCandidateDumper_vtxBDT_bkg)

if customize.vtxProbDumper:
   #vtxProb dumper sequence
   process.path = cms.Path(process.vtxH4GSequence*process.dataRequirements*process.flashggDiPhotonSystematics*process.FlashggH4GCandidate*process.h4gCandidateDumper_vtxProb)

if customize.isCondor:
   customize(process)
