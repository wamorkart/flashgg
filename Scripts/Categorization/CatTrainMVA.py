import ROOT
import argparse
parser =  argparse.ArgumentParser(description='cat MVA')
parser.add_argument('-o', '--output', dest='output', required=True, type=str)
parser.add_argument('-oD', '--outputDir', dest='outputDir', required=True, type=str)
parser.add_argument('-WP','--WP',dest='WP',required = True, type=str)


args = parser.parse_args()
output = args.output
outputDir = args.outputDir
WP = args.WP

## Files used for VBF related training
bkg_file = ROOT.TChain()
bkg_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/data_mix_weight.root/Data_13TeV_H4GTag_0')

sig_file = ROOT.TChain()
sig_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/signal_m_60.root/tagsDumperDumper/trees/SUSYGluGluToHToAA_AToGG_M_60_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
sig_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/signal_m_55.root/tagsDumperDumper/trees/SUSYGluGluToHToAA_AToGG_M_55_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
sig_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/signal_m_50.root/tagsDumperDumper/trees/SUSYGluGluToHToAA_AToGG_M_50_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
sig_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/signal_m_45.root/tagsDumperDumper/trees/SUSYGluGluToHToAA_AToGG_M_45_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
sig_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/signal_m_40.root/tagsDumperDumper/trees/SUSYGluGluToHToAA_AToGG_M_40_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
sig_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/signal_m_35.root/tagsDumperDumper/trees/SUSYGluGluToHToAA_AToGG_M_35_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
sig_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/signal_m_30.root/tagsDumperDumper/trees/SUSYGluGluToHToAA_AToGG_M_30_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
sig_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/signal_m_25.root/tagsDumperDumper/trees/SUSYGluGluToHToAA_AToGG_M_25_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
sig_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/signal_m_20.root/tagsDumperDumper/trees/SUSYGluGluToHToAA_AToGG_M_20_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
sig_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/signal_m_15.root/tagsDumperDumper/trees/SUSYGluGluToHToAA_AToGG_M_15_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
sig_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/signal_m_10.root/tagsDumperDumper/trees/SUSYGluGluToHToAA_AToGG_M_10_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
sig_file.AddFile('/eos/user/t/twamorka/h4g_fullRun2/2016/hadd/signal_m_5.root/tagsDumperDumper/trees/SUSYGluGluToHToAA_AToGG_M_5_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')

f_out = ROOT.TFile(outputDir+output+'.root','RECREATE')

ROOT.TMVA.Tools.Instance()
factory = ROOT.TMVA.Factory("TMVAClassification", f_out,"AnalysisType=Classification")

mvaVars = [
# 'CTStarCS',
# 'CT_a1Pho1',
# 'CT_a2Pho1',
# 'a1_Pt/tp_mass',
# 'a2_Pt/tp_mass',
# 'a1_Pho1PtOvera1Mass',
# 'a2_Pho1PtOvera2Mass',
'pho1_MVA',
'pho2_MVA',
'pho3_MVA',
'pho4_MVA',
#'pairMVAscore'
# 'pho1_pt',
# 'pho2_pt',
# 'pho3_pt',
# 'pho4_pt',
# 'pho1_eta',
# 'pho2_eta',
# 'pho3_eta',
# 'pho4_eta',
# 'tp_pt',
# 'tp_eta'
]

dataloader = ROOT.TMVA.DataLoader("dataset")

for x in mvaVars:
    #factory.AddVariable(x,"F")
    dataloader.AddVariable(x,"F")

#factory.AddSignalTree(sig_file)
#factory.AddBackgroundTree(bkg_file)
dataloader.AddSignalTree(sig_file)
dataloader.AddBackgroundTree(bkg_file)

if (WP == 'veryLoose'):
    Cut_MVA = 'pho1_MVA > -0.9 && pho2_MVA > -0.9 && pho3_MVA > -0.9 && pho4_MVA > -0.9)'
elif (WP == 'Loose'):
    Cut_MVA = 'pho1_MVA > -0.9 && pho2_MVA > -0.9 && pho3_MVA > -0.75 && pho4_MVA > -0.75)'
elif (WP == 'Medium'):
    Cut_MVA = 'pho1_MVA > -0.2 && pho2_MVA > -0.4 && pho3_MVA > -0.75 && pho4_MVA > -0.75)'
else:
    Cut_MVA = 'pho1_MVA > -0.2 && pho2_MVA > -0.4 && pho3_MVA > -0.5 && pho4_MVA > -0.5)'

Cut_Signal = '(pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && abs(pho1_eta) < 2.5 && abs(pho2_eta) < 2.5 && abs(pho3_eta) < 2.5 && abs(pho4_eta) < 2.5 && (abs(pho1_eta) < 1.4442 || abs(pho1_eta) > 1.566) && (abs(pho2_eta) < 1.4442 || abs(pho2_eta) > 1.566) && (abs(pho3_eta) < 1.4442 || abs(pho3_eta) > 1.566) && (abs(pho4_eta) < 1.4442 || abs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 110 && tp_mass <180  &&'
# Cut_Background = '(pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && abs(pho1_eta) < 2.5 && abs(pho2_eta) < 2.5 && abs(pho3_eta) < 2.5 && abs(pho4_eta) < 2.5 && (abs(pho1_eta) < 1.4442 || abs(pho1_eta) > 1.566) && (abs(pho2_eta) < 1.4442 || abs(pho2_eta) > 1.566) && (abs(pho3_eta) < 1.4442 || abs(pho3_eta) > 1.566) && (abs(pho4_eta) < 1.4442 || abs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 110 && tp_mass <180 && !((tp_mass > 115 && tp_mass < 135)) && '
# Cut_Background_weight = '(mix_weight)*'

# Cut_Signal = 'weight_VBF'
# Cut_Background = '(weight_VBF)*(pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && abs(pho1_eta) < 2.5 && abs(pho2_eta) < 2.5 && abs(pho3_eta) < 2.5 && abs(pho4_eta) < 2.5 && (abs(pho1_eta) < 1.4442 || abs(pho1_eta) > 1.566) && (abs(pho2_eta) < 1.4442 || abs(pho2_eta) > 1.566) && (abs(pho3_eta) < 1.4442 || abs(pho3_eta) > 1.566) && (abs(pho4_eta) < 1.4442 || abs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 110 && tp_mass <180 && !((tp_mass > 115 && tp_mass < 135)) &&'

# Cut_Signal = '(1>0 &&'
# Cut_Background = '(mix_weight)*(pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && abs(pho1_eta) < 2.5 && abs(pho2_eta) < 2.5 && abs(pho3_eta) < 2.5 && abs(pho4_eta) < 2.5 && (abs(pho1_eta) < 1.4442 || abs(pho1_eta) > 1.566) && (abs(pho2_eta) < 1.4442 || abs(pho2_eta) > 1.566) && (abs(pho3_eta) < 1.4442 || abs(pho3_eta) > 1.566) && (abs(pho4_eta) < 1.4442 || abs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 110 && tp_mass <180 && !(tp_mass > 115 && tp_mass < 135)  &&  '
Cut_Background = '(pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && abs(pho1_eta) < 2.5 && abs(pho2_eta) < 2.5 && abs(pho3_eta) < 2.5 && abs(pho4_eta) < 2.5 && (abs(pho1_eta) < 1.4442 || abs(pho1_eta) > 1.566) && (abs(pho2_eta) < 1.4442 || abs(pho2_eta) > 1.566) && (abs(pho3_eta) < 1.4442 || abs(pho3_eta) > 1.566) && (abs(pho4_eta) < 1.4442 || abs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 110 && tp_mass <180 && !(tp_mass > 115 && tp_mass < 135)   &&'
# Cut_Background = '(mix_weight)*(!(tp_mass > 115 && tp_mass < 135) &&'
#
Cut_add = 'pho1_MVA > -999. && pho2_MVA > -999. && pho3_MVA > -999. && pho4_MVA > -999)'
sigCut = ROOT.TCut(Cut_Signal+Cut_add)
bkgCut = ROOT.TCut(Cut_Background+Cut_add)

# sigCut = ROOT.TCut(Cut_Signal)
# bkgCut = ROOT.TCut(Cut_Background)


print "S Cut: ", sigCut
print "B Cut: ", bkgCut

dataloader.PrepareTrainingAndTestTree(sigCut,bkgCut,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V")
method = factory.BookMethod( dataloader,ROOT.TMVA.Types.kBDT, "BDT", "!H:!V:NTrees=750:MinNodeSize=3%:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20")#:nCuts=200")

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
