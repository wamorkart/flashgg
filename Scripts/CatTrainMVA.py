import ROOT
import argparse

parser =  argparse.ArgumentParser(description='cat MVA')
parser.add_argument('-t', '--training', dest='training', required=True, type=str)
parser.add_argument('-m', '--mass', dest='mass', required=True, type=str)
parser.add_argument('-y', '--year', dest='year', required=True, type=str)
parser.add_argument('-o', '--output', dest='output', required=True, type=str)

args = parser.parse_args()
training = args.training
sig_list = args.mass.split(',')
year = args.year
output = args.output


Cut_Signal = '(pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && abs(pho1_eta) < 2.5 && abs(pho2_eta) < 2.5 && abs(pho3_eta) < 2.5 && abs(pho4_eta) < 2.5 && (abs(pho1_eta) < 1.4442 || abs(pho1_eta) > 1.566) && (abs(pho2_eta) < 1.4442 || abs(pho2_eta) > 1.566) && (abs(pho3_eta) < 1.4442 || abs(pho3_eta) > 1.566) && (abs(pho4_eta) < 1.4442 || abs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 110 && tp_mass <180 && pho1_MVA > -999. && pho2_MVA > -999. && pho3_MVA > -999. && pho4_MVA > -999)'

Cut_Background = '(pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && abs(pho1_eta) < 2.5 && abs(pho2_eta) < 2.5 && abs(pho3_eta) < 2.5 && abs(pho4_eta) < 2.5 && (abs(pho1_eta) < 1.4442 || abs(pho1_eta) > 1.566) && (abs(pho2_eta) < 1.4442 || abs(pho2_eta) > 1.566) && (abs(pho3_eta) < 1.4442 || abs(pho3_eta) > 1.566) && (abs(pho4_eta) < 1.4442 || abs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 110 && tp_mass <180 && pho1_MVA > -999. && pho2_MVA > -999. && pho3_MVA > -999. && pho4_MVA > -999)'

sig_nums_2016 = []
sig_nums_2017 = []
sig_nums_2018 = []
sig_file_2016 = []
sig_file_2017 = []
sig_file_2018 = []
#
for mass in sig_list:

   sig_file_2016_tmp = ROOT.TChain()
   sig_file_2017_tmp = ROOT.TChain()
   sig_file_2018_tmp = ROOT.TChain()

   if (training=='Standard'):
       tree_name_2016 = '/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2016/hadd/signal_m_'+str(mass)+'.root/tagsDumper/trees/SUSYGluGluToHToAA_AToGG_M_'+str(mass)+'_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0'
       tree_name_2017 = '/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2017/hadd/signal_m_'+str(mass)+'.root/tagsDumper/trees/SUSYGluGluToHToAA_AToGG_M_'+str(mass)+'_TuneCP5_13TeV_pythia8_13TeV_H4GTag_0'
       tree_name_2018 = '/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2018/hadd/signal_m_'+str(mass)+'.root/tagsDumper/trees/HAHMHToAA_AToGG_MA_'+str(mass)+'GeV_TuneCP5_PSweights_13TeV_madgraph_pythia8_13TeV_H4GTag_0'
   elif (training=='GenMass'):
       tree_name_2016 = '/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2016/hadd/signal_m_'+str(mass)+'_genMass.root/SUSYGluGluToHToAA_AToGG_M_'+str(mass)+'_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0'
       tree_name_2017 = '/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2017/hadd/signal_m_'+str(mass)+'_genMass.root/SUSYGluGluToHToAA_AToGG_M_'+str(mass)+'_TuneCP5_13TeV_pythia8_13TeV_H4GTag_0'
       tree_name_2018 = '/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2018/hadd/signal_m_'+str(mass)+'_genMass.root/HAHMHToAA_AToGG_MA_'+str(mass)+'GeV_TuneCP5_PSweights_13TeV_madgraph_pythia8_13TeV_H4GTag_0'

   sig_file_2016_tmp.AddFile(tree_name_2016)
   sig_file_2016.append(sig_file_2016_tmp)
   sig_nums_2016.append(float(sig_file_2016_tmp.GetEntries(Cut_Signal)))

   sig_file_2017_tmp.AddFile(tree_name_2017)
   sig_file_2017.append(sig_file_2017_tmp)
   sig_nums_2017.append(float(sig_file_2017_tmp.GetEntries(Cut_Signal)))

   sig_file_2018_tmp.AddFile(tree_name_2018)
   sig_file_2018.append(sig_file_2018_tmp)
   sig_nums_2018.append(float(sig_file_2018_tmp.GetEntries(Cut_Signal)))

lumi_2016 = 35.9
lumi_2017 = 41.5
lumi_2018 = 54.38

bkg_file_2016 = ROOT.TChain()
bkg_file_2017 = ROOT.TChain()
bkg_file_2018 = ROOT.TChain()

if (training=='Standard'):
    bkg_file_2016.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2016/hadd/data_mix_weight_v4.root/Data_13TeV_H4GTag_0')
    bkg_file_2017.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2017/hadd/data_mix_weight_v4.root/Data_13TeV_H4GTag_0')
    bkg_file_2018.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2018/hadd/data_mix_weight_v4.root/Data_13TeV_H4GTag_0')
elif (training=='GenMass'):
    bkg_file_2016.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2016/hadd/data_mix_weight_v4_genMass.root/Data_13TeV_H4GTag_0')
    bkg_file_2017.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2017/hadd/data_mix_weight_v4_genMass.root/Data_13TeV_H4GTag_0')
    bkg_file_2018.AddFile('/eos/user/t/twamorka/h4g_fullRun2/withSystematics/2018/hadd/data_mix_weight_v4_genMass.root/Data_13TeV_H4GTag_0')

norm_2016 = float(bkg_file_2016.GetEntries(Cut_Background))
norm_2017 = float(bkg_file_2017.GetEntries(Cut_Background))
norm_2018 = float(bkg_file_2018.GetEntries(Cut_Background))

f_out = ROOT.TFile(output+'.root','RECREATE')
ROOT.TMVA.Tools.Instance()
factory = ROOT.TMVA.Factory("TMVAClassification", f_out,"AnalysisType=Classification")

dataloader = ROOT.TMVA.DataLoader("dataset")
mvaVars = []
print "Type of training: ", training
if (training=='Standard'):
    mvaVars = ['a1_mass_dM-a2_mass_dM', 'cosTheta_a1_dM', 'a1_pt_dM', 'a2_pt_dM', 'pho1_MVA', 'pho2_MVA', 'pho3_MVA', 'pho4_MVA']
elif (training=='GenMass'):
    mvaVars = ['genMass', 'a1_mass_dM-a2_mass_dM', 'cosTheta_a1_dM', 'a1_pt_dM', 'a2_pt_dM', 'pho1_MVA', 'pho2_MVA', 'pho3_MVA', 'pho4_MVA']
for x in mvaVars:
   #factory.AddVariable(x,"F")
   print "Variables used for training: ", x
   dataloader.AddVariable(x,"F")

print "Year being considered: ", year
for i,mass in enumerate(sig_list):
   if (year=='2016'):
       print "---> Signal normalization 2016: ",i,mass,float(sig_nums_2016[i])
       dataloader.AddSignalTree(sig_file_2016[i],lumi_2016/float(sig_nums_2016[i]))
   elif (year=='2017'):
       print "---> Signal normalization 2017: ",i,mass,float(sig_nums_2017[i])
       dataloader.AddSignalTree(sig_file_2017[i],lumi_2017/float(sig_nums_2017[i]))
   elif (year=='2018'):
       print "---> Signal normalization 2018: ",i,mass,float(sig_nums_2018[i])
       dataloader.AddSignalTree(sig_file_2018[i],lumi_2018/float(sig_nums_2018[i]))
   elif (year=='Run2'):
       print "---> Signal normalization 2016: ",i,mass,float(sig_nums_2016[i])
       dataloader.AddSignalTree(sig_file_2016[i],lumi_2016/float(sig_nums_2016[i]))
       print "---> Signal normalization 2017: ",i,mass,float(sig_nums_2017[i])
       dataloader.AddSignalTree(sig_file_2017[i],lumi_2017/float(sig_nums_2017[i]))
       print "---> Signal normalization 2018: ",i,mass,float(sig_nums_2018[i])
       dataloader.AddSignalTree(sig_file_2018[i],lumi_2018/float(sig_nums_2018[i]))

if (year=='2016'):
    print "---> Background normalization 2016: ",norm_2016
    dataloader.AddBackgroundTree(bkg_file_2016,lumi_2016/norm_2016)
elif (year=='2017'):
    print "---> Background normalization 2017: ",norm_2017
    dataloader.AddBackgroundTree(bkg_file_2017,lumi_2017/norm_2017)
elif (year=='2018'):
    print "---> Background normalization 2018: ",norm_2018
    dataloader.AddBackgroundTree(bkg_file_2018,lumi_2018/norm_2018)
elif (year=='Run2'):
    print "---> Background normalization 2016: ",norm_2016
    dataloader.AddBackgroundTree(bkg_file_2016,lumi_2016/norm_2016)
    print "---> Background normalization 2017: ",norm_2017
    dataloader.AddBackgroundTree(bkg_file_2017,lumi_2017/norm_2017)
    print "---> Background normalization 2018: ",norm_2018
    dataloader.AddBackgroundTree(bkg_file_2018,lumi_2018/norm_2018)

sigCut = ROOT.TCut(Cut_Signal)
bkgCut = ROOT.TCut(Cut_Background)

print "S Cut: ", sigCut
print "B Cut: ", bkgCut

dataloader.SetWeightExpression("weight","Signal")
dataloader.SetWeightExpression("weight","Background")

dataloader.PrepareTrainingAndTestTree(sigCut,bkgCut,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V")
method = factory.BookMethod( dataloader,ROOT.TMVA.Types.kBDT, "BDTG", "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" )

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
