from ROOT import *

Cut = '1>0'

outputLoc = '/eos/user/t/twamorka/www/H4Gamma/CompareH4GVtxVars_29Jun2020/Signal/'

nbins = 100

Vars = []

Vars.append(['ptAsym','ptAsym',';ptAsym; Normalized Yields',nbins,0, 20])
Vars.append(['ptBal','ptBal',';ptBal; Normalized Yields',nbins,0, 200])
Vars.append(['logSumpt2','logSumpt2',';logSumpt2; Normalized Yields',nbins,0, 20])
Vars.append(['pullConv','pullConv',';pullConv; Normalized Yields',nbins,0, 20])
Vars.append(['nConv','nConv',';nConv; Normalized Yields',nbins,0, 40])

files = []
files.append(['/eos/cms/store/user/bmarzocc/H4G_Analysis/Dumpers_vtxID/SUSYGluGluToHToAA_AToGG_Total_TuneCUETP8M1_13TeV_pythia8.root/h4gCandidateDumper_vtxBDT_sig/trees/SUSYGluGluToHToAA_AToGG_TuneCUETP8M1_13TeV_pythia8_13TeV_4photons_sig','2016',kRed,0,1])
files.append(['/eos/user/t/twamorka/15Nov_2017_VtxBDTTree/signal_m_BDTVtx.root/h4gCandidateDumper_vtxBDT_sig/trees/SUSYGluGluToHToAA_AToGG_TuneCP5_13TeV_pythia8_13TeV_4photons_sig','2017',kBlue,0,1])
