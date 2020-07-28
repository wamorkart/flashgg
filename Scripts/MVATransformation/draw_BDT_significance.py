#!/usr/bin/python
import numpy as n
from ROOT import *
import sys, getopt
from array import array
from optparse import OptionParser
import operator
import math

def drawRatio(sig_file, bkg_file, lumi, label, name, markerStyle, markerColor):

   gStyle.SetOptStat(0000)

   Cut_Signal = 'weight*(pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && abs(pho1_eta) < 2.5 && abs(pho2_eta) < 2.5 && abs(pho3_eta) < 2.5 && abs(pho4_eta) < 2.5 && (abs(pho1_eta) < 1.4442 || abs(pho1_eta) > 1.566) && (abs(pho2_eta) < 1.4442 || abs(pho2_eta) > 1.566) && (abs(pho3_eta) < 1.4442 || abs(pho3_eta) > 1.566) && (abs(pho4_eta) < 1.4442 || abs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 110 && tp_mass <180 && pho1_MVA > -999. && pho2_MVA > -999. && pho3_MVA > -999. && pho4_MVA > -999)'

   Cut_Background = '(pho1_pt > 30 && pho2_pt > 18 && pho3_pt > 15 && pho4_pt > 15 && abs(pho1_eta) < 2.5 && abs(pho2_eta) < 2.5 && abs(pho3_eta) < 2.5 && abs(pho4_eta) < 2.5 && (abs(pho1_eta) < 1.4442 || abs(pho1_eta) > 1.566) && (abs(pho2_eta) < 1.4442 || abs(pho2_eta) > 1.566) && (abs(pho3_eta) < 1.4442 || abs(pho3_eta) > 1.566) && (abs(pho4_eta) < 1.4442 || abs(pho4_eta) > 1.566) && pho1_electronveto==1 && pho2_electronveto==1 && pho3_electronveto==1 && pho4_electronveto==1 && tp_mass > 110 && tp_mass <180 && pho1_MVA > -999. && pho2_MVA > -999. && pho3_MVA > -999. && pho4_MVA > -999)'

   h_bdt_signal = TH1F('h_bdt_signal','h_bdt_signal',100,-1.,1.)
   h_bdt_bkg = TH1F('h_bdt_bkg','h_bdt_bkg',100,-1.,1.)

   sig_file.Draw("bdt>>h_bdt_signal",Cut_Signal)
   bkg_file.Draw("bdt>>h_bdt_bkg",Cut_Background)

   # h_bdt_signal.Scale(lumi)
   # h_bdt_bkg.Scale(1./10.)
   h_bdt_signal.Scale(1./h_bdt_signal.Integral())
   h_bdt_bkg.Scale(1./h_bdt_bkg.Integral())

   h_ratio = h_bdt_signal.Clone()
   h_ratio.Reset()

   for bin in range(1,h_bdt_signal.GetNbinsX()+1):
      if(h_bdt_bkg.GetBinContent(bin)>0.):
         h_ratio.SetBinContent(bin, h_bdt_signal.GetBinContent(bin)/math.sqrt(h_bdt_bkg.GetBinContent(bin)))

   h_ratio.SetMarkerStyle(markerStyle)
   h_ratio.SetMarkerColor(markerColor)
   h_ratio.GetXaxis().SetTitle(label)
   h_ratio.GetYaxis().SetTitle("S/#sqrt{B}")
   h_ratio.SetTitle("")

   h_bdt_signal_cum = h_bdt_signal.GetCumulative(False)
   h_bdt_bkg_cum = h_bdt_bkg.GetCumulative(False)

   h_ratio_cum = h_bdt_signal_cum.Clone()
   h_ratio_cum.Reset()

   for bin in range(1,h_bdt_signal_cum.GetNbinsX()+1):
      if(h_bdt_bkg_cum.GetBinContent(bin)>0.):
         h_ratio_cum.SetBinContent(bin, h_bdt_signal_cum.GetBinContent(bin)/math.sqrt(h_bdt_bkg_cum.GetBinContent(bin)))

   h_ratio_cum.SetMarkerStyle(markerStyle)
   h_ratio_cum.SetMarkerColor(markerColor)
   h_ratio_cum.GetXaxis().SetTitle(label)
   h_ratio_cum.GetYaxis().SetTitle("S/#sqrt{B}")
   h_ratio_cum.SetTitle("")


   return h_ratio, h_ratio_cum

if __name__ == '__main__':

  input =  'CatMVA_PhoMVA_Only'
  output = '/eos/user/t/twamorka/www/H4Gamma/Jul2020_TaggerPlots/CatBDT_ReducedSignal/'

  ## Files used for Signal and Bkgs
  bkg_file_2016 = TChain()
  bkg_file_2016.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2016/hadd/ReducedSignal_CatMVA/'+input+'/data_mix.root/Data_13TeV_H4GTag_0')

  sig_file_2016 = TChain()
  sig_file_2016.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2016/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_60.root/SUSYGluGluToHToAA_AToGG_M_60_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
  sig_file_2016.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2016/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_45.root/SUSYGluGluToHToAA_AToGG_M_45_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
  sig_file_2016.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2016/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_35.root/SUSYGluGluToHToAA_AToGG_M_35_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
  sig_file_2016.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2016/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_25.root/SUSYGluGluToHToAA_AToGG_M_25_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')
  sig_file_2016.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2016/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_15.root/SUSYGluGluToHToAA_AToGG_M_15_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0')

  bkg_file_2017 = TChain()
  bkg_file_2017.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2017/hadd/ReducedSignal_CatMVA/'+input+'/data_mix.root/Data_13TeV_H4GTag_0')

  sig_file_2017 = TChain()
  sig_file_2017.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2017/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_60.root/SUSYGluGluToHToAA_AToGG_M_60_TuneCP5_13TeV_pythia8_13TeV_H4GTag_0')
  sig_file_2017.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2017/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_45.root/SUSYGluGluToHToAA_AToGG_M_45_TuneCP5_13TeV_pythia8_13TeV_H4GTag_0')
  sig_file_2017.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2017/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_35.root/SUSYGluGluToHToAA_AToGG_M_35_TuneCP5_13TeV_pythia8_13TeV_H4GTag_0')
  sig_file_2017.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2017/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_25.root/SUSYGluGluToHToAA_AToGG_M_25_TuneCP5_13TeV_pythia8_13TeV_H4GTag_0')
  sig_file_2017.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2017/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_15.root/SUSYGluGluToHToAA_AToGG_M_15_TuneCP5_13TeV_pythia8_13TeV_H4GTag_0')

  bkg_file_2018 = TChain()
  bkg_file_2018.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2018/hadd/ReducedSignal_CatMVA/'+input+'/data_mix.root/Data_13TeV_H4GTag_0')

  sig_file_2018 = TChain()
  sig_file_2018.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2018/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_60.root/HAHMHToAA_AToGG_MA_60GeV_TuneCP5_PSweights_13TeV_madgraph_pythia8_13TeV_H4GTag_0')
  sig_file_2018.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2018/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_45.root/HAHMHToAA_AToGG_MA_45GeV_TuneCP5_PSweights_13TeV_madgraph_pythia8_13TeV_H4GTag_0')
  sig_file_2018.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2018/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_35.root/HAHMHToAA_AToGG_MA_35GeV_TuneCP5_PSweights_13TeV_madgraph_pythia8_13TeV_H4GTag_0')
  sig_file_2018.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2018/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_25.root/HAHMHToAA_AToGG_MA_25GeV_TuneCP5_PSweights_13TeV_madgraph_pythia8_13TeV_H4GTag_0')
  sig_file_2018.AddFile('/eos/user/t/twamorka/h4g_fullRun2/noSystematics/2018/hadd/ReducedSignal_CatMVA/'+input+'/signal_m_15.root/HAHMHToAA_AToGG_MA_15GeV_TuneCP5_PSweights_13TeV_madgraph_pythia8_13TeV_H4GTag_0')


  h_ratio_2016 = drawRatio(sig_file_2016, bkg_file_2016, 35.9, 'bdt', '2016', 20, kBlack)[0]
  h_ratio_2017 = drawRatio(sig_file_2017, bkg_file_2017, 41.5, 'bdt', '2017', 21, kRed)[0]
  h_ratio_2018 = drawRatio(sig_file_2018, bkg_file_2018, 54.38, 'bdt', '2018', 22, kGreen)[0]

  c = TCanvas()
  h_ratio_2018.Draw('P same')
  h_ratio_2017.Draw('P same')
  h_ratio_2016.Draw('P same')

  leg = TLegend(0.156642, 0.693295, 0.447368, 0.883024)
  leg.SetBorderSize(0)
  leg.SetTextSize(0.04)
  leg.SetFillColor(kWhite)
  leg.SetFillStyle(0)
  leg.AddEntry(h_ratio_2016, '2016','p')
  leg.AddEntry(h_ratio_2017, '2017','p')
  leg.AddEntry(h_ratio_2018, '2018','p')

  leg.Draw('same')
  c.SaveAs(output+input+'_normalized.pdf')
  c.SaveAs(output+input+'_normalized.png')

  h_ratio_2016_cum = drawRatio(sig_file_2016, bkg_file_2016, 35.9, 'bdt', '2016', 20, kBlack)[1]
  h_ratio_2017_cum = drawRatio(sig_file_2017, bkg_file_2017, 41.5, 'bdt', '2017', 21, kRed)[1]
  h_ratio_2018_cum = drawRatio(sig_file_2018, bkg_file_2018, 54.38, 'bdt', '2018', 22, kGreen)[1]

  c = TCanvas()

  h_ratio_2018_cum.Draw('P same')
  h_ratio_2017_cum.Draw('P same')
  h_ratio_2016_cum.Draw('P same')

  leg = TLegend(0.156642, 0.693295, 0.447368, 0.883024)
  leg.SetBorderSize(0)
  leg.SetTextSize(0.04)
  leg.SetFillColor(kWhite)
  leg.SetFillStyle(0)
  leg.AddEntry(h_ratio_2016_cum, '2016','p')
  leg.AddEntry(h_ratio_2017_cum, '2017','p')
  leg.AddEntry(h_ratio_2018_cum, '2018','p')

  leg.Draw('same')
  c.SaveAs(output+input+'_cumulative_normalized.pdf')
  c.SaveAs(output+input+'_cumulative_normalized.png')



  # drawRatio(h_bdt_signal, h_bdt_bkg, 'bdt', 'bdt_significance')
  # drawRatio(h_bdtTransformed_signal, h_bdtTransformed_bkg, 'bdtTransformed', 'bdtTransformed_significance')
  # drawRatio(h_bdt_signal.GetCumulative(False), h_bdt_bkg.GetCumulative(False), 'bdt', 'bdt_significance_cumulative')
  # drawRatio(h_bdtTransformed_signal.GetCumulative(False), h_bdtTransformed_bkg.GetCumulative(False), 'bdtTransformed', 'bdtTransformed_significance_cumulative')
