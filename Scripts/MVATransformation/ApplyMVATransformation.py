#!/usr/bin/python
import numpy as n
from ROOT import *
import sys, getopt
from array import array
import operator
import math
import argparse

if __name__ == '__main__':


  parser =  argparse.ArgumentParser(description='BDT transformation')
  parser.add_argument( "-id", "--inputdir",    dest="inputdir",  required=True, type=str, help="input directory" )
  parser.add_argument( "-iF", "--inputfile",    dest="inputfile",  required=True, type=str, help="input file" )
  parser.add_argument( "-m", "--mass",    dest="mass",  required=False, type=str, help="mass" )
  parser.add_argument( "-y", "--year",    dest="year",  required=True, type=str, help="year" )

  opt = parser.parse_args()
  print opt.inputdir
  print opt.inputfile
  fin = TFile.Open(opt.inputdir+opt.inputfile)
  fTransformed = TFile.Open((opt.inputdir+opt.inputfile).replace(".root","")+"_transformedMVA.root","recreate")

  cumulativeGraph = TFile(opt.inputdir+'cumulativeTransformation.root').Get('cumulativeGraph')
  print "cumulativeGraph: ", opt.inputdir+'cumulativeTransformation.root'

  fileName_trans= opt.inputdir+opt.inputfile
  tree = ''
  if ('signal' in opt.inputfile and opt.year == '2016'):
      tree = 'SUSYGluGluToHToAA_AToGG_M_'+str(opt.mass)+'_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0'
  elif ('signal' in opt.inputfile and opt.year == '2017'):
      tree = 'SUSYGluGluToHToAA_AToGG_M_'+str(opt.mass)+'_TuneCP5_13TeV_pythia8_13TeV_H4GTag_0'
  elif ('signal' in opt.inputfile and opt.year == '2018'):
      tree = 'HAHMHToAA_AToGG_MA_'+str(opt.mass)+'GeV_TuneCP5_PSweights_13TeV_madgraph_pythia8_13TeV_H4GTag_0'
  elif ('data' in opt.inputfile):
      tree = 'Data_13TeV_H4GTag_0'

  print 'Tree: ', tree

  f_trans = TFile(fileName_trans)
  t_trans = f_trans.Get(tree)

  fTransformed = TFile.Open((opt.inputdir+opt.inputfile).replace(".root","")+"_transformedMVA.root","recreate")
  chain = TChain(t_trans.GetName())
  chain.Add(opt.inputdir+opt.inputfile)
  copyTree = chain.CopyTree("")
  copyTree.SetName(tree)
  copyTree.SetTitle(tree)

  transfMVA = array( 'f', [ 0. ] )
  transfBranch = copyTree.Branch("bdtTransformed",transfMVA,"bdtTransformed/F");

  for i,event in enumerate(copyTree):
     if i>copyTree.GetEntries():break

     mva = event.bdt
     transfMVA[0] = cumulativeGraph.Eval(mva)
     # print i, "  " , transfMVA[0]
     transfBranch.Fill()

  fTransformed.Write()
  fTransformed.Close()
