from ROOT import *
from array import array

import argparse
parser =  argparse.ArgumentParser(description='Add Classification BDT weights')
parser.add_argument('-f', '--file', dest='File', required=True, type=str)
parser.add_argument('-m', '--mass', dest='mass', required=False, type=str)
parser.add_argument('-y', '--year', dest='year', required=False, type=str)
parser.add_argument('-t', '--training', dest='training', required=True, type=str)
parser.add_argument('-gm', '--genMass', dest='genMass', required=False, type=float)
parser.add_argument('-W', '--weight', dest='Weight', required=True, type=str)
# parser.add_argument('-T', '--training', dest='training', required=True, type=str)
parser.add_argument('-O', '--output', dest='Out', required=True, type=str)

opt = parser.parse_args()
training = opt.training

reader = TMVA.Reader()

# mvaVars = ['a1_mass_dM-a2_mass_dM', 'cosTheta_a1_dM', 'a1_pt_dM', 'a2_pt_dM', 'pho1_MVA', 'pho2_MVA', 'pho3_MVA', 'pho4_MVA', 'pho1_pixelseed', 'pho2_pixelseed', 'pho3_pixelseed' , 'pho4_pixelseed', 'a1_pt_dM/a1_mass_dM', 'a2_pt_dM/a2_mass_dM', 'a1_pt_dM/tp_mass', 'a2_pt_dM/tp_mass', 'a1_a2_dR_dM']



# genMass = array('f',[0])
diffdM = array('f',[0])
ct1 = array('f',[0])
a1pt = array('f',[0])
a2pt = array('f',[0])
p1mva = array('f',[0])
p2mva = array('f',[0])
p3mva = array('f',[0])
p4mva = array('f',[0])
# p1_PSV = array('f',[0])
# p2_PSV = array('f',[0])
# p3_PSV = array('f',[0])
# p4_PSV = array('f',[0])
# a1_pt_Over_a1_mass = array('f',[0])
# a2_pt_Over_a2_mass = array('f',[0])
# a1_pt_Over_tp_mass = array('f',[0])
# a2_pt_Over_tp_mass = array('f',[0])
a1_a2_dR = array('f',[0])


# if (training=='GenMass'):
    # reader.AddVariable('genMass',genMass)
reader.AddVariable('a1_mass_dM-a2_mass_dM',diffdM)
reader.AddVariable('cosTheta_a1_dM',ct1)
reader.AddVariable('a1_pt_dM',a1pt)
reader.AddVariable('a2_pt_dM',a2pt)
# reader.AddVariable('pho1_MVA',p1mva)
# reader.AddVariable('pho2_MVA',p2mva)
# reader.AddVariable('pho3_MVA',p3mva)
# reader.AddVariable('pho4_MVA',p4mva)
reader.AddVariable('pho1_MVA<-1.? -1.1: pho1_MVA',p1mva)
reader.AddVariable('pho2_MVA<-1.? -1.1: pho2_MVA',p2mva)
reader.AddVariable('pho3_MVA<-1.? -1.1: pho3_MVA',p3mva)
reader.AddVariable('pho4_MVA<-1.? -1.1: pho4_MVA',p4mva)
# reader.AddVariable('pho1_pixelseed',p1_PSV)
# reader.AddVariable('pho2_pixelseed',p2_PSV)
# reader.AddVariable('pho3_pixelseed',p3_PSV)
# reader.AddVariable('pho4_pixelseed',p4_PSV)
# reader.AddVariable('a1_pt_dM/a1_mass_dM',a1_pt_Over_a1_mass)
# reader.AddVariable('a2_pt_dM/a2_mass_dM',a2_pt_Over_a2_mass)
# reader.AddVariable('a1_pt_dM/tp_mass',a1_pt_Over_tp_mass)
# reader.AddVariable('a2_pt_dM/tp_mass',a2_pt_Over_tp_mass)
reader.AddVariable('a1_a2_dR_dM',a1_a2_dR)

# reader.BookMVA("BDT","TrainingWeights/"+opt.Weight)
reader.BookMVA("BDT",opt.Weight)


systLabels = [""]
# for direction in ["Up","Down"]:
#            systLabels.append("MvaShift%s01sigma"%direction)
#            systLabels.append("SigmaEOverEShift%s01sigma"%direction)
#            systLabels.append("MaterialCentralBarrel%s01sigma"%direction)
#            systLabels.append("MaterialOuterBarrel%s01sigma"%direction)
#            systLabels.append("MaterialForward%s01sigma"%direction)
#            systLabels.append("FNUFEB%s01sigma"%direction)
#            systLabels.append("FNUFEE%s01sigma"%direction)
#            systLabels.append("MCScaleGain6EB%s01sigma"%direction)
#            systLabels.append("MCScaleGain1EB%s01sigma"%direction)
#            systLabels.append("MCScaleGain1EB%s01sigma"%direction)
#
#            systLabels.append("MCScaleHighR9EB%s01sigma" % direction)
#            systLabels.append("MCScaleHighR9EE%s01sigma" % direction)
#            systLabels.append("MCScaleLowR9EB%s01sigma" % direction)
#            systLabels.append("MCScaleLowR9EE%s01sigma" % direction)
#            systLabels.append("MCSmearHighR9EBPhi%s01sigma" % direction)
#            systLabels.append("MCSmearHighR9EBRho%s01sigma" % direction)
#            systLabels.append("MCSmearHighR9EEPhi%s01sigma" % direction)
#            systLabels.append("MCSmearHighR9EERho%s01sigma" % direction)
#            systLabels.append("MCSmearLowR9EBPhi%s01sigma" % direction)
#            systLabels.append("MCSmearLowR9EBRho%s01sigma" % direction)
#            systLabels.append("MCSmearLowR9EEPhi%s01sigma" % direction)
#            systLabels.append("MCSmearLowR9EERho%s01sigma" % direction)
#            systLabels.append("ShowerShapeHighR9EB%s01sigma" % direction)
#            systLabels.append("ShowerShapeHighR9EE%s01sigma" % direction)
#            systLabels.append("ShowerShapeLowR9EB%s01sigma" % direction)
#            systLabels.append("ShowerShapeLowR9EE%s01sigma" % direction)
#            systLabels.append("SigmaEOverEShift%s01sigma" % direction)

# print systLabels

FilesToRedo = opt.File.split(',')
for f in FilesToRedo:
    print f
    infilename =  f
    infile = TFile(infilename)
    print "file: ", infile
    tree = ''
    outfile = TFile(opt.Out, "RECREATE")
    if ('signal' in f):
        for sys_i,syst in enumerate(systLabels):
            print('on systematic',sys_i,'/',len(systLabels),':',syst)
            systLabel = ""
            print (syst)
            if syst != "":
               systLabel += '_' + syst
            if (opt.year == "2016"):
                tree = infile.Get("tagsDumper/trees/SUSYGluGluToHToAA_AToGG_M_"+opt.mass+"_TuneCUETP8M1_13TeV_pythia8_13TeV_H4GTag_0" + systLabel)
            elif (opt.year == "2017"):
                tree = infile.Get("tagsDumper/trees/SUSYGluGluToHToAA_AToGG_M_"+opt.mass+"_TuneCP5_13TeV_pythia8_13TeV_H4GTag_0" + systLabel)
            else:
                tree = infile.Get("tagsDumper/trees/HAHMHToAA_AToGG_MA_"+opt.mass+"GeV_TuneCP5_PSweights_13TeV_madgraph_pythia8_13TeV_H4GTag_0" + systLabel)
            nentries = tree.GetEntries()
            print "Signal: ", nentries
            outtree = tree.CloneTree(0)
            bdt = array('f', [0])
            _bdt = outtree.Branch('bdt', bdt, 'bdt/F')
            # print "Gen mass:",float(opt.genMass)
            for i in range(0, nentries):
                if i%1000 == 0: print i
                tree.GetEntry(i)
                # if (training=='GenMass'):
                #     genMass[0] = float(opt.genMass)
                diffdM[0] = float((tree.a1_mass_dM-tree.a2_mass_dM))
                ct1[0] = float(tree.cosTheta_a1_dM)
                a1pt[0] = float(tree.a1_pt_dM)
                a2pt[0] = float(tree.a2_pt_dM)
                # p1mva[0] = float(tree.pho1_MVA)
                # p2mva[0] = float(tree.pho2_MVA)
                # p3mva[0] = float(tree.pho3_MVA)
                # p4mva[0] = float(tree.pho4_MVA)
                p1mva[0] = float(tree.pho1_MVA<-1.? -1.1: tree.pho1_MVA)
                p2mva[0] = float(tree.pho2_MVA<-1.? -1.1: tree.pho2_MVA)
                p3mva[0] = float(tree.pho3_MVA<-1.? -1.1: tree.pho3_MVA)
                p4mva[0] = float(tree.pho4_MVA<-1.? -1.1: tree.pho4_MVA)
                # p1_PSV[0] = float(tree.pho1_pixelseed)
                # p2_PSV[0] = float(tree.pho2_pixelseed)
                # p3_PSV[0] = float(tree.pho3_pixelseed)
                # p4_PSV[0] = float(tree.pho4_pixelseed)
                # a1_pt_Over_a1_mass[0] = float(float(tree.a1_pt_dM)/float(tree.a1_mass_dM))
                # a2_pt_Over_a2_mass[0] = float(float(tree.a2_pt_dM)/float(tree.a2_mass_dM))
                # a1_pt_Over_tp_mass[0] = float(float(tree.a1_pt_dM)/float(tree.tp_mass))
                # a2_pt_Over_tp_mass[0] = float(float(tree.a2_pt_dM)/float(tree.tp_mass))
                a1_a2_dR[0] = float(tree.a1_a2_dR_dM)
                bdt[0] = reader.EvaluateMVA("BDT")

                outtree.Fill()
            outfile.cd()
            if (sys_i == 0):
                outfile.mkdir("tagsDumper/trees")
            outfile.cd("tagsDumper/trees")
            outtree.Write()
    elif ('data_mix' in f):
        tree = infile.Get("Data_13TeV_H4GTag_0")
        nentries = tree.GetEntries()
        print "Data mix: ", nentries
        outtree = tree.CloneTree(0)
        bdt = array('f', [0])
        _bdt = outtree.Branch('bdt', bdt, 'bdt/F')
        # print "Gen mass:",float(opt.genMass)
        for i in range(0, nentries):
            if i%1000 == 0: print i
            tree.GetEntry(i)
            # if (training=='GenMass'):
            #     genMass[0] = float(opt.genMass)
            diffdM[0] = float((tree.a1_mass_dM-tree.a2_mass_dM))
            ct1[0] = float(tree.cosTheta_a1_dM)
            a1pt[0] = float(tree.a1_pt_dM)
            a2pt[0] = float(tree.a2_pt_dM)
            # p1mva[0] = float(tree.pho1_MVA)
            # p2mva[0] = float(tree.pho2_MVA)
            # p3mva[0] = float(tree.pho3_MVA)
            # p4mva[0] = float(tree.pho4_MVA)
            p1mva[0] = float(tree.pho1_MVA<-1.? -1.1: tree.pho1_MVA)
            p2mva[0] = float(tree.pho2_MVA<-1.? -1.1: tree.pho2_MVA)
            p3mva[0] = float(tree.pho3_MVA<-1.? -1.1: tree.pho3_MVA)
            p4mva[0] = float(tree.pho4_MVA<-1.? -1.1: tree.pho4_MVA)
            # p1_PSV[0] = float(tree.pho1_pixelseed)
            # p2_PSV[0] = float(tree.pho2_pixelseed)
            # p3_PSV[0] = float(tree.pho3_pixelseed)
            # p4_PSV[0] = float(tree.pho4_pixelseed)
            # a1_pt_Over_a1_mass[0] = float(float(tree.a1_pt_dM)/float(tree.a1_mass_dM))
            # a2_pt_Over_a2_mass[0] = float(float(tree.a2_pt_dM)/float(tree.a2_mass_dM))
            # a1_pt_Over_tp_mass[0] = float(float(tree.a1_pt_dM)/float(tree.tp_mass))
            # a2_pt_Over_tp_mass[0] = float(float(tree.a2_pt_dM)/float(tree.tp_mass))
            a1_a2_dR[0] = float(tree.a1_a2_dR_dM)
            bdt[0] = reader.EvaluateMVA("BDT")
            # print "bdt: ", bdt[0]

            outtree.Fill()
        outfile.cd()
        outfile.mkdir("tagsDumper/trees")
        outfile.cd("tagsDumper/trees")
        outtree.Write()
    else:
        tree = infile.Get("tagsDumper/trees/Data_13TeV_H4GTag_0")
        # tree = infile.Get("Data_13TeV_H4GTag_0")
        nentries = tree.GetEntries()
        print "Data: ", nentries
        outtree = tree.CloneTree(0)
        bdt = array('f', [0])
        _bdt = outtree.Branch('bdt', bdt, 'bdt/F')
        # print "Gen mass:",float(opt.genMass)
        for i in range(0, nentries):
            if i%1000 == 0: print i
            tree.GetEntry(i)
            # if (training=='GenMass'):
            #     genMass[0] = float(opt.genMass)
            diffdM[0] = float((tree.a1_mass_dM-tree.a2_mass_dM))
            ct1[0] = float(tree.cosTheta_a1_dM)
            a1pt[0] = float(tree.a1_pt_dM)
            a2pt[0] = float(tree.a2_pt_dM)
            # p1mva[0] = float(tree.pho1_MVA)
            # p2mva[0] = float(tree.pho2_MVA)
            # p3mva[0] = float(tree.pho3_MVA)
            # p4mva[0] = float(tree.pho4_MVA)
            p1mva[0] = float(tree.pho1_MVA<-1.? -1.1: tree.pho1_MVA)
            p2mva[0] = float(tree.pho2_MVA<-1.? -1.1: tree.pho2_MVA)
            p3mva[0] = float(tree.pho3_MVA<-1.? -1.1: tree.pho3_MVA)
            p4mva[0] = float(tree.pho4_MVA<-1.? -1.1: tree.pho4_MVA)
            # p1_PSV[0] = float(tree.pho1_pixelseed)
            # p2_PSV[0] = float(tree.pho2_pixelseed)
            # p3_PSV[0] = float(tree.pho3_pixelseed)
            # p4_PSV[0] = float(tree.pho4_pixelseed)
            # a1_pt_Over_a1_mass[0] = float(float(tree.a1_pt_dM)/float(tree.a1_mass_dM))
            # a2_pt_Over_a2_mass[0] = float(float(tree.a2_pt_dM)/float(tree.a2_mass_dM))
            # a1_pt_Over_tp_mass[0] = float(float(tree.a1_pt_dM)/float(tree.tp_mass))
            # a2_pt_Over_tp_mass[0] = float(float(tree.a2_pt_dM)/float(tree.tp_mass))
            a1_a2_dR[0] = float(tree.a1_a2_dR_dM)
            bdt[0] = reader.EvaluateMVA("BDT")

            outtree.Fill()
        outfile.cd()
        outfile.mkdir("tagsDumper/trees")
        outfile.cd("tagsDumper/trees")
        outtree.Write()


    outfile.Write()
