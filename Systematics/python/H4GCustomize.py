import FWCore.ParameterSet.Config as cms

class H4GCustomize():
    """
    HH->WWgg process customizaton class. Started with HH->bbgg customization class.
    """

    def __init__(self, process, customize, metaConditions):
        self.process = process
        self.customize = customize
        self.metaConditions = metaConditions
        # self.tagList = [ ["HHWWggTag",12] ]
        self.tagList = [ ["H4GTag",1] ]
        self.customizeTagSequence()

    def variablesToDump(self):
        var_workspace = [
            # "testVariable[10,0,10] := 5 "
#             "Mjj := dijet().M()"
            # "eventNumber := eventNumber()",
            # "MX := MX()",
            # "HHbbggMVA := MVA()"
        ]
        variables = [
            # "dipho_pt := dipho.pt()"
            # "pho2_MVA := slp_Hgg_MVA"
            # Cut flow variables
            # "passPS[2,0,2] := Cut_Variables[0]",
            # "passPhotonSels[2,0,2] := Cut_Variables[1]",
            # "passbVeto[2,0,2] := Cut_Variables[2]",
            # "ExOneLep[2,0,2] := Cut_Variables[3]",
            # "goodJets[2,0,2] := Cut_Variables[4]",
            # "lp_E[100,0,100] := Leading_Photon.p4().E()",
            # "slp_E[100,0,100] := Subleading_Photon.p4().E()",
            # "lp_initE[100,0,100] := Leading_Photon.energyAtStep('initial')",
            # "slp_initE[100,0,100] := Subleading_Photon.energyAtStep('initial')",
            # "lp_Hgg_MVA[100,-1,1] := lp_Hgg_MVA()",
            # "slp_Hgg_MVA[100,-1,1] := slp_Hgg_MVA()"


            # also want final energies
            # "leadingJet_bDis := leadJet().bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",#FIXME make the btag type configurable?
            # "subleadingJet_bDis := subleadJet().bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",
            # "leadingJet_DeepCSV := leadJet().bDiscriminator('pfDeepCSVJetTags:probb')+leadJet().bDiscriminator('pfDeepCSVJetTags:probbb')",#FIXME make the btag type configurable?
            # "subleadingJet_DeepCSV := subleadJet().bDiscriminator('pfDeepCSVJetTags:probb')+subleadJet().bDiscriminator('pfDeepCSVJetTags:probbb')",
            # "leadingJet_DeepFlavour := leadJet().bDiscriminator('mini_pfDeepFlavourJetTags:probb')+leadJet().bDiscriminator('mini_pfDeepFlavourJetTags:probbb')+leadJet().bDiscriminator('mini_pfDeepFlavourJetTags:problepb')",#FIXME make the btag type configurable?
            # "subleadingJet_DeepFlavour := subleadJet().bDiscriminator('mini_pfDeepFlavourJetTags:probb')+subleadJet().bDiscriminator('mini_pfDeepFlavourJetTags:probbb')+subleadJet().bDiscriminator('mini_pfDeepFlavourJetTags:problepb')",#FIXME make the btag type configurable?
            # "leadingJet_puJetIdMVA := leadJet().puJetIdMVA()",
            # "subleadingJet_puJetIdMVA := subleadJet().puJetIdMVA()",
            # "leadingJet_puJetIdMVA := leadJet().puJetIdMVA()",
            # "subleadingJet_puJetIdMVA := subleadJet().puJetIdMVA()",
            # "absCosThetaStar_CS := abs(getCosThetaStar_CS())",
            # "absCosThetaStar_CS_old := abs(getCosThetaStar_CS_old(6500))",
            # "absCosTheta_bb := abs(CosThetaAngles()[1])",
            # "absCosTheta_gg := abs(CosThetaAngles()[0])",
            # "diphotonCandidatePtOverdiHiggsM := diphotonPtOverM()",
            # "dijetCandidatePtOverdiHiggsM := dijetPtOverM()",
            # "customLeadingPhotonIDMVA := diPhoton.leadingView.phoIdMvaWrtChosenVtx",
            # "customSubLeadingPhotonIDMVA := diPhoton.subLeadingView.phoIdMvaWrtChosenVtx",
            # "EGMLeadingPhotonIDMVA := diPhoton.leadingPhoton.userFloat('EGMPhotonMVA')",
            # "EGMSubLeadingPhotonIDMVA := diPhoton.subLeadingPhoton.userFloat('EGMPhotonMVA')",
            # "leadingPhotonSigOverE := diPhoton.leadingPhoton.sigEOverE",
            # "subleadingPhotonSigOverE := diPhoton.subLeadingPhoton.sigEOverE",
            # "sigmaMOverM := sqrt(0.5*(diPhoton.leadingPhoton.sigEOverE*diPhoton.leadingPhoton.sigEOverE + diPhoton.subLeadingPhoton.sigEOverE*diPhoton.subLeadingPhoton.sigEOverE))",
            # "sigmaMOverMDecorr := getSigmaMDecorr()",
            # "PhoJetMinDr := getPhoJetMinDr()",#up to here input variables to MVA
            # # "HHbbggMVA := MVA()",
            # # "HHbbggMVAprob0 := MVAprob()[0]",
            # "MX := MX()",
            # "genMhh := genMhh()",
            # "Mjj := dijet().M()",
            # "dijet_pt := dijet().pt",
            # "dijet_eta := dijet().eta",
            # "dijet_phi := dijet().phi",
            # "diphoton_pt := diPhoton.pt",
            # "diphoton_eta := diPhoton.eta",
            # "diphoton_phi := diPhoton.phi",
            # "leadingJet_btagWeight := leadJet.weight('JetBTagReshapeWeight') ",
            # "subleadingJet_btagWeight := subleadJet.weight('JetBTagReshapeWeight') ",

            # "diHiggs_pt := getdiHiggsP4().pt()",
            # "diHiggs_mass := getdiHiggsP4().M()",
            # "diHiggs_eta :=  getdiHiggsP4().eta()",
            # "diHiggs_phi := getdiHiggsP4().phi()",
            # "category := categoryNumber()",

            # "leadingPhoton_pt := diPhoton.leadingPhoton.pt",
            # "leadingPhoton_eta := diPhoton.leadingPhoton.eta",
            # "leadingPhoton_phi := diPhoton.leadingPhoton.phi",
            # "subleadingPhoton_pt := diPhoton.subLeadingPhoton.pt",
            # "subleadingPhoton_eta := diPhoton.subLeadingPhoton.eta",
            # "subleadingPhoton_phi := diPhoton.subLeadingPhoton.phi",

            # "leadingJet_pt := leadJet().pt",
            # "leadingJet_eta := leadJet().eta",
            # "leadingJet_phi := leadJet().phi",
            # "leadingJet_mass := leadJet().p4().M()",
            # "leadingJet_hflav := leadJet().hadronFlavour()",
            # "leadingJet_pflav := leadJet().partonFlavour()",

            # "subleadingJet_pt := subleadJet().pt",
            # "subleadingJet_eta := subleadJet().eta",
            # "subleadingJet_phi := subleadJet().phi",
            # "subleadingJet_mass := subleadJet().p4().M()",
            # "subleadingJet_hflav := subleadJet().hadronFlavour()",
            # "subleadingJet_pflav := subleadJet().partonFlavour()",
        ]

        # Save b scores
        # for jeti in range(0,6):
        #     var1 = "jet" + str(jeti) + "_DeepFlavourScore[2,0,2] := ? JetVector.size() >= " + str(jeti + 1) + " ? JetVector[" + str(jeti) + "].bDiscriminator('mini_pfDeepFlavourJetTags:probb') : -99 "
        #     var2 = "jet" + str(jeti) + "_DeepCSVScore[2,0,2] := ? JetVector.size() >= " + str(jeti + 1) + " ? JetVector[" + str(jeti) + "].bDiscriminator('pfDeepCSVJetTags:probb') : -99 "
        #     variables.append(var1)
        #     variables.append(var2)

        if self.customize.doBJetRegression : variables +=[
                # "leadingJet_bRegNNCorr := leadJet().userFloat('bRegNNCorr')",
                # "leadingJet_bRegNNResolution := leadJet().userFloat('bRegNNResolution')",
                # "subleadingJet_bRegNNCorr := subleadJet().userFloat('bRegNNCorr')",
                # "subleadingJet_bRegNNResolution := subleadJet().userFloat('bRegNNResolution')",
                # "sigmaMJets := getSigmaMOverMJets()"
        ]
        # if self.customize.HHWWggReweight > 0:
        #     for num in range(0,12):  #12 benchmarks + 1 SM
        #          variables += ["benchmark_reweight_%d := getBenchmarkReweight(%d)"%(num,num)]
        #          var_workspace += ["benchmark_reweight_%d := getBenchmarkReweight(%d)"%(num,num)]
        #     variables += ["benchmark_reweight_SM := getBenchmarkReweight(12)"]
        #     variables += ["benchmark_reweight_box := getBenchmarkReweight(13)"]
        #     variables += ["benchmark_reweight_2017fake := getBenchmarkReweight(14)"]
        #     var_workspace += ["benchmark_reweight_SM := getBenchmarkReweight(12)"]
        #     var_workspace += ["benchmark_reweight_box := getBenchmarkReweight(13)"]
        #     var_workspace += ["benchmark_reweight_2017fake := getBenchmarkReweight(14)"]

        if self.customize.ttHKillerSaveInputVariables : variables += [
            # "ttH_sumET := sumET()",
            # "ttH_MET := MET()",
            # "ttH_phiMET := phiMET()",
            # "ttH_dPhi1 := dPhi1()",
            # "ttH_dPhi2 := dPhi2()",
            # "ttH_PhoJetMinDr := PhoJetMinDr()",
            # "ttH_njets := njets()",
            # "ttH_Xtt0 := Xtt0()",
            # "ttH_Xtt1 := Xtt1()",
            # "ttH_pte1 := pte1()",
            # "ttH_pte2 := pte2()",
            # "ttH_ptmu1 := ptmu1()",
            # "ttH_ptmu2 := ptmu2()",
            # "ttH_ptdipho := ptdipho()",
            # "ttH_etae1 := etae1()",
            # "ttH_etae2 := etae2()",
            # "ttH_etamu1 := etamu1()",
            # "ttH_etamu2 := etamu2()",
            # "ttH_etadipho := etadipho()",
            # "ttH_phie1 := phie1()",
            # "ttH_phie2 := phie2()",
            # "ttH_phimu1 := phimu1()",
            # "ttH_phimu2 := phimu2()",
            # "ttH_phidipho := phidipho()",
            # "ttH_fabs_CosThetaStar_CS := fabs_CosThetaStar_CS()",
            # "ttH_fabs_CosTheta_bb := fabs_CosTheta_bb()",
            # "ttH_ptjet1 := ptjet1()",
            # "ttH_ptjet2 := ptjet2()",
            # "ttH_etajet1 := etajet1()",
            # "ttH_etajet2 := etajet2()",
            # "ttH_phijet1 := phijet1()",
            # "ttH_phijet2 := phijet2()"
            ]


        # if self.customize.doHHWWggttHKiller :
        #      variables +=[
        #        "ttHScore := ttHScore()",
        #      ]
        #      var_workspace +=[
        #        "ttHScore := ttHScore()",
        #      ]
        return variables

        # if self.customize.dumpWorkspace == False :
        #     return variables
        # else :
        #     return var_workspace


    def systematicVariables(self):
    #   systematicVariables=["CMS_hgg_mass[160,100,180]:=diPhoton().mass","Mjj[120,70,190]:=dijet().M()","HHbbggMVA[100,0,1.]:=MVA()","MX[300,250,5000]:=MX()"]
      systematicVariables=[
          # "CMS_hgg_mass[160,100,180]:=diPhoton().mass",
          # "lp_E[100,0,100] := Leading_Photon.p4().E()",
          # "slp_E[100,0,100] := Subleading_Photon.p4().E()",
          # "lp_initE[100,0,100] := Leading_Photon.energyAtStep('initial')",
          # "slp_initE[100,0,100] := Subleading_Photon.energyAtStep('initial')", # also want final energies
        #  "jet0_btag                       := ? JetVector.size() >= 1 ? JetVector[0].bDiscriminator('mini_pfDeepFlavourJetTags:probb') : -99 ",
      ]
    #   systematicVariables=[]

    #   if self.customize.HHWWggReweight > 0:
    #      for num in range(0,12):  #12 benchmarks
    #         systematicVariables += ["benchmark_reweight_%d[100,0,200] := getBenchmarkReweight(%d)"%(num,num)]
    #      systematicVariables+= ["benchmark_reweight_SM[100,0,200] := getBenchmarkReweight(12)"]
    #      systematicVariables+= ["benchmark_reweight_box[100,0,200] := getBenchmarkReweight(13)"]

    #   if self.customize.doHHWWggttHKiller :
    #          systematicVariables +=["ttHScore := ttHScore()"]

      return systematicVariables


    def variablesToDumpData():
        variables = [
            # "testVariable[100,0,100] := 50"
            # "jet0_btag[2,0,2]                       := ? JetVector.size() >= 1 ? JetVector[0].bDiscriminator('mini_pfDeepFlavourJetTags:probb') : -99 ",
            # "TestVariable:=111"
           #  "leadingJet_DeepCSV := leadJet().bDiscriminator('pfDeepCSVJetTags:probb')+leadJet().bDiscriminator('pfDeepCSVJetTags:probbb')",#FIXME make the btag type configurable?
           #  "subleadingJet_DeepCSV := subleadJet().bDiscriminator('pfDeepCSVJetTags:probb')+subleadJet().bDiscriminator('pfDeepCSVJetTags:probbb')",
           #  "absCosThetaStar_CS := abs(getCosThetaStar_CS())",
           #  "absCosThetaStar_CS_old := abs(getCosThetaStar_CS_old(6500))",
           #  "absCosTheta_bb := abs(CosThetaAngles()[1])",
           #  "absCosTheta_gg := abs(CosThetaAngles()[0])",
           #  "diphotonCandidatePtOverdiHiggsM := diphotonPtOverM()",
           #  "dijetCandidatePtOverdiHiggsM := dijetPtOverM()",
           #  "customLeadingPhotonIDMVA := diPhoton.leadingView.phoIdMvaWrtChosenVtx",
           #  "customSubLeadingPhotonIDMVA := diPhoton.subLeadingView.phoIdMvaWrtChosenVtx",
           #  "leadingPhotonSigOverE := diPhoton.leadingPhoton.sigEOverE",
           #  "subleadingPhotonSigOverE := diPhoton.subLeadingPhoton.sigEOverE",
           #  "sigmaMOverM := sqrt(0.5*(diPhoton.leadingPhoton.sigEOverE*diPhoton.leadingPhoton.sigEOverE + diPhoton.subLeadingPhoton.sigEOverE*diPhoton.subLeadingPhoton.sigEOverE))",
           #  "PhoJetMinDr := getPhoJetMinDr()",#up to here input variables to MVA
           #  "leadingJet_bRegNNResolution := leadJet().userFloat('bRegNNResolution')",
           #  "subleadingJet_bRegNNResolution := subleadJet().userFloat('bRegNNResolution')",
           #  "sigmaMJets := getSigmaMOverMJets()",
            #  "HHbbggMVA := MVA()",
            #  "MX := MX()",
           #  "Mjj := dijet().M()",
           #  "eventNumber := eventNumber()",
             ]

        # for jeti in range(0,6):
        #     var = "jet" + str(jeti) + "_btag[2,0,2] := ? JetVector.size() >= " + str(jeti + 1) + " ? JetVector[" + str(jeti) + "].bDiscriminator('mini_pfDeepFlavourJetTags:probb') : -99 "
        #     variables.append(var)

        # if self.customize.doHHWWggttHKiller : variables +=[
        #     "ttHScore := ttHScore()",
        #    ]
        return variables



    def customizeTagSequence(self):
        self.process.load("flashgg.Taggers.flashggH4GTag_cff")

        if self.customize.doHHWWggTagCutFlow:
            self.process.flashggH4GTag.doHHWWggTagCutFlowAnalysis = cms.bool(True)

        ## customize meta conditions
        # self.process.flashggHHWWggTag.JetIDLevel=cms.string(str(self.metaConditions["doubleHTag"]["jetID"]))
        # self.process.flashggHHWWggTag.MVAConfig.weights=cms.FileInPath(str(self.metaConditions["doubleHTag"]["weightsFile"]))
        # self.process.flashggHHWWggTag.MVAscaling = cms.double(self.metaConditions["doubleHTag"]["MVAscalingValue"])
        # self.process.flashggHHWWggTag.MVAFlatteningFileName = cms.untracked.FileInPath(str(self.metaConditions["doubleHTag"]["MVAFlatteningFileName"]))
        # self.process.flashggHHWWggTag.dottHTagger = cms.bool(self.customize.doHHWWggttHKiller)
        # self.process.flashggHHWWggTag.ttHWeightfile = cms.untracked.FileInPath(str(self.metaConditions["doubleHTag"]["ttHWeightfile"]))
        # self.process.flashggHHWWggTag.ttHKiller_mean = cms.vdouble(self.metaConditions["doubleHTag"]["ttHKiller_mean"])
        # self.process.flashggHHWWggTag.ttHKiller_std = cms.vdouble(self.metaConditions["doubleHTag"]["ttHKiller_std"])
        # self.process.flashggHHWWggTag.ttHKiller_listmean = cms.vdouble(self.metaConditions["doubleHTag"]["ttHKiller_listmean"])
        # self.process.flashggHHWWggTag.ttHKiller_liststd = cms.vdouble(self.metaConditions["doubleHTag"]["ttHKiller_liststd"])

        ## remove single Higgs tags

        print'Removing single Higgs tags'

        if self.customize.H4GTagsOnly:
            self.process.flashggTagSequence.remove(self.process.flashggVBFTag)
            self.process.flashggTagSequence.remove(self.process.flashggTTHLeptonicTag)
            self.process.flashggTagSequence.remove(self.process.flashggTTHHadronicTag)
            self.process.flashggTagSequence.remove(self.process.flashggVHEtTag)
            self.process.flashggTagSequence.remove(self.process.flashggVHLooseTag)
            self.process.flashggTagSequence.remove(self.process.flashggVHTightTag)
            self.process.flashggTagSequence.remove(self.process.flashggVHMetTag)
            self.process.flashggTagSequence.remove(self.process.flashggWHLeptonicTag)
            self.process.flashggTagSequence.remove(self.process.flashggZHLeptonicTag)
            self.process.flashggTagSequence.remove(self.process.flashggVHLeptonicLooseTag)
            self.process.flashggTagSequence.remove(self.process.flashggVHHadronicTag)
            self.process.flashggTagSequence.remove(self.process.flashggVBFMVA)
            self.process.flashggTagSequence.remove(self.process.flashggVBFDiPhoDiJetMVA)
            self.process.flashggTagSequence.remove(self.process.flashggTTHDiLeptonTag)
            self.process.flashggTagSequence.remove(self.process.flashggUntagged)
            self.process.flashggTagSequence.remove(self.process.flashggUntagged)
            self.process.flashggTagSequence.remove(self.process.flashggTHQLeptonicTag)

        self.process.flashggTagSequence.replace(self.process.flashggTagSorter,self.process.flashggH4GTagSequence*self.process.flashggTagSorter)
        self.process.flashggTagSorter.TagPriorityRanges = cms.VPSet( cms.PSet(TagName = cms.InputTag('flashggH4GTag')) )

    def H4GTagMerger(self,systlabels=[]):
        self.process.p.remove(self.process.flashggTagSorter)
        self.process.p.replace(self.process.flashggSystTagMerger,self.process.flashggH4GTagSequence*self.process.flashggTagSorter*self.process.flashggSystTagMerger)
        # print'process.p = ',self.process.p

        ##-- Do I need this part for HHWWgg?
        for systlabel in systlabels:
           if systlabel!='':
             self.process.p.remove(getattr(self.process,'flashggTagSorter'+systlabel))
             self.process.p.replace(self.process.flashggSystTagMerger,getattr(self.process, 'flashggTagSorter'+systlabel)*self.process.flashggSystTagMerger)
           setattr(getattr(self.process, 'flashggTagSorter'+systlabel), 'TagPriorityRanges', cms.VPSet( cms.PSet(TagName = cms.InputTag('flashggH4GTag')) ))
        #    setattr(getattr(self.process, 'flashggTagSorter'+systlabel), 'TagPriorityRanges', cms.VPSet( cms.PSet(TagName = cms.InputTag('flashggHHWWggTag', systlabel)) ))

        # print 'from loop after:',process.flashggSystTagMerger.src


    def H4GTagRunSequence(self,systlabels,jetsystlabels,phosystlabels):
        print'not used'
    #    if self.customize.HHWWggTagsOnly:
        #   print'systlabels = ',systlabels
        #   self.HHWWggTagMerger(systlabels)


        ## Not sure if this is necessary for HHWWgg
    #    if len(systlabels)>1 :
    #       print'[HHWWggTagRunSequence] - Add JetesSuffixes and diphotonsuffices'
    #       getattr(self.process, "flashggHHWWggTag").JetsSuffixes = cms.vstring([systlabels[0]]+jetsystlabels)
    #       getattr(self.process, "flashggHHWWggTag").DiPhotonSuffixes = cms.vstring([systlabels[0]]+phosystlabels)





    #    if self.customize.HHWWggReweight>0:
        #   self.addNodesReweighting()

    #    if self.customize.doHHWWggGenAnalysis:
        #   self.addGenAnalysis()



    def addNodesReweighting(self):
        print'[addNodesReweighting]: Doing Nothing for HHWWgg'
        # if self.customize.HHWWggReweight > 0 :
        #     from flashgg.Taggers.flashggHHWWggReweight_cfi import flashggHHWWggReweight
        #     self.process.flashggHHWWggReweight = flashggHHWWggReweight
        #     self.process.flashggHHWWggReweight.doReweight = self.customize.HHWWggReweight
        #     self.process.flashggHHWWggReweight.weightsFile = cms.untracked.FileInPath(str(self.metaConditions["HHWWggTag"]["NodesReweightingFileName"]))
        #     self.process.p.replace(self.process.flashggHHWWggTag, self.process.flashggHHWWggReweight*self.process.flashggHHWWggTag)


    def addGenAnalysis(self):
        if self.customize.processId == "Data":
            return

        import flashgg.Taggers.dumperConfigTools as cfgTools
        ## load gen-level bbgg
        self.process.load( "flashgg.MicroAOD.flashggGenDiPhotonDiBJetsSequence_cff" )

        ## match gen-level to reco tag
        self.process.load("flashgg.Taggers.flashggTaggedGenDiphotons_cfi")
        self.process.flashggTaggedGenDiphotons.src  = "flashggSelectedGenDiPhotonDiBJets"
        self.process.flashggTaggedGenDiphotons.tags = "flashggTagSorter"
        self.process.flashggTaggedGenDiphotons.remap = self.process.tagsDumper.classifierCfg.remap

        ## prepare gen-level dumper
        self.process.load("flashgg.Taggers.genDiphotonDumper_cfi")
        self.process.genDiphotonDumper.dumpTrees = True
        self.process.genDiphotonDumper.dumpWorkspace = False
        self.process.genDiphotonDumper.src = "flashggTaggedGenDiphotons"

        from flashgg.Taggers.globalVariables_cff import globalVariables
        self.process.genDiphotonDumper.dumpGlobalVariables = True
        self.process.genDiphotonDumper.globalVariables = globalVariables

        genVariables = ["mgg := mass",
                        "mbb := dijet.mass",
                        "mhh := sqrt( pow(energy+dijet.energy,2) - pow(px+dijet.px,2) - pow(py+dijet.py,2) - pow(pz+dijet.pz,2))",


                        "leadPho_px := leadingPhoton.px",
                        "leadPho_py := leadingPhoton.py",
                        "leadPho_pz := leadingPhoton.pz",
                        "leadPho_e  := leadingPhoton.energy",
                        "subleadPho_px := subLeadingPhoton.px",
                        "subleadPho_py := subLeadingPhoton.py",
                        "subleadPho_pz := subLeadingPhoton.pz",
                        "subleadPho_e  := subLeadingPhoton.energy",

                        "leadJet_px := leadingJet.px",
                        "leadJet_py := leadingJet.py",
                        "leadJet_pz := leadingJet.pz",
                        "leadJet_e  := leadingJet.energy",
                        "subleadJet_px := subLeadingJet.px",
                        "subleadJet_py := subLeadingJet.py",
                        "subleadJet_pz := subLeadingJet.pz",
                        "subleadJet_e  := subLeadingJet.energy",

                        ]
        # if self.customize.HHWWggReweight > 0:
        #      for num in range(0,12):
        #            genVariables += ["benchmark_reweight_%d := getHHbbggBenchmarkReweight(%d)"%(num,num)]
        #      genVariables += ["benchmark_reweight_SM := getHHbbggBenchmarkReweight(12)"]
        #      genVariables += ["benchmark_reweight_box := getHHbbggBenchmarkReweight(13)"]
        #      genVariables += ["benchmark_reweight_2017fake := getHHbbggBenchmarkReweight(14)"]

        ## define categories for gen-level dumper
        cfgTools.addCategory(self.process.genDiphotonDumper,  ## events with not reco-level tag
                             "NoTag", 'isTagged("flashggNoTag")',1,
                             variables=genVariables,
                             )

        for tag in self.tagList: ## tagged events
            tagName,subCats = tag
            # need to define all categories explicitely because cut-based classifiers do not look at sub-category number
            for isub in xrange(subCats):
                cfgTools.addCategory(self.process.genDiphotonDumper,
                                     "%s_%d" % ( tagName, isub ),
                                     'isTagged("%s") && categoryNumber == %d' % (tagName, isub),0,
                                     variables=genVariables##+recoVariables
                                     )

        self.process.genp = cms.Path(self.process.flashggGenDiPhotonDiBJetsSequence*self.process.flashggTaggedGenDiphotons*self.process.genDiphotonDumper)
