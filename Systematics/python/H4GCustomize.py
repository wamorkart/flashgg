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
        self.tagList = [ ["H4GTag",3] ]
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
            "pho1_MVA  := pho1_MVA ",
            "pho2_MVA  := pho2_MVA ",
            "pho3_MVA  := pho3_MVA ",
            "pho4_MVA  := pho4_MVA ",
            "pho1_pt   := pho1.pt()",
            "pho2_pt   := pho2.pt()",
            "pho3_pt   := pho3.pt()",
            "pho4_pt   := pho4.pt()",
            # "pho2_pt   := phoP4Corrected_dp[1].pt()",
            # "pho3_pt   := phoP4Corrected_dp[2].pt()",
            # "pho4_pt   := phoP4Corrected_dp[3].pt()",
            "pho1_energy := pho1.energy()",
            # "pho2_energy := phoP4Corrected_dp[1].energy()",
            # "pho3_energy := phoP4Corrected_dp[2].energy()",
            # "pho4_energy := phoP4Corrected_dp[3].energy()",
            "pho1_energy_init := pho1.energyAtStep('initial')",
            "dZ_bdtVtx        := dZ_bdtVtx"
            # "pho2_energy_init := phoP4Corrected_dp[1].energyAtStep('initial')",
            # "pho3_energy_init := phoP4Corrected_dp[2].energyAtStep('initial')",
            # "pho4_energy_init := phoP4Corrected_dp[3].energyAtStep('initial')",
            # "pho1_eta   := phoP4Corrected_dp[0].eta()",
            # "pho2_eta   := phoP4Corrected_dp[1].eta()",
            # "pho3_eta   := phoP4Corrected_dp[2].eta()",
            # "pho4_eta   := phoP4Corrected_dp[3].eta()",
            # "pho1_SC_eta   := phoP4Corrected_dp[0].superCluster().eta()",
            # "pho2_SC_eta   := phoP4Corrected_dp[1].superCluster().eta()",
            # "pho3_SC_eta   := phoP4Corrected_dp[2].superCluster().eta()",
            # "pho4_SC_eta   := phoP4Corrected_dp[3].superCluster().eta()",
            # "pho1_r9   := phoP4Corrected_dp[0].old_r9()",
            # "pho2_r9   := phoP4Corrected_dp[1].old_r9()",
            # "pho3_r9   := phoP4Corrected_dp[2].old_r9()",
            # "pho4_r9   := phoP4Corrected_dp[3].old_r9()",
            # "pho1_full5x5_r9   := phoP4Corrected_dp[0].full5x5_r9()",
            # "pho2_full5x5_r9   := phoP4Corrected_dp[1].full5x5_r9()",
            # "pho3_full5x5_r9   := phoP4Corrected_dp[2].full5x5_r9()",
            # "pho4_full5x5_r9   := phoP4Corrected_dp[3].full5x5_r9()",
            # "pho1_match   := phoP4Corrected_dp[0].genMatchType()",
            # "pho2_match   := phoP4Corrected_dp[1].genMatchType()",
            # "pho3_match   := phoP4Corrected_dp[2].genMatchType()",
            # "pho4_match   := phoP4Corrected_dp[3].genMatchType()",
            # "pho1_PSV   := phoP4Corrected_dp[0].hasPixelSeed()",
            # "pho2_PSV   := phoP4Corrected_dp[1].hasPixelSeed()",
            # "pho3_PSV   := phoP4Corrected_dp[2].hasPixelSeed()",
            # "pho4_PSV   := phoP4Corrected_dp[3].hasPixelSeed()",
            # "pho1_electronVeto   := phoP4Corrected_dp[0].passElectronVeto()",
            # "pho2_electronVeto   := phoP4Corrected_dp[1].passElectronVeto()",
            # "pho3_electronVeto   := phoP4Corrected_dp[2].passElectronVeto()",
            # "pho4_electronVeto   := phoP4Corrected_dp[3].passElectronVeto()",
            # "pho12_dR := deltaR( phoP4Corrected_dp[0].eta, phoP4Corrected_dp[0].phi, phoP4Corrected_dp[1].eta, phoP4Corrected_dp[1].phi )",
            # "pho13_dR := deltaR( phoP4Corrected_dp[0].eta, phoP4Corrected_dp[0].phi, phoP4Corrected_dp[2].eta, phoP4Corrected_dp[2].phi )",
            # "pho14_dR := deltaR( phoP4Corrected_dp[0].eta, phoP4Corrected_dp[0].phi, phoP4Corrected_dp[3].eta, phoP4Corrected_dp[3].phi )",
            # "pho23_dR := deltaR( phoP4Corrected_dp[1].eta, phoP4Corrected_dp[1].phi, phoP4Corrected_dp[2].eta, phoP4Corrected_dp[2].phi )",
            # "pho24_dR := deltaR( phoP4Corrected_dp[1].eta, phoP4Corrected_dp[1].phi, phoP4Corrected_dp[3].eta, phoP4Corrected_dp[3].phi )",
            # "pho34_dR := deltaR( phoP4Corrected_dp[2].eta, phoP4Corrected_dp[2].phi, phoP4Corrected_dp[3].eta, phoP4Corrected_dp[3].phi )",
            # "pho12_m :=  pho12.mass()",
            # "pho13_m :=  pho13.mass()",
            # "pho14_m :=  pho14.mass()",
            # "pho23_m :=  pho23.mass()",
            # "pho24_m :=  pho24.mass()",
            # "pho34_m :=  pho34.mass()",
            # "a1_mass := h4gDiPho1.mass()",
            # "a2_mass := h4gDiPho2.mass()",
            # "avg_a_mass := (h4gDiPho1.mass() + h4gDiPho2.mass())/2",
            # "a1_pt := h4gDiPho1.pt()",
            # "a2_pt := h4gDiPho2.pt()",
            # "a1_eta := h4gDiPho1.eta()",
            # "a2_eta := h4gDiPho2.eta()",
            # "a1_phi := h4gDiPho1.phi()",
            # "a2_phi := h4gDiPho2.phi()",
            # "a1_energy := h4gDiPho1.energy()",
            # "a2_energy := h4gDiPho2.energy()",
            # "a1_dR := deltaR(h4gDiPho1_Pho1.eta(), h4gDiPho1_Pho1.phi(), h4gDiPho1_Pho2.eta(), h4gDiPho1_Pho2.phi())",
            # "a2_dR := deltaR(h4gDiPho2_Pho1.eta(), h4gDiPho2_Pho1.phi(), h4gDiPho2_Pho2.eta(), h4gDiPho2_Pho2.phi())",
            # "a1_a2_dR := deltaR(h4gDiPho1.eta(),h4gDiPho1.phi(),h4gDiPho2.eta(),h4gDiPho2.phi())",
            # "a1_pToverMass := h4gDiPho1.pt()/h4gDiPho1.mass()",
            # "a2_pToverMass := h4gDiPho2.pt()/h4gDiPho2.mass()",
            # "tp_mass := tp.mass()",
            # "tp_pt := tp.pt()",
            # "tp_eta := tp.eta()",
            # "tp_phi := tp.phi()",
            # "tp_energy := tp.energy()",
            # "absCosThetaStar_CS      := getCosThetaStar_CS()",
            # "absCosTheta_pho_a1      := CosThetaAngles()[0]",
            # "absCosTheta_pho_a2      := CosThetaAngles()[1]",
            # "a1_mass_prime := h4gDiPho1_prime.mass()",
            # "a2_mass_prime := h4gDiPho2_prime.mass()",
            # "avg_a_mass_prime := (h4gDiPho1_prime.mass() + h4gDiPho2_prime.mass())/2",
            # "a1_pt_prime := h4gDiPho1_prime.pt()",
            # "a2_pt_prime := h4gDiPho2_prime.pt()",
            # "a1_eta_prime := h4gDiPho1_prime.eta()",
            # "a2_eta_prime := h4gDiPho2_prime.eta()",
            # "a1_phi_prime := h4gDiPho1_prime.phi()",
            # "a2_phi_prime := h4gDiPho2_prime.phi()",
            # "a1_energy_prime := h4gDiPho1_prime.energy()",
            # "a2_energy_prime := h4gDiPho2_prime.energy()",
            # "a1_dR_prime := deltaR(h4gDiPho1_Pho1_prime.eta(), h4gDiPho1_Pho1_prime.phi(), h4gDiPho1_Pho2_prime.eta(), h4gDiPho1_Pho2_prime.phi())",
            # "a2_dR_prime := deltaR(h4gDiPho2_Pho1_prime.eta(), h4gDiPho2_Pho1_prime.phi(), h4gDiPho2_Pho2_prime.eta(), h4gDiPho2_Pho2_prime.phi())",
            # "a1_a2_dR_prime := deltaR(h4gDiPho1_prime.eta(),h4gDiPho1_prime.phi(),h4gDiPho2_prime.eta(),h4gDiPho2_prime.phi())",
            # "absCosThetaStar_CS_prime      := getCosThetaStar_CS_prime()",
            # "absCosTheta_pho_a1_prime      := CosThetaAngles_prime()[0]",
            # "absCosTheta_pho_a2_prime      := CosThetaAngles_prime()[1]",
            # "dZ_bdtVtx               := genVertex.z() - vertex_bdt.z()",
            # "dZ_hggVtx               := genVertex.z() - vertex_diphoton.z()",
            # "dZ_zeroVtx              := genVertex.z() - vertex_zero.z()"
        ]


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
          "pho1_MVA  := pho1_MVA ",
          "pho2_MVA  := pho2_MVA ",
          "pho3_MVA  := pho3_MVA ",
          "pho4_MVA  := pho4_MVA ",
          "pho1_pt   := pho1.pt()"
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
