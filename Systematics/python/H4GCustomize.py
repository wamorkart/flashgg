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
        ]
        variables = [
            "pho1_pt   := pho1.pt()",
            "pho2_pt   := pho2.pt()",
            "pho3_pt   := pho3.pt()",
            "pho4_pt   := pho4.pt()",
            "pho1_eta   := pho1.eta()",
            "pho2_eta   := pho2.eta()",
            "pho3_eta   := pho3.eta()",
            "pho4_eta   := pho4.eta()",
            "pho1_MVA  := pho1_MVA ",
            "pho2_MVA  := pho2_MVA ",
            "pho3_MVA  := pho3_MVA ",
            "pho4_MVA  := pho4_MVA ",
            "pho1_electronveto := pho1.passElectronVeto()",
            "pho2_electronveto := pho2.passElectronVeto()",
            "pho3_electronveto := pho3.passElectronVeto()",
            "pho4_electronveto := pho4.passElectronVeto()",
            "pho1_pixelseed := pho1.hasPixelSeed()",
            "pho2_pixelseed := pho2.hasPixelSeed()",
            "pho3_pixelseed := pho3.hasPixelSeed()",
            "pho4_pixelseed := pho4.hasPixelSeed()",
            "pho1_genType := pho1.genMatchType()",
            "pho2_genType := pho2.genMatchType()",
            "pho3_genType := pho3.genMatchType()",
            "pho4_genType := pho4.genMatchType()",

            "pho1_energy := pho1.energy()",
            "pho1_energy_init := pho1.energyAtStep('initial')",
            "dZ_bdtVtx        := dZ_bdtVtx",
            "tp_pt         := tp.pt()",
            "tp_eta         := tp.eta()",
            "tp_mass         := tp.mass()",
            "a1_mass        := h4gDiPho1_prime.mass()",
            "a2_mass        := h4gDiPho2_prime.mass()",
            "a1_pt        := h4gDiPho1_prime.pt()",
            "a2_pt        := h4gDiPho2_prime.pt()",
            "a1_eta        := h4gDiPho1_prime.eta()",
            "a2_eta        := h4gDiPho2_prime.eta()",
            "a1_dR        := deltaR(h4gDiPho1_Pho1_prime.eta(), h4gDiPho1_Pho1_prime.phi(), h4gDiPho1_Pho2_prime.eta(), h4gDiPho1_Pho2_prime.phi())",
            "a2_dR        := deltaR(h4gDiPho2_Pho1_prime.eta(), h4gDiPho2_Pho1_prime.phi(), h4gDiPho2_Pho2_prime.eta(), h4gDiPho2_Pho2_prime.phi())",
            "a1_a2_dR     := deltaR(h4gDiPho1_prime.eta(),h4gDiPho1_prime.phi(),h4gDiPho2_prime.eta(),h4gDiPho2_prime.phi())",
            "CosThetaStar_CS      := getCosThetaStar_CS()",
            "CosTheta_pho_a1      := CosThetaAngles()[0]",
            "CosTheta_pho_a2      := CosThetaAngles()[1]"

        ]
        return variables


    def dataVariables(self):
        dataVariables=[
        "pho1_pt   := pho1.pt()",
        "pho2_pt   := pho2.pt()",
        "pho3_pt   := pho3.pt()",
        "pho4_pt   := pho4.pt()",
        "pho1_eta   := pho1.eta()",
        "pho2_eta   := pho2.eta()",
        "pho3_eta   := pho3.eta()",
        "pho4_eta   := pho4.eta()",
        "pho1_MVA  := pho1_MVA ",
        "pho2_MVA  := pho2_MVA ",
        "pho3_MVA  := pho3_MVA ",
        "pho4_MVA  := pho4_MVA ",
        "pho1_electronveto := pho1.passElectronVeto()",
        "pho2_electronveto := pho2.passElectronVeto()",
        "pho3_electronveto := pho3.passElectronVeto()",
        "pho4_electronveto := pho4.passElectronVeto()",
        "pho1_pixelseed := pho1.hasPixelSeed()",
        "pho2_pixelseed := pho2.hasPixelSeed()",
        "pho3_pixelseed := pho3.hasPixelSeed()",
        "pho4_pixelseed := pho4.hasPixelSeed()",
        "pho1_genType := pho1.genMatchType()",
        "pho2_genType := pho2.genMatchType()",
        "pho3_genType := pho3.genMatchType()",
        "pho4_genType := pho4.genMatchType()",

        "pho1_energy := pho1.energy()",
        "pho1_energy_init := pho1.energyAtStep('initial')",
        "dZ_bdtVtx        := dZ_bdtVtx",
        "tp_pt         := tp.pt()",
        "tp_eta         := tp.eta()",
        "tp_mass         := tp.mass()",
        "a1_mass        := h4gDiPho1_prime.mass()",
        "a2_mass        := h4gDiPho2_prime.mass()",
        "a1_pt        := h4gDiPho1_prime.pt()",
        "a2_pt        := h4gDiPho2_prime.pt()",
        "a1_eta        := h4gDiPho1_prime.eta()",
        "a2_eta        := h4gDiPho2_prime.eta()",
        "a1_dR        := deltaR(h4gDiPho1_Pho1_prime.eta(), h4gDiPho1_Pho1_prime.phi(), h4gDiPho1_Pho2_prime.eta(), h4gDiPho1_Pho2_prime.phi())",
        "a2_dR        := deltaR(h4gDiPho2_Pho1_prime.eta(), h4gDiPho2_Pho1_prime.phi(), h4gDiPho2_Pho2_prime.eta(), h4gDiPho2_Pho2_prime.phi())",
        "a1_a2_dR     := deltaR(h4gDiPho1_prime.eta(),h4gDiPho1_prime.phi(),h4gDiPho2_prime.eta(),h4gDiPho2_prime.phi())",
        "CosThetaStar_CS      := getCosThetaStar_CS()",
        "CosTheta_pho_a1      := CosThetaAngles()[0]",
        "CosTheta_pho_a2      := CosThetaAngles()[1]"
        ]
        return dataVariables
    def systematicVariables(self):
    #   systematicVariables=["CMS_hgg_mass[160,100,180]:=diPhoton().mass","Mjj[120,70,190]:=dijet().M()","HHbbggMVA[100,0,1.]:=MVA()","MX[300,250,5000]:=MX()"]
      systematicVariables=[
          "pho1_pt   := pho1.pt()",
          "pho2_pt   := pho2.pt()",
          "pho3_pt   := pho3.pt()",
          "pho4_pt   := pho4.pt()",
          "pho1_eta   := pho1.eta()",
          "pho2_eta   := pho2.eta()",
          "pho3_eta   := pho3.eta()",
          "pho4_eta   := pho4.eta()",
          "pho1_MVA  := pho1_MVA ",
          "pho2_MVA  := pho2_MVA ",
          "pho3_MVA  := pho3_MVA ",
          "pho4_MVA  := pho4_MVA ",
          "pho1_electronveto := pho1.passElectronVeto()",
          "pho2_electronveto := pho2.passElectronVeto()",
          "pho3_electronveto := pho3.passElectronVeto()",
          "pho4_electronveto := pho4.passElectronVeto()",
          "pho1_pixelseed := pho1.hasPixelSeed()",
          "pho2_pixelseed := pho2.hasPixelSeed()",
          "pho3_pixelseed := pho3.hasPixelSeed()",
          "pho4_pixelseed := pho4.hasPixelSeed()",
          "pho1_genType := pho1.genMatchType()",
          "pho2_genType := pho2.genMatchType()",
          "pho3_genType := pho3.genMatchType()",
          "pho4_genType := pho4.genMatchType()",

          "pho1_energy := pho1.energy()",
          "pho1_energy_init := pho1.energyAtStep('initial')",
          "dZ_bdtVtx        := dZ_bdtVtx",
          "tp_pt         := tp.pt()",
          "tp_eta         := tp.eta()",
          "tp_mass         := tp.mass()",
          "a1_mass        := h4gDiPho1_prime.mass()",
          "a2_mass        := h4gDiPho2_prime.mass()",
          "a1_pt        := h4gDiPho1_prime.pt()",
          "a2_pt        := h4gDiPho2_prime.pt()",
          "a1_eta        := h4gDiPho1_prime.eta()",
          "a2_eta        := h4gDiPho2_prime.eta()",
          "a1_dR        := deltaR(h4gDiPho1_Pho1_prime.eta(), h4gDiPho1_Pho1_prime.phi(), h4gDiPho1_Pho2_prime.eta(), h4gDiPho1_Pho2_prime.phi())",
          "a2_dR        := deltaR(h4gDiPho2_Pho1_prime.eta(), h4gDiPho2_Pho1_prime.phi(), h4gDiPho2_Pho2_prime.eta(), h4gDiPho2_Pho2_prime.phi())",
          "a1_a2_dR     := deltaR(h4gDiPho1_prime.eta(),h4gDiPho1_prime.phi(),h4gDiPho2_prime.eta(),h4gDiPho2_prime.phi())",
          "CosThetaStar_CS      := getCosThetaStar_CS()",
          "CosTheta_pho_a1      := CosThetaAngles()[0]",
          "CosTheta_pho_a2      := CosThetaAngles()[1]"
           ]

      return systematicVariables





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
