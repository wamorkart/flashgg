import FWCore.ParameterSet.Config as cms

class H4GCustomize():


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
        "pho1_phi   := pho1.phi()",
        "pho2_phi   := pho2.phi()",
        "pho3_phi   := pho3.phi()",
        "pho4_phi   := pho4.phi()",
        "pho1_pfPhoIso03 := pho1_pfPhoIso03",
        "pho2_pfPhoIso03 := pho2_pfPhoIso03",
        "pho3_pfPhoIso03 := pho3_pfPhoIso03",
        "pho4_pfPhoIso03 := pho4_pfPhoIso03",
        "pho1_trkSumPtHollowConeDR03 := pho1_trkSumPtHollowConeDR03",
        "pho2_trkSumPtHollowConeDR03 := pho2_trkSumPtHollowConeDR03",
        "pho3_trkSumPtHollowConeDR03 := pho3_trkSumPtHollowConeDR03",
        "pho4_trkSumPtHollowConeDR03 := pho4_trkSumPtHollowConeDR03",
        "pho1_full5x5_sigmaIetaIeta := pho1_full5x5_sigmaIetaIeta",
        "pho2_full5x5_sigmaIetaIeta := pho2_full5x5_sigmaIetaIeta",
        "pho3_full5x5_sigmaIetaIeta := pho3_full5x5_sigmaIetaIeta",
        "pho4_full5x5_sigmaIetaIeta := pho4_full5x5_sigmaIetaIeta",

        "pho1_MVA  := pho1_MVA ",
        "pho2_MVA  := pho2_MVA ",
        "pho3_MVA  := pho3_MVA ",
        "pho4_MVA  := pho4_MVA ",
        "pho1_energy := pho1.energy()",
        "pho2_energy := pho2.energy()",
        "pho3_energy := pho3.energy()",
        "pho4_energy := pho4.energy()",
        "pho1_electronveto := pho1.passElectronVeto()",
        "pho2_electronveto := pho2.passElectronVeto()",
        "pho3_electronveto := pho3.passElectronVeto()",
        "pho4_electronveto := pho4.passElectronVeto()",
        "pho1_pixelseed := pho1.hasPixelSeed()",
        "pho2_pixelseed := pho2.hasPixelSeed()",
        "pho3_pixelseed := pho3.hasPixelSeed()",
        "pho4_pixelseed := pho4.hasPixelSeed()",
        "pho1_full5x5_r9 := pho1_full5x5_r9",
        "pho2_full5x5_r9 := pho2_full5x5_r9",
        "pho3_full5x5_r9 := pho3_full5x5_r9",
        "pho4_full5x5_r9 := pho4_full5x5_r9",
        "pho1_chHadIso := pho1_egChargedHadronIso",
        "pho2_chHadIso := pho2_egChargedHadronIso",
        "pho3_chHadIso := pho3_egChargedHadronIso",
        "pho4_chHadIso := pho4_egChargedHadronIso",
        "pho1_HoE := pho1_hadronicOverEm",
        "pho2_HoE := pho2_hadronicOverEm",
        "pho3_HoE := pho3_hadronicOverEm",
        "pho4_HoE := pho4_hadronicOverEm",
        "pho1_SC_Eta := pho1_SC_eta",
        "pho2_SC_Eta := pho2_SC_eta",
        "pho3_SC_Eta := pho3_SC_eta",
        "pho4_SC_Eta := pho4_SC_eta",
        "pho1_genType := pho1.genMatchType()",
        "pho1_genType := pho1.genMatchType()",
        "pho2_genType := pho2.genMatchType()",
        "pho3_genType := pho3.genMatchType()",
        "pho4_genType := pho4.genMatchType()",

        "pho1_energy := pho1.energy()",
        "pho1_energy_init := pho1.energyAtStep('initial')",

        "pho12_pt :=  pho12.pt()",
        "pho13_pt :=  pho13.pt()",
        "pho14_pt :=  pho14.pt()",
        "pho23_pt :=  pho23.pt()",
        "pho24_pt :=  pho24.pt()",
        "pho34_pt :=  pho34.pt()",
        "pho12_eta :=  pho12.eta()",
        "pho13_eta :=  pho13.eta()",
        "pho14_eta :=  pho14.eta()",
        "pho23_eta :=  pho23.eta()",
        "pho24_eta :=  pho24.eta()",
        "pho34_eta :=  pho34.eta()",
        "pho12_mass :=  pho12.mass()",
        "pho13_mass :=  pho13.mass()",
        "pho14_mass :=  pho14.mass()",
        "pho23_mass :=  pho23.mass()",
        "pho24_mass :=  pho24.mass()",
        "pho34_mass :=  pho34.mass()",
        "pho12_dR   := deltaR(pho1.eta(), pho1.phi(), pho2.eta(), pho2.phi())",
        "pho13_dR   := deltaR(pho1.eta(), pho1.phi(), pho3.eta(), pho3.phi())",
        "pho14_dR   := deltaR(pho1.eta(), pho1.phi(), pho4.eta(), pho4.phi())",
        "pho23_dR   := deltaR(pho2.eta(), pho2.phi(), pho3.eta(), pho3.phi())",
        "pho24_dR   := deltaR(pho2.eta(), pho2.phi(), pho4.eta(), pho4.phi())",
        "pho34_dR   := deltaR(pho3.eta(), pho3.phi(), pho4.eta(), pho4.phi())",

        "a1_mass_dM      := h4gDiPho1_prime.mass()",
        "a2_mass_dM        := h4gDiPho2_prime.mass()",
        "a1_pt_dM        := h4gDiPho1_prime.pt()",
        "a2_pt_dM        := h4gDiPho2_prime.pt()",
        "a1_eta_dM        := h4gDiPho1_prime.eta()",
        "a2_eta_dM        := h4gDiPho2_prime.eta()",
        "a1_energy_dM     := h4gDiPho1_prime.energy()",
        "a2_energy_dM     := h4gDiPho2_prime.energy()",
        "a1_dR_dM        := deltaR(h4gDiPho1_Pho1_prime.eta(), h4gDiPho1_Pho1_prime.phi(), h4gDiPho1_Pho2_prime.eta(), h4gDiPho1_Pho2_prime.phi())",
        "a2_dR_dM        := deltaR(h4gDiPho2_Pho1_prime.eta(), h4gDiPho2_Pho1_prime.phi(), h4gDiPho2_Pho2_prime.eta(), h4gDiPho2_Pho2_prime.phi())",
        "a1_a2_dR_dM     := deltaR(h4gDiPho1_prime.eta(),h4gDiPho1_prime.phi(),h4gDiPho2_prime.eta(),h4gDiPho2_prime.phi())",
        "cosThetaStarCS_dM := cosThetaStarCS_prime",
        "cosTheta_a1_dM := cosTheta_a1_prime",
        "cosTheta_a2_dM := cosTheta_a2_prime",
        "a1_p1i_dM      := h4gDiPho1_iPho1_prime",
        "a1_p2i_dM      := h4gDiPho1_iPho2_prime",
        "a2_p1i_dM      := h4gDiPho2_iPho1_prime",
        "a2_p2i_dM      := h4gDiPho2_iPho2_prime",


        "dZ_bdtVtx        := dZ_bdtVtx",
        "dZ_ZeroVtx        := dZ_ZeroVtx",
        "dZ_HggVtx        := dZ_HggVtx",
        "tp_pt         := tp.pt()",
        "tp_eta         := tp.eta()",
        "tp_mass         := tp.mass()",
        "gen_pho1_pt     := gen_pho1_pt",
        "gen_pho2_pt     := gen_pho2_pt",
        "gen_pho3_pt     := gen_pho3_pt",
        "gen_pho4_pt     := gen_pho4_pt",
        "gen_pho1_eta     := gen_pho1_eta",
        "gen_pho2_eta     := gen_pho2_eta",
        "gen_pho3_eta     := gen_pho3_eta",
        "gen_pho4_eta     := gen_pho4_eta",
        "gen_pho12_dR     := gen_pho12_dR",
        "gen_pho13_dR     := gen_pho13_dR",
        "gen_pho14_dR     := gen_pho14_dR",
        "gen_pho23_dR     := gen_pho23_dR",
        "gen_pho24_dR     := gen_pho24_dR",
        "gen_pho34_dR     := gen_pho34_dR",
        "gen_pho12_M      := gen_pho12_M",
        "gen_pho13_M      := gen_pho13_M",
        "gen_pho14_M      := gen_pho14_M",
        "gen_pho23_M      := gen_pho23_M",
        "gen_pho24_M      := gen_pho24_M",
        "gen_pho34_M      := gen_pho34_M",
        "gen_a1_pt        := gen_a1_pt",
        "gen_a2_pt        := gen_a2_pt",
        "gen_a1_eta      := gen_a1_eta",
        "gen_a2_eta      := gen_a2_eta",
        "gen_a1a2_dR      := gen_a1a2_dR",
        "gen_h_mass       := gen_h_mass",
        "gen_h_pt         := gen_h_pt",
        "gen_h_eta        := gen_h_eta"

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
        "pho1_phi   := pho1.phi()",
        "pho2_phi   := pho2.phi()",
        "pho3_phi   := pho3.phi()",
        "pho4_phi   := pho4.phi()",
        "pho1_pfPhoIso03 := pho1_pfPhoIso03",
        "pho2_pfPhoIso03 := pho2_pfPhoIso03",
        "pho3_pfPhoIso03 := pho3_pfPhoIso03",
        "pho4_pfPhoIso03 := pho4_pfPhoIso03",
        "pho1_trkSumPtHollowConeDR03 := pho1_trkSumPtHollowConeDR03",
        "pho2_trkSumPtHollowConeDR03 := pho2_trkSumPtHollowConeDR03",
        "pho3_trkSumPtHollowConeDR03 := pho3_trkSumPtHollowConeDR03",
        "pho4_trkSumPtHollowConeDR03 := pho4_trkSumPtHollowConeDR03",
        "pho1_full5x5_sigmaIetaIeta := pho1_full5x5_sigmaIetaIeta",
        "pho2_full5x5_sigmaIetaIeta := pho2_full5x5_sigmaIetaIeta",
        "pho3_full5x5_sigmaIetaIeta := pho3_full5x5_sigmaIetaIeta",
        "pho4_full5x5_sigmaIetaIeta := pho4_full5x5_sigmaIetaIeta",
        "pho1_MVA  := pho1_MVA ",
        "pho2_MVA  := pho2_MVA ",
        "pho3_MVA  := pho3_MVA ",
        "pho4_MVA  := pho4_MVA ",
        "pho1_energy := pho1.energy()",
        "pho2_energy := pho2.energy()",
        "pho3_energy := pho3.energy()",
        "pho4_energy := pho4.energy()",
        "pho1_electronveto := pho1.passElectronVeto()",
        "pho2_electronveto := pho2.passElectronVeto()",
        "pho3_electronveto := pho3.passElectronVeto()",
        "pho4_electronveto := pho4.passElectronVeto()",
        "pho1_pixelseed := pho1.hasPixelSeed()",
        "pho2_pixelseed := pho2.hasPixelSeed()",
        "pho3_pixelseed := pho3.hasPixelSeed()",
        "pho4_pixelseed := pho4.hasPixelSeed()",
        "pho1_full5x5_r9 := pho1_full5x5_r9",
        "pho2_full5x5_r9 := pho2_full5x5_r9",
        "pho3_full5x5_r9 := pho3_full5x5_r9",
        "pho4_full5x5_r9 := pho4_full5x5_r9",
        "pho1_chHadIso := pho1_egChargedHadronIso",
        "pho2_chHadIso := pho2_egChargedHadronIso",
        "pho3_chHadIso := pho3_egChargedHadronIso",
        "pho4_chHadIso := pho4_egChargedHadronIso",
        "pho1_HoE := pho1_hadronicOverEm",
        "pho2_HoE := pho2_hadronicOverEm",
        "pho3_HoE := pho3_hadronicOverEm",
        "pho4_HoE := pho4_hadronicOverEm",
        "pho1_SC_Eta := pho1_SC_eta",
        "pho2_SC_Eta := pho2_SC_eta",
        "pho3_SC_Eta := pho3_SC_eta",
        "pho4_SC_Eta := pho4_SC_eta",
        "pho1_genType := pho1.genMatchType()",
        "pho1_genType := pho1.genMatchType()",
        "pho2_genType := pho2.genMatchType()",
        "pho3_genType := pho3.genMatchType()",
        "pho4_genType := pho4.genMatchType()",

        "pho1_energy := pho1.energy()",
        "pho1_energy_init := pho1.energyAtStep('initial')",

        "pho12_pt :=  pho12.pt()",
        "pho13_pt :=  pho13.pt()",
        "pho14_pt :=  pho14.pt()",
        "pho23_pt :=  pho23.pt()",
        "pho24_pt :=  pho24.pt()",
        "pho34_pt :=  pho34.pt()",
        "pho12_eta :=  pho12.eta()",
        "pho13_eta :=  pho13.eta()",
        "pho14_eta :=  pho14.eta()",
        "pho23_eta :=  pho23.eta()",
        "pho24_eta :=  pho24.eta()",
        "pho34_eta :=  pho34.eta()",
        "pho12_mass :=  pho12.mass()",
        "pho13_mass :=  pho13.mass()",
        "pho14_mass :=  pho14.mass()",
        "pho23_mass :=  pho23.mass()",
        "pho24_mass :=  pho24.mass()",
        "pho34_mass :=  pho34.mass()",
        "pho12_dR   := deltaR(pho1.eta(), pho1.phi(), pho2.eta(), pho2.phi())",
        "pho13_dR   := deltaR(pho1.eta(), pho1.phi(), pho3.eta(), pho3.phi())",
        "pho14_dR   := deltaR(pho1.eta(), pho1.phi(), pho4.eta(), pho4.phi())",
        "pho23_dR   := deltaR(pho2.eta(), pho2.phi(), pho3.eta(), pho3.phi())",
        "pho24_dR   := deltaR(pho2.eta(), pho2.phi(), pho4.eta(), pho4.phi())",
        "pho34_dR   := deltaR(pho3.eta(), pho3.phi(), pho4.eta(), pho4.phi())",

        "a1_mass_dM      := h4gDiPho1_prime.mass()",
        "a2_mass_dM        := h4gDiPho2_prime.mass()",
        "a1_pt_dM        := h4gDiPho1_prime.pt()",
        "a2_pt_dM        := h4gDiPho2_prime.pt()",
        "a1_eta_dM        := h4gDiPho1_prime.eta()",
        "a2_eta_dM        := h4gDiPho2_prime.eta()",
        "a1_energy_dM     := h4gDiPho1_prime.energy()",
        "a2_energy_dM     := h4gDiPho2_prime.energy()",
        "a1_dR_dM        := deltaR(h4gDiPho1_Pho1_prime.eta(), h4gDiPho1_Pho1_prime.phi(), h4gDiPho1_Pho2_prime.eta(), h4gDiPho1_Pho2_prime.phi())",
        "a2_dR_dM        := deltaR(h4gDiPho2_Pho1_prime.eta(), h4gDiPho2_Pho1_prime.phi(), h4gDiPho2_Pho2_prime.eta(), h4gDiPho2_Pho2_prime.phi())",
        "a1_a2_dR_dM     := deltaR(h4gDiPho1_prime.eta(),h4gDiPho1_prime.phi(),h4gDiPho2_prime.eta(),h4gDiPho2_prime.phi())",
        "cosThetaStarCS_dM := cosThetaStarCS_prime",
        "cosTheta_a1_dM := cosTheta_a1_prime",
        "cosTheta_a2_dM := cosTheta_a2_prime",
        "a1_p1i_dM      := h4gDiPho1_iPho1_prime",
        "a1_p2i_dM      := h4gDiPho1_iPho2_prime",
        "a2_p1i_dM      := h4gDiPho2_iPho1_prime",
        "a2_p2i_dM      := h4gDiPho2_iPho2_prime",


        "dZ_bdtVtx        := dZ_bdtVtx",
        "dZ_ZeroVtx        := dZ_ZeroVtx",
        "tp_pt         := tp.pt()",
        "tp_eta         := tp.eta()",
        "tp_mass         := tp.mass()"
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
        "pho1_phi   := pho1.phi()",
        "pho2_phi   := pho2.phi()",
        "pho3_phi   := pho3.phi()",
        "pho4_phi   := pho4.phi()",
        "pho1_pfPhoIso03 := pho1_pfPhoIso03",
        "pho2_pfPhoIso03 := pho2_pfPhoIso03",
        "pho3_pfPhoIso03 := pho3_pfPhoIso03",
        "pho4_pfPhoIso03 := pho4_pfPhoIso03",
        "pho1_trkSumPtHollowConeDR03 := pho1_trkSumPtHollowConeDR03",
        "pho2_trkSumPtHollowConeDR03 := pho2_trkSumPtHollowConeDR03",
        "pho3_trkSumPtHollowConeDR03 := pho3_trkSumPtHollowConeDR03",
        "pho4_trkSumPtHollowConeDR03 := pho4_trkSumPtHollowConeDR03",
        "pho1_full5x5_sigmaIetaIeta := pho1_full5x5_sigmaIetaIeta",
        "pho2_full5x5_sigmaIetaIeta := pho2_full5x5_sigmaIetaIeta",
        "pho3_full5x5_sigmaIetaIeta := pho3_full5x5_sigmaIetaIeta",
        "pho4_full5x5_sigmaIetaIeta := pho4_full5x5_sigmaIetaIeta",
        "pho1_MVA  := pho1_MVA ",
        "pho2_MVA  := pho2_MVA ",
        "pho3_MVA  := pho3_MVA ",
        "pho4_MVA  := pho4_MVA ",
        "pho1_energy := pho1.energy()",
        "pho2_energy := pho2.energy()",
        "pho3_energy := pho3.energy()",
        "pho4_energy := pho4.energy()",
        "pho1_electronveto := pho1.passElectronVeto()",
        "pho2_electronveto := pho2.passElectronVeto()",
        "pho3_electronveto := pho3.passElectronVeto()",
        "pho4_electronveto := pho4.passElectronVeto()",
        "pho1_pixelseed := pho1.hasPixelSeed()",
        "pho2_pixelseed := pho2.hasPixelSeed()",
        "pho3_pixelseed := pho3.hasPixelSeed()",
        "pho4_pixelseed := pho4.hasPixelSeed()",
        "pho1_full5x5_r9 := pho1_full5x5_r9",
        "pho2_full5x5_r9 := pho2_full5x5_r9",
        "pho3_full5x5_r9 := pho3_full5x5_r9",
        "pho4_full5x5_r9 := pho4_full5x5_r9",
        "pho1_chHadIso := pho1_egChargedHadronIso",
        "pho2_chHadIso := pho2_egChargedHadronIso",
        "pho3_chHadIso := pho3_egChargedHadronIso",
        "pho4_chHadIso := pho4_egChargedHadronIso",
        "pho1_HoE := pho1_hadronicOverEm",
        "pho2_HoE := pho2_hadronicOverEm",
        "pho3_HoE := pho3_hadronicOverEm",
        "pho4_HoE := pho4_hadronicOverEm",
        "pho1_SC_Eta := pho1_SC_eta",
        "pho2_SC_Eta := pho2_SC_eta",
        "pho3_SC_Eta := pho3_SC_eta",
        "pho4_SC_Eta := pho4_SC_eta",
        "pho1_genType := pho1.genMatchType()",
        "pho1_genType := pho1.genMatchType()",
        "pho2_genType := pho2.genMatchType()",
        "pho3_genType := pho3.genMatchType()",
        "pho4_genType := pho4.genMatchType()",

        "pho1_energy := pho1.energy()",
        "pho1_energy_init := pho1.energyAtStep('initial')",

        "pho12_pt :=  pho12.pt()",
        "pho13_pt :=  pho13.pt()",
        "pho14_pt :=  pho14.pt()",
        "pho23_pt :=  pho23.pt()",
        "pho24_pt :=  pho24.pt()",
        "pho34_pt :=  pho34.pt()",
        "pho12_eta :=  pho12.eta()",
        "pho13_eta :=  pho13.eta()",
        "pho14_eta :=  pho14.eta()",
        "pho23_eta :=  pho23.eta()",
        "pho24_eta :=  pho24.eta()",
        "pho34_eta :=  pho34.eta()",
        "pho12_mass :=  pho12.mass()",
        "pho13_mass :=  pho13.mass()",
        "pho14_mass :=  pho14.mass()",
        "pho23_mass :=  pho23.mass()",
        "pho24_mass :=  pho24.mass()",
        "pho34_mass :=  pho34.mass()",
        "pho12_dR   := deltaR(pho1.eta(), pho1.phi(), pho2.eta(), pho2.phi())",
        "pho13_dR   := deltaR(pho1.eta(), pho1.phi(), pho3.eta(), pho3.phi())",
        "pho14_dR   := deltaR(pho1.eta(), pho1.phi(), pho4.eta(), pho4.phi())",
        "pho23_dR   := deltaR(pho2.eta(), pho2.phi(), pho3.eta(), pho3.phi())",
        "pho24_dR   := deltaR(pho2.eta(), pho2.phi(), pho4.eta(), pho4.phi())",
        "pho34_dR   := deltaR(pho3.eta(), pho3.phi(), pho4.eta(), pho4.phi())",



        "a1_mass_dM      := h4gDiPho1_prime.mass()",
        "a2_mass_dM        := h4gDiPho2_prime.mass()",
        "a1_pt_dM        := h4gDiPho1_prime.pt()",
        "a2_pt_dM        := h4gDiPho2_prime.pt()",
        "a1_eta_dM        := h4gDiPho1_prime.eta()",
        "a2_eta_dM        := h4gDiPho2_prime.eta()",
        "a1_energy_dM     := h4gDiPho1_prime.energy()",
        "a2_energy_dM     := h4gDiPho2_prime.energy()",
        "a1_dR_dM        := deltaR(h4gDiPho1_Pho1_prime.eta(), h4gDiPho1_Pho1_prime.phi(), h4gDiPho1_Pho2_prime.eta(), h4gDiPho1_Pho2_prime.phi())",
        "a2_dR_dM        := deltaR(h4gDiPho2_Pho1_prime.eta(), h4gDiPho2_Pho1_prime.phi(), h4gDiPho2_Pho2_prime.eta(), h4gDiPho2_Pho2_prime.phi())",
        "a1_a2_dR_dM     := deltaR(h4gDiPho1_prime.eta(),h4gDiPho1_prime.phi(),h4gDiPho2_prime.eta(),h4gDiPho2_prime.phi())",
        "cosThetaStarCS_dM := cosThetaStarCS_prime",
        "cosTheta_a1_dM := cosTheta_a1_prime",
        "cosTheta_a2_dM := cosTheta_a2_prime",
        "a1_p1i_dM      := h4gDiPho1_iPho1_prime",
        "a1_p2i_dM      := h4gDiPho1_iPho2_prime",
        "a2_p1i_dM      := h4gDiPho2_iPho1_prime",
        "a2_p2i_dM      := h4gDiPho2_iPho2_prime",


        "dZ_bdtVtx        := dZ_bdtVtx",
        "dZ_ZeroVtx        := dZ_ZeroVtx",
        "tp_pt         := tp.pt()",
        "tp_eta         := tp.eta()",
        "tp_mass         := tp.mass()",
        "gen_pho1_pt     := gen_pho1_pt",
        "gen_pho2_pt     := gen_pho2_pt",
        "gen_pho3_pt     := gen_pho3_pt",
        "gen_pho4_pt     := gen_pho4_pt",
        "gen_pho1_eta     := gen_pho1_eta",
        "gen_pho2_eta     := gen_pho2_eta",
        "gen_pho3_eta     := gen_pho3_eta",
        "gen_pho4_eta     := gen_pho4_eta",
        "gen_pho12_dR     := gen_pho12_dR",
        "gen_pho13_dR     := gen_pho13_dR",
        "gen_pho14_dR     := gen_pho14_dR",
        "gen_pho23_dR     := gen_pho23_dR",
        "gen_pho24_dR     := gen_pho24_dR",
        "gen_pho34_dR     := gen_pho34_dR",
        "gen_pho12_M      := gen_pho12_M",
        "gen_pho13_M      := gen_pho13_M",
        "gen_pho14_M      := gen_pho14_M",
        "gen_pho23_M      := gen_pho23_M",
        "gen_pho24_M      := gen_pho24_M",
        "gen_pho34_M      := gen_pho34_M",
        "gen_a1_pt        := gen_a1_pt",
        "gen_a2_pt        := gen_a2_pt",
        "gen_a1_eta      := gen_a1_eta",
        "gen_a2_eta      := gen_a2_eta",
        "gen_a1a2_dR      := gen_a1a2_dR",
        "gen_h_mass       := gen_h_mass",
        "gen_h_pt         := gen_h_pt",
        "gen_h_eta        := gen_h_eta"


           ]

      return systematicVariables





    def customizeTagSequence(self):
        self.process.load("flashgg.Taggers.flashggH4GTag_cff")

        self.process.flashggH4GTag.vertexIdMVAweightfileH4G = cms.untracked.FileInPath(str(self.metaConditions["H4GTag"]["vertexIdMVAweightfileH4G"]))
        if (cms.string(str(self.metaConditions["H4GTag"]["year"])) == "2016"):
            from flashgg.Taggers.flashggPreselectedDiPhotons_LowMass16_cfi import flashggPreselectedDiPhotonsLowMass
            #from flashgg.Taggers.flashggPreselectedDiPhotons_LowMass16_NoCuts_cfi import flashggPreselectedDiPhotonsLowMass
            self.process.flashggH4GTag.idSelection = cms.PSet(
                         rho = flashggPreselectedDiPhotonsLowMass.rho,
                         cut = flashggPreselectedDiPhotonsLowMass.cut,
                         variables = flashggPreselectedDiPhotonsLowMass.variables,
                         categories = flashggPreselectedDiPhotonsLowMass.categories)
        elif (cms.string(str(self.metaConditions["H4GTag"]["year"])) == "2017"):
            from flashgg.Taggers.flashggPreselectedDiPhotons_LowMass17_cfi import flashggPreselectedDiPhotonsLowMass
            self.process.flashggH4GTag.idSelection = cms.PSet(
                         rho = flashggPreselectedDiPhotonsLowMass.rho,
                         cut = flashggPreselectedDiPhotonsLowMass.cut,
                         variables = flashggPreselectedDiPhotonsLowMass.variables,
                         categories = flashggPreselectedDiPhotonsLowMass.categories)

        elif (cms.string(str(self.metaConditions["H4GTag"]["year"])) == "2018"):
            from flashgg.Taggers.flashggPreselectedDiPhotons_LowMass18_cfi import flashggPreselectedDiPhotonsLowMass
            self.process.flashggH4GTag.idSelection = cms.PSet(
                         rho = flashggPreselectedDiPhotonsLowMass.rho,
                         cut = flashggPreselectedDiPhotonsLowMass.cut,
                         variables = flashggPreselectedDiPhotonsLowMass.variables,
                         categories = flashggPreselectedDiPhotonsLowMass.categories)


        print'Removing single Higgs tags'

        if self.customize.H4GTagsOnly:
            # self.process.flashggTagSequence.remove(self.process.flashggPrefireDiPhotons)
            #self.process.flashggTagSequence.remove(self.process.flashggPreselectedDiPhotons)
            #self.process.flashggTagSequence.replace(self.process.flashggPreselectedDiPhotons,self.process.flashggDiPhotonSystematics)
            self.process.flashggTagSequence.remove(self.process.flashggDiPhotonMVA)
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
            self.process.flashggTagSequence.remove(self.process.flashggVHhadMVA)
            self.process.flashggTagSequence.remove(self.process.flashggVBFDiPhoDiJetMVA)
            self.process.flashggTagSequence.remove(self.process.flashggTTHDiLeptonTag)
            self.process.flashggTagSequence.remove(self.process.flashggUntagged)
            self.process.flashggTagSequence.remove(self.process.flashggUntagged)
            self.process.flashggTagSequence.remove(self.process.flashggTHQLeptonicTag)

        #self.process.flashggTagSequence.replace(self.process.flashggTagSorter,self.process.flashggDiPhotonSystematics)
        # self.process.flashggTagSequence.replace(self.process.flashggTagSorter,self.process.flashggH4GTagSequence*self.process.flashggTagSorter)
        # self.process.flashggTagSorter.TagPriorityRanges = cms.VPSet( cms.PSet(TagName = cms.InputTag('flashggH4GTag')) )

    def H4GTagMerger(self,systlabels=[]):
        self.process.p.remove(self.process.flashggTagSorter)
        self.process.p.replace(self.process.flashggSystTagMerger,self.process.flashggH4GTagSequence*self.process.flashggTagSorter*self.process.flashggSystTagMerger)
        for systlabel in systlabels:
           if systlabel!='':
             self.process.p.remove(getattr(self.process,'flashggTagSorter'+systlabel))

             self.process.p.replace(self.process.flashggSystTagMerger,getattr(self.process, 'flashggTagSorter'+systlabel)*self.process.flashggSystTagMerger)
           setattr(getattr(self.process, 'flashggTagSorter'+systlabel), 'TagPriorityRanges', cms.VPSet( cms.PSet(TagName = cms.InputTag('flashggH4GTag', systlabel)) ))
        # print 'from loop after:',process.flashggSystTagMerger.src



    def H4GTagRunSequence(self,systlabels,phosystlabels):
        # print'not used'
        if self.customize.H4GTagsOnly:
          self.H4GTagMerger(systlabels)

        if len(systlabels)>1 :
          getattr(self.process, "flashggH4GTag").DiPhotonSuffixes = cms.vstring([systlabels[0]]+phosystlabels)
