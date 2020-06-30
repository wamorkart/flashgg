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


        "diphoPair_MVA := diphoPair_MVA",
        "a1_mass_bdt       := h4gDiPho1.mass()",
        "a2_mass_bdt       := h4gDiPho2.mass()",
        "a1_pt_bdt        := h4gDiPho1.pt()",
        "a2_pt_bdt        := h4gDiPho2.pt()",
        "a1_eta_bdt        := h4gDiPho1.eta()",
        "a2_eta_bdt        := h4gDiPho2.eta()",
        "a1_dR_bdt        := deltaR(h4gDiPho1_Pho1.eta(), h4gDiPho1_Pho1.phi(), h4gDiPho1_Pho2.eta(), h4gDiPho1_Pho2.phi())",
        "a2_dR_bdt        := deltaR(h4gDiPho2_Pho1.eta(), h4gDiPho2_Pho1.phi(), h4gDiPho2_Pho2.eta(), h4gDiPho2_Pho2.phi())",
        "a1_a2_dR_bdt     := deltaR(h4gDiPho1.eta(),h4gDiPho1.phi(),h4gDiPho2.eta(),h4gDiPho2.phi())",
        "cosThetaStarCS_bdt := cosThetaStarCS",
        "cosTheta_a1_bdt := cosTheta_a1",
        "cosTheta_a2_bdt := cosTheta_a2",
        "a1_p1i_bdt      := h4gDiPho1_iPho1",
        "a1_p2i_bdt      := h4gDiPho1_iPho2",
        "a2_p1i_bdt      := h4gDiPho2_iPho1",
        "a2_p2i_bdt      := h4gDiPho2_iPho2",
        "a1_mass_dM      := h4gDiPho1_prime.mass()",
        "a2_mass_dM        := h4gDiPho2_prime.mass()",
        "a1_pt_dM        := h4gDiPho1_prime.pt()",
        "a2_pt_dM        := h4gDiPho2_prime.pt()",
        "a1_eta_dM        := h4gDiPho1_prime.eta()",
        "a2_eta_dM        := h4gDiPho2_prime.eta()",
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
        "tp_pt         := tp.pt()",
        "tp_eta         := tp.eta()",
        "tp_mass         := tp.mass()"

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


        "diphoPair_MVA := diphoPair_MVA",
        "a1_mass_bdt       := h4gDiPho1.mass()",
        "a2_mass_bdt       := h4gDiPho2.mass()",
        "a1_pt_bdt        := h4gDiPho1.pt()",
        "a2_pt_bdt        := h4gDiPho2.pt()",
        "a1_eta_bdt        := h4gDiPho1.eta()",
        "a2_eta_bdt        := h4gDiPho2.eta()",
        "a1_dR_bdt        := deltaR(h4gDiPho1_Pho1.eta(), h4gDiPho1_Pho1.phi(), h4gDiPho1_Pho2.eta(), h4gDiPho1_Pho2.phi())",
        "a2_dR_bdt        := deltaR(h4gDiPho2_Pho1.eta(), h4gDiPho2_Pho1.phi(), h4gDiPho2_Pho2.eta(), h4gDiPho2_Pho2.phi())",
        "a1_a2_dR_bdt     := deltaR(h4gDiPho1.eta(),h4gDiPho1.phi(),h4gDiPho2.eta(),h4gDiPho2.phi())",
        "cosThetaStarCS_bdt := cosThetaStarCS",
        "cosTheta_a1_bdt := cosTheta_a1",
        "cosTheta_a2_bdt := cosTheta_a2",
        "a1_p1i_bdt      := h4gDiPho1_iPho1",
        "a1_p2i_bdt      := h4gDiPho1_iPho2",
        "a2_p1i_bdt      := h4gDiPho2_iPho1",
        "a2_p2i_bdt      := h4gDiPho2_iPho2",
        "a1_mass_dM      := h4gDiPho1_prime.mass()",
        "a2_mass_dM        := h4gDiPho2_prime.mass()",
        "a1_pt_dM        := h4gDiPho1_prime.pt()",
        "a2_pt_dM        := h4gDiPho2_prime.pt()",
        "a1_eta_dM        := h4gDiPho1_prime.eta()",
        "a2_eta_dM        := h4gDiPho2_prime.eta()",
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


        "diphoPair_MVA := diphoPair_MVA",
        "a1_mass_bdt       := h4gDiPho1.mass()",
        "a2_mass_bdt       := h4gDiPho2.mass()",
        "a1_pt_bdt        := h4gDiPho1.pt()",
        "a2_pt_bdt        := h4gDiPho2.pt()",
        "a1_eta_bdt        := h4gDiPho1.eta()",
        "a2_eta_bdt        := h4gDiPho2.eta()",
        "a1_dR_bdt        := deltaR(h4gDiPho1_Pho1.eta(), h4gDiPho1_Pho1.phi(), h4gDiPho1_Pho2.eta(), h4gDiPho1_Pho2.phi())",
        "a2_dR_bdt        := deltaR(h4gDiPho2_Pho1.eta(), h4gDiPho2_Pho1.phi(), h4gDiPho2_Pho2.eta(), h4gDiPho2_Pho2.phi())",
        "a1_a2_dR_bdt     := deltaR(h4gDiPho1.eta(),h4gDiPho1.phi(),h4gDiPho2.eta(),h4gDiPho2.phi())",
        "cosThetaStarCS_bdt := cosThetaStarCS",
        "cosTheta_a1_bdt := cosTheta_a1",
        "cosTheta_a2_bdt := cosTheta_a2",
        "a1_p1i_bdt      := h4gDiPho1_iPho1",
        "a1_p2i_bdt      := h4gDiPho1_iPho2",
        "a2_p1i_bdt      := h4gDiPho2_iPho1",
        "a2_p2i_bdt      := h4gDiPho2_iPho2",
        "a1_mass_dM      := h4gDiPho1_prime.mass()",
        "a2_mass_dM        := h4gDiPho2_prime.mass()",
        "a1_pt_dM        := h4gDiPho1_prime.pt()",
        "a2_pt_dM        := h4gDiPho2_prime.pt()",
        "a1_eta_dM        := h4gDiPho1_prime.eta()",
        "a2_eta_dM        := h4gDiPho2_prime.eta()",
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
        "tp_pt         := tp.pt()",
        "tp_eta         := tp.eta()",
        "tp_mass         := tp.mass()"
           ]

      return systematicVariables





    def customizeTagSequence(self):
        self.process.load("flashgg.Taggers.flashggH4GTag_cff")

        self.process.flashggH4GTag.vertexIdMVAweightfileH4G = cms.untracked.FileInPath(str(self.metaConditions["H4GTag"]["vertexIdMVAweightfileH4G"]))
        if (cms.string(str(self.metaConditions["H4GTag"]["year"])) == "2016"):
            from flashgg.Taggers.flashggPreselectedDiPhotons_LowMass16_cfi import flashggPreselectedDiPhotonsLowMass
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

        # self.process.flashggH4GTag.preselecteddiphotonsLM = cms.
        # teststr = self.metaConditions["H4GTag"]["preSelectedDiPhotons_cfi"]
        # from cms.InputTag(self.metaConditions["H4GTag"]["preSelectedDiPhotons_cfi"]) import flashggPreselectedDiPhotonsLowMass
        # from str(self.metaConditions["H4GTag"]["preSelectedDiPhotons_cfi"]) import flashggPreselectedDiPhotons
        # process.flashggH4GTag.idSelection = cms.PSet(
                    # rho = flashggPreselectedDiPhotonsLowMass.rho,
                    # cut = flashggPreselectedDiPhotonsLowMass.cut,
                    # variables = flashggPreselectedDiPhotonsLowMass.variables,
                    # categories = flashggPreselectedDiPhotonsLowMass.categories)
        # from flashgg.Taggers.flashggPreselectedDiPhotons_LowMass16_cfi import flashggPreselectedDiPhotonsLowMass
        # process.flashggH4GTag.idSelection = cms.PSet(
        #             rho = flashggPreselectedDiPhotonsLowMass.rho,
        #             cut = flashggPreselectedDiPhotonsLowMass.cut,
        #             variables = flashggPreselectedDiPhotonsLowMass.variables,
        #             categories = flashggPreselectedDiPhotonsLowMass.categories)

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
