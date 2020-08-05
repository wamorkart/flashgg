#ifndef flashgg_H4GTag
#define flashgg_H4GTag

// https://root.cern.ch/doc/v608/TLorentzVector_8h_source.html
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"

#include "flashgg/Taggers/interface/FunctionHelpers.h"
#include "TMVA/Reader.h"

namespace flashgg{

  class H4GTag: public DiPhotonTagBase, public reco::LeafCandidate
  {
  public:
    H4GTag();
    ~H4GTag();
     // H4GTag(edm::Ptr<DiPhotonCandidate>);
     H4GTag(edm::Ptr<DiPhotonCandidate>,  flashgg::Photon, flashgg::Photon, flashgg::Photon, flashgg::Photon, edm::Ptr<reco::Vertex>,  float, float);
     H4GTag(edm::Ptr<DiPhotonCandidate>, flashgg::Photon, flashgg::Photon, flashgg::Photon, edm::Ptr<reco::Vertex>, float, float);
     H4GTag(edm::Ptr<DiPhotonCandidate>, flashgg::Photon, flashgg::Photon, edm::Ptr<reco::Vertex>, float, float);


    virtual H4GTag *clone() const override;
    const flashgg::Photon pho1() const {return pho1_;};
    const flashgg::Photon pho2() const {return pho2_;};
    const flashgg::Photon pho3() const {return pho3_;};
    const flashgg::Photon pho4() const {return pho4_;};

    float pho1_MVA() const {return pho1_MVA_;};
    float pho2_MVA() const {return pho2_MVA_;};
    float pho3_MVA() const {return pho3_MVA_;};
    float pho4_MVA() const {return pho4_MVA_;};

    float pho1_full5x5_r9() const {return pho1_full5x5_r9_;};
    float pho2_full5x5_r9() const {return pho2_full5x5_r9_;};
    float pho3_full5x5_r9() const {return pho3_full5x5_r9_;};
    float pho4_full5x5_r9() const {return pho4_full5x5_r9_;};

    float pho1_egChargedHadronIso() const {return pho1_egChargedHadronIso_;};
    float pho2_egChargedHadronIso() const {return pho2_egChargedHadronIso_;};
    float pho3_egChargedHadronIso() const {return pho3_egChargedHadronIso_;};
    float pho4_egChargedHadronIso() const {return pho4_egChargedHadronIso_;};

    float pho1_hadronicOverEm() const {return pho1_hadronicOverEm_;};
    float pho2_hadronicOverEm() const {return pho2_hadronicOverEm_;};
    float pho3_hadronicOverEm() const {return pho3_hadronicOverEm_;};
    float pho4_hadronicOverEm() const {return pho4_hadronicOverEm_;};

    float pho1_SC_eta() const {return pho1_SC_eta_;};
    float pho2_SC_eta() const {return pho2_SC_eta_;};
    float pho3_SC_eta() const {return pho3_SC_eta_;};
    float pho4_SC_eta() const {return pho4_SC_eta_;};


    // float diphoPair_MVA() const {return diphoPair_MVA_;};

    const reco::Candidate::LorentzVector& pho12() const { return pho12_; };
    const reco::Candidate::LorentzVector& pho13() const { return pho13_; };
    const reco::Candidate::LorentzVector& pho14() const { return pho14_; };
    const reco::Candidate::LorentzVector& pho23() const { return pho23_; };
    const reco::Candidate::LorentzVector& pho24() const { return pho24_; };
    const reco::Candidate::LorentzVector& pho34() const { return pho34_; };

    // const reco::Candidate::LorentzVector& h4gDiPho1() const { return dp1_; };
    // const reco::Candidate::LorentzVector& h4gDiPho2() const { return dp2_; };
    // const reco::Candidate::LorentzVector& h4gDiPho1_Pho1() const { return dp1_pho1_; };
    // const reco::Candidate::LorentzVector& h4gDiPho1_Pho2() const { return dp1_pho2_; };
    // const reco::Candidate::LorentzVector& h4gDiPho2_Pho1() const { return dp2_pho1_; };
    // const reco::Candidate::LorentzVector& h4gDiPho2_Pho2() const { return dp2_pho2_; };
    // const int& h4gDiPho1_iPho1() const { return dp1_ipho1_; };
    // const int& h4gDiPho1_iPho2() const { return dp1_ipho2_; };
    // const int& h4gDiPho2_iPho1() const { return dp2_ipho1_; };
    // const int& h4gDiPho2_iPho2() const { return dp2_ipho2_; };
    // const float& cosThetaStarCS() const {return cosThetaStarCS_;};
    // const float& cosTheta_a1() const {return cosTheta_a1_;};
    // const float& cosTheta_a2() const {return cosTheta_a2_;};

    const reco::Candidate::LorentzVector& h4gDiPho1_prime() const { return dp1_prime_; };
    const reco::Candidate::LorentzVector& h4gDiPho2_prime() const { return dp2_prime_; };
    const reco::Candidate::LorentzVector& h4gDiPho1_Pho1_prime() const { return dp1_pho1_prime_; };
    const reco::Candidate::LorentzVector& h4gDiPho1_Pho2_prime() const { return dp1_pho2_prime_; };
    const reco::Candidate::LorentzVector& h4gDiPho2_Pho1_prime() const { return dp2_pho1_prime_; };
    const reco::Candidate::LorentzVector& h4gDiPho2_Pho2_prime() const { return dp2_pho2_prime_; };
    const int& h4gDiPho1_iPho1_prime() const { return dp1_ipho1_prime_; };
    const int& h4gDiPho1_iPho2_prime() const { return dp1_ipho2_prime_; };
    const int& h4gDiPho2_iPho1_prime() const { return dp2_ipho1_prime_; };
    const int& h4gDiPho2_iPho2_prime() const { return dp2_ipho2_prime_; };
    const float& cosThetaStarCS_prime() const {return cosThetaStarCS_prime_;};
    const float& cosTheta_a1_prime() const {return cosTheta_a1_prime_;};
    const float& cosTheta_a2_prime() const {return cosTheta_a2_prime_;};


    const reco::Candidate::LorentzVector& tp() const { return tp_; };

    float dZ_bdtVtx() const {return dZ_bdtVtx_;};
    float dZ_ZeroVtx() const {return dZ_ZeroVtx_;};
    float getCosThetaStar_CS(reco::Candidate::LorentzVector a1, reco::Candidate::LorentzVector a2) const;
    std::vector<float> CosThetaAngles(reco::Candidate::LorentzVector a1, reco::Candidate::LorentzVector a2, reco::Candidate::LorentzVector a1_pho1, reco::Candidate::LorentzVector a2_pho1) const;
    float HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const;






  private:
    flashgg::Photon pho1_;
    flashgg::Photon pho2_;
    flashgg::Photon pho3_;
    flashgg::Photon pho4_;
    float pho1_MVA_;
    float pho2_MVA_;
    float pho3_MVA_;
    float pho4_MVA_;
    float pho1_full5x5_r9_;
    float pho2_full5x5_r9_;
    float pho3_full5x5_r9_;
    float pho4_full5x5_r9_;
    float pho1_egChargedHadronIso_;
    float pho2_egChargedHadronIso_;
    float pho3_egChargedHadronIso_;
    float pho4_egChargedHadronIso_;
    float pho1_hadronicOverEm_;
    float pho2_hadronicOverEm_;
    float pho3_hadronicOverEm_;
    float pho4_hadronicOverEm_;
    float pho1_SC_eta_;
    float pho2_SC_eta_;
    float pho3_SC_eta_;
    float pho4_SC_eta_;
    float diphoPair_MVA_;
    reco::Candidate::LorentzVector pho12_;
    reco::Candidate::LorentzVector pho13_;
    reco::Candidate::LorentzVector pho14_;
    reco::Candidate::LorentzVector pho23_;
    reco::Candidate::LorentzVector pho24_;
    reco::Candidate::LorentzVector pho34_;
    // reco::Candidate::LorentzVector dp1_;
    // reco::Candidate::LorentzVector dp2_;
    // reco::Candidate::LorentzVector dp1_pho1_;
    // reco::Candidate::LorentzVector dp1_pho2_;
    // reco::Candidate::LorentzVector dp2_pho1_;
    // reco::Candidate::LorentzVector dp2_pho2_;
    // int dp1_ipho1_;
    // int dp1_ipho2_;
    // int dp2_ipho1_;
    // int dp2_ipho2_;
    // float cosThetaStarCS_;
    // float cosTheta_a1_;
    // float cosTheta_a2_;
    reco::Candidate::LorentzVector dp1_prime_;
    reco::Candidate::LorentzVector dp2_prime_;
    reco::Candidate::LorentzVector dp1_pho1_prime_;
    reco::Candidate::LorentzVector dp1_pho2_prime_;
    reco::Candidate::LorentzVector dp2_pho1_prime_;
    reco::Candidate::LorentzVector dp2_pho2_prime_;
    int dp1_ipho1_prime_;
    int dp1_ipho2_prime_;
    int dp2_ipho1_prime_;
    int dp2_ipho2_prime_;
    float cosThetaStarCS_prime_;
    float cosTheta_a1_prime_;
    float cosTheta_a2_prime_;
    reco::Candidate::LorentzVector tp_;
    float dZ_bdtVtx_;
    float dZ_ZeroVtx_;


  };
}

#endif
