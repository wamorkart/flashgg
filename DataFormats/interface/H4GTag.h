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
     H4GTag(edm::Ptr<DiPhotonCandidate> ,vector< flashgg::Photon>, vector<edm::Ptr<reco::Vertex>> ,int, vector<int>   );

    // H4GTag(edm::Ptr<DiPhotonCandidate>,vector<edm::Ptr<reco::Vertex>>,edm::Ptr<reco::Vertex>,reco::GenParticle::Point,int,TMVA::Reader *,vector< flashgg::Photon>,vector<int>, float );
    virtual H4GTag *clone() const override;

    const vector<flashgg::Photon> phoP4Corrected_dp() const { return phoP4Corrected_dp_;};
    const float pho1_MVA() const { return pho1_MVA_; };
    const float pho2_MVA() const { return pho2_MVA_; };
    const float pho3_MVA() const { return pho3_MVA_; };
    const float pho4_MVA() const { return pho4_MVA_; };
    const reco::Candidate::LorentzVector& pho12() const { return pho12_; };
    const reco::Candidate::LorentzVector& pho13() const { return pho13_; };
    const reco::Candidate::LorentzVector& pho14() const { return pho14_; };
    const reco::Candidate::LorentzVector& pho23() const { return pho23_; };
    const reco::Candidate::LorentzVector& pho24() const { return pho24_; };
    const reco::Candidate::LorentzVector& pho34() const { return pho34_; };
    const reco::Candidate::LorentzVector& h4gDiPho1() const { return dp1_; };
    const reco::Candidate::LorentzVector& h4gDiPho2() const { return dp2_; };
    const reco::Candidate::LorentzVector& h4gDiPho1_Pho1() const { return dp1_pho1_; };
    const reco::Candidate::LorentzVector& h4gDiPho1_Pho2() const { return dp1_pho2_; };
    const reco::Candidate::LorentzVector& h4gDiPho2_Pho1() const { return dp2_pho1_; };
    const reco::Candidate::LorentzVector& h4gDiPho2_Pho2() const { return dp2_pho2_; };
    const int& h4gDiPho1_iPho1() const { return dp1_ipho1_; };
    const int& h4gDiPho1_iPho2() const { return dp1_ipho2_; };
    const int& h4gDiPho2_iPho1() const { return dp2_ipho1_; };
    const int& h4gDiPho2_iPho2() const { return dp2_ipho2_; };
    const reco::Candidate::LorentzVector& tp() const { return tp_; };
    float getCosThetaStar_CS() const;
    std::vector<float> CosThetaAngles() const;
    float HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const;
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
    float getCosThetaStar_CS_prime() const;
    std::vector<float> CosThetaAngles_prime() const;


  private:
  //   float BS_factor_BDTVtx_;
    vector <flashgg::Photon> phoP4Corrected_dp_;
    float pho1_MVA_;
    float pho2_MVA_;
    float pho3_MVA_;
    float pho4_MVA_;
    reco::Candidate::LorentzVector  pho12_;
    reco::Candidate::LorentzVector  pho13_;
    reco::Candidate::LorentzVector  pho14_;
    reco::Candidate::LorentzVector  pho23_;
    reco::Candidate::LorentzVector  pho24_;
    reco::Candidate::LorentzVector  pho34_;
    reco::Candidate::LorentzVector dp1_;
    reco::Candidate::LorentzVector dp2_;
    reco::Candidate::LorentzVector dp1_pho1_;
    reco::Candidate::LorentzVector dp1_pho2_;
    reco::Candidate::LorentzVector dp2_pho1_;
    reco::Candidate::LorentzVector dp2_pho2_;
    int dp1_ipho1_;
    int dp1_ipho2_;
    int dp2_ipho1_;
    int dp2_ipho2_;
    reco::Candidate::LorentzVector tp_;
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

  };
}

#endif
