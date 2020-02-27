// include path /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_5_0/src/
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/H4GTag.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TMVA/Reader.h"

using namespace flashgg; // makes flashgg sub members visible
H4GTag::H4GTag() : DiPhotonTagBase::DiPhotonTagBase()
// BS_factor_BDTVtx_ (),
// pho1_MVA_ (),
// pho2_MVA_ (),
// pho3_MVA_ (),
// pho4_MVA_ ()
{
}

H4GTag::~H4GTag() {}

// H4GTag::H4GTag(edm::Ptr<DiPhotonCandidate> dipho )
H4GTag::H4GTag(edm::Ptr<DiPhotonCandidate> dipho,vector< flashgg::Photon>phoP4Corrected_dp, vector<edm::Ptr<reco::Vertex>> Vertices,int selected_vertex_index_, vector<int> diphoton_pairing_indices_  )
// H4GTag::H4GTag(edm::Ptr<DiPhotonCandidate> dipho,vector<edm::Ptr<reco::Vertex>> Vertices,edm::Ptr<reco::Vertex>vertex_bdt,reco::GenParticle::Point genVertex, int selected_vertex_index_, TMVA::Reader *VertexProbMva_,vector< flashgg::Photon>phoP4Corrected_dp,vector<int> diphoton_pairing_indices_, float diphoPair_MVA   )
{
  dipho_ = dipho;

  phoP4Corrected_dp_  = phoP4Corrected_dp;
  // cout << "dipho pt " << dipho->pt() << endl;

  // float vtx_X = Vertices[vertex_bdt]->x();
  // float vtx_Y = Vertices[vertex_bdt]->y();
  // float vtx_Z = Vertices[vertex_bdt]->z();

  //--Beam spot reweighting (https://github.com/cms-analysis/flashggFinalFit/blob/e60d53e19ac4f20e7ce187f0a34e483b4fc2a60e/Signal/test/SignalFit.cpp)
  // float mcBeamSpotWidth_=5.14; //cm
  // float dataBeamSpotWidth_=3.5; //cm
  //
  // float dZ_BDTVtx = genVertex.z() - Vertices[selected_vertex_index_]->z();;
  //
  // if (fabs(dZ_BDTVtx) < 0.1 ){
  //   BS_factor_BDTVtx_ =1;
  // } else {
  //   double mcBeamSpot_BDTVtx=TMath::Gaus(dZ_BDTVtx,0,TMath::Sqrt(2)*mcBeamSpotWidth_,true);
  //   double dataBeamSpot_BDTVtx=TMath::Gaus(dZ_BDTVtx,0,TMath::Sqrt(2)*dataBeamSpotWidth_,true);
  //   BS_factor_BDTVtx_ = dataBeamSpot_BDTVtx/mcBeamSpot_BDTVtx;
  // }
  pho1_MVA_ = phoP4Corrected_dp_.size() > 0 ? phoP4Corrected_dp_[0].phoIdMvaDWrtVtx(Vertices[selected_vertex_index_]) : -999;
  pho2_MVA_ = phoP4Corrected_dp_.size() > 0 ? phoP4Corrected_dp_[1].phoIdMvaDWrtVtx(Vertices[selected_vertex_index_]) : -999;
  pho3_MVA_ = phoP4Corrected_dp_.size() > 2 ? phoP4Corrected_dp_[2].phoIdMvaDWrtVtx(Vertices[selected_vertex_index_]) : -999;
  pho4_MVA_ = phoP4Corrected_dp_.size() > 3 ? phoP4Corrected_dp_[3].phoIdMvaDWrtVtx(Vertices[selected_vertex_index_]) : -999;

  pho12_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[1].p4();
  pho13_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[2].p4();
  pho14_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[3].p4();
  pho23_ = phoP4Corrected_dp_[1].p4() + phoP4Corrected_dp_[2].p4();
  pho24_ = phoP4Corrected_dp_[1].p4() + phoP4Corrected_dp_[3].p4();
  pho34_ = phoP4Corrected_dp_[2].p4() + phoP4Corrected_dp_[3].p4();

  flashgg::Photon pho1 = phoP4Corrected_dp_[diphoton_pairing_indices_.at(0)];
  flashgg::Photon pho2 = phoP4Corrected_dp_[diphoton_pairing_indices_.at(1)];
  flashgg::Photon pho3 = phoP4Corrected_dp_[diphoton_pairing_indices_.at(2)];
  flashgg::Photon pho4 = phoP4Corrected_dp_[diphoton_pairing_indices_.at(3)];
  if(pho1.pt() > pho2.pt()){
    dp1_pho1_ = pho1.p4();
    dp1_pho2_ = pho2.p4();
    dp1_ipho1_ = diphoton_pairing_indices_.at(0);
    dp1_ipho2_ = diphoton_pairing_indices_.at(1);
  }else{
    dp1_pho1_ = pho2.p4();
    dp1_pho2_ = pho3.p4();
    dp1_ipho1_ = diphoton_pairing_indices_.at(1);
    dp1_ipho2_ = diphoton_pairing_indices_.at(0);
  }
  if(pho3.pt() > pho4.pt()){
    dp2_pho1_ = pho3.p4();
    dp2_pho2_ = pho4.p4();
    dp2_ipho1_ = diphoton_pairing_indices_.at(2);
    dp2_ipho2_ = diphoton_pairing_indices_.at(3);
  }else{
    dp2_pho1_ = pho4.p4();
    dp2_pho2_ = pho3.p4();
    dp2_ipho1_ = diphoton_pairing_indices_.at(3);
    dp2_ipho2_ = diphoton_pairing_indices_.at(2);
  }
  auto dipho1 = pho1.p4() + pho2.p4();
  auto dipho2 = pho3.p4() + pho4.p4();
  if (dipho1.pt() > dipho2.pt())
  {
    dp1_ = dipho1;
    dp2_ = dipho2;
  }
  else if (dipho1.pt() < dipho2.pt())
  {
    dp1_ = dipho2;
    dp2_ = dipho1;
  }

  tp_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[1].p4() + phoP4Corrected_dp_[2].p4() + phoP4Corrected_dp_[3].p4();
  cout << "higgs mass " << tp_.mass() << endl;

  //--all do alternate diphoton pairing (making pairs with minimum deltaM)
  float minDM = 1000000;
  if (phoP4Corrected_dp_.size() > 3)
  {
    for (int i1=0; i1 < (int) phoP4Corrected_dp_.size(); i1++)
    {
      flashgg::Photon pho1_prime = phoP4Corrected_dp_[i1];
      for (int i2=0; i2 < (int) phoP4Corrected_dp_.size(); i2++)
      {
        if (i2 <= i1 ){continue;}
        flashgg::Photon pho2_prime = phoP4Corrected_dp_[i2];
        for (int i3=0; i3 < (int) phoP4Corrected_dp_.size(); i3++)
        {
          if (i3 == i2 || i3 == i1){continue;}
          flashgg::Photon pho3_prime = phoP4Corrected_dp_[i3];
          for (int i4=0; i4 < (int) phoP4Corrected_dp_.size(); i4++)
          {
            if (i4 <= i3){continue;}
            if (i4 == i1 || i4 == i2){continue;}
            flashgg::Photon pho4_prime = phoP4Corrected_dp_[i4];
            auto dipho1_prime = pho1_prime.p4() + pho2_prime.p4();
            auto dipho2_prime = pho3_prime.p4() + pho4_prime.p4();
            float deltaM = fabs( dipho1_prime.mass() - dipho2_prime.mass());
            if (deltaM < minDM){
              minDM = deltaM;
              dp1_pho1_prime_ = pho1_prime.p4();
              dp1_ipho1_prime_ = i1;
              dp1_pho2_prime_ = pho2_prime.p4();
              dp1_ipho2_prime_ = i2;
              dp2_pho1_prime_ = pho3_prime.p4();
              dp2_ipho1_prime_ = i3;
              dp2_pho2_prime_ = pho4_prime.p4();
              dp2_ipho2_prime_ = i4;
              if (dipho1_prime.pt() > dipho2_prime.pt())
              {
                dp1_prime_ = dipho1_prime;
                dp2_prime_ = dipho2_prime;
              }
              else if (dipho1_prime.pt() < dipho2_prime.pt())
              {
                dp1_prime_ = dipho2_prime;
                dp2_prime_ = dipho1_prime;
              }
            }
          }
        }
      }
    }
  }

}

H4GTag *H4GTag::clone() const
{
  H4GTag *result = new H4GTag(*this);
  return result;
}

float H4GTag::getCosThetaStar_CS() const {

  reco::Candidate::LorentzVector h_lor = dp1_ + dp2_;
  TLorentzVector h;
  h.SetPxPyPzE(h_lor.Px(),h_lor.Py(),h_lor.Pz(),h_lor.E()) ;

  reco::Candidate::LorentzVector a1_lor = dp1_;
  TLorentzVector a_1;
  a_1.SetPxPyPzE(a1_lor.Px(),a1_lor.Py(),a1_lor.Pz(),a1_lor.E()) ;

  a_1.Boost(-h.BoostVector());

  return a_1.CosTheta();
}

std::vector<float> H4GTag::CosThetaAngles() const {
  std::vector<float> helicityThetas;

  TLorentzVector Boosted_a1(0,0,0,0);
  Boosted_a1.SetPxPyPzE(dp1_.px(),dp1_.py(),dp1_.pz(),dp1_.energy()) ;
  TLorentzVector BoostedLeadingPhoton_a1(0,0,0,0);
  BoostedLeadingPhoton_a1.SetPxPyPzE(dp1_pho1_.px(),dp1_pho1_.py(),dp1_pho1_.pz(),dp1_pho1_.energy()) ;

  helicityThetas.push_back( HelicityCosTheta(Boosted_a1, BoostedLeadingPhoton_a1));

  TLorentzVector Boosted_a2(0,0,0,0);
  Boosted_a2.SetPxPyPzE(dp2_.px(),dp2_.py(),dp2_.pz(),dp2_.energy()) ;
  TLorentzVector BoostedLeadingPhoton_a2(0,0,0,0);
  BoostedLeadingPhoton_a2.SetPxPyPzE(dp2_pho1_.px(),dp2_pho1_.py(),dp2_pho1_.pz(),dp2_pho1_.energy()) ;

  helicityThetas.push_back( HelicityCosTheta(Boosted_a2, BoostedLeadingPhoton_a2));

  return helicityThetas;
}
//
// float H4GTag::HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const
// {
//   TVector3 BoostVector = Booster.BoostVector();
//
//   Boosted.Boost( -BoostVector.x(), -BoostVector.y(), -BoostVector.z() );
//   return Boosted.CosTheta();
// }



float H4GTag::getCosThetaStar_CS_prime() const {

  reco::Candidate::LorentzVector h_lor_prime = dp1_prime_ + dp2_prime_;
  TLorentzVector h_prime;
  h_prime.SetPxPyPzE(h_lor_prime.Px(),h_lor_prime.Py(),h_lor_prime.Pz(),h_lor_prime.E()) ;

  reco::Candidate::LorentzVector a1_lor_prime = dp1_prime_;
  TLorentzVector a_1_prime;
  a_1_prime.SetPxPyPzE(a1_lor_prime.Px(),a1_lor_prime.Py(),a1_lor_prime.Pz(),a1_lor_prime.E()) ;

  a_1_prime.Boost(-h_prime.BoostVector());

  return a_1_prime.CosTheta();
}

std::vector<float> H4GTag::CosThetaAngles_prime() const {
  std::vector<float> helicityThetas_prime;

  TLorentzVector Boosted_a1_prime(0,0,0,0);
  Boosted_a1_prime.SetPxPyPzE(dp1_prime_.px(),dp1_prime_.py(),dp1_prime_.pz(),dp1_prime_.energy()) ;
  TLorentzVector BoostedLeadingPhoton_a1_prime(0,0,0,0);
  BoostedLeadingPhoton_a1_prime.SetPxPyPzE(dp1_pho1_prime_.px(),dp1_pho1_prime_.py(),dp1_pho1_prime_.pz(),dp1_pho1_prime_.energy()) ;

  helicityThetas_prime.push_back( HelicityCosTheta(Boosted_a1_prime, BoostedLeadingPhoton_a1_prime));

  TLorentzVector Boosted_a2_prime(0,0,0,0);
  Boosted_a2_prime.SetPxPyPzE(dp2_prime_.px(),dp2_prime_.py(),dp2_prime_.pz(),dp2_prime_.energy()) ;
  TLorentzVector BoostedLeadingPhoton_a2_prime(0,0,0,0);
  BoostedLeadingPhoton_a2_prime.SetPxPyPzE(dp2_pho1_prime_.px(),dp2_pho1_prime_.py(),dp2_pho1_prime_.pz(),dp2_pho1_prime_.energy()) ;

  helicityThetas_prime.push_back( HelicityCosTheta(Boosted_a2_prime, BoostedLeadingPhoton_a2_prime));

  return helicityThetas_prime;
}

float H4GTag::HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const
{
  TVector3 BoostVector = Booster.BoostVector();

  Boosted.Boost( -BoostVector.x(), -BoostVector.y(), -BoostVector.z() );
  return Boosted.CosTheta();
}
