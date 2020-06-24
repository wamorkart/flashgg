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

{
}

H4GTag::~H4GTag() {}

H4GTag::H4GTag(edm::Ptr<DiPhotonCandidate> dipho, flashgg::Photon pho1, flashgg::Photon pho2, flashgg::Photon pho3, flashgg::Photon pho4, edm::Ptr<reco::Vertex> vertex_chosen, float dZ_bdtVtx)
{
  dipho_ = dipho;
  pho1_  = pho1;
  pho2_  = pho2;
  pho3_  = pho3;
  pho4_  = pho4;

  pho1_MVA_ = pho1.phoIdMvaDWrtVtx(vertex_chosen);
  pho2_MVA_ = pho2.phoIdMvaDWrtVtx(vertex_chosen);
  pho3_MVA_ = pho3.phoIdMvaDWrtVtx(vertex_chosen);
  pho4_MVA_ = pho4.phoIdMvaDWrtVtx(vertex_chosen);
  //  cout << pho1_MVA_ << "  " << pho2_MVA_ << "  " << pho3_MVA_ << "  " << pho4_MVA_ << endl;

  tp_ = pho1.p4()+pho2.p4()+pho3.p4()+pho4.p4();

  dZ_bdtVtx_ = dZ_bdtVtx;

  float minDM = 1000000;
  vector <flashgg::Photon> phoVect;
  phoVect.push_back(pho1);
  phoVect.push_back(pho2);
  phoVect.push_back(pho3);
  phoVect.push_back(pho4);

  for (int i1=0; i1 < (int) phoVect.size(); i1++)
  {
    flashgg::Photon pho1_prime = phoVect[i1];
    for (int i2=0; i2 < (int) phoVect.size(); i2++)
    {
      if (i2 <= i1 ){continue;}
      flashgg::Photon pho2_prime = phoVect[i2];
      for (int i3=0; i3 < (int) phoVect.size(); i3++)
      {
        if (i3 == i2 || i3 == i1){continue;}
        flashgg::Photon pho3_prime = phoVect[i3];
        for (int i4=0; i4 < (int) phoVect.size(); i4++)
        {
          if (i4 <= i3){continue;}
          if (i4 == i1 || i4 == i2){continue;}
          flashgg::Photon pho4_prime = phoVect[i4];
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

H4GTag::H4GTag(edm::Ptr<DiPhotonCandidate> dipho, flashgg::Photon pho1, flashgg::Photon pho2, flashgg::Photon pho3, edm::Ptr<reco::Vertex> vertex_chosen, float dZ_bdtVtx)
{
  dipho_ = dipho;
  pho1_  = pho1;
  pho2_  = pho2;
  pho3_  = pho3;

  pho1_MVA_ = pho1.phoIdMvaDWrtVtx(vertex_chosen);
  pho2_MVA_ = pho2.phoIdMvaDWrtVtx(vertex_chosen);
  pho3_MVA_ = pho3.phoIdMvaDWrtVtx(vertex_chosen);

  tp_ = pho1.p4()+pho2.p4()+pho3.p4();

  dZ_bdtVtx_ = dZ_bdtVtx;
}

H4GTag::H4GTag(edm::Ptr<DiPhotonCandidate> dipho, flashgg::Photon pho1, flashgg::Photon pho2, edm::Ptr<reco::Vertex> vertex_chosen, float dZ_bdtVtx)
{
  dipho_ = dipho;
  pho1_  = pho1;
  pho2_  = pho2;

  pho1_MVA_ = pho1.phoIdMvaDWrtVtx(vertex_chosen);
  pho2_MVA_ = pho2.phoIdMvaDWrtVtx(vertex_chosen);

  tp_ = pho1.p4()+pho2.p4();

  dZ_bdtVtx_ = dZ_bdtVtx;
}

float H4GTag::getCosThetaStar_CS() const {

  reco::Candidate::LorentzVector h_lor = dp1_prime_ + dp2_prime_;
  TLorentzVector h;
  h.SetPxPyPzE(h_lor.Px(),h_lor.Py(),h_lor.Pz(),h_lor.E()) ;

  reco::Candidate::LorentzVector a1_lor = dp1_prime_;
  TLorentzVector a_1;
  a_1.SetPxPyPzE(a1_lor.Px(),a1_lor.Py(),a1_lor.Pz(),a1_lor.E()) ;

  a_1.Boost(-h.BoostVector());

  return a_1.CosTheta();
}

std::vector<float> H4GTag::CosThetaAngles() const {
  std::vector<float> helicityThetas;

  TLorentzVector Boosted_a1(0,0,0,0);
  Boosted_a1.SetPxPyPzE(dp1_prime_.px(),dp1_prime_.py(),dp1_prime_.pz(),dp1_prime_.energy()) ;
  TLorentzVector BoostedLeadingPhoton_a1(0,0,0,0);
  BoostedLeadingPhoton_a1.SetPxPyPzE(dp1_pho1_prime_.px(),dp1_pho1_prime_.py(),dp1_pho1_prime_.pz(),dp1_pho1_prime_.energy()) ;

  helicityThetas.push_back( HelicityCosTheta(Boosted_a1, BoostedLeadingPhoton_a1));

  TLorentzVector Boosted_a2(0,0,0,0);
  Boosted_a2.SetPxPyPzE(dp2_prime_.px(),dp2_prime_.py(),dp2_prime_.pz(),dp2_prime_.energy()) ;
  TLorentzVector BoostedLeadingPhoton_a2(0,0,0,0);
  BoostedLeadingPhoton_a2.SetPxPyPzE(dp2_pho1_prime_.px(),dp2_pho1_prime_.py(),dp2_pho1_prime_.pz(),dp2_pho1_prime_.energy()) ;

  helicityThetas.push_back( HelicityCosTheta(Boosted_a2, BoostedLeadingPhoton_a2));

  return helicityThetas;
}

float H4GTag::HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const
{
  TVector3 BoostVector = Booster.BoostVector();

  Boosted.Boost( -BoostVector.x(), -BoostVector.y(), -BoostVector.z() );
  return Boosted.CosTheta();
}

H4GTag *H4GTag::clone() const
{
  H4GTag *result = new H4GTag(*this);
  return result;
}
