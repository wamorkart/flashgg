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

H4GTag::H4GTag(edm::Ptr<DiPhotonCandidate> dipho, flashgg::Photon pho1, flashgg::Photon pho2, flashgg::Photon pho3, flashgg::Photon pho4, edm::Ptr<reco::Vertex> vertex_chosen, float dZ_bdtVtx, float dZ_ZeroVtx, float dZ_HggVtx,   std::vector<reco::Candidate::LorentzVector > genPhos )
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

  pho1_full5x5_r9_ = pho1.full5x5_r9();
  pho2_full5x5_r9_ = pho2.full5x5_r9();
  pho3_full5x5_r9_ = pho3.full5x5_r9();
  pho4_full5x5_r9_ = pho4.full5x5_r9();

  pho1_egChargedHadronIso_ = pho1.egChargedHadronIso();
  pho2_egChargedHadronIso_ = pho2.egChargedHadronIso();
  pho3_egChargedHadronIso_ = pho3.egChargedHadronIso();
  pho4_egChargedHadronIso_ = pho4.egChargedHadronIso();

  pho1_hadronicOverEm_ = pho1.hadronicOverEm();
  pho2_hadronicOverEm_ = pho2.hadronicOverEm();
  pho3_hadronicOverEm_ = pho3.hadronicOverEm();
  pho4_hadronicOverEm_ = pho4.hadronicOverEm();

  pho1_SC_eta_ = pho1.superCluster()->eta();
  pho2_SC_eta_ = pho2.superCluster()->eta();
  pho3_SC_eta_ = pho3.superCluster()->eta();
  pho4_SC_eta_ = pho4.superCluster()->eta();


  // diphoPair_MVA_ = diphoPair_MVA;
  // dp1_ = dp1;
  // dp2_ = dp2;
  // dp1_pho1_ = dp1_pho1;
  // dp1_pho2_ = dp1_pho2;
  // dp2_pho1_ = dp2_pho1;
  // dp2_pho2_ = dp2_pho2;
  //
  // dp1_ipho1_ = dp1_ipho1;
  // dp1_ipho2_ = dp1_ipho2;
  // dp2_ipho1_ = dp2_ipho1;
  // dp2_ipho2_ = dp2_ipho2;



  tp_ = pho1.p4()+pho2.p4()+pho3.p4()+pho4.p4();

  dZ_bdtVtx_ = dZ_bdtVtx;
  dZ_ZeroVtx_ = dZ_ZeroVtx;
  dZ_HggVtx_ = dZ_HggVtx;

  float minDM = 1000000;
  vector <flashgg::Photon> phoVect;
  phoVect.push_back(pho1);
  phoVect.push_back(pho2);
  phoVect.push_back(pho3);
  phoVect.push_back(pho4);

  // -- dM pair method -- //

  //reco::Candidate::LorentzVector dp1_prime_, dp2_prime_, dp1_pho1_prime_, dp1_pho2_prime_, dp2_pho1_prime_,dp2_pho2_prime_   ;

  // min dM pairing begin-------
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

  cosThetaStarCS_prime_ = getCosThetaStar_CS(dp1_prime_, dp2_prime_);
  cosTheta_a1_prime_ = CosThetaAngles(dp1_prime_, dp2_prime_, dp1_pho1_prime_,dp2_pho1_prime_ )[0];
  cosTheta_a2_prime_ = CosThetaAngles(dp1_prime_, dp2_prime_, dp1_pho1_prime_,dp2_pho1_prime_ )[1];


  // -- BDT pair method -- //

  // cosThetaStarCS_ = getCosThetaStar_CS(dp1, dp2);
  // cosTheta_a1_ = CosThetaAngles(dp1, dp2, dp1_pho1, dp2_pho1 )[0];
  // cosTheta_a2_ = CosThetaAngles(dp1, dp2, dp1_pho1, dp2_pho1 )[1];

  pho12_ = pho1.p4() + pho2.p4();
  pho13_ = pho1.p4() + pho3.p4();
  pho14_ = pho1.p4() + pho4.p4();
  pho23_ = pho2.p4() + pho3.p4();
  pho24_ = pho2.p4() + pho4.p4();
  pho34_ = pho3.p4() + pho4.p4();

  if (genPhos.size() == 4)
  {
    gen_pho1_pt_ = genPhos[0].pt();
    gen_pho2_pt_ = genPhos[1].pt();
    gen_pho3_pt_ = genPhos[2].pt();
    gen_pho4_pt_ = genPhos[3].pt();
    gen_pho1_eta_ = genPhos[0].eta();
    gen_pho2_eta_ = genPhos[1].eta();
    gen_pho3_eta_ = genPhos[2].eta();
    gen_pho4_eta_ = genPhos[3].eta();

    gen_pho12_dR_ = deltaR(genPhos[0].eta(), genPhos[0].phi(), genPhos[1].eta(), genPhos[1].phi());
    gen_pho13_dR_ = deltaR(genPhos[0].eta(), genPhos[0].phi(), genPhos[2].eta(), genPhos[2].phi());
    gen_pho14_dR_ = deltaR(genPhos[0].eta(), genPhos[0].phi(), genPhos[3].eta(), genPhos[3].phi());
    gen_pho23_dR_ = deltaR(genPhos[1].eta(), genPhos[1].phi(), genPhos[2].eta(), genPhos[2].phi());
    gen_pho24_dR_ = deltaR(genPhos[1].eta(), genPhos[1].phi(), genPhos[3].eta(), genPhos[3].phi());
    gen_pho34_dR_ = deltaR(genPhos[2].eta(), genPhos[2].phi(), genPhos[3].eta(), genPhos[3].phi());

    gen_pho12_M_ = (genPhos[0]+genPhos[1]).M();
    gen_pho13_M_ = (genPhos[0]+genPhos[2]).M();
    gen_pho14_M_ = (genPhos[0]+genPhos[3]).M();
    gen_pho23_M_ = (genPhos[1]+genPhos[2]).M();
    gen_pho24_M_ = (genPhos[1]+genPhos[3]).M();
    gen_pho34_M_ = (genPhos[2]+genPhos[3]).M();

    gen_a1_pt_ = (genPhos[0]+genPhos[1]).pt();
    gen_a2_pt_ = (genPhos[2]+genPhos[3]).pt();
    gen_a1_eta_ = (genPhos[0]+genPhos[1]).eta();
    gen_a2_eta_ = (genPhos[2]+genPhos[3]).eta();
    gen_a1a2_dR_ = deltaR( (genPhos[0]+genPhos[1]).eta(), (genPhos[0]+genPhos[1]).phi(),(genPhos[2]+genPhos[3]).eta(), (genPhos[2]+genPhos[3]).phi()  );
    gen_h_mass_ = (genPhos[0]+genPhos[1]+genPhos[2]+genPhos[3]).M();
    gen_h_pt_ = (genPhos[0]+genPhos[1]+genPhos[2]+genPhos[3]).pt();
    gen_h_eta_ = (genPhos[0]+genPhos[1]+genPhos[2]+genPhos[3]).eta();

  }




}

H4GTag::H4GTag(edm::Ptr<DiPhotonCandidate> dipho, flashgg::Photon pho1, flashgg::Photon pho2, flashgg::Photon pho3, edm::Ptr<reco::Vertex> vertex_chosen, float dZ_bdtVtx, float dZ_ZeroVtx, float dZ_HggVtx,   std::vector<reco::Candidate::LorentzVector > genPhos)
{
  dipho_ = dipho;
  pho1_  = pho1;
  pho2_  = pho2;
  pho3_  = pho3;

  pho1_MVA_ = pho1.phoIdMvaDWrtVtx(vertex_chosen);
  pho2_MVA_ = pho2.phoIdMvaDWrtVtx(vertex_chosen);
  pho3_MVA_ = pho3.phoIdMvaDWrtVtx(vertex_chosen);

  pho1_full5x5_r9_ = pho1.full5x5_r9();
  pho2_full5x5_r9_ = pho2.full5x5_r9();
  pho3_full5x5_r9_ = pho3.full5x5_r9();

  pho1_egChargedHadronIso_ = pho1.egChargedHadronIso();
  pho2_egChargedHadronIso_ = pho2.egChargedHadronIso();
  pho3_egChargedHadronIso_ = pho3.egChargedHadronIso();

  pho1_hadronicOverEm_ = pho1.hadronicOverEm();
  pho2_hadronicOverEm_ = pho2.hadronicOverEm();
  pho3_hadronicOverEm_ = pho3.hadronicOverEm();

  pho1_SC_eta_ = pho1.superCluster()->eta();
  pho2_SC_eta_ = pho2.superCluster()->eta();
  pho3_SC_eta_ = pho3.superCluster()->eta();

  tp_ = pho1.p4()+pho2.p4()+pho3.p4();

  dZ_bdtVtx_ = dZ_bdtVtx;
  dZ_ZeroVtx_ =dZ_ZeroVtx;

  pho12_ = pho1.p4() + pho2.p4();
  pho13_ = pho1.p4() + pho3.p4();
  pho23_ = pho2.p4() + pho3.p4();
}

H4GTag::H4GTag(edm::Ptr<DiPhotonCandidate> dipho, flashgg::Photon pho1, flashgg::Photon pho2, edm::Ptr<reco::Vertex> vertex_chosen, float dZ_bdtVtx, float dZ_ZeroVtx, float dZ_HggVtx,   std::vector<reco::Candidate::LorentzVector > genPhos)
{
  dipho_ = dipho;
  pho1_  = pho1;
  pho2_  = pho2;

  pho1_MVA_ = pho1.phoIdMvaDWrtVtx(vertex_chosen);
  pho2_MVA_ = pho2.phoIdMvaDWrtVtx(vertex_chosen);

  pho1_full5x5_r9_ = pho1.full5x5_r9();
  pho2_full5x5_r9_ = pho2.full5x5_r9();

  pho1_egChargedHadronIso_ = pho1.egChargedHadronIso();
  pho2_egChargedHadronIso_ = pho2.egChargedHadronIso();

  pho1_hadronicOverEm_ = pho1.hadronicOverEm();
  pho2_hadronicOverEm_ = pho2.hadronicOverEm();

  pho1_SC_eta_ = pho1.superCluster()->eta();
  pho2_SC_eta_ = pho2.superCluster()->eta();

  tp_ = pho1.p4()+pho2.p4();

  dZ_bdtVtx_ = dZ_bdtVtx;
  dZ_ZeroVtx_ = dZ_ZeroVtx;
}

float H4GTag::getCosThetaStar_CS(reco::Candidate::LorentzVector a1, reco::Candidate::LorentzVector a2) const {

  reco::Candidate::LorentzVector h_lor = a1 + a2;
  TLorentzVector h;
  h.SetPxPyPzE(h_lor.Px(),h_lor.Py(),h_lor.Pz(),h_lor.E()) ;

  reco::Candidate::LorentzVector a1_lor = a1;
  TLorentzVector a_1;
  a_1.SetPxPyPzE(a1_lor.Px(),a1_lor.Py(),a1_lor.Pz(),a1_lor.E()) ;

  a_1.Boost(-h.BoostVector());

  return a_1.CosTheta();
}

std::vector<float> H4GTag::CosThetaAngles(reco::Candidate::LorentzVector a1, reco::Candidate::LorentzVector a2, reco::Candidate::LorentzVector a1_pho1, reco::Candidate::LorentzVector a2_pho1) const {
  std::vector<float> helicityThetas;

  TLorentzVector Boosted_a1(0,0,0,0);
  Boosted_a1.SetPxPyPzE(a1.px(),a1.py(),a1.pz(),a1.energy()) ;
  TLorentzVector BoostedLeadingPhoton_a1(0,0,0,0);
  BoostedLeadingPhoton_a1.SetPxPyPzE(a1_pho1.px(),a1_pho1.py(),a1_pho1.pz(),a1_pho1.energy()) ;

  helicityThetas.push_back( HelicityCosTheta(Boosted_a1, BoostedLeadingPhoton_a1));

  TLorentzVector Boosted_a2(0,0,0,0);
  Boosted_a2.SetPxPyPzE(a2.px(),a2.py(),a2.pz(),a2.energy()) ;
  TLorentzVector BoostedLeadingPhoton_a2(0,0,0,0);
  BoostedLeadingPhoton_a2.SetPxPyPzE(a2_pho1.px(),a2_pho1.py(),a2_pho1.pz(),a2_pho1.energy()) ;

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
