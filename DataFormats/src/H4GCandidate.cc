#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "flashgg/DataFormats/interface/VertexCandidateMap.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "flashgg/MicroAOD/interface/CutBasedDiPhotonObjectSelector.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "flashgg/DataFormats/interface/H4GCandidate.h"

#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "flashgg/Taggers/interface/FunctionHelpers.h"

#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "TMVA/Reader.h"
#include <algorithm>

using namespace flashgg;

H4GCandidate::H4GCandidate():
Vertices_(),
Vertex_random_(),
vertex_diphoton_ (),
BS_factor_0Vtx_ (),
BS_factor_HggVtx_(),
BS_factor_RandomVtx_ (),
BS_factor_BDTVtx_ (),
phoP4Corrected_dp_(),
phoP4Corrected1_ (),
phoP4Corrected2_ (),
phoP4Corrected3_ (),
phoP4Corrected4_ (),
pho1_MVA_ (),
pho2_MVA_ (),
pho3_MVA_ (),
pho4_MVA_ (),
pho_MVA_min_(),
pho_MVA_max_(),
dp1_ (),
dp2_ (),
dp1_pho1_ (),
dp1_pho2_ (),
dp2_pho1_ (),
dp2_pho2_ (),
dp1_ipho1_ (),
dp1_ipho2_ (),
dp2_ipho1_ (),
dp2_ipho2_ (),
pho12_ (),
pho13_ (),
pho14_ (),
pho23_ (),
pho24_ (),
pho34_ (),
tp_ (),
isSR_(),
isCR_(),
gen_a1_mass_(),
gen_a2_mass_(),
gen_h_mass_()
// diphoPhotonsVector_()
{}
  H4GCandidate::~H4GCandidate() {}
  H4GCandidate::H4GCandidate( std::vector<edm::Ptr<reco::Vertex>> Vertices, std::vector<edm::Ptr<reco::Vertex>> slim_Vertices, edm::Ptr<reco::Vertex> vertex_diphoton, edm::Ptr<reco::Vertex> vertex_bdt, reco::GenParticle::Point genVertex, math::XYZPoint BSPoint, std::vector<std::vector<float>> Vector, float MVA0, float MVA1, float MVA2, float dZ1, float dZ2, float dZtrue, int hgg_index, int trueVtx_index, int rndVtx_index, int bdtVtx_index, float tp_pt, float nVertices, float nConv, TMVA::Reader *VertexProbMva, double genTotalWeight, std::vector<flashgg::Photon> diphoPhotons, std::vector <edm::Ptr<flashgg::DiPhotonCandidate>> diPhoPtrs,float gen_a1_mass, float gen_a2_mass, float gen_h_mass):
   Vertices_(Vertices), slim_Vertices_(slim_Vertices),vertex_diphoton_(vertex_diphoton), vertex_bdt_(vertex_bdt), genVertex_(genVertex), BSPoint_(BSPoint), Vector_(Vector), MVA0_(MVA0), MVA1_(MVA1), MVA2_(MVA2), dZ1_(dZ1), dZ2_(dZ2), dZtrue_(dZtrue), hgg_index_(hgg_index), trueVtx_index_(trueVtx_index), rndVtx_index_(rndVtx_index), bdtVtx_index_(bdtVtx_index), tp_pt_(tp_pt), nVertices_(nVertices), nConv_(nConv), VertexProbMva_(VertexProbMva), genTotalWeight_(genTotalWeight),diphoPhotons_(diphoPhotons), diPhoPtrs_(diPhoPtrs), gen_a1_mass_(gen_a1_mass), gen_a2_mass_(gen_a2_mass), gen_h_mass_(gen_h_mass)

  {

    isSR_ = 0;
    isCR_ = 0;
    int random_vtx = rand() % slim_Vertices_.size();
    Vertex_random_ = slim_Vertices_[random_vtx];

    float vtx_X = Vertices_[bdtVtx_index]->x();
    float vtx_Y = Vertices_[bdtVtx_index]->y();
    float vtx_Z = Vertices_[bdtVtx_index]->z();
    // float vtx_X = Vertices_[0]->x();
    // float vtx_Y = Vertices_[0]->y();
    // float vtx_Z = Vertices_[0]->z();

    //--Beam spot reweighting (https://github.com/cms-analysis/flashggFinalFit/blob/e60d53e19ac4f20e7ce187f0a34e483b4fc2a60e/Signal/test/SignalFit.cpp)
    float mcBeamSpotWidth_=5.14; //cm
    float dataBeamSpotWidth_=3.5; //cm

    float dZ_HggVtx = genVertex_.z() - vertex_diphoton_->z();
    float dZ_0Vtx = genVertex.z() - Vertices_[0]->z();
    float dZ_RandomVtx = genVertex.z() - Vertex_random_->z();
    float dZ_BDTVtx = genVertex.z() - Vertices_[bdtVtx_index]->z();;


    if (fabs(dZ_HggVtx) < 0.1 ){
      BS_factor_HggVtx_ =1;
    } else {
      double mcBeamSpot_HggVtx=TMath::Gaus(dZ_HggVtx,0,TMath::Sqrt(2)*mcBeamSpotWidth_,true);
      double dataBeamSpot_HggVtx=TMath::Gaus(dZ_HggVtx,0,TMath::Sqrt(2)*dataBeamSpotWidth_,true);
      BS_factor_HggVtx_ = dataBeamSpot_HggVtx/mcBeamSpot_HggVtx;
    }

    if (fabs(dZ_0Vtx) < 0.1 ){
      BS_factor_0Vtx_ =1;
    } else {
      double mcBeamSpot_0Vtx=TMath::Gaus(dZ_0Vtx,0,TMath::Sqrt(2)*mcBeamSpotWidth_,true);
      double dataBeamSpot_0Vtx=TMath::Gaus(dZ_0Vtx,0,TMath::Sqrt(2)*dataBeamSpotWidth_,true);
      BS_factor_0Vtx_ = dataBeamSpot_0Vtx/mcBeamSpot_0Vtx;
    }

    if (fabs(dZ_RandomVtx) < 0.1 ){
      BS_factor_RandomVtx_ =1;
    } else {
      double mcBeamSpot_RandomVtx=TMath::Gaus(dZ_RandomVtx,0,TMath::Sqrt(2)*mcBeamSpotWidth_,true);
      double dataBeamSpot_RandomVtx=TMath::Gaus(dZ_RandomVtx,0,TMath::Sqrt(2)*dataBeamSpotWidth_,true);
      BS_factor_RandomVtx_ = dataBeamSpot_RandomVtx/mcBeamSpot_RandomVtx;
    }

    if (fabs(dZ_BDTVtx) < 0.1 ){
      BS_factor_BDTVtx_ =1;
    } else {
      double mcBeamSpot_BDTVtx=TMath::Gaus(dZ_BDTVtx,0,TMath::Sqrt(2)*mcBeamSpotWidth_,true);
      double dataBeamSpot_BDTVtx=TMath::Gaus(dZ_BDTVtx,0,TMath::Sqrt(2)*dataBeamSpotWidth_,true);
      BS_factor_BDTVtx_ = dataBeamSpot_BDTVtx/mcBeamSpot_BDTVtx;
    }

    math::XYZVector vtx_Pos( vtx_X, vtx_Y, vtx_Z );

    if (diphoPhotons_.size() > 0)
    {
      for (int dp = 0; dp < (int) diphoPhotons_.size(); dp++)
      {
        float sc_X_dp = diphoPhotons_[dp].superCluster()->x();
        float sc_Y_dp = diphoPhotons_[dp].superCluster()->y();
        float sc_Z_dp = diphoPhotons_[dp].superCluster()->z();
        math::XYZVector sc_Pos_dp( sc_X_dp, sc_Y_dp, sc_Z_dp );
        math::XYZVector direction_dp = sc_Pos_dp - vtx_Pos;
        math::XYZVector pho_dp = ( direction_dp.Unit() ) * ( diphoPhotons_[dp].energy() );
        math::XYZTLorentzVector corrected_p4_dp( pho_dp.x(), pho_dp.y(), pho_dp.z(), diphoPhotons_[dp].energy() );
        diphoPhotons_[dp].setP4(corrected_p4_dp);
        phoP4Corrected_dp_.push_back(diphoPhotons_[dp]);
      }
    }

    std::sort(phoP4Corrected_dp_.begin(), phoP4Corrected_dp_.end(), [](const flashgg::Photon a, const flashgg::Photon  b) {return a.pt() > b.pt(); });

    pho1_MVA_ = phoP4Corrected_dp_.size() > 0 ? phoP4Corrected_dp_[0].phoIdMvaDWrtVtx(Vertices_[bdtVtx_index]) : -999;
    pho2_MVA_ = phoP4Corrected_dp_.size() > 0 ? phoP4Corrected_dp_[1].phoIdMvaDWrtVtx(Vertices_[bdtVtx_index]) : -999;
    pho3_MVA_ = phoP4Corrected_dp_.size() > 2 ? phoP4Corrected_dp_[2].phoIdMvaDWrtVtx(Vertices_[bdtVtx_index]) : -999;
    pho4_MVA_ = phoP4Corrected_dp_.size() > 3 ? phoP4Corrected_dp_[3].phoIdMvaDWrtVtx(Vertices_[bdtVtx_index]) : -999;

    std::vector<float> pho_MVA_vec;
    pho_MVA_vec.push_back(pho1_MVA_);
    pho_MVA_vec.push_back(pho2_MVA_);
    pho_MVA_vec.push_back(pho3_MVA_);
    pho_MVA_vec.push_back(pho4_MVA_);

    pho_MVA_min_ = *min_element(pho_MVA_vec.begin(),pho_MVA_vec.end());
    pho_MVA_max_ = *max_element(pho_MVA_vec.begin(),pho_MVA_vec.end());

    float minDM = 1000000;
    if (phoP4Corrected_dp_.size() > 3)
    {
      for (int i1=0; i1 < (int) phoP4Corrected_dp_.size(); i1++)
      {
        flashgg::Photon pho1 = phoP4Corrected_dp_[i1];
        for (int i2=0; i2 < (int) phoP4Corrected_dp_.size(); i2++)
        {
          if (i2 <= i1 ){continue;}
          flashgg::Photon pho2 = phoP4Corrected_dp_[i2];
          for (int i3=0; i3 < (int) phoP4Corrected_dp_.size(); i3++)
          {
            if (i3 == i2 || i3 == i1){continue;}
            flashgg::Photon pho3 = phoP4Corrected_dp_[i3];
            for (int i4=0; i4 < (int) phoP4Corrected_dp_.size(); i4++)
            {
              if (i4 <= i3){continue;}
              if (i4 == i1 || i4 == i2){continue;}
              flashgg::Photon pho4 = phoP4Corrected_dp_[i4];
              auto dipho1 = pho1.p4() + pho2.p4();
              auto dipho2 = pho3.p4() + pho4.p4();
              float deltaM = fabs( dipho1.mass() - dipho2.mass());
              if (deltaM < minDM){
                minDM = deltaM;
                dp1_pho1_ = pho1.p4();
                dp1_ipho1_ = i1;
                dp1_pho2_ = pho2.p4();
                dp1_ipho2_ = i2;
                dp2_pho1_ = pho3.p4();
                dp2_ipho1_ = i3;
                dp2_pho2_ = pho4.p4();
                dp2_ipho2_ = i4;
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
              }
            }
          }
        }
      }

      tp_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[1].p4() + phoP4Corrected_dp_[2].p4() + phoP4Corrected_dp_[3].p4();
      pho12_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[1].p4();
      pho13_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[2].p4();
      pho14_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[3].p4();
      pho23_ = phoP4Corrected_dp_[1].p4() + phoP4Corrected_dp_[2].p4();
      pho24_ = phoP4Corrected_dp_[1].p4() + phoP4Corrected_dp_[3].p4();
      pho34_ = phoP4Corrected_dp_[2].p4() + phoP4Corrected_dp_[3].p4();

      phoP4Corrected1_ = phoP4Corrected_dp_[dp1_ipho1_];
      phoP4Corrected2_ = phoP4Corrected_dp_[dp1_ipho2_];
      phoP4Corrected3_ = phoP4Corrected_dp_[dp2_ipho1_];
      phoP4Corrected4_ = phoP4Corrected_dp_[dp2_ipho2_];

      vector<flashgg::Photon> selPhotons;
      selPhotons = PhotonBasicSelection (phoP4Corrected_dp_);
      vector<flashgg::Photon> IDPhotons;
      std::vector<double> _PhotonID = {-0.9,-0.9,0.6,0.6};
      if (selPhotons.size() == 4)
         {
           std::vector<flashgg::Photon> goodIDPhotons;
           for (unsigned int phoIndex = 0; phoIndex < selPhotons.size(); phoIndex++)
           {
             flashgg::Photon photon = selPhotons[phoIndex];
             if (photon.phoIdMvaDWrtVtx(Vertices_[bdtVtx_index]) < _PhotonID[phoIndex]) continue;
             goodIDPhotons.push_back(photon);
           }
           if (goodIDPhotons.size() == 4)
           {
             isSR_=1;
             isCR_=0;
            }
           if (goodIDPhotons.size() == 2)
           {
            isCR_=1;
            isSR_=0;
            }
          }
    }

    if (phoP4Corrected_dp_.size() == 2)
    {
      dp1_ipho1_ = 0;
      dp1_ipho2_ = 1;
      dp2_ipho1_ = -1;
      dp2_ipho2_ = -1;

      tp_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[1].p4();
      pho12_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[1].p4();


      phoP4Corrected1_ = phoP4Corrected_dp_[dp1_ipho1_];
      phoP4Corrected2_ = phoP4Corrected_dp_[dp1_ipho2_];
    }
    else if (phoP4Corrected_dp_.size() == 3)
    {
      dp1_ipho1_ = 0;
      dp1_ipho2_ = 1;
      dp2_ipho1_ = 2;
      dp2_ipho2_ = -1;

      tp_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[1].p4() + phoP4Corrected_dp_[2].p4();
      pho12_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[1].p4();
      pho13_ = phoP4Corrected_dp_[0].p4() + phoP4Corrected_dp_[2].p4();
      pho23_ = phoP4Corrected_dp_[1].p4() + phoP4Corrected_dp_[2].p4();


      phoP4Corrected1_ = phoP4Corrected_dp_[dp1_ipho1_];
      phoP4Corrected2_ = phoP4Corrected_dp_[dp1_ipho2_];
      phoP4Corrected3_ = phoP4Corrected_dp_[dp2_ipho1_];
    }

    tp_pt =  tp_.pt();
    vtxProbMVA_ = VertexProbMva_->EvaluateMVA( "BDT" );
  }

  // good photon selection

  std::vector<flashgg::Photon> H4GCandidate::PhotonBasicSelection (std::vector<flashgg::Photon> phoCol) const
  {
    std::vector<double> _PhotonPt = {30.,18.,10.,10};
    std::vector<flashgg::Photon> goodPhotons;
    for (unsigned int phoIndex = 0; phoIndex < phoCol.size(); phoIndex++)
    {
      flashgg::Photon photon = phoCol[phoIndex];
      if (photon.pt() < _PhotonPt[phoIndex] || abs(photon.eta() > 2.5 || photon.passElectronVeto() !=1 )) continue;
      goodPhotons.push_back(photon);
    }
    return goodPhotons;
  }

  float H4GCandidate::getCosThetaStar_CS() const {

    reco::Candidate::LorentzVector h_lor = dp1_ + dp2_;
    TLorentzVector h;
    h.SetPxPyPzE(h_lor.Px(),h_lor.Py(),h_lor.Pz(),h_lor.E()) ;

    reco::Candidate::LorentzVector a1_lor = dp1_;
    TLorentzVector a_1;
    a_1.SetPxPyPzE(a1_lor.Px(),a1_lor.Py(),a1_lor.Pz(),a1_lor.E()) ;

    a_1.Boost(-h.BoostVector());

    return a_1.CosTheta();
  }

  float H4GCandidate::getCosThetaStar_CS_old(float ebeam) const {
    TLorentzVector p1, p2;
    p1.SetPxPyPzE(0, 0,  ebeam, ebeam);
    p2.SetPxPyPzE(0, 0, -ebeam, ebeam);

    reco::Candidate::LorentzVector h_lor = dp1_ + dp2_;
    TLorentzVector h;
    h.SetPxPyPzE(h_lor.Px(),h_lor.Py(),h_lor.Pz(),h_lor.E()) ;

    TVector3 boost = - h.BoostVector();
    p1.Boost(boost);
    p2.Boost(boost);
    reco::Candidate::LorentzVector a1_lor = dp1_;
    TLorentzVector a_1;
    a_1.SetPxPyPzE(a1_lor.Px(),a1_lor.Py(),a1_lor.Pz(),a1_lor.E()) ;
    a_1.Boost(boost);

    TVector3 CSaxis = p1.Vect().Unit() - p2.Vect().Unit();
    CSaxis.Unit();
    return cos(   CSaxis.Angle( a_1.Vect().Unit() )    );
  }
  std::vector<float> H4GCandidate::CosThetaAngles() const {
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

  float H4GCandidate::HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const
  {
    TVector3 BoostVector = Booster.BoostVector();

    Boosted.Boost( -BoostVector.x(), -BoostVector.y(), -BoostVector.z() );
    return Boosted.CosTheta();
  }
