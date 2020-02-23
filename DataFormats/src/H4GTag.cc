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

H4GTag::H4GTag(edm::Ptr<DiPhotonCandidate> dipho,vector<edm::Ptr<reco::Vertex>> Vertices,edm::Ptr<reco::Vertex>vertex_bdt,reco::GenParticle::Point genVertex, int selected_vertex_index_, TMVA::Reader *VertexProbMva_,vector< flashgg::Photon>phoP4Corrected_dp,vector<int> diphoton_pairing_indices_, float diphoPair_MVA   )
{
  dipho_ = dipho;
  cout << "dipho pt " << dipho->pt() << endl;

  // float vtx_X = Vertices[vertex_bdt]->x();
  // float vtx_Y = Vertices[vertex_bdt]->y();
  // float vtx_Z = Vertices[vertex_bdt]->z();

  //--Beam spot reweighting (https://github.com/cms-analysis/flashggFinalFit/blob/e60d53e19ac4f20e7ce187f0a34e483b4fc2a60e/Signal/test/SignalFit.cpp)
  float mcBeamSpotWidth_=5.14; //cm
  float dataBeamSpotWidth_=3.5; //cm

  float dZ_BDTVtx = genVertex.z() - Vertices[selected_vertex_index_]->z();;

  if (fabs(dZ_BDTVtx) < 0.1 ){
    BS_factor_BDTVtx_ =1;
  } else {
    double mcBeamSpot_BDTVtx=TMath::Gaus(dZ_BDTVtx,0,TMath::Sqrt(2)*mcBeamSpotWidth_,true);
    double dataBeamSpot_BDTVtx=TMath::Gaus(dZ_BDTVtx,0,TMath::Sqrt(2)*dataBeamSpotWidth_,true);
    BS_factor_BDTVtx_ = dataBeamSpot_BDTVtx/mcBeamSpot_BDTVtx;
  }
  pho1_MVA_ = phoP4Corrected_dp.size() > 0 ? phoP4Corrected_dp[0].phoIdMvaDWrtVtx(Vertices[selected_vertex_index_]) : -999;
  pho2_MVA_ = phoP4Corrected_dp.size() > 0 ? phoP4Corrected_dp[1].phoIdMvaDWrtVtx(Vertices[selected_vertex_index_]) : -999;
  pho3_MVA_ = phoP4Corrected_dp.size() > 2 ? phoP4Corrected_dp[2].phoIdMvaDWrtVtx(Vertices[selected_vertex_index_]) : -999;
  pho4_MVA_ = phoP4Corrected_dp.size() > 3 ? phoP4Corrected_dp[3].phoIdMvaDWrtVtx(Vertices[selected_vertex_index_]) : -999;

}

H4GTag *H4GTag::clone() const
{
  H4GTag *result = new H4GTag(*this);
  return result;
}
