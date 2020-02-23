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

    H4GTag(edm::Ptr<DiPhotonCandidate>,vector<edm::Ptr<reco::Vertex>>,edm::Ptr<reco::Vertex>,reco::GenParticle::Point,int,TMVA::Reader *,vector< flashgg::Photon>,vector<int>, float );
    virtual H4GTag *clone() const override;

    const float pho1_MVA() const { return pho1_MVA_; };
    const float pho2_MVA() const { return pho2_MVA_; };
    const float pho3_MVA() const { return pho3_MVA_; };
    const float pho4_MVA() const { return pho4_MVA_; };

  private:
    float BS_factor_BDTVtx_;
    float pho1_MVA_;
    float pho2_MVA_;
    float pho3_MVA_;
    float pho4_MVA_;
  };
}

#endif
