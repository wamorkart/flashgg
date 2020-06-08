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
H4GTag::H4GTag(edm::Ptr<DiPhotonCandidate> dipho, vector< flashgg::Photon>phoP4Corrected_dp)
{
  dipho_ = dipho;
  phoP4Corrected_dp_  = phoP4Corrected_dp;
}

H4GTag *H4GTag::clone() const
{
  H4GTag *result = new H4GTag(*this);
  return result;
}
