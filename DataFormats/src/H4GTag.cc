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

using namespace flashgg; // makes flashgg sub members visible
H4GTag::H4GTag() : DiPhotonTagBase::DiPhotonTagBase(), mva_(-2.), JetVector_ (), Cut_Variables_ ()
{

}

H4GTag::~H4GTag() {}

void H4GTag::GetPhotons(edm::Ptr<DiPhotonCandidate> dipho)
{
  // dipho->makePhotonsPersistent();
  // auto diphop4 = dipho_.p4();

  // Get photons
  // flashgg::Photon leading_photon = dipho->getLeadingPhoton();
  // flashgg::Photon subleading_photon = dipho->getSubLeadingPhoton();

  // Save as dumper objects
  Leading_Photon_ = dipho->leadingPhoton();
  Subleading_Photon_ = dipho->subLeadingPhoton();

  lp_Hgg_MVA_ = Leading_Photon_->phoIdMvaDWrtVtx( dipho->vtx() );
  slp_Hgg_MVA_ = Subleading_Photon_->phoIdMvaDWrtVtx( dipho->vtx() );
}

// 2 jets, electron
H4GTag::H4GTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron, edm::Ptr<flashgg::Met> MET,
                      edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2, std::vector<flashgg::Jet> tagJets_, std::vector<double> Cut_Variables) : JetVector_(tagJets_), Cut_Variables_(Cut_Variables)
{
  // flashgg::DiPhotonCandidate dipho_ = diphoVector_[0];
  // // if (!dipho_zero_vtx){
  // //   cout << "Diphoton vertex is not zero. Need to recompute" << endl;
  // // }

  // leading_dpho_ = dipho;

  // // Get photons

  dipho_ = dipho;
  GetPhotons(dipho);
}

// 1 jet, electron
// HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron, edm::Ptr<flashgg::Met> MET,
//                       edm::Ptr<flashgg::Jet> jet, std::vector<flashgg::Jet> tagJets_, std::vector<double> Cut_Variables) : JetVector_(tagJets_), Cut_Variables_(Cut_Variables)
// {
//   dipho_ = dipho;
//   GetPhotons(dipho);
// }

// two jets, muon
// HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Muon> muon, edm::Ptr<flashgg::Met> MET,
//                       edm::Ptr<flashgg::Jet> jet1, edm::Ptr<flashgg::Jet> jet2, std::vector<flashgg::Jet> tagJets_, std::vector<double> Cut_Variables) : JetVector_(tagJets_), Cut_Variables_(Cut_Variables)
// {
//   dipho_ = dipho;
//   GetPhotons(dipho);
// }

// 1 jet, muon
// HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Muon> muon, edm::Ptr<flashgg::Met> MET,
//                       edm::Ptr<flashgg::Jet> jet, std::vector<flashgg::Jet> tagJets_, std::vector<double> Cut_Variables) : JetVector_(tagJets_), Cut_Variables_(Cut_Variables)
// {
//   dipho_ = dipho;
//   GetPhotons(dipho);
// }

// Untagged
// HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, std::vector<flashgg::Jet> tagJets_, std::vector<double> Cut_Variables) : JetVector_(tagJets_), Cut_Variables_(Cut_Variables)
// {
//   dipho_ = dipho;
//   GetPhotons(dipho);
// }

// HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Electron> electron, edm::Ptr<flashgg::Met> MET)
// {
//   dipho_ = dipho;
//   GetPhotons(dipho);
// }

// HHWWggTag::HHWWggTag(edm::Ptr<DiPhotonCandidate> dipho, edm::Ptr<flashgg::Muon> muon, edm::Ptr<flashgg::Met> MET)
// {
//   dipho_ = dipho;
//   GetPhotons(dipho);
// }

// You need this because HHWWggTag is derived from another class
H4GTag *H4GTag::clone() const
{
    H4GTag *result = new H4GTag( *this );
    return result;
}
