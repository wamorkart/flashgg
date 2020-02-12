// Abe Tishelman-Charny
// November 2019
// Derived from HH->WWgg event dumper and HH->bbgg tagger

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

#include "FWCore/Common/interface/TriggerNames.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/SinglePhotonView.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/H4GTag.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "flashgg/MicroAOD/interface/CutBasedDiPhotonObjectSelector.h"

#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
#include "flashgg/DataFormats/interface/WHLeptonicTag.h"

#include "flashgg/Taggers/interface/LeptonSelection.h"

#include <vector>
#include <algorithm>
#include "TGraph.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;

namespace flashgg {
  class H4GTagProducer : public EDProducer
  {
  public:
    typedef math::XYZTLorentzVector LorentzVector;

    H4GTagProducer( const ParameterSet &);

  private:
    void produce( Event &, const EventSetup &) override;
    std::vector<edm::EDGetTokenT<edm::View<DiPhotonCandidate> > > diPhotonTokens_;
    std::string inputDiPhotonName_;

    EDGetTokenT<View<DiPhotonCandidate> > diphotonToken_;
    Handle<View<flashgg::DiPhotonCandidate> > diphotons;

    EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
    Handle<View<reco::GenParticle> > genParticle;

    std::vector< std::string > systematicsLabels;
    std::vector<std::string> inputDiPhotonSuffixes_;
    string systLabel_;
  };
  H4GTagProducer::H4GTagProducer( const ParameterSet & pSet):
  diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( pSet.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
  genParticleToken_( consumes<View<reco::GenParticle> >( pSet.getParameter<InputTag> ( "GenParticleTag" ) ) ),
  systLabel_( pSet.getParameter<string> ( "SystLabel" ) )
  {
    inputDiPhotonName_= pSet.getParameter<std::string > ( "DiPhotonName" );
    inputDiPhotonSuffixes_= pSet.getParameter<std::vector<std::string> > ( "DiPhotonSuffixes" );
    std::vector<edm::InputTag>  diPhotonTags;
    for (auto & suffix : inputDiPhotonSuffixes_){
      systematicsLabels.push_back(suffix);
      std::string inputName = inputDiPhotonName_;
      inputName.append(suffix);
      if (!suffix.empty()) diPhotonTags.push_back(edm::InputTag(inputName));
      else  diPhotonTags.push_back(edm::InputTag(inputDiPhotonName_));
    }
    for( auto & tag : diPhotonTags ) { diPhotonTokens_.push_back( consumes<edm::View<flashgg::DiPhotonCandidate> >( tag ) ); }

    produces<vector<H4GTag>>();
    produces<vector<TagTruthBase>>();
  }
  void H4GTagProducer::produce( Event &event, const EventSetup &)
  {
    event.getByToken( diphotonToken_, diphotons );


    std::unique_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );
    edm::RefProd<vector<TagTruthBase> > rTagTruth = event.getRefBeforePut<vector<TagTruthBase> >();

    TagTruthBase truth_obj;

    if( ! event.isRealData() ) {
    Handle<View<reco::GenParticle> > genParticles;
    std::vector<edm::Ptr<reco::GenParticle> > selHiggses;
    event.getByToken( genParticleToken_, genParticles );
    reco::GenParticle::Point higgsVtx(0.,0.,0.);
    for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
        int pdgid = genParticles->ptrAt( genLoop )->pdgId();
        // if( pdgid == 25 || pdgid == 22 ) { // not so sure if this is correct for HHWWgg because of potential photons from hadronization
        if( pdgid == 25 ) {
            higgsVtx = genParticles->ptrAt( genLoop )->vertex();
            // gen_vertex_z = higgsVtx.z();
            break;
        }
    }
    for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
        edm::Ptr<reco::GenParticle> genPar = genParticles->ptrAt(genLoop);
        if (selHiggses.size()>1) break;
      if (genPar->pdgId()==25 && genPar->isHardProcess()){
          selHiggses.push_back(genPar);
      }
    }
    truth_obj.setGenPV( higgsVtx );
    truths->push_back( truth_obj );
}

    Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
    event.getByToken( diphotonToken_, diphotons );


    std::unique_ptr<vector<H4GTag> > H4Gtags( new vector<H4GTag> );
    if (diphotons->size() > 0){
      for( unsigned int diphoIndex = 0; diphoIndex < 1; diphoIndex++ ) {
        edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( diphoIndex );
        H4GTag tag_obj(dipho);
        tag_obj.setSystLabel( systLabel_);
        tag_obj.setDiPhotonIndex( diphoIndex );
        tag_obj.setCategoryNumber( 0 );
        tag_obj.includeWeights( *dipho );
        H4Gtags->push_back( tag_obj );

        if( ! event.isRealData() ) {
                  H4Gtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
                }

      }
    }
    event.put( std::move( H4Gtags ) );
    event.put(std::move(truths));
  }
}

typedef flashgg::H4GTagProducer FlashggH4GTagProducer;
DEFINE_FWK_MODULE( FlashggH4GTagProducer );
