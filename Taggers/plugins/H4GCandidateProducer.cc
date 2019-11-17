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
#include "flashgg/MicroAOD/interface/VertexSelectorBase.h"
#include "flashgg/DataFormats/interface/VertexCandidateMap.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/H4GCandidate.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "flashgg/MicroAOD/interface/CutBasedDiPhotonObjectSelector.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"

#include "TMVA/Reader.h"

#include <vector>
#include <algorithm>
#include "TGraph.h"
#include "TLorentzVector.h"
#include "TTree.h"


using namespace std;
using namespace edm;

int mcTruthVertexIndex( const std::vector<edm::Ptr<reco::GenParticle> > &genParticles ,
        const std::vector<edm::Ptr<reco::Vertex> > &vertices, const double &dzMatch )
{

    reco::Vertex::Point hardVertex( 0, 0, 0 );

    for( unsigned int genLoop = 0 ; genLoop < genParticles.size(); genLoop++ ) {

        if( fabs( genParticles[genLoop]->pdgId() ) < 10 || fabs( genParticles[genLoop]->pdgId() ) == 25 ) {
            hardVertex.SetCoordinates( genParticles[genLoop]->vx(), genParticles[genLoop]->vy(), genParticles[genLoop]->vz() );
            break;
        }
    }
    int  ivMatch = 0;
    double dzMin = 999;
    for( unsigned int iv = 0; iv < vertices.size(); iv++ ) {
        double dz = fabs( vertices[iv]->z() - hardVertex.z() );
        if( dz < dzMin ) {
            ivMatch = iv;
            dzMin   = dz;
        }
    }
    if( dzMin < dzMatch ) { return ivMatch;
    }
     return -1;
    }


namespace flashgg {

  struct Sorter {
        bool operator()( const std::pair<unsigned int, float> pair1, const std::pair<unsigned int, float> pair2 )
        {
            return ( pair1.second > pair2.second );
        };
  };

  class H4GCandidateProducer : public EDProducer
  {
  public:
    //---typedef
    typedef math::XYZTLorentzVector LorentzVector;

    //---ctors
    H4GCandidateProducer();
    H4GCandidateProducer( const ParameterSet & );

    // //Out tree elements:
    edm::Service<TFileService> fs;
    TH1F* vtxHist;
    TH1F* cutFlow;
    TTree* tree;

    // std::vector<LorentzVector> genPhoton_p4;
    std::vector<reco::Candidate::LorentzVector> genPhoton_p4;
  private:
    double genTotalWeight;

    edm::FileInPath vertexIdMVAweightfileH4G_;
    edm::FileInPath vertexProbMVAweightfileH4G_;

    TMVA::Reader *VertexIdMva_;
    TMVA::Reader *VertexProbMva_;

    void produce( Event &, const EventSetup & ) override;

    EDGetTokenT<View<Photon> > photonToken_;
    Handle<View<flashgg::Photon> > photons;

    EDGetTokenT<View<DiPhotonCandidate> > diphotonToken_;
    Handle<View<flashgg::DiPhotonCandidate> > diphotons;

    EDGetTokenT<View<reco::Vertex> > vertexToken_;
    Handle<View<reco::Vertex> > vertex;

    EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
    Handle<View<reco::GenParticle> > genParticle;

    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > genToken2_;

    EDGetTokenT<GenEventInfoProduct > genInfoToken_;
    // Handle<View<GenEventInfoProduct >> genInfo;

    EDGetTokenT<reco::BeamSpot>  beamSpotToken_;
    Handle<reco::BeamSpot>  recoBeamSpotHandle;

    EDGetTokenT<View<reco::Conversion> >  conversionToken_;
    Handle<View<reco::Conversion> > conversionHandle;

    EDGetTokenT<View<reco::Conversion> >  conversionTokenSingleLeg_;
    Handle<View<reco::Conversion> > conversionHandleSingleLeg;

    EDGetTokenT< VertexCandidateMap > vertexCandidateMapToken_;
    unique_ptr<VertexSelectorBase> vertexSelector_;

    std::vector<double> ptCuts_;
    double etaCuts_;
    std::vector<double> mvaCuts_;
    std::vector<double> hMassWindow_;

    bool useSingleLeg_;

    float logSumpt2;
    float ptAsym;
    float ptBal;
    float pullConv;
    float nConv;
    float tp_pt;
    float nVertices;
    float MVA0;
    float MVA1;
    float dZ1;
    float MVA2;
    float dZ2;
    float dZtrue;

    std::vector<std::pair<unsigned int, float> > sorter_;
    unsigned int selected_vertex_index_ = 0;
    unsigned int second_selected_vertex_index_ = 0;
    unsigned int third_selected_vertex_index_ = 0;
    float max_mva_value_ = -999.;
    float second_max_mva_value_ = -999.;
    float third_max_mva_value_ = -999.;

    //---ID selector
    ConsumesCollector cc_;
    CutBasedDiPhotonObjectSelector idSelector_;

    unsigned int isSignal;
    //----output collection
    auto_ptr<vector<H4GCandidate> > H4GColl_;


  };
  //---constructors
  H4GCandidateProducer::H4GCandidateProducer( ):
  photonToken_(),
  diphotonToken_(),
  genParticleToken_(),
  beamSpotToken_(),
  conversionToken_(),
  conversionTokenSingleLeg_(),
  cc_( consumesCollector() ),
  idSelector_( ParameterSet(), cc_ )
  {

  }
  //---standard
  H4GCandidateProducer::H4GCandidateProducer( const ParameterSet & pSet):
  photonToken_( consumes<View<Photon> >( pSet.getParameter<InputTag> ( "PhotonTag" ) ) ),
  diphotonToken_( consumes<View<DiPhotonCandidate> >( pSet.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
  // diphotonToken2_( consumes<View<DiPhotonCandidate> >( pSet.getParameter<InputTag> ( "DiPhotonTag2" ) ) ),
  vertexToken_( consumes<View<reco::Vertex> >( pSet.getParameter<InputTag> ( "VertexTag" ) ) ),
  genParticleToken_( consumes<View<reco::GenParticle> >( pSet.getParameter<InputTag> ( "GenParticleTag" ) ) ),
  genToken2_(consumes<View<pat::PackedGenParticle>> ( pSet.getParameter<InputTag> ( "GenParticleTag2" ) )),
  genInfoToken_(consumes<GenEventInfoProduct> ( pSet.getParameter<InputTag> ( "GenTag" ) )),
  beamSpotToken_( consumes<reco::BeamSpot> ( pSet.getParameter<InputTag> ( "beamSpotTag" ) ) ),
  conversionToken_( consumes <View<reco::Conversion>> ( pSet.getParameter<InputTag> ( "conversionTag" ) ) ),
  conversionTokenSingleLeg_( consumes <View<reco::Conversion>> ( pSet.getParameter<InputTag> ( "conversionTagSingleLeg" ) ) ),
  vertexCandidateMapToken_( consumes<VertexCandidateMap>( pSet.getParameter<InputTag>( "VertexCandidateMapTag" ) ) ),


  cc_( consumesCollector() ),
  idSelector_( pSet.getParameter<ParameterSet> ( "idSelection" ), cc_ )

  {
    const std::string &VertexSelectorName = pSet.getParameter<std::string>( "VertexSelectorName" );
    vertexSelector_.reset( FlashggVertexSelectorFactory::get()->create( VertexSelectorName, pSet ) );

    useSingleLeg_ = pSet.getParameter<bool>( "useSingleLeg" );

    vertexIdMVAweightfileH4G_ = pSet.getParameter<edm::FileInPath>( "vertexIdMVAweightfileH4G_2016" );
    vertexProbMVAweightfileH4G_ = pSet.getParameter<edm::FileInPath>( "vertexProbMVAweightfileH4G_2016" );

    VertexIdMva_ = new TMVA::Reader( "!Color:Silent" );
    VertexIdMva_->AddVariable( "ptAsym", &ptAsym );
    VertexIdMva_->AddVariable( "ptBal", &ptBal );
    VertexIdMva_->AddVariable( "logSumpt2", &logSumpt2 );
    VertexIdMva_->AddVariable( "pullConv", &pullConv );
    VertexIdMva_->AddVariable( "nConv", &nConv );
    VertexIdMva_->BookMVA( "BDT", vertexIdMVAweightfileH4G_.fullPath() );

    VertexProbMva_ = new TMVA::Reader( "!Color:Silent" );
    VertexProbMva_->AddVariable( "tp_pt", &tp_pt);
    VertexProbMva_->AddVariable( "n_vertices", &nVertices );
    VertexProbMva_->AddVariable( "MVA0", &MVA0 );
    VertexProbMva_->AddVariable( "MVA1", &MVA1 );
    VertexProbMva_->AddVariable( "dZ1", &dZ1 );
    VertexProbMva_->AddVariable( "MVA2", &MVA2 );
    VertexProbMva_->AddVariable( "dZ2", &dZ2 );
    VertexProbMva_->AddVariable( "nConv", &nConv );
    VertexProbMva_->BookMVA( "BDT", vertexProbMVAweightfileH4G_.fullPath() );

    vtxHist = fs->make<TH1F> ("vtxHist","vtx Hist",10,0,10);

    ptCuts_ = pSet.getParameter<std::vector<double>>("ptCuts");
    etaCuts_ = pSet.getParameter<double>("etaCuts");
    mvaCuts_ = pSet.getParameter<std::vector<double>>("mvaCuts");
    hMassWindow_ = pSet.getParameter<std::vector<double>>("hMassWindow");

    cutFlow = fs->make<TH1F> ("cutFlow","Cut FLow",10,0,10);

    tree = new TTree("genTree"," gen level H4G Tree");
    tree->Branch("genPhoton_p4",&genPhoton_p4);

    unsigned int def_isSignal = 0;
    isSignal = pSet.getUntrackedParameter<unsigned int>("isSignal", def_isSignal);


    produces<vector<H4GCandidate> > ();
  }
  void H4GCandidateProducer::produce( Event &event, const EventSetup & )
  {
    event.getByToken( photonToken_, photons );
    event.getByToken( diphotonToken_, diphotons );
    event.getByToken( vertexToken_, vertex );
    event.getByToken( genParticleToken_, genParticle );
    event.getByToken( beamSpotToken_, recoBeamSpotHandle );
    event.getByToken( conversionToken_, conversionHandle );
    event.getByToken( conversionTokenSingleLeg_, conversionHandleSingleLeg );

    Handle<VertexCandidateMap> vertexCandidateMap;
    event.getByToken( vertexCandidateMapToken_, vertexCandidateMap );

    int hgg_index = -999;

    genPhoton_p4.clear();

    math::XYZPoint BSPoint;
    if( recoBeamSpotHandle.isValid() ) {
      BSPoint = recoBeamSpotHandle->position();
    }

    //---output collection
    std::unique_ptr<vector<H4GCandidate> > H4GColl_( new vector<H4GCandidate> );
    std::vector <edm::Ptr<reco::Vertex> > Vertices; // Collection of vertices
    std::vector <edm::Ptr<reco::Vertex> > slim_Vertices;
    reco::GenParticle::Point genVertex;

    Handle<View<reco::Vertex> > primaryVertices;
    event.getByToken( vertexToken_, primaryVertices );

    int trueVtxIndexI = -999;
    vector<int>	pvVecNoTrue;
    int irand = -999;
    int randVtxIndexI = -999;
    if( ! event.isRealData() )
    {
      Handle<View<reco::GenParticle> > genParticles;
      event.getByToken( genParticleToken_, genParticles );
      // if (isSignal) {
        trueVtxIndexI = mcTruthVertexIndex( genParticles->ptrs(), primaryVertices->ptrs(), 0.1);
      // }
      for( unsigned int i = 0 ; i < primaryVertices->size() ; i++ ) {
           if( i != (unsigned int)trueVtxIndexI ) { pvVecNoTrue.push_back( i ); }
      }
      if( pvVecNoTrue.size() > 1 ) { irand = rand() % pvVecNoTrue.size(); }
      if( irand != -999 ) { randVtxIndexI = pvVecNoTrue[irand]; }

      for( auto &part : *genParticle ) {
        if( part.pdgId() != 2212 || part.vertex().z() != 0. )
        {
          genVertex = part.vertex();
        }
      }
      for( auto &part : *genParticle )
      {
        if (part.pdgId() == 25 || part.pdgId() == 54)
      {
       genPhoton_p4.push_back(part.daughter(0)->p4());
       genPhoton_p4.push_back(part.daughter(1)->p4());
       }
     }

    }


    edm::Ptr<reco::Vertex> vertex_diphoton;
    edm::Ptr<reco::Vertex> vertex_bdt;
    //---at least one diphoton should pass the low mass hgg pre-selection
    bool atLeastOneDiphoPass = false;
    std::vector<const flashgg::Photon*> phosTemp;
    for( unsigned int dpIndex = 0; dpIndex < diphotons->size(); dpIndex++ )
    {
      edm::Ptr<flashgg::DiPhotonCandidate> thisDPPtr = diphotons->ptrAt( dpIndex );
      vertex_diphoton = diphotons->ptrAt( dpIndex )->vtx();
      flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(thisDPPtr.get());
      atLeastOneDiphoPass |= idSelector_(*thisDPPointer, event);
    }

    int n_photons = -999;
    n_photons = photons->size();
    Handle<GenEventInfoProduct> genInfo;
    if( ! event.isRealData() )
    {
     event.getByToken(genInfoToken_, genInfo);
     genTotalWeight = genInfo->weight();
   } else {
      genTotalWeight = 1;
   }

   cutFlow->Fill(0.0,genTotalWeight); // all events
    if (n_photons == 4 )
    {
     cutFlow->Fill(1.0,genTotalWeight); // events w/ 4 photons
    }

    std::vector<edm::Ptr<flashgg::Photon>> phoPtrVector;
    std::vector< flashgg::Photon> diphoPhotons;
    std::vector<flashgg::Photon*> diphoPhotonsVector;
    std::vector <edm::Ptr<flashgg::DiPhotonCandidate>> diPhoPtrs;

    int n_pho = 0;

    if (atLeastOneDiphoPass)
    {
      for( unsigned int dpIndex = 0; dpIndex < diphotons->size(); dpIndex++ )
      {
        edm::Ptr<flashgg::DiPhotonCandidate> thisDiPho = diphotons->ptrAt( dpIndex );
        flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(thisDiPho.get());
        thisDPPointer->makePhotonsPersistent();
        diPhoPtrs.push_back(thisDiPho);
        auto pho1 = thisDPPointer->getLeadingPhoton();
        auto pho2 = thisDPPointer->getSubLeadingPhoton();

         if (diphoPhotons.size() == 0 ){
          diphoPhotons.push_back(pho1);
          diphoPhotons.push_back(pho2);
          n_pho +=2;
          continue;
        }
        else {
          float minDR1 = 999, minDR2 = 999;
          for (size_t p=0; p < diphoPhotons.size(); p++){
            float deltar1 = sqrt(pow(diphoPhotons[p].superCluster()->eta()-pho1.superCluster()->eta(),2)+pow(diphoPhotons[p].superCluster()->phi()-pho1.superCluster()->phi(),2));
            float deltar2 = sqrt(pow(diphoPhotons[p].superCluster()->eta()-pho2.superCluster()->eta(),2)+pow(diphoPhotons[p].superCluster()->phi()-pho2.superCluster()->phi(),2));
            if (deltar1 < minDR1) minDR1 = deltar1;
            if (deltar2 < minDR2) minDR2 = deltar2;
          }
          if (minDR1 > 0.00001)
          {
            n_pho++;
            diphoPhotons.push_back(pho1);
          }
          if (minDR2 > 0.00001)
          {
            n_pho++;
            diphoPhotons.push_back(pho2);
          }
        }
      }

      std::sort(diphoPhotons.begin(), diphoPhotons.end(), [](const flashgg::Photon a, const flashgg::Photon b) {return a.pt() > b.pt(); });

      for( int phoIndex = 0; phoIndex < n_photons; phoIndex++ )
      {
        edm::Ptr<flashgg::Photon> pho = photons->ptrAt(phoIndex);
        phoPtrVector.push_back(pho);
      }

      std::sort(phoPtrVector.begin(), phoPtrVector.end(), [](const edm::Ptr<flashgg::Photon> a, const edm::Ptr<flashgg::Photon> b) {return a->pt() > b->pt(); });

      for( int v = 0; v < (int) vertex->size(); v++ )
      {

        edm::Ptr<reco::Vertex> vtx = vertex->ptrAt( v );
        Vertices.push_back(vtx);
        if (vertex_diphoton->x() - vtx->x() == 0  && vertex_diphoton->y() - vtx->y() == 0 && vertex_diphoton->z() - vtx->z() == 0 ){
          hgg_index =  v;
        }
        if (fabs(genVertex.z() - vertex_diphoton->z()) > 1 ){
          if (fabs(genVertex.z() - vtx->z()) < 1)
          {
            slim_Vertices.push_back(vtx);
          }
          else{
            slim_Vertices.push_back(vertex->ptrAt( 0 ));
          }
        }
        else {
          slim_Vertices.push_back(vertex_diphoton);
        }
      }
      int trueVtxIndex = trueVtxIndexI;
      int randVtxIndex = randVtxIndexI;

      std::vector<std::vector<float>> vtxVar;
      if (phoPtrVector.size() == 2)
      {
        vtxVar = vertexSelector_->select_h2g(phoPtrVector[0],phoPtrVector[1], Vertices, *vertexCandidateMap,conversionHandle->ptrs(), conversionHandleSingleLeg->ptrs(), BSPoint, useSingleLeg_   );
      }
      if (phoPtrVector.size() == 3)
      {
        vtxVar = vertexSelector_->select_h3g(phoPtrVector[0],phoPtrVector[1],phoPtrVector[2], Vertices, *vertexCandidateMap,conversionHandle->ptrs(), conversionHandleSingleLeg->ptrs(), BSPoint, useSingleLeg_   );
      }
      if (phoPtrVector.size() > 3)
      {
        vtxVar = vertexSelector_->select_h4g(phoPtrVector[0],phoPtrVector[1],phoPtrVector[2],phoPtrVector[3], Vertices, *vertexCandidateMap,conversionHandle->ptrs(), conversionHandleSingleLeg->ptrs(), BSPoint, useSingleLeg_   );
      }

      sorter_.clear();
      selected_vertex_index_ = 0;
      second_selected_vertex_index_ = 0;
      third_selected_vertex_index_ = 0;
      max_mva_value_ = -999.;
      second_max_mva_value_ = -999.;
      third_max_mva_value_ = -999.;
      for( int vtx = 0; vtx < (int) vertex->size(); vtx++ )
      {
           logSumpt2 = vtxVar[0][vtx];
           ptAsym = vtxVar[1][vtx];
           ptBal = vtxVar[2][vtx];
           pullConv = vtxVar[3][vtx];
           nConv = vtxVar[4][vtx];

           float mva_value_ = VertexIdMva_->EvaluateMVA( "BDT" );

           std::pair<unsigned int, float>pairToSort = std::make_pair( (unsigned int)vtx, mva_value_ );
           sorter_.push_back( pairToSort );

           if( mva_value_ > max_mva_value_ ) {
                max_mva_value_ = mva_value_;
                selected_vertex_index_ = vtx;
           }
      }

      std::sort( sorter_.begin(), sorter_.end(), Sorter() );

      if( sorter_.size() > 1 ) {
          second_max_mva_value_ = sorter_[1].second;
          second_selected_vertex_index_ = sorter_[1].first;
      }

      if( sorter_.size() > 2 ) {
          third_max_mva_value_ = sorter_[2].second;
          third_selected_vertex_index_ = sorter_[2].first;
      }

      vertex_bdt = vertex->ptrAt( selected_vertex_index_ );

      MVA0      = max_mva_value_;
      MVA1      = second_max_mva_value_;
      MVA2      = third_max_mva_value_;
      dZ1       = vertex->at(selected_vertex_index_).position().z() - vertex->at(second_selected_vertex_index_).position().z();
      dZ2       = vertex->at(selected_vertex_index_).position().z() - vertex->at(third_selected_vertex_index_).position().z();
      dZtrue    = vertex->at(selected_vertex_index_).position().z() - genVertex.z();
      nVertices = (float) vertex->size();
      nConv = (float)vtxVar[4][selected_vertex_index_];

      if (n_photons == 4)
      {
        cutFlow->Fill(2.0,genTotalWeight); // events passing the online cuts
        if ((photons->ptrAt(0)->p4()+photons->ptrAt(1)->p4()+photons->ptrAt(2)->p4()+photons->ptrAt(3)->p4()).mass() > hMassWindow_[0] && (photons->ptrAt(0)->p4()+photons->ptrAt(1)->p4()+photons->ptrAt(2)->p4()+photons->ptrAt(3)->p4()).mass() < hMassWindow_[1])
        {
          cutFlow->Fill(3.0,genTotalWeight); // events lying within Higgs mass window
          if (photons->ptrAt(0)->pt() > ptCuts_[0] && photons->ptrAt(1)->pt() > ptCuts_[1] && photons->ptrAt(2)->pt() > ptCuts_[2] && photons->ptrAt(3)->pt() > ptCuts_[3])
          {
            if (abs(photons->ptrAt(0)->eta()) < etaCuts_ && abs(photons->ptrAt(1)->eta()) < etaCuts_ && abs(photons->ptrAt(2)->eta()) < etaCuts_ && abs(photons->ptrAt(3)->eta()) < etaCuts_)
            {
            cutFlow->Fill(4.0,genTotalWeight); // events passing kinematic selection
            if (photons->ptrAt(0)->phoIdMvaDWrtVtx(vertex_bdt) > mvaCuts_[0] && photons->ptrAt(1)->phoIdMvaDWrtVtx(vertex_bdt) > mvaCuts_[1] && photons->ptrAt(2)->phoIdMvaDWrtVtx(vertex_bdt) > mvaCuts_[2] && photons->ptrAt(3)->phoIdMvaDWrtVtx(vertex_bdt) > mvaCuts_[3] )
            {
              cutFlow->Fill(5.0,genTotalWeight); // events passing photon ID MVA cut
            }
          }
        }
      }
    }
      H4GCandidate h4g(Vertices, slim_Vertices, vertex_diphoton, vertex_bdt, genVertex, BSPoint, vtxVar, MVA0, MVA1, MVA2, dZ1, dZ2, dZtrue, hgg_index, trueVtxIndex, randVtxIndex, selected_vertex_index_, tp_pt, nVertices, nConv, VertexProbMva_, genTotalWeight,diphoPhotons,diPhoPtrs,genPhoton_p4);
      H4GColl_->push_back(h4g);
    }
    event.put( std::move(H4GColl_) );
    tree->Fill();
  }
}

typedef flashgg::H4GCandidateProducer FlashggH4GCandidateProducer;
DEFINE_FWK_MODULE( FlashggH4GCandidateProducer );
