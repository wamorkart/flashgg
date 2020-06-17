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
#include "flashgg/MicroAOD/interface/VertexSelectorBase.h"
#include "flashgg/DataFormats/interface/VertexCandidateMap.h"
#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/H4GTag.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "flashgg/MicroAOD/interface/CutBasedDiPhotonObjectSelector.h"

#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
#include "flashgg/DataFormats/interface/WHLeptonicTag.h"

#include "flashgg/Taggers/interface/LeptonSelection.h"
#include "TTree.h"
#include "TMVA/Reader.h"
#include <vector>
#include <algorithm>
#include "TGraph.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;

int mcTruthVertexIndex_h4g( const std::vector<edm::Ptr<reco::GenParticle> > &genParticles ,
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
    if( dzMin < dzMatch ) { return ivMatch;}

    return -1;
  }

  namespace flashgg {

    struct Sorter {
      bool operator()( const std::pair<unsigned int, float> pair1, const std::pair<unsigned int, float> pair2 )
      {
        return ( pair1.second > pair2.second );
      };
    };
    struct Sorter_pairs {
      bool operator()( const std::pair<std::vector<int>,float> pair1, std::pair<std::vector<int>,float> pair2 )
      {
        return ( pair1.second > pair2.second );
      };
    };

    class H4GTagProducer : public EDProducer
    {
    public:
      //---typedef
      typedef math::XYZTLorentzVector LorentzVector;

      //---ctors
      H4GTagProducer( const ParameterSet & );



    private:
      void produce( Event &, const EventSetup & ) override;
      std::vector<edm::EDGetTokenT<edm::View<DiPhotonCandidate> > > diPhotonTokens_;
      std::string inputDiPhotonName_;

      EDGetTokenT<View<Photon> > photonToken_;
      Handle<View<flashgg::Photon> > photons;

      EDGetTokenT<View<DiPhotonCandidate> > diphotonToken_;
      Handle<View<flashgg::DiPhotonCandidate> > diphotons;

      EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
      Handle<View<reco::GenParticle> > genParticle;

      EDGetTokenT<View<reco::Vertex> > vertexToken_;
      Handle<View<reco::Vertex> > vertex;


      //---ID selector
      ConsumesCollector cc_;
      CutBasedDiPhotonObjectSelector idSelector_;


      string systLabel_;
      std::vector< std::string > systematicsLabels;
      std::vector<std::string> inputDiPhotonSuffixes_;

      EDGetTokenT<reco::BeamSpot>  beamSpotToken_;
      Handle<reco::BeamSpot>  recoBeamSpotHandle;

      EDGetTokenT<View<reco::Conversion> >  conversionToken_;
      Handle<View<reco::Conversion> > conversionHandle;

      EDGetTokenT<View<reco::Conversion> >  conversionTokenSingleLeg_;
      Handle<View<reco::Conversion> > conversionHandleSingleLeg;

      EDGetTokenT< VertexCandidateMap > vertexCandidateMapToken_;
      unique_ptr<VertexSelectorBase> vertexSelector_;

      edm::InputTag genInfo_;
      edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;

      edm::FileInPath vertexIdMVAweightfileH4G_;
      // edm::FileInPath vertexProbMVAweightfileH4G_;
      TMVA::Reader *VertexIdMva_;
      // TMVA::Reader *VertexProbMva_;

      bool useSingleLeg_;
      float logSumpt2;
      float ptAsym;
      float ptBal;
      float pullConv;
      float nConv;
      // float tp_pt;
      // float nVertices;
      // float MVA0;
      // float MVA1;
      // float dZ1;
      // float MVA2;
      // float dZ2;
      // float dZtrue;

      bool doH4GVertex_;

      std::vector<std::pair<unsigned int, float> > sorter_;
      unsigned int selected_vertex_index_ = 0;
      unsigned int second_selected_vertex_index_ = 0;
      unsigned int third_selected_vertex_index_ = 0;
      float max_mva_value_ = -999.;
      float second_max_mva_value_ = -999.;
      float third_max_mva_value_ = -999.;

    };


    //---standard
    H4GTagProducer::H4GTagProducer( const ParameterSet & pSet):
    photonToken_( consumes<View<Photon> >( pSet.getParameter<InputTag> ( "PhotonTag" ) ) ),
    diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( pSet.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    genParticleToken_( consumes<View<reco::GenParticle> >( pSet.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( pSet.getParameter<InputTag> ( "VertexTag" ) ) ),
    cc_( consumesCollector() ),
    idSelector_( pSet.getParameter<ParameterSet> ( "idSelection" ), cc_ ),
    systLabel_( pSet.getParameter<string> ( "SystLabel" ) ),
    beamSpotToken_( consumes<reco::BeamSpot> ( pSet.getParameter<InputTag> ( "beamSpotTag" ) ) ),
    conversionToken_( consumes <View<reco::Conversion>> ( pSet.getParameter<InputTag> ( "conversionTag" ) ) ),
    conversionTokenSingleLeg_( consumes <View<reco::Conversion>> ( pSet.getParameter<InputTag> ( "conversionTagSingleLeg" ) ) ),
    vertexCandidateMapToken_( consumes<VertexCandidateMap>( pSet.getParameter<InputTag>( "VertexCandidateMapTag" ) ) )


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



      genInfo_ = pSet.getUntrackedParameter<edm::InputTag>( "genInfo", edm::InputTag("generator") );
      genInfoToken_ = consumes<GenEventInfoProduct>( genInfo_ );

      // cout << "**************************** in H4GTagProducer.cc **********************************************" << endl;



      
      useSingleLeg_ = pSet.getParameter<bool>( "useSingleLeg" );
      doH4GVertex_ = pSet.getParameter<bool>("doH4GVertex");
      vertexIdMVAweightfileH4G_ = pSet.getParameter<edm::FileInPath>( "vertexIdMVAweightfileH4G" );
      // vertexProbMVAweightfileH4G_ = pSet.getParameter<edm::FileInPath>( "vertexProbMVAweightfileH4G" );

      if (doH4GVertex_)
      {
      const std::string &VertexSelectorName = pSet.getParameter<std::string>( "VertexSelectorName" );
      vertexSelector_.reset( FlashggVertexSelectorFactory::get()->create( VertexSelectorName, pSet ) );

      VertexIdMva_ = new TMVA::Reader( "!Color:Silent" );
      VertexIdMva_->AddVariable( "ptAsym", &ptAsym );
      VertexIdMva_->AddVariable( "ptBal", &ptBal );
      VertexIdMva_->AddVariable( "logSumpt2", &logSumpt2 );
      VertexIdMva_->AddVariable( "pullConv", &pullConv );
      VertexIdMva_->AddVariable( "nConv", &nConv );
      VertexIdMva_->BookMVA( "BDT", vertexIdMVAweightfileH4G_.fullPath() );
      }
      // VertexProbMva_ = new TMVA::Reader( "!Color:Silent" );
      // VertexProbMva_->AddVariable( "tp_pt", &tp_pt);
      // VertexProbMva_->AddVariable( "n_vertices", &nVertices );
      // VertexProbMva_->AddVariable( "MVA0", &MVA0 );
      // VertexProbMva_->AddVariable( "MVA1", &MVA1 );
      // VertexProbMva_->AddVariable( "dZ1", &dZ1 );
      // VertexProbMva_->AddVariable( "MVA2", &MVA2 );
      // VertexProbMva_->AddVariable( "dZ2", &dZ2 );
      // VertexProbMva_->AddVariable( "nConv", &nConv );
      // VertexProbMva_->BookMVA( "BDT", vertexProbMVAweightfileH4G_.fullPath() );


      produces<vector<H4GTag>>();

      produces<vector<TagTruthBase>>();
    }


    void H4GTagProducer::produce( Event &event, const EventSetup & )
    {

      // cout << "[H4GTagProducer.cc] - Beginning of H4GTagProducer::produce" << endl;

      // update global variables
      // globalVariablesComputer_.update(event);

      // Get particle objects
      event.getByToken( photonToken_, photons );
      event.getByToken( diphotonToken_, diphotons );
      event.getByToken( vertexToken_, vertex );
      event.getByToken( genParticleToken_, genParticle );
      event.getByToken( beamSpotToken_, recoBeamSpotHandle );
      if (doH4GVertex_)
      {
      event.getByToken( conversionToken_, conversionHandle );
      event.getByToken( conversionTokenSingleLeg_, conversionHandleSingleLeg );
      }

      Handle<VertexCandidateMap> vertexCandidateMap;
      event.getByToken( vertexCandidateMapToken_, vertexCandidateMap );

      Handle<View<reco::Vertex> > primaryVertices;
      event.getByToken( vertexToken_, primaryVertices );

      std::unique_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );
      edm::RefProd<vector<TagTruthBase> > rTagTruth = event.getRefBeforePut<vector<TagTruthBase> >();

      //--------vertex output-------//
      vector <edm::Ptr<reco::Vertex> > Vertices; // Collection of vertices


      //---------vertex output end-----//
      //-----------------------------------------------------------------------------------------------------------
      math::XYZPoint BSPoint;
      if( recoBeamSpotHandle.isValid() ) {
        BSPoint = recoBeamSpotHandle->position();
      }

      // MC truth
      TagTruthBase truth_obj;
      // double genMhh=0.;
      int trueVtxIndexI = -999;
      vector<int>	pvVecNoTrue;
      reco::GenParticle::Point higgsVtx(0.,0.,0.);
      reco::Vertex::Point hardVertex( 0, 0, 0 );
      reco::GenParticle::Point genVertex;
      int irand = -999;
      int randVtxIndexI = -999;
      if( ! event.isRealData() ) {

        Handle<View<reco::GenParticle> > genParticles;
        event.getByToken( genParticleToken_, genParticles );
        std::vector<edm::Ptr<reco::GenParticle> > selHiggses;
        trueVtxIndexI = mcTruthVertexIndex_h4g( genParticles->ptrs(), primaryVertices->ptrs(), 0.1);
        for( unsigned int i = 0 ; i < primaryVertices->size() ; i++ ) {
          if( i != (unsigned int)trueVtxIndexI ) { pvVecNoTrue.push_back( i ); }
        }
        if( pvVecNoTrue.size() > 1 ) { irand = rand() % pvVecNoTrue.size(); }
        if( irand != -999 ) { randVtxIndexI = pvVecNoTrue[irand]; }

        for( unsigned int genLoop = 0 ; genLoop < genParticles->ptrs().size(); genLoop++ ) {

        if( fabs( genParticles->ptrs()[genLoop]->pdgId() ) < 10 || fabs( genParticles->ptrs()[genLoop]->pdgId() ) == 25 ) {
            hardVertex.SetCoordinates( genParticles->ptrs()[genLoop]->vx(), genParticles->ptrs()[genLoop]->vy(), genParticles->ptrs()[genLoop]->vz() );
            break;
        }
    }

        for( auto &part : *genParticle ) {
        if( part.pdgId() != 2212 || part.vertex().z() != 0. )
        {
          genVertex = part.vertex();
          //cout << part.pdgId()  <<  "   "  << genVertex.z() << endl;
        }
      }


        // reco::GenParticle::Point higgsVtx(0.,0.,0.);
        for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
          int pdgid = genParticles->ptrAt( genLoop )->pdgId();
          if( fabs( pdgid ) < 10 || fabs( pdgid ) == 25 ) {
            if (pdgid == 25) {  higgsVtx = genParticles->ptrAt( genLoop )->vertex(); }

            // hardVertex.SetCoordinates( genParticles[genLoop]->vx(), genParticles[genLoop]->vy(), genParticles[genLoop]->vz() );
            //cout << "Higgs Vtx " << higgsVtx.z() << endl;
            // gen_vertex_z = higgsVtx.z();
            break;
          }
        }

        truth_obj.setGenPV( higgsVtx );
        truths->push_back( truth_obj );
      }

      Handle<View<flashgg::DiPhotonCandidate> > diPhotons;

      //---
      event.getByToken( diphotonToken_, diphotons ); // without looping over diphoton systematics

      //---
      std::unique_ptr<vector<H4GTag> > H4Gtags( new vector<H4GTag> );
      //
      // edm::Ptr<reco::Vertex> vertex_diphoton;
      edm::Ptr<reco::Vertex> vertex_chosen;
      edm::Ptr<reco::Vertex> vertex_zero;
      bool atLeastOneDiphoPass = false;
      for( unsigned int dpIndex = 0; dpIndex < diphotons->size(); dpIndex++ )
      {
        edm::Ptr<flashgg::DiPhotonCandidate> thisDPPtr = diphotons->ptrAt( dpIndex );
        // vertex_diphoton = diphotons->ptrAt( dpIndex )->vtx();
        flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(thisDPPtr.get());
        atLeastOneDiphoPass |= idSelector_(*thisDPPointer, event);
      }

      // cout << "atLeastOneDiphoPass: " << atLeastOneDiphoPass << endl;

      //if (photons->size()>3){ // without systematics (?)
      std::vector<edm::Ptr<flashgg::Photon>> phoPtrVector;
      // if (atLeastOneDiphoPass && photons->size() > 3){
      if (atLeastOneDiphoPass){
        edm::Ptr<flashgg::DiPhotonCandidate> dipho;

        vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoVec;
        int n_pho = 0;

        vector<flashgg::Photon> phoVector;
        // std::vector<edm::Ptr<flashgg::Photon>> phoPtrVector;

        // cout << "# of diphotons " << diphotons->size() << endl;

        for( unsigned int diphoIndex = 0; diphoIndex < diphotons->size(); diphoIndex++ ) {

          dipho = diphotons->ptrAt( diphoIndex ); // without systematic look (?)
          diphoVec.push_back(dipho);
          // create a vector of unique photons from all diphotons
          edm::Ptr<flashgg::DiPhotonCandidate> thisDiPho = diphotons->ptrAt( diphoIndex );
          flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(thisDiPho.get());
          thisDPPointer->makePhotonsPersistent();
          auto pho1 = thisDPPointer->getLeadingPhoton();
          auto pho2 = thisDPPointer->getSubLeadingPhoton();

          if (phoVector.size() == 0)
          {
            phoVector.push_back(pho1);
            phoVector.push_back(pho2);
            n_pho +=2;
            continue;
          }
          else{
            float minDR1 = 999, minDR2 = 999;
            for (size_t p=0; p < phoVector.size(); p++){
              float deltar1 = sqrt(pow(phoVector[p].superCluster()->eta()-pho1.superCluster()->eta(),2)+pow(phoVector[p].superCluster()->phi()-pho1.superCluster()->phi(),2));
              float deltar2 = sqrt(pow(phoVector[p].superCluster()->eta()-pho2.superCluster()->eta(),2)+pow(phoVector[p].superCluster()->phi()-pho2.superCluster()->phi(),2));
              if (deltar1 < minDR1) minDR1 = deltar1;
              if (deltar2 < minDR2) minDR2 = deltar2;
            }
            if (minDR1 > 0.00001)
            {
              n_pho++;
              phoVector.push_back(pho1);
            }
            if (minDR2 > 0.00001)
            {
              n_pho++;
              phoVector.push_back(pho2);
            }
          }
        }
        // sort the photons in Pt
        std::sort(phoVector.begin(), phoVector.end(), [](const flashgg::Photon a, const flashgg::Photon b) {return a.pt() > b.pt(); });

        //-- prepare a vector of ptr to photons, to be used for vertex selection
        for( int phoIndex = 0; phoIndex < (int) photons->size(); phoIndex++ )
        {
          edm::Ptr<flashgg::Photon> pho = photons->ptrAt(phoIndex);
          phoPtrVector.push_back(pho);
        }
        std::sort(phoPtrVector.begin(), phoPtrVector.end(), [](const edm::Ptr<flashgg::Photon> a, const edm::Ptr<flashgg::Photon> b) {return a->pt() > b->pt(); });

        for( int v = 0; v < (int) vertex->size(); v++ )
        {
          edm::Ptr<reco::Vertex> vtx = vertex->ptrAt( v );
          Vertices.push_back(vtx);
        }
        int trueVtxIndex = trueVtxIndexI;
        int randVtxIndex = randVtxIndexI;


        std::vector<std::vector<float>> vtxVar;
        sorter_.clear();
        selected_vertex_index_ = 0;
        second_selected_vertex_index_ = 0;
        third_selected_vertex_index_ = 0;
        max_mva_value_ = -999.;
        second_max_mva_value_ = -999.;
        third_max_mva_value_ = -999.;

        if (doH4GVertex_)
        {
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

          // sorter_.clear();
          // selected_vertex_index_ = 0;
          // second_selected_vertex_index_ = 0;
          // third_selected_vertex_index_ = 0;
          // max_mva_value_ = -999.;
          // second_max_mva_value_ = -999.;
          // third_max_mva_value_ = -999.;

          for( int vtx = 0; vtx < (int) vertex->size(); vtx++ )
          {
            logSumpt2 = vtxVar[0][vtx];
            ptAsym = vtxVar[1][vtx];
            ptBal = vtxVar[2][vtx];
            pullConv = vtxVar[3][vtx];
            nConv = vtxVar[4][vtx];

            float mva_value_bdt_ = VertexIdMva_->EvaluateMVA( "BDT" );

            // cout << "mva value " << mva_value_bdt_ << endl;

            std::pair<unsigned int, float>pairToSort = std::make_pair( (unsigned int)vtx, mva_value_bdt_ );
            sorter_.push_back( pairToSort );

            if( mva_value_bdt_ > max_mva_value_ ) {
              max_mva_value_ = mva_value_bdt_;
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

        }

        vertex_chosen = vertex->ptrAt( selected_vertex_index_ );
        vertex_zero = vertex->ptrAt( 0 );
        //cout << "# of photons: " << phoVector.size() << endl;
        //cout << "selected vertex index " << selected_vertex_index_ << endl;

        float dZ_bdtVtx = -999;

        if( ! event.isRealData() ){ dZ_bdtVtx =   hardVertex.z() - vertex_chosen->z(); }

        //cout << dZ_bdtVtx << endl;

        // MVA0      = max_mva_value_;
        // MVA1      = second_max_mva_value_;
        // MVA2      = third_max_mva_value_;
        // dZ1       = vertex->at(selected_vertex_index_).position().z() - vertex->at(second_selected_vertex_index_).position().z();
        // dZ2       = vertex->at(selected_vertex_index_).position().z() - vertex->at(third_selected_vertex_index_).position().z();
        // dZtrue    = vertex->at(selected_vertex_index_).position().z() - genVertex.z();
        // nVertices = (float) vertex->size();
        // nConv = (float)vtxVar[4][selected_vertex_index_];

        float vtx_X = Vertices[selected_vertex_index_]->x();
        float vtx_Y = Vertices[selected_vertex_index_]->y();
        float vtx_Z = Vertices[selected_vertex_index_]->z();
        math::XYZVector vtx_Pos( vtx_X, vtx_Y, vtx_Z );
        // cout << "vtx X " << vtx_X << "vtx Y" << vtx_Y << "vtx Z " << vtx_Z << "selected_vertex_index_: " << selected_vertex_index_ << endl;

        vector<flashgg::Photon> phoP4Corrected_dp;
        if (phoVector.size() > 0)
        {
          for (int dp = 0; dp < (int) phoVector.size(); dp++)
          {
            float sc_X_dp = phoVector[dp].superCluster()->x();
            float sc_Y_dp = phoVector[dp].superCluster()->y();
            float sc_Z_dp = phoVector[dp].superCluster()->z();
            math::XYZVector sc_Pos_dp( sc_X_dp, sc_Y_dp, sc_Z_dp );
            math::XYZVector direction_dp = sc_Pos_dp - vtx_Pos;
            math::XYZVector pho_dp = ( direction_dp.Unit() ) * ( phoVector[dp].energy() );
            math::XYZTLorentzVector corrected_p4_dp( pho_dp.x(), pho_dp.y(), pho_dp.z(), phoVector[dp].energy() );
            phoVector[dp].setP4(corrected_p4_dp);
            phoP4Corrected_dp.push_back(phoVector[dp]);
          }
        }
        sorter_.clear();
        std::sort(phoP4Corrected_dp.begin(), phoP4Corrected_dp.end(), [](const flashgg::Photon a, const flashgg::Photon  b) {return a.pt() > b.pt(); });

        // for (int p=0; p < (int) phoP4Corrected_dp.size(); p++)
        // {
        //   cout << phoP4Corrected_dp[p].pt() <<endl;
        // }

        if (phoP4Corrected_dp.size() > 3)
        {
                   //cout  << phoP4Corrected_dp[0].phoIdMvaDWrtVtx(Vertices[selected_vertex_index_]) << "  " << phoP4Corrected_dp[1].phoIdMvaDWrtVtx(Vertices[selected_vertex_index_]) << "  " << phoP4Corrected_dp[2].phoIdMvaDWrtVtx(Vertices[selected_vertex_index_]) << "  " << phoP4Corrected_dp[3].phoIdMvaDWrtVtx(Vertices[selected_vertex_index_]) << endl;
                   //cout <<"HERE Vtx Zero " << phoP4Corrected_dp[0].phoIdMvaDWrtVtx(vertex_zero) << "  " << phoP4Corrected_dp[1].phoIdMvaDWrtVtx(vertex_zero) << "  " << phoP4Corrected_dp[2].phoIdMvaDWrtVtx(vertex_zero) << "  " << phoP4Corrected_dp[3].phoIdMvaDWrtVtx(vertex_zero) << endl;

                   //cout  << phoP4Corrected_dp[0].pt() << "  " << phoP4Corrected_dp[1].pt() << "  " << phoP4Corrected_dp[2].pt() << "  " << phoP4Corrected_dp[3].pt() << endl;
                   //cout  << phoP4Corrected_dp[0].eta() << "  " << phoP4Corrected_dp[1].eta() << "  " << phoP4Corrected_dp[2].eta() << "  " << phoP4Corrected_dp[3].eta() << endl;

          // int catnum = 0;
          //cout << "4 photon vtx "<< selected_vertex_index_ << endl;
          H4GTag tag_obj (dipho, phoP4Corrected_dp[0], phoP4Corrected_dp[1], phoP4Corrected_dp[2], phoP4Corrected_dp[3], vertex_chosen, dZ_bdtVtx);
          // cout << "systlabel: " << systLabel_ << endl;
          tag_obj.setSystLabel( systLabel_);
          // tag_obj.setDiPhotonIndex( diphoIndex );
          // tag_obj.setMVA( -0.9 );
          tag_obj.setCategoryNumber( 0 );
          for (int i1=0; i1 < (int) diphoVec.size() ; i1++)
          {
            auto diphotemp = diphoVec[i1];
            tag_obj.includeWeights(* diphotemp);
          }

          H4Gtags->push_back( tag_obj );
          if( ! event.isRealData() ) {
            H4Gtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
          }

        }

        if (phoP4Corrected_dp.size() == 3)
        {
          // int catnum = 1;
          //cout << "3 photon vtx "<< selected_vertex_index_ << endl;
          H4GTag tag_obj (dipho, phoP4Corrected_dp[0], phoP4Corrected_dp[1], phoP4Corrected_dp[2],vertex_chosen, dZ_bdtVtx);
          // cout << "systlabel: " << systLabel_ << endl;
          tag_obj.setSystLabel( systLabel_);
          // tag_obj.setDiPhotonIndex( diphoIndex );
          // tag_obj.setMVA( -0.9 );
          tag_obj.setCategoryNumber( 1 );
          for (int i1=0; i1 < (int) diphoVec.size() ; i1++)
          {
            auto diphotemp = diphoVec[i1];
            tag_obj.includeWeights(* diphotemp);
          }

          H4Gtags->push_back( tag_obj );
          if( ! event.isRealData() ) {
            H4Gtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
          }

        }

        if (phoP4Corrected_dp.size() == 2)
        {
          // int catnum = 2;

          //cout << "2 photon vtx "<< selected_vertex_index_ << endl;
          H4GTag tag_obj (dipho, phoP4Corrected_dp[0], phoP4Corrected_dp[1], vertex_chosen, dZ_bdtVtx);
          //cout << "breaks here" << endl;
          // cout << "systlabel: " << systLabel_ << endl;
          tag_obj.setSystLabel( systLabel_);
          // tag_obj.setDiPhotonIndex( diphoIndex );
          // tag_obj.setMVA( -0.9 );
          tag_obj.setCategoryNumber( 2 );
          for (int i1=0; i1 < (int) diphoVec.size() ; i1++)
          {
            auto diphotemp = diphoVec[i1];
            tag_obj.includeWeights(* diphotemp);
          }

          H4Gtags->push_back( tag_obj );
          if( ! event.isRealData() ) {
             //cout << "should not enter here " << endl;
            H4Gtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
          }

        }


      }
      event.put( std::move( H4Gtags ) );
      event.put( std::move( truths ) );

    } // H4GTagProducer::produce


  } // namespace flashgg

  typedef flashgg::H4GTagProducer FlashggH4GTagProducer;
  DEFINE_FWK_MODULE( FlashggH4GTagProducer );
