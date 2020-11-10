#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"

#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/H4GTag.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"
#include "flashgg/MicroAOD/interface/MVAComputer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "flashgg/MicroAOD/interface/CutBasedDiPhotonObjectSelector.h"
#include "flashgg/MicroAOD/interface/VertexSelectorBase.h"
#include "flashgg/DataFormats/interface/VertexCandidateMap.h"



#include <vector>
#include <algorithm>
#include "TGraph.h"

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

std::pair<int,int> photonIndices(Handle<View<flashgg::DiPhotonCandidate> > diPhotons, vector<flashgg::Photon> phoP4Corrected_dp)
  {
    std::pair<int,int> pairs;
    int index_dipho1 = -999;
    int index_dipho2 = -999;
    for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ ) {

          double minLeading = 999.;
          double minSubleading = 999.;
          int leadIndex = -1;
          int subleadIndex = -1;

          auto dipho = diPhotons->ptrAt( diphoIndex );


          for( unsigned int phoIndex = 0; phoIndex < phoP4Corrected_dp.size(); phoIndex++ ) {
              double dR_leading = deltaR(dipho->leadingPhoton()->superCluster()->eta(),dipho->leadingPhoton()->superCluster()->phi(),phoP4Corrected_dp[phoIndex].superCluster()->eta(),phoP4Corrected_dp[phoIndex].superCluster()->phi());
              double dR_subleading = deltaR(dipho->subLeadingPhoton()->superCluster()->eta(),dipho->subLeadingPhoton()->superCluster()->phi(),phoP4Corrected_dp[phoIndex].superCluster()->eta(),phoP4Corrected_dp[phoIndex].superCluster()->phi());
               if(dR_leading<minLeading){
                 leadIndex = phoIndex;
                 minLeading = dR_leading;
              }
              if(dR_subleading<minSubleading){
                 subleadIndex = phoIndex;
                 minSubleading = dR_subleading;
              }
          }
          if((leadIndex == 0 && subleadIndex == 1 ) || (leadIndex == 1 && subleadIndex == 0 ) )
          {
            index_dipho1 = diphoIndex;
          }
          if((leadIndex == 2 && subleadIndex == 3) || (leadIndex == 3 && subleadIndex == 2 ))
          {
            index_dipho2 = diphoIndex;
          }

    }
    pairs = std::make_pair(index_dipho1,index_dipho2);
    return pairs;
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
      H4GTagProducer( const ParameterSet & );
      // std::vector<reco::Candidate::LorentzVector> genPhoton_p4;

    private:
      void produce( Event &, const EventSetup & ) override;
      vector<flashgg::Photon> getPhotons(Handle<View<flashgg::DiPhotonCandidate> > diPhotons );
      std::string inputDiPhotonName_;
      std::vector<std::string> inputDiPhotonSuffixes_;
      std::vector<edm::EDGetTokenT<edm::View<DiPhotonCandidate> > > diPhotonTokens_;

      EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
      Handle<View<reco::GenParticle> > genParticle;

      EDGetTokenT<View<Photon> > photonToken_;
      Handle<View<flashgg::Photon> > photons;

      EDGetTokenT<View<reco::Vertex> > vertexToken_;
      Handle<View<reco::Vertex> > vertex;


      //---ID selector
      ConsumesCollector cc_;
      CutBasedDiPhotonObjectSelector idSelector_;

      std::vector< std::string > systematicsLabels;

      EDGetTokenT<reco::BeamSpot>  beamSpotToken_;
      Handle<reco::BeamSpot>  recoBeamSpotHandle;

      EDGetTokenT<View<reco::Conversion> >  conversionToken_;
      Handle<View<reco::Conversion> > conversionHandle;

      EDGetTokenT<View<reco::Conversion> >  conversionTokenSingleLeg_;
      Handle<View<reco::Conversion> > conversionHandleSingleLeg;

      EDGetTokenT< VertexCandidateMap > vertexCandidateMapToken_;
      unique_ptr<VertexSelectorBase> vertexSelector_;

      edm::FileInPath vertexIdMVAweightfileH4G_;

      TMVA::Reader *VertexIdMva_;
      bool useSingleLeg_;
      float logSumpt2;
      float ptAsym;
      float ptBal;
      float pullConv;
      float nConv;
      bool doH4GVertex_;

      std::vector<std::pair<unsigned int, float> > sorter_;
      unsigned int selected_vertex_index_ = 0 ;
      float max_mva_value_  = -999;


    };

    H4GTagProducer::H4GTagProducer( const ParameterSet &iConfig ) :
    genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    photonToken_( consumes<View<Photon> >( iConfig.getParameter<InputTag> ( "PhotonTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
    cc_( consumesCollector() ),
    idSelector_( iConfig.getParameter<ParameterSet> ( "idSelection" ), cc_ ),
    beamSpotToken_( consumes<reco::BeamSpot> ( iConfig.getParameter<InputTag> ( "beamSpotTag" ) ) ),
    conversionToken_( consumes <View<reco::Conversion>> ( iConfig.getParameter<InputTag> ( "conversionTag" ) ) ),
    conversionTokenSingleLeg_( consumes <View<reco::Conversion>> ( iConfig.getParameter<InputTag> ( "conversionTagSingleLeg" ) ) ),
    vertexCandidateMapToken_( consumes<VertexCandidateMap>( iConfig.getParameter<InputTag>( "VertexCandidateMapTag" ) ) )


    {
      inputDiPhotonName_= iConfig.getParameter<std::string > ( "DiPhotonName" );
      inputDiPhotonSuffixes_= iConfig.getParameter<std::vector<std::string> > ( "DiPhotonSuffixes" );
      std::vector<edm::InputTag>  diPhotonTags;
      for (auto & suffix : inputDiPhotonSuffixes_){
        systematicsLabels.push_back(suffix);
        std::string inputName = inputDiPhotonName_;
        inputName.append(suffix);
        if (!suffix.empty()) diPhotonTags.push_back(edm::InputTag(inputName));
        else  diPhotonTags.push_back(edm::InputTag(inputDiPhotonName_));
      }
      for( auto & tag : diPhotonTags ) { diPhotonTokens_.push_back( consumes<edm::View<flashgg::DiPhotonCandidate> >( tag ) ); }

      doH4GVertex_ = iConfig.getParameter<bool>("doH4GVertex");

      if (doH4GVertex_)
      {
        useSingleLeg_ = iConfig.getParameter<bool>( "useSingleLeg" );
        vertexIdMVAweightfileH4G_ = iConfig.getUntrackedParameter<edm::FileInPath>( "vertexIdMVAweightfileH4G" );
        const std::string &VertexSelectorName = iConfig.getParameter<std::string>( "VertexSelectorName" );
        vertexSelector_.reset( FlashggVertexSelectorFactory::get()->create( VertexSelectorName, iConfig ) );

        VertexIdMva_ = new TMVA::Reader( "!Color:Silent" );
        VertexIdMva_->AddVariable( "ptAsym", &ptAsym );
        VertexIdMva_->AddVariable( "ptBal", &ptBal );
        VertexIdMva_->AddVariable( "logSumpt2", &logSumpt2 );
        VertexIdMva_->AddVariable( "pullConv", &pullConv );
        VertexIdMva_->AddVariable( "nConv", &nConv );
        VertexIdMva_->BookMVA( "BDT", vertexIdMVAweightfileH4G_.fullPath() );
      }

      for (auto & systname : systematicsLabels) {
        produces<vector<H4GTag>>(systname);
      }
      produces<vector<TagTruthBase>>();
    }




    void H4GTagProducer::produce( Event &evt, const EventSetup & )
    {

      // prepare output
      std::unique_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );
      edm::RefProd<vector<TagTruthBase> > rTagTruth = evt.getRefBeforePut<vector<TagTruthBase> >();


      // MC truth
      TagTruthBase truth_obj;

      int trueVtxIndexI = -999;
      vector<int>	pvVecNoTrue;
      reco::GenParticle::Point higgsVtx(0.,0.,0.);
      reco::Vertex::Point hardVertex( 0, 0, 0 );
      reco::GenParticle::Point genVertex;

      // std::vector<reco::Candidate::LorentzVector> genPhoton_p4;
      // genPhoton_p4.clear();

      std::vector<reco::Candidate::LorentzVector > genPhos;
      // std::vector<edm::Ptr<reco::GenParticle>> genPhos;
      if( ! evt.isRealData() ) {
        Handle<View<reco::GenParticle> > genParticles;
        evt.getByToken( genParticleToken_, genParticles );

        Handle<View<reco::Vertex> > primaryVertices;
        evt.getByToken( vertexToken_, primaryVertices );

        // std::vector<reco::Candidate::LorentzVector > genPhos;
        trueVtxIndexI = mcTruthVertexIndex_h4g( genParticles->ptrs(), primaryVertices->ptrs(), 0.1);
        for( unsigned int i = 0 ; i < primaryVertices->size() ; i++ ) {
          if( i != (unsigned int)trueVtxIndexI ) { pvVecNoTrue.push_back( i ); }
        }
        for( unsigned int genLoop = 0 ; genLoop < genParticles->ptrs().size(); genLoop++ ) {

          if( fabs( genParticles->ptrs()[genLoop]->pdgId() ) < 10 || fabs( genParticles->ptrs()[genLoop]->pdgId() ) == 25 ) {
            hardVertex.SetCoordinates( genParticles->ptrs()[genLoop]->vx(), genParticles->ptrs()[genLoop]->vy(), genParticles->ptrs()[genLoop]->vz() );
            break;
          }
        }


        for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
          edm::Ptr<reco::GenParticle> part = genParticles->ptrAt(genLoop);
          if (part->pdgId() == 25) {
            genPhos.push_back(part->daughter(0)->p4());
          }
          if( part->pdgId() != 2212 || part->vertex().z() != 0. )
          {
            genVertex = part->vertex();
          }
        }

        for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
          int pdgid = genParticles->ptrAt( genLoop )->pdgId();
          edm::Ptr<reco::GenParticle> part = genParticles->ptrAt(genLoop);
          if( fabs( pdgid ) < 10 || fabs( pdgid ) == 25 ) {
            if (pdgid == 25) {
              higgsVtx = genParticles->ptrAt( genLoop )->vertex();
            }

            break;
          }
        }

        truth_obj.setGenPV(higgsVtx);
        truths->push_back( truth_obj );
      }

      selected_vertex_index_ = 0;
      edm::Ptr<reco::Vertex> vertex_chosen;
      vector <edm::Ptr<reco::Vertex> > Vertices;

      float vtx_X = -999;
      float vtx_Y = -999;
      float vtx_Z = -999;

      float dZ_bdtVtx = -999;
      float dZ_ZeroVtx = -999;
      float dZ_HggVtx = -999;

      double e_veto_dipho_0 = -999;
      double tmpCentralObjectWeight = -999;

      bool atLeastOneDiphoPass_nominal = false;

      for (unsigned int diphoton_idx = 0; diphoton_idx < diPhotonTokens_.size(); diphoton_idx++)
      {
        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        evt.getByToken( diPhotonTokens_[diphoton_idx], diPhotons );
        evt.getByToken(photonToken_, photons);

        evt.getByToken(vertexToken_, vertex );

        evt.getByToken( beamSpotToken_, recoBeamSpotHandle );
        evt.getByToken( conversionToken_, conversionHandle );
        evt.getByToken( conversionTokenSingleLeg_, conversionHandleSingleLeg );
        Handle<VertexCandidateMap> vertexCandidateMap;
        evt.getByToken( vertexCandidateMapToken_, vertexCandidateMap );

        math::XYZPoint BSPoint;
        if( recoBeamSpotHandle.isValid() ) {
          BSPoint = recoBeamSpotHandle->position();
        }

        edm::Ptr<reco::Vertex> vertex_diphoton;
        for( unsigned int dpIndex = 0; dpIndex < diPhotons->size(); dpIndex++ )
        {
          edm::Ptr<flashgg::DiPhotonCandidate> thisDPPtr = diPhotons->ptrAt( dpIndex );
          flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(thisDPPtr.get());
          if (dpIndex == 0) {vertex_diphoton = thisDPPtr->vtx();}
          atLeastOneDiphoPass_nominal |= idSelector_(*thisDPPointer, evt);
        }

        if (atLeastOneDiphoPass_nominal){
          std::vector<edm::Ptr<flashgg::Photon>> phoPtrVector;

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

          max_mva_value_ = -999.;

          std::vector<std::vector<float>> vtxVar;
          sorter_.clear();

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

          for( int vtx = 0; vtx < (int) vertex->size(); vtx++ )
          {
            logSumpt2 = vtxVar[0][vtx];
            ptAsym = vtxVar[1][vtx];
            ptBal = vtxVar[2][vtx];
            pullConv = vtxVar[3][vtx];
            nConv = vtxVar[4][vtx];

            float mva_value_bdt_ = VertexIdMva_->EvaluateMVA( "BDT" );


            std::pair<unsigned int, float>pairToSort = std::make_pair( (unsigned int)vtx, mva_value_bdt_ );
            sorter_.push_back( pairToSort );

            if( mva_value_bdt_ > max_mva_value_ ) {
              max_mva_value_ = mva_value_bdt_;
              selected_vertex_index_ = vtx;
            }
          }

          if (doH4GVertex_)
          {
            vertex_chosen = vertex->ptrAt( selected_vertex_index_ );
            vtx_X = vertex_chosen->x();
            vtx_Y = vertex_chosen->y();
            vtx_Z = vertex_chosen->z();

          }
          else
          {
            vertex_chosen = vertex->ptrAt(0);
          }
          if( ! evt.isRealData() ){
            dZ_bdtVtx =   hardVertex.z() - vertex_chosen->z();
            dZ_ZeroVtx = hardVertex.z() - vertex->ptrAt(0)->z();
            dZ_HggVtx = hardVertex.z() - vertex_diphoton->z();

          }
        }

        if (diphoton_idx > 0) break;
      }

      math::XYZVector vtx_Pos( vtx_X, vtx_Y, vtx_Z );


      // read diphotons
      for (unsigned int diphoton_idx = 0; diphoton_idx < diPhotonTokens_.size(); diphoton_idx++) {//looping over all diphoton systematics

        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        evt.getByToken( diPhotonTokens_[diphoton_idx], diPhotons );

        evt.getByToken(vertexToken_, vertex );


        std::unique_ptr<vector<H4GTag> > tags( new vector<H4GTag> );

        bool atLeastOneDiphoPass = false;
        vector<int> preselDiPhoIndex;


        for( unsigned int dpIndex = 0; dpIndex < diPhotons->size(); dpIndex++ )
        {
          edm::Ptr<flashgg::DiPhotonCandidate> thisDPPtr = diPhotons->ptrAt( dpIndex );
          flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(thisDPPtr.get());
          atLeastOneDiphoPass |= idSelector_(*thisDPPointer, evt);
          if (atLeastOneDiphoPass)
          {
            preselDiPhoIndex.push_back(dpIndex);
          }
        }


        if (atLeastOneDiphoPass_nominal)
        {
          if (atLeastOneDiphoPass){

            edm::Ptr<flashgg::DiPhotonCandidate> dipho;

            vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoVec;
            // vector<flashgg::DiPhotonCandidate> diphoPtrVec;
            vector <flashgg::Photon> vecPho = getPhotons(diPhotons);

            for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ ) {

              dipho = diPhotons->ptrAt( diphoIndex ); // without systematic look (?)
              diphoVec.push_back(dipho);

            }

            vector<flashgg::Photon> phoP4Corrected_dp;
            vector<const flashgg::Photon *> phoP4Corrected_dp_Ptr;
            if (vecPho.size() > 0)
            {
              for (int dp = 0; dp < (int) vecPho.size(); dp++)
              {
                float sc_X_dp = vecPho[dp].superCluster()->x();
                float sc_Y_dp = vecPho[dp].superCluster()->y();
                float sc_Z_dp = vecPho[dp].superCluster()->z();
                math::XYZVector sc_Pos_dp( sc_X_dp, sc_Y_dp, sc_Z_dp );
                math::XYZVector direction_dp = sc_Pos_dp - vtx_Pos;
                math::XYZVector pho_dp = ( direction_dp.Unit() ) * ( vecPho[dp].energy() );
                math::XYZTLorentzVector corrected_p4_dp( pho_dp.x(), pho_dp.y(), pho_dp.z(), vecPho[dp].energy() );
                vecPho[dp].setP4(corrected_p4_dp);
                phoP4Corrected_dp.push_back(vecPho[dp]);
                // phoP4Corrected_dp_Ptr.push_back(&vecPho[dp]);
              }
            }

            sorter_.clear();
            std::sort(phoP4Corrected_dp.begin(), phoP4Corrected_dp.end(), [](const flashgg::Photon a, const flashgg::Photon  b) {return a.pt() > b.pt(); });





            if (phoP4Corrected_dp.size() > 3){
              H4GTag tag_obj (dipho, phoP4Corrected_dp[0], phoP4Corrected_dp[1], phoP4Corrected_dp[2], phoP4Corrected_dp[3], vertex_chosen, dZ_bdtVtx, dZ_ZeroVtx, dZ_HggVtx,genPhos );

              vector<flashgg::Photon> phoP4Corrected_dp_reduced;

              for (int p=0; p < 4; p++)
              {
                phoP4Corrected_dp_reduced.push_back(phoP4Corrected_dp[p]);
              }
             std::pair<int,int> pairs = photonIndices(diPhotons,phoP4Corrected_dp_reduced);
              tag_obj.setCategoryNumber( 0 );
              tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
              //cout << "***********************" << endl;
              //cout << "Weights before including weights: " << endl;
              //for (auto it = tag_obj.weightListBegin() ; it != tag_obj.weightListEnd(); it++){
              //  cout << "weight " << *it << " " << tag_obj.weight(*it) << endl;
              //}
cout << "********************************************" << endl;
              e_veto_dipho_0 = -999;
              tmpCentralObjectWeight = -999;
              tag_obj.includeWeights( *diphoVec[preselDiPhoIndex[0]] ); // Add highest pT diphoton weights to central object weight
              tmpCentralObjectWeight = tag_obj.centralWeight(); // get central object weight
              cout << "central weight before removing electron veto SF: " << tmpCentralObjectWeight << endl;

              e_veto_dipho_0 = tag_obj.weight("electronVetoSFCentral"); // get electron veto SF from 0th diphoton
              cout << "e veto dipho 0 value : " << e_veto_dipho_0 << endl;
              tag_obj.setCentralWeight(tmpCentralObjectWeight / e_veto_dipho_0); // remove electron veto of 0th pT contribution from central weight

              cout << "central event weight after dividing out e veto dipho 0 : " << tag_obj.centralWeight() << endl;

              // Reset electronVetoSF weights because we want weights from 4 photons here not from two photons from leading diphoton
              tag_obj.setWeight("electronVetoSFDown01sigma",1);
              tag_obj.setWeight("electronVetoSFCentral",1);
              tag_obj.setWeight("electronVetoSFUp01sigma",1);


              cout << "printing out electronVetoSF* weights after setting explictly to 1 " << endl;
              cout << "electronVetoSFDown01Sigma: " << tag_obj.weight("electronVetoSFDown01sigma") << endl;
              cout << "electronVetoSFCentral : " << tag_obj.weight("electronVetoSFCentral") << endl;
              cout << "electronVetoSFUp01Sigma: " << tag_obj.weight("electronVetoSFUp01sigma") << endl;
              cout << "central event weight after setting electronVetoSF* weights to 1: " << tag_obj.centralWeight() << endl;
              //tag_obj.setCentralWeight(1);
              //cout << "electronVetoSFDown01Sigma: " << tag_obj.weight("electronVetoSFDown01Sigma") << endl;
              //cout << "electronVetoSFCentral : " << tag_obj.weight("electronVetoSFCentral") << endl;
              //cout << "electronVetoSFUp01Sigma: " << tag_obj.weight("electronVetoSFUp01Sigma") << endl;
              //cout << "central weight after setting it to 1: "  << tag_obj.centralWeight() << endl;
              // Include electron veto SF of photons used by tagger, while simultaneously updating electronVetoSF* branch values
              tag_obj.includeWeightsByLabel(*diphoVec[std::get<0>(pairs)],"electronVetoSFCentral");
              tag_obj.includeWeightsByLabel(*diphoVec[std::get<1>(pairs)],"electronVetoSFCentral");

              tag_obj.includeWeightsByLabel(*diphoVec[std::get<0>(pairs)],"electronVetoSFDown01sigma");
              tag_obj.includeWeightsByLabel(*diphoVec[std::get<1>(pairs)],"electronVetoSFDown01sigma");

              tag_obj.includeWeightsByLabel(*diphoVec[std::get<0>(pairs)],"electronVetoSFUp01sigma");
              tag_obj.includeWeightsByLabel(*diphoVec[std::get<1>(pairs)],"electronVetoSFUp01sigma");

              cout << "trying to get electronVetoSF up and down 1 sigmas for paired diphotons" << endl;


              // Set tag object's electronVetoSF weight equal to the four photons used in the H4G event
              //tag_obj.setWeight("electronVetoSFCentral",1);
              // check to see what setweight is doing
              cout << "electronVetoSFDown01Sigma: " << tag_obj.weight("electronVetoSFDown01sigma") << endl;
              cout << "electronVetoSFCentral : " << tag_obj.weight("electronVetoSFCentral") << endl;
              cout << "electronVetoSFUp01Sigma: " << tag_obj.weight("electronVetoSFUp01sigma") << endl;
              // check that central weight is updated
              cout << "central event weight : " << tag_obj.centralWeight() << endl;
              // check the prefire weight
              cout << "prefireWeightUp01Sigma : " << tag_obj.weight("prefireWeightUp01sigma") << endl;
              cout << "prefireWeightDown01Sigma : " << tag_obj.weight("prefireWeightDown01sigma") << endl;
cout << "******************************************************* " << endl;
//tag_obj.includeWeightsByLabel( *diphoVec[preselDiPhoIndex[0]] ,"PreselSF");
              //tag_obj.includeWeightsByLabel( *diphoVec[preselDiPhoIndex[0]] ,"TriggerWeight");
              //tag_obj.includeWeightsByLabel( *diphoVec[preselDiPhoIndex[0]] ,"prefireWeight",false);
              //tag_obj.includeWeightsByLabel(*diphoVec[std::get<0>(pairs)],"electronVetoSF");
              //tag_obj.includeWeightsByLabel(*diphoVec[std::get<1>(pairs)],"electronVetoSF");

              //tag_obj.setCentralWeight(1);
              //cout << "tag_obj.weight(PreselSF):" << tag_obj.weight("PreselSF") << endl;
              //cout << "tag_obj.centralWeight(PreselSFCentral):" << tag_obj.centralWeight("PreselSFCentral")
              //cout << "tag_obj.weight('PreselSFCentral'):" << tag_obj.weight("PreselSFCentral") << endl;
              //cout << "tag_obj.centralWeight():" << tag_obj.centralWeight() << endl;

              //cout << "************************" << endl;
              //cout << "After including weights" << endl;

              //for (auto it = tag_obj.weightListBegin() ; it != tag_obj.weightListEnd(); it++) {
              //          std::cout << " Weight Debug " << *it << " " << tag_obj.weight(*it) << std::endl;

               //     }


              //cout << "now trying to set central weight" << endl;
              //cout << "tag_obj.setCentralWeight(1)" << endl;
              //tag_obj.setCentralWeight(1);
              //cout << "get central weight: " << tag_obj.weight("Central") << endl;

              //cout << "now trying to set central weight" << endl;
              //cout << "tag_obj.setCentralWeight(1)" << endl;
              //tag_obj.setCentralWeight(1);
              //cout << "get central weight: " << tag_obj.weight("Central") << endl;

              //for (auto it = tag_obj.weightListBegin() ; it != tag_obj.weightListEnd(); it++){
              //      cout << "Weight Debug " << *it << " " << tag_obj.weight(*it) << endl;
              //}

              //cout << "checking if there is a PreselSF weight now " << endl;
              //cout << "tag_obj.hasWeight(PreselSF) " << tag_obj.hasWeight("PreselSF") << endl;


              tags->push_back(tag_obj);
              if (!evt.isRealData())
              {
                tags->back().setTagTruth(edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ));
              }
            }

            if (phoP4Corrected_dp.size() == 3 )
            {
              H4GTag tag_obj (dipho, phoP4Corrected_dp[0], phoP4Corrected_dp[1], phoP4Corrected_dp[2],vertex_chosen, dZ_bdtVtx, dZ_ZeroVtx, dZ_HggVtx,genPhos);
              tag_obj.setCategoryNumber( 1 );
              tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
              tag_obj.includeWeightsByLabel( *diphoVec[preselDiPhoIndex[0]] ,"PreselSF");
              tag_obj.includeWeightsByLabel( *diphoVec[preselDiPhoIndex[0]] ,"TriggerWeight");


              tags->push_back( tag_obj );
              if( ! evt.isRealData() ) {
                tags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
              }

            }

            if (phoP4Corrected_dp.size() == 2 )
            {
              H4GTag tag_obj (dipho, phoP4Corrected_dp[0], phoP4Corrected_dp[1], vertex_chosen, dZ_bdtVtx, dZ_ZeroVtx, dZ_HggVtx,genPhos);
              tag_obj.setCategoryNumber( 2 );
              tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
              tag_obj.includeWeightsByLabel( *diphoVec[preselDiPhoIndex[0]] ,"PreselSF");
              tag_obj.includeWeightsByLabel( *diphoVec[preselDiPhoIndex[0]] ,"TriggerWeight");

              tags->push_back( tag_obj );
              if( ! evt.isRealData() ) {
                tags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
              }
            }
          }
        }
        evt.put( std::move( tags ),inputDiPhotonSuffixes_[diphoton_idx] );

      }
      evt.put( std::move( truths ) );
    }

/*vector<flashgg::Photon> H4GTagProducer::getPhotons(Handle<View<flashgg::DiPhotonCandidate> > diPhotons )
    {
      vector<flashgg::Photon> phoVector;
      int n_pho = 0;
      for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ ) {
          edm::Ptr<flashgg::DiPhotonCandidate> thisDiPho = diPhotons->ptrAt( diphoIndex );
        flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(thisDiPho.get());
        thisDPPointer->makePhotonsPersistent();
        auto pho1 = thisDPPointer->getLeadingPhoton();
        auto pho2 = thisDPPointer->getSubLeadingPhoton();

       // if (phoVector.size() == 0)
        //  {
       //   phoVector.push_back(pho1);
      //    std::sort(phoVector.begin(), phoVector.end(), [](const flashgg::Photon a, const flashgg::Photon b) {return a.pt() > b.pt(); });
      std::sort(phoVector.begin(), phoVector.end(), [](const flashgg::Photon a, const flashgg::Photon b) {return a.pt() > b.pt(); });

      return phoVector;

    }
  }
*/

    vector<flashgg::Photon> H4GTagProducer::getPhotons(Handle<View<flashgg::DiPhotonCandidate> > diPhotons )
    {
      vector<flashgg::Photon> phoVector;
      int n_pho = 0;
      for( unsigned int diphoIndex = 0; diphoIndex < diPhotons->size(); diphoIndex++ ) {

        // create a vector of unique photons from all diphotons
        edm::Ptr<flashgg::DiPhotonCandidate> thisDiPho = diPhotons->ptrAt( diphoIndex );
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
            float deltar1 = deltaR(phoVector[p].superCluster()->eta(),phoVector[p].superCluster()->phi(),pho1.superCluster()->eta(),pho1.superCluster()->phi());
            float deltar2 = deltaR(phoVector[p].superCluster()->eta(),phoVector[p].superCluster()->phi(),pho2.superCluster()->eta(),pho2.superCluster()->phi());
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

      return phoVector;

    }
  }

  typedef flashgg::H4GTagProducer FlashggH4GTagProducer;
  DEFINE_FWK_MODULE( FlashggH4GTagProducer );
  // Local Variables:
  // mode:c++
  // indent-tabs-mode:nil
  // tab-width:4
  // c-basic-offset:4
  // End:
  // vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
