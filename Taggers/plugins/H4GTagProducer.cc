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

#include "flashgg/MicroAOD/interface/CutBasedDiPhotonObjectSelector.h"

#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
#include "flashgg/DataFormats/interface/WHLeptonicTag.h"

#include "flashgg/Taggers/interface/LeptonSelection.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include <vector>
#include <algorithm>
#include "TGraph.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;

int mcTrueVertexIndex( const std::vector<edm::Ptr<reco::GenParticle> > &genParticles ,
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
    class H4GTagProducer : public EDProducer
    {
    public:
      //---typedef
      typedef math::XYZTLorentzVector LorentzVector;

      //---ctors
      // HHWWggTagProducer();
      H4GTagProducer( const ParameterSet & );

      //---Outtree
      edm::Service<TFileService> fs;

      std::vector<reco::Candidate::LorentzVector> genPhoton_p4;

    private:
      double genTotalWeight;
      bool checkPassMVAs(const flashgg::Photon*& leading_photon, const flashgg::Photon*& subleading_photon, edm::Ptr<reco::Vertex>& diphoton_vertex);
      void produce( Event &, const EventSetup & ) override;
      std::vector<edm::EDGetTokenT<edm::View<DiPhotonCandidate> > > diPhotonTokens_;
      std::string inputDiPhotonName_;

      std::string inputJetsName_;
      std::vector<std::string> inputJetsSuffixes_;
      unsigned int inputJetsCollSize_;

      // Adding Jets
      std::vector<edm::EDGetTokenT<edm::View<flashgg::Jet> > > jetTokens_;

      EDGetTokenT<View<Photon> > photonToken_;
      Handle<View<flashgg::Photon> > photons;

      EDGetTokenT<View<DiPhotonCandidate> > diphotonToken_;
      Handle<View<flashgg::DiPhotonCandidate> > diphotons;

      EDGetTokenT<View<reco::Vertex> > vertexToken_;
      Handle<View<reco::Vertex> > vertex;

      EDGetTokenT< VertexCandidateMap > vertexCandidateMapToken_;
      unique_ptr<VertexSelectorBase> vertexSelector_;

      EDGetTokenT<View<reco::Conversion> >  conversionToken_;
      Handle<View<reco::Conversion> > conversionHandle;
      
      EDGetTokenT<View<reco::Conversion> >  conversionTokenSingleLeg_;
      Handle<View<reco::Conversion> > conversionHandleSingleLeg;

      EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
      Handle<View<reco::GenParticle> > genParticle;

      EDGetTokenT<reco::BeamSpot>  beamSpotToken_;
      Handle<reco::BeamSpot>  recoBeamSpotHandle;

      EDGetTokenT<View<Electron> > electronToken_;
      Handle<View<flashgg::Electron> > electrons;

      EDGetTokenT<View<Muon> > muonToken_;
      Handle<View<flashgg::Muon> > muons;

      EDGetTokenT<View<flashgg::Met> > METToken_;
      Handle<View<flashgg::Met> > METs;

      EDGetTokenT<View<DiPhotonMVAResult> > mvaResultToken_;
      Handle<View<flashgg::DiPhotonMVAResult> > mvaResults;

      Handle<View<reco::Vertex> > vertices;

      EDGetTokenT<double> rhoTag_;
      edm::EDGetTokenT<edm::TriggerResults> triggerRECO_;
      edm::EDGetTokenT<edm::TriggerResults> triggerPAT_;
      edm::EDGetTokenT<edm::TriggerResults> triggerFLASHggMicroAOD_;
      string systLabel_;
      edm::Handle<double>  rho;

      // std::vector<edm::EDGetTokenT<View<flashgg::Jet> > > JetToken_;


      std::vector< std::string > systematicsLabels;
      std::vector<std::string> inputDiPhotonSuffixes_;

      //---ID selector
      ConsumesCollector cc_;
      GlobalVariablesComputer globalVariablesComputer_;
      CutBasedDiPhotonObjectSelector idSelector_;


      bool useSingleLeg_;
      //----output collection
      // auto_ptr<vector<HHWWggCandidate> > HHWWggColl_;

      // variables from WHLeptonicTagProducer
      double leptonPtThreshold_;
      double muonEtaThreshold_;
      double leadPhoOverMassThreshold_;
      double subleadPhoOverMassThreshold_;
      double MVAThreshold_;
      double deltaRMuonPhoThreshold_;
      double jetsNumberThreshold_;
      double jetPtThreshold_;
      double jetEtaThreshold_;
      double muPFIsoSumRelThreshold_;
      double PhoMVAThreshold_;
      double METThreshold_;
      bool useVertex0only_;
      double deltaRJetMuonThreshold_;
      double deltaRPhoLeadJet_;
      double deltaRPhoSubLeadJet_;

      double DeltaRTrkElec_;
      double TransverseImpactParam_;
      double LongitudinalImpactParam_;

      double deltaRPhoElectronThreshold_;
      double deltaMassElectronZThreshold_;

      bool hasGoodElec = false;
      bool hasGoodMuons = false;

      vector<double> nonTrigMVAThresholds_;
      vector<double> nonTrigMVAEtaCuts_;

      double electronIsoThreshold_;
      double electronNumOfHitsThreshold_;
      vector<double> electronEtaThresholds_;
      bool useElectronMVARecipe_;
      bool useElectronLooseID_;
      // string bTag_;
      double btagThresh_;
      bool doHHWWggTagCutFlowAnalysis_;


      edm::InputTag genInfo_;
      edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
    };

    //---constructors
    // HHWWggTagProducer::HHWWggTagProducer( ):
    // photonToken_(),
    // diphotonToken_()
    // genParticleToken_(),
    // electronToken_(),
    // muonToken_(),
    // METToken_(),
    // cc_( consumesCollector() )
    // // idSelector_( ParameterSet(), cc_ )

    // {}

    //---standard
    H4GTagProducer::H4GTagProducer( const ParameterSet & pSet):
    photonToken_( consumes<View<Photon> >( pSet.getParameter<InputTag> ( "PhotonTag" ) ) ),
    diphotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( pSet.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( pSet.getParameter<InputTag> ( "VertexTag" ) ) ),
    vertexCandidateMapToken_( consumes<VertexCandidateMap>( pSet.getParameter<InputTag>( "VertexCandidateMapTag" ) ) ),
    conversionToken_( consumes <View<reco::Conversion>> ( pSet.getParameter<InputTag> ( "conversionTag" ) ) ),
    conversionTokenSingleLeg_( consumes <View<reco::Conversion>> ( pSet.getParameter<InputTag> ( "conversionTagSingleLeg" ) ) ),
    genParticleToken_( consumes<View<reco::GenParticle> >( pSet.getParameter<InputTag> ( "GenParticleTag" ) ) ),
    beamSpotToken_( consumes<reco::BeamSpot> ( pSet.getParameter<InputTag> ( "beamSpotTag" ) ) ),
    electronToken_( consumes<View<Electron> >( pSet.getParameter<InputTag> ( "ElectronTag" ) ) ),
    muonToken_( consumes<View<Muon> >( pSet.getParameter<InputTag> ( "MuonTag" ) ) ),
    METToken_( consumes<View<Met> >( pSet.getParameter<InputTag> ( "METTag" ) ) ),
    mvaResultToken_( consumes<View<flashgg::DiPhotonMVAResult> >( pSet.getParameter<InputTag> ( "MVAResultTag" ) ) ),
    rhoTag_( consumes<double>( pSet.getParameter<InputTag>( "rhoTag" ) ) ),
    triggerRECO_( consumes<edm::TriggerResults>(pSet.getParameter<InputTag>("RECOfilters") ) ),
    triggerPAT_( consumes<edm::TriggerResults>(pSet.getParameter<InputTag>("PATfilters") ) ),
    triggerFLASHggMicroAOD_( consumes<edm::TriggerResults>( pSet.getParameter<InputTag>("FLASHfilters") ) ),
    systLabel_( pSet.getParameter<string> ( "SystLabel" ) ),
    cc_( consumesCollector() ),
    globalVariablesComputer_(pSet.getParameter<edm::ParameterSet>("globalVariables"), cc_),
    idSelector_( pSet.getParameter<ParameterSet> ( "idSelection" ), cc_ )

    {
      const std::string &VertexSelectorName = pSet.getParameter<std::string>( "VertexSelectorName" );
      vertexSelector_.reset( FlashggVertexSelectorFactory::get()->create( VertexSelectorName, pSet ) );
      useSingleLeg_ = pSet.getParameter<bool>( "useSingleLeg" );

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

      bool breg = 0;

      inputJetsName_= pSet.getParameter<std::string> ( "JetsName" );
      inputJetsCollSize_= pSet.getParameter<unsigned int> ( "JetsCollSize" );
      inputJetsSuffixes_= pSet.getParameter<std::vector<std::string> > ( "JetsSuffixes" );
      // cout << "inputJetsCollSize_ = " << inputJetsCollSize_ << endl;
      if (breg){
        std::vector<edm::InputTag>  jetTags; // With bregression on
        for (auto & suffix : inputJetsSuffixes_) {
          if (!suffix.empty()) systematicsLabels.push_back(suffix);  //nominal is already put in the diphoton loop
          for (unsigned int i = 0; i < inputJetsCollSize_ ; i++) {
            std::string bregtag = suffix;
            bregtag.append(std::to_string(i));
            if (breg) jetTags.push_back(edm::InputTag(inputJetsName_,bregtag)); // With bregression on
          }
        }

        for( auto & tag : jetTags ) { jetTokens_.push_back( consumes<edm::View<flashgg::Jet> >( tag ) ); } // With bregression on
      }

      // Jets without bregression
      if (!breg){
        auto jetTags = pSet.getParameter<std::vector<edm::InputTag> > ( "JetTags" );
        for( auto & tag : jetTags ) { jetTokens_.push_back( consumes<edm::View<flashgg::Jet> >( tag ) ); }
      }


      genInfo_ = pSet.getUntrackedParameter<edm::InputTag>( "genInfo", edm::InputTag("generator") );
      genInfoToken_ = consumes<GenEventInfoProduct>( genInfo_ );
      // numDiphoCand = fs->make<TH1F> ("numDiphoCand","numDiphoCand",10,0,10);
      // diphoton_idx_h = fs->make<TH1F> ("diphoton_idx_h","diphoton_idx_h",20,0,20);
      // diPhotons_size_h = fs->make<TH1F> ("diPhotons_size_h","diPhotons_size_h",20,0,20);

      // indexes = fs->make<TH1F> ("indexes","indexes",5,0,5);
      // btags = fs->make<TH1F> ("btags","btags",100,0,1);

      // numEvents = fs->make<TH1F> ("numEvents","numEvents",1,0,10);

      // gen_weights = fs->make<TH1F> ("gen_weights","gen_weights",1000,-2,2);
      // vars = fs->make<TH1F> ("vars","vars",10,0,10);
      // cutFlow = fs->make<TH1F> ("cutFlow","Cut Flow",10,0,10);
      // WTags = fs->make<TH1F> ("WTags","W Tags",3,0,3);

      leptonPtThreshold_ = pSet.getParameter<double>( "leptonPtThreshold");
      muonEtaThreshold_ = pSet.getParameter<double>( "muonEtaThreshold");
      leadPhoOverMassThreshold_ = pSet.getParameter<double>( "leadPhoOverMassThreshold");
      subleadPhoOverMassThreshold_ = pSet.getParameter<double>( "subleadPhoOverMassThreshold");
      MVAThreshold_ = pSet.getParameter<double>( "MVAThreshold");
      deltaRMuonPhoThreshold_ = pSet.getParameter<double>( "deltaRMuonPhoThreshold");
      jetsNumberThreshold_ = pSet.getParameter<double>( "jetsNumberThreshold");
      jetPtThreshold_ = pSet.getParameter<double>( "jetPtThreshold");
      jetEtaThreshold_ = pSet.getParameter<double>( "jetEtaThreshold");
      muPFIsoSumRelThreshold_ = pSet.getParameter<double>( "muPFIsoSumRelThreshold");
      PhoMVAThreshold_ = pSet.getParameter<double>( "PhoMVAThreshold");
      METThreshold_ = pSet.getParameter<double>( "METThreshold");
      useVertex0only_              = pSet.getParameter<bool>("useVertex0only");
      deltaRJetMuonThreshold_ = pSet.getParameter<double>( "deltaRJetMuonThreshold");
      deltaRPhoLeadJet_ = pSet.getParameter<double>( "deltaRPhoLeadJet");
      deltaRPhoSubLeadJet_ = pSet.getParameter<double>( "deltaRPhoSubLeadJet");

      DeltaRTrkElec_ = pSet.getParameter<double>( "DeltaRTrkElec");
      TransverseImpactParam_ = pSet.getParameter<double>( "TransverseImpactParam");
      LongitudinalImpactParam_ = pSet.getParameter<double>( "LongitudinalImpactParam");

      deltaRPhoElectronThreshold_ = pSet.getParameter<double>( "deltaRPhoElectronThreshold");
      deltaMassElectronZThreshold_ = pSet.getParameter<double>( "deltaMassElectronZThreshold");

      nonTrigMVAThresholds_ =  pSet.getParameter<vector<double > >( "nonTrigMVAThresholds");
      nonTrigMVAEtaCuts_ =  pSet.getParameter<vector<double > >( "nonTrigMVAEtaCuts");
      electronIsoThreshold_ = pSet.getParameter<double>( "electronIsoThreshold");
      electronNumOfHitsThreshold_ = pSet.getParameter<double>( "electronNumOfHitsThreshold");
      electronEtaThresholds_ = pSet.getParameter<vector<double > >( "electronEtaThresholds");
      useElectronMVARecipe_=pSet.getParameter<bool>("useElectronMVARecipe");
      useElectronLooseID_=pSet.getParameter<bool>("useElectronLooseID");
      // bTag_ = pSet.getParameter<string> ( "bTag");
      btagThresh_ = pSet.getParameter<double>( "btagThresh");
      doHHWWggTagCutFlowAnalysis_ = pSet.getParameter<bool>( "doHHWWggTagCutFlowAnalysis");

      produces<vector<H4GTag>>();
      // for (auto & systname : systematicsLabels) { // to deal with systematics in producer
      //     produces<vector<HHWWggTag>>(systname);
      // }
      produces<vector<TagTruthBase>>();
      // cout << "**************************** in HHWWggTagProducer.cc **********************************************" << endl;
    }

    // bool HHWWggTagProducer::PassMVASelections()
    // {
    //
    // }

    bool H4GTagProducer::checkPassMVAs( const flashgg::Photon*& leading_photon, const flashgg::Photon*& subleading_photon, edm::Ptr<reco::Vertex>& diphoton_vertex){

      bool debug_mva = 0;

      // MVA Check variables
      double lp_mva_thresh = 0.07;
      double slp_mva_thresh = -0.03;

      bool lead_pass_TightPhoID = 0, sublead_pass_TightPhoID = 0;
      double lp_Hgg_MVA = -99, slp_Hgg_MVA = -99;
      double leading_pho_eta = -99, sub_leading_pho_eta = -99;

      // Get MVA values wrt diphoton vertex
      lp_Hgg_MVA = leading_photon->phoIdMvaDWrtVtx( diphoton_vertex );
      slp_Hgg_MVA = subleading_photon->phoIdMvaDWrtVtx( diphoton_vertex );

      // Get eta values
      leading_pho_eta = leading_photon->p4().eta();
      sub_leading_pho_eta = subleading_photon->p4().eta();

      // Debug
      if(debug_mva){
        cout << "leading mva: " << lp_Hgg_MVA << endl;
        cout << "subleading mva: " << slp_Hgg_MVA << endl;
      }
      // leading photon
      // EB
      if (( abs(leading_pho_eta) > 0) && ( abs(leading_pho_eta) < 1.4442)){
        // if (lead_pho_EG_MVA_ > 0.42) lead_pass_TightPhoID = 1;
        if (lp_Hgg_MVA > lp_mva_thresh) lead_pass_TightPhoID = 1;
      }

      // EE
      else if (( abs(leading_pho_eta) > 1.566) && ( abs(leading_pho_eta) < 2.5)){
        // if (lead_pho_EG_MVA_ > 0.14) lead_pass_TightPhoID = 1;
        if (lp_Hgg_MVA > slp_mva_thresh) lead_pass_TightPhoID = 1;
      }

      // SubLeading Photon
      // EB
      if (( abs(sub_leading_pho_eta) > 0) && ( abs(sub_leading_pho_eta) < 1.4442)){
        // if (sublead_pho_EG_MVA_ > 0.42) sublead_pass_TightPhoID = 1;
        if (slp_Hgg_MVA > lp_mva_thresh) sublead_pass_TightPhoID = 1;
      }

      // EE
      else if (( abs(sub_leading_pho_eta) > 1.566) && ( abs(sub_leading_pho_eta) < 2.5)){
        // if (sublead_pho_EG_MVA_ > 0.14) sublead_pass_TightPhoID = 1;
        if (slp_Hgg_MVA > slp_mva_thresh) sublead_pass_TightPhoID = 1;
      }

      if (lead_pass_TightPhoID && sublead_pass_TightPhoID){
        if(debug_mva) cout << "PASS MVA Selections" << endl;
        return 1;
      }

      else{
        if(debug_mva) cout << "FAIL MVA Selections" << endl;
        return 0;
      }

    }

    void H4GTagProducer::produce( Event &event, const EventSetup & )
    {

      // cout << "[HHWWggTagProducer.cc] - Beginning of HHWWggTagProducer::produce" << endl;

      // update global variables
      // globalVariablesComputer_.update(event);

      // Get particle objects
      event.getByToken( photonToken_, photons );
      event.getByToken( diphotonToken_, diphotons );
      event.getByToken( genParticleToken_, genParticle );
      event.getByToken( electronToken_, electrons );
      event.getByToken( muonToken_, muons );
      event.getByToken( METToken_, METs );
      event.getByToken( mvaResultToken_, mvaResults );
      event.getByToken( vertexToken_, vertex );
      event.getByToken( rhoTag_, rho);
      event.getByToken( conversionToken_, conversionHandle );
      event.getByToken( conversionTokenSingleLeg_, conversionHandleSingleLeg );
      Handle<VertexCandidateMap> vertexCandidateMap;
      event.getByToken( vertexCandidateMapToken_, vertexCandidateMap );
      event.getByToken( beamSpotToken_, recoBeamSpotHandle );

      double rho_    = *rho;
      int hgg_index = -999;
      genPhoton_p4.clear();

      math::XYZPoint BSPoint;
      if( recoBeamSpotHandle.isValid() ) {
        BSPoint = recoBeamSpotHandle->position();
      }
      // Set cut booleans
      // std::vector<double> Cut_Results = {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // Cut_Results[i] = 1: Event Passed Cut i
      std::vector<double> Cut_Variables(20,0.0); // Cut_Results[i] = 1.0: Event Passed Cut i
      // std::vector<double> Vertex_Variables(20,0.0); // Cut_Results[i] = 1.0: Event Passed Cut i

      // Cut Variables
      // double has_PS_Dipho = 0, pass_METfilters = 0, dipho_vertex_is_zero = 0, pass_leadPhoOverMassThreshold = 0, pass_subleadPhoOverMassThreshold = 0,
      //   pass_LeadPhoton_MVA = 0, pass_SubLeadPhoton_MVA = 0, pass_dipho_MVA = 0, number_passed_jetid = 0;
      // double dipho_vertex_is_zero = -999;
      // double SLW_Tag = 0.; // Semi-Leptonic W Tag
      // double FLW_Tag = 0.; // Fully-Leptonic W Tag
      // double FHW_Tag = 0.; // Fully-Hadronic W Tag
      // bool PS_dipho_tag = 0; // preselected diphoton

      //---output collection
      // std::unique_ptr<vector<HHWWggCandidate> > HHWWggColl_( new vector<HHWWggCandidate> );
      // std::unique_ptr<vector<HHWWggTag> > tags( new vector<HHWWggTag> );
      // int n_METs = METs->size(); // Should be 1, but using as a way to obtain met four vector
      int n_good_electrons = 0;
      int n_good_muons = 0;
      int n_good_leptons = 0;
      int n_good_jets = 0;
      bool hasHighbTag = 0;
      float btagVal = 0;
      // double dipho_MVA = -99;
      // double lead_pho_Hgg_MVA = -99, sublead_pho_Hgg_MVA = -99;
      // double CMS_hgg_mass = -99;
      // float bDiscriminatorValue = -2.;

      bool passMVAs = 0; // True if leading and subleading photons pass MVA selections

      // Saved Objects after selections
      std::vector<flashgg::Jet> tagJets_;
      std::vector<flashgg::Muon> goodMuons_;
      std::vector<flashgg::Electron> goodElectrons_;
      std::vector<flashgg::Met> theMET_;
      std::vector<flashgg::DiPhotonCandidate> diphoVector_;
      reco::GenParticle::Point genVertex;
      Handle<View<reco::Vertex> > primaryVertices;
      event.getByToken( vertexToken_, primaryVertices );

      int trueVtxIndexI = -999;
      vector<int>	pvVecNoTrue;
      int irand = -999;
      int randVtxIndexI = -999;
      //--------vertex output-------//
      vector <edm::Ptr<reco::Vertex> > Vertices; // Collection of vertices
      vector <edm::Ptr<reco::Vertex> > slim_Vertices;

      std::vector<edm::Ptr<reco::GenParticle>> genPhos;
      std::unique_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );
      edm::RefProd<vector<TagTruthBase> > rTagTruth = event.getRefBeforePut<vector<TagTruthBase> >();


      //-----------------------------------------------------------------------------------------------------------

      // Vertex variables
      // double gen_vertex_z = -999;
      // double hgg_vertex_z = -999;
      // double zero_vertex_z = -999;
      // double vertex_diff_zeroeth = -999;
      // double vertex_diff_hgg = -999;
      // double num_vertices = -999;

      // int diphoton_vertex_index = -99;
      // const edm::Ptr<reco::Vertex> dipho_vertex;
      edm::Ptr<reco::Vertex> diphoton_vertex;
      // edm::Ptr<reco::Vertex> zero_vertex;
      // num_vertices = (double)vertices->size();
      // cout << "vertices->size() = " << vertices->size() << endl;
      // if (vertices->size() > 0){
      //   zero_vertex = vertices->ptrAt( 0 );
      // }

      // MC truth
      TagTruthBase truth_obj;
      // double genMhh=0.;
      if( ! event.isRealData() ) {
        Handle<View<reco::GenParticle> > genParticles;
        std::vector<edm::Ptr<reco::GenParticle> > selHiggses;
        event.getByToken( genParticleToken_, genParticles );
        reco::GenParticle::Point higgsVtx(0.,0.,0.);
        trueVtxIndexI = mcTrueVertexIndex( genParticles->ptrs(), primaryVertices->ptrs(), 0.1);

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

        for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
          int pdgid = genParticles->ptrAt( genLoop )->pdgId();
          if( pdgid == 25 || pdgid == 22 ) { // not so sure if this is correct for HHWWgg because of potential photons from hadronization
            higgsVtx = genParticles->ptrAt( genLoop )->vertex();
            // gen_vertex_z = higgsVtx.z();
            break;
          }
        }
        for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
          edm::Ptr<reco::GenParticle> genPar = genParticles->ptrAt(genLoop);
          // if (selHiggses.size()>1) break;
          // if (genPar->pdgId()==25 && genPar->isHardProcess()){
          //   selHiggses.push_back(genPar);
          // }
          edm::Ptr<reco::GenParticle> part = genParticles->ptrAt(genLoop);
          if (genPar->pdgId() ==  25 || genPar->pdgId() == 54)
          {
            if (genPar->daughter(0)->pdgId() == 22 && genPar->daughter(1)->pdgId() == 22)
            {
              genPhos.push_back(part);
              genPhoton_p4.push_back(genPar->daughter(0)->p4());
              genPhoton_p4.push_back(genPar->daughter(1)->p4());
            }
          }
        }

        truth_obj.setGenPV( higgsVtx );
        truths->push_back( truth_obj );
      }

      // // Get Gen vertex
      // bool got_gen_vertex = 0;
      // if (! event.isRealData()){
      //   for( auto &part : *genParticle ) {
      //     if( ( part.pdgId() != 2212 || part.vertex().z() != 0.) && (!got_gen_vertex) ){
      //       genVertex = part.vertex();
      //       gen_vertex_z = genVertex.z();
      //       got_gen_vertex = 1;
      //       // cout << "Gen vertex z: " << genVertex.z() << endl;
      //     }
      //   }
      //   if (!got_gen_vertex){
      //     cout << "**********WARNING: Did not obtain non-zero GEN vertex from GEN particles" << endl;
      //   }
      // }

      // METfilters
      // bool passMETfilters=1;
      //Get trigger results relevant to MET filters

      edm::Handle<edm::TriggerResults> triggerBits;
      if(! event.isRealData() )
      event.getByToken( triggerPAT_, triggerBits );
      else
      event.getByToken( triggerRECO_, triggerBits );

      edm::Handle<edm::TriggerResults> triggerFLASHggMicroAOD;
      event.getByToken( triggerFLASHggMicroAOD_, triggerFLASHggMicroAOD );
      const edm::TriggerNames &triggerNames = event.triggerNames( *triggerBits );

      //check if passMETfilters
      std::vector<std::string> flagList {"Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_goodVertices","Flag_eeBadScFilter"};
      for( unsigned int i = 0; i < triggerNames.triggerNames().size(); i++ )
      {
        if(!triggerBits->accept(i))
        for(size_t j=0;j<flagList.size();j++)
        {
          if(flagList[j]==triggerNames.triggerName(i))
          {
            // passMETfilters=0;
            break;
          }
        }
      }

      std::vector<std::string> flashggFlagList {"flag_BadChargedCandidateFilter","flag_BadPFMuonFilter","flag_globalTightHalo2016Filter"};
      const edm::TriggerNames &flashggtriggerNames = event.triggerNames( *triggerFLASHggMicroAOD );
      for( unsigned int i = 0; i < flashggtriggerNames.triggerNames().size(); i++ )
      {
        if(!triggerFLASHggMicroAOD->accept(i))
        for(size_t j=0;j<flagList.size();j++)
        {
          if(flagList[j]==flashggtriggerNames.triggerName(i))
          {
            // passMETfilters=0;
            break;
          }
        }
      }


      if(doHHWWggTagCutFlowAnalysis_) Cut_Variables[0] = 1.0; // passed diphoton preselection (?) use value to check
      std::unique_ptr<vector<H4GTag> > H4Gtags( new vector<H4GTag> );

      edm::Ptr<reco::Vertex> vertex_diphoton;
      bool atLeastOneDiphoPass = false;
      for( unsigned int dpIndex = 0; dpIndex < diphotons->size(); dpIndex++ )
       {
         edm::Ptr<flashgg::DiPhotonCandidate> thisDPPtr = diphotons->ptrAt( dpIndex );
         vertex_diphoton = diphotons->ptrAt( dpIndex )->vtx();
         flashgg::DiPhotonCandidate * thisDPPointer = const_cast<flashgg::DiPhotonCandidate *>(thisDPPtr.get());
         atLeastOneDiphoPass |= idSelector_(*thisDPPointer, event);
       }
      
       cout << "atLeastOneDiphoPass: " << atLeastOneDiphoPass << endl;
      
       if (atLeastOneDiphoPass)
       {
         edm::Ptr<flashgg::DiPhotonCandidate> dipho;
         vector<edm::Ptr<flashgg::DiPhotonCandidate>> diphoVec;
         int n_pho = 0;
      
         vector<flashgg::Photon> phoVector;
         std::vector<edm::Ptr<flashgg::Photon>> phoPtrVector;
         cout << "# of diphotons " << diphotons->size() << endl;
      
         for( unsigned int diphoIndex = 0; diphoIndex < diphotons->size(); diphoIndex++ ) {
     
           dipho = diphotons->ptrAt( diphoIndex );
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

  cout << "# of photons " << phoVector.size() << endl;
  for( int phoIndex = 0; phoIndex < (int) phoVector.size(); phoIndex++ )
  {
    // edm::Ptr<flashgg::Photon> pho = photons->ptrAt(phoIndex);
    cout << "photon pt " << phoVector[phoIndex].pt() << endl;
    // phoPtrVector.push_back(pho);
  }
  //-- prepare a vector of ptr to photons, to be used for vertex selection
  for( int phoIndex = 0; phoIndex < (int) photons->size(); phoIndex++ )
  {
    edm::Ptr<flashgg::Photon> pho = photons->ptrAt(phoIndex);
    // cout << "photon pt " << pho->pt() << endl;
    phoPtrVector.push_back(pho);
  }
  for( int v = 0; v < (int) vertex->size(); v++ )
  {
    edm::Ptr<reco::Vertex> vtx = vertex->ptrAt( v );
    Vertices.push_back(vtx);
    if (vertex_diphoton->x() - vtx->x() == 0  && vertex_diphoton->y() - vtx->y() == 0 && vertex_diphoton->z() - vtx->z() == 0 )
    {
      hgg_index =  v;
    }
    if (fabs(genVertex.z() - vertex_diphoton->z()) > 1 )
    {
      if (fabs(genVertex.z() - vtx->z()) < 1)
      {
        slim_Vertices.push_back(vtx);
      }
      else{
        slim_Vertices.push_back(vertex->ptrAt( 0 ));
      }
    }
    else{
      slim_Vertices.push_back(vertex_diphoton);
    }
  }
  cout << hgg_index << randVtxIndexI << endl;
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



  if (photons->size() == 2)
  {
    cout << "2 photon category " << endl;
  }
  if (photons->size() == 3)
  {
    cout << "3 photon category " << endl;
  }
  if (photons->size() > 3)
  {
    cout << "4 photon category " << endl;
  }
}


      // if(diphotons->size() > 1){
      //   cout << "diphotons->size(): " << diphotons->size() << endl;
      //   for (unsigned int diphoIndex = 0; diphoIndex < diphotons->size(); diphoIndex++){
      //     edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( diphoIndex );
      //     cout << "dipho->pt(): " << dipho->pt() << endl;
      //     cout << "diphoton->sumPt(): " << dipho->sumPt() << endl;
      //   }
      // }

      if (diphotons->size() > 0){ // without systematics (?)
        for( unsigned int diphoIndex = 0; diphoIndex < 1; diphoIndex++ ) { // only look at highest pt dipho
          edm::Ptr<flashgg::DiPhotonCandidate> dipho = diphotons->ptrAt( diphoIndex ); // without systematic look (?)
          edm::Ptr<flashgg::DiPhotonMVAResult> mvares = mvaResults->ptrAt( diphoIndex );
          diphoton_vertex = dipho->vtx();

          // MVA selections
          // kinematic cuts on diphotons
          const flashgg::Photon* leadPho = dipho->leadingPhoton();
          const flashgg::Photon* subleadPho = dipho->subLeadingPhoton();

          diphoton_vertex = dipho->vtx();

          passMVAs = 0;
          passMVAs = checkPassMVAs(leadPho, subleadPho, diphoton_vertex);

          if(doHHWWggTagCutFlowAnalysis_){
            if(!passMVAs) Cut_Variables[1] = 0.0;
            else Cut_Variables[1] = 1.0; // passed photon MVAs (and all photon selections)
          }

          else{
            if(!passMVAs) continue; // Do not save event if leading and subleading photons don't pass MVA cuts
          }
          if(!passMVAs && !doHHWWggTagCutFlowAnalysis_) cout << "*********************************************problem" << endl;

          hasGoodElec = false;
          hasGoodMuons = false;

          // Electrons
          std::vector<edm::Ptr<Electron> > goodElectrons = selectStdElectrons( electrons->ptrs(), dipho, vertices->ptrs(), leptonPtThreshold_, electronEtaThresholds_,
          useElectronMVARecipe_,useElectronLooseID_,
          deltaRPhoElectronThreshold_,DeltaRTrkElec_,deltaMassElectronZThreshold_,
          rho_, event.isRealData() );

          // Muons
          std::vector<edm::Ptr<flashgg::Muon> > goodMuons = selectMuons( muons->ptrs(), dipho, vertex->ptrs(), muonEtaThreshold_, leptonPtThreshold_,
          muPFIsoSumRelThreshold_, deltaRMuonPhoThreshold_, deltaRMuonPhoThreshold_ );

          n_good_electrons = goodElectrons.size();
          n_good_muons = goodMuons.size();
          n_good_leptons = n_good_electrons + n_good_muons;
          hasGoodElec = ( goodElectrons.size() > 0 );
          hasGoodMuons = ( goodMuons.size() > 0 );

          // Jets
          unsigned int jetCollectionIndex = diphotons->at( diphoIndex ).jetCollectionIndex(); // not looping over systematics
          edm::Handle<edm::View<flashgg::Jet> > Jets_;


          event.getByToken( jetTokens_[jetCollectionIndex], Jets_); // testing


          std::vector<edm::Ptr<Jet> > tagJets;

          // Jet Selections
          for( unsigned int candIndex_outer = 0; candIndex_outer <  Jets_->size() ; candIndex_outer++ )
          {
            bool keepJet=true;
            edm::Ptr<flashgg::Jet> thejet = Jets_->ptrAt( candIndex_outer );


            if( fabs( thejet->eta() ) > jetEtaThreshold_ ) { keepJet=false; }

            if( thejet->pt() < jetPtThreshold_ ) { keepJet=false; }
            float dRPhoLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->leadingPhoton()->superCluster()->eta(), dipho->leadingPhoton()->superCluster()->phi() ) ;
            float dRPhoSubLeadJet = deltaR( thejet->eta(), thejet->phi(), dipho->subLeadingPhoton()->superCluster()->eta(),
            dipho->subLeadingPhoton()->superCluster()->phi() );

            if( dRPhoLeadJet < deltaRPhoLeadJet_ || dRPhoSubLeadJet < deltaRPhoSubLeadJet_ ) { keepJet=false; }
            if( hasGoodElec )
            for( unsigned int electronIndex = 0; electronIndex < goodElectrons.size(); electronIndex++ )
            {
              Ptr<flashgg::Electron> electron = goodElectrons[electronIndex];
              float dRJetElectron = deltaR( thejet->eta(), thejet->phi(), electron->eta(), electron->phi() ) ;
              if( dRJetElectron < deltaRJetMuonThreshold_ ) { keepJet=false; }
            }
            if( hasGoodMuons )
            for( unsigned int muonIndex = 0; muonIndex < goodMuons.size(); muonIndex++ )
            {
              Ptr<flashgg::Muon> muon = goodMuons[muonIndex];
              float dRJetMuon = deltaR( thejet->eta(), thejet->phi(), muon->eta(), muon->phi() ) ;
              if( dRJetMuon < deltaRJetMuonThreshold_ ) { keepJet=false; }
            }


            if(keepJet)
            tagJets.push_back( thejet );

          }

          // If jet collection has a jet suspected to be a b jet, don't save the event
          hasHighbTag = 0;
          for (unsigned int j = 0; j < tagJets.size(); j++){
            Ptr<flashgg::Jet> jet_ = tagJets[j];
            btagVal = jet_->bDiscriminator("mini_pfDeepFlavourJetTags:probb");
            if (btagVal > btagThresh_) hasHighbTag = 1;
          }

          // If doing cut flow analysis, don't continue
          if(doHHWWggTagCutFlowAnalysis_){
            if(hasHighbTag)
            Cut_Variables[2] = 0.0; // does not pass bveto
            else
            Cut_Variables[2] = 1.0; // passes bveto
          }
          else if(hasHighbTag) continue; // Skip event if it has at least one jet with a btag above threshold

          n_good_jets = tagJets.size();


          // MET
          if( METs->size() != 1 ) { std::cout << "WARNING - #MET is not 1" << std::endl;}
          Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );


          for (unsigned int i = 0; i < tagJets.size(); i++){
            auto tag_jet = tagJets[i];
            flashgg::Jet * thisJetPointer = const_cast<flashgg::Jet *>(tag_jet.get());
            tagJets_.push_back(*thisJetPointer);
          }


          if (doHHWWggTagCutFlowAnalysis_){
            if (n_good_leptons == 1) Cut_Variables[3] = 1.0; // exactly one good lepton
            else Cut_Variables[3] = 0.0;
            if (n_good_jets >= 2) Cut_Variables[4] = 1.0; // at least 2 good jets
            else Cut_Variables[4] = 0.0;
          }
          else
          if ((n_good_leptons != 1) || (n_good_jets <= 1)) continue;

          //-- Tag object
          if ( (n_good_leptons == 1) && (n_good_jets >= 2)){

            int catnum = 0;
            Ptr<flashgg::Jet> jet1 = tagJets[0];
            Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );

            if (n_good_electrons == 1){
              Ptr<flashgg::Electron> tag_electron = goodElectrons[0];
              catnum = 1;

              // if (n_good_jets == 2){
              Ptr<flashgg::Jet> jet2 = tagJets[1];
              H4GTag tag_obj;
              // HHWWggTag tag_obj_0;
              if (doHHWWggTagCutFlowAnalysis_){
                H4GTag tag_obj_(dipho, tag_electron, theMET, jet1, jet2, tagJets_, Cut_Variables); // electron, MET, jet1, jet2
                tag_obj = tag_obj_;
              }
              else{
                H4GTag tag_obj_(dipho, tag_electron, theMET, jet1, jet2); // diphoton, electron, MET, jet1, jet2
                tag_obj = tag_obj_;
              }

              tag_obj.setSystLabel( systLabel_);
              tag_obj.setDiPhotonIndex( diphoIndex );
              // tag_obj.setMVA( -0.9 );
              tag_obj.setCategoryNumber( catnum );
              tag_obj.includeWeights( *dipho );

              H4Gtags->push_back( tag_obj );


              if( ! event.isRealData() ) {
                H4Gtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
              }

            }

            if (n_good_muons == 1){
              Ptr<flashgg::Muon> tag_muon = goodMuons[0];
              catnum = 2;


              Ptr<flashgg::Jet> jet2 = tagJets[1];


              H4GTag tag_obj;
              // HHWWggTag tag_obj_0;
              if (doHHWWggTagCutFlowAnalysis_){
                H4GTag tag_obj_(dipho, tag_muon, theMET, jet1, jet2, tagJets_, Cut_Variables); // muon, MET, jet1, jet2
                tag_obj = tag_obj_;
              }
              else{
                H4GTag tag_obj_(dipho, tag_muon, theMET, jet1, jet2); // diphoton, muon, MET, jet1, jet2
                tag_obj = tag_obj_;
              }


              tag_obj.setSystLabel(systLabel_);

              tag_obj.setDiPhotonIndex( diphoIndex );
              // tag_obj.setMVA( -0.9 );
              tag_obj.setCategoryNumber( catnum ); // 2 for muon events
              tag_obj.includeWeights( *dipho );

              // tag_obj.setEventNumber(event.id().event() );
              // cout << "Pushing back tag object w/ muon" << endl;
              H4Gtags->push_back( tag_obj );


              if( ! event.isRealData() ) {
                H4Gtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
              }

            }

          } // if ( (n_good_leptons == 1) && (n_good_jets >= 2) )

          // Untagged category
          // Don't have semileptonic selections. Only results in further analysis if doHHWWggTagCutFlowAnalysis_ == True
          else {
            if(doHHWWggTagCutFlowAnalysis_){
              H4GTag tag_obj(dipho, tagJets_, Cut_Variables);
              // if (loopOverJets == 1) tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
              // else tag_obj.setSystLabel( inputJetsSuffixes_[jet_col_idx]);
              tag_obj.setSystLabel(systLabel_);

              tag_obj.setDiPhotonIndex( diphoIndex );
              // tag_obj.setMVA( -0.9 );
              tag_obj.setCategoryNumber( 0 );
              tag_obj.includeWeights( *dipho );

              // tag_obj.setEventNumber(event.id().event() );
              // cout << "Pushing back tag object w/ electron" << endl;
              H4Gtags->push_back( tag_obj );

              if( ! event.isRealData() ) {
                H4Gtags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
              }
            }
          }

        } // only look at highest pt dipho

      } // if at least 1 PS diphoton

      event.put( std::move( H4Gtags ) );
      event.put( std::move( truths ) );

      // cout << "[HHWWggTagProducer.cc] - End of HHWWggTagProducer::produce" << endl;

    } // H4GTagProducer::produce

  } // namespace flashgg

  typedef flashgg::H4GTagProducer FlashggH4GTagProducer;
  DEFINE_FWK_MODULE( FlashggH4GTagProducer );
