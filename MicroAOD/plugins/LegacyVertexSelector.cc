#include "flashgg/MicroAOD/interface/VertexSelectorBase.h"
#include "DataFormats/Common/interface/Handle.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TMVA/Reader.h"
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//--H4G includes
// #include "flashgg/DataFormats/interface/H4GCandidate.h"


namespace flashgg {

    struct Sorter {
        bool operator()( const std::pair<unsigned int, float> pair1, const std::pair<unsigned int, float> pair2 )
        {
            return ( pair1.second > pair2.second );
        };
    };

    class LegacyVertexSelector : public VertexSelectorBase
    {

    public:
        LegacyVertexSelector( const edm::ParameterSet & );
        ~LegacyVertexSelector();
        edm::Ptr<reco::Vertex> select( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &, const std::vector<edm::Ptr<reco::Vertex> > &,
                                       const VertexCandidateMap &vertexCandidateMap,
                                       const std::vector<edm::Ptr<reco::Conversion> > &,
                                       const std::vector<edm::Ptr<reco::Conversion> > &,
                                       const math::XYZPoint &,
                                       bool
                                     ) override;

        // for gamma+jets studies
        edm::Ptr<reco::Vertex> select( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<pat::Jet> &, const std::vector<edm::Ptr<reco::Vertex> > &,
                                       const VertexCandidateMap &vertexCandidateMap,
                                       const std::vector<edm::Ptr<reco::Conversion> > &,
                                       const std::vector<edm::Ptr<reco::Conversion> > &,
                                       const math::XYZPoint &,
                                       bool
                                       ) override;

        // for H4G studies
        std::vector<std::vector<float>> select_h2g( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &, const std::vector<edm::Ptr<reco::Vertex> > &,
                                      const VertexCandidateMap &vertexCandidateMap,
                                      const std::vector<edm::Ptr<reco::Conversion> > &,
                                      const std::vector<edm::Ptr<reco::Conversion> > &,
                                      const math::XYZPoint &,
                                      bool
                                    ) ;

       std::vector<std::vector<float>> select_h3g( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &,const edm::Ptr<flashgg::Photon> &, const std::vector<edm::Ptr<reco::Vertex> > &,
                                     const VertexCandidateMap &vertexCandidateMap,
                                     const std::vector<edm::Ptr<reco::Conversion> > &,
                                     const std::vector<edm::Ptr<reco::Conversion> > &,
                                     const math::XYZPoint &,
                                     bool
                                   ) ;

        std::vector<std::vector<float>> select_h4g( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &,const edm::Ptr<flashgg::Photon> &,const edm::Ptr<flashgg::Photon> &, const std::vector<edm::Ptr<reco::Vertex> > &,
                                      const VertexCandidateMap &vertexCandidateMap,
                                      const std::vector<edm::Ptr<reco::Conversion> > &,
                                      const std::vector<edm::Ptr<reco::Conversion> > &,
                                      const math::XYZPoint &,
                                      bool
                                    ) override;


        void writeInfoFromLastSelectionTo( flashgg::DiPhotonCandidate & ) override;
        void writeInfoFromLastSelectionTo( flashgg::PhotonJetCandidate & ) override;
        void writeInfoFromLastSelectionTo( flashgg::H4GCandidate & ) override;

        double vtxZFromConvOnly( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<reco::Conversion> &, const math::XYZPoint & ) const;
        double vtxZFromConvSuperCluster( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<reco::Conversion> &, const math::XYZPoint & ) const;
        double vtxZFromConv( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<reco::Conversion> &, const math::XYZPoint & ) const;
        double vtxDZFromConv( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<reco::Conversion> & ) const;

        double zFromConvPair( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &,
                              const int , const int, const int, const int,
                              const std::vector<edm::Ptr<reco::Conversion> > &,
                              const std::vector<edm::Ptr<reco::Conversion> > &,
                              const math::XYZPoint & ) const;

        double sZFromConvPair( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &,
                               const int , const int, const int, const int,
                               const std::vector<edm::Ptr<reco::Conversion> > &,
                               const std::vector<edm::Ptr<reco::Conversion> > & ) const;


        double zFromConvPair_h3g( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &,
                              const int , const int, const int, const int, const int, const int,
                              const std::vector<edm::Ptr<reco::Conversion> > &,
                              const std::vector<edm::Ptr<reco::Conversion> > &,
                              const math::XYZPoint & ) const;

        double sZFromConvPair_h3g( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &,
                                const int , const int, const int, const int, const int, const int,
                                const std::vector<edm::Ptr<reco::Conversion> > &,
                                const std::vector<edm::Ptr<reco::Conversion> > & ) const;

        double zFromConvPair_h4g( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &,
                              const int , const int, const int, const int, const int, const int, const int, const int,
                              const std::vector<edm::Ptr<reco::Conversion> > &,
                              const std::vector<edm::Ptr<reco::Conversion> > &,
                              const math::XYZPoint & ) const;

        double sZFromConvPair_h4g( const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &, const edm::Ptr<flashgg::Photon> &,
                                const int , const int, const int, const int, const int, const int, const int, const int,
                                const std::vector<edm::Ptr<reco::Conversion> > &,
                                const std::vector<edm::Ptr<reco::Conversion> > & ) const;

        std::vector<int> IndexMatchedConversion( const edm::Ptr<flashgg::Photon> &, const std::vector<edm::Ptr<reco::Conversion> > &,
                const std::vector<edm::Ptr<reco::Conversion> > &, bool ) const;

        void Initialize();

    private:

        edm::FileInPath vertexIdMVAweightfile_;
        edm::FileInPath vertexProbMVAweightfile_;

        edm::FileInPath vertexIdMVAweightfileH4G_;

        unsigned int nVtxSaveInfo;
        bool trackHighPurity;
        bool pureGeomConvMatching;
        double dRexclude;
        double sigma1Pix;
        double sigma1Tib;
        double sigma1Tob;
        double sigma1PixFwd;
        double sigma1Tid;
        double sigma1Tec;
        double sigma2Pix;
        double sigma2Tib;
        double sigma2Tob;
        double sigma2PixFwd;
        double sigma2Tid;
        double sigma2Tec;
        double singlelegsigma1Pix;
        double singlelegsigma1Tib;
        double singlelegsigma1Tob;
        double singlelegsigma1PixFwd;
        double singlelegsigma1Tid;
        double singlelegsigma1Tec;
        double singlelegsigma2Pix;
        double singlelegsigma2Tib;
        double singlelegsigma2Tob;
        double singlelegsigma2PixFwd;
        double singlelegsigma2Tid;
        double singlelegsigma2Tec;

    protected:

        TMVA::Reader *VertexIdMva_;
        bool initialized_;
        float logsumpt2_;
        float ptbal_;
        float ptasym_;
        float nConv_;
        float pull_conv_;

        TMVA::Reader *VertexIdMvaH4G_;

        float logsumpt2selected_;
        float ptbalselected_;
        float ptasymselected_;

        TMVA::Reader *VertexProbMva_;
        float dipho_pt_;
        float nVert_;
        float MVA0_;
        float MVA1_;
        float dZ1_;
        float MVA2_;
        float dZ2_;
        float vtxprobmva_;

        std::vector<float> vlogsumpt2_;
        std::vector<float> vptbal_;
        std::vector<float> vptasym_;
        std::vector<float> vpull_conv_;
        std::vector<float> vnConv_;
        std::vector<float> vmva_value_;
        std::vector<unsigned int> vmva_sortedindex_;
        std::vector<edm::Ptr<reco::Vertex> >  vVtxPtr_;

    };

    LegacyVertexSelector::LegacyVertexSelector( const edm::ParameterSet &iConfig ) :
        VertexSelectorBase( iConfig )
    {
        vertexIdMVAweightfile_ = iConfig.getParameter<edm::FileInPath>( "vertexIdMVAweightfile" );
        vertexProbMVAweightfile_ = iConfig.getParameter<edm::FileInPath>( "vertexProbMVAweightfile" );
        vertexIdMVAweightfileH4G_ = iConfig.getParameter<edm::FileInPath>( "vertexIdMVAweightfile" );

        nVtxSaveInfo          = iConfig.getUntrackedParameter<unsigned int>( "nVtxSaveInfo" );
        trackHighPurity       = iConfig.getParameter<bool>( "trackHighPurity" );
        pureGeomConvMatching  = iConfig.getParameter<bool>( "pureGeomConvMatching" );
        dRexclude             = iConfig.getParameter<double>( "dRexclude" );
        sigma1Pix             = iConfig.getParameter<double>( "sigma1Pix" );
        sigma1Tib             = iConfig.getParameter<double>( "sigma1Tib" );
        sigma1Tob             = iConfig.getParameter<double>( "sigma1Tob" );
        sigma1PixFwd          = iConfig.getParameter<double>( "sigma1PixFwd" );
        sigma1Tid             = iConfig.getParameter<double>( "sigma1Tid" );
        sigma1Tec             = iConfig.getParameter<double>( "sigma1Tec" );
        sigma2Pix             = iConfig.getParameter<double>( "sigma2Pix" );
        sigma2Tib             = iConfig.getParameter<double>( "sigma2Tib" );
        sigma2Tob             = iConfig.getParameter<double>( "sigma2Tob" );
        sigma2PixFwd          = iConfig.getParameter<double>( "sigma2PixFwd" );
        sigma2Tid             = iConfig.getParameter<double>( "sigma2Tid" );
        sigma2Tec             = iConfig.getParameter<double>( "sigma2Tec" );
        singlelegsigma1Pix    = iConfig.getParameter<double>( "singlelegsigma1Pix" );
        singlelegsigma1Tib    = iConfig.getParameter<double>( "singlelegsigma1Tib" );
        singlelegsigma1Tob    = iConfig.getParameter<double>( "singlelegsigma1Tob" );
        singlelegsigma1PixFwd = iConfig.getParameter<double>( "singlelegsigma1PixFwd" );
        singlelegsigma1Tid    = iConfig.getParameter<double>( "singlelegsigma1Tid" );
        singlelegsigma1Tec    = iConfig.getParameter<double>( "singlelegsigma1Tec" );
        singlelegsigma2Pix    = iConfig.getParameter<double>( "singlelegsigma2Pix" );
        singlelegsigma2Tib    = iConfig.getParameter<double>( "singlelegsigma2Tib" );
        singlelegsigma2Tob    = iConfig.getParameter<double>( "singlelegsigma2Tob" );
        singlelegsigma2PixFwd = iConfig.getParameter<double>( "singlelegsigma2PixFwd" );
        singlelegsigma2Tid    = iConfig.getParameter<double>( "singlelegsigma2Tid" );
        singlelegsigma2Tec    = iConfig.getParameter<double>( "singlelegsigma2Tec" );

        initialized_ = false;
    }

    void LegacyVertexSelector::Initialize()
    {
        VertexIdMva_ = new TMVA::Reader( "!Color:Silent" );
        VertexIdMva_->AddVariable( "ptasym", &ptasym_ );
        VertexIdMva_->AddVariable( "ptbal", &ptbal_ );
        VertexIdMva_->AddVariable( "logsumpt2", &logsumpt2_ );
        VertexIdMva_->AddVariable( "limPullToConv", &pull_conv_ );
        VertexIdMva_->AddVariable( "nConv", &nConv_ );

        VertexIdMva_->BookMVA( "BDT", vertexIdMVAweightfile_.fullPath() );

        VertexProbMva_ = new TMVA::Reader( "!Color:Silent" );

        VertexProbMva_->AddVariable( "pt", &dipho_pt_ );
        VertexProbMva_->AddVariable( "NVert", &nVert_ );
        VertexProbMva_->AddVariable( "MVA0", &MVA0_ );
        VertexProbMva_->AddVariable( "MVA1", &MVA1_ );
        VertexProbMva_->AddVariable( "DZ1", &dZ1_ );
        VertexProbMva_->AddVariable( "MVA2", &MVA2_ );
        VertexProbMva_->AddVariable( "DZ2", &dZ2_ );
        VertexProbMva_->AddVariable( "NConv", &nConv_ );

        VertexProbMva_->BookMVA( "BDT", vertexProbMVAweightfile_.fullPath() );

        VertexIdMvaH4G_ = new TMVA::Reader( "!Color:Silent" );
        // VertexIdMvaH4G_->AddVariable( "ptasym", &ptasym_ );
        // VertexIdMvaH4G_->AddVariable( "ptbal", &ptbal_ );
        VertexIdMvaH4G_->AddVariable( "logsumpt2", &logsumpt2_ );

        VertexIdMvaH4G_->BookMVA( "BDT", vertexIdMVAweightfileH4G_.fullPath() );

        initialized_ = true;
    }

    LegacyVertexSelector::~LegacyVertexSelector()
    {
        delete VertexIdMva_;
        delete VertexProbMva_;
        delete VertexIdMvaH4G_;
    }

    double LegacyVertexSelector::vtxZFromConvOnly( const edm::Ptr<flashgg::Photon> &pho, const edm::Ptr<reco:: Conversion> &conversion,
            const math::XYZPoint &beamSpot ) const
    {
        double dz = 0;
        if( conversion->nTracks() == 2 ) {
            double r = sqrt( conversion->refittedPairMomentum().perp2() );
            dz = ( conversion->conversionVertex().z() - beamSpot.z() )
                 -
                 ( ( conversion->conversionVertex().x() - beamSpot.x() ) * conversion->refittedPair4Momentum().x() + ( conversion->conversionVertex().y() - beamSpot.y() ) *
                   conversion->refittedPair4Momentum().y() ) / r * conversion->refittedPair4Momentum().z() / r;
        }
        if( conversion->nTracks() == 1 ) {
            double r = sqrt( conversion->tracksPin()[0].x() * conversion->tracksPin()[0].x() + conversion->tracksPin()[0].y() * conversion->tracksPin()[0].y() );
            dz = ( conversion->conversionVertex().z() - beamSpot.z() )
                 -
                 ( ( conversion->conversionVertex().x() - beamSpot.x() ) * conversion->tracksPin()[0].x() + ( conversion->conversionVertex().y() - beamSpot.y() ) *
                   conversion->tracksPin()[0].y() ) / r * conversion->tracksPin()[0].z() / r;
        }
        return dz + beamSpot.z();
    }

    double LegacyVertexSelector::vtxZFromConvSuperCluster( const edm::Ptr<flashgg::Photon> &pho, const edm::Ptr<reco:: Conversion> &conversion,
            const math::XYZPoint &beamSpot ) const
    {
        // get the z from conversion plus SuperCluster
        double deltaX1 =  pho->caloPosition().x() - conversion->conversionVertex().x();
        double deltaY1 =  pho->caloPosition().y() - conversion->conversionVertex().y();
        double deltaZ1 =  pho->caloPosition().z() - conversion->conversionVertex().z();
        double R1 = sqrt( deltaX1 * deltaX1 + deltaY1 * deltaY1 );
        double tantheta = R1 / deltaZ1;

        double deltaX2 = conversion->conversionVertex().x() - beamSpot.x();
        double deltaY2 = conversion->conversionVertex().y() - beamSpot.y();
        double R2 = sqrt( deltaX2 * deltaX2 + deltaY2 * deltaY2 );
        double deltaZ2 = R2 / tantheta;
        double higgsZ =  pho->caloPosition().z() - deltaZ1 - deltaZ2;
        return higgsZ;
    }

    double LegacyVertexSelector::vtxZFromConv( const edm::Ptr<flashgg::Photon> &pho, const edm::Ptr<reco::Conversion> &conversion,
            const math::XYZPoint &beamSpot ) const
    {
        double ReturnValue = 0;
        double perp = sqrt( conversion->conversionVertex().x() * conversion->conversionVertex().x() + conversion->conversionVertex().y() *
                            conversion->conversionVertex().y() );

        float nTracksConv = conversion->nTracks();

        if( nTracksConv == 2 ) {
            if( fabs( pho->superCluster()->eta() ) < 1.5 ) {
                if( perp <= 15.0 ) {

                    if( sigma1Pix < sigma2Pix )
                    { ReturnValue = vtxZFromConvOnly( pho, conversion, beamSpot ); }
                    else
                    { ReturnValue = vtxZFromConvSuperCluster( pho, conversion, beamSpot ); }
                } else if( perp > 15 && perp <= 60.0 ) {

                    if( sigma1Tib < sigma2Tib )
                    { ReturnValue = vtxZFromConvOnly( pho, conversion, beamSpot ); }
                    else
                    { ReturnValue = vtxZFromConvSuperCluster( pho, conversion, beamSpot ); }
                } else {

                    if( sigma1Tob < sigma2Tob )
                    { ReturnValue = vtxZFromConvOnly( pho, conversion, beamSpot ); }
                    else
                    { ReturnValue = vtxZFromConvSuperCluster( pho, conversion, beamSpot ); }
                }
            } else {
                if( fabs( conversion->conversionVertex().z() ) <= 50.0 ) {

                    if( sigma1PixFwd < sigma2PixFwd )
                    { ReturnValue = vtxZFromConvOnly( pho, conversion, beamSpot ); }
                    else
                    { ReturnValue = vtxZFromConvSuperCluster( pho, conversion, beamSpot ); }
                } else if( fabs( conversion->conversionVertex().z() ) > 50.0 && fabs( conversion->conversionVertex().z() ) <= 100.0 ) {
                    if( sigma1Tid < sigma2Tid )
                    { ReturnValue = vtxZFromConvOnly( pho, conversion, beamSpot ); }
                    else
                    { ReturnValue = vtxZFromConvSuperCluster( pho, conversion, beamSpot ); }
                } else {

                    if( sigma1Tec < sigma2Tec )
                    { ReturnValue = vtxZFromConvOnly( pho, conversion, beamSpot ); }
                    else
                    { ReturnValue = vtxZFromConvSuperCluster( pho, conversion, beamSpot ); }
                }
            }
        }
        if( nTracksConv == 1 ) {
            if( fabs( pho->superCluster()->eta() ) < 1.5 ) {
                if( perp <= 15.0 ) {

                    if( singlelegsigma1Pix < singlelegsigma2Pix )
                    { ReturnValue = vtxZFromConvOnly( pho, conversion, beamSpot ); }
                    else
                    { ReturnValue = vtxZFromConvSuperCluster( pho, conversion, beamSpot ); }
                } else if( perp > 15 && perp <= 60.0 ) {

                    if( singlelegsigma1Tib < singlelegsigma2Tib )
                    { ReturnValue = vtxZFromConvOnly( pho, conversion, beamSpot ); }
                    else
                    { ReturnValue = vtxZFromConvSuperCluster( pho, conversion, beamSpot ); }
                } else {

                    if( singlelegsigma1Tob < singlelegsigma2Tob )
                    { ReturnValue = vtxZFromConvOnly( pho, conversion, beamSpot ); }
                    else
                    { ReturnValue = vtxZFromConvSuperCluster( pho, conversion, beamSpot ); }
                }
            } else {
                if( fabs( conversion->conversionVertex().z() ) <= 50.0 ) {

                    if( singlelegsigma1PixFwd < singlelegsigma2PixFwd )
                    { ReturnValue = vtxZFromConvOnly( pho, conversion, beamSpot ); }
                    else
                    { ReturnValue = vtxZFromConvSuperCluster( pho, conversion, beamSpot ); }
                } else if( fabs( conversion->conversionVertex().z() ) > 50.0 && fabs( conversion->conversionVertex().z() ) <= 100.0 ) {

                    if( singlelegsigma1Tid < singlelegsigma2Tid )
                    { ReturnValue = vtxZFromConvOnly( pho, conversion, beamSpot ); }
                    else
                    { ReturnValue = vtxZFromConvSuperCluster( pho, conversion, beamSpot ); }
                } else {

                    if( singlelegsigma1Tec < singlelegsigma2Tec )
                    { ReturnValue = vtxZFromConvOnly( pho, conversion, beamSpot ); }
                    else
                    { ReturnValue = vtxZFromConvSuperCluster( pho, conversion, beamSpot ); }
                }
            }
        }
        return ReturnValue;
    }



    double LegacyVertexSelector::vtxDZFromConv( const edm::Ptr<flashgg::Photon> &pho, const edm::Ptr<reco::Conversion> &conversion ) const
    {
        double dz = -99999;
        double perp = sqrt( conversion->conversionVertex().x() * conversion->conversionVertex().x() + conversion->conversionVertex().y() *
                            conversion->conversionVertex().y() );
        if( conversion->nTracks() == 2 ) {
            if( fabs( pho->superCluster()->eta() ) < 1.5 ) { // barrel
                if( perp <= 15 ) {
                    dz = sigma1Pix;
                } else if( perp > 15 && perp <= 60 ) {
                    dz = sigma2Tib;
                } else {
                    dz = sigma2Tob;
                }
            } else { // endcap
                if( fabs( conversion->conversionVertex().z() ) <= 50 ) {
                    dz = sigma1PixFwd;
                } else if( fabs( conversion->conversionVertex().z() ) > 50 && fabs( conversion->conversionVertex().z() ) <= 100 ) {
                    dz = sigma1Tid;
                } else {
                    dz = sigma2Tec;
                }
            }
        } else if( conversion->nTracks() == 1 ) {
            if( fabs( pho->superCluster()->eta() ) < 1.5 ) { // barrel
                if( perp <= 15 ) {
                    dz = singlelegsigma1Pix;
                } else if( perp > 15 && perp <= 60 ) {
                    dz = singlelegsigma2Tib;
                } else {
                    dz = singlelegsigma2Tob;
                }

            } else { // endcap
                if( fabs( conversion->conversionVertex().z() ) <= 50 ) {
                    dz = singlelegsigma1PixFwd;
                } else if( fabs( conversion->conversionVertex().z() ) > 50 && fabs( conversion->conversionVertex().z() ) <= 100 ) {
                    dz = singlelegsigma1Tid;
                } else {
                    dz = singlelegsigma2Tec;
                }
            }
        }
        return dz;
    }

    double LegacyVertexSelector::zFromConvPair( const edm::Ptr<flashgg::Photon> &p1,
            const edm::Ptr<flashgg::Photon> &p2,
            const int index_conversionLead,
            const int index_conversionTrail,
            const int nConvLegs_LeadPhoton,
            const int nConvLegs_TrailPhoton,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg,
            const math::XYZPoint &beamSpot ) const
    {
        double zconv = 0;
        float z1 = 0;
        float sz1 = 0;
        float z2 = 0;
        float sz2 = 0;

        if( index_conversionLead != -1  && index_conversionTrail == -1 ) {
            if( nConvLegs_LeadPhoton == 2 ) { zconv = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot ); }
            if( nConvLegs_LeadPhoton == 1 ) { zconv = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );}
        }

        if( index_conversionLead == -1 && index_conversionTrail != -1 ) {
            if( nConvLegs_TrailPhoton == 2 ) { zconv = vtxZFromConv( p2, conversionsVector[index_conversionTrail], beamSpot ); }
            if( nConvLegs_TrailPhoton == 1 ) { zconv = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionTrail], beamSpot ); }
        }

        if( index_conversionLead != -1 && index_conversionTrail != -1 ) {

            if( nConvLegs_LeadPhoton == 2 && nConvLegs_TrailPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionTrail], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionTrail] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_TrailPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionTrail], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionTrail] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_TrailPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionTrail], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionTrail] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_TrailPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionTrail], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionTrail] );
            }
            if( sz1 != 0 && sz2 != 0 ) {
                zconv  = ( z1 / sz1 / sz1 + z2 / sz2 / sz2 ) / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 );
            }
        }
        return zconv;
    }

    double LegacyVertexSelector::zFromConvPair_h3g( const edm::Ptr<flashgg::Photon> &p1,
            const edm::Ptr<flashgg::Photon> &p2, const edm::Ptr<flashgg::Photon> &p3,
            const int index_conversionLead,
            const int index_conversionSecond,
            const int index_conversionThird,
            const int nConvLegs_LeadPhoton,
            const int nConvLegs_SecondPhoton,
            const int nConvLegs_ThirdPhoton,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg,
            const math::XYZPoint &beamSpot ) const
    {
        double zconv = 0;
        float z1 = 0;
        float sz1 = 0;
        float z2 = 0;
        float sz2 = 0;
        float z3 = 0;
        float sz3 = 0;

        if( index_conversionLead != -1  && index_conversionSecond == -1  && index_conversionThird == -1) {
            if( nConvLegs_LeadPhoton == 2 ) { zconv = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot ); }
            if( nConvLegs_LeadPhoton == 1 ) { zconv = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );}
        }

        if( index_conversionLead == -1 && index_conversionSecond != -1 && index_conversionThird == -1 ) {
            if( nConvLegs_SecondPhoton == 2 ) { zconv = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot ); }
            if( nConvLegs_SecondPhoton == 1 ) { zconv = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot ); }
        }

        if( index_conversionLead == -1 && index_conversionSecond == -1 && index_conversionThird != -1 ) {
            if( nConvLegs_ThirdPhoton == 2 ) { zconv = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot ); }
            if( nConvLegs_ThirdPhoton == 1 ) { zconv = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot ); }
        }

        if( index_conversionLead != -1 && index_conversionSecond != -1 && index_conversionThird == -1 ) {

            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            }
            if( sz1 != 0 && sz2 != 0 ) {
                zconv  = ( z1 / sz1 / sz1 + z2 / sz2 / sz2 ) / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 );
            }
        }

        if( index_conversionLead != -1 && index_conversionSecond == -1 && index_conversionThird != -1 ) {

            if( nConvLegs_LeadPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( sz1 != 0 && sz3 != 0 ) {
                zconv  = ( z1 / sz1 / sz1 + z3 / sz3 / sz3 ) / ( 1. / sz1 / sz1 + 1. / sz3 / sz3 );
            }
        }
        if( index_conversionLead == -1 && index_conversionSecond != -1 && index_conversionThird != -1 ) {

            if( nConvLegs_ThirdPhoton == 2 && nConvLegs_SecondPhoton == 2 ) {
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            }
            if( nConvLegs_ThirdPhoton == 1 && nConvLegs_SecondPhoton == 1 ) {
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            }
            if( nConvLegs_ThirdPhoton == 2 && nConvLegs_SecondPhoton == 1 ) {
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            }
            if( nConvLegs_ThirdPhoton == 1 && nConvLegs_SecondPhoton == 2 ) {
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            }
            if( sz3 != 0 && sz2 != 0 ) {
                zconv  = ( z3 / sz3 / sz3 + z2 / sz2 / sz2 ) / ( 1. / sz3 / sz3 + 1. / sz2 / sz2 );
            }
        }
        if( index_conversionLead != -1 && index_conversionSecond != -1 && index_conversionThird != -1 ) {
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
            z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
            z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
            z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
            z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
            z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
            z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
            z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
            z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
          }
          if( sz1 != 0 && sz2 != 0  && sz3 != 0 ) {
              zconv  = ( z1 / sz1 / sz1 + z2 / sz2 / sz2 + z3 / sz3 / sz3 ) / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 + 1. / sz3 / sz3 );
          }
        }
        return zconv;
    }

    double LegacyVertexSelector::zFromConvPair_h4g( const edm::Ptr<flashgg::Photon> &p1,
            const edm::Ptr<flashgg::Photon> &p2, const edm::Ptr<flashgg::Photon> &p3, const edm::Ptr<flashgg::Photon> &p4,
            const int index_conversionLead,
            const int index_conversionSecond,
            const int index_conversionThird,
            const int index_conversionFourth,
            const int nConvLegs_LeadPhoton,
            const int nConvLegs_SecondPhoton,
            const int nConvLegs_ThirdPhoton,
            const int nConvLegs_FourthPhoton,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg,
            const math::XYZPoint &beamSpot ) const
    {
        double zconv = 0;
        float z1 = 0;
        float sz1 = 0;
        float z2 = 0;
        float sz2 = 0;
        float z3 = 0;
        float sz3 = 0;
        float z4 = 0;
        float sz4 = 0;

        if( index_conversionLead != -1  && index_conversionSecond == -1 && index_conversionThird == -1 && index_conversionFourth == -1) {
            if( nConvLegs_LeadPhoton == 2 ) { zconv = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot ); }
            if( nConvLegs_LeadPhoton == 1 ) { zconv = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );}
        }

        if( index_conversionLead == -1 && index_conversionSecond != -1 && index_conversionThird == -1 && index_conversionFourth == -1 ) {

            if( nConvLegs_SecondPhoton == 2 ) { zconv = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot ); }
            if( nConvLegs_SecondPhoton == 1 ) { zconv = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot ); }
        }

        if( index_conversionLead == -1 && index_conversionSecond == -1 && index_conversionThird != -1 && index_conversionFourth == -1 ) {

            if( nConvLegs_ThirdPhoton == 2 ) { zconv = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot ); }
            if( nConvLegs_ThirdPhoton == 1 ) { zconv = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot ); }
        }

        if( index_conversionLead == -1 && index_conversionSecond == -1 && index_conversionThird == -1 && index_conversionFourth != -1 ) {

            if( nConvLegs_FourthPhoton == 2 ) { zconv = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot ); }
            if( nConvLegs_FourthPhoton == 1 ) { zconv = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot ); }
        }

        if( index_conversionLead != -1 && index_conversionSecond != -1 && index_conversionThird == -1 && index_conversionFourth == -1) {

            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            }
            if( sz1 != 0 && sz2 != 0 ) {
                zconv  = ( z1 / sz1 / sz1 + z2 / sz2 / sz2 ) / ( 1. / sz1 / sz1 + 1. / sz2 / sz2  );
            }
        }

        if( index_conversionLead == -1 && index_conversionSecond == -1 && index_conversionThird != -1 && index_conversionFourth != -1) {

            if( nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 2 ) {
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
            }
            if( nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 1 ) {
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
            }
            if( nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 1 ) {
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
            }
            if( nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 2 ) {
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
            }
            if( sz3 != 0 && sz4 != 0 ) {
                zconv  = ( z3 / sz3 / sz3 + z4 / sz4 / sz4 ) / ( 1. / sz3 / sz3 + 1. / sz4 / sz4  );
            }
        }

        if( index_conversionLead != -1 && index_conversionSecond == -1 && index_conversionThird == -1 && index_conversionFourth != -1) {

            if( nConvLegs_LeadPhoton == 2 && nConvLegs_FourthPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_FourthPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_FourthPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_FourthPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
            }
            if( sz1 != 0 && sz4 != 0 ) {
                zconv  = ( z1 / sz1 / sz1 + z4 / sz4 / sz4 ) / ( 1. / sz1 / sz1 + 1. / sz4 / sz4  );
            }
        }

        if( index_conversionLead == -1 && index_conversionSecond != -1 && index_conversionThird != -1 && index_conversionFourth == -1) {

            if( nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_SecondPhoton == 1 && nConvLegs_SecondPhoton == 2 ) {
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( sz2 != 0 && sz3!= 0 ) {
                zconv  = ( z2 / sz2 / sz2 + z3 / sz3 / sz3 ) / ( 1. / sz2 / sz2 + 1. / sz3 / sz3  );
            }
        }

        if( index_conversionLead != -1 && index_conversionSecond != -1 && index_conversionThird != -1 && index_conversionFourth == -1) {

            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( sz1 != 0 && sz2!= 0 && sz3 != 0) {
                zconv  = ( z1 / sz1 / sz1 + z2 / sz2 / sz2 + z3 / sz3 / sz3 ) / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 + 1. / sz3 / sz3  );
            }
        }

        if( index_conversionLead != -1 && index_conversionSecond != -1 && index_conversionThird == -1 && index_conversionFourth != -1) {

            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_FourthPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_FourthPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_FourthPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_FourthPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_FourthPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_FourthPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_FourthPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_FourthPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
            }
            if( sz1 != 0 && sz2!= 0 && sz4 != 0) {
                zconv  = ( z1 / sz1 / sz1 + z2 / sz2 / sz2 + z4 / sz4 / sz4 ) / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 + 1. / sz4 / sz4  );
            }
        }

        if( index_conversionLead != -1 && index_conversionSecond == -1 && index_conversionThird != -1 && index_conversionFourth != -1) {

            if( nConvLegs_LeadPhoton == 1 && nConvLegs_FourthPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_FourthPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_FourthPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_FourthPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_FourthPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_FourthPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_FourthPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_FourthPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                z1  = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( sz1 != 0 && sz4!= 0 && sz3 != 0) {
                zconv  = ( z1 / sz1 / sz1 + z4 / sz4 / sz4 + z3 / sz3 / sz3 ) / ( 1. / sz1 / sz1 + 1. / sz4 / sz4 + 1. / sz3 / sz3  );
            }
        }
        if( index_conversionLead == -1 && index_conversionSecond != -1 && index_conversionThird != -1 && index_conversionFourth != -1) {

            if( nConvLegs_FourthPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_FourthPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( nConvLegs_FourthPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_FourthPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_FourthPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                z4  = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( nConvLegs_FourthPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_FourthPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                z2  = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( nConvLegs_FourthPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                z4  = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                z2  = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                z3  = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( sz4 != 0 && sz2!= 0 && sz3 != 0) {
                zconv  = ( z4 / sz4 / sz4 + z2 / sz2 / sz2 + z3 / sz3 / sz3 ) / ( 1. / sz4 / sz4 + 1. / sz2 / sz2 + 1. / sz3 / sz3  );
            }
        }
        if( index_conversionLead != -1 && index_conversionSecond != -1 && index_conversionThird != -1 && index_conversionFourth != -1) {

          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 1) {

            z1 = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 2) {

            z1 = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 1) {

            z1 = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 1) {

            z1 = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 1) {

            z1 = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 2) {

            z1 = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 1) {

            z1 = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 1) {

            z1 = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 2) {

            z1 = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 2) {

            z1 = vtxZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 1) {

            z1 = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 2) {

            z1 = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 2) {

            z1 = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 2) {

            z1 = vtxZFromConv( p1, conversionsVector[index_conversionLead], beamSpot );
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            z2 = vtxZFromConv( p2, conversionsVector[index_conversionSecond], beamSpot );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            z3 = vtxZFromConv( p3, conversionsVector[index_conversionThird], beamSpot );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            z4 = vtxZFromConv( p4, conversionsVector[index_conversionFourth], beamSpot );
            sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
          }
          if( sz1 != 0 && sz2 != 0 && sz3 != 0 && sz4 != 0 ) {
                zconv  = ( z1 / sz1 / sz1 + z2 / sz2 / sz2 + z3 / sz3 / sz3 + z4 / sz4 / sz4 ) / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 + 1. / sz3 / sz3 + 1. / sz4 / sz4);
            }
        }
        return zconv;
    }

    double LegacyVertexSelector::sZFromConvPair( const edm::Ptr<flashgg::Photon> &p1,
            const edm::Ptr<flashgg::Photon> &p2,
            int index_conversionLead,
            int index_conversionTrail,
            const int nConvLegs_LeadPhoton,
            const int nConvLegs_TrailPhoton,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg ) const
    {
        double szconv = 0;
        float sz1 = 0;
        float sz2 = 0;

        if( index_conversionLead != -1  && index_conversionTrail == -1 ) {
            if( nConvLegs_LeadPhoton == 2 ) { szconv = vtxDZFromConv( p1, conversionsVector[index_conversionLead] ); }
            if( nConvLegs_LeadPhoton == 1 ) {  szconv = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] ); }
        }

        if( index_conversionLead == -1 && index_conversionTrail != -1 ) {
            if( nConvLegs_TrailPhoton == 2 ) { szconv = vtxDZFromConv( p2, conversionsVector[index_conversionTrail] ); }
            if( nConvLegs_TrailPhoton == 1 ) { szconv = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionTrail] ); }
        }

        if( index_conversionLead != -1 && index_conversionTrail != -1 ) {
            if( nConvLegs_LeadPhoton == 2 &&  nConvLegs_TrailPhoton == 2 ) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionTrail] );
            }
            if( nConvLegs_LeadPhoton == 1 &&  nConvLegs_TrailPhoton == 1 ) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionTrail] );
            }
            if( nConvLegs_LeadPhoton == 2 &&  nConvLegs_TrailPhoton == 1 ) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionTrail] );
            }
            if( nConvLegs_LeadPhoton == 1 &&  nConvLegs_TrailPhoton == 2 ) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionTrail] );
            }
            if( sz1 != 0 && sz2 != 0 ) {
                szconv = sqrt( 1. / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 ) ) ;
            }
        }
        return szconv;
    }

    double LegacyVertexSelector::sZFromConvPair_h3g( const edm::Ptr<flashgg::Photon> &p1,
            const edm::Ptr<flashgg::Photon> &p2, const edm::Ptr<flashgg::Photon> &p3,
            const int index_conversionLead,
            const int index_conversionSecond,
            const int index_conversionThird,
            const int nConvLegs_LeadPhoton,
            const int nConvLegs_SecondPhoton,
            const int nConvLegs_ThirdPhoton,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg) const
    {
        double szconv = 0;
        float sz1 = 0;
        float sz2 = 0;
        float sz3 = 0;

        if( index_conversionLead != -1  && index_conversionSecond == -1  && index_conversionThird == -1) {
            if( nConvLegs_LeadPhoton == 2 ) { szconv = vtxDZFromConv( p1, conversionsVector[index_conversionLead] ); }
            if( nConvLegs_LeadPhoton == 1 ) { szconv = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );}
        }

        if( index_conversionLead == -1 && index_conversionSecond != -1 && index_conversionThird == -1 ) {
            if( nConvLegs_SecondPhoton == 2 ) { szconv = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] ); }
            if( nConvLegs_SecondPhoton == 1 ) { szconv = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] ); }
        }

        if( index_conversionLead == -1 && index_conversionSecond == -1 && index_conversionThird != -1 ) {
            if( nConvLegs_ThirdPhoton == 2 ) { szconv = vtxDZFromConv( p3, conversionsVector[index_conversionThird] ); }
            if( nConvLegs_ThirdPhoton == 1 ) { szconv = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] ); }
        }

        if( index_conversionLead != -1 && index_conversionSecond != -1 && index_conversionThird == -1 ) {

            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 ) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 ) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 ) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 ) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            }
            if( sz1 != 0 && sz2 != 0 ) {
                szconv  = sqrt( 1. / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 ) ) ;
            }
        }

        if( index_conversionLead != -1 && index_conversionSecond == -1 && index_conversionThird != -1 ) {

            if( nConvLegs_LeadPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
            }
            if( nConvLegs_LeadPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
            }
            if( sz1 != 0 && sz3 != 0 ) {
                szconv  = sqrt( 1. / ( 1. / sz1 / sz1 + 1. / sz3 / sz3 ) ) ;
            }
        }
        if( index_conversionLead == -1 && index_conversionSecond != -1 && index_conversionThird != -1 ) {

            if( nConvLegs_ThirdPhoton == 2 && nConvLegs_SecondPhoton == 2 ) {
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            }
            if( nConvLegs_ThirdPhoton == 1 && nConvLegs_SecondPhoton == 1 ) {
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            }
            if( nConvLegs_ThirdPhoton == 2 && nConvLegs_SecondPhoton == 1 ) {
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            }
            if( nConvLegs_ThirdPhoton == 1 && nConvLegs_SecondPhoton == 2 ) {
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            }
            if( sz3 != 0 && sz2 != 0 ) {
                szconv  = sqrt( 1. / ( 1. / sz2 / sz2 + 1. / sz3 / sz3 ) ) ;
            }
        }
        if( index_conversionLead != -1 && index_conversionSecond != -1 && index_conversionThird != -1 ) {
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
            sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
          }
          if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
            sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
            sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
            sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
          }
          if( sz1 != 0 && sz2 != 0  && sz3 != 0 ) {
            szconv  = sqrt( 1. / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 + 1. / sz3 / sz3) ) ;
          }
        }
        return szconv;
    }


        double LegacyVertexSelector::sZFromConvPair_h4g( const edm::Ptr<flashgg::Photon> &p1,
                const edm::Ptr<flashgg::Photon> &p2, const edm::Ptr<flashgg::Photon> &p3, const edm::Ptr<flashgg::Photon> &p4,
                const int index_conversionLead,
                const int index_conversionSecond,
                const int index_conversionThird,
                const int index_conversionFourth,
                const int nConvLegs_LeadPhoton,
                const int nConvLegs_SecondPhoton,
                const int nConvLegs_ThirdPhoton,
                const int nConvLegs_FourthPhoton,
                const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,
                const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg
                 ) const
        {
            double szconv = 0;
            float sz1 = 0;
            float sz2 = 0;
            float sz3 = 0;
            float sz4 = 0;

            if( index_conversionLead != -1  && index_conversionSecond == -1 && index_conversionThird == -1 && index_conversionFourth == -1) {
                if( nConvLegs_LeadPhoton == 2 ) { szconv = vtxDZFromConv( p1, conversionsVector[index_conversionLead] ); }
                if( nConvLegs_LeadPhoton == 1 ) { szconv = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );}
            }

            if( index_conversionLead == -1 && index_conversionSecond != -1 && index_conversionThird == -1 && index_conversionFourth == -1 ) {
                if( nConvLegs_SecondPhoton == 2 ) { szconv = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] ); }
                if( nConvLegs_SecondPhoton == 1 ) { szconv = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] ); }
            }

            if( index_conversionLead == -1 && index_conversionSecond == -1 && index_conversionThird != -1 && index_conversionFourth == -1 ) {
                if( nConvLegs_ThirdPhoton == 2 ) { szconv = vtxDZFromConv( p2, conversionsVector[index_conversionThird] ); }
                if( nConvLegs_ThirdPhoton == 1 ) { szconv = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionThird] ); }
            }

            if( index_conversionLead == -1 && index_conversionSecond == -1 && index_conversionThird == -1 && index_conversionFourth != -1 ) {
                if( nConvLegs_FourthPhoton == 2 ) { szconv = vtxDZFromConv( p2, conversionsVector[index_conversionFourth] ); }
                if( nConvLegs_FourthPhoton == 1 ) { szconv = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionFourth] ); }
            }

            if( index_conversionLead != -1 && index_conversionSecond != -1 && index_conversionThird == -1 && index_conversionFourth == -1) {

                if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                }
                if( sz1 != 0 && sz2 != 0 ) {
                    szconv  = ( 1. / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 )  );
                }
            }

            if( index_conversionLead == -1 && index_conversionSecond == -1 && index_conversionThird != -1 && index_conversionFourth != -1) {

                if( nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 2 ) {
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                }
                if( nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 1 ) {
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                }
                if( nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 1 ) {
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                }
                if( nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 2 ) {
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                }
                if( sz3 != 0 && sz4 != 0 ) {
                    szconv  = ( 1. / ( 1. / sz3 / sz3 + 1. / sz4 / sz4 )  );
                }
            }

            if( index_conversionLead != -1 && index_conversionSecond == -1 && index_conversionThird == -1 && index_conversionFourth != -1) {

                if( nConvLegs_LeadPhoton == 2 && nConvLegs_FourthPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_FourthPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_FourthPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_FourthPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                }
                if( sz1 != 0 && sz4 != 0 ) {
                    szconv  = ( 1. / ( 1. / sz1 / sz1 + 1. / sz4 / sz4 )  );
                }
            }

            if( index_conversionLead == -1 && index_conversionSecond != -1 && index_conversionThird != -1 && index_conversionFourth == -1) {

                if( nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_SecondPhoton == 1 && nConvLegs_SecondPhoton == 2 ) {
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( sz2 != 0 && sz3!= 0 ) {
                    szconv  = ( 1. / ( 1. / sz2 / sz2 + 1. / sz3 / sz3 )  );
                }
            }

            if( index_conversionLead != -1 && index_conversionSecond != -1 && index_conversionThird != -1 && index_conversionFourth == -1) {

                if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( sz1 != 0 && sz2!= 0 && sz3 != 0) {
                    szconv  = ( 1. / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 + 1. / sz3 / sz3 )  );
                }
            }

            if( index_conversionLead != -1 && index_conversionSecond != -1 && index_conversionThird == -1 && index_conversionFourth != -1) {

                if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_FourthPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_FourthPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_FourthPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_FourthPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_FourthPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_FourthPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_FourthPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_FourthPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                }
                if( sz1 != 0 && sz2!= 0 && sz4 != 0) {
                    szconv  = ( 1. / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 + 1. / sz4 / sz4 )  );
                }
            }

            if( index_conversionLead != -1 && index_conversionSecond == -1 && index_conversionThird != -1 && index_conversionFourth != -1) {

                if( nConvLegs_LeadPhoton == 1 && nConvLegs_FourthPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_FourthPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_FourthPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_FourthPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 1 && nConvLegs_FourthPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_FourthPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_FourthPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( nConvLegs_LeadPhoton == 2 && nConvLegs_FourthPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                    sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( sz1 != 0 && sz4!= 0 && sz3 != 0) {
                    szconv  = ( 1. / ( 1. / sz1 / sz1 + 1. / sz4 / sz4 + 1. / sz3 / sz3 )  );
                }
            }
            if( index_conversionLead == -1 && index_conversionSecond != -1 && index_conversionThird != -1 && index_conversionFourth != -1) {

                if( nConvLegs_FourthPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_FourthPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( nConvLegs_FourthPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_FourthPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 ) {
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_FourthPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                    sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( nConvLegs_FourthPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 ) {
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                }
                if( nConvLegs_FourthPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 ) {
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                    sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( nConvLegs_FourthPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 ) {
                    sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
                    sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                    sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                }
                if( sz4 != 0 && sz2!= 0 && sz3 != 0) {
                    szconv  = ( 1. / ( 1. / sz4 / sz4 + 1. / sz2 / sz2 + 1. / sz3 / sz3 )  );
                }
            }

            if( index_conversionLead != -1 && index_conversionSecond != -1 && index_conversionThird != -1 && index_conversionFourth != -1) {

              if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 1) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 2) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 1) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 1) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 1) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 2) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 1) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 1) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 2) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 1 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 2) {
                sz1 = vtxDZFromConv( p1, conversionsVectorSingleLeg[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 1) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVectorSingleLeg[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 1 && nConvLegs_FourthPhoton == 2) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVectorSingleLeg[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 1 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 2) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVectorSingleLeg[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
              }
              if( nConvLegs_LeadPhoton == 2 && nConvLegs_SecondPhoton == 2 && nConvLegs_ThirdPhoton == 2 && nConvLegs_FourthPhoton == 2) {
                sz1 = vtxDZFromConv( p1, conversionsVector[index_conversionLead] );
                sz2 = vtxDZFromConv( p2, conversionsVector[index_conversionSecond] );
                sz3 = vtxDZFromConv( p3, conversionsVector[index_conversionThird] );
                sz4 = vtxDZFromConv( p4, conversionsVector[index_conversionFourth] );
              }
              if( sz1 != 0 && sz2 != 0 && sz3 != 0 && sz4 != 0 ) {
                    szconv  = sqrt( 1. / ( 1. / sz1 / sz1 + 1. / sz2 / sz2 + 1. / sz3 / sz3 + 1. / sz4 / sz4 ) );
                }
            }
            return szconv;
        }

    std::vector<int> LegacyVertexSelector::IndexMatchedConversion( const edm::Ptr<flashgg::Photon> &g,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,  const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg,
            bool useSingleLeg ) const
    {
        double mindR = 999;
        int nConvLegs = 0;
        bool doOneLeg = true;

        std::vector<int> result;

        if(!pureGeomConvMatching) assert( g->hasConversionTracks() );
        int selected_conversion_index = -1;

        if( (g->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching){

            for( unsigned int i = 0; i < conversionsVector.size(); i++ ) {
                edm::Ptr<reco::Conversion> conv = conversionsVector[i];
                if( conv->nTracks() == 2 ) {
                    if( !conv->isConverted() ) { continue; }
                    if( conv->refittedPair4Momentum().pt() < 10. ) { continue; }
                    if( TMath::Prob( conv->conversionVertex().chi2(), conv->conversionVertex().ndof() ) < 1e-6 ) { continue; }

                    TVector3 VtxtoSC;
                    VtxtoSC.SetXYZ( g->superCluster()->position().x() - conv->conversionVertex().x(),
                                    g->superCluster()->position().y() - conv->conversionVertex().y(),
                                    g->superCluster()->position().z() - conv->conversionVertex().z() );
                    TVector3 RefPairMo;
                    RefPairMo.SetXYZ( conv->refittedPairMomentum().x(), conv->refittedPairMomentum().y(), conv->refittedPairMomentum().z() );
                    double dR = 0;
                    dR = VtxtoSC.DeltaR( RefPairMo );
                    if( dR < mindR ) {
                        mindR = dR;
                        selected_conversion_index = i;
                    }
                }
            }
            if( mindR < 0.1 ) {
                result.push_back( selected_conversion_index );
                nConvLegs = 2;
                result.push_back( nConvLegs );
                doOneLeg = false;
            }
            if( doOneLeg && useSingleLeg ) {
                mindR = 999;
                for( unsigned int j = 0; j < conversionsVectorSingleLeg.size(); j++ ) {
                    edm::Ptr<reco::Conversion> conv = conversionsVectorSingleLeg[j];
                    if( conv->nTracks() == 1 ) {
                        TVector3 VtxtoSC;
                        VtxtoSC.SetXYZ( g->superCluster()->position().x() - conv->conversionVertex().x(),
                                        g->superCluster()->position().y() - conv->conversionVertex().y(),
                                        g->superCluster()->position().z() - conv->conversionVertex().z() );
                        TVector3 RefPairMo;
                        float oneLegTrack_X = conv->tracksPin()[0].x();
                        float oneLegTrack_Y = conv->tracksPin()[0].y();
                        float oneLegTrack_Z = conv->tracksPin()[0].z();

                        RefPairMo.SetXYZ( oneLegTrack_X, oneLegTrack_Y, oneLegTrack_Z );
                        double dR = 0;
                        dR = VtxtoSC.DeltaR( RefPairMo );
                        if( dR < mindR ) {
                            mindR = dR;
                            selected_conversion_index = j;
                        }
                    }
                }
                if( mindR < 0.1 ) {
                    result.push_back( selected_conversion_index );
                    nConvLegs = 1;
                    result.push_back( nConvLegs );
                }
            }
        }

        if( mindR < 0.1 )
            {return result;}
        else {
            result.push_back( -1 );
            result.push_back( -1 );
            return result;
        }
    }

    edm::Ptr<reco::Vertex> LegacyVertexSelector::select( const edm::Ptr<flashgg::Photon> &g1, const edm::Ptr<flashgg::Photon> &g2,
            const std::vector<edm::Ptr<reco::Vertex> > &vtxs,
            const VertexCandidateMap &vertexCandidateMap,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,
            const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg,
            const math::XYZPoint &beamSpot,
            bool useSingleLeg
            //						      const std::map<std::string,double> & param,
            //                                                      const float & beamsig
                                                       )
    {

        vlogsumpt2_.clear();
        vptbal_.clear();
        vptasym_.clear();
        vpull_conv_.clear();
        vnConv_.clear();
        vmva_value_.clear();
        vVtxPtr_.clear();
        vmva_sortedindex_.clear();

        std::vector<std::pair<unsigned int, float> > sorter;

        int IndexMatchedConversionLeadPhoton = -1;
        int IndexMatchedConversionTrailPhoton = -1;

        std::vector<int> vIndexMatchedConversionLeadPhoton;
        std::vector<int> vIndexMatchedConversionTrailPhoton;

        float nConvLegs_TrailPhoton = 0;
        float nConvLegs_LeadPhoton = 0;

        if( conversionsVector.size() > 0 || conversionsVectorSingleLeg.size() > 0 ) {
            if( (g1->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching) {
                vIndexMatchedConversionLeadPhoton = IndexMatchedConversion( g1, conversionsVector, conversionsVectorSingleLeg, useSingleLeg );
                IndexMatchedConversionLeadPhoton = vIndexMatchedConversionLeadPhoton[0];
                nConvLegs_LeadPhoton = vIndexMatchedConversionLeadPhoton[1];
            }
            if( (g2->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching) {
                vIndexMatchedConversionTrailPhoton = IndexMatchedConversion( g2, conversionsVector, conversionsVectorSingleLeg, useSingleLeg );
                IndexMatchedConversionTrailPhoton = vIndexMatchedConversionTrailPhoton[0];
                nConvLegs_TrailPhoton = vIndexMatchedConversionTrailPhoton[1];
            }
        }

        unsigned int vertex_index;

        unsigned int selected_vertex_index = 0;
        unsigned int second_selected_vertex_index = 0;
        unsigned int third_selected_vertex_index = 0;
        float max_mva_value = -999;
        float second_max_mva_value = -999;
        float third_max_mva_value = -999;

        if( !initialized_ ) {
            Initialize();
        }

        std::vector<float> vlogsumpt2;
        std::vector<float> vptbal;
        std::vector<float> vptasym;
        std::vector<float> vpull_conv;
        std::vector<float> vnConv;
        std::vector<float> vmva_value;
        std::vector<edm::Ptr<reco::Vertex> >  vVtxPtr;

        for( vertex_index = 0 ; vertex_index < vtxs.size() ; vertex_index++ ) {
            edm::Ptr<reco::Vertex> vtx = vtxs[vertex_index];

            TVector3 Photon1Dir;
            TVector3 Photon1Dir_uv;
            TVector3 Photon2Dir;
            TVector3 Photon2Dir_uv;
            TLorentzVector p14;
            TLorentzVector p24;
            Photon1Dir.SetXYZ( g1->superCluster()->position().x() - vtx->position().x(), g1->superCluster()->position().y() - vtx->position().y(),
                               g1->superCluster()->position().z() - vtx->position().z() );
            Photon2Dir.SetXYZ( g2->superCluster()->position().x() - vtx->position().x(), g2->superCluster()->position().y() - vtx->position().y(),
                               g2->superCluster()->position().z() - vtx->position().z() );
            Photon1Dir_uv = Photon1Dir.Unit() * g1->superCluster()->rawEnergy();
            Photon2Dir_uv = Photon2Dir.Unit() * g2->superCluster()->rawEnergy();
            p14.SetPxPyPzE( Photon1Dir_uv.x(), Photon1Dir_uv.y(), Photon1Dir_uv.z(), g1->superCluster()->rawEnergy() );
            p24.SetPxPyPzE( Photon2Dir_uv.x(), Photon2Dir_uv.y(), Photon2Dir_uv.z(), g2->superCluster()->rawEnergy() );

            TVector2 sumpt;
            double sumpt2_out = 0;
            double sumpt2_in = 0;
            double ptbal = 0;
            double ptasym = 0;

            sumpt.Set( 0., 0. );


            auto mapRange = std::equal_range( vertexCandidateMap.begin(), vertexCandidateMap.end(), vtx, flashgg::compare_with_vtx() );
            if( mapRange.first == mapRange.second ) { continue; }
            for( auto pair_iter = mapRange.first ; pair_iter != mapRange.second ; pair_iter++ ) {
                const edm::Ptr<pat::PackedCandidate> cand = pair_iter->second;
                TVector3 tk;
                TVector2 tkXY;
                double dr1 = 0;
                double dr2 = 0;
                tk.SetXYZ( cand->px(), cand->py(), cand->pz() );
                tkXY = tk.XYvector();
                dr1 = tk.DeltaR( p14.Vect() );
                dr2 = tk.DeltaR( p24.Vect() );
                bool isPure = cand->trackHighPurity();
                if( !isPure && trackHighPurity ) { continue; }

                if( dr1 < dRexclude || dr2 < dRexclude ) {
                    sumpt2_in += tkXY.Mod2();
                    continue;
                }
                sumpt += tkXY;
                sumpt2_out += tkXY.Mod2();
                ptbal -= tkXY * ( p14 + p24 ).Vect().XYvector().Unit();
            }

            ptasym = ( sumpt.Mod() - ( p14 + p24 ).Vect().XYvector().Mod() ) / ( sumpt.Mod() + ( p14 + p24 ).Vect().XYvector().Mod() );
            ptasym_ = ptasym;

            float nConv = 0;
            if( IndexMatchedConversionLeadPhoton != -1 ) { ++nConv; }
            if( IndexMatchedConversionTrailPhoton != -1 ) { ++nConv; }

            float zconv = 0;
            float szconv = 0;
            float pull_conv = 999;

            if( nConv != 0 ) {
                zconv = zFromConvPair( g1, g2, IndexMatchedConversionLeadPhoton, IndexMatchedConversionTrailPhoton, nConvLegs_LeadPhoton, nConvLegs_TrailPhoton,
                                       conversionsVector, conversionsVectorSingleLeg, beamSpot );

                szconv = sZFromConvPair( g1, g2, IndexMatchedConversionLeadPhoton, IndexMatchedConversionTrailPhoton, nConvLegs_LeadPhoton, nConvLegs_TrailPhoton,
                                         conversionsVector, conversionsVectorSingleLeg );
                if( szconv != 0 ) {
                    pull_conv = fabs( vtx->position().z() - zconv ) / szconv;
                } else {
                    pull_conv = 10.;
                }
            }

            if( pull_conv > 10. ) { pull_conv = 10.; }

            logsumpt2_ = log( sumpt2_in + sumpt2_out );
            ptbal_ = ptbal;
            pull_conv_ = pull_conv;
            nConv_ = nConv;
            float mva_value = VertexIdMva_->EvaluateMVA( "BDT" );

            vlogsumpt2.push_back( logsumpt2_ );
            vptbal.push_back( ptbal_ );
            vptasym.push_back( ptasym_ );
            vpull_conv.push_back( pull_conv_ );
            vnConv.push_back( nConv_ );
            vmva_value.push_back( mva_value );
            vVtxPtr.push_back( vtx );

            std::pair<unsigned int, float>pairToSort = std::make_pair( vmva_value.size() - 1, mva_value );
            sorter.push_back( pairToSort );

            if( mva_value > max_mva_value ) {
                max_mva_value = mva_value;
                selected_vertex_index = vertex_index;
                logsumpt2selected_ = logsumpt2_;
                ptbalselected_ = ptbal_;
                ptasymselected_ = ptasym_;
            }
        }

        std::sort( sorter.begin(), sorter.end(), Sorter() );

        if( sorter.size() > 1 ) {
            second_max_mva_value = sorter[1].second;
            second_selected_vertex_index = sorter[1].first;
        } else{
            second_max_mva_value = -2;
        }
        if( sorter.size() > 2 ) {
            third_max_mva_value = sorter[2].second;
            third_selected_vertex_index = sorter[2].first;
        }else{
            third_max_mva_value = -2 ;
        }

        for( unsigned int jj = 0; jj < sorter.size(); jj++ ) {

            if( vlogsumpt2_.size() < nVtxSaveInfo ) {

                vmva_sortedindex_.push_back( sorter[jj].first );
                vlogsumpt2_.push_back( vlogsumpt2[sorter[jj].first] );
                vptbal_.push_back( vptbal[sorter[jj].first] );
                vptasym_.push_back( vptasym[sorter[jj].first] );
                vpull_conv_.push_back( vpull_conv[sorter[jj].first] );
                vnConv_.push_back( vnConv[sorter[jj].first] );
                vmva_value_.push_back( vmva_value[sorter[jj].first] );
                vVtxPtr_.push_back( vVtxPtr[sorter[jj].first] );
            }

        }

        dipho_pt_ = ( g1->p4() + g2->p4() ).pt();
        nVert_    = vtxs.size();
        MVA0_     = max_mva_value;
        MVA1_     = second_max_mva_value;
        dZ1_      = vtxs[selected_vertex_index]->position().z() - vtxs[second_selected_vertex_index]->position().z();
        MVA2_     = third_max_mva_value;
        dZ2_      = vtxs[selected_vertex_index]->position().z() - vtxs[third_selected_vertex_index]->position().z();

        if( sorter.size() < 2 ) dZ1_=100;
        if( sorter.size() < 3 ) dZ2_=100;


        vtxprobmva_ = VertexProbMva_->EvaluateMVA( "BDT" );

        return vtxs[selected_vertex_index];
    }


    edm::Ptr<reco::Vertex> LegacyVertexSelector::select( const edm::Ptr<flashgg::Photon> &g1, const edm::Ptr<pat::Jet> &g2,
                                                         const std::vector<edm::Ptr<reco::Vertex> > &vtxs,
                                                         const VertexCandidateMap &vertexCandidateMap,
                                                         const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,
                                                         const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg,
                                                         const math::XYZPoint &beamSpot,
                                                         bool useSingleLeg
                                                         )
    {

        vlogsumpt2_.clear();
        vptbal_.clear();
        vptasym_.clear();
        vpull_conv_.clear();
        vnConv_.clear();
        vmva_value_.clear();
        vVtxPtr_.clear();
        vmva_sortedindex_.clear();

        std::vector<std::pair<unsigned int, float> > sorter;

        int IndexMatchedConversionLeadPhoton = -1;
        int IndexMatchedConversionTrailPhoton = -1;

        std::vector<int> vIndexMatchedConversionLeadPhoton;
        std::vector<int> vIndexMatchedConversionTrailPhoton;

        float nConvLegs_LeadPhoton = 0;
        float nConvLegs_TrailPhoton = 0;



        if( conversionsVector.size() > 0 ) {
            if( (g1->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching ) {
                vIndexMatchedConversionLeadPhoton = IndexMatchedConversion( g1, conversionsVector, conversionsVectorSingleLeg, useSingleLeg );
                IndexMatchedConversionLeadPhoton = vIndexMatchedConversionLeadPhoton[0];
                nConvLegs_LeadPhoton = vIndexMatchedConversionLeadPhoton[1];
            }
        }

        unsigned int vertex_index;

        unsigned int selected_vertex_index = 0;
        unsigned int second_selected_vertex_index = 0;
        unsigned int third_selected_vertex_index = 0;
        float max_mva_value = -999;
        float second_max_mva_value = -999;
        float third_max_mva_value = -999;

        if( !initialized_ ) {
            Initialize();
        }

        std::vector<float> vlogsumpt2;
        std::vector<float> vptbal;
        std::vector<float> vptasym;
        std::vector<float> vpull_conv;
        std::vector<float> vnConv;
        std::vector<float> vmva_value;
        std::vector<edm::Ptr<reco::Vertex> >  vVtxPtr;

        for( vertex_index = 0 ; vertex_index < vtxs.size() ; vertex_index++ ) {
            edm::Ptr<reco::Vertex> vtx = vtxs[vertex_index];

            TVector3 Photon1Dir;
            TVector3 Photon1Dir_uv;
            TLorentzVector p14;
            Photon1Dir.SetXYZ( g1->superCluster()->position().x() - vtx->position().x(), g1->superCluster()->position().y() - vtx->position().y(),
                               g1->superCluster()->position().z() - vtx->position().z() );
            Photon1Dir_uv = Photon1Dir.Unit() * g1->superCluster()->rawEnergy();
            p14.SetPxPyPzE( Photon1Dir_uv.x(), Photon1Dir_uv.y(), Photon1Dir_uv.z(), g1->superCluster()->rawEnergy() );

            TLorentzVector p24;
            p24.SetPxPyPzE(g2->px(), g2->py(), g2->pz(), g2->energy());

            TVector2 sumpt;
            double sumpt2_out = 0;
            double sumpt2_in = 0;
            double ptbal = 0;
            double ptasym = 0;

            sumpt.Set( 0., 0. );

            auto mapRange = std::equal_range( vertexCandidateMap.begin(), vertexCandidateMap.end(), vtx, flashgg::compare_with_vtx() );
            if( mapRange.first == mapRange.second ) { continue; }

            for( auto pair_iter = mapRange.first ; pair_iter != mapRange.second ; pair_iter++ ) {
                const edm::Ptr<pat::PackedCandidate> cand = pair_iter->second;

                TVector3 tk;
                TVector2 tkXY;
                double dr1 = 0;
                double dr2 = 0;
                tk.SetXYZ( cand->px(), cand->py(), cand->pz() );
                tkXY = tk.XYvector();
                dr1 = tk.DeltaR( p14.Vect() );
                dr2 = tk.DeltaR( p24.Vect() );

                // gamma+jet: skip tracks around the jet direction
                if ( dr2 < 0.4 ) { continue; }

                bool isPure = cand->trackHighPurity();
                if( !isPure && trackHighPurity ) { continue; }

                if( dr1 < dRexclude || dr2 < dRexclude ) {
                    sumpt2_in += tkXY.Mod2();
                    continue;
                }
                sumpt += tkXY;
                sumpt2_out += tkXY.Mod2();
                ptbal -= tkXY * ( p14 + p24 ).Vect().XYvector().Unit();
            }

            ptasym = ( sumpt.Mod() - ( p14 + p24 ).Vect().XYvector().Mod() ) / ( sumpt.Mod() + ( p14 + p24 ).Vect().XYvector().Mod() );
            ptasym_ = ptasym;

            float nConv = 0;
            if( IndexMatchedConversionLeadPhoton != -1 ) { ++nConv; }

            float zconv = 0;
            float szconv = 0;
            float pull_conv = 999;

            if( nConv != 0 ) {

                const edm::Ptr<flashgg::Photon> &dummyg2 = g1;

                zconv = zFromConvPair( g1, dummyg2, IndexMatchedConversionLeadPhoton, IndexMatchedConversionTrailPhoton, nConvLegs_LeadPhoton, nConvLegs_TrailPhoton,
                                       conversionsVector, conversionsVectorSingleLeg, beamSpot );

                szconv = sZFromConvPair( g1, dummyg2, IndexMatchedConversionLeadPhoton, IndexMatchedConversionTrailPhoton, nConvLegs_LeadPhoton, nConvLegs_TrailPhoton,
                                         conversionsVector, conversionsVectorSingleLeg );
                if( szconv != 0 ) {
                    pull_conv = fabs( vtx->position().z() - zconv ) / szconv;
                } else {
                    pull_conv = 10.;
                }
            }

            if( pull_conv > 10. ) { pull_conv = 10.; }


            logsumpt2_ = log( sumpt2_in + sumpt2_out );
            ptbal_ = ptbal;
            pull_conv_ = pull_conv;
            nConv_ = nConv;
            float mva_value = VertexIdMva_->EvaluateMVA( "BDT" );

            vlogsumpt2.push_back( logsumpt2_ );
            vptbal.push_back( ptbal_ );
            vptasym.push_back( ptasym_ );
            vpull_conv.push_back( pull_conv_ );
            vnConv.push_back( nConv_ );
            vmva_value.push_back( mva_value );
            vVtxPtr.push_back( vtx );

            std::pair<unsigned int, float>pairToSort = std::make_pair( vmva_value.size() - 1, mva_value );
            sorter.push_back( pairToSort );

            if( mva_value > max_mva_value ) {
                max_mva_value = mva_value;
                selected_vertex_index = vertex_index;
                logsumpt2selected_ = logsumpt2_;
                ptbalselected_ = ptbal_;
                ptasymselected_ = ptasym_;
            }
        }

        std::sort( sorter.begin(), sorter.end(), Sorter() );

        if( sorter.size() > 1 ) {
            second_max_mva_value = sorter[1].second;
            second_selected_vertex_index = sorter[1].first;
        }
        if( sorter.size() > 2 ) {
            third_max_mva_value = sorter[2].second;
            third_selected_vertex_index = sorter[2].first;
        }

        for( unsigned int jj = 0; jj < sorter.size(); jj++ ) {

            if( vlogsumpt2_.size() < nVtxSaveInfo ) {

                vmva_sortedindex_.push_back( sorter[jj].first );
                vlogsumpt2_.push_back( vlogsumpt2[sorter[jj].first] );
                vptbal_.push_back( vptbal[sorter[jj].first] );
                vptasym_.push_back( vptasym[sorter[jj].first] );
                vpull_conv_.push_back( vpull_conv[sorter[jj].first] );
                vnConv_.push_back( vnConv[sorter[jj].first] );
                vmva_value_.push_back( vmva_value[sorter[jj].first] );
                vVtxPtr_.push_back( vVtxPtr[sorter[jj].first] );
            }

        }
        dipho_pt_ = ( g1->p4() + g2->p4() ).pt();
        nVert_    = vtxs.size();
        MVA0_     = max_mva_value;
        MVA1_     = second_max_mva_value;
        dZ1_      = vtxs[selected_vertex_index]->position().z() - vtxs[second_selected_vertex_index]->position().z();
        MVA2_     = third_max_mva_value;
        dZ2_      = vtxs[selected_vertex_index]->position().z() - vtxs[third_selected_vertex_index]->position().z();

        vtxprobmva_ = VertexProbMva_->EvaluateMVA( "BDT" );

        return vtxs[selected_vertex_index];

    }

    std::vector<std::vector<float>> LegacyVertexSelector::select_h2g( const edm::Ptr<flashgg::Photon> &g1, const edm::Ptr<flashgg::Photon> &g2,
                                                         const std::vector<edm::Ptr<reco::Vertex> > &vtxs,
                                                         const VertexCandidateMap &vertexCandidateMap,
                                                         const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,
                                                         const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg,
                                                         const math::XYZPoint &beamSpot,
                                                         bool useSingleLeg
                                                         )
    {

      vlogsumpt2_.clear();
      vptbal_.clear();
      vptasym_.clear();
      vpull_conv_.clear();
      vnConv_.clear();

      int IndexMatchedConversionLeadPhoton = -1;
      int IndexMatchedConversionSecondPhoton = -1;

      std::vector<int> vIndexMatchedConversionLeadPhoton;
      std::vector<int> vIndexMatchedConversionSecondPhoton;

      float nConvLegs_LeadPhoton = 0;
      float nConvLegs_SecondPhoton = 0;

      if( conversionsVector.size() > 0 || conversionsVectorSingleLeg.size() > 0 ) {
      if( (g1->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching) {
           vIndexMatchedConversionLeadPhoton = IndexMatchedConversion( g1, conversionsVector, conversionsVectorSingleLeg, useSingleLeg );
           IndexMatchedConversionLeadPhoton = vIndexMatchedConversionLeadPhoton[0];
           nConvLegs_LeadPhoton = vIndexMatchedConversionLeadPhoton[1];
          }
      if( (g2->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching) {
           vIndexMatchedConversionSecondPhoton = IndexMatchedConversion( g2, conversionsVector, conversionsVectorSingleLeg, useSingleLeg );
           IndexMatchedConversionSecondPhoton = vIndexMatchedConversionSecondPhoton[0];
           nConvLegs_SecondPhoton = vIndexMatchedConversionSecondPhoton[1];
          }
      }

      std::vector<std::vector<float>> megaVec;

      unsigned int vertex_index;
      for( vertex_index = 0 ; vertex_index < vtxs.size() ; vertex_index++ ) {
            edm::Ptr<reco::Vertex> vtx = vtxs[vertex_index];

            TVector3 Photon1Dir;
            TVector3 Photon1Dir_uv;
            TVector3 Photon2Dir;
            TVector3 Photon2Dir_uv;

            TLorentzVector p14;
            TLorentzVector p24;

            Photon1Dir.SetXYZ( g1->superCluster()->position().x() - vtx->position().x(), g1->superCluster()->position().y() - vtx->position().y(),
                               g1->superCluster()->position().z() - vtx->position().z() );
            Photon2Dir.SetXYZ( g2->superCluster()->position().x() - vtx->position().x(), g2->superCluster()->position().y() - vtx->position().y(),
                               g2->superCluster()->position().z() - vtx->position().z() );

            Photon1Dir_uv = Photon1Dir.Unit() * g1->superCluster()->rawEnergy();
            Photon2Dir_uv = Photon2Dir.Unit() * g2->superCluster()->rawEnergy();
            p14.SetPxPyPzE( Photon1Dir_uv.x(), Photon1Dir_uv.y(), Photon1Dir_uv.z(), g1->superCluster()->rawEnergy() );
            p24.SetPxPyPzE( Photon2Dir_uv.x(), Photon2Dir_uv.y(), Photon2Dir_uv.z(), g2->superCluster()->rawEnergy() );

            TVector2 sumpt;
            double sumpt2_out = 0;
            double sumpt2_in = 0;
            double ptbal = 0;
            double ptasym = 0;

            sumpt.Set( 0., 0. );

            auto mapRange = std::equal_range( vertexCandidateMap.begin(), vertexCandidateMap.end(), vtx, flashgg::compare_with_vtx() );
            if( mapRange.first == mapRange.second ) { continue; }
            for( auto pair_iter = mapRange.first ; pair_iter != mapRange.second ; pair_iter++ ) {
                const edm::Ptr<pat::PackedCandidate> cand = pair_iter->second;
                TVector3 tk;
                TVector2 tkXY;
                double dr1 = 0.;
                double dr2 = 0.;

                tk.SetXYZ( cand->px(), cand->py(), cand->pz() );
                tkXY = tk.XYvector();

                dr1 = tk.DeltaR(p14.Vect());
                dr2 = tk.DeltaR(p24.Vect());

                bool isPure = cand->trackHighPurity();
                if( !isPure && trackHighPurity ) { continue; }
                if( dr1 < dRexclude || dr2 < dRexclude ) {
                    sumpt2_in += tkXY.Mod2();
                    continue;
                }
                sumpt += tkXY;
                sumpt2_out += tkXY.Mod2();
                ptbal -= tkXY * ( p14 + p24 ).Vect().XYvector().Unit();
              }

              ptasym = ( sumpt.Mod() - ( p14 + p24  ).Vect().XYvector().Mod() ) / ( sumpt.Mod() + ( p14 + p24  ).Vect().XYvector().Mod() );
              ptasym_ = ptasym;

              float nConv = 0;
              if( IndexMatchedConversionLeadPhoton != -1 ) { ++nConv; }
              if( IndexMatchedConversionSecondPhoton != -1 ) { ++nConv; }

              float zconv = 0;
              float szconv = 0;
              float pull_conv = 999;

              if( nConv != 0 ) {
              zconv = zFromConvPair( g1, g2, IndexMatchedConversionLeadPhoton, IndexMatchedConversionSecondPhoton, nConvLegs_LeadPhoton, nConvLegs_SecondPhoton,
                                      conversionsVector, conversionsVectorSingleLeg, beamSpot );

              szconv = sZFromConvPair( g1, g2, IndexMatchedConversionLeadPhoton, IndexMatchedConversionSecondPhoton, nConvLegs_LeadPhoton, nConvLegs_SecondPhoton,
                                       conversionsVector, conversionsVectorSingleLeg );
              if( szconv != 0 ) {
                  pull_conv = fabs( vtx->position().z() - zconv ) / szconv;
              } else {
                  pull_conv = 10.;
                }
              }

              if( pull_conv > 10. ) { pull_conv = 10.; }

              logsumpt2_ = log( sumpt2_in + sumpt2_out );
              ptbal_ = ptbal;
              pull_conv_ = pull_conv;
              nConv_ = nConv;

              vlogsumpt2_.push_back( logsumpt2_ );
              vptasym_.push_back(ptasym_);
              vptbal_.push_back(ptbal_);
              vpull_conv_.push_back(pull_conv_);
              vnConv_.push_back(nConv_);
      }
      megaVec.push_back(vlogsumpt2_);
      megaVec.push_back(vptasym_);
      megaVec.push_back(vptbal_);
      megaVec.push_back(vpull_conv_);
      megaVec.push_back(vnConv_);
      return megaVec;
    }

    std::vector<std::vector<float>> LegacyVertexSelector::select_h3g( const edm::Ptr<flashgg::Photon> &g1, const edm::Ptr<flashgg::Photon> &g2, const edm::Ptr<flashgg::Photon> &g3,
                                                         const std::vector<edm::Ptr<reco::Vertex> > &vtxs,
                                                         const VertexCandidateMap &vertexCandidateMap,
                                                         const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,
                                                         const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg,
                                                         const math::XYZPoint &beamSpot,
                                                         bool useSingleLeg
                                                         )
    {

      vlogsumpt2_.clear();
      vptbal_.clear();
      vptasym_.clear();
      vpull_conv_.clear();
      vnConv_.clear();

      int IndexMatchedConversionLeadPhoton = -1;
      int IndexMatchedConversionSecondPhoton = -1;
      int IndexMatchedConversionThirdPhoton = -1;

      std::vector<int> vIndexMatchedConversionLeadPhoton;
      std::vector<int> vIndexMatchedConversionSecondPhoton;
      std::vector<int> vIndexMatchedConversionThirdPhoton;

      float nConvLegs_LeadPhoton = 0;
      float nConvLegs_SecondPhoton = 0;
      float nConvLegs_ThirdPhoton = 0;

      if( conversionsVector.size() > 0 || conversionsVectorSingleLeg.size() > 0 ) {
      if( (g1->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching) {
           vIndexMatchedConversionLeadPhoton = IndexMatchedConversion( g1, conversionsVector, conversionsVectorSingleLeg, useSingleLeg );
           IndexMatchedConversionLeadPhoton = vIndexMatchedConversionLeadPhoton[0];
           nConvLegs_LeadPhoton = vIndexMatchedConversionLeadPhoton[1];
          }
      if( (g2->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching) {
           vIndexMatchedConversionSecondPhoton = IndexMatchedConversion( g2, conversionsVector, conversionsVectorSingleLeg, useSingleLeg );
           IndexMatchedConversionSecondPhoton = vIndexMatchedConversionSecondPhoton[0];
           nConvLegs_SecondPhoton = vIndexMatchedConversionSecondPhoton[1];
          }
     if( (g3->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching) {
           vIndexMatchedConversionThirdPhoton = IndexMatchedConversion( g3, conversionsVector, conversionsVectorSingleLeg, useSingleLeg );
           IndexMatchedConversionThirdPhoton = vIndexMatchedConversionThirdPhoton[0];
           nConvLegs_ThirdPhoton = vIndexMatchedConversionThirdPhoton[1];
         }
      }


      std::vector<std::vector<float>> megaVec;

      unsigned int vertex_index;
      for( vertex_index = 0 ; vertex_index < vtxs.size() ; vertex_index++ ) {
            edm::Ptr<reco::Vertex> vtx = vtxs[vertex_index];

            TVector3 Photon1Dir;
            TVector3 Photon1Dir_uv;
            TVector3 Photon2Dir;
            TVector3 Photon2Dir_uv;
            TVector3 Photon3Dir;
            TVector3 Photon3Dir_uv;

            TLorentzVector p14;
            TLorentzVector p24;
            TLorentzVector p34;
            Photon1Dir.SetXYZ( g1->superCluster()->position().x() - vtx->position().x(), g1->superCluster()->position().y() - vtx->position().y(),
                               g1->superCluster()->position().z() - vtx->position().z() );
            Photon2Dir.SetXYZ( g2->superCluster()->position().x() - vtx->position().x(), g2->superCluster()->position().y() - vtx->position().y(),
                               g2->superCluster()->position().z() - vtx->position().z() );
            Photon3Dir.SetXYZ( g3->superCluster()->position().x() - vtx->position().x(), g3->superCluster()->position().y() - vtx->position().y(),
                               g3->superCluster()->position().z() - vtx->position().z() );

            Photon1Dir_uv = Photon1Dir.Unit() * g1->superCluster()->rawEnergy();
            Photon2Dir_uv = Photon2Dir.Unit() * g2->superCluster()->rawEnergy();
            Photon3Dir_uv = Photon3Dir.Unit() * g3->superCluster()->rawEnergy();
            p14.SetPxPyPzE( Photon1Dir_uv.x(), Photon1Dir_uv.y(), Photon1Dir_uv.z(), g1->superCluster()->rawEnergy() );
            p24.SetPxPyPzE( Photon2Dir_uv.x(), Photon2Dir_uv.y(), Photon2Dir_uv.z(), g2->superCluster()->rawEnergy() );
            p34.SetPxPyPzE( Photon3Dir_uv.x(), Photon3Dir_uv.y(), Photon3Dir_uv.z(), g3->superCluster()->rawEnergy() );

            TVector2 sumpt;
            double sumpt2_out = 0;
            double sumpt2_in = 0;
            double ptbal = 0;
            double ptasym = 0;

            sumpt.Set( 0., 0. );

            auto mapRange = std::equal_range( vertexCandidateMap.begin(), vertexCandidateMap.end(), vtx, flashgg::compare_with_vtx() );
            if( mapRange.first == mapRange.second ) { continue; }
            for( auto pair_iter = mapRange.first ; pair_iter != mapRange.second ; pair_iter++ ) {
                const edm::Ptr<pat::PackedCandidate> cand = pair_iter->second;
                TVector3 tk;
                TVector2 tkXY;
                double dr1 = 0.;
                double dr2 = 0.;
                double dr3 = 0.;

                tk.SetXYZ( cand->px(), cand->py(), cand->pz() );
                tkXY = tk.XYvector();

                dr1 = tk.DeltaR(p14.Vect());
                dr2 = tk.DeltaR(p24.Vect());
                dr3 = tk.DeltaR(p34.Vect());

                bool isPure = cand->trackHighPurity();
                if( !isPure && trackHighPurity ) { continue; }
                if( dr1 < dRexclude || dr2 < dRexclude || dr3 < dRexclude ) {
                    sumpt2_in += tkXY.Mod2();
                    continue;
                }
                sumpt += tkXY;
                sumpt2_out += tkXY.Mod2();
                ptbal -= tkXY * ( p14 + p24 + p34 ).Vect().XYvector().Unit();
              }

              ptasym = ( sumpt.Mod() - ( p14 + p24 + p34  ).Vect().XYvector().Mod() ) / ( sumpt.Mod() + ( p14 + p24 +p34  ).Vect().XYvector().Mod() );
              ptasym_ = ptasym;

              float nConv = 0;
              if( IndexMatchedConversionLeadPhoton != -1 ) { ++nConv; }
              if( IndexMatchedConversionSecondPhoton != -1 ) { ++nConv; }
              if( IndexMatchedConversionThirdPhoton != -1 ) { ++nConv; }

              float zconv = 0;
              float szconv = 0;
              float pull_conv = 999;

              if( nConv != 0 ) {
              zconv = zFromConvPair_h3g( g1, g2, g3, IndexMatchedConversionLeadPhoton, IndexMatchedConversionSecondPhoton,IndexMatchedConversionThirdPhoton, nConvLegs_LeadPhoton, nConvLegs_SecondPhoton,
                                     nConvLegs_ThirdPhoton, conversionsVector, conversionsVectorSingleLeg, beamSpot );

              szconv = sZFromConvPair_h3g( g1, g2, g3, IndexMatchedConversionLeadPhoton, IndexMatchedConversionSecondPhoton, IndexMatchedConversionThirdPhoton, nConvLegs_LeadPhoton, nConvLegs_SecondPhoton, nConvLegs_ThirdPhoton,
                                       conversionsVector, conversionsVectorSingleLeg );
              if( szconv != 0 ) {
                  pull_conv = fabs( vtx->position().z() - zconv ) / szconv;
              } else {
                  pull_conv = 10.;
                }
              }

              if( pull_conv > 10. ) { pull_conv = 10.; }

              logsumpt2_ = log( sumpt2_in + sumpt2_out );
              ptbal_ = ptbal;
              pull_conv_ = pull_conv;
              nConv_ = nConv;

              vlogsumpt2_.push_back( logsumpt2_ );
              vptasym_.push_back(ptasym_);
              vptbal_.push_back(ptbal_);
              vpull_conv_.push_back(pull_conv_);
              vnConv_.push_back(nConv_);

      }
      megaVec.push_back(vlogsumpt2_);
      megaVec.push_back(vptasym_);
      megaVec.push_back(vptbal_);
      megaVec.push_back(vpull_conv_);
      megaVec.push_back(vnConv_);
      return megaVec;
    }

    std::vector<std::vector<float>> LegacyVertexSelector::select_h4g( const edm::Ptr<flashgg::Photon> &g1, const edm::Ptr<flashgg::Photon> &g2, const edm::Ptr<flashgg::Photon> &g3, const edm::Ptr<flashgg::Photon> &g4,
                                                         const std::vector<edm::Ptr<reco::Vertex> > &vtxs,
                                                         const VertexCandidateMap &vertexCandidateMap,
                                                         const std::vector<edm::Ptr<reco::Conversion> > &conversionsVector,
                                                         const std::vector<edm::Ptr<reco::Conversion> > &conversionsVectorSingleLeg,
                                                         const math::XYZPoint &beamSpot,
                                                         bool useSingleLeg
                                                         )
    {

      vlogsumpt2_.clear();
      vptbal_.clear();
      vptasym_.clear();
      vpull_conv_.clear();
      vnConv_.clear();

      int IndexMatchedConversionLeadPhoton = -1;
      int IndexMatchedConversionSecondPhoton = -1;
      int IndexMatchedConversionThirdPhoton = -1;
      int IndexMatchedConversionFourthPhoton = -1;

      std::vector<int> vIndexMatchedConversionLeadPhoton;
      std::vector<int> vIndexMatchedConversionSecondPhoton;
      std::vector<int> vIndexMatchedConversionThirdPhoton;
      std::vector<int> vIndexMatchedConversionFourthPhoton;

      float nConvLegs_LeadPhoton = 0;
      float nConvLegs_SecondPhoton = 0;
      float nConvLegs_ThirdPhoton = 0;
      float nConvLegs_FourthPhoton = 0;

      if( conversionsVector.size() > 0 || conversionsVectorSingleLeg.size() > 0 ) {
      if( (g1->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching) {
           vIndexMatchedConversionLeadPhoton = IndexMatchedConversion( g1, conversionsVector, conversionsVectorSingleLeg, useSingleLeg );
           IndexMatchedConversionLeadPhoton = vIndexMatchedConversionLeadPhoton[0];
           nConvLegs_LeadPhoton = vIndexMatchedConversionLeadPhoton[1];
          }
      if( (g2->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching) {
           vIndexMatchedConversionSecondPhoton = IndexMatchedConversion( g2, conversionsVector, conversionsVectorSingleLeg, useSingleLeg );
           IndexMatchedConversionSecondPhoton = vIndexMatchedConversionSecondPhoton[0];
           nConvLegs_SecondPhoton = vIndexMatchedConversionSecondPhoton[1];
          }
     if( (g3->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching) {
           vIndexMatchedConversionThirdPhoton = IndexMatchedConversion( g3, conversionsVector, conversionsVectorSingleLeg, useSingleLeg );
           IndexMatchedConversionThirdPhoton = vIndexMatchedConversionThirdPhoton[0];
           nConvLegs_ThirdPhoton = vIndexMatchedConversionThirdPhoton[1];
         }
     if( (g4->hasConversionTracks() && !pureGeomConvMatching) || pureGeomConvMatching) {
         vIndexMatchedConversionFourthPhoton = IndexMatchedConversion( g4, conversionsVector, conversionsVectorSingleLeg, useSingleLeg );
         IndexMatchedConversionFourthPhoton = vIndexMatchedConversionFourthPhoton[0];
         nConvLegs_FourthPhoton = vIndexMatchedConversionFourthPhoton[1];
         }
      }


      std::vector<std::vector<float>> megaVec;

      unsigned int vertex_index;
      // if ( !initialized_ ){
      //   Initialize();
      // }
      for( vertex_index = 0 ; vertex_index < vtxs.size() ; vertex_index++ ) {
            edm::Ptr<reco::Vertex> vtx = vtxs[vertex_index];

            TVector3 Photon1Dir;
            TVector3 Photon1Dir_uv;
            TVector3 Photon2Dir;
            TVector3 Photon2Dir_uv;
            TVector3 Photon3Dir;
            TVector3 Photon3Dir_uv;
            TVector3 Photon4Dir;
            TVector3 Photon4Dir_uv;
            TLorentzVector p14;
            TLorentzVector p24;
            TLorentzVector p34;
            TLorentzVector p44;
            Photon1Dir.SetXYZ( g1->superCluster()->position().x() - vtx->position().x(), g1->superCluster()->position().y() - vtx->position().y(),
                               g1->superCluster()->position().z() - vtx->position().z() );
            Photon2Dir.SetXYZ( g2->superCluster()->position().x() - vtx->position().x(), g2->superCluster()->position().y() - vtx->position().y(),
                               g2->superCluster()->position().z() - vtx->position().z() );
            Photon3Dir.SetXYZ( g3->superCluster()->position().x() - vtx->position().x(), g3->superCluster()->position().y() - vtx->position().y(),
                               g3->superCluster()->position().z() - vtx->position().z() );
            Photon4Dir.SetXYZ( g4->superCluster()->position().x() - vtx->position().x(), g4->superCluster()->position().y() - vtx->position().y(),
                               g4->superCluster()->position().z() - vtx->position().z() );
            Photon1Dir_uv = Photon1Dir.Unit() * g1->superCluster()->rawEnergy();
            Photon2Dir_uv = Photon2Dir.Unit() * g2->superCluster()->rawEnergy();
            Photon3Dir_uv = Photon3Dir.Unit() * g3->superCluster()->rawEnergy();
            Photon4Dir_uv = Photon4Dir.Unit() * g4->superCluster()->rawEnergy();
            p14.SetPxPyPzE( Photon1Dir_uv.x(), Photon1Dir_uv.y(), Photon1Dir_uv.z(), g1->superCluster()->rawEnergy() );
            p24.SetPxPyPzE( Photon2Dir_uv.x(), Photon2Dir_uv.y(), Photon2Dir_uv.z(), g2->superCluster()->rawEnergy() );
            p34.SetPxPyPzE( Photon3Dir_uv.x(), Photon3Dir_uv.y(), Photon3Dir_uv.z(), g3->superCluster()->rawEnergy() );
            p44.SetPxPyPzE( Photon4Dir_uv.x(), Photon4Dir_uv.y(), Photon4Dir_uv.z(), g4->superCluster()->rawEnergy() );

            TVector2 sumpt;
            double sumpt2_out = 0;
            double sumpt2_in = 0;
            double ptbal = 0;
            double ptasym = 0;

            sumpt.Set( 0., 0. );

            auto mapRange = std::equal_range( vertexCandidateMap.begin(), vertexCandidateMap.end(), vtx, flashgg::compare_with_vtx() );
            if( mapRange.first == mapRange.second ) { continue; }
            for( auto pair_iter = mapRange.first ; pair_iter != mapRange.second ; pair_iter++ ) {
                const edm::Ptr<pat::PackedCandidate> cand = pair_iter->second;
                TVector3 tk;
                TVector2 tkXY;
                double dr1 = 0.;
                double dr2 = 0.;
                double dr3 = 0.;
                double dr4 = 0.;

                tk.SetXYZ( cand->px(), cand->py(), cand->pz() );
                tkXY = tk.XYvector();

                dr1 = tk.DeltaR(p14.Vect());
                dr2 = tk.DeltaR(p24.Vect());
                dr3 = tk.DeltaR(p34.Vect());
                dr4 = tk.DeltaR(p44.Vect());

                bool isPure = cand->trackHighPurity();
                if( !isPure && trackHighPurity ) { continue; }
                if( dr1 < dRexclude || dr2 < dRexclude || dr3 < dRexclude || dr4 < dRexclude ) {
                    sumpt2_in += tkXY.Mod2();
                    continue;
                }
                sumpt += tkXY;
                sumpt2_out += tkXY.Mod2();
                ptbal -= tkXY * ( p14 + p24 + p34 + p44).Vect().XYvector().Unit();
              }

              ptasym = ( sumpt.Mod() - ( p14 + p24 + p34 + p44 ).Vect().XYvector().Mod() ) / ( sumpt.Mod() + ( p14 + p24 +p34 + p44 ).Vect().XYvector().Mod() );
              ptasym_ = ptasym;

              float nConv = 0;
              if( IndexMatchedConversionLeadPhoton != -1 ) { ++nConv; }
              if( IndexMatchedConversionSecondPhoton != -1 ) { ++nConv; }
              if( IndexMatchedConversionThirdPhoton != -1 ) { ++nConv; }
              if( IndexMatchedConversionFourthPhoton != -1 ) { ++nConv; }

              float zconv = 0;
              float szconv = 0;
              float pull_conv = 999;
              if( nConv != 0 ) {
              zconv = zFromConvPair_h4g( g1, g2, g3, g4, IndexMatchedConversionLeadPhoton, IndexMatchedConversionSecondPhoton,IndexMatchedConversionThirdPhoton, IndexMatchedConversionFourthPhoton, nConvLegs_LeadPhoton, nConvLegs_SecondPhoton,
                                     nConvLegs_ThirdPhoton, nConvLegs_FourthPhoton, conversionsVector, conversionsVectorSingleLeg, beamSpot );
              szconv = sZFromConvPair_h4g( g1, g2, g3, g4, IndexMatchedConversionLeadPhoton, IndexMatchedConversionSecondPhoton, IndexMatchedConversionThirdPhoton, IndexMatchedConversionFourthPhoton, nConvLegs_LeadPhoton, nConvLegs_SecondPhoton, nConvLegs_ThirdPhoton, nConvLegs_FourthPhoton,
                                       conversionsVector, conversionsVectorSingleLeg );

              if( szconv != 0 ) {
                  pull_conv = fabs( vtx->position().z() - zconv ) / szconv;
              } else {
                  pull_conv = 10.;
                }
              }

              if( pull_conv > 10. ) { pull_conv = 10.; }

              logsumpt2_ = log( sumpt2_in + sumpt2_out );
              ptbal_ = ptbal;
              pull_conv_ = pull_conv;
              nConv_ = nConv;
              // float mva_value = VertexIdMva_->EvaluateMVA("BDT");

              // cout << "mva value = " << mva_value;

              vlogsumpt2_.push_back( logsumpt2_ );
              vptasym_.push_back(ptasym_);
              vptbal_.push_back(ptbal_);
              vpull_conv_.push_back(pull_conv_);
              vnConv_.push_back(nConv_);

      }
      megaVec.push_back(vlogsumpt2_);
      megaVec.push_back(vptasym_);
      megaVec.push_back(vptbal_);
      megaVec.push_back(vpull_conv_);
      megaVec.push_back(vnConv_);
      return megaVec;
    }

    void LegacyVertexSelector::writeInfoFromLastSelectionTo( flashgg::DiPhotonCandidate &dipho )
    {

        dipho.setLogSumPt2( logsumpt2selected_ );
        dipho.setPtBal( ptbalselected_ );
        dipho.setPtAsym( ptasymselected_ );

        dipho.setNConv( nConv_ );
        dipho.setPullConv( pull_conv_ );
        dipho.setNVert( nVert_ );
        dipho.setMVA0( MVA0_ );
        dipho.setMVA1( MVA1_ );
        dipho.setMVA2( MVA2_ );
        dipho.setDZ1( dZ1_ );
        dipho.setDZ2( dZ2_ );

        dipho.setVNConv( vnConv_ );
        dipho.setVPullConv( vpull_conv_ );
        dipho.setVPtBal( vptbal_ );
        dipho.setVPtAsym( vptasym_ );
        dipho.setVLogSumPt2( vlogsumpt2_ );
        dipho.setVMVA( vmva_value_ );
        dipho.setVVtxPtr( vVtxPtr_ );
        dipho.setVMVASortedIndex( vmva_sortedindex_ );

        dipho.setVtxProbMVA( vtxprobmva_ );
    }


    void LegacyVertexSelector::writeInfoFromLastSelectionTo( flashgg::PhotonJetCandidate &phojet )
    {

        phojet.setLogSumPt2( logsumpt2selected_ );
        phojet.setPtBal( ptbalselected_ );
        phojet.setPtAsym( ptasymselected_ );

        phojet.setNConv( nConv_ );
        phojet.setPullConv( pull_conv_ );
        phojet.setNVert( nVert_ );
        phojet.setMVA0( MVA0_ );
        phojet.setMVA1( MVA1_ );
        phojet.setMVA2( MVA2_ );
        phojet.setDZ1( dZ1_ );
        phojet.setDZ2( dZ2_ );

        phojet.setVNConv( vnConv_ );
        phojet.setVPullConv( vpull_conv_ );
        phojet.setVPtBal( vptbal_ );
        phojet.setVPtAsym( vptasym_ );
        phojet.setVLogSumPt2( vlogsumpt2_ );
        phojet.setVMVA( vmva_value_ );
        phojet.setVVtxPtr( vVtxPtr_ );
        phojet.setVMVASortedIndex( vmva_sortedindex_ );

        phojet.setVtxProbMVA( vtxprobmva_ );
    }

    void LegacyVertexSelector::writeInfoFromLastSelectionTo( flashgg::H4GCandidate &h4g )
    {
      // h4g.setLogSumPt2( logsumpt2_);
      // h4g.setPtBal( ptbal_ );
      // h4g.setPtAsym( ptasym_);
      //
      // h4g.setVLogSumPt2( vlogsumpt2_ );
      // h4g.setVPtBal(vptbal_ );
      // h4g.setVPtAsym( vptasym_ );
      //
      // h4g.setLogSumPt2_Zero( vlogsumpt2_[0]);
    }

} // namespace flashgg


DEFINE_EDM_PLUGIN( FlashggVertexSelectorFactory,
                   flashgg::LegacyVertexSelector,
                   "FlashggLegacyVertexSelector" );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
