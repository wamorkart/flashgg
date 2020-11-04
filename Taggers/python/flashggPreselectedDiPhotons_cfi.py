import FWCore.ParameterSet.Config as cms

phoEffArea=cms.PSet( var=cms.string("abs(superCluster.eta)"), bins=cms.vdouble(0.,0.9,1.5,2,2.2,3), vals=cms.vdouble(0.16544,0.16544,0.13212,0.13212,0.13212) )

rediscoveryHLTvariables = cms.vstring(
    "pfPhoIso03",
    "trkSumPtHollowConeDR03",
    "full5x5_sigmaIetaIeta",
    "full5x5_r9",
    "passElectronVeto"
    )

#cuts to mimic category trigger cuts
rediscoveryHLTcutsV1 = cms.VPSet(
    cms.PSet(cut=cms.string("isEB && full5x5_r9>0.85"), ##EB high R9
             selection = cms.VPSet(
            cms.PSet(max=cms.string("10000000000"),
                     rhocorr=phoEffArea,
                     ),
            cms.PSet(max=cms.string("10000000000")),
            cms.PSet(max=cms.string("10000000000")),
            cms.PSet(min=cms.string("-10000000000")),
            cms.PSet(min=cms.string("-10000000000"))
            ),
             ),

    cms.PSet(cut=cms.string("isEE && full5x5_r9>0.90"),  ##EE high R9
             selection = cms.VPSet(
            cms.PSet(max=cms.string("10000000000"),
                     rhocorr=phoEffArea,
                     ),
            cms.PSet(max=cms.string("10000000000")),
            cms.PSet(max=cms.string("10000000000")),
            cms.PSet(min=cms.string("-10000000000")),
            cms.PSet(min=cms.string("-10000000000"))
            ),
             ),
    cms.PSet(cut=cms.string("isEB && full5x5_r9<=0.85"),  #EB low R9
             selection = cms.VPSet(
            cms.PSet(max=cms.string("10000000000"),
                     rhocorr=phoEffArea,
                     ),
            cms.PSet(max=cms.string("10000000000")),
            cms.PSet(max=cms.string("10000000000")),
            cms.PSet(min=cms.string("-10000000000")),
            cms.PSet(min=cms.string("-10000000000"))
            ),
             ),
    cms.PSet(cut=cms.string("isEE && full5x5_r9<=0.90"),  ##EE low R9
             selection = cms.VPSet(
            cms.PSet(max=cms.string("10000000000"),
                     rhocorr=phoEffArea,
                     ),
            cms.PSet(max=cms.string("10000000000")),
            cms.PSet(max=cms.string("10000000000")),
            cms.PSet(min=cms.string("-10000000000")),
            cms.PSet(min=cms.string("-10000000000"))
            ),
             )
    )

#cuts here mimic the miniAOD photon cuts and the non-category based trigger cuts
#Also included: the super-loose ID MVA cuts
flashggPreselectedDiPhotons = cms.EDFilter(
    "GenericDiPhotonCandidateSelector",
    src = cms.InputTag(""),
    rho = cms.InputTag("fixedGridRhoAll"),
    cut = cms.string(
        "1"
        ),
    variables = rediscoveryHLTvariables,
    categories = rediscoveryHLTcutsV1
    )
