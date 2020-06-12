import FWCore.ParameterSet.Config as cms

# No default. Latest is only in the EGM tool
photonSmearBins = cms.PSet()
photonScaleUncertBins = cms.PSet()

mvaShiftBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)"),
    bins = cms.VPSet(
                     cms.PSet( lowBounds = cms.vdouble(0.000), upBounds = cms.vdouble(999.), values = cms.vdouble( 0.0 ), uncertainties = cms.vdouble( 0.00 ))
                     )
    )

# Lowmass Preselection SF and uncertainties of 2018 : Prasant https://indico.cern.ch/event/890825/contributions/3762465/attachments/1992801/3323372/Update_Preselection_LoosePhotonIDMVA_Autumn18_SF_24022020.pdf 
preselBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.0 ) , upBounds = cms.vdouble( 1.5, 0.85 ) , values = cms.vdouble( 1.0068 ) , uncertainties = cms.vdouble( 0.0588 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85 ), upBounds = cms.vdouble( 1.5, 999.0 ) , values = cms.vdouble( 1.0053 ) , uncertainties = cms.vdouble( 0.0201 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.0 ) , upBounds = cms.vdouble( 6.0, 0.9) , values = cms.vdouble( 1.0 ) , uncertainties = cms.vdouble( 0.0 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.9 ) , upBounds = cms.vdouble( 6.0, 999.0 ) , values = cms.vdouble( 1.0288 ) , uncertainties = cms.vdouble( 0.0060 )  )
        )
    )


# JTao: Low mass 2018 case https://indico.cern.ch/event/873822/contributions/3712516/attachments/1972421/3281504/202001_LM2018PreliminaryResults_present.pdf
electronVetoBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.00 ) , upBounds = cms.vdouble( 1.5, 0.85 ) , values = cms.vdouble( 0.9683 ) , uncertainties = cms.vdouble( 0.0053 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85 ) , upBounds = cms.vdouble( 1.5, 999. ) , values = cms.vdouble( 0.9676 ) , uncertainties = cms.vdouble( 0.0015 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.00 ) , upBounds = cms.vdouble( 6.0, 0.90 ) , values = cms.vdouble( 1.0 ) , uncertainties = cms.vdouble( 0.0 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.90 ) , upBounds = cms.vdouble( 6.0, 999. ) , values = cms.vdouble( 0.8844 ) , uncertainties = cms.vdouble( 0.0073 )  )
        )
    )


# JMalcles - based on JTao SF + ttH efficiencies. calculated to preserve nTot ttH. 

leadPixelSeedBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9", "hasPixelSeed"),
    bins = cms.VPSet(
        # No Pixel Seed
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.00, -0.1 ) , upBounds = cms.vdouble( 1.5, 0.85 , 0.1) , values = cms.vdouble(0.978 ) , uncertainties = cms.vdouble( -0.00401807 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85, -0.1 ) , upBounds = cms.vdouble( 1.5, 999. , 0.1) , values = cms.vdouble(0.9824) , uncertainties = cms.vdouble( -0.00200421 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.00, -0.1 ) , upBounds = cms.vdouble( 6.0, 0.90 , 0.1) , values = cms.vdouble(0.9168 ) , uncertainties = cms.vdouble( -0.0224756 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.90 ,-0.1) , upBounds = cms.vdouble( 6.0, 999., 0.1) , values = cms.vdouble( 0.9403 ) , uncertainties = cms.vdouble(  -0.00631264  )  ),        
        # Yes Pixel Seed
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.00 ,0.9) , upBounds = cms.vdouble( 1.5, 0.85, 1.1 ) , values = cms.vdouble( 1.08876 ) , uncertainties = cms.vdouble( 0.0162106 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85 ,0.9) , upBounds = cms.vdouble( 1.5, 999., 1.1 ) , values = cms.vdouble( 1.5961) , uncertainties = cms.vdouble( 0.0678807 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.00 ,0.9) , upBounds = cms.vdouble( 6.0, 0.90, 1.1 ) , values = cms.vdouble( 1.09763 ) , uncertainties = cms.vdouble( 0.0263745 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.90 ,0.9) , upBounds = cms.vdouble( 6.0, 999., 1.1 ) , values = cms.vdouble( 1.20264 ) , uncertainties = cms.vdouble(  0.0214274 )  )
        )
    )


subleadPixelSeedBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9", "hasPixelSeed"),
    bins = cms.VPSet(
        # No Pixel Seed
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.00, -0.1 ) , upBounds = cms.vdouble( 1.5, 0.85 , 0.1) , values = cms.vdouble( 0.978  ) , uncertainties = cms.vdouble(  -0.00415083)  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85, -0.1 ) , upBounds = cms.vdouble( 1.5, 999. , 0.1) , values = cms.vdouble( 0.9824 ) , uncertainties = cms.vdouble( -0.00280026 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.00, -0.1 ) , upBounds = cms.vdouble( 6.0, 0.90 , 0.1) , values = cms.vdouble( 0.9168 ) , uncertainties = cms.vdouble( -0.0225538  )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.90 ,-0.1) , upBounds = cms.vdouble( 6.0, 999., 0.1) , values = cms.vdouble( 0.9403 ) , uncertainties = cms.vdouble( -0.00655045 )  ),        
        # Yes Pixel Seed
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.00 ,0.9) , upBounds = cms.vdouble( 1.5, 0.85, 1.1 ) , values = cms.vdouble( 1.13196) , uncertainties = cms.vdouble( 0.0248967 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85 ,0.9) , upBounds = cms.vdouble( 1.5, 999., 1.1 ) , values = cms.vdouble( 1.61512 ) , uncertainties = cms.vdouble(  0.0978689 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.00 ,0.9) , upBounds = cms.vdouble( 6.0, 0.90, 1.1 ) , values = cms.vdouble( 1.10623 ) , uncertainties = cms.vdouble(  0.0287957 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.90 ,0.9) , upBounds = cms.vdouble( 6.0, 999., 1.1 ) , values = cms.vdouble( 1.20311 ) , uncertainties = cms.vdouble(  0.0222861 )  )
        )
    )

FNUFBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.00 ) , upBounds = cms.vdouble( 1.5, 0.94 ) , values = cms.vdouble( 0. ) , uncertainties = cms.vdouble(  0.0007 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.94 ) , upBounds = cms.vdouble( 1.5, 999. ) , values = cms.vdouble( 0. ) , uncertainties = cms.vdouble( -0.0007 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.00 ) , upBounds = cms.vdouble( 6.0, 0.94 ) , values = cms.vdouble( 0. ) , uncertainties = cms.vdouble(  0.0007 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.94 ) , upBounds = cms.vdouble( 6.0, 999. ) , values = cms.vdouble( 0. ) , uncertainties = cms.vdouble( -0.0007 )  )
        )
    )

# from Martina: https://indico.cern.ch/event/628676/contributions/2546615/attachments/1440085/2216643/20170405_martina_regrEchecksUpdate.pdf
showerShapeBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        #EB low R9
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.00 ) , upBounds = cms.vdouble( 1.0, 0.94 ) , values = cms.vdouble( 0. ) , uncertainties = cms.vdouble( -0.0001 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.0, 0.00 ) , upBounds = cms.vdouble( 1.5, 0.94 ) , values = cms.vdouble( 0. ) , uncertainties = cms.vdouble(  0.0002 )  ) ,
        #EB high R9
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.94 ) , upBounds = cms.vdouble( 1.0, 999. ) , values = cms.vdouble( 0. ) , uncertainties = cms.vdouble( -0.0006 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.0, 0.94 ) , upBounds = cms.vdouble( 1.5, 999. ) , values = cms.vdouble( 0. ) , uncertainties = cms.vdouble( -0.0011 )  ) ,
        #EE low R9
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.00 ) , upBounds = cms.vdouble( 2.0, 0.94 ) , values = cms.vdouble( 0. ) , uncertainties = cms.vdouble( 0.0015 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 2.0, 0.00 ) , upBounds = cms.vdouble( 6.0, 0.94 ) , values = cms.vdouble( 0. ) , uncertainties = cms.vdouble( 0.0004 )  ) ,
        #EE high R9
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.94 ) , upBounds = cms.vdouble( 2.0, 999. ) , values = cms.vdouble( 0. ) , uncertainties = cms.vdouble( 0.0002 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 2.0, 0.94 ) , upBounds = cms.vdouble( 6.0, 999. ) , values = cms.vdouble( 0. ) , uncertainties = cms.vdouble( 0.0003 )  ) 
        )
    )



#with full 2017 dataset from Linda on April 30 for preapproval.  Trigger scale factors for use without HLT applied in MC
leadTriggerScaleBins = cms.PSet(
    variables = cms.vstring("full5x5_r9","abs(superCluster.eta)","pt"),
    bins = cms.VPSet(


	#pt binning seeded leg: 0.,  35., 37., 40., 45., 50., 60., 70., 90., 300	
        
	###BARREL	
	#0 <eta < 1.5, R9<0.50 ==> No photons with r9 < 0.5 in the barrel
        cms.PSet(lowBounds = cms.vdouble(0.0,0.,0.), upBounds = cms.vdouble(0.50,1.5,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),


        #0 <eta < 1.5, 0.50<R9<0.54 
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,0.), upBounds = cms.vdouble(0.54,1.5,35), values = cms.vdouble(0.7809388692), uncertainties = cms.vdouble(0.0056916545,0.0056916545)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,35.), upBounds = cms.vdouble(0.54,1.5,37.), values = cms.vdouble(0.9191352822), uncertainties = cms.vdouble(0.0035817647,0.0035817647)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,37.), upBounds = cms.vdouble(0.54,1.5,40.), values = cms.vdouble(0.9309284430), uncertainties = cms.vdouble(0.0043363573,0.0043363573)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,40.), upBounds = cms.vdouble(0.54,1.5,45.), values = cms.vdouble(0.9283272209), uncertainties = cms.vdouble(0.0023562783,0.0023562783)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,45.), upBounds = cms.vdouble(0.54,1.5,50.), values = cms.vdouble(0.9316272382), uncertainties = cms.vdouble(0.0040759286,0.0040759286)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,50.), upBounds = cms.vdouble(0.54,1.5,60.), values = cms.vdouble(0.9448619060), uncertainties = cms.vdouble(0.0151644118,0.0151644118)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,60.), upBounds = cms.vdouble(0.54,1.5,70.), values = cms.vdouble(0.9477220486), uncertainties = cms.vdouble(0.0258457214,0.0258457214)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,70.), upBounds = cms.vdouble(0.54,1.5,90.), values = cms.vdouble(0.9455777786), uncertainties = cms.vdouble(0.0131709562,0.0131709562)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,90.), upBounds = cms.vdouble(0.54,1.5,999999.), values = cms.vdouble(0.9282674041), uncertainties = cms.vdouble(0.0350558270,0.0350558270)), 
                 
                 
        #0 <eta < 1.5, 0.54<R9<0.85 
        cms.PSet(lowBounds = cms.vdouble(0.54,0.,0.), upBounds = cms.vdouble(0.85,1.5,35.), values = cms.vdouble(0.8209142931), uncertainties = cms.vdouble(0.0010115643,0.0010115643)),
        cms.PSet(lowBounds = cms.vdouble(0.54,0.,35.), upBounds = cms.vdouble(0.85,1.5,37.), values = cms.vdouble(0.9473472272), uncertainties = cms.vdouble(0.0015686238,0.0015686238)), 
        cms.PSet(lowBounds = cms.vdouble(0.54,0.,37.), upBounds = cms.vdouble(0.85,1.5,40.), values = cms.vdouble(0.9545920641), uncertainties = cms.vdouble(0.0031514496,0.0031514496)),
        cms.PSet(lowBounds = cms.vdouble(0.54,0.,40.), upBounds = cms.vdouble(0.85,1.5,45.), values = cms.vdouble(0.9563820089), uncertainties = cms.vdouble(0.0010016443,0.0010016443)),
        cms.PSet(lowBounds = cms.vdouble(0.54,0.,45.), upBounds = cms.vdouble(0.85,1.5,50.), values = cms.vdouble(0.9596097798), uncertainties = cms.vdouble(0.0011668878,0.0011668878)),
        cms.PSet(lowBounds = cms.vdouble(0.54,0.,50.), upBounds = cms.vdouble(0.85,1.5,60.), values = cms.vdouble(0.9613614182), uncertainties = cms.vdouble(0.0023494340,0.0023494340)),
        cms.PSet(lowBounds = cms.vdouble(0.54,0.,60.), upBounds = cms.vdouble(0.85,1.5,70.), values = cms.vdouble(0.9683269990), uncertainties = cms.vdouble(0.0016207766,0.0016207766)),
        cms.PSet(lowBounds = cms.vdouble(0.54,0.,70.), upBounds = cms.vdouble(0.85,1.5,90.), values = cms.vdouble(0.9755505038), uncertainties = cms.vdouble(0.0022129914,0.0022129914)),
        cms.PSet(lowBounds = cms.vdouble(0.54,0.,90.), upBounds = cms.vdouble(0.85,1.5,999999.), values = cms.vdouble(0.9728120369), uncertainties = cms.vdouble(0.0064619954,0.0064619954)),


        
        #0 <eta < 1.5, 0.85<R9<999 
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,0.), upBounds = cms.vdouble(999,1.5,35.), values = cms.vdouble(0.8777756549), uncertainties = cms.vdouble(0.0078958838,0.0078958838)),

        cms.PSet(lowBounds = cms.vdouble(0.85,0.,35.), upBounds = cms.vdouble(999,1.5,37.), values = cms.vdouble(0.9591566003), uncertainties = cms.vdouble(0.0011344872,0.0011344872)),

        cms.PSet(lowBounds = cms.vdouble(0.85,0.,37.), upBounds = cms.vdouble(999,1.5,40.), values = cms.vdouble(0.9661504394), uncertainties = cms.vdouble(0.0011211070,0.0011211070)),

        cms.PSet(lowBounds = cms.vdouble(0.85,0.,40.), upBounds = cms.vdouble(999,1.5,45.), values = cms.vdouble(0.9713429409), uncertainties = cms.vdouble(0.0010005257,0.0010005257)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,45.), upBounds = cms.vdouble(999,1.5,50.), values = cms.vdouble(0.9766069398), uncertainties = cms.vdouble(0.0010598956,0.0010598956)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,50.), upBounds = cms.vdouble(999,1.5,60.), values = cms.vdouble(0.9785793731), uncertainties = cms.vdouble(0.0019282260,0.0019282260)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,60.), upBounds = cms.vdouble(999,1.5,70.), values = cms.vdouble(0.9784872329), uncertainties = cms.vdouble(0.0033155470,0.0033155470)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,70.), upBounds = cms.vdouble(999,1.5,90.), values = cms.vdouble(0.9807043393), uncertainties = cms.vdouble(0.0011278073,0.0011278073)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,90.), upBounds = cms.vdouble(999,1.5,999999.), values = cms.vdouble(0.9872168099), uncertainties = cms.vdouble(0.0010175510,0.0010175510)),


	###ENDCAPS
 	#1.5 <eta < 3, R9<0.9 ==> No low-R9 photons for te seeded leg/leading photon in the Endcaps
        cms.PSet(lowBounds = cms.vdouble(0.,1.5,0.), upBounds = cms.vdouble(0.9,3.,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),

        #1.5 <eta < 3, 0.9<R9<0.93 ==> R9 turn-on bin - Endcaps
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,0), upBounds = cms.vdouble(0.93,3.,35.), values = cms.vdouble(0.8336632168), uncertainties = cms.vdouble(0.0063368601,0.0063368601)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,35.), upBounds = cms.vdouble(0.93,3.,37.), values = cms.vdouble(0.9593305114), uncertainties = cms.vdouble(0.0016819998,0.0016819998)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,37.), upBounds = cms.vdouble(0.93,3.,40.), values = cms.vdouble(0.9639324466), uncertainties = cms.vdouble(0.0017621054,0.0017621054)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,40.), upBounds = cms.vdouble(0.93,3.,45.), values = cms.vdouble(0.9673347931), uncertainties = cms.vdouble(0.0010168368,0.0010168368)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,45.), upBounds = cms.vdouble(0.93,3.,50.), values = cms.vdouble(0.9730845967), uncertainties = cms.vdouble(0.0012127740,0.0012127740)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,50.), upBounds = cms.vdouble(0.93,3.,60.), values = cms.vdouble(0.9700713622), uncertainties = cms.vdouble(0.0018934304,0.0018934304)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,60.), upBounds = cms.vdouble(0.93,3.,70.), values = cms.vdouble(0.9753541182), uncertainties = cms.vdouble(0.0024945324,0.0024945324)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,70.), upBounds = cms.vdouble(0.93,3.,90.), values = cms.vdouble(0.9778381112), uncertainties = cms.vdouble(0.0027074031,0.0027074031)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,90.), upBounds = cms.vdouble(0.93,3.,9999.), values = cms.vdouble(0.9840133922), uncertainties = cms.vdouble(0.0032652893,0.0032652893)),

              
        #1.5 <eta < 3, 0.93<R9<999 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,0), upBounds = cms.vdouble(999,3.,35.), values = cms.vdouble(0.8453711721), uncertainties = cms.vdouble(0.0020634406,0.0020634406)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,35.), upBounds = cms.vdouble(999,3.,37.), values = cms.vdouble(0.9712187003), uncertainties = cms.vdouble(0.0011020605,0.0011020605)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,37.), upBounds = cms.vdouble(999,3.,40.), values = cms.vdouble(0.9768360010), uncertainties = cms.vdouble(0.0011907161,0.0011907161)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,40.), upBounds = cms.vdouble(999,3.,45.), values = cms.vdouble(0.9803619699), uncertainties = cms.vdouble(0.0012222053,0.0012222053)), 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,45.), upBounds = cms.vdouble(999,3.,50.), values = cms.vdouble(0.9838732468), uncertainties = cms.vdouble(0.0043378749,0.0043378749)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,50.), upBounds = cms.vdouble(999,3.,60.), values = cms.vdouble(0.9849675199), uncertainties = cms.vdouble(0.0098590040,0.0098590040)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,60.), upBounds = cms.vdouble(999,3.,70.), values = cms.vdouble(0.9905641943), uncertainties = cms.vdouble(0.0029688155,0.0029688155)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,70.), upBounds = cms.vdouble(999,3.,90.), values = cms.vdouble(0.9891408232), uncertainties = cms.vdouble(0.0072279420,0.0072279420)), 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,90.), upBounds = cms.vdouble(999,3.,9999.), values = cms.vdouble(0.9931558912), uncertainties = cms.vdouble(0.0021266845,0.0021266845)), 


        #eta > 3. not used 
        cms.PSet(lowBounds = cms.vdouble(0.0,3.,0.), upBounds = cms.vdouble(999.,999.,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.))
        )
    )


subleadTriggerScaleBins = cms.PSet(
    variables = cms.vstring("full5x5_r9","abs(superCluster.eta)","pt"),
    bins = cms.VPSet(

	#pt binning for unseeded leg 0., 25., 28., 31., 35.,40.,45.,50.,60.,70.,90.,300.	
        

	###BARREL
	#0. <eta < 1.5, R9<0.5 ==> No photons with r9 < 0.5 in the barrel
        cms.PSet(lowBounds = cms.vdouble(0.0,0,0.), upBounds = cms.vdouble(0.5,1.5,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),

        #eta < 1.5, 0.5<R9<0.54 ==> R9 Turn-on bin - low R9 - Barrel 
        cms.PSet(lowBounds = cms.vdouble(0.50,0,0.), upBounds = cms.vdouble(0.54,1.5,25.), values = cms.vdouble(0.9138844940), uncertainties = cms.vdouble(0.0095172890,0.0095172890)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,25.), upBounds = cms.vdouble(0.54,1.5,28.), values = cms.vdouble(0.9391151967), uncertainties = cms.vdouble(0.0091460702,0.0091460702)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0,28.), upBounds = cms.vdouble(0.54,1.5,31.), values = cms.vdouble(0.9471744557), uncertainties = cms.vdouble(0.0094074668,0.0094074668)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0,31.), upBounds = cms.vdouble(0.54,1.5,35.), values = cms.vdouble(0.9657520875), uncertainties = cms.vdouble(0.0039026683,0.0039026683)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,35.), upBounds = cms.vdouble(0.54,1.5,40.), values = cms.vdouble(0.9726946474), uncertainties = cms.vdouble(0.0046277278,0.0046277278)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,40.), upBounds = cms.vdouble(0.54,1.5,45.), values = cms.vdouble(0.9735944404), uncertainties = cms.vdouble(0.0018454424,0.0018454424)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,45.), upBounds = cms.vdouble(0.54,1.5,50.), values = cms.vdouble(0.9736234735), uncertainties = cms.vdouble(0.0019540762,0.0019540762)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0,50.), upBounds = cms.vdouble(0.54,1.5,60.), values = cms.vdouble(0.9786018408), uncertainties = cms.vdouble(0.0143690034,0.0143690034)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,60.), upBounds = cms.vdouble(0.54,1.5,70.), values = cms.vdouble(0.9704059294), uncertainties = cms.vdouble(0.0184516481,0.0184516481)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,70.), upBounds = cms.vdouble(0.54,1.5,90.), values = cms.vdouble(0.9666644891), uncertainties = cms.vdouble(0.0170072811,0.0170072811)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,90.), upBounds = cms.vdouble(0.54,1.5,999999.), values = cms.vdouble(0.9473277036), uncertainties = cms.vdouble(0.0479759614,0.0479759614)),



        #eta < 1.5, 0.54<R9<0.85 ==> low R9 - Barrel
        cms.PSet(lowBounds = cms.vdouble(0.54,0,0.), upBounds = cms.vdouble(0.85,1.5,25.), values = cms.vdouble(0.9437983952), uncertainties = cms.vdouble(0.0033166891,0.0033166891)), 
        cms.PSet(lowBounds = cms.vdouble(0.54,0,25.), upBounds = cms.vdouble(0.85,1.5,28.), values = cms.vdouble(0.9758766263), uncertainties = cms.vdouble(0.0037443283,0.0037443283)),
        cms.PSet(lowBounds = cms.vdouble(0.54,0,28.), upBounds = cms.vdouble(0.85,1.5,31.), values = cms.vdouble(0.9828895074), uncertainties = cms.vdouble(0.0014681851,0.0014681851)), 
        cms.PSet(lowBounds = cms.vdouble(0.54,0,31.), upBounds = cms.vdouble(0.85,1.5,35.), values = cms.vdouble(0.9876097195), uncertainties = cms.vdouble(0.0018252738,0.0018252738)), 
        cms.PSet(lowBounds = cms.vdouble(0.54,0,35.), upBounds = cms.vdouble(0.85,1.5,40.), values = cms.vdouble(0.9892959590), uncertainties = cms.vdouble(0.0010888097,0.0010888097)),
        cms.PSet(lowBounds = cms.vdouble(0.54,0,40.), upBounds = cms.vdouble(0.85,1.5,45.), values = cms.vdouble(0.9897007026), uncertainties = cms.vdouble(0.0010931961,0.0010931961)),
        cms.PSet(lowBounds = cms.vdouble(0.54,0,45.), upBounds = cms.vdouble(0.85,1.5,50.), values = cms.vdouble(0.9882825409), uncertainties = cms.vdouble(0.0010451187,0.0010451187)), 
        cms.PSet(lowBounds = cms.vdouble(0.54,0,50.), upBounds = cms.vdouble(0.85,1.5,60.), values = cms.vdouble(0.9879121093), uncertainties = cms.vdouble(0.0010589355,0.0010589355)),
        cms.PSet(lowBounds = cms.vdouble(0.54,0,60.), upBounds = cms.vdouble(0.85,1.5,70.), values = cms.vdouble(0.9919510847), uncertainties = cms.vdouble(0.0027610542,0.0027610542)), 
        cms.PSet(lowBounds = cms.vdouble(0.54,0,70.), upBounds = cms.vdouble(0.85,1.5,90.), values = cms.vdouble(0.9913742928), uncertainties = cms.vdouble(0.0015724445,0.0015724445)), 
        cms.PSet(lowBounds = cms.vdouble(0.54,0,90.), upBounds = cms.vdouble(0.85,1.5,999999.), values = cms.vdouble(0.9852282167), uncertainties = cms.vdouble(0.0046772429,0.0046772429)),

	#eta < 1.5, 0.85<R9<999 ==> High R9 - Barrel	       
          cms.PSet(lowBounds = cms.vdouble(0.85,0,0.), upBounds = cms.vdouble(999,1.5,25.), values = cms.vdouble(0.9672297042), uncertainties = cms.vdouble(0.0015837834,0.0015837834)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0,25.), upBounds = cms.vdouble(999,1.5,28.), values = cms.vdouble(0.9805682395), uncertainties = cms.vdouble(0.0032401849,0.0032401849)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,28.), upBounds = cms.vdouble(999,1.5,31.), values = cms.vdouble(0.9857864878), uncertainties = cms.vdouble(0.0035391237,0.0035391237)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,31.), upBounds = cms.vdouble(999,1.5,35.), values = cms.vdouble(0.9899655333), uncertainties = cms.vdouble(0.0014305149,0.0014305149)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,35.), upBounds = cms.vdouble(999,1.5,40.), values = cms.vdouble(0.9926847217), uncertainties = cms.vdouble(0.0010036853,0.0010036853)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,40.), upBounds = cms.vdouble(999,1.5,45.), values = cms.vdouble(0.9936363464), uncertainties = cms.vdouble(0.0010018311,0.0010018311)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,45.), upBounds = cms.vdouble(999,1.5,50.), values = cms.vdouble(0.9938550364), uncertainties = cms.vdouble(0.0022661436,0.0022661436)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,50.), upBounds = cms.vdouble(999,1.5,60.), values = cms.vdouble(0.9942535185), uncertainties = cms.vdouble(0.0010869726,0.0010869726)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,60.), upBounds = cms.vdouble(999,1.5,70.), values = cms.vdouble(0.9949263283), uncertainties = cms.vdouble(0.0011394169,0.0011394169)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,70.), upBounds = cms.vdouble(999,1.5,90.), values = cms.vdouble(0.9948820767), uncertainties = cms.vdouble(0.0017256019,0.0017256019)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,90.), upBounds = cms.vdouble(999,1.5,999999.), values = cms.vdouble(0.9949257401), uncertainties = cms.vdouble(0.0012735262,0.0012735262)), 

	###ENDCAPS	
	#1.5 <eta < 3., R9<0.9 ==> No photons with r9 < 0.9 in the endcaps
        cms.PSet(lowBounds = cms.vdouble(0.0,1.5,0.), upBounds = cms.vdouble(0.9,3.,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),

        #1.5 <eta < 3, 0.9<R9<0.93 ==> R9 Turn-on bin - low R9 - Endcaps
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,0.), upBounds = cms.vdouble(0.93,3.,25.), values = cms.vdouble(0.9121171624), uncertainties = cms.vdouble(0.0022336702,0.0022336702)),
 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,25.), upBounds = cms.vdouble(0.93,3.,28.), values = cms.vdouble(0.9543055879), uncertainties = cms.vdouble(0.0065733116,0.0065733116)),

        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,28.), upBounds = cms.vdouble(0.93,3.,31.), values = cms.vdouble(0.9584879018), uncertainties = cms.vdouble(0.0038973106,0.0038973106)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,31.), upBounds = cms.vdouble(0.93,3.,35.), values = cms.vdouble(0.9665292836), uncertainties = cms.vdouble(0.0013507542,0.0013507542)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,35.), upBounds = cms.vdouble(0.93,3.,40.), values = cms.vdouble(0.9724877347), uncertainties = cms.vdouble(0.0011959782,0.0011959782)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,40.), upBounds = cms.vdouble(0.93,3.,45.), values = cms.vdouble(0.9756256530), uncertainties = cms.vdouble(0.0010215966,0.0010215966)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,45.), upBounds = cms.vdouble(0.93,3.,50.), values = cms.vdouble(0.9794385124), uncertainties = cms.vdouble(0.0010991265,0.0010991265)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,50.), upBounds = cms.vdouble(0.93,3.,60.), values = cms.vdouble(0.9785902688), uncertainties = cms.vdouble(0.0013014346,0.0013014346)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,60.), upBounds = cms.vdouble(0.93,3.,70.), values = cms.vdouble(0.9771428485), uncertainties = cms.vdouble(0.0023350841,0.0023350841)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,70.), upBounds = cms.vdouble(0.93,3.,90.), values = cms.vdouble(0.9780830492), uncertainties = cms.vdouble(0.0026833573,0.0026833573)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,90.), upBounds = cms.vdouble(0.93,3.,999999.), values = cms.vdouble(0.9811520203), uncertainties = cms.vdouble(0.0026233174,0.0026233174)),

      
	##1.5 <eta < 3, 0.93<R9<999 ==> high R9 - Endcaps	
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,0.), upBounds = cms.vdouble(999,3.,25.), values = cms.vdouble(0.9433319000), uncertainties = cms.vdouble(0.0036178299,0.0036178299)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,25.), upBounds = cms.vdouble(999,3.,28.), values = cms.vdouble(0.9836486488), uncertainties = cms.vdouble(0.0022606561,0.0022606561)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,28.), upBounds = cms.vdouble(999,3.,31.), values = cms.vdouble(0.9849810016), uncertainties = cms.vdouble(0.0037905292,0.0037905292)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,31.), upBounds = cms.vdouble(999,3.,35.), values = cms.vdouble(0.9898516787),  uncertainties = cms.vdouble(0.0012461983,0.0012461983)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,35.), upBounds = cms.vdouble(999,3.,40.), values = cms.vdouble(0.9925685115), uncertainties = cms.vdouble(0.0010564254,0.0010564254)), 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,40.), upBounds = cms.vdouble(999,3.,45.), values = cms.vdouble(0.9934337616), uncertainties = cms.vdouble(0.0010008334,0.0010008334)), 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,45.), upBounds = cms.vdouble(999,3.,50.), values = cms.vdouble(0.9932323204), uncertainties = cms.vdouble(0.0010568528,0.0010568528)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,50.), upBounds = cms.vdouble(999,3.,60.), values = cms.vdouble(0.9933956724), uncertainties = cms.vdouble(0.0010399275,0.0010399275)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,60.), upBounds = cms.vdouble(999,3.,70.), values = cms.vdouble(0.9955477196),  uncertainties = cms.vdouble(0.0010488087,0.0010488087)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,70.), upBounds = cms.vdouble(999,3.,90.), values = cms.vdouble(0.9935966416), uncertainties = cms.vdouble(0.0010357135,0.0010357135)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,90.), upBounds = cms.vdouble(999,3.,999999.), values = cms.vdouble(0.9966331942), uncertainties = cms.vdouble(0.0022287400,0.0022287400)),

         #eta > 3. not used
        cms.PSet(lowBounds = cms.vdouble(0.0,3.,0), upBounds = cms.vdouble(999,999,999999), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.))
        )
    )

# Lowmass  LooseIDMVA cut SF and uncertainities Update 2018: Prasant https://indico.cern.ch/event/890825/contributions/3762465/attachments/1992801/3323372/Update_Preselection_LoosePhotonIDMVA_Autumn18_SF_24022020.pdf
looseMvaBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.0  ) , upBounds = cms.vdouble( 1.5, 0.85  ) , values = cms.vdouble( 1.0012 ) , uncertainties = cms.vdouble( 0.0029 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85 ) , upBounds = cms.vdouble( 1.5, 999.0 ) , values = cms.vdouble( 0.9998 ) , uncertainties = cms.vdouble( 0.0007 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.0  ) , upBounds = cms.vdouble( 6.0, 0.9   ) , values = cms.vdouble( 1.0 ) , uncertainties = cms.vdouble( 0.0 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.9  ) , upBounds = cms.vdouble( 6.0, 999.0 ) , values = cms.vdouble( 1.0009 ) , uncertainties = cms.vdouble( 0.0006 ) ) ) )

# RELATIVE shift of sigmaE/E --> 0.05 corresponds to a shift of 5%
sigmaEOverEShiftBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)"),
    bins = cms.VPSet(
                     cms.PSet( lowBounds = cms.vdouble(0.000), upBounds = cms.vdouble(999.), values = cms.vdouble( 0.0 ), uncertainties = cms.vdouble( 0.05 ))
                     )
    )


RVBins = cms.PSet(
    variables = cms.vstring("pt"),
    bins = cms.VPSet(
        cms.PSet(lowBounds = cms.vdouble(0), upBounds = cms.vdouble(5), values = cms.vdouble(1.00046,0.999332), uncertainties = cms.vdouble(0.00167952,0.00167952,0.00243437,0.00243437)),
        cms.PSet(lowBounds = cms.vdouble(5), upBounds = cms.vdouble(10), values = cms.vdouble(1.01403,0.978954), uncertainties = cms.vdouble(0.00132967,0.00132967,0.00199484,0.00199484)),
        cms.PSet(lowBounds = cms.vdouble(10), upBounds = cms.vdouble(15), values = cms.vdouble(1.0031,0.99458), uncertainties = cms.vdouble(0.00127284,0.00127284,0.00222552,0.00222552)),
        cms.PSet(lowBounds = cms.vdouble(15), upBounds = cms.vdouble(20), values = cms.vdouble(0.992237,1.01713), uncertainties = cms.vdouble(0.00122729,0.00122729,0.00270787,0.00270787)),
        cms.PSet(lowBounds = cms.vdouble(20), upBounds = cms.vdouble(30), values = cms.vdouble(0.990433,1.02985), uncertainties = cms.vdouble(0.000854325,0.000854325,0.00266531,0.00266531)),
        cms.PSet(lowBounds = cms.vdouble(30), upBounds = cms.vdouble(40), values = cms.vdouble(0.988515,1.05637), uncertainties = cms.vdouble(0.000847473,0.000847473,0.00415923,0.00415923)),
        cms.PSet(lowBounds = cms.vdouble(40), upBounds = cms.vdouble(50), values = cms.vdouble(0.988526,1.07976), uncertainties = cms.vdouble(0.000864982,0.000864982,0.00601261,0.00601261)),
        cms.PSet(lowBounds = cms.vdouble(50), upBounds = cms.vdouble(60), values = cms.vdouble(0.988509,1.11643), uncertainties = cms.vdouble(0.000909363,0.000909363,0.00921419,0.00921419)),
        cms.PSet(lowBounds = cms.vdouble(60), upBounds = cms.vdouble(80), values = cms.vdouble(0.989606,1.14786), uncertainties = cms.vdouble(0.000690743,0.000690743,0.00982573,0.00982573)),
        cms.PSet(lowBounds = cms.vdouble(80), upBounds = cms.vdouble(100), values = cms.vdouble(0.991492,1.16885), uncertainties = cms.vdouble(0.000759541,0.000759541,0.0150743,0.0150743)),
        cms.PSet(lowBounds = cms.vdouble(100), upBounds = cms.vdouble(140), values = cms.vdouble(0.997022,1.07771), uncertainties = cms.vdouble(0.00066297,0.00066297,0.0173001,0.0173001)),
        cms.PSet(lowBounds = cms.vdouble(140), upBounds = cms.vdouble(200), values = cms.vdouble(0.999255,1.02942), uncertainties = cms.vdouble(0.000738493,0.000738493,0.0291629,0.0291629)),
        cms.PSet(lowBounds = cms.vdouble(200), upBounds = cms.vdouble(400), values = cms.vdouble(1.00079,0.943138), uncertainties = cms.vdouble(0.000985164,0.000985164,0.0710487,0.0710487)),
        # maximum energy beyond 7000 because of wonky DiPhotons in data
        cms.PSet(lowBounds = cms.vdouble(400 ) , upBounds = cms.vdouble(999999999 ) , values = cms.vdouble(1,1              ) , uncertainties = cms.vdouble(0.,0.,0.,0.))       
        )
    )      
 
RVBinsNvtx = cms.PSet(
    variables = cms.vstring("nVert"),
     bins = cms.VPSet(
        cms.PSet(lowBounds = cms.vdouble(-0.5), upBounds = cms.vdouble(10.5), values = cms.vdouble(1.02898,0.828452), uncertainties = cms.vdouble(0.00155068,0.00155068,0.00918045,0.00918045)),
        cms.PSet(lowBounds = cms.vdouble(10.5), upBounds = cms.vdouble(12.5), values = cms.vdouble(1.00775,0.960156), uncertainties = cms.vdouble(0.0013211,0.0013211,0.00679271,0.00679271)),
        cms.PSet(lowBounds = cms.vdouble(12.5), upBounds = cms.vdouble(14.5), values = cms.vdouble(1.00406,0.980929), uncertainties = cms.vdouble(0.00113947,0.00113947,0.00535269,0.00535269)),
        cms.PSet(lowBounds = cms.vdouble(14.5), upBounds = cms.vdouble(16.5), values = cms.vdouble(1.00159,0.992869), uncertainties = cms.vdouble(0.00109956,0.00109956,0.00493888,0.00493888)),
        cms.PSet(lowBounds = cms.vdouble(16.5), upBounds = cms.vdouble(18.5), values = cms.vdouble(0.993201,1.02899), uncertainties = cms.vdouble(0.00112887,0.00112887,0.00481407,0.00481407)),
        cms.PSet(lowBounds = cms.vdouble(18.5), upBounds = cms.vdouble(20.5), values = cms.vdouble(0.991425,1.03468), uncertainties = cms.vdouble(0.0012414,0.0012414,0.00502105,0.00502105)),
        cms.PSet(lowBounds = cms.vdouble(20.5), upBounds = cms.vdouble(22.5), values = cms.vdouble(0.989716,1.03941), uncertainties = cms.vdouble(0.00142369,0.00142369,0.00545553,0.00545553)),
        cms.PSet(lowBounds = cms.vdouble(22.5), upBounds = cms.vdouble(25.5), values = cms.vdouble(0.98674,1.04837), uncertainties = cms.vdouble(0.00147513,0.00147513,0.00538112,0.00538112)),
        cms.PSet(lowBounds = cms.vdouble(25.5), upBounds = cms.vdouble(30.5), values = cms.vdouble(0.976922,1.07893), uncertainties = cms.vdouble(0.00188024,0.00188024,0.00643049,0.00643049)),
        cms.PSet(lowBounds = cms.vdouble(30.5), upBounds = cms.vdouble(100.5), values = cms.vdouble(0.959731,1.13018), uncertainties = cms.vdouble(0.00440431,0.00440431,0.0142389,0.0142389)),
        # just in case
        cms.PSet(lowBounds = cms.vdouble(100.5 ) , upBounds = cms.vdouble(999999999 ) , values = cms.vdouble(1,1              ) , uncertainties = cms.vdouble(0.,0.,0.,0.))       
        )
    )     




# Photon categoryscale [$\times 10^{-2}$]
# central EB (eta<0.8) low r9  0.035
# central EB (eta<0.8) high r9  0.033
# intermediate EB (0.8 < eta < 1.0) low r9  0.058
# intermediate EB (0.8 < eta < 1.0) high r9  0.12
# outer EB (eta>1.0) low r9  0.22
# outer EB (eta>1.0) high r9  0.34
# EE low r9  0.22
# EE high r9  0.34
#
# Copied from Run 1 values
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/HggLegacySystematicUncertainties#Upstream_material_electron_to_ph
# NB these are comparable to other scale uncertainties
materialBinsRun1 = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.0  ), upBounds = cms.vdouble( 0.8, 0.94  ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.00035 ) ),
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.94 ), upBounds = cms.vdouble( 0.8, 999.0 ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.00033 ) ),
        cms.PSet( lowBounds = cms.vdouble( 0.8, 0.0  ), upBounds = cms.vdouble( 1.0, 0.94  ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.00058 ) ),
        cms.PSet( lowBounds = cms.vdouble( 0.8, 0.94 ), upBounds = cms.vdouble( 1.0, 999.0 ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.0012 ) ),
        cms.PSet( lowBounds = cms.vdouble( 1.0, 0.0  ), upBounds = cms.vdouble( 1.5, 0.94  ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.0022 ) ),
        cms.PSet( lowBounds = cms.vdouble( 1.0, 0.94 ), upBounds = cms.vdouble( 1.5, 999.0 ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.0034 ) ),
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.0  ), upBounds = cms.vdouble( 999., 0.94  ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.0022 ) ),
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.94 ), upBounds = cms.vdouble( 999., 999.0 ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.0034 ) ),

        )
    )

materialBinsICHEP = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.0 ), upBounds = cms.vdouble( 1.5, 0.94 ), values = cms.vdouble( 0. ) , uncertainties = cms.vdouble( 0.00070 ) ),
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.94 ), upBounds = cms.vdouble( 1.5, 999. ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.00036 ) ),
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.0 ), upBounds = cms.vdouble( 6.0, 0.94  ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.00089 ) ),
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.9 ), upBounds = cms.vdouble( 6.0, 999. ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.00170 ) )
        )
    )

# Moriond17 studies from https://indico.cern.ch/event/628619/contributions/2552389/attachments/1442921/2222173/MaterialStudy_10042017_v2.pdf
#
# All numbers for 5% material uncertainty (x10^-3)
#
# Category1:CentralBarrelandlowr9 (|eta|<1.0andr9<0.94)     0.455
# Category2:CentralBarrelandhighr9(|eta|<1.0andr9>0.94)     0.233
# Category3:OuterBarrelandlowr9 (1.0<|eta|<1.5andr9<0.94)   2.089
# Category4:OuterBarrelandhighr9 (1.0<|eta|<1.5andr9>0.94)   N/A   (taken to be equal to cat3)
# Category5:Endcapandlowr9 (|eta|>1.5andr9<0.94)            1.090
# Category6:Endcapandhighr9 (|eta|>1.5andr9>0.94)           2.377

materialBinsMoriond17 = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.0  ), upBounds = cms.vdouble( 1.0, 0.94  ), values = cms.vdouble( 0. ),  uncertainties = cms.vdouble( 0.000455 ) ),
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.94 ), upBounds = cms.vdouble( 1.0, 999.0 ), values = cms.vdouble( 0. ),  uncertainties = cms.vdouble( 0.000233 ) ),
        cms.PSet( lowBounds = cms.vdouble( 1.0, 0.0  ), upBounds = cms.vdouble( 1.5, 0.94  ), values = cms.vdouble( 0. ),  uncertainties = cms.vdouble( 0.002089 ) ),
        cms.PSet( lowBounds = cms.vdouble( 1.0, 0.94 ), upBounds = cms.vdouble( 1.5, 999.0 ), values = cms.vdouble( 0. ),  uncertainties = cms.vdouble( 0.002089 ) ),
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.0  ), upBounds = cms.vdouble( 999., 0.94  ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.001090 ) ),
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.94 ), upBounds = cms.vdouble( 999., 999.0 ), values = cms.vdouble( 0. ), uncertainties = cms.vdouble( 0.002377 ) ),

        )
    )

emptyBins = cms.PSet(
    variables = cms.vstring("1"),
    bins = cms.VPSet()
    )

emptySigma = cms.PSet(
    firstVar = cms.vint32(),
    secondVar = cms.vint32()
)


scalesAndSmearingsPrefix = cms.string("EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2018_Step2Closure_CoarseEtaR9Gain_v2")
#scalesAndSmearingsPrefixForSigmaEOverE = scalesAndSmearingsPrefix

MCScaleHighR9EB = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("MCScaleHighR9EB"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)<1.5"),
          BinList = photonScaleUncertBins,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )

MCScaleLowR9EB = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("MCScaleLowR9EB"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("full5x5_r9<0.94&&abs(superCluster.eta)<1.5"),
          BinList = photonScaleUncertBins,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )

MCScaleHighR9EE = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("MCScaleHighR9EE"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)>=1.5"),
          BinList = photonScaleUncertBins,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )

MCScaleLowR9EE = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("MCScaleLowR9EE"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("full5x5_r9<0.94&&abs(superCluster.eta)>=1.5"),
          BinList = photonScaleUncertBins,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )

MaterialCentralBarrel = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("MaterialCentralBarrel"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("abs(superCluster.eta)<1.0"),
          BinList = materialBinsMoriond17,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )

MaterialOuterBarrel = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("MaterialOuterBarrel"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("abs(superCluster.eta)>=1.0&&abs(superCluster.eta)<1.5"),
          BinList = materialBinsMoriond17,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )

MaterialForward = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("MaterialForward"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("abs(superCluster.eta)>=1.5"),
          BinList = materialBinsMoriond17,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )

FNUFEB = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("FNUFEB"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("abs(superCluster.eta)<1.5"),
          BinList = FNUFBins,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )

FNUFEE = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("FNUFEE"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("abs(superCluster.eta)>=1.5"),
          BinList = FNUFBins,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )

ShowerShapeHighR9EB = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("ShowerShapeHighR9EB"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)<1.5"),
          BinList = showerShapeBins,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )

ShowerShapeHighR9EE = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("ShowerShapeHighR9EE"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)>=1.5"),
          BinList = showerShapeBins,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )

ShowerShapeLowR9EB = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("ShowerShapeLowR9EB"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("full5x5_r9<=0.94&&abs(superCluster.eta)<1.5"),
          BinList = showerShapeBins,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )

ShowerShapeLowR9EE = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScale"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("ShowerShapeLowR9EE"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("full5x5_r9<=0.94&&abs(superCluster.eta)>=1.5"),
          BinList = showerShapeBins,
          ApplyCentralValue = cms.bool(False),
          Debug = cms.untracked.bool(False)
          )


MCSmearHighR9EE = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSmearConstant"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("MCSmearHighR9EE"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)>=1.5"),
          BinList = photonSmearBins,
          # has to match the labels embedded in the photon object as
          # defined e.g. in flashgg/MicroAOD/python/flashggRandomizedPerPhotonDiPhotonProducer_cff.py
          #           or in flashgg/MicroAOD/python/flashggRandomizedPhotonProducer_cff.py (if at MicroAOD prod.)
          RandomLabel = cms.string("rnd_g_E"),
          Debug = cms.untracked.bool(False),
          ExaggerateShiftUp = cms.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

MCSmearLowR9EE = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSmearConstant"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("MCSmearLowR9EE"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("full5x5_r9<0.94&&abs(superCluster.eta)>=1.5"),
          BinList = photonSmearBins,
          # has to match the labels embedded in the photon object as
          # defined e.g. in flashgg/MicroAOD/python/flashggRandomizedPerPhotonDiPhotonProducer_cff.py
          #           or in flashgg/MicroAOD/python/flashggRandomizedPhotonProducer_cff.py (if at MicroAOD prod.)
          RandomLabel = cms.string("rnd_g_E"),
          Debug = cms.untracked.bool(False),
          ExaggerateShiftUp = cms.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

MCSmearHighR9EB = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSmearConstant"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("MCSmearHighR9EB"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)<1.5"),
          BinList = photonSmearBins,
          # has to match the labels embedded in the photon object as
          # defined e.g. in flashgg/MicroAOD/python/flashggRandomizedPerPhotonDiPhotonProducer_cff.py
          #           or in flashgg/MicroAOD/python/flashggRandomizedPhotonProducer_cff.py (if at MicroAOD prod.)
          RandomLabel = cms.string("rnd_g_E"),
          Debug = cms.untracked.bool(False),
          ExaggerateShiftUp = cms.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

MCSmearLowR9EB = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSmearConstant"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("MCSmearLowR9EB"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("full5x5_r9<=0.94&&abs(superCluster.eta)<1.5"),
          BinList = photonSmearBins,
          # has to match the labels embedded in the photon object as
          # defined e.g. in flashgg/MicroAOD/python/flashggRandomizedPerPhotonDiPhotonProducer_cff.py
          #           or in flashgg/MicroAOD/python/flashggRandomizedPhotonProducer_cff.py (if at MicroAOD prod.)
          RandomLabel = cms.string("rnd_g_E"),
          Debug = cms.untracked.bool(False),
          ExaggerateShiftUp = cms.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

MvaShift = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonMvaTransform"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("MvaShift"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("1"),
#          CorrectionFile = cms.FileInPath("flashgg/Systematics/data/SystematicsIDMVA_LegRunII_v1_2018.root"),
          CorrectionFile = cms.FileInPath("flashgg/Systematics/data/SystematicsIDMVA_LegRunII_LMDoubleEE_2018.root"),
          BinList = mvaShiftBins,
          Debug = cms.untracked.bool(False),
          ApplyCentralValue = cms.bool(False)
          )

PreselSF = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonWeight"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("PreselSF"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("1"),
          BinList = preselBins,
          Debug = cms.untracked.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

electronVetoSF = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonWeight"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("electronVetoSF"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("1"),
          BinList = electronVetoBins,
          Debug = cms.untracked.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

TriggerWeight = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonWeight"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("TriggerWeight"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("pt<99999"),
          BinList = leadTriggerScaleBins,
          BinList2 = subleadTriggerScaleBins,
          Debug = cms.untracked.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

PixelSeedWeight = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonWeight"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("PixelSeedWeight"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("pt<99999"),
          BinList = leadPixelSeedBins,
          BinList2 = subleadPixelSeedBins,
          Debug = cms.untracked.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

LooseMvaSF = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonWeight"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("LooseMvaSF"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("1"),
          BinList = looseMvaBins,
          Debug = cms.untracked.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

SigmaEOverEShift = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSigEOverEShift"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("SigmaEOverEShift"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("1"),
          BinList = sigmaEOverEShiftBins,
          Debug = cms.untracked.bool(False),
          ApplyCentralValue = cms.bool(False)
          )

SigmaEOverESmearing = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSigEoverESmearing"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("SigmaEOverESmearing"),
          NSigmas = cms.vint32(),
          OverallRange = cms.string("1"),
          BinList = photonSmearBins,
          Debug = cms.untracked.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

SigmaEOverESmearing_EGM = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSigEoverESmearingEGMTool"),
          MethodName = cms.string("FlashggDiPhotonFromPhoton"),
          Label = cms.string("SigmaEOverESmearing"),
          CorrectionFile = scalesAndSmearingsPrefix,
          NSigmas = cms.vint32(),
          OverallRange = cms.string("1"),
          BinList = emptyBins,
          Debug = cms.untracked.bool(False),
          ExaggerateShiftUp = cms.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

FracRVWeight = cms.PSet( MethodName = cms.string("FlashggDiPhotonWeightFromFracRV"),
          Label = cms.string("FracRVWeight"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("pt<99999"),
          BinList = RVBins,
          Debug = cms.untracked.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

FracRVNvtxWeight = cms.PSet( MethodName = cms.string("FlashggDiPhotonWeightFromFracRV"),
          Label = cms.string("FracRVNvtxWeight"),
          NSigmas = cms.vint32(-1,1),
          OverallRange = cms.string("nVert<99999"),
          BinList = RVBinsNvtx,
          Debug = cms.untracked.bool(False),
          ApplyCentralValue = cms.bool(True)
          )

MCSmearHighR9EE_EGM = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSmearStochasticEGMTool"),
         MethodName = cms.string("FlashggDiPhotonFromPhoton2D"),
         Label = cms.string("MCSmearHighR9EE"),
         FirstParameterName = cms.string("Rho"),
         SecondParameterName = cms.string("Phi"),
         CorrectionFile = scalesAndSmearingsPrefix,
         NSigmas = cms.PSet( firstVar = cms.vint32(1,-1,0,0),
                            secondVar = cms.vint32(0,0,1,-1)),
         OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)>=1.5"),
         BinList = emptyBins,
         # has to match the labels embedded in the photon object as
         # defined e.g. in flashgg/MicroAOD/python/flashggRandomizedPerPhotonDiPhotonProducer_cff.py
         #           or in flashgg/MicroAOD/python/flashggRandomizedPhotonProducer_cff.py (if at MicroAOD prod.)
         RandomLabel = cms.string("rnd_g_E"),
         Debug = cms.untracked.bool(False),
         ExaggerateShiftUp = cms.bool(False),
         ApplyCentralValue = cms.bool(True)
         )

MCSmearLowR9EE_EGM = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSmearStochasticEGMTool"),
         MethodName = cms.string("FlashggDiPhotonFromPhoton2D"),
         Label = cms.string("MCSmearLowR9EE"),
         FirstParameterName = cms.string("Rho"),
         SecondParameterName = cms.string("Phi"),
         CorrectionFile = scalesAndSmearingsPrefix,
         NSigmas = cms.PSet( firstVar = cms.vint32(1,-1,0,0),
                            secondVar = cms.vint32(0,0,1,-1)),
         OverallRange = cms.string("full5x5_r9<=0.94&&abs(superCluster.eta)>=1.5"),
         BinList = emptyBins,
         # has to match the labels embedded in the photon object as
         # defined e.g. in flashgg/MicroAOD/python/flashggRandomizedPerPhotonDiPhotonProducer_cff.py
         #           or in flashgg/MicroAOD/python/flashggRandomizedPhotonProducer_cff.py (if at MicroAOD prod.)
         RandomLabel = cms.string("rnd_g_E"),
         Debug = cms.untracked.bool(False),
         ExaggerateShiftUp = cms.bool(False),
         ApplyCentralValue = cms.bool(True)
         )

MCSmearHighR9EB_EGM = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSmearStochasticEGMTool"),
         MethodName = cms.string("FlashggDiPhotonFromPhoton2D"),
         Label = cms.string("MCSmearHighR9EB"),
         FirstParameterName = cms.string("Rho"),
         SecondParameterName = cms.string("Phi"),
         CorrectionFile = scalesAndSmearingsPrefix,
         NSigmas = cms.PSet( firstVar = cms.vint32(1,-1,0,0),
                             secondVar = cms.vint32(0,0,1,-1)),
         OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)<1.5"),
         BinList = emptyBins,
         # has to match the labels embedded in the photon object as
         # defined e.g. in flashgg/MicroAOD/python/flashggRandomizedPerPhotonDiPhotonProducer_cff.py
         #           or in flashgg/MicroAOD/python/flashggRandomizedPhotonProducer_cff.py (if at MicroAOD prod.)
         RandomLabel = cms.string("rnd_g_E"),
         Debug = cms.untracked.bool(False),
         ExaggerateShiftUp = cms.bool(False),
         ApplyCentralValue = cms.bool(True)
         )

MCSmearLowR9EB_EGM = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonSmearStochasticEGMTool"),
         MethodName = cms.string("FlashggDiPhotonFromPhoton2D"),
         Label = cms.string("MCSmearLowR9EB"),
         FirstParameterName = cms.string("Rho"),
         SecondParameterName = cms.string("Phi"),
         CorrectionFile = scalesAndSmearingsPrefix,
         NSigmas = cms.PSet( firstVar = cms.vint32(1,-1,0,0),
                            secondVar = cms.vint32(0,0,1,-1)),
         OverallRange = cms.string("full5x5_r9<=0.94&&abs(superCluster.eta)<1.5"),
         BinList = emptyBins,
         # has to match the labels embedded in the photon object as
         # defined e.g. in flashgg/MicroAOD/python/flashggRandomizedPerPhotonDiPhotonProducer_cff.py
         #           or in flashgg/MicroAOD/python/flashggRandomizedPhotonProducer_cff.py (if at MicroAOD prod.)
         RandomLabel = cms.string("rnd_g_E"),
         Debug = cms.untracked.bool(False),
         ExaggerateShiftUp = cms.bool(False),
         ApplyCentralValue = cms.bool(True)
         )

MCScaleHighR9EB_EGM = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScaleEGMTool"),
         MethodName = cms.string("FlashggDiPhotonFromPhoton"),
         Label = cms.string("MCScaleHighR9EB"),
         NSigmas = cms.vint32(-1,1),
         OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)<1.5"),
         BinList = emptyBins,
         CorrectionFile = scalesAndSmearingsPrefix,
         ApplyCentralValue = cms.bool(False),
         UncertaintyBitMask = cms.string("011"),#cms.string("110"),
         ExaggerateShiftUp = cms.bool(False),
         Debug = cms.untracked.bool(False)
         )

MCScaleLowR9EB_EGM = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScaleEGMTool"),
         MethodName = cms.string("FlashggDiPhotonFromPhoton"),
         Label = cms.string("MCScaleLowR9EB"),
         NSigmas = cms.vint32(-1,1),
         OverallRange = cms.string("full5x5_r9<0.94&&abs(superCluster.eta)<1.5"),
         BinList = emptyBins,
         CorrectionFile = scalesAndSmearingsPrefix,
         ApplyCentralValue = cms.bool(False),
         UncertaintyBitMask = cms.string("011"),#cms.string("110"),
         ExaggerateShiftUp = cms.bool(False),
         Debug = cms.untracked.bool(False)
         )

MCScaleHighR9EE_EGM = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScaleEGMTool"),
         MethodName = cms.string("FlashggDiPhotonFromPhoton"),
         Label = cms.string("MCScaleHighR9EE"),
         NSigmas = cms.vint32(-1,1),
         OverallRange = cms.string("full5x5_r9>0.94&&abs(superCluster.eta)>=1.5"),
         BinList = emptyBins,
         CorrectionFile = scalesAndSmearingsPrefix,
         ApplyCentralValue = cms.bool(False),
         UncertaintyBitMask = cms.string("011"),#cms.string("110"),
         ExaggerateShiftUp = cms.bool(False),
         Debug = cms.untracked.bool(False)
         )

MCScaleLowR9EE_EGM = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScaleEGMTool"),
         MethodName = cms.string("FlashggDiPhotonFromPhoton"),
         Label = cms.string("MCScaleLowR9EE"),
         NSigmas = cms.vint32(-1,1),
         OverallRange = cms.string("full5x5_r9<0.94&&abs(superCluster.eta)>=1.5"),
         BinList = emptyBins,
         CorrectionFile = scalesAndSmearingsPrefix,
         ApplyCentralValue = cms.bool(False),
         UncertaintyBitMask = cms.string("011"),#cms.string("110"),
         ExaggerateShiftUp = cms.bool(False),
         Debug = cms.untracked.bool(False)
         )

MCScaleGain6EB_EGM = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScaleEGMTool"),
         MethodName = cms.string("FlashggDiPhotonFromPhoton"),
         Label = cms.string("MCScaleGain6EB"),
         NSigmas = cms.vint32(-1,1),
         OverallRange = cms.string("hasSwitchToGain6&&abs(superCluster.eta)<1.5"),
         BinList = emptyBins,
         CorrectionFile = scalesAndSmearingsPrefix,
         ApplyCentralValue = cms.bool(False),
         ExaggerateShiftUp = cms.bool(False),
         UncertaintyBitMask = cms.string("100"),#cms.string("001"), # this should be a bit mask, don't know how to make it in python now
         Debug = cms.untracked.bool(False)
         )

MCScaleGain1EB_EGM = cms.PSet( PhotonMethodName = cms.string("FlashggPhotonScaleEGMTool"),
         MethodName = cms.string("FlashggDiPhotonFromPhoton"),
         Label = cms.string("MCScaleGain1EB"),
         NSigmas = cms.vint32(-1,1),
         OverallRange = cms.string("hasSwitchToGain1&&abs(superCluster.eta)<1.5"),
         BinList = emptyBins,
         CorrectionFile = scalesAndSmearingsPrefix,
         ApplyCentralValue = cms.bool(False),
         ExaggerateShiftUp = cms.bool(False),
         UncertaintyBitMask = cms.string("100"),#cms.string("001"), # this should be a bit mask, don't know how to make it in python now
         Debug = cms.untracked.bool(False)
         )

