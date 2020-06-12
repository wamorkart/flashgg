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

# Lowmass Preselection SF and uncertainties Update : Prasant (18122019)
preselBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.0 ) , upBounds = cms.vdouble( 1.5, 0.85 ) , values = cms.vdouble( 1.0050 ) , uncertainties = cms.vdouble( 0.0333 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85 ), upBounds = cms.vdouble( 1.5, 999.0 ) , values = cms.vdouble( 0.9934 ) , uncertainties = cms.vdouble( 0.0139 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.0 ) , upBounds = cms.vdouble( 6.0, 0.9) , values = cms.vdouble( 0.9946 ) , uncertainties = cms.vdouble( 0.0155 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.9 ) , upBounds = cms.vdouble( 6.0, 999.0 ) , values = cms.vdouble( 0.9977 ) , uncertainties = cms.vdouble( 0.0070 )  )
        )
    )


# JTao: Low mass case
electronVetoBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.00 ) , upBounds = cms.vdouble( 1.5, 0.85 ) , values = cms.vdouble( 0.9610 ) , uncertainties = cms.vdouble( 0.0054 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85 ) , upBounds = cms.vdouble( 1.5, 999. ) , values = cms.vdouble( 0.9754 ) , uncertainties = cms.vdouble( 0.0017 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.00 ) , upBounds = cms.vdouble( 6.0, 0.90 ) , values = cms.vdouble( 0.9161 ) , uncertainties = cms.vdouble( 0.0293 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.90 ) , upBounds = cms.vdouble( 6.0, 999. ) , values = cms.vdouble( 0.9202 ) , uncertainties = cms.vdouble( 0.0078 )  )
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


	#pt binning seeded leg: 0.,  33, 35., 37., 40., 45., 50., 60., 70., 90., 300	
        
	###BARREL	
	#0 <eta < 1.5, R9<0.85 ==> No low-R9 photons for te seeded leg/leading photon in the barrel
        cms.PSet(lowBounds = cms.vdouble(0.0,0.,0.), upBounds = cms.vdouble(0.85,1.5,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),


        #0 <eta < 1.5, 0.85<R9<0.88 ==> R9 turn-on bin - Barrel
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,0.), upBounds = cms.vdouble(0.88,1.5,33), values = cms.vdouble(0.697427), uncertainties = cms.vdouble(0.00137832 , 0.00137864)),  
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,33), upBounds = cms.vdouble(0.88,1.5,35.), values = cms.vdouble(0.882715), uncertainties = cms.vdouble(0.00107557 , 0.00107622)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,35.), upBounds = cms.vdouble(0.88,1.5,37.), values = cms.vdouble(0.894138), uncertainties = cms.vdouble(0.00100788 , 0.00100859)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,37.), upBounds = cms.vdouble(0.88,1.5,40.), values = cms.vdouble(0.902827), uncertainties = cms.vdouble(0.00100849 , 0.00100921)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,40.), upBounds = cms.vdouble(0.88,1.5,45.), values = cms.vdouble(0.912114), uncertainties = cms.vdouble(0.00100914 , 0.00100987)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,45.), upBounds = cms.vdouble(0.88,1.5,50.), values = cms.vdouble(0.915718), uncertainties = cms.vdouble(0.00100939 , 0.00101013)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,50.), upBounds = cms.vdouble(0.88,1.5,60.), values = cms.vdouble(0.920406), uncertainties = cms.vdouble(0.00100972 , 0.00101047)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,60.), upBounds = cms.vdouble(0.88,1.5,70.), values = cms.vdouble(0.924396), uncertainties = cms.vdouble(0.00133935 , 0.00133992)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,70.), upBounds = cms.vdouble(0.88,1.5,90.), values = cms.vdouble(0.92879), uncertainties = cms.vdouble(0.00155606 , 0.00155655)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,90.), upBounds = cms.vdouble(0.88,1.5,999999.), values = cms.vdouble(0.939501), uncertainties = cms.vdouble(0.00167275 , 0.00167322)),


        #0 <eta < 1.5, 0.88<R9<999 
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,0.), upBounds = cms.vdouble(999,1.5,33), values = cms.vdouble(0.764902), uncertainties = cms.vdouble(0.000999564 , 0.00100009)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,33), upBounds = cms.vdouble(999,1.5,35.), values = cms.vdouble(0.930524), uncertainties = cms.vdouble(0.00101045 , 0.00101121)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,35.), upBounds = cms.vdouble(999,1.5,37.), values = cms.vdouble(0.938207), uncertainties = cms.vdouble(0.001011 , 0.00101178)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,37.), upBounds = cms.vdouble(999,1.5,40.), values = cms.vdouble(0.943564), uncertainties = cms.vdouble(0.00101139 , 0.00101217)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,40.), upBounds = cms.vdouble(999,1.5,45.), values = cms.vdouble(0.947708), uncertainties = cms.vdouble(0.00101169 , 0.00101248)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,45.), upBounds = cms.vdouble(999,1.5,50.), values = cms.vdouble(0.949471), uncertainties = cms.vdouble(0.00101182 , 0.00101261)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,50.), upBounds = cms.vdouble(999,1.5,60.), values = cms.vdouble(0.949001), uncertainties = cms.vdouble(0.00101178 , 0.00101258)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,60.), upBounds = cms.vdouble(999,1.5,70.), values = cms.vdouble(0.95082), uncertainties = cms.vdouble(0.00101192 , 0.00101271)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,70.), upBounds = cms.vdouble(999,1.5,90.), values = cms.vdouble(0.957201), uncertainties = cms.vdouble(0.00101239 , 0.00101319)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,90.), upBounds = cms.vdouble(999,1.5,999999.), values = cms.vdouble(0.963203), uncertainties = cms.vdouble(0.00101283 , 0.00101364)),


	###ENDCAPS
 	#1.5 <eta < 3, R9<0.9 ==> No low-R9 photons for te seeded leg/leading photon in the Endcaps
        cms.PSet(lowBounds = cms.vdouble(0.,1.5,0.), upBounds = cms.vdouble(0.9,3.,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),

        #1.5 <eta < 3, 0.9<R9<0.93 ==> R9 turn-on bin - Endcaps
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,0), upBounds = cms.vdouble(0.93,3.,33.), values = cms.vdouble(0.533179), uncertainties = cms.vdouble(0.00166349 , 0.00166417)),   
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,33.), upBounds = cms.vdouble(0.93,3.,35.), values = cms.vdouble(0.78167), uncertainties = cms.vdouble(0.0015198 , 0.0015214)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,35.), upBounds = cms.vdouble(0.93,3.,37.), values = cms.vdouble(0.853343), uncertainties = cms.vdouble(0.00121737 , 0.00121976)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,37.), upBounds = cms.vdouble(0.93,3.,40.), values = cms.vdouble(0.882869), uncertainties = cms.vdouble(0.00104966 , 0.00105263)),       
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,40.), upBounds = cms.vdouble(0.93,3.,45.), values = cms.vdouble(0.901054), uncertainties = cms.vdouble(0.00105239 , 0.00105548)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,45.), upBounds = cms.vdouble(0.93,3.,50.), values = cms.vdouble(0.914075), uncertainties = cms.vdouble(0.00105438 , 0.00105755)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,50.), upBounds = cms.vdouble(0.93,3.,60.), values = cms.vdouble(0.921084), uncertainties = cms.vdouble(0.00105546 , 0.00105868)),  
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,60.), upBounds = cms.vdouble(0.93,3.,70.), values = cms.vdouble(0.933343), uncertainties = cms.vdouble(0.00161384 , 0.001616)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,70.), upBounds = cms.vdouble(0.93,3.,90.), values = cms.vdouble(0.943441), uncertainties = cms.vdouble(0.00177307 , 0.00177507)),  
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,90.), upBounds = cms.vdouble(0.93,3.,9999.), values = cms.vdouble(0.958013), uncertainties = cms.vdouble(0.00183042 , 0.00183242)), 
        
        #1.5 <eta < 3, 0.93<R9<999 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,0), upBounds = cms.vdouble(999,3.,33.), values = cms.vdouble(0.587245), uncertainties = cms.vdouble(0.00129654 , 0.0012976)),  
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,33.), upBounds = cms.vdouble(999,3.,35.), values = cms.vdouble(0.860594), uncertainties = cms.vdouble(0.00104777 , 0.0010506)),  
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,35.), upBounds = cms.vdouble(999,3.,37.), values = cms.vdouble(0.904747), uncertainties = cms.vdouble(0.00105296 , 0.00105606)), 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,37.), upBounds = cms.vdouble(999,3.,40.), values = cms.vdouble(0.926744), uncertainties = cms.vdouble(0.00105634 , 0.00105959)),       
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,40.), upBounds = cms.vdouble(999,3.,45.), values = cms.vdouble(0.942814), uncertainties = cms.vdouble(0.00105886 , 0.00106221)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,45.), upBounds = cms.vdouble(999,3.,50.), values = cms.vdouble(0.951324), uncertainties = cms.vdouble(0.00106021 , 0.00106362)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,50.), upBounds = cms.vdouble(999,3.,60.), values = cms.vdouble(0.953806), uncertainties = cms.vdouble(0.0010606 , 0.00106403)), 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,60.), upBounds = cms.vdouble(999,3.,70.), values = cms.vdouble(0.958042), uncertainties = cms.vdouble(0.00106809 , 0.00107152)), 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,70.), upBounds = cms.vdouble(999,3.,90.), values = cms.vdouble(0.96411), uncertainties = cms.vdouble(0.00108331 , 0.00108674)), 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,90.), upBounds = cms.vdouble(999,3.,9999.), values = cms.vdouble(0.966957), uncertainties = cms.vdouble(0.00123276 , 0.00123579)),

        #eta > 3. not used 
        cms.PSet(lowBounds = cms.vdouble(0.0,3.,0.), upBounds = cms.vdouble(999.,999.,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.))
        )
    )


subleadTriggerScaleBins = cms.PSet(
    variables = cms.vstring("full5x5_r9","abs(superCluster.eta)","pt"),
    bins = cms.VPSet(

	#pt binning for unseeded leg 0., 25., 28., 31., 40., 70., 300.	
        

	###BARREL
	#0. <eta < 1.5, R9<0.5 ==> No photons with r9 < 0.5 in the barrel
        cms.PSet(lowBounds = cms.vdouble(0.0,0,0.), upBounds = cms.vdouble(0.5,1.5,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),

        #eta < 1.5, 0.5<R9<0.53 ==> R9 Turn-on bin - low R9 - Barrel 
        cms.PSet(lowBounds = cms.vdouble(0.50,0,0.), upBounds = cms.vdouble(0.53,1.5,25.), values = cms.vdouble(0.762829), uncertainties = cms.vdouble(0.0270595 , 0.0270596)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0,25.), upBounds = cms.vdouble(0.53,1.5,28.), values = cms.vdouble(0.881282), uncertainties = cms.vdouble(0.0435529 , 0.043553)),   
        cms.PSet(lowBounds = cms.vdouble(0.50,0,28.), upBounds = cms.vdouble(0.53,1.5,31.), values = cms.vdouble(0.926675), uncertainties = cms.vdouble(0.0414056 , 0.0414057)),  
        cms.PSet(lowBounds = cms.vdouble(0.50,0,31.), upBounds = cms.vdouble(0.53,1.5,40.), values = cms.vdouble(0.902538), uncertainties = cms.vdouble(0.0226052 , 0.0226054)),  
        cms.PSet(lowBounds = cms.vdouble(0.50,0,40.), upBounds = cms.vdouble(0.53,1.5,70.), values = cms.vdouble(0.881184), uncertainties = cms.vdouble(0.0254469 , 0.0254471)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,70.), upBounds = cms.vdouble(0.53,1.5,999999.), values = cms.vdouble(0.821656), uncertainties = cms.vdouble(0.178344 , 0.248831)),

        #eta < 1.5, 0.53<R9<0.85 ==> low R9 - Barrel
        cms.PSet(lowBounds = cms.vdouble(0.53,0,0.), upBounds = cms.vdouble(0.85,1.5,25.), values = cms.vdouble(0.841898), uncertainties = cms.vdouble(0.0115231 , 0.0115234)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0,25.), upBounds = cms.vdouble(0.85,1.5,28.), values = cms.vdouble(0.925735), uncertainties = cms.vdouble(0.014695 , 0.0146953)),  
        cms.PSet(lowBounds = cms.vdouble(0.53,0,28.), upBounds = cms.vdouble(0.85,1.5,31.), values = cms.vdouble(0.929883), uncertainties = cms.vdouble(0.0138927 , 0.013893)),   
        cms.PSet(lowBounds = cms.vdouble(0.53,0,31.), upBounds = cms.vdouble(0.85,1.5,40.), values = cms.vdouble(0.935106), uncertainties = cms.vdouble(0.0112811 , 0.0112815)), 
        cms.PSet(lowBounds = cms.vdouble(0.53,0,40.), upBounds = cms.vdouble(0.85,1.5,70.), values = cms.vdouble(0.940301), uncertainties = cms.vdouble(0.0113313 , 0.0113318)), 
        cms.PSet(lowBounds = cms.vdouble(0.53,0,70.), upBounds = cms.vdouble(0.85,1.5,999999.), values = cms.vdouble(0.938715), uncertainties = cms.vdouble(0.0255649 , 0.0255651)),  

 	#eta < 1.5, 0.85<R9<999 ==> High R9 - Barrel	       
        cms.PSet(lowBounds = cms.vdouble(0.85,0,0.), upBounds = cms.vdouble(999,1.5,25.), values = cms.vdouble(0.931975), uncertainties = cms.vdouble(0.0152747 , 0.0152748)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,25.), upBounds = cms.vdouble(999,1.5,28.), values = cms.vdouble(0.959115), uncertainties = cms.vdouble(0.0171234 , 0.0171235)),   
        cms.PSet(lowBounds = cms.vdouble(0.85,0,28.), upBounds = cms.vdouble(999,1.5,31.), values = cms.vdouble(0.96881), uncertainties = cms.vdouble(0.0129254 , 0.0156261)),  
        cms.PSet(lowBounds = cms.vdouble(0.85,0,31.), upBounds = cms.vdouble(999,1.5,40.), values = cms.vdouble(0.975143), uncertainties = cms.vdouble(0.0114152 , 0.012184)),  
        cms.PSet(lowBounds = cms.vdouble(0.85,0,40.), upBounds = cms.vdouble(999,1.5,70.), values = cms.vdouble(0.970157), uncertainties = cms.vdouble(0.011705 , 0.0117051)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0,70.), upBounds = cms.vdouble(999,1.5,999999.), values = cms.vdouble(0.975143), uncertainties = cms.vdouble(0.0114152 , 0.0168243)), 


	###ENDCAPS	
	#1.5 <eta < 3., R9<0.8 ==> No photons with r9 < 0.8 in the endcaps
        cms.PSet(lowBounds = cms.vdouble(0.0,1.5,0.), upBounds = cms.vdouble(0.8,3.,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),

        #1.5 <eta < 3, 0.8<R9<0.83 ==> R9 Turn-on bin - low R9 - Endcaps
        cms.PSet(lowBounds = cms.vdouble(0.8,1.5,0.), upBounds = cms.vdouble(0.83,3.,25.), values = cms.vdouble(0.709933), uncertainties = cms.vdouble(0.0303323 , 0.0303343)),
        cms.PSet(lowBounds = cms.vdouble(0.8,1.5,25.), upBounds = cms.vdouble(0.83,3.,28.), values = cms.vdouble(0.815994), uncertainties = cms.vdouble(0.0450792 , 0.0450809)), 
        cms.PSet(lowBounds = cms.vdouble(0.8,1.5,28.), upBounds = cms.vdouble(0.83,3.,31.), values = cms.vdouble(0.832926), uncertainties = cms.vdouble(0.0387947 , 0.0387969)),
        cms.PSet(lowBounds = cms.vdouble(0.8,1.5,31.), upBounds = cms.vdouble(0.83,3.,40.), values = cms.vdouble(0.850107), uncertainties = cms.vdouble(0.021142 , 0.0211461)),
        cms.PSet(lowBounds = cms.vdouble(0.8,1.5,40.), upBounds = cms.vdouble(0.83,3.,70.), values = cms.vdouble(0.882189), uncertainties = cms.vdouble(0.0240573 , 0.0240612)),
        cms.PSet(lowBounds = cms.vdouble(0.8,1.5,70.), upBounds = cms.vdouble(0.83,3.,999999.), values = cms.vdouble(0.695921), uncertainties = cms.vdouble(0.0863995 , 0.0864002)),

        ##1.5 <eta < 3, 0.83<R9<0.9 ==> low R9 - Endcaps
	cms.PSet(lowBounds = cms.vdouble(0.83,1.5,0.), upBounds = cms.vdouble(0.9,3.,25.), values = cms.vdouble(0.803149), uncertainties = cms.vdouble(0.0184513 , 0.0184555)),
        cms.PSet(lowBounds = cms.vdouble(0.83,1.5,25.), upBounds = cms.vdouble(0.9,3.,28.), values = cms.vdouble(0.934628), uncertainties = cms.vdouble(0.0249662 , 0.0249704)),
        cms.PSet(lowBounds = cms.vdouble(0.83,1.5,28.), upBounds = cms.vdouble(0.9,3.,31.), values = cms.vdouble(0.932669), uncertainties = cms.vdouble(0.0229007 , 0.0229053)),
        cms.PSet(lowBounds = cms.vdouble(0.83,1.5,31.), upBounds = cms.vdouble(0.9,3.,40.), values = cms.vdouble(0.939581), uncertainties = cms.vdouble(0.0147877 , 0.0147949)),
        cms.PSet(lowBounds = cms.vdouble(0.83,1.5,40.), upBounds = cms.vdouble(0.9,3.,70.), values = cms.vdouble(0.933654), uncertainties = cms.vdouble(0.0148203 , 0.0148274)), 
        cms.PSet(lowBounds = cms.vdouble(0.83,1.5,70.), upBounds = cms.vdouble(0.9,3.,999999.), values = cms.vdouble(0.965191), uncertainties = cms.vdouble(0.034809 , 0.0574877)), 

	##1.5 <eta < 3, 0.9<R9<999 ==> high R9 - Endcaps	
 	cms.PSet(lowBounds = cms.vdouble(0.9,1.5,0.), upBounds = cms.vdouble(999,3.,25.), values = cms.vdouble(0.842708), uncertainties = cms.vdouble(0.0177746 , 0.017775)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,25.), upBounds = cms.vdouble(999,3.,28.), values = cms.vdouble(0.945462), uncertainties = cms.vdouble(0.0233602 , 0.0233606)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,28.), upBounds = cms.vdouble(999,3.,31.), values = cms.vdouble(0.962701), uncertainties = cms.vdouble(0.0188846 , 0.0205055)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,31.), upBounds = cms.vdouble(999,3.,40.), values = cms.vdouble(0.965889), uncertainties = cms.vdouble(0.0139393 , 0.01394)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,40.), upBounds = cms.vdouble(999,3.,70.), values = cms.vdouble(0.969626), uncertainties = cms.vdouble(0.0134006 , 0.0134013)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,70.), upBounds = cms.vdouble(999,3.,999999.), values = cms.vdouble(0.977952), uncertainties = cms.vdouble(0.0114919 , 0.031207)),


         #eta > 3. not used
        cms.PSet(lowBounds = cms.vdouble(0.0,3.,0), upBounds = cms.vdouble(999,999,999999), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.))
        )
    )

# Lowmass  LooseIDMVA cut SF and uncertainities Update : Prasant (18/12/2019)
looseMvaBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.0  ) , upBounds = cms.vdouble( 1.5, 0.85  ) , values = cms.vdouble( 0.9950 ) , uncertainties = cms.vdouble( 0.0026 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85 ) , upBounds = cms.vdouble( 1.5, 999.0 ) , values = cms.vdouble( 0.9979 ) , uncertainties = cms.vdouble( 0.0017 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.0  ) , upBounds = cms.vdouble( 6.0, 0.9   ) , values = cms.vdouble( 1.0007 ) , uncertainties = cms.vdouble( 0.0022 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.9  ) , upBounds = cms.vdouble( 6.0, 999.0 ) , values = cms.vdouble( 1.0001 ) , uncertainties = cms.vdouble( 0.0012 ) ) ) )

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


scalesAndSmearingsPrefix = cms.string("EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc")
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
          CorrectionFile = cms.FileInPath("flashgg/Systematics/data/SystematicsIDMVA_LegRunII_v1_2017.root"),
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

