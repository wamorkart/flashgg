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

# from Arnab 
preselBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.0 ) , upBounds = cms.vdouble( 1.5, 0.85 ) , values = cms.vdouble( 1.0399 ) , uncertainties = cms.vdouble( 0.0305 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85 ), upBounds = cms.vdouble( 1.5, 999. ) , values = cms.vdouble( 1.0152 ) , uncertainties = cms.vdouble( 0.0193 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.0 ) , upBounds = cms.vdouble( 6.0, 0.9  ) , values = cms.vdouble( 1.0 ) , uncertainties = cms.vdouble( 0.0 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.9 ) , upBounds = cms.vdouble( 6.0, 999. ) , values = cms.vdouble( 0.9717 ) , uncertainties = cms.vdouble( 0.0335 )  )
        )
    )


# slide  ...
electronVetoBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.00 ) , upBounds = cms.vdouble( 1.5, 0.85 ) , values = cms.vdouble( 1.0281 ) , uncertainties = cms.vdouble( 0.0123 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85 ) , upBounds = cms.vdouble( 1.5, 999. ) , values = cms.vdouble( 1.0090 ) , uncertainties = cms.vdouble( 0.0059 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.00 ) , upBounds = cms.vdouble( 6.0, 0.90 ) , values = cms.vdouble( 1.0 ) , uncertainties = cms.vdouble( 0.0 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.90 ) , upBounds = cms.vdouble( 6.0, 999. ) , values = cms.vdouble( 1.0266 ) , uncertainties = cms.vdouble( 0.0162 )  )
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


#with full 2016  dataset from Linda
leadTriggerScaleBins = cms.PSet(
    variables = cms.vstring("full5x5_r9","abs(superCluster.eta)","pt"),
    bins = cms.VPSet(

        #0 <eta < 1.5, R9<0.50 ==> Low mass: in fact not used
        cms.PSet(lowBounds = cms.vdouble(0.0,0.,0.), upBounds = cms.vdouble(0.50,1.5,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),

        #0 <eta < 1.5, 0.50<R9<0.53 ==> Low mass: turn-on bin ==> "OR" efficiencies 
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,0.), upBounds = cms.vdouble(0.53,1.5,33.3333), values = cms.vdouble(0.594793), uncertainties = cms.vdouble(0.00612789,0.00609896)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,33.3333), upBounds = cms.vdouble(0.53,1.5,35.), values = cms.vdouble(0.662945), uncertainties = cms.vdouble(0.00783236,0.00774567)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,35.), upBounds = cms.vdouble(0.53,1.5,40.), values = cms.vdouble(0.692324), uncertainties = cms.vdouble(0.00398642,0.00395829)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,40.), upBounds = cms.vdouble(0.53,1.5,45.), values = cms.vdouble(0.708503), uncertainties = cms.vdouble(0.00343288,0.00340952)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,45.), upBounds = cms.vdouble(0.53,1.5,50.), values = cms.vdouble(0.732292), uncertainties = cms.vdouble(0.00438748,0.00434297)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,50.), upBounds = cms.vdouble(0.53,1.5,60.), values = cms.vdouble(0.744355), uncertainties = cms.vdouble(0.00646887,0.0063652)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,60.), upBounds = cms.vdouble(0.53,1.5,70.), values = cms.vdouble(0.792566), uncertainties = cms.vdouble(0.0150806,0.0143441)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,70.), upBounds = cms.vdouble(0.53,1.5,90.), values = cms.vdouble(0.808), uncertainties = cms.vdouble(0.0283497,0.025675)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0.,90.), upBounds = cms.vdouble(0.53,1.5,999999.), values = cms.vdouble(0.783784), uncertainties = cms.vdouble(0.0896274,0.0714411)),


        #eta < 1.5, 0.53<R9<0.85 ==> "OR" efficiencies
        cms.PSet(lowBounds = cms.vdouble(0.53,0.,0.), upBounds = cms.vdouble(0.85,1.5,33.3333), values = cms.vdouble(0.880773), uncertainties = cms.vdouble(0.00100007,0.000992926)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0.,33.3333), upBounds = cms.vdouble(0.85,1.5,35.), values = cms.vdouble(0.951792), uncertainties = cms.vdouble(0.000842411,0.000828832)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0.,35.), upBounds = cms.vdouble(0.85,1.5,40.), values = cms.vdouble(0.961186), uncertainties = cms.vdouble(0.000371599,0.000368245)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0.,40.), upBounds = cms.vdouble(0.85,1.5,45.), values = cms.vdouble(0.967805), uncertainties = cms.vdouble(0.000267569,0.000265453)), 
        cms.PSet(lowBounds = cms.vdouble(0.53,0.,45.), upBounds = cms.vdouble(0.85,1.5,50.), values = cms.vdouble(0.971748), uncertainties = cms.vdouble(0.00030314,0.000300042)), 
        cms.PSet(lowBounds = cms.vdouble(0.53,0.,50.), upBounds = cms.vdouble(0.85,1.5,60.), values = cms.vdouble(0.97485), uncertainties = cms.vdouble(0.000389931,0.000384194)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0.,60.), upBounds = cms.vdouble(0.85,1.5,70.), values = cms.vdouble(0.979271), uncertainties = cms.vdouble(0.000769622,0.000743186)), 
        cms.PSet(lowBounds = cms.vdouble(0.53,0.,70.), upBounds = cms.vdouble(0.85,1.5,90.), values = cms.vdouble(0.983438), uncertainties = cms.vdouble(0.00107967,0.00101667)), 
        cms.PSet(lowBounds = cms.vdouble(0.53,0.,90.), upBounds = cms.vdouble(0.85,1.5,999999.), values = cms.vdouble(0.990463), uncertainties = cms.vdouble(0.00143767,0.00126133)),


        #eta < 1.5, 0.85<R9<0.88 ==> "AND" efficiencies
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,0.), upBounds = cms.vdouble(0.88,1.5,33.3333), values = cms.vdouble(0.480261), uncertainties = cms.vdouble(0.0047909,0.00479449)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,33.3333), upBounds = cms.vdouble(0.88,1.5,35.), values = cms.vdouble(0.544853), uncertainties = cms.vdouble(0.00595334,0.00594074)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,35.), upBounds = cms.vdouble(0.88,1.5,40.), values = cms.vdouble(0.564528), uncertainties = cms.vdouble(0.00285718,0.00285294)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,40.), upBounds = cms.vdouble(0.88,1.5,45.), values = cms.vdouble(0.584602), uncertainties = cms.vdouble(0.00218262,0.00217934)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,45.), upBounds = cms.vdouble(0.88,1.5,50.), values = cms.vdouble(0.610664), uncertainties = cms.vdouble(0.00255131,0.00254531)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,50.), upBounds = cms.vdouble(0.88,1.5,60.), values = cms.vdouble(0.638541), uncertainties = cms.vdouble(0.00321898,0.00320672)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,60.), upBounds = cms.vdouble(0.88,1.5,70.), values = cms.vdouble(0.679178), uncertainties = cms.vdouble(0.0061075,0.00604779)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,70.), upBounds = cms.vdouble(0.88,1.5,90.), values = cms.vdouble(0.719895), uncertainties = cms.vdouble(0.0079413,0.00780903)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,90.), upBounds = cms.vdouble(0.88,1.5,9999.), values = cms.vdouble(0.782112), uncertainties = cms.vdouble(0.00942959,0.0091521)),


         #eta < 1.5, R9>0.88 ==> "AND" efficiencies
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,0.), upBounds = cms.vdouble(999.,1.5,33.3333), values = cms.vdouble(0.911612), uncertainties = cms.vdouble(0.000551594,0.000548525)), 
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,33.3333), upBounds = cms.vdouble(999.,1.5,35.), values = cms.vdouble(0.964408), uncertainties = cms.vdouble(0.000451792,0.00044639)), 
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,35.), upBounds = cms.vdouble(999.,1.5,40.), values = cms.vdouble(0.973082), uncertainties = cms.vdouble(0.000187761,0.000186505)),  
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,40.), upBounds = cms.vdouble(999.,1.5,45.), values = cms.vdouble(0.978494), uncertainties = cms.vdouble(0.000126388,0.000125671)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,45.), upBounds = cms.vdouble(999.,1.5,50.), values = cms.vdouble(0.981905), uncertainties = cms.vdouble(0.000137401,0.000136392)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,50.), upBounds = cms.vdouble(999.,1.5,60.), values = cms.vdouble(0.982553), uncertainties = cms.vdouble(0.000174033,0.000172359)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,60.), upBounds = cms.vdouble(999.,1.5,70.), values = cms.vdouble(0.983744), uncertainties = cms.vdouble(0.000329198,0.000322855)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,70.), upBounds = cms.vdouble(999.,1.5,90.), values = cms.vdouble(0.986766), uncertainties = cms.vdouble(0.000405361,0.000393681)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0.,90.), upBounds = cms.vdouble(999.,1.5,9999.), values = cms.vdouble(0.990753), uncertainties = cms.vdouble(0.000407389,0.000390745)),

       #1.5 <eta < 3, 0.9<R9<0.93 ==> "AND" efficiencies (turn-on bin)
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,0), upBounds = cms.vdouble(0.93,3,33.3333), values = cms.vdouble(0.595211), uncertainties = cms.vdouble(0.00560117,0.00557685)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,33.3333), upBounds = cms.vdouble(0.93,3.,35.), values = cms.vdouble(0.702697), uncertainties = cms.vdouble(0.00674031,0.00665491)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,35.), upBounds = cms.vdouble(0.93,3.,40.), values = cms.vdouble(0.734526), uncertainties = cms.vdouble(0.00328452,0.00325908)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,40.), upBounds = cms.vdouble(0.93,3.,45.), values = cms.vdouble(0.760539), uncertainties = cms.vdouble(0.00252922,0.00251124)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,45.), upBounds = cms.vdouble(0.93,3.,50.), values = cms.vdouble(0.778439), uncertainties = cms.vdouble(0.00304564,0.00301634)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,50.), upBounds = cms.vdouble(0.93,3.,60.), values = cms.vdouble(0.791388), uncertainties = cms.vdouble(0.00380274,0.0037531)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,60.), upBounds = cms.vdouble(0.93,3.,70.), values = cms.vdouble(0.81457), uncertainties = cms.vdouble(0.00719971,0.00699529)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,70.), upBounds = cms.vdouble(0.93,3.,90.), values = cms.vdouble(0.837738), uncertainties = cms.vdouble(0.00897237,0.00860159)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,90.), upBounds = cms.vdouble(0.93,3.,9999.), values = cms.vdouble(0.884554), uncertainties = cms.vdouble(0.0097871,0.00914326)),
        
        #1.5 <eta < 3, R9>0.93 ==> "AND" efficiencies
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,0.), upBounds = cms.vdouble(999.,3.,33.3333), values = cms.vdouble(0.853852), uncertainties = cms.vdouble(0.00154649,0.00153317)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,33.3333), upBounds = cms.vdouble(999.,3.,35.), values = cms.vdouble(0.950491), uncertainties = cms.vdouble(0.00127802,0.00124798)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,35.), upBounds = cms.vdouble(999.,3.,40.), values = cms.vdouble(0.965209), uncertainties = cms.vdouble(0.000558641,0.000550223)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,40.), upBounds = cms.vdouble(999.,3.,45.), values = cms.vdouble(0.970964), uncertainties = cms.vdouble(0.000416357,0.000410705)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,45.), upBounds = cms.vdouble(999.,3.,50.), values = cms.vdouble(0.973743), uncertainties = cms.vdouble(0.000499808,0.000490835)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,50.), upBounds = cms.vdouble(999.,3.,60.), values = cms.vdouble(0.977689), uncertainties = cms.vdouble(0.000597825,0.000582818)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,60.), upBounds = cms.vdouble(999.,3.,70.), values = cms.vdouble(0.981832), uncertainties = cms.vdouble(0.00107608,0.00101868)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,70.), upBounds = cms.vdouble(999.,3.,90.), values = cms.vdouble(0.989488), uncertainties = cms.vdouble(0.00110048,0.0010017)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,90.), upBounds = cms.vdouble(999,3.,9999.), values = cms.vdouble(0.993663), uncertainties = cms.vdouble(0.00105488,0.00091404)),

        #1.5 <eta < 3, R9<0.9 ==> Low mass : in fact not used
        cms.PSet(lowBounds = cms.vdouble(0.,1.5,0.), upBounds = cms.vdouble(0.9,3.,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),

        #eta > 3. ==> Low mass : in fact not used
        cms.PSet(lowBounds = cms.vdouble(0.0,3.,0.), upBounds = cms.vdouble(999.,999.,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.))
        )
    )

leadTriggerScaleBinsEBhiR9OR= cms.PSet(
    variables = cms.vstring("full5x5_r9","abs(superCluster.eta)","pt"),
    bins = cms.VPSet(
        #everything a part from eta < 1.5, R9>0.85 ==> "OR" efficiency is 0, but needed because the preselection is applied after 
        #0 <eta < 1.5, R9<0.85 ==> Low mass: in fact not used
        cms.PSet(lowBounds = cms.vdouble(0.0,0.,0.), upBounds = cms.vdouble(0.85,1.5,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),

        #eta < 1.5, R9 >0.85 (OR path has its turn-on at 0.50-0.53, no need to put a turn-on bin between 0.85-0.88)
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,0.), upBounds = cms.vdouble(999,1.5,33.3333), values = cms.vdouble(0.917161), uncertainties = cms.vdouble(0.000453874,0.000451638)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,33.3333), upBounds = cms.vdouble(999,1.5,35.), values = cms.vdouble(0.969697), uncertainties = cms.vdouble(0.000353493,0.00034958)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,35.), upBounds = cms.vdouble(999,1.5,40.), values = cms.vdouble(0.977014), uncertainties = cms.vdouble(0.000146821,0.000145917)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,40.), upBounds = cms.vdouble(999,1.5,45.), values = cms.vdouble(0.980632), uncertainties = cms.vdouble(0.000100407,9.99024e-05)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,45.), upBounds = cms.vdouble(999,1.5,50.), values = cms.vdouble(0.983649), uncertainties = cms.vdouble(0.00010883,0.000108127)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,50.), upBounds = cms.vdouble(999,1.5,60.), values = cms.vdouble(0.984346), uncertainties = cms.vdouble(0.00013642,0.000135269)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,60.), upBounds = cms.vdouble(999,1.5,70.), values = cms.vdouble(0.98573), uncertainties = cms.vdouble(0.000254647,0.0002503)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,70.), upBounds = cms.vdouble(999,1.5,90.), values = cms.vdouble(0.98907), uncertainties = cms.vdouble(0.00030857,0.000300329)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0.,90.), upBounds = cms.vdouble(999,1.5,9999.), values = cms.vdouble(0.992703), uncertainties = cms.vdouble(0.000311862,0.000299456)), 

        #eta >1.5 , R9>0 ==> Low mass : in fact not used
        cms.PSet(lowBounds = cms.vdouble(0.,1.5,0), upBounds = cms.vdouble(999,999,999999), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.))

        )
    )


subleadTriggerScaleBins = cms.PSet(
    variables = cms.vstring("full5x5_r9","abs(superCluster.eta)","pt"),
    bins = cms.VPSet(
        #0. <eta < 1.5, R9<0.5 ==> Low mass: in fact not used
        cms.PSet(lowBounds = cms.vdouble(0.0,0,0.), upBounds = cms.vdouble(0.5,1.5,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),

         #eta < 1.5, 0.5<R9<0.53 ==> Low mass: turn-on bin ==> "OR" efficiencies 
        cms.PSet(lowBounds = cms.vdouble(0.50,0,0.), upBounds = cms.vdouble(0.53,1.5,32.), values = cms.vdouble(0.641895), uncertainties = cms.vdouble(0.00565514,0.00563762)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0,32.), upBounds = cms.vdouble(0.53,1.5,34.), values = cms.vdouble(0.682347), uncertainties = cms.vdouble(0.00769889,0.00761554)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0,34.), upBounds = cms.vdouble(0.53,1.5,36.), values = cms.vdouble(0.695856), uncertainties = cms.vdouble(0.00718414,0.00710465)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0,36.), upBounds = cms.vdouble(0.53,1.5,38.), values = cms.vdouble(0.695441), uncertainties = cms.vdouble(0.00682556,0.00675585)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,38.), upBounds = cms.vdouble(0.53,1.5,40.), values = cms.vdouble(0.714604), uncertainties = cms.vdouble(0.00652138,0.00645205)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,40.), upBounds = cms.vdouble(0.53,1.5,45.), values = cms.vdouble(0.719518), uncertainties = cms.vdouble(0.00443017,0.00441113)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,45.), upBounds = cms.vdouble(0.53,1.5,50.), values = cms.vdouble(0.74371), uncertainties = cms.vdouble(0.00516489,0.00512556)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,50.), upBounds = cms.vdouble(0.53,1.5,55.), values = cms.vdouble(0.754574), uncertainties = cms.vdouble(0.00792411,0.00777959)), 
        cms.PSet(lowBounds = cms.vdouble(0.50,0,55.), upBounds = cms.vdouble(0.53,1.5,65.), values = cms.vdouble(0.763377), uncertainties = cms.vdouble(0.0102859,0.0100044)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,65.), upBounds = cms.vdouble(0.53,1.5,70.), values = cms.vdouble(0.796185), uncertainties = cms.vdouble(0.0276435,0.025054)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,70.), upBounds = cms.vdouble(0.53,1.5,90.), values = cms.vdouble(0.802318), uncertainties = cms.vdouble(0.0278713,0.0251395)),
        cms.PSet(lowBounds = cms.vdouble(0.50,0,90.), upBounds = cms.vdouble(0.53,1.5,999999.), values = cms.vdouble(0.769673), uncertainties = cms.vdouble(0.0880707,0.0702263)),

        #eta < 1.5, 0.53<R9<0.85 ==> "OR" efficiencies
        cms.PSet(lowBounds = cms.vdouble(0.53,0,0.), upBounds = cms.vdouble(0.85,1.5,32.), values = cms.vdouble(0.963248), uncertainties = cms.vdouble(0.00581956,0.00581934)), 
        cms.PSet(lowBounds = cms.vdouble(0.53,0,32.), upBounds = cms.vdouble(0.85,1.5,34.), values = cms.vdouble(0.970865), uncertainties = cms.vdouble(0.0045656,0.00456476)),  
        cms.PSet(lowBounds = cms.vdouble(0.53,0,34.), upBounds = cms.vdouble(0.85,1.5,36.), values = cms.vdouble(0.973962), uncertainties = cms.vdouble(0.00412894,0.00412828)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0,36.), upBounds = cms.vdouble(0.85,1.5,38.), values = cms.vdouble(0.975698), uncertainties = cms.vdouble(0.0040755,0.00407505)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0,38.), upBounds = cms.vdouble(0.85,1.5,40.), values = cms.vdouble(0.97649), uncertainties = cms.vdouble(0.00402688,0.00402655)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0,40.), upBounds = cms.vdouble(0.85,1.5,45.), values = cms.vdouble(0.978155), uncertainties = cms.vdouble(0.00404844,0.0040484)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0,45.), upBounds = cms.vdouble(0.85,1.5,50.), values = cms.vdouble(0.978768), uncertainties = cms.vdouble(0.00400031,0.00400024)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0,50.), upBounds = cms.vdouble(0.85,1.5,55.), values = cms.vdouble(0.979168), uncertainties = cms.vdouble(0.00399705,0.00399673)),  
        cms.PSet(lowBounds = cms.vdouble(0.53,0,55.), upBounds = cms.vdouble(0.85,1.5,65.), values = cms.vdouble(0.979229), uncertainties = cms.vdouble(0.0039998,0.00399918)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0,65.), upBounds = cms.vdouble(0.85,1.5,70.), values = cms.vdouble(0.978977), uncertainties = cms.vdouble(0.00403552,0.00402446)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0,70.), upBounds = cms.vdouble(0.85,1.5,90.), values = cms.vdouble(0.978874), uncertainties = cms.vdouble(0.00404716,0.00403931)),
        cms.PSet(lowBounds = cms.vdouble(0.53,0,90.), upBounds = cms.vdouble(0.85,1.5,999999.), values = cms.vdouble(0.979177), uncertainties = cms.vdouble(0.00412128,0.00408621)),

       #eta < 1.5, 0.85< R9<0.88 ==> "AND" efficiencies
        cms.PSet(lowBounds = cms.vdouble(0.85,0,0), upBounds = cms.vdouble(0.88,1.5,32.), values = cms.vdouble(0.520538), uncertainties = cms.vdouble(0.00431893,0.00431599)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,32.), upBounds = cms.vdouble(0.88,1.5,34.), values = cms.vdouble(0.544783), uncertainties = cms.vdouble(0.00622536,0.0062114)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,34.), upBounds = cms.vdouble(0.88,1.5,36.), values = cms.vdouble(0.552613), uncertainties = cms.vdouble(0.00568908,0.00567638)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,36.), upBounds = cms.vdouble(0.88,1.5,38.), values = cms.vdouble(0.569686), uncertainties = cms.vdouble(0.00521165,0.00519895)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0,38.), upBounds = cms.vdouble(0.88,1.5,40.), values = cms.vdouble(0.57511), uncertainties = cms.vdouble(0.0048133,0.00480255)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,40.), upBounds = cms.vdouble(0.88,1.5,45.), values = cms.vdouble(0.586177), uncertainties = cms.vdouble(0.00323418,0.00323159)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,45.), upBounds = cms.vdouble(0.88,1.5,50.), values = cms.vdouble(0.612741), uncertainties = cms.vdouble(0.00356864,0.00356375)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0,50.), upBounds = cms.vdouble(0.88,1.5,55.), values = cms.vdouble(0.635098), uncertainties = cms.vdouble(0.00472043,0.004704)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,55.), upBounds = cms.vdouble(0.88,1.5,65.), values = cms.vdouble(0.656744), uncertainties = cms.vdouble(0.00528625,0.0052596)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0,65.), upBounds = cms.vdouble(0.88,1.5,70.), values = cms.vdouble(0.686169), uncertainties = cms.vdouble(0.0107875,0.0106001)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,70.), upBounds = cms.vdouble(0.88,1.5,90.), values = cms.vdouble(0.716083), uncertainties = cms.vdouble(0.00852188,0.00838884)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,90.), upBounds = cms.vdouble(0.88,1.5,999999.), values = cms.vdouble(0.771437), uncertainties = cms.vdouble(0.0100283,0.00975703)),  

        #eta < 1.5, R9>0.88 ==> "AND" efficiencies
        cms.PSet(lowBounds = cms.vdouble(0.88,0,0.), upBounds = cms.vdouble(999,1.5,32.), values = cms.vdouble(0.957582), uncertainties = cms.vdouble(0.00408733,0.00408722)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0,32.), upBounds = cms.vdouble(999,1.5,34.), values = cms.vdouble(0.96933), uncertainties = cms.vdouble(0.00398626,0.00398589)),  
        cms.PSet(lowBounds = cms.vdouble(0.88,0,34.), upBounds = cms.vdouble(999,1.5,36.), values = cms.vdouble(0.971038), uncertainties = cms.vdouble(0.00395399,0.00395373)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0,36.), upBounds = cms.vdouble(999,1.5,38.), values = cms.vdouble(0.973016), uncertainties = cms.vdouble(0.00395445,0.00395427)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0,38.), upBounds = cms.vdouble(999,1.5,40.), values = cms.vdouble(0.975327), uncertainties = cms.vdouble(0.00395841,0.00395829)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0,40.), upBounds = cms.vdouble(999,1.5,45.), values = cms.vdouble(0.977557), uncertainties = cms.vdouble(0.00396522,0.00396521)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0,45.), upBounds = cms.vdouble(999,1.5,50.), values = cms.vdouble(0.978469), uncertainties = cms.vdouble(0.00396684,0.00396682)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0,50.), upBounds = cms.vdouble(999,1.5,55.), values = cms.vdouble(0.979083), uncertainties = cms.vdouble(0.00397048,0.00397039)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0,55.), upBounds = cms.vdouble(999,1.5,65.), values = cms.vdouble(0.979269), uncertainties = cms.vdouble(0.00397217,0.00397202)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0,65.), upBounds = cms.vdouble(999,1.5,70.), values = cms.vdouble(0.979605), uncertainties = cms.vdouble(0.00398867,0.00398687)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0,70.), upBounds = cms.vdouble(999,1.5,90.), values = cms.vdouble(0.978773), uncertainties = cms.vdouble(0.00398014,0.00397916)),
        cms.PSet(lowBounds = cms.vdouble(0.88,0,90.), upBounds = cms.vdouble(999,1.5,999999.), values = cms.vdouble(0.980216), uncertainties = cms.vdouble(0.0039884,0.0039869)),

        #1.5 <eta < 3, 0.9<R9<0.93 ==> "AND" efficiencies
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,0.), upBounds = cms.vdouble(0.93,3.,32.), values = cms.vdouble(0.685109), uncertainties = cms.vdouble(0.014049,0.0140405)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,32.), upBounds = cms.vdouble(0.93,3.,34.), values = cms.vdouble(0.697586), uncertainties = cms.vdouble(0.0150491,0.0150125)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,34.), upBounds = cms.vdouble(0.93,3.,36.), values = cms.vdouble(0.708077), uncertainties = cms.vdouble(0.0150137,0.014983)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,36.), upBounds = cms.vdouble(0.93,3.,38.), values = cms.vdouble(0.71455), uncertainties = cms.vdouble(0.0149684,0.0149427)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,38.), upBounds = cms.vdouble(0.93,3.,40.), values = cms.vdouble(0.728204), uncertainties = cms.vdouble(0.0150431,0.0150221)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,40.), upBounds = cms.vdouble(0.93,3.,45.), values = cms.vdouble(0.7459), uncertainties = cms.vdouble(0.014774,0.0147706)),
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,45.), upBounds = cms.vdouble(0.93,3.,50.), values = cms.vdouble(0.761913), uncertainties = cms.vdouble(0.0151802,0.0151738)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,50.), upBounds = cms.vdouble(0.93,3.,55.), values = cms.vdouble(0.769585), uncertainties = cms.vdouble(0.0157417,0.0157184)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,55.), upBounds = cms.vdouble(0.93,3.,65.), values = cms.vdouble(0.784002), uncertainties = cms.vdouble(0.0162511,0.0162123)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,65.), upBounds = cms.vdouble(0.93,3.,70.), values = cms.vdouble(0.796937), uncertainties = cms.vdouble(0.0198298,0.0194583)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,70.), upBounds = cms.vdouble(0.93,3.,90.), values = cms.vdouble(0.815029), uncertainties = cms.vdouble(0.0183423,0.0181501)), 
        cms.PSet(lowBounds = cms.vdouble(0.9,1.5,90.), upBounds = cms.vdouble(0.93,3.,999999.), values = cms.vdouble(0.858506), uncertainties = cms.vdouble(0.0195013,0.0191683)), 

        #1.5 <eta < 3, R9>0.93 ==> "AND" efficiencies
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,0.), upBounds = cms.vdouble(999,3.,32), values = cms.vdouble(0.952053), uncertainties = cms.vdouble(0.0187433,0.0187431)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,32), upBounds = cms.vdouble(999,3.,34), values = cms.vdouble(0.960736), uncertainties = cms.vdouble(0.0187688,0.0187678)), 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,34), upBounds = cms.vdouble(999,3.,36), values = cms.vdouble(0.963011), uncertainties = cms.vdouble(0.0187951,0.0187943)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,36), upBounds = cms.vdouble(999,3.,38), values = cms.vdouble(0.962529), uncertainties = cms.vdouble(0.0187849,0.0187842)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,38), upBounds = cms.vdouble(999,3.,40), values = cms.vdouble(0.965013), uncertainties = cms.vdouble(0.0188305,0.01883)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,40), upBounds = cms.vdouble(999,3.,45),  values = cms.vdouble(0.965534), uncertainties = cms.vdouble(0.0188371,0.018837)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,45), upBounds = cms.vdouble(999,3.,50), values = cms.vdouble(0.966361), uncertainties = cms.vdouble(0.0188532,0.0188531)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,50), upBounds = cms.vdouble(999,3.,55), values = cms.vdouble(0.965761), uncertainties = cms.vdouble(0.0188452,0.0188446)), 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,55), upBounds = cms.vdouble(999,3.,65), values = cms.vdouble(0.966527), uncertainties = cms.vdouble(0.0188619,0.0188609)), 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,65), upBounds = cms.vdouble(999,3.,70), values = cms.vdouble(0.964668), uncertainties = cms.vdouble(0.0188709,0.0188586)), 
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,70), upBounds = cms.vdouble(999,3.,90), values = cms.vdouble(0.966477), uncertainties = cms.vdouble(0.0188786,0.0188733)),
        cms.PSet(lowBounds = cms.vdouble(0.93,1.5,90), upBounds = cms.vdouble(999,3.,999999), values = cms.vdouble(0.96697), uncertainties = cms.vdouble(0.0188979,0.0188896)),

         #1.5 <eta < 3, R9<0.9 ==> Low mass: in fact not used
        cms.PSet(lowBounds = cms.vdouble(0.,1.5,0.), upBounds = cms.vdouble(0.9,3.,999999), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),

        cms.PSet(lowBounds = cms.vdouble(0.0,3.,0), upBounds = cms.vdouble(999,999,999999), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.))
        )
    )

subleadTriggerScaleBinsEBhiR9OR = cms.PSet(
    variables = cms.vstring("full5x5_r9","abs(superCluster.eta)","pt"),
    bins = cms.VPSet(
        #0. < eta < 1.5, R9<0.85 ==> Low mass: in fact not used
        cms.PSet(lowBounds = cms.vdouble(0.0,0,0.), upBounds = cms.vdouble(0.85,1.5,999999.), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.)),

       #0. < eta < 1.5, R9>0.85 (OR path has its turn-on at 0.50-0.53, no need of a turn-on bin between 0.85-0.88)
        cms.PSet(lowBounds = cms.vdouble(0.85,0,0), upBounds = cms.vdouble(999,1.5,32.), values = cms.vdouble(0.973687), uncertainties = cms.vdouble(0.00507977,0.00507975)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0,32.), upBounds = cms.vdouble(999,1.5,34.), values = cms.vdouble(0.978161), uncertainties = cms.vdouble(0.00499689,0.00499683)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0,34.), upBounds = cms.vdouble(999,1.5,36.), values = cms.vdouble(0.979632), uncertainties = cms.vdouble(0.00498149,0.00498146)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0,36.), upBounds = cms.vdouble(999,1.5,38.), values = cms.vdouble(0.980581), uncertainties = cms.vdouble(0.00498197,0.00498194)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0,38.), upBounds = cms.vdouble(999,1.5,40.), values = cms.vdouble(0.981697), uncertainties = cms.vdouble(0.00498496,0.00498495)),  
        cms.PSet(lowBounds = cms.vdouble(0.85,0,40.), upBounds = cms.vdouble(999,1.5,45.), values = cms.vdouble(0.982505), uncertainties = cms.vdouble(0.00498911,0.00498911)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,45.), upBounds = cms.vdouble(999,1.5,50.), values = cms.vdouble(0.983331), uncertainties = cms.vdouble(0.00499181,0.00499181)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,50.), upBounds = cms.vdouble(999,1.5,55.), values = cms.vdouble(0.984008), uncertainties = cms.vdouble(0.00499507,0.00499506)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,55.), upBounds = cms.vdouble(999,1.5,65.), values = cms.vdouble(0.984377), uncertainties = cms.vdouble(0.00499693,0.00499692)),
        cms.PSet(lowBounds = cms.vdouble(0.85,0,65), upBounds = cms.vdouble(999,1.5,70), values = cms.vdouble(0.984618), uncertainties = cms.vdouble(0.00499824,0.00499815)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0,70), upBounds = cms.vdouble(999,1.5,90), values = cms.vdouble(0.984768), uncertainties = cms.vdouble(0.00499902,0.00499899)), 
        cms.PSet(lowBounds = cms.vdouble(0.85,0,90), upBounds = cms.vdouble(999,1.5,999999), values = cms.vdouble(0.98487), uncertainties = cms.vdouble(0.00499951,0.00499945)),

        #eta >1.5, R9>0. ==> Low mass: in fact not used
        cms.PSet(lowBounds = cms.vdouble(0.,1.5,0.), upBounds = cms.vdouble(999,999,999999), values = cms.vdouble(1.), uncertainties = cms.vdouble(0.,0.))

        )
    )


# from Arnab via                                                                                                                                
looseMvaBins = cms.PSet(
    variables = cms.vstring("abs(superCluster.eta)","full5x5_r9"),
    bins = cms.VPSet(
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.0  ) , upBounds = cms.vdouble( 1.5, 0.85  ) , values = cms.vdouble( 1.0144 ) , uncertainties = cms.vdouble( 0.0088 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 0.0, 0.85 ) , upBounds = cms.vdouble( 1.5, 999.0 ) , values = cms.vdouble( 1.0074 ) , uncertainties = cms.vdouble( 0.0034 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.0  ) , upBounds = cms.vdouble( 6.0, 0.9   ) , values = cms.vdouble( 1.0 ) , uncertainties = cms.vdouble( 0.0 )  ) ,
        cms.PSet( lowBounds = cms.vdouble( 1.5, 0.9  ) , upBounds = cms.vdouble( 6.0, 999.0 ) , values = cms.vdouble( 1.0331 ) , uncertainties = cms.vdouble( 0.0167 )  ) 
        )
    )



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

#scalesAndSmearingsPrefix = cms.string("EgammaAnalysis/ElectronTools/data/80X_ichepV2_2016_pho")
#scalesAndSmearingsPrefixForSigmaEOverE = cms.string("EgammaAnalysis/ElectronTools/data/Golden22June")
##scalesAndSmearingsPrefix = cms.string("EgammaAnalysis/ElectronTools/data/Winter_2016_reReco_v1_ele")
##scalesAndSmearingsPrefixForSigmaEOverE = cms.string("EgammaAnalysis/ElectronTools/data/Winter_2016_reReco_v1_ele")
scalesAndSmearingsPrefix = cms.string("EgammaAnalysis/ElectronTools/data/Moriond17_74x_pho")
scalesAndSmearingsPrefixForSigmaEOverE = cms.string("EgammaAnalysis/ElectronTools/data/Winter_2016_reReco_v1_ele")

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
          CorrectionFile = cms.FileInPath("flashgg/MicroAOD/data/transformationIDMVA_v2.root"),
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
          #new for low mass
          BinList3 = leadTriggerScaleBinsEBhiR9OR,
          BinList4 = subleadTriggerScaleBinsEBhiR9OR,
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
          CorrectionFile = scalesAndSmearingsPrefixForSigmaEOverE,
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

flashggDiPhotonSystematics = cms.EDProducer('FlashggDiPhotonSystematicProducer',
		src = cms.InputTag("flashggUpdatedIdMVADiPhotons"),
                SystMethods2D = cms.VPSet(),
                # the number of syst methods matches the number of nuisance parameters
                # assumed for a given systematic uncertainty and is NOT required
                # to match 1-to-1 the number of bins above.
                SystMethods = cms.VPSet(
                    MCScaleHighR9EB,
                    MCScaleLowR9EB,
                    MCScaleHighR9EE,
                    MCScaleLowR9EE,
                    MCScaleGain6EB_EGM,
                    MCScaleGain1EB_EGM,
                    MaterialCentralBarrel,
                    MaterialOuterBarrel,
                    MaterialForward,
                    ShowerShapeHighR9EB,
                    ShowerShapeHighR9EE,
                    ShowerShapeLowR9EB,
                    ShowerShapeLowR9EE,
                    FNUFEB,
                    FNUFEE,
                    MCSmearHighR9EE,
                    MCSmearLowR9EE,
                    MCSmearHighR9EB,
                    MCSmearLowR9EB,
                    MvaShift,
                    PreselSF,
                    electronVetoSF,
                    TriggerWeight,
                    LooseMvaSF,
                    SigmaEOverEShift,
                    SigmaEOverESmearing,
                    FracRVWeight,
                    FracRVNvtxWeight
                )
)
