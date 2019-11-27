#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
//#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

/*#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#endif*/

void diphoPairMva_H4G()
{
   
    TString outfileName( "outputTMVA_diphoPair.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    TFile* inputFile = TFile::Open("/eos/cms/store/user/bmarzocc/H4G_Analysis/Dumpers_pairBDT/SUSYGluGluToHToAA_AToGG_M-60_TuneCUETP8M1_13TeV_pythia8_diPairMVA_final.root");

    TTree* signal_ggf  = (TTree*)inputFile->Get("diphotonPair_BDT_sig");
    TTree* bkg_ggf     = (TTree*)inputFile->Get("diphotonPair_BDT_bkg");
    
    TMVA::Tools::Instance();

    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    factory->AddVariable("dipho1_energy", "dipho1_energy", "", 'F'); 
    factory->AddVariable("dipho1_pt", "dipho1_pt", "", 'F');
    factory->AddVariable("dipho1_eta", "dipho1_eta", "", 'F'); 
    //factory->AddVariable("dipho1_phi", "dipho1_phi", "", 'F');
    factory->AddVariable("dipho1_dR", "dipho1_dR", "", 'F');
    factory->AddVariable("deltaM1_gen1", "deltaM1_gen1", "", 'F');
    factory->AddVariable("dipho2_energy", "dipho2_energy", "", 'F'); 
    factory->AddVariable("dipho2_pt", "dipho2_pt", "", 'F');
    factory->AddVariable("dipho2_eta", "dipho2_eta", "", 'F'); 
    //factory->AddVariable("dipho2_phi", "dipho2_phi", "", 'F');
    factory->AddVariable("dipho2_dR", "dipho2_dR", "", 'F');
    factory->AddVariable("deltaM2_gen1", "deltaM2_gen1", "", 'F');
    factory->AddVariable("dipair_dR", "dipair_dR", "", 'F');
   
    Double_t signalWeight = 1.0;
    Double_t backgroundWeight = 1.0;
  
    factory->AddSignalTree( signal_ggf, signalWeight );
    factory->AddBackgroundTree( bkg_ggf, backgroundWeight );

    //factory->SetWeightExpression( "genweight" );

    outputFile->cd();

    TCut mycuts = "fabs(dipho1_energy)<999. && fabs(dipho1_pt)<999. && fabs(dipho1_eta)<999. && fabs(dipho1_dR)<999. && fabs(deltaM1_gen1)<999. && fabs(dipho2_energy)<999. && fabs(dipho2_pt)<999. && fabs(dipho2_eta)<999. && fabs(dipho2_dR)<999. && fabs(deltaM2_gen1)<999. && fabs(dipair_dR)<999.";
    TCut mycutb = "fabs(dipho1_energy)<999. && fabs(dipho1_pt)<999. && fabs(dipho1_eta)<999. && fabs(dipho1_dR)<999. && fabs(deltaM1_gen1)<999. && fabs(dipho2_energy)<999. && fabs(dipho2_pt)<999. && fabs(dipho2_eta)<999. && fabs(dipho2_dR)<999. && fabs(deltaM2_gen1)<999. && fabs(dipair_dR)<999.";

    factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                    "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

    TMVA::MethodCategory* mcat = 0;

    TString theCat1Vars = "dipho1_energy:dipho1_pt:dipho1_eta:dipho1_dR:deltaM1_gen1:dipho2_energy:dipho2_pt:dipho2_eta:dipho2_dR:deltaM2_gen1:dipair_dR";
    TString theCat2Vars = "dipho1_energy:dipho1_pt:dipho1_eta:dipho1_dR:deltaM1_gen1:dipho2_energy:dipho2_pt:dipho2_eta:dipho2_dR:deltaM2_gen1:dipair_dR";
 
    TMVA::MethodBase* BDT_Cat = factory->BookMethod( TMVA::Types::kCategory, "diphoPairTMVA","" );
    mcat = dynamic_cast<TMVA::MethodCategory*>(BDT_Cat);

    mcat->AddMethod( "(deltaM1_gen1>-999. & deltaM1_gen1<=5.) && (deltaM2_gen1>-999. & deltaM2_gen1<=5.)", theCat1Vars, TMVA::Types::kBDT, "diphoPairTMVA_highRes","!H:!V:!CreateMVAPdfs:NTrees=1000:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:MaxDepth=3:MinNodeSize=2:NegWeightTreatment=ignorenegweightsintraining");
    mcat->AddMethod( "(deltaM1_gen1>5. && deltaM1_gen1<999.) || (deltaM2_gen1>5. && deltaM2_gen1<999.)",  theCat2Vars, TMVA::Types::kBDT, "diphoPairTMVA_lowRes","!H:!V:!CreateMVAPdfs:NTrees=1000:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:MaxDepth=3:MinNodeSize=2:NegWeightTreatment=ignorenegweightsintraining");

    factory->TrainAllMethods();

    factory->TestAllMethods();

    factory->EvaluateAllMethods();

   // --------------------------------------------------------------

    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;
}

