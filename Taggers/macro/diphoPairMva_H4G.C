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

    TFile* inputFile = TFile::Open("/eos/cms/store/user/bmarzocc/H4G_Analysis/Dumpers_pairBDT/SUSYGluGluToHToAA_AToGG_M-60_TuneCUETP8M1_13TeV_pythia8.root");

    TTree* signal_ggf  = (TTree*)inputFile->Get("FlashggH4GCandidate/diphotonPair_BDT_sig");
    TTree* bkg_ggf     = (TTree*)inputFile->Get("FlashggH4GCandidate/diphotonPair_BDT_bkg");
    
    TMVA::Tools::Instance();

    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    factory->AddVariable("dipho_energy", "dipho_energy", "", 'F'); 
    factory->AddVariable("dipho_pt", "dipho_pt", "", 'F');
    factory->AddVariable("dipho_eta", "dipho_eta", "", 'F'); 
    //factory->AddVariable("dipho_phi", "dipho_phi", "", 'F');
    factory->AddVariable("dipho_dR", "dipho_dR", "", 'F');
    factory->AddVariable("deltaM_gen1", "deltaM_gen1", "", 'F');
   
    Double_t signalWeight = 1.0;
    Double_t backgroundWeight = 1.0;
  
    factory->AddSignalTree( signal_ggf, signalWeight );
    factory->AddBackgroundTree( bkg_ggf, backgroundWeight );

    //factory->SetWeightExpression( "genweight" );

    outputFile->cd();

    TCut mycuts = "fabs(dipho_energy)<999. && fabs(dipho_pt)<999. && fabs(dipho_eta)<999. && fabs(dipho_dR)<999. && fabs(deltaM_gen1)<999.";
    TCut mycutb = "fabs(dipho_energy)<999. && fabs(dipho_pt)<999. && fabs(dipho_eta)<999. && fabs(dipho_dR)<999. && fabs(deltaM_gen1)<999.";

    factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                    "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

    TMVA::MethodCategory* mcat = 0;

    TString theCat1Vars = "dipho_energy:dipho_pt:dipho_eta:dipho_dR:deltaM_gen1";
    TString theCat2Vars = "dipho_energy:dipho_pt:dipho_eta:dipho_dR:deltaM_gen1";
 
    TMVA::MethodBase* BDT_Cat = factory->BookMethod( TMVA::Types::kCategory, "diphoPairTMVA","" );
    mcat = dynamic_cast<TMVA::MethodCategory*>(BDT_Cat);

    mcat->AddMethod( "deltaM_gen1>-999. & deltaM_gen1<=5.", theCat1Vars, TMVA::Types::kBDT, "diphoPairTMVA_highRes","!H:!V:!CreateMVAPdfs:NTrees=1000:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:MaxDepth=3:MinNodeSize=2:NegWeightTreatment=ignorenegweightsintraining");
    mcat->AddMethod( "deltaM_gen1>5. && deltaM_gen1<999.",  theCat2Vars, TMVA::Types::kBDT, "diphoPairTMVA_lowRes","!H:!V:!CreateMVAPdfs:NTrees=1000:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:MaxDepth=3:MinNodeSize=2:NegWeightTreatment=ignorenegweightsintraining");

    factory->TrainAllMethods();

    factory->TestAllMethods();

    factory->EvaluateAllMethods();

   // --------------------------------------------------------------

    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;
}

