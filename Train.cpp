#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/Classification.h"

// The default output roc curve is (1-eff_bkg)-eff_sig, this function is used to convert the (1-eff_bkg) to rej_bkg (i.e. 1/eff_bkg) 
void ResetPoints(TGraph* graph) {
  Int_t nPoints = graph->GetN();
  for (Int_t i = 0; i < nPoints; ++i) {
    if(graph->GetPointY(i)!=1) graph->SetPointY(i, 1/(1-graph->GetPointY(i)));
    else graph->SetPointY(i, graph->GetPointY(i-1));
  }
}

void Train(UInt_t jobs = 1)
{
   TMVA::Tools::Instance();

   TFile* input{nullptr};
   TString fname = "/Users/diegomac/Documents/HEP/VBF-Analysis/TauTau/TauhadTaulep/High-Mass/MVATrainingData/TrainingData.root";
   if (!gSystem->AccessPathName(fname)) {
      input = TFile::Open(fname); // check if file in local directory exists
   } else {
      TFile::SetCacheFileDir(".");
      input = TFile::Open("TrainingData.root", "CACHEREAD");
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }

   // Register the training and test trees
   TTree *signalTree = (TTree *)input->Get("SIGNAL");
   TTree *background = (TTree *)input->Get("BACKGROUND");

   TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   dataloader->AddVariable("mjj", "mjj", "GeV", 'D');
   dataloader->AddVariable("deltaRapidity", "deltaRapidity", "Units", 'D');
   //dataloader->AddVariable("deltaPhiLT", "deltaPhiLT", "Units", 'D');
   //dataloader->AddVariable("jetRNNScore", "jetRNNScore", "Units", 'D');
   dataloader->AddVariable("ptBalance", "ptBalance", "Units", 'D');
   dataloader->AddVariable("zCentrality", "zCentrality", "Units", 'D');
   dataloader->AddVariable("omega", "omega", "Units", 'D');
   //dataloader->AddVariable("reco_mass", "reco_mass", "GeV", 'D');
   //dataloader->AddVariable("lepNuPt", "lepNuPt", "GeV", 'D');
   //dataloader->AddVariable("transverseMassLep", "transverseMassLep", "GeV", 'D');
   //dataloader->AddVariable("massTauLep", "massTauLep", "GeV", 'D');
   //dataloader->AddVariable("nLightJets", "nLightJets", "Units", 'I');
   //dataloader->AddVariable("tau_p4->Pt()", "tauPt", "GeV", 'D');
   //dataloader->AddVariable("lep_p4->Pt()", "lepPt", "GeV", 'D');
   //dataloader->AddVariable("jet0_p4->Pt()", "jet0Pt", "GeV", 'D');
   //dataloader->AddVariable("jet1_p4->Pt()", "jet1Pt", "GeV", 'D');
   //dataloader->AddVariable("met_p4->Pt()", "met", "GeV", 'D');
   
   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   //dataloader->AddSpectator("spec1 := var1*2", "Spectator 1", "units", 'F');
   //dataloader->AddSpectator("spec2 := var1*3", "Spectator 2", "units", 'F');

   // global event weights per tree (see below for setting event-wise weights)
   //Double_t signalWeight = 1.0;
   //Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree(signalTree);
   dataloader->AddBackgroundTree(background);

   // Set individual event weights (the variables must exist in the original TTree)
   // -  for signal    : `dataloader->SetSignalWeightExpression    ("weight1*weight2");`
   // -  for background: `dataloader->SetBackgroundWeightExpression("weight1*weight2");`
   dataloader->SetSignalWeightExpression("mcWeight");
   dataloader->SetBackgroundWeightExpression("mcWeight");
   
   // specify the number of events used for training, the rest will be used for evaluation
   TCut preCut = "reco_mass>135 && mcWeight < 1";
   dataloader->PrepareTrainingAndTestTree(preCut,preCut, "nTrain_Signal=11300:nTrain_Background=23000:nTest_Signal=0:nTest_Background=0:SplitMode=Random:SplitSeed=70:NormMode=EqualNumEvents:V");

   // open output file
   TFile *outputFile = TFile::Open("TrainResultsDefault.root", "RECREATE");
   (TMVA::gConfig().GetVariablePlotting()).fNbinsMVAoutput=10;

   // initialise the TMVA classification workflow
   TMVA::Factory factory{"", outputFile, "AnalysisType=Classification"};
   //book classifiers to be trained with the hyperparameters setup
   factory.BookMethod(dataloader,TMVA::Types::kBDT, "BDT-with-mT", "!H:!V:NTrees=125:MinNodeSize=3%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost=True:BaggedSampleFraction=0.5:UseNvars=2:nCuts=80:MaxDepth=2:SeparationType=SDivSqrtSPlusB");
   factory.BookMethod(dataloader,TMVA::Types::kBDT, "BDT-no-mT", "!H:!V:NTrees=200:MinNodeSize=3%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost=True:BaggedSampleFraction=0.5:UseNvars=2:nCuts=80:MaxDepth=2:SeparationType=SDivSqrtSPlusB");
                                                                                                                     
   // Train and Test all methods
   // Training and Evaluation
   factory.TrainAllMethods();
   factory.TestAllMethods();
   factory.EvaluateAllMethods();
   
   // Clean up
   outputFile->Close();
}
