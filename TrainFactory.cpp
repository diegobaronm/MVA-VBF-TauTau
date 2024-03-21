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

void AverageAndStdDeviation(const std::vector<double>& data, double& average, double& stdDeviation) {
  average = 0;
  stdDeviation = 0;
  for (auto& value : data) {
    average += value;
  }
  average /= data.size();
  for (auto& value : data) {
    stdDeviation += (value - average) * (value - average);
  }
  stdDeviation /= (data.size()-1);
  stdDeviation = sqrt(stdDeviation);
}


void TrainFactory(UInt_t jobs = 1)
{
   TMVA::Tools::Instance();

   TFile* input{nullptr};
   TString fname = "/Users/user/Documents/HEP/VBF-Analysis/TauTau/TauhadTaulep/High-Mass/TrainingBDTSamples/TrainingData.root";
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
   dataloader->AddVariable("mjj", "m_{jj}", "GeV", 'D');
   dataloader->AddVariable("deltaRapidity", "#Delta y_{jj}", "", 'D');
   //dataloader->AddVariable("deltaPhiLT", "deltaPhiLT", "Units", 'D');
   //dataloader->AddVariable("jetRNNScore", "jetRNNScore", "Units", 'D');
   dataloader->AddVariable("ptBalance", "pT_{balance}", "", 'D');
   dataloader->AddVariable("zCentrality", "#xi(Z)", "", 'D');
   dataloader->AddVariable("omega", "#Omega", "", 'D');
   //dataloader->AddVariable("reco_mass", "reco_mass", "GeV", 'D');
   //dataloader->AddVariable("lepNuPt", "lepNuPt", "GeV", 'D');
   //dataloader->AddVariable("transverseMassLep", "m_{T}(l)", "GeV", 'D');
   dataloader->AddVariable("transverseRecoMassVariable", "m_{T}(l)/f(m_{reco})", "GeV", 'D');
   //dataloader->AddVariable( "transRecoMassRatio := transverseMassLep/pow(reco_mass,0.3)", "m_{T}(l)/(m_{reco}^{0.3})", "", 'D' );
   //dataloader->AddVariable("massTauLep", "massTauLep", "GeV", 'D');
   //dataloader->AddVariable("nLightJets", "nLightJets", "Units", 'I');
   //dataloader->AddVariable("tau_p4->Pt()", "tauPt", "GeV", 'D');
   //dataloader->AddVariable("lep_p4->Pt()", "lepPt", "GeV", 'D');
   //dataloader->AddVariable("jet0_p4->Pt()", "jet0Pt", "GeV", 'D');
   //dataloader->AddVariable("jet1_p4->Pt()", "jet1Pt", "GeV", 'D');
   //dataloader->AddVariable("met_p4->Pt()", "met", "GeV", 'D');
   dataloader->AddSpectator("eventNumber",'I');
   
   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   //dataloader->AddSpectator("spec1 := var1*2", "Spectator 1", "units", 'F');
   //dataloader->AddSpectator("spec2 := var1*3", "Spectator 2", "units", 'F');

   // global event weights per tree (see below for setting event-wise weights)
   //Double_t signalWeight = 1.0;
   //Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree(signalTree,1.0);
   dataloader->AddBackgroundTree(background,1.0);

   // Set individual event weights (the variables must exist in the original TTree)
   // -  for signal    : `dataloader->SetSignalWeightExpression    ("weight1*weight2");`
   // -  for background: `dataloader->SetBackgroundWeightExpression("weight1*weight2");`
   dataloader->SetSignalWeightExpression("mcWeight");
   dataloader->SetBackgroundWeightExpression("mcWeight");
   
   // specify the number of events used for training, the rest will be used for evaluation
   TCut preCut = "reco_mass>135 && mcWeight < 1"; // SplitSeed=70 // nTrain_Signal=11750:nTrain_Background=25690
   dataloader->PrepareTrainingAndTestTree(preCut,preCut, "nTrain_Signal=11750:nTrain_Background=25710:nTest_Signal=0:nTest_Background=0:SplitMode=Random:SplitSeed=70:NormMode=EqualNumEvents:V");
   int NumFolds = 10;
   TString splitExpr = ":SplitType=Deterministic:SplitExpr=int([eventNumber])\%int(10)";
   // open output file
   TFile *outputFile = TFile::Open("TrainResults.root", "RECREATE");
   (TMVA::gConfig().GetVariablePlotting()).fNbinsMVAoutput=10;

   TString opt = "!V:!Silent:!Correlations:AnalysisType=Classification:NumFolds=10:FoldFileOutput=True"+splitExpr;
   TMVA::CrossValidation cv {"validateBDT_Test", dataloader, outputFile, opt};

   // initialise the TMVA classification workflow
   //TMVA::Factory factory{"", outputFile, "AnalysisType=Classification"};
   //book classifiers to be trained with the hyperparameters setup
   //factory.BookMethod(dataloader,TMVA::Types::kBDT, "BDT1", "!H:!V:NTrees=125:MinNodeSize=3%:BoostType=Grad:Shrinkage=0.10:"
   //                                          "UseBaggedBoost=True:BaggedSampleFraction=0.5:UseNvars=2:nCuts=80:MaxDepth=3:SeparationType=SDivSqrtSPlusB");                                                                                  
   //factory.BookMethod(dataloader,TMVA::Types::kBDT, "BDT2", "!H:!V:NTrees=125:MinNodeSize=3%:BoostType=Grad:Shrinkage=0.10:"
   //                                          "UseBaggedBoost=True:BaggedSampleFraction=0.5:UseNvars=2:nCuts=80:MaxDepth=2");
   //factory.BookMethod(dataloader,TMVA::Types::kBDT, "BDTStandard", "!H:!V:NTrees=125:MinNodeSize=3%:BoostType=Grad:Shrinkage=0.10:"
   //                                         "UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=80:MaxDepth=3:SeparationType=CrossEntropy");                                                                                                                      
   // Train and Test all methods
   // Training and Evaluation
   //factory.TrainAllMethods();
   //factory.TestAllMethods();
   //factory.EvaluateAllMethods();

   // BDT with mT
   cv.BookMethod(TMVA::Types::kBDT, "BDT-HM-10Folds", "!H:!V:NTrees=125:MinNodeSize=3%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost=True:BaggedSampleFraction=0.5:UseNvars=2:nCuts=80:MaxDepth=2:SeparationType=SDivSqrtSPlusB");
   //BDT without mT
   //cv.BookMethod(TMVA::Types::kBDT, "BDT-HM-NOmT-10Folds", "!H:!V:NTrees=125:MinNodeSize=3%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost=True:BaggedSampleFraction=0.5:UseNvars=2:nCuts=80:MaxDepth=2:SeparationType=SDivSqrtSPlusB");
   //BDT without mT and cut
   //cv.BookMethod(TMVA::Types::kBDT, "BDT-HM-NOmTCut-10Folds", "!H:!V:NTrees=125:MinNodeSize=3%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost=True:BaggedSampleFraction=0.5:UseNvars=2:nCuts=80:MaxDepth=2:SeparationType=SDivSqrtSPlusB");
   cv.Evaluate();

   size_t iMethod = 0;
   for (auto && result : cv.GetResults()) {
      std::cout << "Summary for method "
               << cv.GetMethods()[iMethod++].GetValue<TString>("MethodTitle")
               << std::endl;
      std::vector<double> ROCValues{};
      for (UInt_t iFold = 0; iFold<cv.GetNumFolds(); ++iFold) {
         ROCValues.push_back(result.GetROCValues()[iFold]);
         std::cout << "\tFold " << iFold << ": "
         << "ROC int: " << result.GetROCValues()[iFold]
         << ", "
         << "SignalEff @ BkgEff=0.1: " << result.GetEff10Values()[iFold]
         << ", "
         << "Separation: " << result.GetSepValues()[iFold]
         << std::endl;
      } 
      double mean{};
      double stddev{};
      std::cout << "SUMMARY: " << std::endl;
      AverageAndStdDeviation(ROCValues, mean, stddev);
      std::cout << "ROC Integral: " << mean << " +- " << stddev << std::endl;
      AverageAndStdDeviation(result.GetSepValues(), mean, stddev);
      std::cout << "Separation: " << mean << " +- " << stddev << std::endl << std::endl;
   }
   // Clean up
   outputFile->Close();
}
