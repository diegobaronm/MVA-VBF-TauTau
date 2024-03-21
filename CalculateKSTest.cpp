std::vector<std::string> methodNames = {"BDT-HM-10Folds"}; //,"BDT2Levels+Gini","BDT3Levels","BDTMoreTrees", "BDTEvenMoreTrees"};
int numberOfFolds = 10;

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

void CalculateKSTestFunction(const int nFolds,const std::vector<std::string>& methodNames){

    // Loop over all methods
    for (size_t m{0};m<methodNames.size();m++){
        std::cout << "SUMMARY" << methodNames[m] << std::endl;
        std::vector<double> signalKS{};
        std::vector<double> backgroundKS{};
        // Loop over all folds
        for(int fold{1};fold<=nFolds;fold++){
            // Open the file
            TFile* inputFile = TFile::Open((methodNames[m]+"_fold"+std::to_string(fold)+".root").c_str(),"READ");

            // Navigate to the correct directory
            std::string directoryName = "dataset/Method_BDT/"+methodNames[m]+"_fold"+std::to_string(fold)+"/";
            std::string histogramName = "MVA_"+methodNames[m]+"_fold"+std::to_string(fold);

            // Get Signal histograms
            TH1* sTest = (TH1*)inputFile->Get((directoryName+histogramName+"_S").c_str());
            TH1* sTrain = (TH1*)inputFile->Get((directoryName+histogramName+"_Train_S").c_str());
            // Get Background histograms
            TH1* bTest = (TH1*)inputFile->Get((directoryName+histogramName+"_B").c_str());
            TH1* bTrain = (TH1*)inputFile->Get((directoryName+histogramName+"_Train_B").c_str());

            // Normalise all histograms 
            // Signal test
            double integral = sTest->Integral(1,-1);
            sTest->Scale(1./integral);
            // Signal train
            integral = sTrain->Integral(1,-1);
            sTrain->Scale(1./integral);
            // Background test
            integral = bTest->Integral(1,-1);
            bTest->Scale(1./integral);
            // Background train
            integral = bTrain->Integral(1,-1);
            bTrain->Scale(1./integral);

            double kolS = sTrain->KolmogorovTest( sTest, "X" );
            double kolB = bTrain->KolmogorovTest( bTest, "X" );

            signalKS.push_back(kolS);
            backgroundKS.push_back(kolB);
            
            std::cout << "KS Fold "+std::to_string(fold)+" ------ " << "BG :" << kolB << " Signal: "<< kolS << std::endl;

            inputFile->Close();
        }
        // Signal average and std deviation
        double averageSignalKS{0};
        double stdDeviationSignalKS{0};
        // Background average and std deviation
        double averageBackgroundKS{0};
        double stdDeviationBackgroundKS{0};

        AverageAndStdDeviation(signalKS,averageSignalKS,stdDeviationSignalKS);
        AverageAndStdDeviation(backgroundKS,averageBackgroundKS,stdDeviationBackgroundKS);

        std::cout<< "Average KS Signal: " << averageSignalKS << " +/- " << stdDeviationSignalKS << std::endl;
        std::cout<< "Average KS Background: " << averageBackgroundKS << " +/- " << stdDeviationBackgroundKS << std::endl;

    }
}


void CalculateKSTest(){
    CalculateKSTestFunction(numberOfFolds, methodNames);
}