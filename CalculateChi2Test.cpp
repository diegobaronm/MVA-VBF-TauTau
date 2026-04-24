#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"

// This script calculates the KS test for the signal and background distributions for multiple MVAs and folds.

std::vector<std::string> methodNames = {"BDT_NewCuts_NewQCD_NoRegSplits"}; //,"BDT2Levels+Gini","BDT3Levels","BDTMoreTrees", "BDTEvenMoreTrees"};
int numberOfFolds = 5;

// Calculate the mean and standard deviation of a vector of doubles.
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

double ComputeHistogramChi2(const TH1* referenceHistogram, const TH1* comparisonHistogram, bool ignoreReferenceVariance = false) {
  if (referenceHistogram == nullptr || comparisonHistogram == nullptr) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const int numberOfBins = referenceHistogram->GetNbinsX();
  if (numberOfBins != comparisonHistogram->GetNbinsX()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double chi2 = 0.0;
  for (int bin = 1; bin <= numberOfBins; ++bin) {
    const double referenceCount = referenceHistogram->GetBinContent(bin);
    const double comparisonCount = comparisonHistogram->GetBinContent(bin);
    const double referenceError = referenceHistogram->GetBinError(bin);
    const double comparisonError = comparisonHistogram->GetBinError(bin);
    double variance;
    if (ignoreReferenceVariance) {
      variance = comparisonError * comparisonError;
    } else {
      variance = referenceError * referenceError + comparisonError * comparisonError;
    }

    if (variance <= 0.0) {
      continue;
    }

    const double difference = referenceCount - comparisonCount;
    chi2 += difference * difference / variance;
  }

  return chi2;
}

TH1* CloneHistogramDetached(const TH1* histogram, const std::string& cloneName) {
  if (histogram == nullptr) {
    return nullptr;
  }

  TH1* clone = (TH1*)histogram->Clone(cloneName.c_str());
  if (clone != nullptr) {
    clone->SetDirectory(nullptr);
  }
  return clone;
}

void AddHistogramToSum(TH1*& accumulatedHistogram, const TH1* histogram, const std::string& histogramName) {
  if (histogram == nullptr) {
    return;
  }

  if (accumulatedHistogram == nullptr) {
    accumulatedHistogram = CloneHistogramDetached(histogram, histogramName);
    return;
  }

  if (accumulatedHistogram->GetNbinsX() != histogram->GetNbinsX()) {
    std::cerr << "Cannot sum histograms with different bin counts: "
              << accumulatedHistogram->GetName() << " and " << histogram->GetName() << std::endl;
    return;
  }

  for (int bin = 0; bin <= accumulatedHistogram->GetNbinsX() + 1; ++bin) {
    const double accumulatedContent = accumulatedHistogram->GetBinContent(bin);
    const double accumulatedError = accumulatedHistogram->GetBinError(bin);
    const double histogramContent = histogram->GetBinContent(bin);
    const double histogramError = histogram->GetBinError(bin);

    accumulatedHistogram->SetBinContent(bin, accumulatedContent + histogramContent);
    accumulatedHistogram->SetBinError(
      bin,
      std::sqrt(accumulatedError * accumulatedError + histogramError * histogramError)
    );
  }
}

void SaveFoldComparisonPlot(
  const std::string& methodName,
  const int fold,
  TH1* signalTrain,
  TH1* signalTest,
  TH1* backgroundTrain,
  TH1* backgroundTest,
  const double chi2Signal,
  const double chi2Background
) {
  if (signalTrain == nullptr || signalTest == nullptr || backgroundTrain == nullptr || backgroundTest == nullptr) {
    return;
  }

  gStyle->SetOptStat(0);

  TCanvas canvas((methodName + "_fold" + std::to_string(fold) + "_canvas").c_str(), "Train/Test comparison", 900, 700);
  canvas.cd();
  canvas.SetTicks(1, 1);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.05);
  canvas.SetTopMargin(0.10);
  canvas.SetBottomMargin(0.14);

  TH1* signalTrainBand = (TH1*)signalTrain->Clone((std::string(signalTrain->GetName()) + "_band").c_str());
  TH1* signalTrainLine = (TH1*)signalTrain->Clone((std::string(signalTrain->GetName()) + "_line").c_str());
  TH1* signalTestPoints = (TH1*)signalTest->Clone((std::string(signalTest->GetName()) + "_points").c_str());
  TH1* backgroundTrainBand = (TH1*)backgroundTrain->Clone((std::string(backgroundTrain->GetName()) + "_band").c_str());
  TH1* backgroundTrainLine = (TH1*)backgroundTrain->Clone((std::string(backgroundTrain->GetName()) + "_line").c_str());
  TH1* backgroundTestPoints = (TH1*)backgroundTest->Clone((std::string(backgroundTest->GetName()) + "_points").c_str());

  signalTrainBand->SetDirectory(nullptr);
  signalTrainLine->SetDirectory(nullptr);
  signalTestPoints->SetDirectory(nullptr);
  backgroundTrainBand->SetDirectory(nullptr);
  backgroundTrainLine->SetDirectory(nullptr);
  backgroundTestPoints->SetDirectory(nullptr);

  signalTrainBand->SetTitle(((std::string("Signal and Background - Fold ") + std::to_string(fold)) + ";BDT score;Events").c_str());
  signalTrainBand->SetLineColor(kRed + 1);
  signalTrainBand->SetFillColorAlpha(kRed + 1, 0.35);
  signalTrainBand->SetMarkerSize(0.0);
  signalTrainLine->SetLineColor(kRed + 1);
  signalTrainLine->SetLineWidth(2);
  signalTrainLine->SetFillStyle(0);
  signalTestPoints->SetLineColor(kBlack);
  signalTestPoints->SetMarkerColor(kRed);
  signalTestPoints->SetMarkerStyle(20);
  signalTestPoints->SetMarkerSize(0.9);

  backgroundTrainBand->SetLineColor(kBlue + 1);
  backgroundTrainBand->SetFillColorAlpha(kBlue + 1, 0.35);
  backgroundTrainBand->SetMarkerSize(0.0);
  backgroundTrainLine->SetLineColor(kBlue + 1);
  backgroundTrainLine->SetLineWidth(2);
  backgroundTrainLine->SetFillStyle(0);
  backgroundTestPoints->SetLineColor(kBlack);
  backgroundTestPoints->SetMarkerColor(kBlue);
  backgroundTestPoints->SetMarkerStyle(21);
  backgroundTestPoints->SetMarkerSize(0.9);

  double maximum = signalTrainBand->GetMaximum();
  maximum = std::max(maximum, signalTestPoints->GetMaximum());
  maximum = std::max(maximum, backgroundTrainBand->GetMaximum());
  maximum = std::max(maximum, backgroundTestPoints->GetMaximum());

  signalTrainBand->SetMaximum(maximum * 1.35);
  signalTrainBand->SetMinimum(0.0);
  signalTrainBand->GetXaxis()->SetTitleSize(0.060);
  signalTrainBand->GetYaxis()->SetTitleSize(0.060);
  signalTrainBand->GetXaxis()->SetLabelSize(0.050);
  signalTrainBand->GetYaxis()->SetLabelSize(0.050);
  signalTrainBand->GetXaxis()->SetTickLength(0.035);
  signalTrainBand->GetYaxis()->SetTickLength(0.020);
  signalTrainBand->GetXaxis()->SetTitleOffset(1.05);
  signalTrainBand->GetYaxis()->SetTitleOffset(0.9);

  signalTrainBand->Draw("E2");
  signalTrainLine->Draw("HIST SAME");
  backgroundTrainBand->Draw("E2 SAME");
  backgroundTrainLine->Draw("HIST SAME");
  signalTestPoints->Draw("E1 SAME");
  backgroundTestPoints->Draw("E1 SAME");

  TLegend legend(0.54, 0.64, 0.90, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(signalTrainBand, "Signal train", "lf");
  legend.AddEntry(signalTestPoints, "Signal test", "lep");
  legend.AddEntry(backgroundTrainBand, "Background train", "lf");
  legend.AddEntry(backgroundTestPoints, "Background test", "lep");
  legend.Draw();

  std::ostringstream signalChi2Label;
  signalChi2Label.setf(std::ios::fixed);
  signalChi2Label.precision(3);
  signalChi2Label << "Signal #chi^{2}/dof = " << chi2Signal;

  std::ostringstream backgroundChi2Label;
  backgroundChi2Label.setf(std::ios::fixed);
  backgroundChi2Label.precision(3);
  backgroundChi2Label << "Background #chi^{2}/dof = " << chi2Background;

  TLatex latex;
  latex.SetNDC(true);
  latex.SetTextSize(0.040);
  latex.DrawLatex(0.14, 0.84, signalChi2Label.str().c_str());
  latex.DrawLatex(0.14, 0.78, backgroundChi2Label.str().c_str());

  canvas.SaveAs((methodName + "_fold" + std::to_string(fold) + "_train_test_comparison.pdf").c_str());

  delete signalTrainBand;
  delete signalTrainLine;
  delete signalTestPoints;
  delete backgroundTrainBand;
  delete backgroundTrainLine;
  delete backgroundTestPoints;
}

void SaveCombinedComparisonPlot(
  const std::string& methodName,
  TH1* signalTrain,
  TH1* signalTest,
  TH1* backgroundTrain,
  TH1* backgroundTest,
  const double chi2Signal,
  const double chi2Background
) {
  if (signalTrain == nullptr || signalTest == nullptr || backgroundTrain == nullptr || backgroundTest == nullptr) {
    return;
  }

  gStyle->SetOptStat(0);

  TCanvas canvas((methodName + "_combined_canvas").c_str(), "Combined train/test comparison", 900, 700);
  canvas.cd();
  canvas.SetTicks(1, 1);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.05);
  canvas.SetTopMargin(0.10);
  canvas.SetBottomMargin(0.14);

  //TH1* signalTrainBand = CloneHistogramDetached(signalTrain, methodName + "_combined_signal_train_band");
  TH1* signalTrainLine = CloneHistogramDetached(signalTrain, methodName + "_combined_signal_train_line");
  TH1* signalTestPoints = CloneHistogramDetached(signalTest, methodName + "_combined_signal_test_points");
  //TH1* backgroundTrainBand = CloneHistogramDetached(backgroundTrain, methodName + "_combined_background_train_band");
  TH1* backgroundTrainLine = CloneHistogramDetached(backgroundTrain, methodName + "_combined_background_train_line");
  TH1* backgroundTestPoints = CloneHistogramDetached(backgroundTest, methodName + "_combined_background_test_points");

  //signalTrainBand->SetTitle((std::string("Combined signal and background;BDT score;Events") ).c_str());
  //signalTrainBand->SetLineColor(kRed + 1);
  //signalTrainBand->SetFillColorAlpha(kRed + 1, 0.35);
  //signalTrainBand->SetMarkerSize(0.0);
  signalTrainLine->SetTitle((std::string("Combined signal and background;BDT score;Events") ).c_str());
  signalTrainLine->SetLineColor(kRed + 1);
  signalTrainLine->SetLineWidth(2);
  signalTrainLine->SetFillStyle(0);
  signalTestPoints->SetLineColor(kBlack);
  signalTestPoints->SetMarkerColor(kRed);
  signalTestPoints->SetMarkerStyle(20);
  signalTestPoints->SetMarkerSize(0.9);

  //backgroundTrainBand->SetLineColor(kBlue + 1);
  //backgroundTrainBand->SetFillColorAlpha(kBlue + 1, 0.35);
  //backgroundTrainBand->SetMarkerSize(0.0);
  backgroundTrainLine->SetLineColor(kBlue + 1);
  backgroundTrainLine->SetLineWidth(2);
  backgroundTrainLine->SetFillStyle(0);
  backgroundTestPoints->SetLineColor(kBlack);
  backgroundTestPoints->SetMarkerColor(kBlue);
  backgroundTestPoints->SetMarkerStyle(21);
  backgroundTestPoints->SetMarkerSize(0.9);

  double maximum = backgroundTestPoints->GetMaximum();
  maximum = std::max(maximum, signalTestPoints->GetMaximum());
  //maximum = std::max(maximum, backgroundTrainBand->GetMaximum());
  maximum = std::max(maximum, backgroundTestPoints->GetMaximum());

  signalTrainLine->SetMaximum(maximum * 1.35);
  signalTrainLine->SetMinimum(0.0);
  signalTrainLine->GetXaxis()->SetTitleSize(0.060);
  signalTrainLine->GetYaxis()->SetTitleSize(0.060);
  signalTrainLine->GetXaxis()->SetLabelSize(0.050);
  signalTrainLine->GetYaxis()->SetLabelSize(0.050);
  signalTrainLine->GetXaxis()->SetTickLength(0.035);
  signalTrainLine->GetYaxis()->SetTickLength(0.020);
  signalTrainLine->GetXaxis()->SetTitleOffset(1.05);
  signalTrainLine->GetYaxis()->SetTitleOffset(0.9);

  //signalTrainBand->Draw("E2");
  signalTrainLine->Draw("HIST");
  //backgroundTrainBand->Draw("E2 SAME");
  backgroundTrainLine->Draw("HIST SAME");
  signalTestPoints->Draw("E1 SAME");
  backgroundTestPoints->Draw("E1 SAME");

  TLegend legend(0.54, 0.64, 0.90, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  //legend.AddEntry(signalTrainBand, "Signal train sum", "lf");
  legend.AddEntry(signalTestPoints, "Signal test sum", "lep");
  //legend.AddEntry(backgroundTrainBand, "Background train sum", "lf");
  legend.AddEntry(backgroundTestPoints, "Background test sum", "lep");
  legend.Draw();

  std::ostringstream signalChi2Label;
  signalChi2Label.setf(std::ios::fixed);
  signalChi2Label.precision(3);
  signalChi2Label << "Signal #chi^{2}/dof = " << chi2Signal;

  std::ostringstream backgroundChi2Label;
  backgroundChi2Label.setf(std::ios::fixed);
  backgroundChi2Label.precision(3);
  backgroundChi2Label << "Background #chi^{2}/dof = " << chi2Background;

  TLatex latex;
  latex.SetNDC(true);
  latex.SetTextSize(0.040);
  latex.DrawLatex(0.14, 0.84, signalChi2Label.str().c_str());
  latex.DrawLatex(0.14, 0.78, backgroundChi2Label.str().c_str());

  canvas.SaveAs((methodName + "_combined_train_test_comparison.pdf").c_str());

  //delete signalTrainBand;
  delete signalTrainLine;
  delete signalTestPoints;
  //delete backgroundTrainBand;
  delete backgroundTrainLine;
  delete backgroundTestPoints;
}

// Main function to calculate the KS test results.
void CalculateKSTestFunction(const int nFolds,const std::vector<std::string>& methodNames){

    // Loop over all methods
    for (size_t m{0};m<methodNames.size();m++){
        std::cout << "SUMMARY" << methodNames[m] << std::endl;
        std::vector<double> signalKS{};
        std::vector<double> backgroundKS{};
        std::vector<double> chi2Signal{};
        std::vector<double> chi2Background{};
        TH1* summedSignalTrain{nullptr};
        TH1* summedSignalTest{nullptr};
        TH1* summedBackgroundTrain{nullptr};
        TH1* summedBackgroundTest{nullptr};
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

            AddHistogramToSum(summedSignalTrain, sTrain, methodNames[m] + "_summed_signal_train");
            AddHistogramToSum(summedSignalTest, sTest, methodNames[m] + "_summed_signal_test");
            AddHistogramToSum(summedBackgroundTrain, bTrain, methodNames[m] + "_summed_background_train");
            AddHistogramToSum(summedBackgroundTest, bTest, methodNames[m] + "_summed_background_test");

            double kolS = sTrain->KolmogorovTest( sTest, "X=10000" );
            double kolB = bTrain->KolmogorovTest( bTest, "X=10000" );

            double chi2S = ComputeHistogramChi2(sTrain, sTest)/sTrain->GetNbinsX();
            double chi2B = ComputeHistogramChi2(bTrain, bTest)/bTrain->GetNbinsX();
            chi2Signal.push_back(chi2S);
            chi2Background.push_back(chi2B);
            std::cout << "Chi2 Fold "+std::to_string(fold)+" ------ " << "BG :" << chi2B << " Signal: "<< chi2S << std::endl;

            SaveFoldComparisonPlot(methodNames[m], fold, sTrain, sTest, bTrain, bTest, chi2S, chi2B);

            signalKS.push_back(kolS);
            backgroundKS.push_back(kolB);
            
            //std::cout << "KS Fold "+std::to_string(fold)+" ------ " << "BG :" << kolB << " Signal: "<< kolS << std::endl;

            inputFile->Close();
        }
        // Signal average and std deviation
        double averageSignalKS{0};
        double stdDeviationSignalKS{0};
        // Background average and std deviation
        double averageBackgroundKS{0};
        double stdDeviationBackgroundKS{0};

        AverageAndStdDeviation(signalKS, averageSignalKS, stdDeviationSignalKS);
        AverageAndStdDeviation(backgroundKS, averageBackgroundKS, stdDeviationBackgroundKS);

        double averageChi2Signal{0};
        double stdDeviationChi2Signal{0};
        double averageChi2Background{0};
        double stdDeviationChi2Background{0};
        AverageAndStdDeviation(chi2Signal, averageChi2Signal, stdDeviationChi2Signal);
        AverageAndStdDeviation(chi2Background, averageChi2Background, stdDeviationChi2Background);

        double combinedChi2Signal = std::numeric_limits<double>::quiet_NaN();
        double combinedChi2Background = std::numeric_limits<double>::quiet_NaN();
        if (summedSignalTrain != nullptr && summedSignalTest != nullptr) {
          combinedChi2Signal = ComputeHistogramChi2(summedSignalTrain, summedSignalTest, true) / summedSignalTrain->GetNbinsX();
        }
        if (summedBackgroundTrain != nullptr && summedBackgroundTest != nullptr) {
          combinedChi2Background = ComputeHistogramChi2(summedBackgroundTrain, summedBackgroundTest, true) / summedBackgroundTrain->GetNbinsX();
        }
        SaveCombinedComparisonPlot(
          methodNames[m],
          summedSignalTrain,
          summedSignalTest,
          summedBackgroundTrain,
          summedBackgroundTest,
          combinedChi2Signal,
          combinedChi2Background
        );

        std::cout<< "Average Chi2 Signal: " << averageChi2Signal << " +/- " << stdDeviationChi2Signal << std::endl;
        std::cout<< "Average Chi2 Background: " << averageChi2Background << " +/- " << stdDeviationChi2Background << std::endl;
        std::cout<< "Combined Chi2 Signal: " << combinedChi2Signal << std::endl;
        std::cout<< "Combined Chi2 Background: " << combinedChi2Background << std::endl;

        std::cout<< "Average KS Signal: " << averageSignalKS << " +/- " << stdDeviationSignalKS << std::endl;
        std::cout<< "Average KS Background: " << averageBackgroundKS << " +/- " << stdDeviationBackgroundKS << std::endl;

        delete summedSignalTrain;
        delete summedSignalTest;
        delete summedBackgroundTrain;
        delete summedBackgroundTest;

    }
}


void CalculateKSTest(){
    CalculateKSTestFunction(numberOfFolds, methodNames);
}