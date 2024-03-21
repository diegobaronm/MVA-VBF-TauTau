std::map<std::string,std::vector<double>> binsMap {
{"mjj",{0.0,250.0,500.,750.0,1000.,1250.0,1500.0,2000.,3000.,4000.,5000.}},
{"deltaRapidity",{0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.}},
{"ptBalance",{0.0,0.01,0.02,0.03,0.05,0.1,0.15,0.2,0.5}},
{"zCentrality",{0.0,0.1,0.2,0.3,0.5,1.0,2.0}},
//{"jetRNNScore",{0.0,0.15,0.25,0.4,0.5,0.6,0.7,0.8,0.9,1.0}},
{"omega",{-2.0,-1.0,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,2.0}},
//{"reco_mass",{0,60,80,90,100,110,130,150,200,300,500,1000}},
{"transverseMassLep",{0,20,40,60,80,100,120,150,200,250,500}},
//{"tau_p4",{0,20,40,60,80,100,120,150,200}}*/
};

std::map<std::string,std::string> unitsMap {
{"mjj","GeV"},
{"deltaRapidity",""},
{"ptBalance",""},
{"zCentrality",""},
{"jetRNNScore",""},
{"omega",""},
{"reco_mass","GeV"},
{"transverseMassLep","GeV"},
{"tau_p4","GeV"}
};

std::map<std::string,std::vector<double>> initialBiningMap {
{"mjj",{5000,0,5000}},
{"deltaRapidity",{20,0,10}},
{"ptBalance",{50,0,0.5}},
{"zCentrality",{20,0,2}},
{"jetRNNScore",{20,0,1}},
{"omega",{40,-2,2}},
{"reco_mass",{100,0,1000}},
{"transverseMassLep",{50,0,500}},
{"tau_p4",{20,0,200}}
};

std::map<std::string,double> normMap {
{"mjj",250},
{"deltaRapidity",0.5},
{"ptBalance",0.03},
{"zCentrality",0.1},
{"jetRNNScore",0.15},
{"omega",0.2},
{"reco_mass",20},
{"transverseMassLep",20},
{"tau_p4",20}
};

std::map<std::string,std::string> xTitleMap {
{"mjj","m_{jj}"},
{"deltaRapidity","#Delta y_{jj}"},
{"ptBalance","pT_{balance}"},
{"zCentrality","#xi(Z)"},
{"jetRNNScore","jetRNNScore"},
{"omega","#Omega"},
{"reco_mass","reco_mass"},
{"transverseMassLep","m_{T}(l)"},
{"tau_p4","tauPt"}
};

bool isBackground = true;

void makePlot(){
    TH1::AddDirectory(kFALSE);
    // Set the pads
    TCanvas canvas("canvas","Variable comparisons",1000,600);
    canvas.cd();
    
    TPad* superPad = new TPad("superPad","Super Pad",0,0,1,1,kWhite,10);
    superPad->Divide(3,2);
    superPad->Draw();
    for (int i{1};i<=6;i++) {
        superPad->GetPad(i)->SetFillColor(kWhite);
        superPad->GetPad(i)->SetLogy();
    };
    
    std::string directoryPath ="/Users/diegomac/Documents/HEP/VBF-Analysis/TauTau/TauhadTaulep/High-Mass/MVATrainingData/";
    std::string directoryPath2 ="/Users/diegomac/Documents/HEP/VBF-Analysis/TauTau/TauhadTaulep/Z-peak/MVATrainingData/";
    std::string file1Name = "TrainingData.root";
    std::string file2Name = "Ztautau_1GapJet.root";
    std::string signalName = "TrainingData.root";
    //std::string signal2Name = "TrainingData.root";

    std::string NameOne = "QCDjj no-gap-jets";
    std::string NameTwo = "QCDjj 1-gap-jet";
    std::string NameSignal = "Signal";
    //std::string NameSignal2 = "Signal Z-peak";

    TFile* file1 = (TFile*)TFile::Open((directoryPath+file1Name).c_str(),"READ");
    TFile* file2 = (TFile*)TFile::Open((directoryPath+file2Name).c_str(),"READ");
    TFile* fileSignal = (TFile*)TFile::Open((directoryPath+signalName).c_str(),"READ");
    //TFile* fileSignal2 = (TFile*)TFile::Open((directoryPath2+signal2Name).c_str(),"READ");

    std::string treeName = isBackground ? "BACKGROUND" : "SIGNAL";

    TTree * tree1  = (TTree*)file1->Get(treeName.c_str());
    TTree * tree2  = (TTree*)file2->Get(treeName.c_str());
    TTree * signalTree  = (TTree*)fileSignal->Get("SIGNAL");
    //TTree * signalTree2  = (TTree*)fileSignal2->Get("SIGNAL");

    int i{1};
    std::string drawString;
    for (const auto& it : binsMap) {
        std::string histoName = it.first;
        std::cout<< histoName << std::endl;
        int nBins = initialBiningMap[histoName][0];
        double firstBin = initialBiningMap[histoName][1];
        double lastBin = initialBiningMap[histoName][2];
        double * bins = binsMap[histoName].data();
        drawString = histoName+">>h1("+std::to_string(nBins)+","+std::to_string(firstBin)+","+std::to_string(lastBin)+")";

        TLegend* legend = new TLegend(0.75,0.75,0.95,0.95);

        std::string cut = "reco_mass>135 && mcWeight < 1";
        std::string cutZ = "reco_mass>66 && mcWeight < 1";

        // Get Pad and Fill histogram and calculate integral
        superPad->GetPad(i)->cd();
        tree1->Draw(drawString.c_str(),cut.c_str(),"hist");
        TH1F* h1 = (TH1F*)superPad->GetPad(i)->GetPrimitive("h1");
        double integral = h1->Integral(1,-1);
        // Rebin and scale
        TH1F* h1new = (TH1F*)h1->Rebin(binsMap[histoName].size()-1,(histoName+"1").c_str(),bins);
        h1new->Scale(normMap[histoName]/integral,"width");
        h1new->SetLineColorAlpha(kBlue,0.6);
        // Remove old histogram and draw new one
        legend->AddEntry(h1new,NameOne.c_str(),"l");
        h1new->GetYaxis()->SetRangeUser(0.0001,4.5);
        
        // Second sample
        drawString = histoName+">>h2("+std::to_string(nBins)+","+std::to_string(firstBin)+","+std::to_string(lastBin)+")";
        // Get Pad and Fill histogram and calculate integral
        superPad->GetPad(i)->cd();
        tree2->Draw(drawString.c_str(),cut.c_str(),"hist same");
        TH1F* h2 = (TH1F*)superPad->GetPad(i)->GetPrimitive("h2");
        integral = h2->Integral(1,-1);
        // Rebin and scale
        TH1F* h2new = (TH1F*)h2->Rebin(binsMap[histoName].size()-1,(histoName+"2").c_str(),bins);
        h2new->Scale(normMap[histoName]/integral,"width");
        h2new->SetLineColor(kBlue+1);
        // Remove old histogram and draw new one
        legend->AddEntry(h2new,NameTwo.c_str(),"l");

        // Signal sample
        drawString = histoName+">>h3("+std::to_string(nBins)+","+std::to_string(firstBin)+","+std::to_string(lastBin)+")";
        // Get Pad and Fill histogram and calculate integral
        superPad->GetPad(i)->cd();
        signalTree->Draw(drawString.c_str(),cut.c_str(),"hist same");
        TH1F* h3 = (TH1F*)superPad->GetPad(i)->GetPrimitive("h3");
        integral = h3->Integral(1,-1);
        // Rebin and scale
        TH1F* h3new = (TH1F*)h3->Rebin(binsMap[histoName].size()-1,(histoName+"3").c_str(),bins);
        h3new->Scale(normMap[histoName]/integral,"width");
        h3new->SetLineColorAlpha(kRed,0.60);
        // Remove old histogram and draw new one
        legend->AddEntry(h3new,NameSignal.c_str(),"l");

        /*// Signal sample 2
        drawString = histoName+">>h4("+std::to_string(nBins)+","+std::to_string(firstBin)+","+std::to_string(lastBin)+")";
        // Get Pad and Fill histogram and calculate integral
        superPad->GetPad(i)->cd();
        signalTree2->Draw(drawString.c_str(),cutZ.c_str(),"hist same");
        TH1F* h4 = (TH1F*)superPad->GetPad(i)->GetPrimitive("h4");
        integral = h4->Integral(1,-1);
        // Rebin and scale
        TH1F* h4new = (TH1F*)h4->Rebin(binsMap[histoName].size()-1,(histoName+"4").c_str(),bins);
        h4new->Scale(normMap[histoName]/integral,"width");
        h4new->SetLineColor(kRed+1);
        // Remove old histogram and draw new one
        legend->AddEntry(h4new,NameSignal2.c_str(),"l");*/

        // Style
        gStyle->SetOptStat(0);
        std::string xTitle = xTitleMap[histoName]+" "+unitsMap[histoName];
        h1new->GetXaxis()->SetTitle(xTitle.c_str());
        std::string yTitle = "1/N Events/ "+std::to_string(std::round(100*normMap[histoName])/100)+" "+unitsMap[histoName];
        h1new->GetYaxis()->SetTitle(yTitle.c_str());
        h1new->SetTitle("");
        h1new->SetTitleSize(0.06,"X");
        h1new->SetTitleSize(0.05,"Y");
        h1new->SetTitleOffset(0.96,"X");
        h1new->GetYaxis()->SetLabelSize(0.056);
        h1new->GetXaxis()->SetLabelSize(0.056);
        if (histoName=="mjj" || histoName== "transverseMassLep") h1new->GetXaxis()->SetNdivisions(5,2,0,false);
        superPad->GetPad(i)->SetMargin(0.15,0.05,0.15,0.05);

        superPad->GetPad(i)->GetListOfPrimitives()->Remove(h1);
        superPad->GetPad(i)->GetListOfPrimitives()->Remove(h2);
        superPad->GetPad(i)->GetListOfPrimitives()->Remove(h3);
        //superPad->GetPad(i)->GetListOfPrimitives()->Remove(h4);
        h1new->Draw("hist e1 X0");
        h2new->Draw("hist e1 X0 same");
        h3new->Draw("hist e1 X0 same");
        //h4new->Draw("hist e1 X0 same");
        legend->Draw("same");

        i++;
    }


    file1->Close();
    file2->Close();
    fileSignal->Close();
    //fileSignal2->Close();

    canvas.Print("testing.pdf");
}