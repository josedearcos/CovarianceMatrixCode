void LoganCrossCheck()
{
    
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);
    
    //Load spectra
    TFile* LoganSpectrumF = new TFile("roofile_rawE.root");
    LoganPromptH = (TH1D*)LoganSpectrumF->Get("hEp");//prompt spectrum
    LoganDelayH = (TH1D*)LoganSpectrumF->Get("hEd");//delay spectrum
    delete LoganSpectrumF;
    
    TFile* MySpectrumF = new TFile("../ResponseMatrices/Hydrogen/NominalResponseMatrix.root");
    
    TH2D* NominalResponseMatrixH = (TH2D*)MySpectrumF->Get("FineEvisEnu1");
    
    delete MySpectrumF;
    TH1D* TrueADH;
    
    const bool UseADAntineutrinoSpectrum;
    TFile* MyTrueADF;
    if(UseADAntineutrinoSpectrum){
        //this includes oscillation effect
        MyTrueADF = TFile::Open("../RootOutputs/Hydrogen/NominalOutputs/Oscillation.root");
        MyTrueADF->cd("Total AD Spectra after oscillation");
  
        TrueADH = (TH1D*)gDirectory->Get("Total spectrum after oscillation at AD1");
    }
    else
    {
        //reactor
        MyTrueADF = TFile::Open("../RootOutputs/Reactor/NominalOutputs/ReactorSpectrum.root");
        TrueADH = (TH1D*)gDirectory->Get("SpectrumFromReactor1");
    }
    MyTrueADF->Close();;
    
    TH1D* VisibleSpectrumH = new TH1D("VisibleH","VisibleH",NominalResponseMatrixH->GetYaxis()->GetNbins(),0,12);
    
    for(Int_t i=0; i<NominalResponseMatrixH->GetYaxis()->GetNbins(); i++)//visible
    {
        for(Int_t j=0; j<NominalResponseMatrixH->GetXaxis()->GetNbins(); j++)//true
        {
            //[True_0 ... True_n] = [0,0]  [...] [0,m]  * [Vis_0]
            //                      [0,n]  [...] [n,m]     [...]
            //                                            [Vis_m]
            
            // [nx1] = [n x m] x [mx1]
            
            std::cout << "Visible: " <<  VisibleSpectrumH->GetBinContent(i+1) << " - True: " << TrueADH->GetBinContent(j+1) << " - Visible bin: " << i <<  " - True bin: " << j << std::endl;

            VisibleSpectrumH->SetBinContent(i+1,VisibleSpectrumH->GetBinContent(i+1)+NominalResponseMatrixH->GetBinContent(j+1+(1.8/0.05),i+1)*TrueADH->GetBinContent(j+1));
            
        }
    }
    
    //Compare spectra
    
    ComparisonH = (TH1D*)LoganPromptH->Clone("Comparison");
    
    //normalize area
    VisibleSpectrumH->Scale(LoganPromptH->Integral()/VisibleSpectrumH->Integral());
    ComparisonH->Add(VisibleSpectrumH,-1);
    
    ComparisonH->Divide(LoganPromptH);
    
    TCanvas* ComparisonC = new TCanvas("PromptComparison","Prompt Comparison");
    ComparisonC->Divide(3,1);
    
    ComparisonC->cd(1);
    LoganPromptH->Draw();
    
    ComparisonC->cd(2);
    VisibleSpectrumH->Draw();
    
    ComparisonC->cd(3);
    ComparisonH->Draw();
    
    ComparisonC->Print("../Images/CrossChecks/LoganPrompt.eps");
    delete ComparisonC;
}