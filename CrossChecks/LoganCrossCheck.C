void LoganCrossCheck()
{
    
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);
    
    const Int_t volumesX = 2;//or 10
    const Int_t volumesY = 1;//or 10

    //Load spectra
    TFile* LoganSpectrumF = new TFile("roofile_rawE.root");
        LoganPromptH = (TH1D*)LoganSpectrumF->Get("hEp");//prompt spectrum
        LoganDelayH = (TH1D*)LoganSpectrumF->Get("hEd");//delay spectrum
    delete LoganSpectrumF;
    
    TH2D* NominalResponseMatrixH[volumesX][volumesY];//One per cell
    TFile* MySpectrumF = new TFile("../ResponseMatrices/Hydrogen/NominalResponseMatrix.root");
   
    for(Int_t x=0; x<volumesX; x++)
    {
        for(Int_t y=0; y<volumesY; y++)
        {
            NominalResponseMatrixH[x][y] = (TH2D*)MySpectrumF->Get(Form("FineEvisEnu1,Cell%d,%d",x,y));
        }
    }
    
    std::cout << "Nominal response matrixH dimensions: x " << NominalResponseMatrixH[0][0]->GetXaxis()->GetNbins() << " y " << NominalResponseMatrixH[0][0]->GetYaxis()->GetNbins() << std::endl;
    
    delete MySpectrumF;
    TH1D* TrueADH;
    
    const bool UseADAntineutrinoSpectrum = 0;
    const bool LoganInput = 1;
    
    const bool AverageEvents = 1;
    
    TFile* MyTrueADF;
    if(UseADAntineutrinoSpectrum){
        
        std::cout << " USING AD OSCILLATED SPECTRUM" << std::endl;

        //this includes oscillation effect
        MyTrueADF = TFile::Open("../RootOutputs/Hydrogen/NominalOutputs/Oscillation.root");
        MyTrueADF->cd("Total AD Spectra after oscillation");
  
        TrueADH = (TH1D*)gDirectory->Get("Total spectrum after oscillation at AD1");
    }
    else if(LoganInput==0)
    {
        std::cout << " USING NOMINAL REACTOR SPECTRUM" << std::endl;

        //reactor
        MyTrueADF = TFile::Open("../RootOutputs/Hydrogen/NominalOutputs/AntineutrinoSpectrum.root");
        TrueADH = (TH1D*)gDirectory->Get("AntineutrinoSpectrumFromReactor1");
        
    }
    else
    {
        std::cout << " USING LOGAN'S H_EV " << std::endl;
        
        MyTrueADF = TFile::Open("./file_h_Ev.root");
        TrueADH = (TH1D*)gDirectory->Get("h_Ev");
    }
    
    MyTrueADF->Close();;
    
    std::cout << "True Energy Spectrum dimension: x " << TrueADH->GetXaxis()->GetNbins() << std::endl;
    
    TH1D* VisibleSpectrumH[volumesX][volumesY];
    
    TH2D* AverageH;

    if(AverageEvents)
    {
        //Read normalization map
        H2MapF = TFile::Open("../Inputs/HInputs/MapEventRatio2Total.root");
        
        AverageH = (TH2D*)gDirectory->Get("MapEventsRatio2total");
        
        delete H2MapF;
        
        std::cout << "Integral should be 1 : " << AverageH->Integral() << std::endl;
        
    }
    else
    {
        AverageH = new TH2D("MapEventsRatio2center","MapEventsRatio2center",volumesX,0,4,volumesY,-2,2);
        
        for(Int_t x=0; x<volumesX; x++)
        {
            for(Int_t y=0; y<volumesY; y++)
            {
                AverageH->SetBinContent(x+1,y+1,1);
            }
        }
    }
    
    for(Int_t x=0; x<volumesX; x++)
    {
        for(Int_t y=0; y<volumesY; y++)
        {
            VisibleSpectrumH[x][y] = new TH1D(Form("VisibleHCell%d,%d",x,y),Form("VisibleHCell%d,%d",x,y),NominalResponseMatrixH[0][0]->GetYaxis()->GetNbins(),0,12);

            for(Int_t i=0; i<NominalResponseMatrixH[0][0]->GetYaxis()->GetNbins(); i++)//visible
            {
                for(Int_t j=0; j<NominalResponseMatrixH[0][0]->GetXaxis()->GetNbins()-(1.8/0.05); j++)//true
                {
                    //[True_0 ... True_n] = [0,0]  [...] [0,m]  * [Vis_0]
                    //                      [0,n]  [...] [n,m]     [...]
                    //                                            [Vis_m]
                    
                    // [nx1] = [n x m] x [mx1]
                    
                    //std::cout << "Visible: " <<  VisibleSpectrumH[x][y]->GetBinContent(i+1) << " - True: " << TrueADH->GetBinContent(j+1) << " - Visible bin: " << i <<  " - True bin: " << j << std::endl;
                    
                    VisibleSpectrumH[x][y]->SetBinContent(i+1,VisibleSpectrumH[x][y]->GetBinContent(i+1)+NominalResponseMatrixH[x][y]->GetBinContent(j+1+(1.8/0.05),i+1)*TrueADH->Interpolate(NominalResponseMatrixH[x][y]->GetXaxis()->GetBinCenter(j+1+(1.8/0.05))));
                   
                    if(x==0&&y==0&&i==0)
                    {
                        std::cout << "energy : " << NominalResponseMatrixH[x][y]->GetXaxis()->GetBinCenter(j+1+(1.8/0.05)) << " in bin " << j+1 << std::endl;
                    }
                   // std::cout << " Energy bin j " << j << " - " << NominalResponseMatrixH[x][y]->GetXaxis()->GetBinCenter(j+1+(1.8/0.05)) << std::endl;
                
                }
            }
        }
    }
    
    TH1D* TotalVisibleSpectrumH = new TH1D("VisibleH","VisibleH",NominalResponseMatrixH[0][0]->GetYaxis()->GetNbins(),0,12);

    for(Int_t x=0; x<volumesX; x++)
    {
        for(Int_t y=0; y<volumesY; y++)
        {
            //Add all cells:
            TotalVisibleSpectrumH->Add(VisibleSpectrumH[x][y],AverageH->GetBinContent(x+1,y+1));
        }
    }
    //Compare spectra
    
    ComparisonH = (TH1D*)LoganPromptH->Clone("Comparison");
    
    std::cout << " My file integral: " << TotalVisibleSpectrumH->Integral() << " vs Logan's Integral: " << LoganPromptH->Integral() << std::endl;
                                           
    //normalize area
    TotalVisibleSpectrumH->Scale(LoganPromptH->Integral()/TotalVisibleSpectrumH->Integral());
    
    ComparisonH->Add(TotalVisibleSpectrumH,-1);
    
    ComparisonH->Divide(LoganPromptH);
    
    TCanvas* ComparisonC = new TCanvas("PromptComparison","Prompt Comparison");
    ComparisonC->Divide(3,1);
    
    ComparisonC->cd(1);
    LoganPromptH->Draw();
    
    ComparisonC->cd(2);
    TotalVisibleSpectrumH->Draw();
    
    ComparisonC->cd(3);
    ComparisonH->Draw();
    
    if(UseADAntineutrinoSpectrum){
        
        ComparisonC->Print("../Images/CrossChecks/LoganPromptUsingAD1.eps");

    }
    else if(LoganInput==0)
    {
        ComparisonC->Print("../Images/CrossChecks/LoganPromptUsingReactor1.eps");

    }
    else
    {
        ComparisonC->Print("../Images/CrossChecks/LoganPromptUsingLoganEv.eps");
    }
    
    delete ComparisonC;
}