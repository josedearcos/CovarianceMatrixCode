void LoganCrossCheck()
{
    
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);
    
    const Int_t NADs = 6;
    const Int_t volumesX = 2;//or 10
    const Int_t volumesY = 1;//or 10
    const Int_t TrueBins = 39;
    const Int_t VisBins = 34;
    //Load spectra
    TFile* LoganSpectrumF = new TFile("roofile_rawE.root");
        LoganPromptH = (TH1D*)LoganSpectrumF->Get("hEp");//prompt spectrum
        LoganDelayH = (TH1D*)LoganSpectrumF->Get("hEd");//delay spectrum
    delete LoganSpectrumF;
    
    TH2D* NominalResponseMatrixH[volumesX][volumesY];//One per cell
    TFile* MySpectrumF = new TFile(Form("../ResponseMatrices/Hydrogen/NominalResponseMatrix%i_%i.root",VisBins,TrueBins));
   
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
    
    TH1D* VisibleSpectrumH[15*volumesX*volumesY][NADs][volumesX][volumesY];
    
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
            VisibleSpectrumH[0][0][x][y] = new TH1D(Form("VisibleHCell%d,%d",x,y),Form("VisibleHCell%d,%d",x,y),NominalResponseMatrixH[0][0]->GetYaxis()->GetNbins(),0,12);

            for(Int_t i=0; i<NominalResponseMatrixH[0][0]->GetYaxis()->GetNbins(); i++)//visible
            {
                for(Int_t j=0; j<NominalResponseMatrixH[0][0]->GetXaxis()->GetNbins()-(1.8/0.05); j++)//true
                {
                    //[True_0 ... True_n] = [0,0]  [...] [0,m]  * [Vis_0]
                    //                      [0,n]  [...] [n,m]     [...]
                    //                                            [Vis_m]
                    
                    // [nx1] = [n x m] x [mx1]
                    
                    //std::cout << "Visible: " <<  VisibleSpectrumH[0][x][y]->GetBinContent(i+1) << " - True: " << TrueADH->GetBinContent(j+1) << " - Visible bin: " << i <<  " - True bin: " << j << std::endl;
                    
                    VisibleSpectrumH[0][0][x][y]->SetBinContent(i+1,VisibleSpectrumH[0][0][x][y]->GetBinContent(i+1)+NominalResponseMatrixH[x][y]->GetBinContent(j+1+(1.8/0.05),i+1)*TrueADH->Interpolate(NominalResponseMatrixH[x][y]->GetXaxis()->GetBinCenter(j+1+(1.8/0.05))));
                   
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
            TotalVisibleSpectrumH->Add(VisibleSpectrumH[0][0][x][y],AverageH->GetBinContent(x+1,y+1));
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
    
    
    //Show all the outputs from each respective systematic:
    
    TH2D* Histo[15][NADs][volumesX][volumesY];
    
    Int_t num = 1;
    TString fname[100];
    
    const char *dirname;

    dirname= "../ResponseMatrices/Hydrogen/";
    char *ext= (Form("34_39.root"));
    
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    
    if(files)
    {
        TSystemFile *file;
        TIter next(files);
        
        while ((file=(TSystemFile*)next()))
        {
            fname[num] = file->GetName();
            
            if (!file->IsDirectory() && fname[num].EndsWith(ext))
            {
                std::cout << fname[num] << std::endl;
                
                TFile* FileF = TFile::Open(dirname+fname[num]);
                for(Int_t AD=0; AD<NADs; AD++)
                {
                    for(Int_t x=0; x<volumesX; x++)
                    {
                        for(Int_t y=0; y<volumesY; y++)
                        {
                            Histo[num][AD][x][y] = (TH2D*)FileF->Get(Form("FineEvisEnu%i,Cell%d,%d",AD+1,x,y));
                        }
                    }
                }
                FileF->Close();
                num++;
            }
        }
    }
    
    
    TCanvas* CanvasSys[15];
    
    for (Int_t systematic = 1; systematic<num; systematic++)
    {
        CanvasSys[systematic] = new TCanvas(Form("CanvasSys%d",systematic),Form("CanvasSys%d",systematic));
        
        CanvasSys[systematic]->Divide(NADs,volumesX);

        for(Int_t AD=0; AD<NADs; AD++)
        {
            for(Int_t x=0; x<volumesX; x++)
            {
                for(Int_t y=0; y<volumesY; y++)
                {
                    VisibleSpectrumH[systematic][AD][x][y] = new TH1D(Form("VisibleHAD%dCell%d,%d,%d",AD,x,y,systematic),Form("VisibleHAD%dCell%d,%d,%d",AD,x,y,systematic),NominalResponseMatrixH[0][0]->GetYaxis()->GetNbins(),0,12);
                    
                    for(Int_t i=0; i<NominalResponseMatrixH[0][0]->GetYaxis()->GetNbins(); i++)//visible
                    {
                        for(Int_t j=0; j<NominalResponseMatrixH[0][0]->GetXaxis()->GetNbins()-(1.8/0.05); j++)//true
                        {
                            //[True_0 ... True_n] = [0,0]  [...] [0,m]  * [Vis_0]
                            //                      [0,n]  [...] [n,m]     [...]
                            //                                            [Vis_m]
                            
                            // [nx1] = [n x m] x [mx1]
                            
                            
                            VisibleSpectrumH[systematic][AD][x][y]->SetBinContent(i+1,VisibleSpectrumH[systematic][AD][x][y]->GetBinContent(i+1)+Histo[systematic][AD][x][y]->GetBinContent(j+1+(1.8/0.05),i+1)*TrueADH->Interpolate(Histo[systematic][AD][x][y]->GetXaxis()->GetBinCenter(j+1+(1.8/0.05))));
                            
                            //if(x==0&&y==0&&i==0)
                            //{
                            //                            std::cout << "energy : " << Histo[systematic][AD][x][y]->GetXaxis()->GetBinCenter(j+1+(1.8/0.05)) << " in bin " << j+1 << std::endl;
                            //}
                            // std::cout << " Energy bin j " << j << " - " << NominalResponseMatrixH[x][y]->GetXaxis()->GetBinCenter(j+1+(1.8/0.05)) << std::endl;
                            
                        }
                    }
                    
                    CanvasSys[systematic]->cd(AD+NADs*x+1);
                    
                    VisibleSpectrumH[systematic][AD][x][y]->Draw();
                }
            }
        }
        
        TString FileNameSave;
        
        switch (systematic) {
            case 1:
                FileNameSave = "AllSystematics";
                break;
            case 2:
                FileNameSave = "IAV";
                break;
            case 3:
                FileNameSave = "NL";
                break;
            case 4:
                FileNameSave = "Nominal";
                break;
            case 5:
                FileNameSave = "OAV";
                break;
            case 6:
                FileNameSave = "RelativeEnergyScale";
                break;
            case 7:
                FileNameSave = "Resolution";
                break;
                
            default:
                break;
        }
 
        CanvasSys[systematic]->Print(FileNameSave+".eps");
        CanvasSys[systematic]->Close();
    }
    

    TH1D* F_obs[NADs][15][volumesX][volumesY];
    TH1D* RebinedF_obs[NADs][15][volumesX][volumesY];
    TH1D* F_pred[NADs][15][volumesX][volumesY];
    TH1D* RebinedF_pred[NADs][15][volumesX][volumesY];

    Double_t Cov;

    TH2D* CovH[NADs][15][volumesX];

    //Calculate covariance matrix:
    TCanvas* CanvasCovSys[15];
    std::cout << "Get binning " << std::endl;

    const Int_t n_evis_bins=34;
    Double_t evis_bins[n_evis_bins];
    
    evis_bins[0] = 1.5;
    
    for (Int_t i = 0; i < n_evis_bins-1; i++)
    {
        evis_bins[i+1] = 0.2 * i + 1.6;
        
        std::cout << "EVIS: " << evis_bins[i+1] << std::endl;
    }
    evis_bins[n_evis_bins] = 12;
    
    std::cout << "Calculate cov " << std::endl;

    for (Int_t systematic = 1; systematic<num; systematic++)
    {
        CanvasCovSys[systematic] = new TCanvas(Form("CanvasCovSys%d",systematic),Form("CanvasCovSys%d",systematic));
        
        CanvasCovSys[systematic]->Divide(NADs,volumesX);
        
        for(Int_t AD=0; AD<NADs; AD++)
        {
            for(Int_t x=0; x<volumesX; x++)
            {
                CovH[AD][systematic][x] = new TH2D(Form("ConvarianceMatrixAD%d_%i_volume%i",AD,systematic,x),Form("ConvarianceMatrixAD%d_%i_volume%i",AD,systematic,x),n_evis_bins,evis_bins[0],evis_bins[n_evis_bins],n_evis_bins,evis_bins[0],evis_bins[n_evis_bins]);
                
                for(Int_t y=0; y<volumesY; y++)
                {
                    std::cout << "CLONING REBINNING " << systematic << std::endl;

                    F_pred[AD][systematic][x][y] = (TH1D*)VisibleSpectrumH[4][AD][x][y]->Clone(Form("PRED_Nominal%d_%i_volume%i",AD,systematic,x));
                    RebinedF_pred[AD][systematic][x][y] = (TH1D*)F_pred[AD][systematic][x][y]->Rebin(n_evis_bins,Form("Rebined_Nominal%d_%i_volume%i",AD,systematic,x),evis_bins);
                    F_obs[AD][systematic][x][y] = (TH1D*)VisibleSpectrumH[systematic][AD][x][y]->Clone(Form("OBS_Nominal%d_%i_volume%i",AD,systematic,x));
                    RebinedF_obs[AD][systematic][x][y] = (TH1D*)F_obs[AD][systematic][x][y]->Rebin(n_evis_bins,Form("Rebined_Varied%d_%i_volume%i",AD,systematic,x),evis_bins);

                    for(Int_t a=0; a<n_evis_bins; a++)
                    {
                        for(Int_t b=0; b<n_evis_bins; b++)
                        {
                            Cov=(RebinedF_pred[AD][systematic][x][y]->GetBinContent(a+1) - RebinedF_obs[AD][systematic][x][y]->GetBinContent(a+1))*(RebinedF_pred[AD][systematic][x][y]->GetBinContent(b+1) - RebinedF_obs[AD][systematic][x][y]->GetBinContent(b+1));
                            
                            CovH[AD][systematic][x]->SetBinContent(a+1,b+1,Cov);
                        }
                    }
                }
                
                CanvasCovSys[systematic]->cd(AD+NADs*x+1);
                
                CovH[AD][systematic][x]->Draw("colz");
            }
        }
        
        TString FileNameSave2;

        switch (systematic) {
            case 1:
                FileNameSave2 = "AllSystematics";
                break;
            case 2:
                FileNameSave2 = "IAV";
                break;
            case 3:
                FileNameSave2 = "NL";
                break;
            case 4:
                FileNameSave2 = "Nominal";
                break;
            case 5:
                FileNameSave2 = "OAV";
                break;
            case 6:
                FileNameSave2 = "RelativeEnergyScale";
                break;
            case 7:
                FileNameSave2 = "Resolution";
                break;
                
            default:
                break;
        }
        
        CanvasCovSys[systematic]->Print("CovarianceMatrix"+FileNameSave2+".eps");
        CanvasCovSys[systematic]->Close();
    }
    
    TCanvas* CanvasRebin[15];
    
    std::cout << "PRINT REBINNED " << std::endl;
    
    for (Int_t systematic = 1; systematic<num; systematic++)
    {
        CanvasRebin[systematic] = new TCanvas(Form("CanvasSysRebin%d",systematic),Form("CanvasSysRebin%d",systematic));
        
        CanvasRebin[systematic]->Divide(NADs,volumesX);
        
        for(Int_t AD=0; AD<NADs; AD++)
        {
            for(Int_t x=0; x<volumesX; x++)
            {
                for(Int_t y=0; y<volumesY; y++)
                {
                    CanvasRebin[systematic]->cd(AD+NADs*x+1);
                    
                    RebinedF_pred[AD][systematic][x][y]->Draw();
                }
            }
        }
        
        TString FileNameSave3;
        
        switch (systematic) {
            case 1:
                FileNameSave3 = "RebinAllSystematics";
                break;
            case 2:
                FileNameSave3 = "RebinIAV";
                break;
            case 3:
                FileNameSave3 = "RebinNL";
                break;
            case 4:
                FileNameSave3 = "RebinNominal";
                break;
            case 5:
                FileNameSave3 = "RebinOAV";
                break;
            case 6:
                FileNameSave3 = "RebinRelativeEnergyScale";
                break;
            case 7:
                FileNameSave3 = "RebinResolution";
                break;
                
            default:
                break;
        }
        
        CanvasRebin[systematic]->Print(FileNameSave3+".eps");
        CanvasRebin[systematic]->Close();
    }
}