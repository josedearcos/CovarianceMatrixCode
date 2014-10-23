void plot_covariance_matrix()
{
    TH1::AddDirectory(kFALSE);
    
    gStyle->SetErrorX(0.0001);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.05,"x");
    gStyle->SetTitleSize(0.05,"y");
    
    //To draw using a better palette:
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    
    TH2D* MatrixTotalH[1];
    TH2D* MatrixIAVH[1];
    TH2D* MatrixNLH[1];
    TH2D* MatrixResoH[1];
    TH2D* MatrixRelativeEnergyH[1];
    TH2D* MatrixPowerH[1];
    TH2D* MatrixIsotopeH[1];
    TH2D* MatrixSystematicH[1];
    TH2D* MatrixBackgroundH[1];
    TH2D* MatrixDistortAmC[1];
    TH2D* MatrixDistortLiHe[1];
    TH2D* MatrixDistortFN[1];
    TH2D* MatrixVaryAmC[1];
    TH2D* MatrixVaryLiHe[1];
    TH2D* MatrixVaryFN[1];
    TH2D* MatrixVaryAcc[1];

    Char_t CorrelationName[100];
//
//    /CovarianceMatrices/Combine1/CovarianceMatricesRoot/DistortAmC.root				IAV.root				ReactorPower.root			StatisticalCovarianceMatrix.root	VaryFastNeutrons.root
//    DistortFastNeutrons.root		Isotope.root				RelativeEnergyScale.root		VaryAccidental.root			VaryLiHe.root
//    DistortLiHe.root			NL.root					Resolution.root				VaryAmC.root
//

//    sprintf(CorrelationName,"Vary Accidental Correlation Matrix%d",week);
//    sprintf(CorrelationName,"Vary LiHe Correlation Matrix%d",week);
//    sprintf(CorrelationName,"Vary FN Correlation Matrix%d",week);
//    sprintf(CorrelationName,"Distort LiHe Correlation Matrix%d",week);
//    sprintf(CorrelationName,"Distort FN Correlation Matrix%d",week);
//    sprintf(CorrelationName,"Isotope Correlation Matrix%d",week);
//    sprintf(CorrelationName,"Absolute Energy Scale Correlation Matrix%d",week);
    
    TFile* MatrixF = TFile::Open("../CovarianceMatrices/FitterCovarianceMatrixResultsPeriod0.root");
    
    TCanvas *c1 = new TCanvas("Total Cov Matrix","Total Cov Matrix",600,600);
    
    for(Int_t week=0;week<1;++week)
    {
            MatrixTotalH[week]=(TH2D*)gDirectory->Get("Total Covariance Matrix");
        
            MatrixTotalH[week]->SetTitle("Total Covariance Matrix");
        
            MatrixTotalH[week]->Draw("colz");
    }
    
    c1->Print("../Images/CovarianceMatrices/TotalCovMatrix.eps", "png");
    
    MatrixF->Close();
    
    TFile* MatrixF1 = TFile::Open("../CovarianceMatrices/Combine1/CovarianceMatricesRoot/DistortAmC.root");
    
    TCanvas *c2 = new TCanvas("Distort AMC Matrix","Distort AMC Matrix",600,600);

    for(Int_t week=0;week<1;++week)
    {
        sprintf(CorrelationName,"Distort AmC Correlation Matrix%d",week);

        MatrixDistortAmC[week]=(TH2D*)gDirectory->Get(CorrelationName);
        
        MatrixDistortAmC[week]->SetTitle(CorrelationName);
        
        MatrixDistortAmC[week]->Draw("colz");
        
    }

    c2->Print("../Images/CovarianceMatrices/DistortAmC.eps", "png");
    
    MatrixF1->Close();
    
    TFile* MatrixF2 = TFile::Open("../CovarianceMatrices/Combine1/CovarianceMatricesRoot/RelativeEnergyScale.root");
    
    TCanvas *c3 = new TCanvas("Relative Energy Scale Matrix","Relative Energy Scale Matrix",600,600);
    
    for(Int_t week=0;week<1;++week)
    {
        sprintf(CorrelationName,"Relative Energy Scale Correlation Matrix%d",week);

        MatrixRelativeEnergyH[week]=(TH2D*)gDirectory->Get(CorrelationName);
        
        MatrixRelativeEnergyH[week]->SetTitle(CorrelationName);
        
        MatrixRelativeEnergyH[week]->Draw("colz");
    }
    
    c3->Print("../Images/CovarianceMatrices/RelativeEnergyMatrix.eps", "png");
    
    MatrixF2->Close();
    
    TFile* MatrixF3 = TFile::Open("../CovarianceMatrices/Combine1/CovarianceMatricesRoot/IAV.root");
    
    TCanvas *c4 = new TCanvas("IAV Matrix","IAV Matrix",600,600);
    
    for(Int_t week=0;week<1;++week)
    {
        sprintf(CorrelationName,"IAV Correlation Matrix%d",week);
        
        MatrixIAVH[week]=(TH2D*)gDirectory->Get(CorrelationName);
        
        MatrixIAVH[week]->SetTitle(CorrelationName);
        
        MatrixIAVH[week]->Draw("colz");
    }
    
    c4->Print("../Images/CovarianceMatrices/IAVMatrix.eps", "png");
    
    MatrixF3->Close();
    
    TFile* MatrixF4 = TFile::Open("../CovarianceMatrices/Combine1/CovarianceMatricesRoot/Resolution.root");
    
    TCanvas *c5 = new TCanvas("Resolution Matrix","Resolution Matrix",600,600);
    
    for(Int_t week=0;week<1;++week)
    {
        sprintf(CorrelationName,"Resolution Correlation Matrix%d",week);
        
        MatrixResoH[week]=(TH2D*)gDirectory->Get(CorrelationName);
        
        MatrixResoH[week]->SetTitle(CorrelationName);
        
        MatrixResoH[week]->Draw("colz");
    }
    
    c5->Print("../Images/CovarianceMatrices/ResolutionMatrix.eps", "png");
    
    MatrixF4->Close();
    
    TFile* MatrixF5 = TFile::Open("../CovarianceMatrices/Combine1/CovarianceMatricesRoot/ReactorPower.root");
    
    TCanvas *c6 = new TCanvas("Power Matrix","Power Matrix",600,600);
    
    for(Int_t week=0;week<1;++week)
    {
        sprintf(CorrelationName,"Reactor Power Correlation Matrix%d",week);
        
        MatrixPowerH[week]=(TH2D*)gDirectory->Get(CorrelationName);
        
        MatrixPowerH[week]->SetTitle(CorrelationName);
        
        MatrixPowerH[week]->Draw("colz");
    }
    
    c6->Print("../Images/CovarianceMatrices/ReactorPower.eps", "png");
    
    MatrixF5->Close();
    
    TFile* MatrixF6 = TFile::Open("../CovarianceMatrices/Combine1/CovarianceMatricesRoot/NL.root");
    
    TCanvas *c7 = new TCanvas("NL Matrix","NL Matrix",600,600);
    
    for(Int_t week=0;week<1;++week)
    {
        sprintf(CorrelationName,"NL Correlation Matrix%d",week);
        
        MatrixPowerH[week]=(TH2D*)gDirectory->Get(CorrelationName);
        
        MatrixPowerH[week]->SetTitle(CorrelationName);
        
        MatrixPowerH[week]->Draw("colz");
    }
    
    c7->Print("../Images/CovarianceMatrices/NL.eps", "png");
    
    MatrixF6->Close();
    
    TFile* MatrixF7 = TFile::Open("../CovarianceMatrices/Combine1/CovarianceMatricesRoot/VaryAmC.root");
    
    TCanvas *c8 = new TCanvas("NL Matrix","NL Matrix",600,600);
    
    for(Int_t week=0;week<1;++week)
    {
        sprintf(CorrelationName,"Vary AmC Correlation Matrix%d",week);
        
        MatrixPowerH[week]=(TH2D*)gDirectory->Get(CorrelationName);
        
        MatrixPowerH[week]->SetTitle(CorrelationName);
        
        MatrixPowerH[week]->Draw("colz");
    }
    
    c8->Print("../Images/CovarianceMatrices/NL.eps", "png");
    
    MatrixF7->Close();
    
    
    

}