#include <iostream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include <string>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TF1.h"
#include "../NominalData.h"

//  Needs improvement (the fits to the chi2 aren't optimized and produce problematic results)
//
int plot_range()
{
    TH1::AddDirectory(kFALSE);
    
    Double_t MaxAccuracy[4];
    const Int_t NExperiments = 21;
    Double_t Accuracy[4][NExperiments];
    Int_t markerStyle=24;
    Int_t markerColor=33;
    Int_t lineStyle=2;
    NominalData* Data = new NominalData(0,2);//nGd,DataSet = P12E
    
    Double_t BestDeltaM = Data->GetDm2ee();
    Double_t BestSin22t13 = Data->GetSin22t13();
    
    Int_t NWeeks = Data->GetWeeks();
    
    TH1D* DeltaSinChiSquare[NExperiments];
    TH1D* DeltaDeltaChiSquare[NExperiments];
    TH1D* SinDeltaChiSquare[NExperiments];
    TH1D* SinSinChiSquare[NExperiments];
    
    TGraphAsymmErrors* SinSin;
    TGraphAsymmErrors* SinDelta;
    TGraphAsymmErrors* DeltaSin;
    TGraphAsymmErrors* DeltaDelta;
    
    Double_t YSinSin[NExperiments];
    Double_t YSinSinError1[NExperiments];
    Double_t YSinSinError2[NExperiments];
    
    Double_t YSinDelta[NExperiments];
    Double_t YSinDeltaError1[NExperiments];
    Double_t YSinDeltaError2[NExperiments];
    
    Double_t YDeltaSin[NExperiments];
    Double_t YDeltaSinError1[NExperiments];
    Double_t YDeltaSinError2[NExperiments];
    
    Double_t YDeltaDelta[NExperiments];
    Double_t YDeltaDeltaError1[NExperiments];
    Double_t YDeltaDeltaError2[NExperiments];
    
    Double_t MinChi2[NExperiments][4];
    
    TFile* PlotF[4];
    
    Double_t XSin22t13[NExperiments];
    Double_t XDm2ee[NExperiments];
    
    TCanvas* c[4];
    
    for (Int_t i = 0; i<NExperiments; i++)
    {
        XSin22t13[i] = i*0.2/(NExperiments-1);
        
        XDm2ee[i] = 0.0015+i*(0.002/(NExperiments-1));
    }
    
    Double_t minX[NExperiments];
    Double_t maxX[NExperiments];
    
    for (Int_t count = 0; count<4; count++)
    {
        MaxAccuracy[count] = 0;

        if (count==0)//  Fit Sin22t13 varying Sin22t13 for a fixed DeltaM
        {
            TH1D* FitSinSinChiSquareResult;
            
            for (Int_t week = 0; week<NWeeks; week++)
            {
                for (Int_t i = 0; i<NExperiments; i++)
                {
                    PlotF[count] = TFile::Open("../ChiSquare/Gadolinium/Combine2/Sin22t13Sin22t13Range.root");
                    
                    SinSinChiSquare[i]=(TH1D*)gDirectory->Get(Form("Sin%f_Distribution_DeltaM%f_Period%d",XSin22t13[i],BestDeltaM,week));
                    
                    Double_t Minimum = 100000000.;
                    
                    Double_t Originalwidth = (SinSinChiSquare[i]->GetXaxis()->GetXmax()-SinSinChiSquare[i]->GetXaxis()->GetXmin())/((SinSinChiSquare[i]->GetXaxis()->GetNbins()-1));
                    
                    for (Int_t j = 0; j < SinSinChiSquare[i]->GetXaxis()->GetNbins(); j++)
                    {
                        if(SinSinChiSquare[i]->GetBinContent(j+1)<Minimum)
                        {
                            Minimum = SinSinChiSquare[i]->GetBinContent(j+1);
                            YSinSin[i] = SinSinChiSquare[i]->GetXaxis()->GetXmin()+j*Originalwidth;
                        }
                    }
                    
//                    std::cout << "MIN CHI" << Minimum << std::endl;
                    
                    minX[i] = 20.;
                    maxX[i] = 0.2;
                    Double_t xvalue;
                    Double_t CHISINMIN=1000000;
                    Double_t SINMIN = 0.2;
                    FitSinSinChiSquareResult = new TH1D(Form("Interpolated Sin%f_Distribution_DeltaM%f_Period%d",XSin22t13[i],BestDeltaM,week),Form("Interpolated Sin%f_Distribution_DeltaM%f_Period%d",XSin22t13[i],BestDeltaM,week),SinSinChiSquare[i]->GetXaxis()->GetNbins()*100,0,0.2);
                    
                    Double_t width = (SinSinChiSquare[i]->GetXaxis()->GetXmax()-SinSinChiSquare[i]->GetXaxis()->GetXmin())/((SinSinChiSquare[i]->GetXaxis()->GetNbins()-1)*100);
                    
                    //                    std::cout << "Width" << SinSinChiSquare[i]->GetXaxis()->GetXmax() << " "  << SinSinChiSquare[i]->GetXaxis()->GetXmin() <<  " " << SinSinChiSquare[i]->GetXaxis()->GetNbins() << " " << width << std::endl;
                    
                    for (Int_t j = 0; j<=((SinSinChiSquare[i]->GetXaxis()->GetNbins()-1)*100); j++)
                    {
                        xvalue = SinSinChiSquare[i]->GetXaxis()->GetXmin()+j*width;
                        
                        FitSinSinChiSquareResult->SetBinContent(j+1,SinSinChiSquare[i]->Interpolate(xvalue)-Minimum);
                        
                        if(FitSinSinChiSquareResult->GetBinContent(j+1)<CHISINMIN)
                        {
                            CHISINMIN = FitSinSinChiSquareResult->GetBinContent(j+1);
                            SINMIN=xvalue;
                            YSinSin[i] = xvalue;
                        }
                    }
                    
                    //                    bool minFound = 1;
                    
                    for (Int_t j = 0; j<=((SinSinChiSquare[i]->GetXaxis()->GetNbins()-1)*100); j++)
                    {
                        xvalue = SinSinChiSquare[i]->GetXaxis()->GetXmin()+j*width;
                        
                        //                        std::cout << xvalue << std::endl;
                        
                        //                        std::cout << xvalue << " " << ((TMath::Floor(FitSinSinChiSquareResult->GetBinContent(j+1))==0)) << std::endl;
                        
                        if((TMath::Floor(FitSinSinChiSquareResult->GetBinContent(j+1))==0)&&(xvalue<=YSinSin[i])&&(minX[i]>=xvalue))//Look for minimum X where Δχ2 == 1
                        {
                            
                            //                            std::cout << xvalue << " " << ((TMath::Floor(FitSinSinChiSquareResult->GetBinContent(j+1))==0)&&(minX<=xvalue)) << std::endl;
                            
                            
                            minX[i] = xvalue;
                            //                            minFound = 0;
                            //                            std::cout << "ERROR 1 " << minX[i] << std::endl;
                            
                        }
                        else if((TMath::Floor(FitSinSinChiSquareResult->GetBinContent(j+1))==0)&&(xvalue>YSinSin[i]))//Look for maximum X where Δχ2 == 1
                        {
                            maxX[i] = xvalue;
                            
                            //                            std::cout << maxX << std::endl;
                            
                        }
                        
                    }
                    
                    YSinSinError1[i] = TMath::Abs(YSinSin[i]-minX[i]);// First crossing with 1
                    //                    std::cout << "YSINERROR1 " << YSinSinError1[i] << std::endl;
                    
                    YSinSinError2[i] = TMath::Abs(maxX[i]-YSinSin[i]);// Second crossing with 1
                    //                    std::cout << "YSINERROR2 " << YSinSinError2[i] << std::endl;
                    
                    delete FitSinSinChiSquareResult;
                    
                    Accuracy[count][i] = TMath::Abs((YSinSin[i] - XSin22t13[i])/XSin22t13[i]);
                    
                    if(MaxAccuracy[count] < Accuracy[count][i])
                    {
                        MaxAccuracy[count] = Accuracy[count][i];
                    }
                }
            }
            
            SinSin = new TGraphAsymmErrors(NExperiments,XSin22t13,YSinSin,0,0,YSinSinError1,YSinSinError2);
            
            c[count] = new TCanvas("Sin22t13 Sin22t13 Range Fit","Sin22t13 Sin22t13 Range Fit",500,500);
            
            gPad->SetLeftMargin(0.2);
            
            SinSin->SetTitle("");
            SinSin->SetLineColor(2);
            SinSin->SetLineWidth(1);
            SinSin->SetLineStyle(lineStyle);
            SinSin->SetMarkerColor(markerColor);
            SinSin->SetMarkerStyle(markerStyle);
            SinSin->GetXaxis()->SetTitle("true sin^{2}(2#theta_{13})");
            SinSin->GetYaxis()->SetTitle("best fit sin^{2}(2#theta_{13})");
            SinSin->GetYaxis()->SetTitleOffset(1.7);
            SinSin->Draw("ALP>");
            
            c[count]->Print("../Images/FitterTests/SinSinRange.eps");
            c[count]->Close();
            PlotF[count]->Close();
            
        }
        else if(count==1)// Fit DeltaM varying Sin22t13 for a fixed DeltaM
        {
            TH1D* FitDeltaSinChiSquareResult;
            
            for (Int_t week = 0; week<NWeeks; week++)
            {
                for (Int_t i = 0; i<NExperiments; i++)
                {
                    
                    PlotF[count] = TFile::Open("../ChiSquare/Gadolinium/Combine2/DeltaMSin22t13Range.root");
                    
                    DeltaSinChiSquare[i]=(TH1D*)gDirectory->Get(Form("Sin%f_Distribution_DeltaM%f_Period%d",XSin22t13[i],BestDeltaM,week));
                    
                    Double_t Minimum = 100000000.;
                    
                    Double_t Originalwidth = (DeltaSinChiSquare[i]->GetXaxis()->GetXmax()-DeltaSinChiSquare[i]->GetXaxis()->GetXmin())/((DeltaSinChiSquare[i]->GetXaxis()->GetNbins()-1));
                    
                    
                    for (Int_t j = 0; j < DeltaSinChiSquare[i]->GetXaxis()->GetNbins(); j++)
                    {
                        if(DeltaSinChiSquare[i]->GetBinContent(j+1)<Minimum)
                        {
                            Minimum = DeltaSinChiSquare[i]->GetBinContent(j+1);
                            YDeltaSin[i] = DeltaSinChiSquare[i]->GetXaxis()->GetXmin()+j*Originalwidth;
                        }
                    }
                    
//                    std::cout << "MIN CHI" << Minimum << std::endl;
                    
                    minX[i] = 20.;
                    maxX[i] = 0.2;
                    Double_t xvalue;
                    Double_t CHISINMIN=1000000;
                    Double_t DELTAMIN = 0.0035;
                    FitDeltaSinChiSquareResult = new TH1D(Form("Interpolated Sin%f_Distribution_DeltaM%f_Period%d",XSin22t13[i],BestDeltaM,week),Form("Interpolated Sin%f_Distribution_DeltaM%f_Period%d",XSin22t13[i],BestDeltaM,week),DeltaSinChiSquare[i]->GetXaxis()->GetNbins()*100,0,0.2);
                    
                    Double_t width = (DeltaSinChiSquare[i]->GetXaxis()->GetXmax()-DeltaSinChiSquare[i]->GetXaxis()->GetXmin())/((DeltaSinChiSquare[i]->GetXaxis()->GetNbins()-1)*100);
                    
                    //                    std::cout << "Width" << DeltaSinChiSquare[i]->GetXaxis()->GetXmax() << " "  << DeltaSinChiSquare[i]->GetXaxis()->GetXmin() <<  " " << DeltaSinChiSquare[i]->GetXaxis()->GetNbins() << " " << width << std::endl;
                    
                    for (Int_t j = 0; j<=((DeltaSinChiSquare[i]->GetXaxis()->GetNbins()-1)*100); j++)
                    {
                        xvalue = DeltaSinChiSquare[i]->GetXaxis()->GetXmin()+j*width;
                        
                        FitDeltaSinChiSquareResult->SetBinContent(j+1,DeltaSinChiSquare[i]->Interpolate(xvalue)-Minimum);
                        
                        if(FitDeltaSinChiSquareResult->GetBinContent(j+1)<CHISINMIN)
                        {
                            CHISINMIN = FitDeltaSinChiSquareResult->GetBinContent(j+1);
                            DELTAMIN=xvalue;
                            YDeltaSin[i] = xvalue;
                        }
                    }
                    
                    //                    bool minFound = 1;
                    
                    for (Int_t j = 0; j<=((DeltaSinChiSquare[i]->GetXaxis()->GetNbins()-1)*100); j++)
                    {
                        xvalue = DeltaSinChiSquare[i]->GetXaxis()->GetXmin()+j*width;
                        
                        //                        std::cout << xvalue << std::endl;
                        
                        //                        std::cout << xvalue << " " << ((TMath::Floor(FitDeltaSinChiSquare->GetBinContent(j+1))==0)) << std::endl;
                        
                        if((TMath::Floor(FitDeltaSinChiSquareResult->GetBinContent(j+1))==0)&&(xvalue<=YDeltaSin[i])&&(minX[i]>=xvalue))//Look for minimum X where Δχ2 == 1
                        {
                            
                            //                            std::cout << xvalue << " " << ((TMath::Floor(FitDeltaSinChiSquare->GetBinContent(j+1))==0)&&(minX<=xvalue)) << std::endl;
                            
                            
                            minX[i] = xvalue;
                            //                            minFound = 0;
                            //                            std::cout << "ERROR 1 " << minX[i] << std::endl;
                            
                        }
                        else if((TMath::Floor(FitDeltaSinChiSquareResult->GetBinContent(j+1))==0)&&(xvalue>YDeltaSin[i]))//Look for maximum X where Δχ2 == 1
                        {
                            maxX[i] = xvalue;
                            
                            //                            std::cout << maxX << std::endl;
                            
                        }
                        
                    }
                    
                    YDeltaSinError1[i] = TMath::Abs(YDeltaSin[i]-minX[i]);// First crossing with 1
                    //                    std::cout << "YSINERROR1 " << YSinSinError1[i] << std::endl;
                    
                    YDeltaSinError2[i] = TMath::Abs(maxX[i]-YDeltaSin[i]);// Second crossing with 1
                    //                    std::cout << "YSINERROR2 " << YSinSinError2[i] << std::endl;
                    
                    delete FitDeltaSinChiSquareResult;
                    
                    Accuracy[count][i] = TMath::Abs((YDeltaSin[i] - BestDeltaM)/BestDeltaM);
                    
                    if(MaxAccuracy[count] < Accuracy[count][i])
                    {
                        MaxAccuracy[count] = Accuracy[count][i];
                    }
                }
            }
            
            c[count] = new TCanvas("Delta Sin22t13 Range Fit","Delta Sin22t13 Range Fit",500,500);
            gPad->SetLeftMargin(0.2);
            DeltaSin = new TGraphAsymmErrors(NExperiments,XSin22t13,YDeltaSin,0,0,YDeltaSinError1,YDeltaSinError2);
            
            DeltaSin->SetTitle("");
            DeltaSin->SetLineColor(2);
            DeltaSin->SetLineWidth(1);
            DeltaSin->SetLineStyle(lineStyle);
            DeltaSin->SetMarkerColor(markerColor);
            DeltaSin->SetMarkerStyle(markerStyle);
            DeltaSin->GetYaxis()->SetTitleOffset(2.2);
            DeltaSin->GetXaxis()->SetTitle("true sin^{2}(2#theta_{13})");
            DeltaSin->GetYaxis()->SetTitle("best fit #Deltam^{2}_{ee} (eV^{2})");
            DeltaSin->Draw("ALP>");
            
            c[count]->Print("../Images/FitterTests/DeltaSinRange.eps");
            c[count]->Close();
            
            PlotF[count]->Close();
        }
        else if(count==2)//  Fit Sin22t13 varying DeltaM for a fixed Sin22t13
        {
            TH1D* FitSinDeltaChiSquareResult;
            
            for (Int_t week = 0; week<NWeeks; week++)
            {
                for (Int_t i = 0; i<NExperiments; i++)
                {
                    PlotF[count] = TFile::Open("../ChiSquare/Gadolinium/Combine2/Sin22t13Delta2MeeRange.root");
                    
                    SinDeltaChiSquare[i]=(TH1D*)gDirectory->Get(Form("Sin%f_Distribution_DeltaM%f_Period%d",BestSin22t13,XDm2ee[i],week));
                    
                    Double_t Minimum = 100000000.;
                    
                    Double_t Originalwidth = (SinDeltaChiSquare[i]->GetXaxis()->GetXmax()-SinDeltaChiSquare[i]->GetXaxis()->GetXmin())/((SinDeltaChiSquare[i]->GetXaxis()->GetNbins()-1));
                    
                    
                    for (Int_t j = 0; j < SinDeltaChiSquare[i]->GetXaxis()->GetNbins(); j++)
                    {
                        if(SinDeltaChiSquare[i]->GetBinContent(j+1)<Minimum)
                        {
                            Minimum = SinDeltaChiSquare[i]->GetBinContent(j+1);
                            YSinDelta[i] = SinDeltaChiSquare[i]->GetXaxis()->GetXmin()+j*Originalwidth;
                        }
                    }
                    
//                    std::cout << "MIN CHI" << Minimum << std::endl;
                    
                    minX[i] = 20.;
                    maxX[i] = 0.2;
                    Double_t xvalue;
                    Double_t CHISINMIN=1000000;
                    Double_t SINMIN = 0.2;
                    FitSinDeltaChiSquareResult = new TH1D(Form("Interpolated Sin%f_Distribution_DeltaM%f_Period%d",BestSin22t13,XDm2ee[i],week),Form("Interpolated Sin%f_Distribution_DeltaM%f_Period%d",BestSin22t13,XDm2ee[i],week),SinDeltaChiSquare[i]->GetXaxis()->GetNbins()*100,0,0.2);
                    
                    Double_t width = (SinDeltaChiSquare[i]->GetXaxis()->GetXmax()-SinDeltaChiSquare[i]->GetXaxis()->GetXmin())/((SinDeltaChiSquare[i]->GetXaxis()->GetNbins()-1)*100);
                    
                    //                    std::cout << "Width" << SinDeltaChiSquare[i]->GetXaxis()->GetXmax() << " "  << SinDeltaChiSquare[i]->GetXaxis()->GetXmin() <<  " " << SinDeltaChiSquare[i]->GetXaxis()->GetNbins() << " " << width << std::endl;
                    
                    for (Int_t j = 0; j<=((SinDeltaChiSquare[i]->GetXaxis()->GetNbins()-1)*100); j++)
                    {
                        xvalue = SinDeltaChiSquare[i]->GetXaxis()->GetXmin()+j*width;
                        
                        FitSinDeltaChiSquareResult->SetBinContent(j+1,SinDeltaChiSquare[i]->Interpolate(xvalue)-Minimum);
                        
                        if(FitSinDeltaChiSquareResult->GetBinContent(j+1)<CHISINMIN)
                        {
                            CHISINMIN = FitSinDeltaChiSquareResult->GetBinContent(j+1);
                            SINMIN=xvalue;
                            YSinDelta[i] = xvalue;
                        }
                    }
                    
                    //                    bool minFound = 1;
                    
                    for (Int_t j = 0; j<=((SinDeltaChiSquare[i]->GetXaxis()->GetNbins()-1)*100); j++)
                    {
                        xvalue = SinDeltaChiSquare[i]->GetXaxis()->GetXmin()+j*width;
                        
                        //                        std::cout << xvalue << std::endl;
                        
                        //                        std::cout << xvalue << " " << ((TMath::Floor(FitSinDeltaChiSquareResult->GetBinContent(j+1))==0)) << std::endl;
                        
                        if((TMath::Floor(FitSinDeltaChiSquareResult->GetBinContent(j+1))==0)&&(xvalue<=YSinDelta[i])&&(minX[i]>=xvalue))//Look for minimum X where Δχ2 == 1
                        {
                            
                            //                            std::cout << xvalue << " " << ((TMath::Floor(FitSinDeltaChiSquare->GetBinContent(j+1))==0)&&(minX<=xvalue)) << std::endl;
                            
                            
                            minX[i] = xvalue;
                            //                            minFound = 0;
                            //                            std::cout << "ERROR 1 " << minX[i] << std::endl;
                            
                        }
                        else if((TMath::Floor(FitSinDeltaChiSquareResult->GetBinContent(j+1))==0)&&(xvalue>YSinDelta[i]))//Look for maximum X where Δχ2 == 1
                        {
                            maxX[i] = xvalue;
                            
                            //                            std::cout << maxX << std::endl;
                            
                        }
                        
                    }
                    
                    YSinDeltaError1[i] = TMath::Abs(YSinDelta[i]-minX[i]);// First crossing with 1
                    //                    std::cout << "YSINERROR1 " << YSinDeltaError1[i] << std::endl;
                    
                    YSinDeltaError2[i] = TMath::Abs(maxX[i]-YSinDelta[i]);// Second crossing with 1
                    //                    std::cout << "YSINERROR2 " << YSinDeltaError2[i] << std::endl;
                    
                    delete FitSinDeltaChiSquareResult;
                    
                    Accuracy[count][i] = TMath::Abs((YSinDelta[i] - BestSin22t13)/BestSin22t13);
                    
                    if(MaxAccuracy[count] < Accuracy[count][i])
                    {
                        MaxAccuracy[count] = Accuracy[count][i];
                    }
                }
            }
            
            SinDelta = new TGraphAsymmErrors(NExperiments,XDm2ee,YSinDelta,0,0,YSinDeltaError1,YSinDeltaError2);
            
            c[count] = new TCanvas("Sin22t13 Delta Range Fit","Sin22t13 Delta Range Fit",500,500);
            
            gPad->SetLeftMargin(0.2);
            
            SinDelta->SetTitle("");
            SinDelta->SetLineColor(2);
            SinDelta->SetLineWidth(1);
            SinDelta->SetLineStyle(lineStyle);
            SinDelta->SetMarkerColor(markerColor);
            SinDelta->SetMarkerStyle(markerStyle);
            SinDelta->GetYaxis()->SetTitleOffset(1.7);
            SinDelta->GetXaxis()->SetTitle("true #Deltam^{2}_{ee} (eV^{2})");
            SinDelta->GetYaxis()->SetTitle("best fit sin^{2}(2#theta_{13})");
            SinDelta->Draw("ALP>");
            
            c[count]->Print("../Images/FitterTests/SinDeltaRange.eps");
            c[count]->Close();
            
            PlotF[count]->Close();
            
        }
        else if(count==3)//  Fit DeltaM varying DeltaM for a fixed Sin22t13
        {
            
            TH1D* FitDeltaDeltaChiSquareResult;
            
            for (Int_t week = 0; week<NWeeks; week++)
            {
                for (Int_t i = 0; i<NExperiments; i++)
                {
                    PlotF[count] = TFile::Open("../ChiSquare/Gadolinium/Combine2/DeltaMDelta2MeeRange.root");
                    
                    DeltaDeltaChiSquare[i]=(TH1D*)gDirectory->Get(Form("Sin%f_Distribution_DeltaM%f_Period%d",BestSin22t13,XDm2ee[i],week));
                    
                    Double_t Minimum = 100000000.;
                    
                    Double_t Originalwidth = (DeltaDeltaChiSquare[i]->GetXaxis()->GetXmax()-DeltaDeltaChiSquare[i]->GetXaxis()->GetXmin())/((DeltaDeltaChiSquare[i]->GetXaxis()->GetNbins()-1));
                    
                    
                    for (Int_t j = 0; j < DeltaDeltaChiSquare[i]->GetXaxis()->GetNbins(); j++)
                    {
                        if(DeltaDeltaChiSquare[i]->GetBinContent(j+1)<Minimum)
                        {
                            Minimum = DeltaDeltaChiSquare[i]->GetBinContent(j+1);
                            YDeltaDelta[i] = DeltaDeltaChiSquare[i]->GetXaxis()->GetXmin()+j*Originalwidth;
                        }
                    }
                    
//                    std::cout << "MIN CHI" << Minimum << std::endl;
                    
                    minX[i] = 20.;
                    maxX[i] = 0.2;
                    Double_t xvalue;
                    Double_t CHISINMIN=1000000;
                    Double_t DELTAMIN = 0.0035;
                    FitDeltaDeltaChiSquareResult = new TH1D(Form("Interpolated Sin%f_Distribution_DeltaM%f_Period%d",XSin22t13[i],BestDeltaM,week),Form("Interpolated Sin%f_Distribution_DeltaM%f_Period%d",XSin22t13[i],BestDeltaM,week),DeltaDeltaChiSquare[i]->GetXaxis()->GetNbins()*100,0,0.2);
                    
                    Double_t width = (DeltaDeltaChiSquare[i]->GetXaxis()->GetXmax()-DeltaDeltaChiSquare[i]->GetXaxis()->GetXmin())/((DeltaDeltaChiSquare[i]->GetXaxis()->GetNbins()-1)*100);
                    
                    //                    std::cout << "Width" << SinDeltaChiSquare[i]->GetXaxis()->GetXmax() << " "  << SinDeltaChiSquare[i]->GetXaxis()->GetXmin() <<  " " << SinDeltaChiSquare[i]->GetXaxis()->GetNbins() << " " << width << std::endl;
                    
                    for (Int_t j = 0; j<=((DeltaDeltaChiSquare[i]->GetXaxis()->GetNbins()-1)*100); j++)
                    {
                        xvalue = DeltaDeltaChiSquare[i]->GetXaxis()->GetXmin()+j*width;
                        
                        FitDeltaDeltaChiSquareResult->SetBinContent(j+1,DeltaDeltaChiSquare[i]->Interpolate(xvalue)-Minimum);
                        
                        if(FitDeltaDeltaChiSquareResult->GetBinContent(j+1)<CHISINMIN)
                        {
                            CHISINMIN = FitDeltaDeltaChiSquareResult->GetBinContent(j+1);
                            DELTAMIN=xvalue;
                            YDeltaDelta[i] = xvalue;
                        }
                    }
                    
                    //                    bool minFound = 1;
                    
                    for (Int_t j = 0; j<=((DeltaDeltaChiSquare[i]->GetXaxis()->GetNbins()-1)*100); j++)
                    {
                        xvalue = DeltaDeltaChiSquare[i]->GetXaxis()->GetXmin()+j*width;
                        
                        //                        std::cout << xvalue << std::endl;
                        
                        //                        std::cout << xvalue << " " << ((TMath::Floor(FitSinDeltaChiSquareResult->GetBinContent(j+1))==0)) << std::endl;
                        
                        if((TMath::Floor(FitDeltaDeltaChiSquareResult->GetBinContent(j+1))==0)&&(xvalue<=YDeltaDelta[i])&&(minX[i]>=xvalue))//Look for minimum X where Δχ2 == 1
                        {
                            
                            //                            std::cout << xvalue << " " << ((TMath::Floor(FitSinDeltaChiSquare->GetBinContent(j+1))==0)&&(minX<=xvalue)) << std::endl;
                            
                            
                            minX[i] = xvalue;
                            //                            minFound = 0;
                            //                            std::cout << "ERROR 1 " << minX[i] << std::endl;
                            
                        }
                        else if((TMath::Floor(FitDeltaDeltaChiSquareResult->GetBinContent(j+1))==0)&&(xvalue>YDeltaDelta[i]))//Look for maximum X where Δχ2 == 1
                        {
                            maxX[i] = xvalue;
                            
                            //                            std::cout << maxX << std::endl;
                            
                        }
                        
                    }
                    
                    YDeltaDeltaError1[i] = TMath::Abs(YDeltaDelta[i]-minX[i]);// First crossing with 1
                    //                    std::cout << "YSINERROR1 " << YSinDeltaError1[i] << std::endl;
                    
                    YDeltaDeltaError2[i] = TMath::Abs(maxX[i]-YDeltaDelta[i]);// Second crossing with 1
                    //                    std::cout << "YSINERROR2 " << YSinDeltaError2[i] << std::endl;
                    
                    delete FitDeltaDeltaChiSquareResult;
                    
                    Accuracy[count][i] = TMath::Abs((YDeltaDelta[i] - XDm2ee[i])/XDm2ee[i]);
                    
                    if(MaxAccuracy[count] < Accuracy[count][i])
                    {
                        MaxAccuracy[count] = Accuracy[count][i];
                    }
                }
            }
            
            DeltaDelta = new TGraphAsymmErrors(NExperiments,XDm2ee,YDeltaDelta,0,0,YDeltaDeltaError1,YDeltaDeltaError2);
            
            c[count] = new TCanvas("Sin22t13 Delta Range Fit","Sin22t13 Delta Range Fit",500,500);
            
            gPad->SetLeftMargin(0.2);
            
            DeltaDelta->SetTitle("");
            DeltaDelta->SetLineColor(2);
            DeltaDelta->SetLineWidth(1);
            DeltaDelta->SetLineStyle(lineStyle);
            DeltaDelta->SetMarkerColor(markerColor);
            DeltaDelta->SetMarkerStyle(markerStyle);
            DeltaDelta->GetYaxis()->SetTitleOffset(2.2);
            DeltaDelta->GetXaxis()->SetTitle("true #Deltam^{2}_{ee} (eV^{2})");
            DeltaDelta->GetYaxis()->SetTitle("best fit #Deltam^{2}_{ee} (eV^{2})");
            DeltaDelta->Draw("ALP>");
            
            c[count]->Print("../Images/FitterTests/DeltaDeltaRange.eps");
            c[count]->Close();
            
            PlotF[count]->Close();
        }
    }
    
    for(Int_t count = 0; count<4; count++)
    {
        for(Int_t i = 0; i<NExperiments; i++)
        {
            std::cout << " Count " << count << " Experiment " << i << " Accuracy : " <<  Accuracy[count][i] << std::endl;
        }
    }
    
    std::cout << "Accuracy of the fitter for" << std::endl << "    SinSin:" << MaxAccuracy[0] << std::endl << "    DeltaSin: " << MaxAccuracy[1] << std::endl << "    SinDelta: "  << MaxAccuracy[2] << std::endl << "   DeltaDelta: " << MaxAccuracy[3] << std::endl;
    
    return 0;
}

