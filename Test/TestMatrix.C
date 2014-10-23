#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2.h"
#include "TMath.h"
#include <vector>
#include <math.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TCanvas.h"
#include "TAxis.h"
#include "TArrayD.h"

const Int_t MatrixBins = 240;
const Int_t n_etrue_bins=39;
const Int_t n_evis_bins=37;

const NormalizeRow=0;
// In this class I'm going to make a small simulation of the whole code with flat inputs to look for bugs.
void TestMatrix()
{
    
    Double_t Norma[MatrixBins];
    Double_t InitialVisibleEnergy=0.7;
    Double_t InitialEnergy=1.8;
    Double_t FinalVisibleEnergy=12;
    Double_t FinalEnergy=9.6;
    Double_t evis_bins[n_evis_bins+1];
    Double_t enu_bins[n_etrue_bins+1];
    
    TH1F* OscDeltaVisibleSpectrumH[MatrixBins];
    TH1F* PredictionTrueH;
    TH1F* RebinnedPredictionTrueH;
    TH1F* PredictionVisH[n_etrue_bins];
    TH1F* VisibleHisto;
    TH1F* FineVisibleHisto;
    TH1F* RebinnedFineVisibleHisto;
    TH1F* FinePredictionTrueH;
    
    for (Int_t i = 0; i <= n_etrue_bins; i++)
    {
        enu_bins[i] = 0.2 * i + InitialEnergy;
    }
    
    evis_bins[0] = 0.7;
    for (Int_t i = 0; i < n_evis_bins-1; i++)
    {
        evis_bins[i+1] = 0.2 * i + 1.0;
    }
    evis_bins[n_evis_bins] = FinalVisibleEnergy;
    
    // // // // // // // // // // // // // // // // // // // // // // // // //
    //                               Flat Spectrum
    // // // // // // // // // // // // // // // // // // // // // // // // //
    
    PredictionTrueH = new TH1F(Form("True Spectrum"),Form("True Spectrum"),MatrixBins,0,FinalVisibleEnergy);
 
    for (Int_t j = 0; j < MatrixBins; j++)
    {
        PredictionTrueH->SetBinContent(j+1,1);
    }
    
    RebinnedPredictionTrueH=(TH1F*)PredictionTrueH->Rebin(n_etrue_bins,Form("Rebinned True Spectrum"),enu_bins);
    FinePredictionTrueH = new TH1F(Form("Fine True Spectrum"),Form("Fine True Spectrum"),MatrixBins,0,FinalVisibleEnergy);
  
    for (Int_t j = 0; j < MatrixBins; j++)
    {
        FinePredictionTrueH->SetBinContent(j+1,0);
    }
    for (Int_t j = 36; j < 192; j++)
    {
        FinePredictionTrueH->SetBinContent(j+1,1);
    }
    // // // // // // // // // // // // // // // // // // // // // // // // //
    //                               Energy Matrix
    // // // // // // // // // // // // // // // // // // // // // // // // //
    
    for (Int_t i = 0; i < MatrixBins; i++)
    {
        OscDeltaVisibleSpectrumH[i] = new TH1F(Form("Spectrum%d",i),Form("Spectrum%d",i),MatrixBins,0,FinalVisibleEnergy);
        for (Int_t j = 0; j < MatrixBins; j++)
        {
            OscDeltaVisibleSpectrumH[i]->SetBinContent(i+1,1);//    Energy Deltas with 1's in each true energy.
        }
    }
    
    TH2F* FineMatrixH = new TH2F("FineMatrix","FineMatrix",MatrixBins,0,FinalVisibleEnergy,MatrixBins,0,FinalVisibleEnergy);
    TH2F* MatrixH = new TH2F("Matrix","Matrix",n_etrue_bins,enu_bins,n_evis_bins,evis_bins);
    FineMatrixH->Reset();
    MatrixH->Reset();
    
    //Fine version
    for (Int_t i = 0; i < MatrixBins; i++)
    {
        for (Int_t j = 0; j < MatrixBins; j++)
        {
            FineMatrixH->SetBinContent(i+1,j+1, OscDeltaVisibleSpectrumH[i]->GetBinContent(j+1));
        }
    }

    //Rebin version
    TH1F* OscDeltaVisibleSpectrumSumH[n_etrue_bins];

    for (Int_t i = 0; i < n_etrue_bins; i++)
    {
        OscDeltaVisibleSpectrumSumH[i] = new TH1F(Form("Fine Spectrum Vis for True Energy Index %d",i),Form("Fine Spectrum Vis i %d",i),MatrixBins,0,FinalVisibleEnergy);
        
        for (Int_t TrueEnergyIndex = Int_t(enu_bins[i]*MatrixBins/FinalVisibleEnergy); TrueEnergyIndex < Int_t(enu_bins[i+1]*MatrixBins/FinalVisibleEnergy); TrueEnergyIndex++)
        {
            std::cout << "enu_bins[i]*MatrixBins/FinalVisibleEnergy" << enu_bins[i]*MatrixBins/FinalVisibleEnergy << std::endl;
            std::cout << "enu_bins[i+1]*MatrixBins/FinalVisibleEnergy" <<enu_bins[i+1]*MatrixBins/FinalVisibleEnergy<< std::endl;
            std::cout <<TrueEnergyIndex << std::endl;

            OscDeltaVisibleSpectrumSumH[i]->Add(OscDeltaVisibleSpectrumH[TrueEnergyIndex]);
        }
        
        PredictionVisH[i]=(TH1F*)OscDeltaVisibleSpectrumSumH[i]->Rebin(n_evis_bins,Form("Rebinned Spectrum Vis for True Energy Index %d",i),evis_bins);

        for (Int_t j = 0; j < n_evis_bins; j++)
        {
            MatrixH->SetBinContent(i+1,j+1, PredictionVisH[i]->GetBinContent(j+1));
        }
    }
    
    if(NormalizeRow)
    {
        for(Int_t j=0;j<n_evis_bins;j++)
        {
            Norma[j]=0;
            
            for(Int_t i=0;i<n_etrue_bins;i++)
            {
                Norma[j] = Norma[j]+MatrixH->GetBinContent(i+1,j+1);
            }
        }
        
        for (Int_t i = 0; i < n_etrue_bins; i++)
        {
            for (Int_t j = 0; j < n_evis_bins; j++)
            {
                if(Norma[j]!=0)
                {
                    MatrixH->SetBinContent(i+1,j+1,MatrixH->GetBinContent(i+1,j+1)/Norma[j]);//Normalization so Σi E(i,j) = 1; (Σ(x axis) =1)
                }
            }
        }
    }
    else
    {
        for(Int_t j=0;j<n_etrue_bins;j++)
        {
            Norma[j]=0;
            
            for(Int_t i=0;i<n_evis_bins;i++)
            {
                Norma[j] = Norma[j]+MatrixH->GetBinContent(j+1,i+1);
            }
        }
        
        for (Int_t i = 0; i < n_evis_bins; i++)
        {
            for (Int_t j = 0; j < n_etrue_bins; j++)
            {
                if(Norma[j]!=0)
                {
                    MatrixH->SetBinContent(j+1,i+1,MatrixH->GetBinContent(j+1,i+1)/Norma[j]);//Normalization so Σj E(i,j) = 1; (Σ(y axis) =1)
                }
            }
        }
    }

    //Code to check correct normalization:
    Double_t Normx[MatrixBins]={};
    Double_t Normy[MatrixBins]={};
    Double_t FineNormx[MatrixBins]={};
    Double_t FineNormy[MatrixBins]={};
    
    for(Int_t i=0;i<n_etrue_bins;i++)
    {
        for(Int_t j=0;j<n_evis_bins;j++)
        {
            //Check normalization
            Normy[i]=Normy[i]+MatrixH->GetBinContent(i+1,j+1);
        }
        
    }
    for(Int_t j=0;j<n_evis_bins;j++)
    {
        for(Int_t i=0;i<n_etrue_bins;i++)
        {
            Normx[j]=Normx[j]+MatrixH->GetBinContent(i+1,j+1);
        }
    }
    
    for(Int_t j=0;j<n_evis_bins;j++)
    {
        TH1D *px = MatrixH->ProjectionX("x",j,j);
        std::cout << "Norma in X" << Normx[j]<<std::endl;
        std::cout << "Norma in X " << px->Integral() <<std::endl;;
    }
    for(Int_t i=0;i<n_etrue_bins;i++)
    {
        TH1D *py = MatrixH->ProjectionY("y",i,i);
        std::cout << "Norma in Y" << Normy[i]<<std::endl;
        std::cout << "Norma in Y " << py->Integral() <<std::endl;;
    }
    for(Int_t i=0;i<MatrixBins;i++)
    {
        for(Int_t j=0;j<MatrixBins;j++)
        {
            //Check normalization
            FineNormy[i]=FineNormy[i]+FineMatrixH->GetBinContent(i+1,j+1);
        }
        
    }
    for(Int_t j=0;j<MatrixBins;j++)
    {
        for(Int_t i=0;i<MatrixBins;i++)
        {
            FineNormx[j]=FineNormx[j]+FineMatrixH->GetBinContent(i+1,j+1);
        }
    }
    
    for(Int_t j=0;j<MatrixBins;j++)
    {
        TH1D *Finepx = FineMatrixH->ProjectionX("x",j+1,j+1);
        std::cout << "Fine Norma in X" << FineNormx[j]<<std::endl;
        std::cout << "Fine Norma in X " << Finepx->Integral() <<std::endl;;
    }
    for(Int_t i=0;i<MatrixBins;i++)
    {
        TH1D *Finepy = FineMatrixH->ProjectionY("y",i+1,i+1);
        std::cout << "Fine Norma in Y" << FineNormy[i]<<std::endl;
        std::cout << "Fine Norma in Y " << Finepy->Integral() <<std::endl;;
    }
    
    
    
    // Here multiply Matrix and flat spectrum, both of them rebinned.
   
    VisibleHisto = new TH1F(Form("VisibleHisto Spectrum"),Form("VisibleHisto Spectrum"),n_evis_bins,evis_bins);
    FineVisibleHisto = new TH1F(Form("FineVisibleHisto Spectrum"),Form("FineVisibleHisto Spectrum"),MatrixBins,0,FinalVisibleEnergy);

    for(Int_t i=1; i<=n_evis_bins; i++)
    {
        for(Int_t j=1; j<=n_etrue_bins; j++)
        {
            VisibleHisto->SetBinContent(i,VisibleHisto->GetBinContent(i) + MatrixH->GetBinContent(j,i) * RebinnedPredictionTrueH->GetBinContent(j));
        }
    }
    
    // Here fine matrix, and then rebin visiblehisto spectrum
    
    for(Int_t i=1; i<=MatrixBins; i++)
    {
        for(Int_t j=1; j<=MatrixBins; j++)
        {
            FineVisibleHisto->SetBinContent(i,FineVisibleHisto->GetBinContent(i) + FineMatrixH->GetBinContent(i,j) * FinePredictionTrueH->GetBinContent(j));
        }
    }
    
    RebinnedFineVisibleHisto=(TH1F*)FineVisibleHisto->Rebin(n_evis_bins,Form("Rebinned Visible Histo"),evis_bins);

    // Check both are the same, if they differ find out why.
    // Plot
    TCanvas* c1 = new TCanvas("a","a",500,500);
    c1->Divide(2,2);
    
    c1->cd(1);
    
    PredictionVisH[0]->Draw();
    c1->cd(2);
    
    PredictionVisH[1]->Draw();
    c1->cd(3);
    OscDeltaVisibleSpectrumSumH[0]->Draw();
    c1->cd(4);
    OscDeltaVisibleSpectrumSumH[1]->Draw();
    
    TCanvas* c2 = new TCanvas("Matrix","Matrix");

    MatrixH->Draw("colz");
    
    TCanvas* c3 = new TCanvas("RebinnedResults","RebinnedResults");
    c3->Divide(2,2);

    c3->cd(1);
    PredictionTrueH->Draw();
    
    c3->cd(2);
    RebinnedPredictionTrueH->Draw();
    
    c3->cd(3);
    VisibleHisto->Draw();
    
    TCanvas* c4 = new TCanvas("FineMatrix","FineMatrix");
    
    FineMatrixH->Draw("colz");
    
    TCanvas* c5 = new TCanvas("FineResults","FineResults");
    c5->Divide(2,2);
    
    c5->cd(1);
    PredictionTrueH->Draw();
    
    c5->cd(2);
    FineVisibleHisto->Draw();
    
    c5->cd(3);
    RebinnedFineVisibleHisto->Draw();

}