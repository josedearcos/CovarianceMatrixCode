#pragma once
#include "TH2F.h"
#include "TFile.h"
#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2.h"
#include "TMath.h"
#include <math.h>
#include <TArrayD.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TCanvas.h"    
#include "NominalData.h"
#include "CrossSection.h"
#include "CovarianceMatrix3.h"
#include "Prediction.h"

const Int_t n_evis_bins =37;
const Int_t nsteps = 101;

class Fitter
{
private:

    NominalData* Nom;

    Prediction* Pred;
    TRandom3* rand;
    
    bool StatisticalFluctuation;
    bool SimpleReactorModel;
    //AD configuration parameters:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    Int_t x,y;

    Int_t Combine;
    Int_t MaxNear;
    Int_t MaxFar;
    Int_t MaxBins;
    Int_t Nweeks;
    Int_t hall;
    
    Double_t chi2;
    Double_t Sin22t13;
    Double_t sin22t13[nsteps];
    
    //  Background rates:
    Double_t ScaleAcc[MaxDetectors][MaxPeriods];
    Double_t ScaleLiHe[MaxDetectors][MaxPeriods];
    Double_t ScaleFN[MaxDetectors][MaxPeriods];
    Double_t ScaleAmC[MaxDetectors][MaxPeriods];
    TH1F* BackgroundSpectrumH[MaxDetectors][MaxPeriods];
    TH1F* NearBackgroundSpectrumH[MaxNearDetectors][MaxPeriods];
    TH1F* FarBackgroundSpectrumH[MaxFarDetectors][MaxPeriods];
    TH1F* RandomBackgroundSpectrumH[MaxDetectors][MaxPeriods];
    TH1F* NearRandomBackgroundSpectrumH[MaxNearDetectors][MaxPeriods];
    TH1F* FarRandomBackgroundSpectrumH[MaxFarDetectors][MaxPeriods];
    TH1F* AccidentalsH[MaxDetectors][MaxPeriods];
    TH1F* LiHeH[MaxDetectors][MaxPeriods];
    TH1F* FastNeutronsH[MaxDetectors][MaxPeriods];
    TH1F* AmCH[MaxDetectors][MaxPeriods];
    
    // Data and prediction histograms
    TH1F* PredictionVisH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];//  Prediction with no alterations, nominal backgrounds may be added. Vis Binning (after detector response is applied)
    TH1F* ADSpectrumVisH[MaxDetectors][MaxPeriods];
    
    TH1F* FarDataH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];//  Prediction with no alterations, nominal backgrounds may be added. Vis Binning (after detector response is applied)
    TH1F* NearDataH[MaxDetectors][MaxPeriods];
    
    //  Root Histograms
    TH2F* DAmCCovarianceMatrixH[MaxPeriods];
    TH2F* VAmCCovarianceMatrixH[MaxPeriods];
    TH2F* DFNCovarianceMatrixH[MaxPeriods];
    TH2F* VFNCovarianceMatrixH[MaxPeriods];
    TH2F* DLiHeCovarianceMatrixH[MaxPeriods];
    TH2F* VLiHeCovarianceMatrixH[MaxPeriods];
    TH2F* VAccCovarianceMatrixH[MaxPeriods];
    TH2F* IsotopeCovarianceMatrixH[MaxPeriods];
    TH2F* ReactorPowerCovarianceMatrixH[MaxPeriods];
    TH2F* RelativeEnergyCovarianceMatrixH[MaxPeriods];
    TH2F* IAVCovarianceMatrixH[MaxPeriods];
    TH2F* NLCovarianceMatrixH[MaxPeriods];
    TH2F* ResolutionCovarianceMatrixH[MaxPeriods];
    TH2F* StatisticalCovarianceMatrixH[MaxPeriods];
    TH2F* BackgroundsCovarianceMatrixH[MaxPeriods];
    TH2F* SystematicCovarianceMatrixH[MaxPeriods];
    TH1F* ChiSquareH[MaxPeriods];
    TH2F* TotalCovarianceMatrixH[MaxPeriods];
    TH2F* InvTotalCovarianceMatrixH[MaxPeriods];
    TH2F* StatCov2H[MaxPeriods];
    TH2F* UnityH[MaxPeriods];
    
    //  Matrices
    Double_t DAmCCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t VAmCCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t DFNCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t VFNCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t DLiHeCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t VLiHeCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t VAccCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t IsotopeCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t ReactorPowerCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t RelativeEnergyCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t IAVCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t NLCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t ResolutionCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t StatisticalCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t BackgroundsCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t SystematicCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t InvTotalCovarianceMatrixM[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t CovStat[9*MaxNbins][9*MaxNbins][MaxPeriods];

    Double_t UnityM[9*MaxNbins][9*MaxNbins][MaxPeriods];

    TArrayD TotalCovarianceMatrixM;

    Double_t Sigma_Near[MaxFarDetectors][MaxNearDetectors][MaxPeriods][MaxNbins];
    Double_t Sigma_Far[MaxFarDetectors][MaxNearDetectors][MaxPeriods][MaxNbins];
    
    void GenerateStatisticalCovarianceMatrix(Int_t);
    void SaveStatisticalCovarianceMatrix(Int_t);
    void NormCov(TH2F*);
    void LoadBackgrounds();
//    void LoadRootCovarianceMatrices();
    void LoadTxtCovarianceMatrices();
    void CombineMatrices(Int_t);
    void InvertMatrix();
    void LoadData();
    void LoadPrediction();
    void SaveChiSquare(Int_t);
    void SaveCovarianceMatrices(Int_t);
    void ScanParameters();
    void ChiSquare();
    
    void ApplyStatisticalFluctuation(TH1F*);

public:
    Fitter();
    Fitter(NominalData*);
    void MainFitter();
};

Fitter :: Fitter()
{
    Nom = new NominalData();
    Pred = new Prediction(Nom);
    rand = new TRandom3();
    
    Combine = Nom->GetCombineMode();
    Nweeks = Nom->GetWeeks();
    
    SimpleReactorModel = Nom->GetSimpleReactorModel();
    StatisticalFluctuation = Nom->GetStatisticalFluctuation();
    Sin22t13 = Nom->GetSin22t13();
    NADs = Nom->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    if(NADs == 8)//    ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }

    if (Combine == 1)
    {
        MaxNear=1;
        MaxFar=1;
        MaxBins=n_evis_bins;
    }
    else if(Combine == 2)
    {
        MaxNear=2;
        MaxFar=1;
        MaxBins=2*n_evis_bins;
    }
    else
    {
        MaxNear = ADsEH1+ADsEH2;
        MaxFar = ADsEH3;
        MaxBins=9*n_evis_bins;
    }
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            ScaleAcc[AD][week]=Nom->GetAccidentalEvents(AD,week);
            ScaleLiHe[AD][week]=Nom->GetLiHeEvents(hall,week);
            ScaleFN[AD][week]=Nom->GetFNEvents(hall,week);
            ScaleAmC[AD][week]=Nom->GetAmCEvents(hall,week);
        }
    }
}

Fitter :: Fitter(NominalData* Data)
{
    Pred = new Prediction(Data);
    rand = new TRandom3();

    Combine = Data->GetCombineMode();
    Nweeks = Data->GetWeeks();
    
    SimpleReactorModel = Data->GetSimpleReactorModel();
    StatisticalFluctuation = Data->GetStatisticalFluctuation();
    Sin22t13 = Data->GetSin22t13();
    NADs = Data->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    if(NADs == 8)//    ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }

    if (Combine == 1)
    {
        MaxNear=1;
        MaxFar=1;
        MaxBins=n_evis_bins;
    }
    else if(Combine == 2)
    {
        MaxNear=2;
        MaxFar=1;
        MaxBins=2*n_evis_bins;
    }
    else
    {
        MaxNear = ADsEH1+ADsEH2;
        MaxFar = ADsEH3;
        MaxBins=9*n_evis_bins;
    }
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            ScaleAcc[AD][week]=Data->GetAccidentalEvents(AD,week);
            ScaleLiHe[AD][week]=Data->GetLiHeEvents(hall,week);
            ScaleFN[AD][week]=Data->GetFNEvents(hall,week);
            ScaleAmC[AD][week]=Data->GetAmCEvents(hall,week);
        }
    }
}

void Fitter :: MainFitter()
{
//    LoadRootCovarianceMatrices();
    LoadTxtCovarianceMatrices();
    LoadData();

    ScanParameters();
    //for theta13 Pred->Set(theta13)
    //Multiply SystematicCovVariance times Fi,Fj
    //Sum Final Cov Matrix and Invert
    //Use same Fi,Fj to Fit the Data
    //Keep minimum chi2 value
    
    //Repeat for DeltaM
}

void Fitter :: ScanParameters()
{
    //  Reset file.
    TFile* ResetStatisticalCovMatrixF = TFile::Open("./CovarianceMatrices/StatisticalCovarianceMatrix.root","recreate");
    ResetStatisticalCovMatrixF->Close();
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        const Double_t s22t13start=0.05;
        const Double_t s22t13end=0.15;
        Double_t chi2_min = 1e10;
        Double_t s2t_min = 1e10;
        Int_t degenerancy;
        
         ChiSquareH[week] = new TH1F("ChiSquare Distribution","ChiSquare Distribution",nsteps,s22t13start,s22t13end);

        for(Int_t step=0;step<nsteps;++step)
        {
            sin22t13[step]=(s22t13end-s22t13start)*step*1./(nsteps-1)+s22t13start;
            std::cout << "sin22t13 " << sin22t13[step] << std::endl;
            
            Pred->MakePrediction(sin22t13[step]);
            LoadPrediction();
            InvertMatrix();
            ChiSquare();
            ChiSquareH[week]->SetBinContent(step+1,TMath::Abs(chi2));

            if (TMath::Abs(chi2)<chi2_min)
            {
                s2t_min=sin22t13[step];
                chi2_min=TMath::Abs(chi2);
            }
            else if(TMath::Abs(chi2)==chi2_min)
            {
                degenerancy++;
            }
        }
        
        std::cout << "Chisquare min" << chi2_min << std::endl;
        std::cout << "Sin22t13 min" << s2t_min << std::endl;
        std::cout << "ChiMin Square degenerated " << degenerancy << " times" << std::endl;
        
    SaveChiSquare(week);
    delete ChiSquareH[week];
    }
}

void Fitter :: ChiSquare()
{
    chi2=0;
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        x =0;
        y =0;
        
        for (Int_t neari=0; neari<MaxNear; neari++)
        {
            Int_t Ni1,Ni2,Ni3,Ni4;
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Ni1=1; Ni2=0; Ni3=0; Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<MaxFar; fari++)
            {
                Int_t Fi1,Fi2,Fi3,Fi4;
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Fi1=fari+1;Fi2 = 0;Fi3 = 0;Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFar; Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFar; Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=MaxFar; Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<MaxNear; nearj++)
                {
                    Int_t Nj1,Nj2,Nj3,Nj4;
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<MaxFar; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        Int_t Fj1,Fj2,Fj3,Fj4;
                        if(Nj1!=Nj2){Fj1=farj+1;Fj2 = 0;Fj3 = 0;Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFar; Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFar; Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=MaxFar; Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;

                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                //                chi2+= ((FarDataH[0][0][0]->GetBinContent(i+1)-PredictionVisH[0][0][0]->GetBinContent(i+1))*InvTotalCovarianceMatrixH[week]->GetBinContent(i+1,j+1)*(FarDataH[0][0][0]->GetBinContent(j+1)-PredictionVisH[0][0][0]->GetBinContent(j+1)));
                                chi2+= ((FarDataH[fari][neari][week]->GetBinContent(i+1)-PredictionVisH[fari][neari][week]->GetBinContent(i+1))*InvTotalCovarianceMatrixM[x][y][week]*(FarDataH[farj][nearj][week]->GetBinContent(j+1)-PredictionVisH[farj][nearj][week]->GetBinContent(j+1)));
                                
                                
                                std::cout << "InvTotalCovarianceMatrixH " << InvTotalCovarianceMatrixM[x][y][week]<< std::endl;
                                std::cout << "FarDataHi " << FarDataH[fari][neari][week]->GetBinContent(i+1) << std::endl;
                                std::cout << "PredictionVisH i" << PredictionVisH[fari][neari][week]->GetBinContent(i+1) << std::endl;
                                std::cout << "FarDataH - PredictionVisH i " <<(FarDataH[fari][neari][week]->GetBinContent(i+1)-PredictionVisH[fari][neari][week]->GetBinContent(i+1)) << std::endl;
                                std::cout << "FarDataHj " << FarDataH[farj][nearj][week]->GetBinContent(j+1) << std::endl;
                                std::cout << "PredictionVisHj " << PredictionVisH[farj][nearj][week]->GetBinContent(j+1) << std::endl;
                                std::cout << "FarDataH - PredictionVisH j " <<(FarDataH[farj][nearj][week]->GetBinContent(j+1)-PredictionVisH[farj][nearj][week]->GetBinContent(j+1))<< std::endl;
                                std::cout << "portion of chi2 " <<((FarDataH[fari][neari][week]->GetBinContent(i+1)-PredictionVisH[fari][neari][week]->GetBinContent(i+1))*InvTotalCovarianceMatrixM[x][y][week]*(FarDataH[farj][nearj][week]->GetBinContent(j+1)-PredictionVisH[farj][nearj][week]->GetBinContent(j+1))) << std::endl;
                                
                            }
                        }
                    }
                }
            }
        }
        std::cout << " Total Chisquare" << chi2 << std::endl;
    }
}

void Fitter :: LoadData()
{
    if(SimpleReactorModel)
    {
        Pred->MakePrediction(Sin22t13);

        // Let's use the nominal data for a fixed real sin22t13
        TFile* DataF = TFile::Open("./RootOutputs/PredictedSpectrum.root");
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for (Int_t near = 0; near<MaxNear; near++)
            {
                NearDataH[near][week]=(TH1F*)gDirectory->Get(Form("Near Prediction AD%d",near));
            }
            for (Int_t near =0; near<MaxNear; near++)
            {
                for (Int_t far =0; far<MaxFar; far++)
                {
                    FarDataH[far][near][week]=(TH1F*)gDirectory->Get(Form("Far Prediction AD%d from AD%d",far,near));//
                    if (StatisticalFluctuation)
                    {
                        ApplyStatisticalFluctuation(FarDataH[far][near][week]);
                    }
                }
            }
        }
        DataF->Close();
    }
}

void Fitter :: LoadPrediction()
{
    TFile* PredictionF = TFile::Open("./RootOutputs/PredictedSpectrum.root");
    for (Int_t week = 0; week<Nweeks; week++)
    {
//        for (Int_t near = 0; near<MaxNear; near++)
//        {
//            ADSpectrumVisH[near][week]=(TH1F*)gDirectory->Get(Form("Near Prediction AD%d",near));
//        }
        for (Int_t near =0; near<MaxNear; near++)
        {
            for (Int_t far =0; far<MaxFar; far++)
            {
                PredictionVisH[far][near][week]=(TH1F*)gDirectory->Get(Form("Far Prediction AD%d from AD%d",far,near));
            }
        }
    }
    PredictionF->Close();
}
//
//void Fitter :: LoadRootCovarianceMatrices()
//{
//    for (Int_t week = 0; week<Nweeks; week++)
//    {
//        TFile* VAccCovarianceMatrixF = new TFile("./CovarianceMatrices/VaryAccidentalCovarianceMatrix.root");
//        VAccCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("Vary Accidental Covariance Matrix%d",week));
//        VAccCovarianceMatrixF->Close();
//        
//        TFile* VFNCovarianceMatrixF = new TFile("./CovarianceMatrices/VaryFastNeutronsCovarianceMatrix.root");
//        VFNCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("Vary FN Covariance Matrix%d",week));
//        VFNCovarianceMatrixF->Close();
//        
//        TFile* VLiHeCovarianceMatrixF = new TFile("./CovarianceMatrices/VaryLiHeCovarianceMatrix.root");
//        VLiHeCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("Vary LiHe Covariance Matrix%d",week));
//        VLiHeCovarianceMatrixF->Close();
//        
//        TFile* VAmCCovarianceMatrixF = new TFile("./CovarianceMatrices/VaryAmCCovarianceMatrix.root");
//        VAmCCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("Vary AmC Covariance Matrix%d",week));
//        VAmCCovarianceMatrixF->Close();
//        
//        TFile* DFNCovarianceMatrixF = new TFile("./CovarianceMatrices/DistortFastNeutronsCovarianceMatrix.root");
//        DFNCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("Distort FN Covariance Matrix%d",week));
//        DFNCovarianceMatrixF->Close();
//        
//        TFile* DLiHeCovarianceMatrixF = new TFile("./CovarianceMatrices/DistortLiHeCovarianceMatrix.root");
//        DLiHeCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("Distort LiHe Covariance Matrix%d",week));
//        DLiHeCovarianceMatrixF->Close();
//        
//        TFile* DAmCCovarianceMatrixF = new TFile("./CovarianceMatrices/DistortAmCCovarianceMatrix.root");
//        DAmCCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("Distort AmC Covariance Matrix%d",week));
//        DAmCCovarianceMatrixF->Close();
//        
//        TFile* IsotopeCovarianceMatrixF = new TFile("./CovarianceMatrices/IsotopeCovarianceMatrix.root");
//        IsotopeCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("Isotope Covariance Matrix%d",week));
//        IsotopeCovarianceMatrixF->Close();
//        
//        TFile* ReactorPowerCovarianceMatrixF = new TFile("./CovarianceMatrices/ReactorPowerCovarianceMatrix.root");
//        ReactorPowerCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("Reactor Power Covariance Matrix%d",week));
//        ReactorPowerCovarianceMatrixF->Close();
//        //    TFile* AbsoluteEnergyCovarianceMatrixF = new TFile("./CovarianceMatrices/AbsoluteEnergyCovarianceMatrix.root");    //    AbsoluteEnergyCovarianceMatrixF->Close();
//        //    TFile* RelativeEnergyCovarianceMatrixF = new TFile("./CovarianceMatrices/RelativeEnergyCovarianceMatrix.root");
//        //    RelativeEnergyCovarianceMatrixF->Close();
//        
//        TFile* IAVCovarianceMatrixF = new TFile("./CovarianceMatrices/IAVCovarianceMatrix.root");
//        IAVCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("IAV Covariance Matrix%d",week));
//        IAVCovarianceMatrixF->Close();
//        
//        TFile* NLCovarianceMatrixF = new TFile("./CovarianceMatrices/NLCovarianceMatrix.root");
//        NLCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("NL Covariance Matrix%d",week));
//        NLCovarianceMatrixF->Close();
//        
//        TFile* ResolutionCovarianceMatrixF = new TFile("./CovarianceMatrices/ResolutionCovarianceMatrix.root");
//        ResolutionCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("Resolution Covariance Matrix%d",week));
//        ResolutionCovarianceMatrixF->Close();
//        
////        TFile* StatisticalCovarianceMatrixF = new TFile("./CovarianceMatrices/StatisticalCovarianceMatrix.root");
////        StatisticalCovarianceMatrixH[week]=(TH2F*)gDirectory->Get(Form("Statistical Covariance Matrix%d",week));
////        StatisticalCovarianceMatrixF->Close();
//        
//        //    TFile* Efficiency CovarianceMatrixF = new TFile("./CovarianceMatrices/EfficiencyCovarianceMatrix.root");
//        //    EfficiencyCovarianceMatrixF->Close();
//    }
//}


void Fitter :: LoadTxtCovarianceMatrices()
{
    Char_t filenameCov[100];

    for (Int_t week = 0; week<Nweeks; week++)
    {
        sprintf(filenameCov,"./CovarianceMatrices/VaryAccidentalCovarianceMatrix.txt");

        ifstream covfile_vacc(filenameCov);
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j <MaxBins; j++)
            {
                covfile_vacc >>  VAccCovarianceMatrixM[i][j][week];
            }
        }
        
        sprintf(filenameCov,"./CovarianceMatrices/VaryFastNeutronsCovarianceMatrix.txt");
        ifstream covfile_vfn(filenameCov);
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j <MaxBins; j++)
            {
                covfile_vfn >>  VFNCovarianceMatrixM[i][j][week];
            }
        }
        
        sprintf(filenameCov,"./CovarianceMatrices/VaryLiHeCovarianceMatrix.txt");
        ifstream covfile_vlihe(filenameCov);
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j <MaxBins; j++)
            {
                covfile_vlihe >>  VLiHeCovarianceMatrixM[i][j][week];
            }
        }
        sprintf(filenameCov,"./CovarianceMatrices/VaryAmCCovarianceMatrix.txt");
        ifstream covfile_vamc(filenameCov);
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j <MaxBins; j++)
            {
                covfile_vamc >>  VAmCCovarianceMatrixM[i][j][week];
            }
        }
        
        sprintf(filenameCov,"./CovarianceMatrices/DistortFastNeutronsCovarianceMatrix.txt");
        ifstream covfile_dfn(filenameCov);
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j <MaxBins; j++)
            {
                covfile_dfn >>  DFNCovarianceMatrixM[i][j][week];
            }
        }
        
        sprintf(filenameCov,"./CovarianceMatrices/DistortLiHeCovarianceMatrix.txt");
        ifstream covfile_dlihe(filenameCov);
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j <MaxBins; j++)
            {
                covfile_dlihe >>  DLiHeCovarianceMatrixM[i][j][week];
            }
        }
        
        sprintf(filenameCov,"./CovarianceMatrices/DistortAmCCovarianceMatrix.txt");
        ifstream covfile_damc(filenameCov);
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j <MaxBins; j++)
            {
                covfile_damc >>  DAmCCovarianceMatrixM[i][j][week];
            }
        }
        
        sprintf(filenameCov,"./CovarianceMatrices/IsotopeCovarianceMatrix.txt");
        ifstream covfile_iso(filenameCov);
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j <MaxBins; j++)
            {
                covfile_iso >>  IsotopeCovarianceMatrixM[i][j][week];
            }
        }
        
        sprintf(filenameCov,"./CovarianceMatrices/ReactorPowerCovarianceMatrix.txt");
        ifstream covfile_pow(filenameCov);
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j <MaxBins; j++)
            {
                covfile_pow >> ReactorPowerCovarianceMatrixM[i][j][week];
            }
        }
        
        sprintf(filenameCov,"./CovarianceMatrices/IAVCovarianceMatrix.txt");
        ifstream covfile_iav(filenameCov);
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j <MaxBins; j++)
            {
                covfile_iav >> IAVCovarianceMatrixM[i][j][week];
            }
        }
        
        sprintf(filenameCov,"./CovarianceMatrices/NLCovarianceMatrix.txt");
        ifstream covfile_nl(filenameCov);
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j <MaxBins; j++)
            {
                covfile_nl >> NLCovarianceMatrixM[i][j][week];
            }
        }
        
        sprintf(filenameCov,"./CovarianceMatrices/ResolutionCovarianceMatrix.txt");
        ifstream covfile_reso(filenameCov);
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j <MaxBins; j++)
            {
                covfile_reso >> ResolutionCovarianceMatrixM[i][j][week];
            }
        }
//        sprintf(filenameCov,"./CovarianceMatrices/AbsoluteEnergyCovarianceMatrix.txt");
//        ifstream covfile_abs(filenameCov);
//        for (Int_t i = 0; i < MaxBins; i++)
//        {
//            for (Int_t j = 0; j <MaxBins; j++)
//            {
//                covfile_abs >>  AbsoluteEnergyCovarianceMatrixM[i][j][week];
//            }
//        }
//        sprintf(filenameCov,"./CovarianceMatrices/RelativeEnergyCovarianceMatrix.txt");
//        ifstream covfile_rel(filenameCov);
//        for (Int_t i = 0; i < MaxBins; i++)
//        {
//            for (Int_t j = 0; j <MaxBins; j++)
//            {
//                covfile_rel >> RelativeEnergyCovarianceMatrixM[i][j][week];
//            }
//        }
//        sprintf(filenameCov,"./CovarianceMatrices/EfficiencyCovarianceMatrix.txt");
//        ifstream covfile_eff(filenameCov);
//        for (Int_t i = 0; i < MaxBins; i++)
//        {
//            for (Int_t j = 0; j <MaxBins; j++)
//            {
//                covfile_eff >> EfficiencyCovarianceMatrixM[i][j][week];
//            }
//        }
    }
}

void Fitter :: SaveCovarianceMatrices(Int_t week)
{
        Int_t x =0;
        Int_t y =0;
        
        for (Int_t neari=0; neari<MaxNear; neari++)
        {
            Int_t Ni1,Ni2,Ni3,Ni4;
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Ni1=1; Ni2=0; Ni3=0; Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<MaxFar; fari++)
            {
                Int_t Fi1,Fi2,Fi3,Fi4;
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Fi1=fari+1;Fi2 = 0;Fi3 = 0;Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFar; Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFar; Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=MaxFar; Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<MaxNear; nearj++)
                {
                    Int_t Nj1,Nj2,Nj3,Nj4;
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<MaxFar; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        Int_t Fj1,Fj2,Fj3,Fj4;
                        if(Nj1!=Nj2){Fj1=farj+1;Fj2 = 0;Fj3 = 0;Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFar; Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFar; Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=MaxFar; Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;

                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                std::cout << "INDEX Y" << y << std::endl;
                                IsotopeCovarianceMatrixH[week]->SetBinContent(x+1,y+1,IsotopeCovarianceMatrixM[x][y][week]);
                                ReactorPowerCovarianceMatrixH[week]->SetBinContent(x+1,y+1,ReactorPowerCovarianceMatrixM[x][y][week]);
                                IAVCovarianceMatrixH[week]->SetBinContent(x+1,y+1,IAVCovarianceMatrixM[x][y][week]);
                                NLCovarianceMatrixH[week]->SetBinContent(x+1,y+1,NLCovarianceMatrixM[x][y][week]);
                                ResolutionCovarianceMatrixH[week]->SetBinContent(x+1,y+1,ResolutionCovarianceMatrixM[x][y][week]);

                                VAccCovarianceMatrixH[week]->SetBinContent(x+1,y+1,VAccCovarianceMatrixM[x][y][week]);
                                VLiHeCovarianceMatrixH[week]->SetBinContent(x+1,y+1,VLiHeCovarianceMatrixM[x][y][week]);
                                VFNCovarianceMatrixH[week]->SetBinContent(x+1,y+1,VFNCovarianceMatrixM[x][y][week]);
                                VAmCCovarianceMatrixH[week]->SetBinContent(x+1,y+1,VAmCCovarianceMatrixM[x][y][week]);
                                DLiHeCovarianceMatrixH[week]->SetBinContent(x+1,y+1,DLiHeCovarianceMatrixM[x][y][week]);
                                DFNCovarianceMatrixH[week]->SetBinContent(x+1,y+1,DFNCovarianceMatrixM[x][y][week]);
                                DAmCCovarianceMatrixH[week]->SetBinContent(x+1,y+1,DAmCCovarianceMatrixM[x][y][week]);
                                SystematicCovarianceMatrixH[week]->SetBinContent(x+1,y+1,SystematicCovarianceMatrixM[x][y][week]);
                                BackgroundsCovarianceMatrixH[week]->SetBinContent(x+1,y+1,BackgroundsCovarianceMatrixM[x][y][week]);
                                StatisticalCovarianceMatrixH[week]->SetBinContent(x+1,y+1,CovStat[x][y][week]);
                                TotalCovarianceMatrixH[week]->SetBinContent(x+1,y+1,TotalCovarianceMatrixM[y+x*MaxBins]);
                                InvTotalCovarianceMatrixH[week]->SetBinContent(x+1,y+1,InvTotalCovarianceMatrixM[x][y][week]);
                            }
                        }
                    }
                }
            }
        }
    
    //Check unity matrix:
    UnityH[week]=(TH2F*)TotalCovarianceMatrixH[week]->Clone("Unity");
    UnityH[week]->Reset();
    
    for(Int_t i=0; i<MaxBins; i++)
    {
        for(Int_t j=0; j<MaxBins; j++)
        {
            UnityM[i][j][week] = 0;
            for(Int_t k=0; k<MaxBins; k++)
            {
                UnityM[i][j][week] += TotalCovarianceMatrixM[i+MaxBins*k] * InvTotalCovarianceMatrixM[k][j][week];
            }
            UnityH[week]->SetBinContent(i+1,j+1,UnityM[i][j][week]);
        }
    }

        TFile* SaveCovarianceMatricesF = new TFile(Form("./CovarianceMatrices/FitterCovarianceMatrixResultsPeriod%d.root",week),"recreate");
        VAccCovarianceMatrixH[week]->Write();
        VFNCovarianceMatrixH[week]->Write();
        VLiHeCovarianceMatrixH[week]->Write();
        VAmCCovarianceMatrixH[week]->Write();
        DFNCovarianceMatrixH[week]->Write();
        DLiHeCovarianceMatrixH[week]->Write();
        DAmCCovarianceMatrixH[week]->Write();
        
        IsotopeCovarianceMatrixH[week]->Write("Isotope Matrix");
        ReactorPowerCovarianceMatrixH[week]->Write("Power Matrix");
        IAVCovarianceMatrixH[week]->Write("IAV Matrix");
        NLCovarianceMatrixH[week]->Write("NL Matrix");
        ResolutionCovarianceMatrixH[week]->Write("Reso Matrix");
        
        BackgroundsCovarianceMatrixH[week]->Write("Background Covariance Matrix");
        SystematicCovarianceMatrixH[week]->Write("Systematic Covariance Matrix");
        TotalCovarianceMatrixH[week]->Write("Total Covariance Matrix");
        InvTotalCovarianceMatrixH[week]->Write("Inv Matrix");
        UnityH[week]->Write("Unity Matrix");
        SaveCovarianceMatricesF->Close();
}

void Fitter :: CombineMatrices(Int_t week)
{
    BackgroundsCovarianceMatrixH[week]= new TH2F("Background Covariance Matrix","Background Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    SystematicCovarianceMatrixH[week]=(TH2F*)BackgroundsCovarianceMatrixH[week]->Clone("Systematic Covariance Matrix");
    SystematicCovarianceMatrixH[week]->SetTitle("Systematic Covariance Matrix");
    StatisticalCovarianceMatrixH[week]=(TH2F*)BackgroundsCovarianceMatrixH[week]->Clone("Statistical Covariance Matrix");
    StatisticalCovarianceMatrixH[week]->SetTitle("Statistical Covariance Matrix");
    TotalCovarianceMatrixH[week]=(TH2F*)BackgroundsCovarianceMatrixH[week]->Clone("Total Covariance Matrix");
    TotalCovarianceMatrixH[week]->SetTitle("Total Covariance Matrix");
    InvTotalCovarianceMatrixH[week]=(TH2F*)TotalCovarianceMatrixH[week]->Clone("Inv Total Covariance Matrix");
    InvTotalCovarianceMatrixH[week]->SetTitle("Inv Total Covariance Matrix");
    
    VAccCovarianceMatrixH[week]=(TH2F*)BackgroundsCovarianceMatrixH[week]->Clone("Vary Acc Covariance Matrix");
    VAccCovarianceMatrixH[week]->SetTitle("Vary Accidentals Covariance Matrix");
    VFNCovarianceMatrixH[week]=(TH2F*)BackgroundsCovarianceMatrixH[week]->Clone("Vary FN Covariance Matrix");
    VFNCovarianceMatrixH[week]->SetTitle("Vary FN Covariance Matrix");
    VLiHeCovarianceMatrixH[week]=(TH2F*)BackgroundsCovarianceMatrixH[week]->Clone("Vary LiHe Covariance Matrix");
    VLiHeCovarianceMatrixH[week]->SetTitle("Vary LiHe Covariance Matrix");
    VAmCCovarianceMatrixH[week]=(TH2F*)BackgroundsCovarianceMatrixH[week]->Clone("Vary AmC Covariance Matrix");
    VAmCCovarianceMatrixH[week]->SetTitle("Vary AmC Covariance Matrix");

    DFNCovarianceMatrixH[week]=(TH2F*)BackgroundsCovarianceMatrixH[week]->Clone("Distort FN Covariance Matrix");
    DFNCovarianceMatrixH[week]->SetTitle("Distort FN Covariance Matrix");
    DLiHeCovarianceMatrixH[week]=(TH2F*)BackgroundsCovarianceMatrixH[week]->Clone("Distort LiHe Covariance Matrix");
    DLiHeCovarianceMatrixH[week]->SetTitle("Distort LiHe Covariance Matrix");
    DAmCCovarianceMatrixH[week]=(TH2F*)BackgroundsCovarianceMatrixH[week]->Clone("Distort AmC Covariance Matrix");
    DAmCCovarianceMatrixH[week]->SetTitle("Distort AmC Covariance Matrix");
    
    IsotopeCovarianceMatrixH[week]=(TH2F*)SystematicCovarianceMatrixH[week]->Clone("Isotope Covariance Matrix");
    IsotopeCovarianceMatrixH[week]->SetTitle("Reactor Spectrum Covariance Matrix");
    ReactorPowerCovarianceMatrixH[week]=(TH2F*)SystematicCovarianceMatrixH[week]->Clone("Reactor Power Covariance Matrix");
    ReactorPowerCovarianceMatrixH[week]->SetTitle("Reactor Power Covariance Matrix");
    IAVCovarianceMatrixH[week]=(TH2F*)SystematicCovarianceMatrixH[week]->Clone("IAV Covariance Matrix");
    IAVCovarianceMatrixH[week]->SetTitle("IAV Covariance Matrix");
    NLCovarianceMatrixH[week]=(TH2F*)SystematicCovarianceMatrixH[week]->Clone("NL Covariance Matrix");
    NLCovarianceMatrixH[week]->SetTitle("NL Covariance Matrix");
    ResolutionCovarianceMatrixH[week]=(TH2F*)SystematicCovarianceMatrixH[week]->Clone("Resolution Covariance Matrix");
    ResolutionCovarianceMatrixH[week]->SetTitle("Resolution Covariance Matrix");

    GenerateStatisticalCovarianceMatrix(week);
    SaveStatisticalCovarianceMatrix(week);

//        TFile* SaveStatisticalCovMatrixF = TFile::Open("./CovarianceMatrices/SystematicCovarianceMatrixBeforeNormalize.root","recreate");
//        IsotopeCovarianceMatrixH[week]->Write();
//        ReactorPowerCovarianceMatrixH[week]->Write();
//        IAVCovarianceMatrixH[week]->Write();
//        NLCovarianceMatrixH[week]->Write();
//        ResolutionCovarianceMatrixH[week]->Write();
//        SaveStatisticalCovMatrixF->Close();
    
        x =0;
        y =0;
        
        for (Int_t neari=0; neari<MaxNear; neari++)
        {
            Int_t Ni1,Ni2,Ni3,Ni4;
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Ni1=1; Ni2=0; Ni3=0; Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<MaxFar; fari++)
            {
                Int_t Fi1,Fi2,Fi3,Fi4;
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Fi1=fari+1;Fi2 = 0;Fi3 = 0;Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFar; Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFar; Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=MaxFar; Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<MaxNear; nearj++)
                {
                    Int_t Nj1,Nj2,Nj3,Nj4;
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<MaxFar; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        Int_t Fj1,Fj2,Fj3,Fj4;
                        if(Nj1!=Nj2){Fj1=farj+1;Fj2 = 0;Fj3 = 0;Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFar; Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFar; Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=MaxFar; Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;

                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                
                                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                //                                                              Denormalize Oscillation in Systematic Covariance Matrix
                                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//                                IsotopeCovarianceMatrixH[week]->SetBinContent(x+1,y+1,IsotopeCovarianceMatrixH[week]->GetBinContent(x+1,y+1)*(PredictionVisH[fari][neari][week]->GetBinContent(i+1)*PredictionVisH[farj][nearj][week]->GetBinContent(j+1)));
                                IsotopeCovarianceMatrixM[x][y][week]=IsotopeCovarianceMatrixM[x][y][week]*(PredictionVisH[fari][neari][week]->GetBinContent(i+1)*PredictionVisH[farj][nearj][week]->GetBinContent(j+1));

//                                ReactorPowerCovarianceMatrixH[week]->SetBinContent(x+1,y+1,ReactorPowerCovarianceMatrixH[week]->GetBinContent(x+1,y+1)*(PredictionVisH[fari][neari][week]->GetBinContent(i+1)*PredictionVisH[farj][nearj][week]->GetBinContent(j+1)));
                                ReactorPowerCovarianceMatrixM[x][y][week]=ReactorPowerCovarianceMatrixM[x][y][week]*(PredictionVisH[fari][neari][week]->GetBinContent(i+1)*PredictionVisH[farj][nearj][week]->GetBinContent(j+1));
                                
//                                IAVCovarianceMatrixH[week]->SetBinContent(x+1,y+1,IAVCovarianceMatrixH[week]->GetBinContent(x+1,y+1)*(PredictionVisH[fari][neari][week]->GetBinContent(i+1)*PredictionVisH[farj][nearj][week]->GetBinContent(j+1)));
                                IAVCovarianceMatrixM[x][y][week]=IAVCovarianceMatrixM[x][y][week]*(PredictionVisH[fari][neari][week]->GetBinContent(i+1)*PredictionVisH[farj][nearj][week]->GetBinContent(j+1));
                                
//                                NLCovarianceMatrixH[week]->SetBinContent(x+1,y+1,NLCovarianceMatrixH[week]->GetBinContent(x+1,y+1)*(PredictionVisH[fari][neari][week]->GetBinContent(i+1)*PredictionVisH[farj][nearj][week]->GetBinContent(j+1)));
                                NLCovarianceMatrixM[x][y][week]=NLCovarianceMatrixM[x][y][week]*(PredictionVisH[fari][neari][week]->GetBinContent(i+1)*PredictionVisH[farj][nearj][week]->GetBinContent(j+1));
                                
//                                ResolutionCovarianceMatrixH[week]->SetBinContent(x+1,y+1,ResolutionCovarianceMatrixH[week]->GetBinContent(x+1,y+1)*(PredictionVisH[fari][neari][week]->GetBinContent(i+1)*PredictionVisH[farj][nearj][week]->GetBinContent(j+1)));
                                ResolutionCovarianceMatrixM[x][y][week]=ResolutionCovarianceMatrixM[x][y][week]*(PredictionVisH[fari][neari][week]->GetBinContent(i+1)*PredictionVisH[farj][nearj][week]->GetBinContent(j+1));
                                
                                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                //                                                             Add Systematic Covariance Matrix
                                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                SystematicCovarianceMatrixM[x][y][week]=
//                                IsotopeCovarianceMatrixM[x][y][week]+
//                                ReactorPowerCovarianceMatrixM[x][y][week]+
                                IAVCovarianceMatrixM[x][y][week]+
                                NLCovarianceMatrixM[x][y][week]+
                                ResolutionCovarianceMatrixM[x][y][week];
                                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                //                                                             Add Background Covariance Matrix
                                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                BackgroundsCovarianceMatrixM[x][y][week]=
                                VAccCovarianceMatrixM[x][y][week]+
                                VLiHeCovarianceMatrixM[x][y][week]+
                                VFNCovarianceMatrixM[x][y][week]+
                                VAmCCovarianceMatrixM[x][y][week]+
                                DLiHeCovarianceMatrixM[x][y][week]+
                                DFNCovarianceMatrixM[x][y][week]+
                                DAmCCovarianceMatrixM[x][y][week];
                                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                //                                                      Add all matrices into a Total Covariance Matrix
                                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                TotalCovarianceMatrixM.Set(MaxBins*MaxBins);

                                TotalCovarianceMatrixM[x*MaxBins+y]=
                                CovStat[x][y][week];
                                //+BackgroundsCovarianceMatrixM[x][y][week]+SystematicCovarianceMatrixM[x][y][week];
                            }
                        }
                    }
                }
            }
        }
    
//        TFile* SaveSystematicF = TFile::Open("./CovarianceMatrices/SystematicCovarianceMatrixNormalized.root","recreate");
//        IsotopeCovarianceMatrixH[week]->Write();
//        ReactorPowerCovarianceMatrixH[week]->Write();
//        IAVCovarianceMatrixH[week]->Write();
//        NLCovarianceMatrixH[week]->Write();
//        ResolutionCovarianceMatrixH[week]->Write();
//        SaveSystematicF->Close();
//    
//        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        //                                                                            Add all matrices into a Total Covariance Matrix
//        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//        BackgroundsCovarianceMatrixH[week]=(TH2F*)VAccCovarianceMatrixH[week]->Clone("Background Covariance Matrix");
//        BackgroundsCovarianceMatrixH[week]->Add(VFNCovarianceMatrixH[week]);
//        BackgroundsCovarianceMatrixH[week]->Add(DFNCovarianceMatrixH[week]);
//        BackgroundsCovarianceMatrixH[week]->Add(VLiHeCovarianceMatrixH[week]);
//        BackgroundsCovarianceMatrixH[week]->Add(DLiHeCovarianceMatrixH[week]);
//        BackgroundsCovarianceMatrixH[week]->Add(VAmCCovarianceMatrixH[week]);
//        BackgroundsCovarianceMatrixH[week]->Add(DAmCCovarianceMatrixH[week]);

//        SystematicCovarianceMatrixH[week]=(TH2F*)IsotopeCovarianceMatrixH[week]->Clone();// Don't include yet, need to apply cholevski to get varied spectra from Christine model.
//        SystematicCovarianceMatrixH[week]->Add(ReactorPowerCovarianceMatrixH[week]);
//        SystematicCovarianceMatrixH[week]->Add(IAVCovarianceMatrixH[week]);
//        SystematicCovarianceMatrixH[week]=(TH2F*)(IAVCovarianceMatrixH[week]->Clone("Systematic Covariance Matrix"));
//        SystematicCovarianceMatrixH[week]->Add(NLCovarianceMatrixH[week]);
//        SystematicCovarianceMatrixH[week]->Add(ResolutionCovarianceMatrixH[week]);
    
//
//     // TotalCovarianceMatrixH[week]=(TH2F*)BackgroundsCovarianceMatrixH[week]->Clone();
//      //   TotalCovarianceMatrixH[week]=(TH2F*)SystematicCovarianceMatrixH[week]->Clone();
//
//     // TotalCovarianceMatrixH[week]->Add(SystematicCovarianceMatrixH[week]);//This should be multiplied in the fitter process!
//        TotalCovarianceMatrixH[week]=(TH2F*)StatCov2H[week]->Clone("Total Covariance Matrix");
//       // TotalCovarianceMatrixH[week]->Add(StatCov2H[week]);
}

void Fitter :: NormCov(TH2F* Histo)
{
    TH1F* CopyHisto = (TH1F*)Histo->Clone();
    x =0;
    y =0;
    
    for (Int_t neari=0; neari<MaxNear; neari++)
    {
        Int_t Ni1,Ni2,Ni3,Ni4;
        //Logic for the 2D matrix index done up to 8 ADs
        if(neari==0){Ni1=1; Ni2=0; Ni3=0; Ni4=0;}
        if(neari==1){Ni2++;}
        if(neari==2){Ni3++;}
        if(neari==3){Ni4++;}
        
        for (Int_t fari=0; fari<MaxFar; fari++)
        {
            Int_t Fi1,Fi2,Fi3,Fi4;
            //Logic for the 2D matrix index done up to 8 ADs
            if(Ni1!=Ni2){Fi1=fari+1;Fi2 = 0;Fi3 = 0;Fi4 = 0;}
            if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFar; Fi2=fari+1;}
            if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFar; Fi3=fari+1;}
            if(Ni3==Ni4&&Ni4==1){Fi3=MaxFar; Fi4=fari+1;}
            
            for (Int_t nearj=0; nearj<MaxNear; nearj++)
            {
                Int_t Nj1,Nj2,Nj3,Nj4;
                //Logic for the 2D matrix index done up to 8 ADs
                if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                if(nearj==1){Nj2++;}
                if(nearj==2){Nj3++;}
                if(nearj==3){Nj4++;}
                
                for (Int_t farj=0; farj<MaxFar; farj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    Int_t Fj1,Fj2,Fj3,Fj4;
                    if(Nj1!=Nj2){Fj1=farj+1;Fj2 = 0;Fj3 = 0;Fj4 = 0;}
                    if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFar; Fj2=farj+1;}
                    if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFar; Fj3=farj+1;}
                    if(Nj3==Nj4&&Nj4==1){Fj3=MaxFar; Fj4=farj+1;}
                    
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {//columns
                        x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;

                        for (Int_t j = 0; j<n_evis_bins; j++)
                        {//rows
                            y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                            
                            if((Histo->GetBinContent(x+1,x+1)*Histo->GetBinContent(y+1,y+1))==0)
                            {
                                std::cout << "NAN IN COV CALCULATION" << std::endl;
                                Histo->SetBinContent(x+1,y+1,1);//To avoid (0/0) nans when bin contents are practically the same; (this happens when bins are empty, in this case the correlation should be 1)
                                //                                    std::cout << "Norm cov tried to be inf" << std::endl;
                            }
                            else
                            {
                               Histo->SetBinContent(x+1,y+1,CopyHisto->GetBinContent(x+1,y+1)/(sqrt(CopyHisto->GetBinContent(x+1,x+1)*CopyHisto->GetBinContent(y+1,y+1))));
                            }
                        }
                    }
                }
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Produces the Statistical Covariance Matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Fitter :: GenerateStatisticalCovarianceMatrix(Int_t week)
{
    StatCov2H[week] = new TH2F(Form("Statistical Covariance Matrix%d",week),Form("Statistical Covariance Matrix%d",week),MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    LoadBackgrounds();
    //Nominal backgrounds:
    
    for (Int_t AD=0; AD<NADs; AD++)
    {
        //Nominal background spectrum
        BackgroundSpectrumH[AD][week]=(TH1F*)AccidentalsH[0][0]->Clone();
        BackgroundSpectrumH[AD][week]->Reset();
        BackgroundSpectrumH[AD][week]->Add(AccidentalsH[AD][week]);
        BackgroundSpectrumH[AD][week]->Add(LiHeH[AD][week]);
        BackgroundSpectrumH[AD][week]->Add(FastNeutronsH[AD][week]);
        BackgroundSpectrumH[AD][week]->Add(AmCH[AD][week]);
        //            BackgroundSpectrumH[AD][week]->Write();
    }
    
    //Combine matrices in 9x9, 2x2 or 1x1 prediction
    if (Combine == 1)
    {
        NearBackgroundSpectrumH[0][week]=(TH1F*)BackgroundSpectrumH[0][week]->Clone();
        NearBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[1][week]);//All near hall detectors together
        NearBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[2][week]);
        NearBackgroundSpectrumH[0][week]->Scale(1/3);
        
        FarBackgroundSpectrumH[0][week]=(TH1F*)BackgroundSpectrumH[3][week]->Clone();
        FarBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[4][week]);//All far hall detectors together
        FarBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[5][week]);
        FarBackgroundSpectrumH[0][week]->Scale(1/3);
    }
    else if(Combine == 2)
    {
        
        FarBackgroundSpectrumH[0][week]=(TH1F*)BackgroundSpectrumH[3][week]->Clone();//All far hall detectors together
        FarBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[4][week]);
        FarBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[5][week]);
        FarBackgroundSpectrumH[0][week]->Scale(1/3);

        NearBackgroundSpectrumH[0][week]=(TH1F*)BackgroundSpectrumH[0][week]->Clone();//DB
        NearBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[1][week]);
        NearBackgroundSpectrumH[0][week]->Scale(0.5);
        NearBackgroundSpectrumH[1][week]=(TH1F*)BackgroundSpectrumH[2][week]->Clone();//LO
    }
    else
    {
        NearBackgroundSpectrumH[0][week]=(TH1F*)BackgroundSpectrumH[0][week]->Clone();
        NearBackgroundSpectrumH[1][week]=(TH1F*)BackgroundSpectrumH[1][week]->Clone();
        NearBackgroundSpectrumH[2][week]=(TH1F*)BackgroundSpectrumH[2][week]->Clone();
        
        FarBackgroundSpectrumH[0][week]=(TH1F*)BackgroundSpectrumH[3][week]->Clone();
        FarBackgroundSpectrumH[1][week]=(TH1F*)BackgroundSpectrumH[4][week]->Clone();
        FarBackgroundSpectrumH[2][week]=(TH1F*)BackgroundSpectrumH[5][week]->Clone();
    }
    
    for (Int_t far=0; far<MaxFar; far++)
    {
        for (Int_t near=0; near<MaxNear; near++)
        {
            for (Int_t pts = 0; pts < n_evis_bins; pts++)
            {
                Sigma_Far[far][near][week][pts]=sqrt(PredictionVisH[far][near][week]->GetBinContent(pts+1)+FarBackgroundSpectrumH[far][week]->GetBinContent(pts+1));
                
                //N OBSERVED, F PREDICTED FROM N OBSERVED, BACKGROUNDS PREDICTED (WITHOUT ANY FLUCTUATION)
                Sigma_Near[far][near][week][pts]=(PredictionVisH[far][near][week]->GetBinContent(pts+1)/ NearDataH[near][week]->GetBinContent(pts+1))*sqrt(NearDataH[near][week]->GetBinContent(pts+1)+NearBackgroundSpectrumH[near][week]->GetBinContent(pts+1));
            }
        }
    }
    
    Int_t x =0;
    Int_t y =0;
    
    for (Int_t neari=0; neari<MaxNear; neari++)
    {
        Int_t Ni1,Ni2,Ni3,Ni4;
        //Logic for the 2D matrix index done up to 8 ADs
        if(neari==0){Ni1=1; Ni2=0; Ni3=0; Ni4=0;}
        if(neari==1){Ni2++;}
        if(neari==2){Ni3++;}
        if(neari==3){Ni4++;}
        
        for (Int_t fari=0; fari<MaxFar; fari++)
        {
            Int_t Fi1,Fi2,Fi3,Fi4;
            //Logic for the 2D matrix index done up to 8 ADs
            if(Ni1!=Ni2){Fi1=fari+1;Fi2 = 0;Fi3 = 0;Fi4 = 0;}
            if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFar; Fi2=fari+1;}
            if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFar; Fi3=fari+1;}
            if(Ni3==Ni4&&Ni4==1){Fi3=MaxFar; Fi4=fari+1;}
            
            for (Int_t nearj=0; nearj<MaxNear; nearj++)
            {
                Int_t Nj1,Nj2,Nj3,Nj4;
                //Logic for the 2D matrix index done up to 8 ADs
                if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                if(nearj==1){Nj2++;}
                if(nearj==2){Nj3++;}
                if(nearj==3){Nj4++;}
                
                for (Int_t farj=0; farj<MaxFar; farj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    Int_t Fj1,Fj2,Fj3,Fj4;
                    if(Nj1!=Nj2){Fj1=farj+1;Fj2 = 0;Fj3 = 0;Fj4 = 0;}
                    if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFar; Fj2=farj+1;}
                    if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFar; Fj3=farj+1;}
                    if(Nj3==Nj4&&Nj4==1){Fj3=MaxFar; Fj4=farj+1;}
                    
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {//columns
                        x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;

                        for (Int_t j = 0; j<n_evis_bins; j++)
                        {//rows
                            y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                            
                            //Near component correlated
                            if(neari==nearj && fari!=farj)
                            {
                                CovStat[x][y][week]=Sigma_Near[fari][neari][week][i]*Sigma_Near[farj][nearj][week][j];
                            }
                            //Far component correlated
                            if(fari==farj && neari!=nearj)
                            {
                                CovStat[x][y][week]=Sigma_Far[fari][neari][week][i]*Sigma_Far[farj][nearj][week][j];
                            }
                            if(neari==nearj && fari==farj)
                            {
                                //General covariance
                                CovStat[x][y][week]=(Sigma_Near[fari][neari][week][i]*Sigma_Near[fari][neari][week][j])+(Sigma_Far[farj][nearj][week][i]*Sigma_Far[farj][nearj][week][j]);
                            }
                            //Uncorrelated terms
                            if(neari!=nearj && fari!=farj)
                            {
                                CovStat[x][y][week]=0;
                            }
                            StatCov2H[week]->SetBinContent(x+1,y+1,CovStat[x][y][week]);
                        }
                    }
                }
            }
        }
    }
    TFile* SaveNoNormStatF = TFile::Open("./CovarianceMatrices/StatisticalCovarianceMatrix.root","update");
    StatCov2H[week]->Write();
    SaveNoNormStatF->Close();
    
//    NormCov(StatCov2H[week]);
    // StatCov2H[week]->Draw("colz");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Save Statistical Covariance Matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Fitter :: SaveStatisticalCovarianceMatrix(Int_t week)
{
    //Save statistical matrix
    TFile* SaveStatisticalCovMatrixF = TFile::Open("./CovarianceMatrices/StatisticalCovarianceMatrix.root","update");
    StatCov2H[week]->Write();
    
    for (Int_t near = 0; near<MaxNear; near++)
    {
        NearDataH[near][week]->Write();
        for (Int_t far =0; far<MaxFar; far++)
        {
            FarDataH[far][near][week]->Write();
            PredictionVisH[far][near][week]->Write();
        }
    }
    
    SaveStatisticalCovMatrixF->Close();
    
    //Save in a txt file
    if(WriteOutput)
    {
        ofstream statf("CovarianceMatrices/StatisticalCovarianceMatrix.txt");
        Int_t x =0;
        Int_t y =0;
        
        for (Int_t neari=0; neari<MaxNear; neari++)
        {
            Int_t Ni1,Ni2,Ni3,Ni4;
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Ni1=1; Ni2=0; Ni3=0; Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<MaxFar; fari++)
            {
                Int_t Fi1,Fi2,Fi3,Fi4;
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Fi1=fari+1;Fi2 = 0;Fi3 = 0;Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFar;Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFar;Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=MaxFar;Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<MaxNear; nearj++)
                {
                    Int_t Nj1,Nj2,Nj3,Nj4;
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<MaxFar; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        Int_t Fj1,Fj2,Fj3,Fj4;
                        if(Nj1!=Nj2){Fj1=farj+1;Fj2 = 0;Fj3 = 0;Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFar; Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFar; Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=MaxFar; Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;

                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                
                                statf << StatCov2H[week]->GetBinContent(x+1,y+1) << " ";
                            }
                        }
                    }
                }
            }
        }
        statf << std::endl;
        statf.close();
    }
    delete StatCov2H[week];
}

void Fitter :: InvertMatrix()
{
    for (Int_t week = 0; week<Nweeks; week++)
    {
//        InvTotalCovarianceMatrixH[week] = (TH2F*)TotalCovarianceMatrixH[week]->Clone();
        
//        InvTotalCovarianceMatrixH[week]->Reset();
        
//        for (Int_t i = 0; i < MaxBins; i++)
//        {
//            for (Int_t j = 0; j < MaxBins; j++)
//            {
//                TotalCovarianceMatrixArray[i*MaxBins+j] = TotalCovarianceMatrixH[week]->GetBinContent(i+1,j+1);
//                //            std::cout << "ELEMENTS TOTAL COV" << TotalCovarianceMatrixArray[i*MaxBins+j] << std::endl;
//            }
//        }

        CombineMatrices(week);
        TMatrixD* TotalCovarianceMatrix = new TMatrixD(MaxBins,MaxBins);
//        TotalCovarianceMatrix->SetMatrixArray(TotalCovarianceMatrixArray);
        TotalCovarianceMatrix->SetMatrixArray(TotalCovarianceMatrixM.GetArray());

        TotalCovarianceMatrix->Invert();
        
        Double_t* InvTotalCovarianceMatrixArray = TotalCovarianceMatrix->GetMatrixArray();
        
//        for (Int_t i = 0; i < MaxBins; i++)
//        {
//            for (Int_t j = 0; j < MaxBins; j++)
//            {
//                InvTotalCovarianceMatrixH[week]->SetBinContent(i+1,j+1,InvTotalCovarianceMatrixArray[i*MaxBins+j]);
//            }
//        }
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j < MaxBins; j++)
            {
                InvTotalCovarianceMatrixM[i][j][week]=InvTotalCovarianceMatrixArray[(j+i*MaxBins)];
            }
        }

        delete TotalCovarianceMatrix;
        
        SaveCovarianceMatrices(week);
    }
}

void Fitter :: SaveChiSquare(Int_t week)
{
    std::cout << "SAVING HISTOGRAMS" << std::endl;
    TFile* SaveF = TFile::Open("./RootOutputs/ChiSquare.root","recreate");
    ChiSquareH[week]->Write();
    SaveF->Close();
}

void Fitter :: LoadBackgrounds()
{
    TFile* BackgroundsF = TFile::Open("./BackgroundSpectrum/Backgrounds.root");
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t AD =0; AD<NADs; AD++)
        {
            AccidentalsH[AD][week]= (TH1F*)gDirectory->Get(Form("Accidentals_AD%i",AD+1));
            AccidentalsH[AD][week]->Scale(ScaleAcc[AD][week]/AccidentalsH[AD][week]->Integral());
            LiHeH[AD][week]=(TH1F*)gDirectory->Get("LiHe");//Missing LiHe inputs so far in Hydrogen Analysis
            LiHeH[AD][week]->Scale(ScaleLiHe[AD][week]/LiHeH[AD][week]->Integral());
            FastNeutronsH[AD][week]=(TH1F*)gDirectory->Get("FN");
            FastNeutronsH[AD][week]->Scale(ScaleFN[AD][week]/FastNeutronsH[AD][week]->Integral());
            AmCH[AD][week]=(TH1F*)gDirectory->Get("AmC");
            AmCH[AD][week]->Scale(ScaleAmC[AD][week]/AmCH[AD][week]->Integral());
        }
    }
    BackgroundsF->Close();
}

void Fitter :: ApplyStatisticalFluctuation(TH1F* Histo)
{
    for(Int_t VisibleEnergyIndex=1;VisibleEnergyIndex<=Histo->GetXaxis()->GetNbins();VisibleEnergyIndex++)
    {
        Histo->SetBinContent(VisibleEnergyIndex,(Double_t)(rand->Poisson(Histo->GetBinContent(VisibleEnergyIndex))));
    }
}

