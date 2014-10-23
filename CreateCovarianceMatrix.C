#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH2.h"
#include "TMath.h"
#include <vector>
#include "NominalData.h"
#include <math.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TCanvas.h"
#include "NominalData.h"
#include "TAxis.h"
#include "TArrayD.h"

class CreateCovarianceMatrix()
{
    Int_t n_evis_bins;
    Int_t n_etrue_bins;
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.
    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    Int_t MaxBins;
    Int_t Combine;
    
private:
    NominalData* Data;
public:
    CreateCovarianceMatrix();
    CreateCovarianceMatrix(NominalData*);
    ~CreateCovarianceMatrix();
}

CreateCovarianceMatrix :: CreateCovarianceMatrix()
{
}
CreateCovarianceMatrix :: CreateCovarianceMatrix(NominalData* Data)
{
    Combine = Data->GetCombineMode();
    //  Linear binning
    if(LinearBinning)
    {
        n_evis_bins = Data->GetNbins();
    }
    //  Non-linear binning
    else
    {
        n_evis_bins=37;
    }
    
    MaxBins = 9*n_evis_bins;
}

CreateCovarianceMatrix :: ~CreateCovarianceMatrix()
{
}

void CreateCovarianceMatrix :: LoadRootCovarianceMatrices()
{
    
}

void CreateCovarianceMatrix :: LoadTxTCovarianceMatrices()
{
    
}

void CreateCovarianceMatrix :: CreateCovarianceMatrixMain()
{
    std::cout <<  "***********************************************************************************************" << std::endl;
    std::cout <<  "\t Summing Matrices " << std::endl;
    
    TH2D* BackgroundsCovarianceMatrixH = new TH2D("Background Covariance Matrix","Background Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    TH2D* SystematicCovarianceMatrixH = new TH2D("Systematic Covariance Matrix","Systematic Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);

    TH2D* IsotopeCovarianceMatrixH = new TH2D(" Isotope Covariance Matrix"," Isotope Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    TH2D* ReactorPowerCovarianceMatrixH = new TH2D(" Reactor Power Covariance Matrix"," Reactor Power Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    TH2D* RelativeEnergyScaleCovarianceMatrixH = new TH2D(" Relative Energy Scale Covariance Matrix"," Relative Energy Scale Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    TH2D* RelativeEnergyOffsetCovarianceMatrixH = new TH2D(" Relative Energy Offset Covariance Matrix"," Relative Energy Offset Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    TH2D* AbsoluteEnergyScaleCovarianceMatrixH = new TH2D(" Absolute Energy Scale Covariance Matrix"," Absolute Energy Scale Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    TH2D* AbsoluteEnergyOffsetCovarianceMatrixH = new TH2D(" Absolute Energy Offset Covariance Matrix"," Absolute Energy Offset Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    TH2D* IAVCovarianceMatrixH = new TH2D(" IAV Covariance Matrix"," IAV Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    TH2D* NLCovarianceMatrixH = new TH2D(" NL Covariance Matrix"," NL Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    TH2D* ResolutionCovarianceMatrixH = new TH2D(" Resolution Covariance Matrix"," Resolution Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    TH2D* Sin22t12CovarianceMatrixH = new TH2D(" Sin22t12 Covariance Matrix"," Sin22t12 Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    
    std::cout <<  "\t Matrices Created " << std::endl;
    
    x =0;
    y =0;
    
    for (Int_t neari=0; neari<MaxNearOscModel; neari++)
    {
        Int_t Ni1,Ni2,Ni3,Ni4;
        //Logic for the 2D matrix index done up to 8 ADs
        if(neari==0){Ni1=1; Ni2=0; Ni3=0; Ni4=0;}
        if(neari==1){Ni2++;}
        if(neari==2){Ni3++;}
        if(neari==3){Ni4++;}
        
        for (Int_t fari=0; fari<MaxFarOscModel; fari++)
        {
            Int_t Fi1,Fi2,Fi3,Fi4;
            //Logic for the 2D matrix index done up to 8 ADs
            if(Ni1!=Ni2){Fi1=fari+1;Fi2 = 0;Fi3 = 0;Fi4 = 0;}
            if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFarOscModel; Fi2=fari+1;}
            if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFarOscModel; Fi3=fari+1;}
            if(Ni3==Ni4&&Ni4==1){Fi3=MaxFarOscModel; Fi4=fari+1;}
            
            for (Int_t nearj=0; nearj<MaxNearOscModel; nearj++)
            {
                Int_t Nj1,Nj2,Nj3,Nj4;
                //Logic for the 2D matrix index done up to 8 ADs
                if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                if(nearj==1){Nj2++;}
                if(nearj==2){Nj3++;}
                if(nearj==3){Nj4++;}
                
                for (Int_t farj=0; farj<MaxFarOscModel; farj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    Int_t Fj1,Fj2,Fj3,Fj4;
                    if(Nj1!=Nj2){Fj1=farj+1;Fj2 = 0;Fj3 = 0;Fj4 = 0;}
                    if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFarOscModel; Fj2=farj+1;}
                    if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFarOscModel; Fj3=farj+1;}
                    if(Nj3==Nj4&&Nj4==1){Fj3=MaxFarOscModel; Fj4=farj+1;}
                    
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {//columns
                        x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                        
                        for (Int_t j = 0; j<n_evis_bins; j++)
                        {//rows
                            y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                            
                            if(ReadTxt)
                            {
                                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                //                                                             Add Systematic Covariance Matrix
                                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                SystematicCovarianceMatrixM[x+y*MaxBins]=
                                IsotopeCovarianceMatrixM[x+y*MaxBins]+
                                ReactorPowerCovarianceMatrixM[x+y*MaxBins]+
                                RelativeEnergyScaleCovarianceMatrixM[x+y*MaxBins]+
                                RelativeEnergyOffsetCovarianceMatrixM[x+y*MaxBins]+
                                AbsoluteEnergyScaleCovarianceMatrixM[x+y*MaxBins]+
                                AbsoluteEnergyOffsetCovarianceMatrixM[x+y*MaxBins]+
                                IAVCovarianceMatrixM[x+y*MaxBins]+
                                NLCovarianceMatrixM[x+y*MaxBins]+
                                ResolutionCovarianceMatrixM[x+y*MaxBins]+
                                Sin22t12CovarianceMatrixM[x+y*MaxBins];
                                
                                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                //                                                             Add Background Covariance Matrix
                                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                BackgroundsCovarianceMatrixM[x+y*MaxBins]=
                                VAccCovarianceMatrixM[x+y*MaxBins]+
                                VLiHeCovarianceMatrixM[x+y*MaxBins]+
                                VFNCovarianceMatrixM[x+y*MaxBins]+
                                VAmCCovarianceMatrixM[x+y*MaxBins]+
                                DLiHeCovarianceMatrixM[x+y*MaxBins]+
                                DFNCovarianceMatrixM[x+y*MaxBins]+
                                DAmCCovarianceMatrixM[x+y*MaxBins];
                                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                //                                                      Add all matrices into a Total Covariance Matrix
                                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                
                                TotalCovarianceMatrixM[x*MaxBins+y]=CovStat[x+y*MaxBins]+BackgroundsCovarianceMatrixM[x+y*MaxBins]+SystematicCovarianceMatrixM[x+y*MaxBins];
                                
                                //  Fill Histograms
                                TotalCovarianceMatrixH->SetBinContent(x+1,y+1,TotalCovarianceMatrixM[x*MaxBins+y]);
                                BackgroundsCovarianceMatrixH->SetBinContent(x+1,y+1,BackgroundsCovarianceMatrixM[x+y*MaxBins]);
                                SystematicCovarianceMatrixH->SetBinContent(x+1,y+1,SystematicCovarianceMatrixM[x+y*MaxBins]);
                                
                                VAccCovarianceMatrixH->SetBinContent(x+1,y+1,VAccCovarianceMatrixM[x+y*MaxBins]);
                                VLiHeCovarianceMatrixH->SetBinContent(x+1,y+1,VLiHeCovarianceMatrixM[x+y*MaxBins]);
                                VFNCovarianceMatrixH->SetBinContent(x+1,y+1,VFNCovarianceMatrixM[x+y*MaxBins]);
                                VAmCCovarianceMatrixH->SetBinContent(x+1,y+1,VAmCCovarianceMatrixM[x+y*MaxBins]);
                                DLiHeCovarianceMatrixH->SetBinContent(x+1,y+1,DLiHeCovarianceMatrixM[x+y*MaxBins]);
                                DFNCovarianceMatrixH->SetBinContent(x+1,y+1,DFNCovarianceMatrixM[x+y*MaxBins]);
                                DAmCCovarianceMatrixH->SetBinContent(x+1,y+1,DAmCCovarianceMatrixM[x+y*MaxBins]);
                                
                                IsotopeCovarianceMatrixH->SetBinContent(x+1,y+1,IsotopeCovarianceMatrixM[x+y*MaxBins]);
                                ReactorPowerCovarianceMatrixH->SetBinContent(x+1,y+1,ReactorPowerCovarianceMatrixM[x+y*MaxBins]);
                                RelativeEnergyScaleCovarianceMatrixH->SetBinContent(x+1,y+1,RelativeEnergyScaleCovarianceMatrixM[x+y*MaxBins]);
                                RelativeEnergyOffsetCovarianceMatrixH->SetBinContent(x+1,y+1,RelativeEnergyOffsetCovarianceMatrixM[x+y*MaxBins]);
                                AbsoluteEnergyScaleCovarianceMatrixH->SetBinContent(x+1,y+1,AbsoluteEnergyScaleCovarianceMatrixM[x+y*MaxBins]);
                                AbsoluteEnergyOffsetCovarianceMatrixH->SetBinContent(x+1,y+1,AbsoluteEnergyOffsetCovarianceMatrixM[x+y*MaxBins]);
                                IAVCovarianceMatrixH->SetBinContent(x+1,y+1,IAVCovarianceMatrixM[x+y*MaxBins]);
                                NLCovarianceMatrixH->SetBinContent(x+1,y+1,NLCovarianceMatrixM[x+y*MaxBins]);
                                ResolutionCovarianceMatrixH->SetBinContent(x+1,y+1,ResolutionCovarianceMatrixM[x+y*MaxBins]);
                                Sin22t12CovarianceMatrixH->SetBinContent(x+1,y+1,Sin22t12CovarianceMatrixM[x+y*MaxBins]);
                            }
                            else
                            {
                                IsotopeCovarianceMatrixH->SetBinContent(x+1,y+1,IsotopeCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari+MaxFarOscModel*neari]->GetBinContent(i+1)*PredictionH[farj+MaxFarOscModel*nearj]->GetBinContent(j+1)));
                                
                                ReactorPowerCovarianceMatrixH->SetBinContent(x+1,y+1,ReactorPowerCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari+MaxFarOscModel*neari]->GetBinContent(i+1)*PredictionH[farj+MaxFarOscModel*nearj]->GetBinContent(j+1)));
                                
                                RelativeEnergyScaleCovarianceMatrixH->SetBinContent(x+1,y+1,RelativeEnergyScaleCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari+MaxFarOscModel*neari]->GetBinContent(i+1)*PredictionH[farj+MaxFarOscModel*nearj]->GetBinContent(j+1)));
                                
                                RelativeEnergyOffsetCovarianceMatrixH->SetBinContent(x+1,y+1,RelativeEnergyOffsetCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari+MaxFarOscModel*neari]->GetBinContent(i+1)*PredictionH[farj+MaxFarOscModel*nearj]->GetBinContent(j+1)));
                                
                                AbsoluteEnergyScaleCovarianceMatrixH->SetBinContent(x+1,y+1,AbsoluteEnergyScaleCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari+MaxFarOscModel*neari]->GetBinContent(i+1)*PredictionH[farj+MaxFarOscModel*nearj]->GetBinContent(j+1)));
                                
                                AbsoluteEnergyOffsetCovarianceMatrixH->SetBinContent(x+1,y+1,AbsoluteEnergyOffsetCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari+MaxFarOscModel*neari]->GetBinContent(i+1)*PredictionH[farj+MaxFarOscModel*nearj]->GetBinContent(j+1)));
                                
                                IAVCovarianceMatrixH->SetBinContent(x+1,y+1,IAVCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari+MaxFarOscModel*neari]->GetBinContent(i+1)*PredictionH[farj+MaxFarOscModel*nearj]->GetBinContent(j+1)));
                                
                                NLCovarianceMatrixH->SetBinContent(x+1,y+1,NLCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari+MaxFarOscModel*neari]->GetBinContent(i+1)*PredictionH[farj+MaxFarOscModel*nearj]->GetBinContent(j+1)));
                                
                                ResolutionCovarianceMatrixH->SetBinContent(x+1,y+1,ResolutionCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari+MaxFarOscModel*neari]->GetBinContent(i+1)*PredictionH[farj+MaxFarOscModel*nearj]->GetBinContent(j+1)));
                                
                                Sin22t12CovarianceMatrixH->SetBinContent(x+1,y+1,Sin22t12CovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari+MaxFarOscModel*neari]->GetBinContent(i+1)*PredictionH[farj+MaxFarOscModel*nearj]->GetBinContent(j+1)));
                                
                            }
                        }
                    }
                }
            }
        }
}