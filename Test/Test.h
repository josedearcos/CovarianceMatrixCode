#define ActivateTest

#pragma once
#include "TColor.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2.h"
#include "TMath.h"
#include <math.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TCanvas.h"
#include "TF1.h"
#include <vector>
#include "TRandom3.h"
#include "TTree.h"
#include <sstream>
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "FitBackgrounds2.h"
#include "NominalData.h"
#include "CrossSection.h"
#include "ReactorSpectrumMultiple.h"
#include "AntineutrinoSpectrum.h"
#include "Oscillation.h"
#include "OscillationReactor.h"
#include "Prediction.h"
#include "CovarianceMatrix3.h"


const bool TestExternalInputs = 0;

const Int_t Combine = 2;
const Int_t P12E=0;

const Int_t Nweeks = 1;
const Int_t ADsEH1=2;
const Int_t ADsEH2=1;
const Int_t ADsEH3=3;
const Int_t NADs = 6;
const Int_t MaxNearOscModel = ADsEH1+ADsEH2;
const Int_t MaxFarOscModel = ADsEH3;
const Int_t n_evis_bins=37;
const Int_t n_etrue_bins =39;
const Int_t MaxBins = 9*n_evis_bins;
const Double_t InitialEnergy = 1.8;
const Double_t FinalVisibleEnergy = 12.0;

class Test
{
public:
    void TestCovarianceMatrixConstruction();//  Tested, works
    void TestCovarianceMatrixSamples();
    void TestCombineMain();//   Tested, works
    Test();
    void TestInputData();// Tested, works
    void TestRebinMatrix();//Tested, rebin works.
    void TestTree();
    void CompareWithHenoch();
    void VectorOnes();
    void LUDecomp();
    void LBNLPredictions();
    void FitterMainTests();
    void TestAllPredictions(bool);
    void TestSuperHistograms();//Superhistograms in LBNL are the same between different ADs, I don't understand this! It doesn't agree with their definition
    void LoganCrossCheck();//Test logan Toy MC output and mine.
    void Near_To_Far_Prediction_VS_Far_Prediction();
    
private:
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.
    
    std::vector<TH1D*> CombineData(std::vector<TH1D*>,std::vector<TH1D*>);
    Int_t Period;
    Int_t NReactorPeriods;
    Int_t NADs;
    Int_t NSamples;
    Int_t CombineMode;
    Int_t NFits;
    bool Analysis;
    bool flagCombine;
    bool flagNear;
    bool flagFar;
    bool Binning;
    Int_t PlotBin;
    bool Automatic;
    bool ToyMC;
    bool Minuit;
    bool deleteFlag;
    bool deleteFlagSpec;
    Int_t MaxFar;
    Int_t MaxNear;
    Int_t NearTrueIndex;
    
    Double_t Sin22t13SliderValue;
    Double_t DeltaMSilderValue;
    
    bool UseToyMCTree;
    
    bool Fit2D;
    bool FitSin22t13;
    //Manual Control:
    //To calculate Covariance Matrices set to 1. Activate only one at a time.
    //Backgrounds
    bool VaryAccidentalMatrix;
    bool VaryLiHeMatrix;
    bool VaryFastNeutronsMatrix;
    bool VaryAmCMatrix;
    bool DistortLiHeMatrix;
    bool DistortFastNeutronsMatrix;
    bool DistortAmCMatrix;
    //Systematics
    bool IsotopeMatrix;
    bool ReactorPowerMatrix;
    bool EnergyScaleMatrix;
    //        bool EnergyOffsetMatrix;
    //        bool AbsoluteScaleMatrix;
    //        bool AbsoluteOffsetMatrix;
    bool IAVMatrix;
    bool NLMatrix;
    bool ResolutionMatrix;
    bool Sin22t12Matrix;
    bool EfficiencyMatrix;
    
    bool PlotCovariance;
    bool NL[3];
    bool FitterMode[2];
    bool Fitter1DMode[2];
    
    bool AutomaticBudget; //  To loop the fitter over all possible systematics
    bool TurnOnBudget; // Generates Turn On Error Budget
    bool TurnOffBudget;// Generates Turn Off Error Budget
    bool StatisticalFluctuation;
    
    std::string ResponseMatrixDirectory;
    std::string ToyMCSampleDirectory;
    std::string NominalPredictionsDirectory;
    std::string SysCovarianceMatrixDirectory;
    std::string BkgCovarianceMatrixDirectory;
};
Test :: Test()
{
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);
    
    for (Int_t i = 0; i <= n_etrue_bins; i++)
    {
        enu_bins[i] = 0.2 * i + InitialEnergy;
    }
    
    evis_bins[0] = 0.7;
    
    for (Int_t i = 0; i < n_evis_bins-1; i++)
    {
        evis_bins[i+1] = 0.2 * i + 1.0;
        
        //        std::cout << evis_bins[i+1] << std::endl;
    }
    evis_bins[n_evis_bins] = FinalVisibleEnergy;
}

void Test::TestCombineMain()
{
    Int_t CombineMaxBins;
    Double_t TotalCovarianceMatrixM[9*n_evis_bins*9*n_evis_bins];
    Double_t final_covmatrix[9*n_evis_bins*9*n_evis_bins];
    
    Int_t CombineMaxFar;
    Int_t CombineMaxNear;
    
    if (Combine == 0)
    {
        CombineMaxFar = ADsEH3;
        CombineMaxNear = ADsEH1+ADsEH2;
        CombineMaxBins = 9*n_evis_bins;
    }
    if (Combine == 1)
    {
        CombineMaxFar = 1;
        CombineMaxNear = 1;
        CombineMaxBins = n_evis_bins;
    }
    if (Combine == 2)
    {
        CombineMaxFar = 1;
        CombineMaxNear = 2;
        CombineMaxBins = 2*n_evis_bins;
    }
    
    Int_t x =0;
    Int_t y =0;
    
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
            if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFarOscModel;Fi2=fari+1;}
            if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFarOscModel;Fi3=fari+1;}
            if(Ni3==Ni4&&Ni4==1){Fi3=MaxFarOscModel;Fi4=fari+1;}
            
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
                            
                            TotalCovarianceMatrixM[y+MaxBins*x]=1;
                            
                            if(Combine == 0)
                            {
                                final_covmatrix[y+CombineMaxBins*x] = TotalCovarianceMatrixM[y+MaxBins*x];
                            }
                            else if(Combine == 1)
                            {//CombineMaxBins = n_evis_bins
                                if(x%n_evis_bins == i && y%n_evis_bins ==j)
                                {
                                    //                                    if(x%n_evis_bins == 0 && y%n_evis_bins ==0)
                                    //                                    {
                                    ////                                        std::cout << count++ << std:: endl;
                                    //                                    }
                                    final_covmatrix[j+CombineMaxBins*i] = final_covmatrix[j+CombineMaxBins*i] + 1./(ADsEH1+ADsEH2)*1./(ADsEH1+ADsEH2)*TotalCovarianceMatrixM[y+MaxBins*x];
                                }
                            }
                            else if(Combine == 2)
                            {
                                if(neari<ADsEH1 && nearj<ADsEH1)
                                {
                                    if(x%n_evis_bins == i && y%n_evis_bins==j)
                                    {
                                        final_covmatrix[j+CombineMaxBins*i] = final_covmatrix[j+CombineMaxBins*i] + 1./ADsEH1*1./ADsEH1*TotalCovarianceMatrixM[y+MaxBins*x];
                                    }
                                }
                                else
                                {
                                    if(x%n_evis_bins == i && y%n_evis_bins==j)
                                    {
                                        final_covmatrix[j+CombineMaxBins*i] = final_covmatrix[j+CombineMaxBins*i] + 1./ADsEH2*1./ADsEH2*TotalCovarianceMatrixM[y+MaxBins*x];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    TH2D* PlotH = new TH2D("plot","plot",CombineMaxBins,0,CombineMaxBins,CombineMaxBins,0,CombineMaxBins);
    
    for (Int_t i = 0; i<CombineMaxBins; i++)
    {//columns
        
        for (Int_t j = 0; j<CombineMaxBins; j++)
        {//rows
            PlotH->SetBinContent(i+1,j+1,final_covmatrix[j+CombineMaxBins*i]);
        }
    }
    TCanvas* c1 = new TCanvas("c1");
    
    PlotH->Draw("colz");
    
    std::vector<TH1D*> Pred(ADsEH3*(ADsEH1+ADsEH2));
    std::vector<TH1D*> Far(ADsEH3*(ADsEH1+ADsEH2));
    std::vector<TH1D*> CombinedPred(ADsEH3*(ADsEH1+ADsEH2));
    std::vector<TH1D*> CombinedFar(ADsEH3*(ADsEH1+ADsEH2));
    
    TH1D* UnitH = new TH1D("unit","unit",n_evis_bins,0,n_evis_bins);
    for (Int_t i = 0; i<n_evis_bins; i++)
    {
        UnitH->SetBinContent(i+1,1);
    }
    for (Int_t far = 0; far < ADsEH3; far++)
    {
        for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            Pred[near+far*(ADsEH1+ADsEH2)]=(TH1D*)UnitH->Clone();
            
            if(near<ADsEH1)
            {
                //                Pred[near+far*(ADsEH1+ADsEH2)]->Add(UnitH);
            }
            
            Far[near+far*(ADsEH1+ADsEH2)]=(TH1D*)UnitH->Clone();
        }
    }
    
    CombinedPred = CombineData(Pred,CombinedPred);
    CombinedFar = CombineData(Far,CombinedFar);
    
    TCanvas* c2 = new TCanvas("c2");
    c2->Divide(3,3);
    Int_t counti = 1;
    for (Int_t far = 0; far < CombineMaxFar; far++)
    {
        for (Int_t near = 0; near < CombineMaxNear; near++)
        {
            c2->cd(counti++);
            CombinedPred[near+far*(ADsEH1+ADsEH2)]->Draw();
        }
    }
}

std::vector<TH1D*> Test :: CombineData(std::vector<TH1D*> Histo,std::vector<TH1D*> CombinedHisto)
{    //    Combine matrices in 9x9, 2x2 or 1x1 prediction
    
    Int_t CombineMaxFar, CombineMaxNear, CombineMaxBins;
    
    if (Combine==0)
    {
        CombineMaxFar = ADsEH3;
        CombineMaxNear = ADsEH1+ADsEH2;
        CombineMaxBins = MaxBins;
    }
    else if(Combine == 1)
    {
        CombineMaxFar = 1;
        CombineMaxNear = 1;
        CombineMaxBins=n_evis_bins;
    }
    else if(Combine == 2)
    {
        CombineMaxFar = 1;
        CombineMaxNear = 2;
        CombineMaxBins=2*n_evis_bins;
    }
    
    for (Int_t far = 0; far < CombineMaxFar; far++)
    {
        for (Int_t near = 0; near < CombineMaxNear; near++)
        {
            
            CombinedHisto[far+CombineMaxFar*near] = (TH1D*)Histo[far+CombineMaxFar*near]->Clone();
            
            if(Combine!=0)
            {
                CombinedHisto[far+CombineMaxFar*near]->Reset();
            }
        }
    }
    
    //Combine matrices in 9x9, 2x2 or 1x1 prediction
    if (Combine == 1)
    {
        for (Int_t far = 0; far < ADsEH3; far++)
        {
            for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
            {
                //1x1
                CombinedHisto[0+0]->Add(Histo[far+ADsEH3*near]);
            }
        }
        
        CombinedHisto[0+0]->Scale(1./(ADsEH1+ADsEH2));
        
    }
    else if(Combine == 2)
    {
        for (Int_t far = 0; far < ADsEH3; far++)
        {
            for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
            {
                //2x2
                
                if(near<ADsEH1)
                {
                    CombinedHisto[0+0]->Add(Histo[far+ADsEH3*near]);
                }
                else
                {
                    CombinedHisto[0+1*CombineMaxFar]->Add(Histo[far+ADsEH3*near]);
                }
            }
        }
        CombinedHisto[0+0]->Scale(1./ADsEH1);
        CombinedHisto[0+1*CombineMaxFar]->Scale(1./ADsEH2);
    }
    
    return CombinedHisto;
}

void Test :: TestCovarianceMatrixConstruction()
{
    TRandom3* rand = new TRandom3(0);
    
    TH2D* Cov2H = new TH2D("Cov","Cov",9*n_evis_bins,0,9*n_evis_bins,9*n_evis_bins,0,9*n_evis_bins);
    Double_t Cov[9*n_evis_bins*9*n_evis_bins];
    TH1D* AlteredPredictionVisH[(ADsEH1+ADsEH2)*ADsEH3];
    TH1D* PredictionVisH[(ADsEH1+ADsEH2)*ADsEH3];
    Double_t randvec[(ADsEH1+ADsEH2)*ADsEH3];
    
    Double_t Error = 0.5;//relative error
    Double_t Nominal = 1;
    
    for (Int_t near=0; near<(ADsEH1+ADsEH2); near++)
    {
        for (Int_t far=0; far<ADsEH3; far++)
        {
            AlteredPredictionVisH[near*ADsEH3+far] = new TH1D(Form("Altered Pred %d %d",near,far),"Altered Pred",n_evis_bins,0,n_evis_bins);
            
            PredictionVisH[near*ADsEH3+far] = new TH1D(Form("Pred %d %d",near,far),"Pred",n_evis_bins,0,n_evis_bins);
            
            for(Int_t i = 1; i<=n_evis_bins; i++)
            {
                AlteredPredictionVisH[near*ADsEH3+far]->SetBinContent(i,1);
                PredictionVisH[near*ADsEH3+far]->SetBinContent(i,1);
            }
        }
    }
    
    for (Int_t near=0; near<(ADsEH1+ADsEH2); near++)
    {
        for (Int_t far=0; far<ADsEH3; far++)
        {
            rand->SetSeed(0);
            randvec[near*ADsEH3+far] = Nominal*(1+Error*rand->Gaus());
            if(near < ADsEH1)
            {
                randvec[near*ADsEH3+far]=randvec[0];
            }
            AlteredPredictionVisH[near*ADsEH3+far]->Scale(randvec[near*ADsEH3+far]);
        }
    }
    
    Int_t x =0;
    Int_t y =0;
    
    for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
    {
        Int_t Ni1,Ni2,Ni3,Ni4;
        //Logic for the 2D matrix index done up to 8 ADs
        if(neari==0){Ni1=1; Ni2=0; Ni3=0; Ni4=0;}
        if(neari==1){Ni2++;}
        if(neari==2){Ni3++;}
        if(neari==3){Ni4++;}
        
        for (Int_t fari=0; fari<ADsEH3; fari++)
        {
            Int_t Fi1,Fi2,Fi3,Fi4;
            //Logic for the 2D matrix index done up to 8 ADs
            if(Ni1!=Ni2){Fi1=fari+1;Fi2 = 0;Fi3 = 0;Fi4 = 0;}
            if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3; Fi2=fari+1;}
            if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3; Fi3=fari+1;}
            if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3; Fi4=fari+1;}
            //                std::cout <<"\n"<< "Ni" << Ni1 << Ni2 << Ni3 << Ni4 << "\n";
            //                std::cout << "Fi"<< Fi1 << Fi2 << Fi3 << Fi4 << "\n";
            //                std::cout << "====================================" << icounter++ <<"\n";
            //                std::cout << Ni1*Fi1 << Ni2*Fi2 << Ni3*Fi3 << Ni4*Fi4;
            
            
            for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
            {
                Int_t Nj1,Nj2,Nj3,Nj4;
                //Logic for the 2D matrix index done up to 8 ADs
                if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                if(nearj==1){Nj2++;}
                if(nearj==2){Nj3++;}
                if(nearj==3){Nj4++;}
                
                for (Int_t farj=0; farj<ADsEH3; farj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    Int_t Fj1,Fj2,Fj3,Fj4;
                    if(Nj1!=Nj2){Fj1=farj+1;Fj2 = 0;Fj3 = 0;Fj4 = 0;}
                    if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3; Fj2=farj+1;}
                    if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3; Fj3=farj+1;}
                    if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3; Fj4=farj+1;}
                    //                        std::cout <<"\n"<< "Nj"<< Nj1 << Nj2 << Nj3 << Nj4 << "\n";
                    //                        std::cout << "Fj"<< Fj1 << Fj2 << Fj3 << Fj4 << "\n";
                    //                        std::cout << "===================================="<< jcounter++ <<"\n";
                    //                        std::cout << Nj1*Fj1 << Nj2*Fj2 << Nj3*Fj3 << Nj4*Fj4 <<"\n";
                    //                        std::cout <<"\n";
                    
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {//columns
                        x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                        
                        for (Int_t j = 0; j<n_evis_bins; j++)
                        {//rows
                            y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                            
                            //                            if(BackgroundE != (CovarianceMatrix3::BackgroundType)(-1))// Background has been varied
                            //                            {
                            //                                Cov[x*MaxBins+y]=(FarRandomBackgroundSpectrumH[fari]->GetBinContent(i+1)-FarBackgroundSpectrumH[fari]->GetBinContent(i+1))*(FarRandomBackgroundSpectrumH[farj]->GetBinContent(j+1)-FarBackgroundSpectrumH[farj]->GetBinContent(j+1));
                            //                            }
                            //                            else // Systematic has been varied, or nothing, in the latter case this should be 0. (Check)
                            //                            {
                            Cov[x*MaxBins+y]=(AlteredPredictionVisH[neari*ADsEH3+fari]->GetBinContent(i+1)-PredictionVisH[neari*ADsEH3+fari]->GetBinContent(i+1))*(AlteredPredictionVisH[nearj*ADsEH3+farj]->GetBinContent(j+1)-PredictionVisH[nearj*ADsEH3+farj]->GetBinContent(j+1));
                            //                            }
                            
                            
                            Cov2H->SetBinContent(x+1,y+1,Cov[x*MaxBins+y]);
                            
                        }
                    }
                }
            }
        }
    }
    
    TFile* TestF = new TFile("TestConstruction.root","recreate");
    for (Int_t near=0; near<(ADsEH1+ADsEH2); near++)
    {
        for (Int_t far=0; far<ADsEH3; far++)
        {
            AlteredPredictionVisH[near*ADsEH3+far]->Write();
            PredictionVisH[near*ADsEH3+far]->Write();
            Cov2H->Write();
        }
    }
    delete TestF;
}

void Test::TestInputData()
{
    Char_t NearData[100];
    Char_t NearDataSpec[100];
    
    if(P12E==0)//LBNL
    {
        sprintf(NearData,"./Inputs/GdInputs/ibd_eprompt_shapes.root");
    }
    else if(P12E==1)//P12E
    {
        sprintf(NearData,"./Inputs/GdInputs/file_P12E_nGd_IBD_Acc_spectrum.root");
    }
    
    std::vector<TH1D*> ADSpectrumVisH;
    std::vector<TH1D*> ADSpectrumH;
    
    TH2D* TransEnergyMatrix;
    TH2D* EnergyMatrixH;
    TH2D* InvEnergyMatrixH;
    
    ADSpectrumVisH.resize((ADsEH1+ADsEH2)*Nweeks);
    ADSpectrumH.resize((ADsEH1+ADsEH2)*Nweeks);
    
    TArrayD EnergyMatrixM;
    EnergyMatrixM.Set(MatrixBins*MatrixBins);
    
    Double_t InvEnergyMatrixM[MatrixBins][MatrixBins];
    
    InvEnergyMatrixH= new TH2D("Inv Energy Matrix","Inv Energy Matrix", MatrixBins,0,MatrixBins,MatrixBins,0,MatrixBins);
    
    TFile* TransEnergyMatrixDataF = new TFile("./ResponseMatrices/NominalResponseMatrix.root");
    
    TransEnergyMatrix = (TH2D*)gDirectory->Get("EnuEvis11");
    EnergyMatrixH = (TH2D*)gDirectory->Get("EvisEnu11");
    
    delete TransEnergyMatrixDataF;
    
    //    TMatrixD* InvertEnergyMatrix = new TMatrixD(MatrixBins,MatrixBins);
    //
    //    for(Int_t i = 0; i<MatrixBins; i++)
    //    {
    //        for(Int_t j = 0; j<MatrixBins; j++)
    //        {
    //            EnergyMatrixM[(j+MatrixBins*i)]=EnergyMatrixH->GetBinContent(i+1,j+1);
    //        }
    //    }
    //
    //    InvertEnergyMatrix->SetMatrixArray(EnergyMatrixM.GetArray());
    //
    //    InvertEnergyMatrix->Invert();
    //
    //    Double_t* InvertEnergyMatrixArray = InvertEnergyMatrix->GetMatrixArray();
    //
    //    for (Int_t i = 0; i < MatrixBins; i++)
    //    {
    //        for (Int_t j = 0; j < MatrixBins; j++)
    //        {
    //            InvEnergyMatrixM[i][j]=InvertEnergyMatrixArray[(i*MatrixBins+j)];
    //            InvEnergyMatrixH->SetBinContent(i+1,j+1,InvEnergyMatrixM[i][j]);
    //        }
    //    }
    //Rebin matrix to match near data:
    //
    //    TH2D* RebinInvEnergyMatrix = new TH2D("Rebin Inv Energy","Rebin Inv Energy",n_evis_bins,evis_bins,n_etrue_bins,enu_bins);
    //    TH2D* NormRebinInvEnergyMatrix = new TH2D("Norm Rebin Inv Energy","Norm Rebin Inv Energy",n_evis_bins,evis_bins,n_etrue_bins,enu_bins);
    //
    //    for (Int_t i = 0; i < n_evis_bins; i++)
    //    {
    //        for (Int_t VisibleEnergyIndex = Int_t(evis_bins[i]*MatrixBins/FinalVisibleEnergy)+1; VisibleEnergyIndex <= Int_t(evis_bins[i+1]*MatrixBins/FinalVisibleEnergy); VisibleEnergyIndex++)//Although the results are integers I need to add Int_t() if I don't want the precision errors to mess up the rebinning.
    //        {
    //            for (Int_t j = 0; j < n_etrue_bins; j++)
    //            {
    //                for (Int_t TrueEnergyIndex = Int_t(enu_bins[j]*MatrixBins/FinalVisibleEnergy)+1; TrueEnergyIndex <= Int_t(enu_bins[j+1]*MatrixBins/FinalVisibleEnergy); TrueEnergyIndex++)//Although the results are integers I need to add Int_t() if I don't want the precision errors to mess up the rebinning.
    //                {
    //
    //                    RebinInvEnergyMatrix->SetBinContent(j+1,i+1, RebinInvEnergyMatrix->GetBinContent(j+1,i+1)+TransEnergyMatrix->GetBinContent(TrueEnergyIndex,VisibleEnergyIndex));
    //
    //                }
    //            }
    //        }
    //
    //    }
    //
    //    //Normalize
    //    const Int_t LimitVis = n_evis_bins;
    //    const Int_t LimitTrue = n_etrue_bins;
    //
    //    Double_t NormaInv[LimitVis];
    //
    //    for(Int_t j=0;j<LimitVis;j++)
    //    {
    //        NormaInv[j]=0;
    //
    //        for(Int_t i=0;i<LimitTrue;i++)
    //        {
    //            NormaInv[j] = NormaInv[j]+RebinInvEnergyMatrix->GetBinContent(j+1,i+1);
    //        }
    //    }
    //
    //
    //
    //    for (Int_t i = 0; i < LimitVis; i++)
    //    {
    //        for (Int_t j = 0; j < LimitTrue; j++)
    //        {
    //
    ////            std::cout << NormaInv[i] << " " << i << std::endl;
    //
    //            if(NormaInv[i]!=0)
    //            {
    //                NormRebinInvEnergyMatrix->SetBinContent(i+1,j+1,RebinInvEnergyMatrix->GetBinContent(i+1,j+1)/NormaInv[i]);//Normalization so Σj E(i,j) = 1; (Σ(y axis) =1)
    //            }
    //        }
    //    }
    
    //Invert EnergyMatrix, rebin in rows and check if this produces a better spectrum
    
    for(Int_t week = 0; week < Nweeks; week++)
    {
        for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            if(Nweeks == 1)
            {
                if(P12E==0)//LBNL
                {
                    sprintf(NearDataSpec,"h_ibd_eprompt_inclusive_ad%i",near+1);
                }
                else if(P12E==1)//P12E
                {
                    if(near<ADsEH1)
                    {
                        sprintf(NearDataSpec,"hist_AccSub_%d_EH1",near+1);
                    }
                    if(near>ADsEH1)
                    {
                        sprintf(NearDataSpec,"hist_AccSub_%d_EH2",near-ADsEH1+1);
                    }
                }
            }
            else
            {
                if(P12E==0)//LBNL
                {
                    sprintf(NearDataSpec,"h_ibd_eprompt_week%d_ad%d", week, near+1);
                }
            }
            
            TFile* NearDataF = new TFile(NearData);
            
            ADSpectrumVisH[week+near*Nweeks] = (TH1D*)gDirectory->Get(NearDataSpec);
            
            delete NearDataF;
            
            //            for(Int_t j=1; j<=n_evis_bins; j++)
            //            {
            //                ADSpectrumVisH[week+near*Nweeks]->SetBinContent(j,1);
            //            }
            
            ADSpectrumH[week+near*Nweeks] = new TH1D("TrueEnergyData","TrueEnergyData",n_etrue_bins,enu_bins);
            
            for(Int_t i=1; i<=n_etrue_bins; i++)
            {
                for(Int_t j=1; j<=n_evis_bins; j++)
                {
                    //  Fluctuations in this one:
                    ADSpectrumH[week+near*Nweeks]->SetBinContent(i,ADSpectrumH[week+near*Nweeks]->GetBinContent(i)+TransEnergyMatrix->GetBinContent(j,i)*ADSpectrumVisH[week+near*Nweeks]->GetBinContent(j));
                }
            }
            ////////////////////////////////////////I need to multiply by the inverse of the corresponding AD response matrix here (Visible energy)
            //            ADSpectrumH[week+near*Nweeks]->Scale(1./(DetectorEfficiency[week+near*Nweeks]*LiveTime[week+near*Nweeks]*DetectorProtons[near]));// Remove AD effects
        }
    }
    
    
    TFile* DataAfterEvis = new TFile("DataAfterEvis.root","recreate");
    for(Int_t week = 0; week < Nweeks; week++)
    {
        TransEnergyMatrix->Write();
        EnergyMatrixH->Write();
        //        RebinInvEnergyMatrix->Write();
        //        NormRebinInvEnergyMatrix->Write();
        
        InvEnergyMatrixH->Write();
        delete TransEnergyMatrix;
        delete EnergyMatrixH;
        //        delete InvEnergyMatrixH;
        //        delete RebinInvEnergyMatrix;
        for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            ADSpectrumVisH[week+near*Nweeks]->Write();
            ADSpectrumH[week+near*Nweeks]->Write();
            delete ADSpectrumH[week+near*Nweeks];
            delete ADSpectrumVisH[week+near*Nweeks];
        }
    }
    DataAfterEvis->Close();
}

void Test :: TestRebinMatrix()
{
    TH1D* OscDeltaVisibleSpectrumH[MatrixBins];
    TH1D* OscDeltaVisibleSpectrumSumH[n_etrue_bins];
    TH1D* RebinnedOscDeltaVisibleSpectrumSumH[MatrixBins];
    
    for(Int_t TrueEnergyIndex = 0; TrueEnergyIndex<MatrixBins; TrueEnergyIndex++)
    {
        OscDeltaVisibleSpectrumH[TrueEnergyIndex]= new TH1D(Form("Fine Visible Spectrum Vis for True Energy Index %d",TrueEnergyIndex),Form("Fine Visible Spectrum Vis for True Energy Index %d",TrueEnergyIndex),MatrixBins,0,FinalVisibleEnergy);
        for(Int_t k = 0; k<MatrixBins; k++)
        {
            if(k==TrueEnergyIndex)
            {
                OscDeltaVisibleSpectrumH[TrueEnergyIndex]->SetBinContent(k+1,1);//Build diagonal test matrix.
                //                OscDeltaVisibleSpectrumH[TrueEnergyIndex]->SetBinContent(k+2,2);//Build diagonal test matrix.
            }
        }
    }
    
    std::cout << "\t Rebin Matrix from Fine to Coarse binning" << std::endl;
    
    Int_t LimitVis = n_evis_bins;
    Int_t LimitTrue = n_etrue_bins;
    
    TH2D* EnergyMatrixH = new TH2D("EvisEnu","EvisEnu",LimitTrue,enu_bins,LimitVis,evis_bins);
    
    TH2D* RowEnergyMatrixH = new TH2D("EnuEvis","EnuEvis", LimitTrue,enu_bins,LimitVis,evis_bins);
    
    for (Int_t i = 0; i < n_etrue_bins; i++)
    {
        OscDeltaVisibleSpectrumSumH[i]= new TH1D(Form("Fine Visible Spectrum Sum Vis for True Energy Index %d",i),Form("Fine Visible Spectrum Sum Vis for True Energy Index %d",i),MatrixBins,0,FinalVisibleEnergy);
        
        for (Int_t TrueEnergyIndex = Int_t(enu_bins[i]*MatrixBins/FinalVisibleEnergy); TrueEnergyIndex < Int_t(enu_bins[i+1]*MatrixBins/FinalVisibleEnergy); TrueEnergyIndex++)//Although the results are integers I need to add Int_t() if I don't want the precision errors to mess up the rebinning.
        {
            OscDeltaVisibleSpectrumSumH[i]->Add(OscDeltaVisibleSpectrumH[TrueEnergyIndex]);
        }
        
        RebinnedOscDeltaVisibleSpectrumSumH[i]=(TH1D*)OscDeltaVisibleSpectrumSumH[i]->Rebin(n_evis_bins,Form("Rebinned Visible Spectrum Vis for True Energy Index %d",i),evis_bins);
        
        for (Int_t j = 0; j < n_evis_bins; j++)
        {
            EnergyMatrixH->SetBinContent(i+1,j+1, RebinnedOscDeltaVisibleSpectrumSumH[i]->GetBinContent(j+1));
        }
    }
    for (Int_t i = 0; i < n_etrue_bins; i++)
    {
        for (Int_t j = 0; j < n_evis_bins; j++)
        {
            RowEnergyMatrixH->SetBinContent(i+1,j+1,EnergyMatrixH->GetBinContent(i+1,j+1));
        }
    }
    
    TFile* FileRebiNF = new TFile("TestRebin.root","recreate");
    EnergyMatrixH->Write();
    RowEnergyMatrixH->Write();
    
    for(Int_t TrueEnergyIndex = 0; TrueEnergyIndex<MatrixBins; TrueEnergyIndex++)
    {
        OscDeltaVisibleSpectrumH[TrueEnergyIndex]->Write();
    }
    for (Int_t i = 0; i < n_etrue_bins; i++)
    {
        OscDeltaVisibleSpectrumSumH[i]->Write();
    }
    
    delete FileRebiNF;
}


void Test :: TestTree()
{
    
    TH1D* PredictionH;
    Double_t sen22t13 = 0.09;
    Double_t dm2_31 = 0.0024;
    
    TTree *T;
    
    TFile *f = new TFile("./ToyMCTrees/ToyMCTreeCombined2.root");
    T = (TTree*)f->Get("TVar");
    
    TH1D* TNominalHisto = 0;
    TH1D* TDataHisto = 0;
    TH1D* TVarHisto = 0;
    
    Long_t entries= T->GetEntries();
    std::cout << entries << std::endl;
    
    Int_t far = 0;
    Int_t near = 0;
    Int_t MaxFarCombine = 1;
    Int_t MaxNearCombine = 2;
    
    Double_t s22t13start=0.0;
    Double_t s22t13end=0.2;
    Double_t dm2_31start=0.0015;
    Double_t dm2_31end=0.0035;
    
    Int_t Nsteps = 101;
    Int_t NstepsM = Nsteps;
    
    Double_t Sin22t13[Nsteps];
    Double_t Dm2_31[Nsteps];
    
    Double_t DeltaWidth =(dm2_31end-dm2_31start)*1./(NstepsM-1);
    Double_t SinWidth = (s22t13end-s22t13start)*1./(Nsteps-1);
    
    for(Int_t step=0;step < Nsteps;++step)
    {
        Sin22t13[step] = SinWidth*step + s22t13start;
        TCanvas *c1 = new TCanvas("c1","test",11000,11000);
        c1->Divide(11,11);
        
        for(Int_t stepM = 0; stepM < NstepsM ; stepM++)
        {
            Dm2_31[stepM] = DeltaWidth*stepM + dm2_31start;
            
            //                    T->GetEntry(step*NstepsM+stepM);//Calculate entry number for choosen sin22t13 and dm2_31. The tree contains 100*100 grid points, with MaxFarCombine*MaxNearCombine predictions in each point.
            c1->cd(1+stepM);
            T->Draw("NominalHistoDayaBay.Draw()","","goff",1,(step*NstepsM+stepM));
            c1->Update();
            
            //            c1->cd(2+3*(step*NstepsM+stepM));
            //            T->Draw("DataHistoDayaBay.Draw()","","goff",1,step*NstepsM+stepM);
            //            c1->cd(3+3*(step*NstepsM+stepM));
            //            T->Draw("VariationHistoDayaBay.Draw()","","goff",1,step*NstepsM+stepM);
            
            //            std::cout << step*NstepsM+stepM << std::endl;
            
        }
        c1->Print(Form("./Images/TreeSamples/Tree_%f.eps",Sin22t13[step]),".eps");
        
        delete c1;
        
    }
    
    delete f;
    
}

void Test :: CompareWithHenoch()
{
    TFile* HenochPredictionsF = new TFile("./Inputs/GdInputs/IHEP_data_lbnlbin_6AD.root");
    
    std::cout << "TEST HENOCH INPUTS " << std::endl;
    
    
    std::string lines;
    
    Double_t dNdE[NADs];
    TH1D* NominalSpectra[NADs];
    
    for (Int_t AD = 0; AD < NADs; AD++)
    {
        NominalSpectra[AD] = new TH1D("Spectrum","Spectrum",240,0,240);
        dNdE[AD]=0;
    }
    
    Int_t curSample = 0;
    Double_t eNu;
    
    std::string line;
    
    std::ifstream fileData("./Test/lbnl_toyspectra_p12e_blinded/toySpectra_nominal_p12e_blinded_s2t_0.100_dm2_0.00235.txt");

    while(!fileData.eof())
    {
        getline(fileData,line);
        std::string firstchar = line.substr(0,1);
        
        if(firstchar=="#") continue;//<-- ignore lines with comments
        
        std::istringstream iss(line);
        
        Int_t column=0;
        while(iss)
        {
            std::string sub; iss >> sub;
            if(column==0) eNu=atof(sub.c_str());
            
            std::cout << "energy " << eNu << std::endl;
            
            if(column>1 && sub!="")
            {
                dNdE[column-1]=atof(sub.c_str());
                NominalSpectra[column-1]->SetBinContent(curSample+1,dNdE[column-1]);
                std::cout << "spectra " << dNdE[column-1] << std::endl;

            }
            
            column+=1;
        }//looping over columns
        std::cout << " sample " << curSample << std::endl;

        curSample++;
    }

    fileData.close();
    
    std::cout << "TOTAL SAMPLES: " << curSample-- << std::endl;
    
    for (Int_t AD = 0; AD < NADs; AD++)
    {
        NominalSpectra[AD]->Rebin(n_evis_bins,Form("SpectrumAD%d",AD),evis_bins);
    }
    
    char *dirname;
    char *ext= Form(".txt");
   
    dirname = "/Users/royal/DayaBay/OthersDYBCode/Fits/LBNLShapeFit/covariance_matrices_unified_p12e_unblinded/";

    Int_t num = 0;
    TString fname[30];
    
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    if (files)
    {
        TSystemFile *file;
        TIter next(files);
        
        while ((file=(TSystemFile*)next()))
        {
            fname[num] = file->GetName();
            
            if (!file->IsDirectory() && fname[num].EndsWith(ext))
            {
                std::string line2;
                Int_t number_of_lines;
                Int_t spaces;
                std::ifstream fileCov(fname[num]);
                Int_t x =0;
                Int_t y =0;
                Double_t value;
                
                TH2D* MatrixH = new TH2D("MatrixH","MatrixH",n_evis_bins*9,0,n_evis_bins*9,n_evis_bins*9,0,n_evis_bins*9);
                
                while (std::getline(fileCov, line2))
                {
                    std::string firstchar = line2.substr(0,1);
                    if(line2.empty())
                    {
                        spaces++;
                        y++;
                        x=0;
                        
                        continue;//<-- ignore lines with comments
                    }
                    
                    value=atof(line2.c_str());
                    
                    std::cout << "ATOF OF 0 " << atof("");
                    
                    std::cout << "COV MATRIX VALUE : " << value << " x: " << x << " y: " << y << std::endl;
                    MatrixH->SetBinContent(x+1,y+1, value);
                    
                    x++;
                    ++number_of_lines;
                }
                
                std::cout << "Number of lines in text file: " << number_of_lines << " - spaces :" << spaces << std::endl;
                
                TCanvas* CovMatrixC = new TCanvas("CovMatrix","CovMatrix");
                MatrixH->Draw("colz");
//                CovMatrixC->Print(("./Images/Test/"+fname[num].operator(0,fname[num].length() - 4)+".eps"));
                delete CovMatrixC;
                
                num++;
            }
        }
    }
}

//Input a flat vector to the fitter, use a diagonal covariance matrix as well and check that the fitter works for this simple case
void Test :: VectorOnes()
{
    TH1D* OnesH[2][1];

    for(Int_t near = 0; near < 2; near++)
    {
        for(Int_t far = 0; far < 1; far++)
        {
            OnesH[near][far] = new TH1D(Form("OnesHFar%d_Near%d",far,near),Form("OnesHFar%d_Near%d",far,near),n_evis_bins,evis_bins);
            
            for(Int_t i = 0; i<2*n_evis_bins; i++)
            {
                OnesH[near][far]->SetBinContent(i+1,1);
            }
        }
    }
    
    TTree* OnesT = new TTree("TVar","TVar");
    
    OnesT->Branch("NominalHistoDayaBay", "TH1D", &OnesH[0][0]);
    OnesT->Branch("NominalHistoLingAo","TH1D", &OnesH[1][0]);
    OnesT->Branch("DataHistoDayaBay","TH1D", &OnesH[0][0]);
    OnesT->Branch("DataHistoLingAo","TH1D", &OnesH[1][0]);
    
    for(Int_t sin22t13 = 0; sin22t13 <101; sin22t13++)
    {
        for(Int_t dm2_13 = 0; dm2_13 <101; dm2_13++)
        {
            for(Int_t near = 0; near < 2; near++)
            {
                for(Int_t far = 0; far < 1; far++)
                {
                    OnesH[near][far]->Scale(1.*(sin22t13+1)*(dm2_13+1));
                }
            }
            OnesT->Fill();
            for(Int_t near = 0; near < 2; near++)
            {
                for(Int_t far = 0; far < 1; far++)
                {
                    OnesH[near][far]->Scale(1./(sin22t13+1)*(dm2_13+1));
                }
            }
        }
    }
    
    TFile* TreeOnesF = new TFile("./ToyMCTrees/TreeOnesCombine2.root","recreate");
    
    OnesT->Write();
    
    delete TreeOnesF;
}

void Test :: LUDecomp()
{
//    Double_t Example[3][3];
//    Example[0][0]=1;
//    Example[1][0]=2;
//    Example[2][0]=3;
//    Example[0][1]=4;
//    Example[1][1]=5;
//    Example[2][1]=7;
//    Example[0][2]=6;
//    Example[1][2]=8;
//    Example[2][2]=9;
    
    Double_t Example[9];
    Example[0]=25;
    Example[1]=15;
    Example[2]=-5;
    Example[3]=15;
    Example[4]=18;
    Example[5]=0;
    Example[6]=-5;
    Example[7]=0;
    Example[8]=11;
    
    Int_t MaxMatrix = 3;
    
    TMatrixD renormmatrix(3,3,&Example[0]);
    
    TDecompChol cholrenorm(renormmatrix);//  M = L*U
    cholrenorm.Decompose();
    TMatrixD cmat(cholrenorm.GetU());//   U
    TMatrixD tcmat(cmat.Transpose(cmat));// L
    
    Double_t* tmp_matrix = tcmat.GetMatrixArray();
    
    TH2D* TestM = new TH2D("","",MaxMatrix,0,MaxMatrix,MaxMatrix,0,MaxMatrix);
    gStyle->SetOptStat(kFALSE);
    
    for (Int_t i = 0; i < MaxMatrix; i++)
    {
        for (Int_t j = 0; j < MaxMatrix; j++)
        {
            TestM->SetBinContent(i+1,j+1,tmp_matrix[i*MaxMatrix+j]);
        }
    }
    TCanvas* c1 = new TCanvas("c1","c1");
    c1->Divide(3,1);
    c1->cd(1);
        renormmatrix.Draw("colz");
    c1->cd(2);
    cmat.Draw("colz");
    c1->cd(3);
    TestM->Draw("colz");
    c1->Print("./Images/Test/TestLU.eps");
    delete c1;
}

void Test :: LBNLPredictions()
{
    
    TFile* LoadPredictionF;
    TFile* LoadPredictionF1;

    TH1D* LBNLHisto[NADs];
    TH1D* MyHisto[NADs];
    TH1D* AccidentalsH[NADs];
    TH1D* LiHeH[NADs];
    TH1D* FastNeutronsH[NADs];
    TH1D* AmCH[NADs];
    TH1D* BackgroundSpectrumH[NADs];
    LoadPredictionF = new TFile("/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/ToyMCTrees/nGdInputs/toySpectra_allsys_and_stat.root");
    
    //Get random sample, need to see how the tree is built, this is nominal sample so far.
    for (Int_t AD = 0; AD < NADs; AD++)
    {
        LBNLHisto[AD]=(TH1D*)gDirectory->Get(Form("h_nominal_ad%i",AD+1));//For 6ADs (Yasu's inputs)
    }
    
    delete LoadPredictionF;
    
    TFile* BackgroundsF = new TFile("./BackgroundSpectrum/GDBackground/Backgrounds.root");
    for (Int_t AD = 0; AD < NADs; AD++)
    {
        AccidentalsH[AD]=(TH1D*)gDirectory->Get(Form("Accidentals_AD%i",AD));
        LiHeH[AD]=(TH1D*)gDirectory->Get(Form("LiHe_AD%i",AD));//Missing LiHe inputs so far in Hydrogen Analysis
        FastNeutronsH[AD]=(TH1D*)gDirectory->Get(Form("FN_AD%i",AD));
        AmCH[AD]=(TH1D*)gDirectory->Get(Form("AmC_AD%i",AD));
        
        //Nominal background spectrum
        BackgroundSpectrumH[AD] = new TH1D(Form("Background Spectrum_%d",AD),Form("Background Spectrum_%d",AD),n_evis_bins,evis_bins);
        BackgroundSpectrumH[AD]->Add(AccidentalsH[AD]);
        BackgroundSpectrumH[AD]->Add(LiHeH[AD]);
        BackgroundSpectrumH[AD]->Add(FastNeutronsH[AD]);
        BackgroundSpectrumH[AD]->Add(AmCH[AD]);
        
        //Substract Backgrounds

        LBNLHisto[AD]->Add(BackgroundSpectrumH[AD],-1);
        
        std::cout << " LBNLHisto HISTO AD " << AD << " is " << LBNLHisto[AD]->Integral() << std::endl;

        LBNLHisto[AD]->Scale(1./LBNLHisto[AD]->Integral());
        delete BackgroundSpectrumH[AD];
    }
    delete BackgroundsF;
    
    LoadPredictionF1 = new TFile("/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/RootOutputs/Gadolinium/NominalOutputs/Oscillation.root");

    LoadPredictionF1->cd("Total AD Spectra after oscillation");
    
    //Get random sample, need to see how the tree is built, this is nominal sample so far.
    for (Int_t AD = 0; AD < NADs; AD++)
    {
        MyHisto[AD]= (TH1D*)gDirectory->Get(Form("Oscillation Prediction AD%i, week0",AD+1));
        
        std::cout << " MY HISTO AD " << AD << " is " << MyHisto[AD]->Integral() << std::endl;
        
        MyHisto[AD]->Scale(1./MyHisto[AD]->Integral());
    }
    
    TCanvas* ShapeDifferenceC = new TCanvas("","");
    
    ShapeDifferenceC->Divide(NADs/2,2);
    
    for (Int_t AD = 0; AD < NADs; AD++)
    {
        LBNLHisto[AD]->Add(MyHisto[AD],-1);
        LBNLHisto[AD]->Divide(MyHisto[AD]);
        
        ShapeDifferenceC->cd(AD+1);
        
        LBNLHisto[AD]->Draw("hist");
    }
    
    ShapeDifferenceC->Print("./Images/DifferenceBetweenLBNLpredictiosnAndMine.eps");
    delete ShapeDifferenceC;
    
    delete LoadPredictionF1;

    
}

void Test::FitterMainTests()
{
    //test to check the process of creating tree and recovering tree samples properly
    
    Int_t Combine=1;
    Double_t Tsin22t13;
    Double_t Tdm2_ee;
    Int_t Experiment;

    Int_t TNominalHisto;
    Int_t TNominalHisto1;
    
    TTree* T;
    Int_t FakeExperiments = 100;
    Int_t GridSteps = 21;
    Int_t steps = 6;
    Double_t dm2_eeend = 10;
    Double_t dm2_eestart = 5;
    Double_t s22t13end = 5;
    Double_t s22t13start = 0;
    
    TTree *TTest= new TTree("TTest","TestExperiments");
    
    TTest->Branch("VariationHistoDayaBay",&TNominalHisto,"TNominalHisto/I");
    TTest->Branch("VariationHistoLingAo",&TNominalHisto1,"TNominalHisto1/I");

    TTest->Branch("sin22t13",&Tsin22t13,"sin22t13/D");
    TTest->Branch("dm2_ee",&Tdm2_ee,"dm2_ee/D");
    TTest->Branch("Experiment",&Experiment,"Experiment/I");
    
    
    Double_t DeltaWidth = (dm2_eeend-dm2_eestart)*1./(steps-1);
    Double_t SinWidth = (s22t13end-s22t13start)*1./(steps-1);
    
    for(Int_t step=0;step < steps;step++)
    {
        Tsin22t13 = SinWidth*step + s22t13start;
        
        for(Int_t stepM = 0; stepM < steps ; stepM++)
        {
            Tdm2_ee = DeltaWidth*stepM + dm2_eestart;
            
            for(Int_t i =0;i<FakeExperiments;i++)
            {
                Experiment = i;
                
                for (Int_t far = 0; far<1; far++)
                {
                    for (Int_t near = 0; near<2; near++)
                    {
                        if(near==0)
                        {
                            TNominalHisto = 100+Experiment;
                        }
                        else
                        {
                            TNominalHisto1 = 200+Experiment;
                        }
                    }
                }
                
                TTest->Fill();
            }
        }
    }

    TFile* f = new TFile("./ToyMCTrees/TestTree.root","recreate");
    
    TTest->Write();//Save file for additional checks;
    
    delete f;

    TFile* f1 = new TFile("./ToyMCTrees/TestTree.root");

    //Proceed to recovery
    T = (TTree*)f1->Get("TTest");
    
    Double_t sin22t13,dm2_ee;
    Int_t experiment;

    T->SetBranchAddress("sin22t13",&Tsin22t13);
    T->SetBranchAddress("dm2_ee",&Tdm2_ee);

    T->SetBranchAddress("VariationHistoDayaBay",&TNominalHisto);
    T->SetBranchAddress("VariationHistoLingAo",&TNominalHisto1);
    
    T->SetBranchAddress("Experiment",&experiment);

    for(Int_t step=0;step < steps;step++)
    {
        sin22t13 = SinWidth*step + s22t13start;
        
        for(Int_t stepM = 0; stepM < steps ; stepM++)
        {
            dm2_ee = DeltaWidth*stepM + dm2_eestart;
            
            for(Int_t i =0;i<FakeExperiments;i++)
            {
                Experiment = i;
                
                T->GetEntry(Experiment+FakeExperiments*((steps*(sin22t13-s22t13start)/SinWidth)+((dm2_ee-dm2_eestart)/(DeltaWidth))));//Calculate entry number for choosen sin22t13 and dm2_ee. The tree contains 100*100 grid points, with MaxFarCombine*MaxNearCombine predictions in each point.
           
                std::cout << "Data TREE ENTRY: " << experiment << " " << experiment+FakeExperiments*(steps*(Tsin22t13/SinWidth)+((Tdm2_ee-dm2_eestart)/(DeltaWidth))) << " " << (Tdm2_ee-dm2_eestart)/(DeltaWidth) << " " << (Tsin22t13/SinWidth) << std::endl;
                
                std::cout << "Data TREE ENTRY should be: " << Experiment << " " << Experiment+FakeExperiments*(steps*(sin22t13/SinWidth)+((dm2_ee-dm2_eestart)/(DeltaWidth))) << " " << (dm2_ee-dm2_eestart)/(DeltaWidth) << " " << (sin22t13/SinWidth) << std::endl;

            }
        }
    }
    delete f1;
}

void Test::TestAllPredictions(bool Analysis)
{
    //Antineutrino
    
    //Reactor
    
    //Reactor Prediction
    
    bool mode = 0;//Nominal
    Int_t week = 0;
    bool ToyMC = 1;
    bool CovMatrix=1;
    
    NominalData* Data = new NominalData(Analysis,2);//Gadollinium
    
    Data->SetCombineMode(0); //0 is 9x9, 1 is 1x1 and 2 is 2x2f

    //    Data->SetSin22t12(0);
    //Parameters of the model
    Data->SetAnalysis(Analysis);//  Gd or H data
    Data->SetBinning(0);//  0 for LBNL binning or 1 for Linear binning
    Data->SetWeeks(1); //   1 or 20 so far
    
    Data->SetUnifiedModel(1);
    
    Data->SetToyMC(ToyMC);
    
    Data->SetAllRandomSystematics(0);
    Data->SetStatisticalFluctuation(0);
    
    Double_t sin22t13=0.09;
    Double_t dm2_32 = 2.41e-3;//eV2
    
    Int_t hierarchy=1;//-1 for inverted           // Set as external input to be set in the GUI   !!!!
    
    Double_t dm2_ee = dm2_32+hierarchy*5.21e-5;
    
    Prediction* Pred = new Prediction(Data);
    Pred->MakePrediction(sin22t13, dm2_ee,  mode,  week,  ToyMC, CovMatrix);
    
    TH1D* ReactorPred[NADs];
    TH1D* Prediction[MaxFarADs][MaxNearADs];
    TH1D* DifferencePred[MaxFarADs][MaxNearADs];
    
    //Difference between reactorpred and reactor should be the background spectrum:
    //compare to the one produced in makepred to see if the background is correctly added

    TCanvas* CheckReactorPred = new TCanvas("2","2");

    CheckReactorPred->Divide(NADs,2);
    
    for (Int_t AD=0; AD<NADs; AD++)
    {
        ReactorPred[AD] = Pred->GetReactorPrediction(AD);
        
        ReactorPred[AD]->SetStats(kTRUE);

        CheckReactorPred->cd(AD+NADs+1);
        ReactorPred[AD]->Draw("HIST");
    }

    CheckReactorPred->Print("./Images/Test/ReactorPredictionWithAndWithoutBkg.eps");
    
    delete CheckReactorPred;
    
    //Combined-Reactor prediction difference
    //combine=0 is default;

    TCanvas* CheckPredDifferenceC = new TCanvas("11","11");
    TCanvas* CheckPredDifferenceC2 = new TCanvas("22","22");
    CheckPredDifferenceC->Divide(NADs/2,NADs/2);
    CheckPredDifferenceC2->Divide(NADs/2,2);
    for (Int_t far=0; far<NADs/2; far++)
    {
        for (Int_t near=0; near<NADs/2; near++)
        {
            Prediction[far][near]= Pred->GetPrediction(far,near);//Has backgrounds included in oscillation.h
            Prediction[far][near]->SetStats(kTRUE);

            DifferencePred[far][near] = (TH1D*)Prediction[far][near]->Clone();
            DifferencePred[far][near]->Add(ReactorPred[far+NADs/2],-1);
            
            CheckPredDifferenceC->cd(near+(far*NADs/2)+1);
            DifferencePred[far][near]->Draw("HIST");

            CheckPredDifferenceC2->cd(far+1);
            Prediction[far][near]->Draw("HIST");

            CheckPredDifferenceC2->cd(far+NADs/2+1);
            ReactorPred[far+NADs/2]->Draw("HIST");

        }
    }
    CheckPredDifferenceC->Print("./Images/Test/PredictionsDifferenceWithBackground.eps");
    CheckPredDifferenceC2->Print("./Images/Test/ReactorAndRelativePredictionsDifference.eps");

    delete CheckPredDifferenceC;
    delete CheckPredDifferenceC2;
    //Combined Prediction
    
    
    
    delete Pred;
}

void Test::TestSuperHistograms()
{
    TFile* LBNLSuperHistoF = new TFile("/Users/royal/DayaBay/OthersDYBCode/Berkeley/Theta13Analysis2/ShapeFit/Flux/SuperHistograms_20week.root");
    
    TH1F* LBNLSuperHisto[6][6];
    TH1F* LBNLSuperHistoRatio[6][6];

    for(Int_t reactor = 0; reactor<6; reactor++)
    {
        for(Int_t ad = 0; ad<6; ad++)
        {
            LBNLSuperHisto[reactor][ad]=(TH1F*)gDirectory->Get(Form("h_super_idet%d_icore%d",ad,reactor));
            LBNLSuperHistoRatio[reactor][ad]=(TH1F*)LBNLSuperHisto[reactor][ad]->Clone();
        }
    }
    
    for(Int_t reactor = 0; reactor<6; reactor++)
    {
        for(Int_t ad = 0; ad<6; ad++)
        {
            LBNLSuperHistoRatio[reactor][ad]->Divide(LBNLSuperHisto[0][0]);
        }
    }
    
    delete LBNLSuperHistoF;
    
    TFile* MySuperHistoF = new TFile("/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/Inputs/Gadolinium/SuperFlux20.root");
    
    TH1F* MySuperHisto[6][6];
    TH1F* MySuperHistoRatio[6][6];
    
    for(Int_t reactor = 0; reactor<6; reactor++)
    {
        for(Int_t ad = 0; ad<NADs; ad++)
        {
            MySuperHisto[reactor][ad]=(TH1F*)gDirectory->Get(Form("Rebinned Flux Spectrum from reactor%d with ad%d efficiencies",reactor+1,ad+1));
            MySuperHistoRatio[reactor][ad]=(TH1F*)MySuperHisto[reactor][ad]->Clone();
        }
    }
    
    for(Int_t reactor = 0; reactor<6; reactor++)
    {
        for(Int_t ad = 0; ad<6; ad++)
        {
            MySuperHistoRatio[reactor][ad]->Divide(MySuperHisto[0][0]);
        }
    }
    
    delete MySuperHistoF;
    
    TCanvas* SuperHistoResultsC = new TCanvas("","");
    TCanvas* SuperHistoResultsC2 = new TCanvas("2","2");

    SuperHistoResultsC->Divide(6,2);
    SuperHistoResultsC2->Divide(6,2);

        for(Int_t ad = 0; ad<NADs; ad++)
        {
            SuperHistoResultsC->cd(ad+1);

            MySuperHistoRatio[0][ad]->Draw();
            
            SuperHistoResultsC2->cd(ad+1);
            MySuperHisto[0][ad]->Draw();

            SuperHistoResultsC->cd(ad+1+NADs);
            LBNLSuperHistoRatio[0][ad]->Draw();
            
            SuperHistoResultsC2->cd(ad+1+NADs);
            LBNLSuperHisto[0][ad]->Draw();
        }
    SuperHistoResultsC->Print("./Images/Test/SuperHistogramTest.eps");
    
    SuperHistoResultsC2->Print("./Images/Test/SuperHistogramTest2.eps");
}//Results to be understood, the comparison between superhistograms from LBNL are flat, meaning that they are the same for all ADs, but different for each reactor. In my case since AD efficiencies are applied all 36 histograms are different

void Test :: TestCovarianceMatrixSamples()//A test to decide how many samples are necessary to build the covariance matrices.
{
    bool Generate = 0;//1 to generate 500 covariance matrices in a file, 0 to do the test

    if(TestExternalInputs)
    {
        ToyMCSampleDirectory = "/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/ToyMCTrees/nGdInputs/toySpectra_allsys_and_stat.root";
        NominalPredictionsDirectory = "/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/ToyMCTrees/nGdInputs/toySpectra_allsys_and_stat.root";
        ResponseMatrixDirectory = "/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/matrix_evis_to_enu_unified_p12e_unblinded.root";
        SysCovarianceMatrixDirectory = "/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/CovarianceMatrices/Berkeley/matrix_sigsys_full.root";
        BkgCovarianceMatrixDirectory = "/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/CovarianceMatrices/Berkeley/matrix_bgsys_full.root";
    }
    else
    {
        ToyMCSampleDirectory =  "";
        NominalPredictionsDirectory = "";
        ResponseMatrixDirectory = "";
        SysCovarianceMatrixDirectory = "";
        BkgCovarianceMatrixDirectory = "";
    }
    
    Minuit = 0;//Choose between applying Minuit or manual grid fitting.
    Fit2D = 0; // 1 for 2D Fit, 0 for 1D Fit
    FitSin22t13 = 1; //  1 for Sin22t13 Fit, 0 for DM2ee Fit.
    NFits = 101;//101 in the final version
    Period=1;
    NReactorPeriods=20;
    PlotBin=0;//Initial bin
    Binning = 0;
    NADs=6;
    ToyMC=1;//1 ToyMC, 0 Data. To test the fitter use 1, to fit real data use 0.
    deleteFlag=0;
    deleteFlagSpec=0;
    NearTrueIndex = 0;
    flagCombine = 1;//Show Combine Plot as default
    Analysis = 0; //0 for Gd, 1 for H

    NSamples = 500; //change default to 500 once debug is completed
    CombineMode = 2; // change default to 2 once debug is completed
    NL[0]=0;//BCW
    NL[1]=0;//LBNL
    NL[2]=1;//Unified (Default)

    //         EnergyOffsetMatrix=0;
    //         AbsoluteScaleMatrix=0;
    //         AbsoluteOffsetMatrix=0;
    IAVMatrix=0;
    NLMatrix=0;
    ResolutionMatrix=0;
    Sin22t12Matrix=0;
    EfficiencyMatrix=0;
    
    UseToyMCTree = 0;
    PlotCovariance= 1;
    
    //To draw using a better palette:
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    
    Int_t DataSet=2;//0 is Simulation, 2 is P12E
    
    NominalData* Data = new NominalData(Analysis,DataSet);
    
    Data->SetToyMCSamplesDirectory(ToyMCSampleDirectory);
    Data->SetPredictionDirectory(NominalPredictionsDirectory);
    Data->SetResponseDirectory(ResponseMatrixDirectory);
    Data->SetBkgCovDirectory(BkgCovarianceMatrixDirectory);
    Data->SetSysCovDirectory(SysCovarianceMatrixDirectory);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // CrossCalc takes a lot of time and it needs to be run only once to produce the root file. Uncomment only if binning is different from standards and you want to produce a new cross section for this new binning.
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(Analysis)//nGd cross section is already tabulated, this is to calculate nH
    {
        CrossSection* Cross = new CrossSection(Data);
        Cross->CrossCalc();
        delete Cross;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Parameters of the model
    
    Data->SetNSamples(NSamples);// 500 in the final version.
    
    Data->SetToyMC(1);//  1 for Toy MC, 0 for data. To produce covariance matrices we use ToyMC.
    
    Data->SetCombineMode(CombineMode); //0 is 9x9, 1 is 1x1 and 2 is 2x2
    Data->SetUseToyMCTree(UseToyMCTree);
    Data->SetBinning(Binning);//  0 for LBNL binning or 1 for Linear binning
    Data->SetWeeks(Period);
    Data->SetNReactorPeriods(NReactorPeriods);
    std::cout <<  "********************************************************************************************************" << std::endl;
    std::cout <<  "**********************************            MAIN             *****************************************" << std::endl;
    std::cout <<  "********************************************************************************************************" << std::endl;
    
    //Chose Data Set
    if(Data->GetAnalysis())//   Hydrogen data
    {
        switch (DataSet)
        {
            case 0://  Simple reactor model used as input data
                std::cout << "\t Loading simple reactor model" << std::endl;
                break;
            case 1://   P12E
                if(1==Data->GetWeeks())
                {
                    std::cout << "\t Loading nH P12E Data" << std::endl;
                    Data->LoadMainData("./Inputs/HInputs/P12E_Inclusive.txt");
                }
                else
                {
                    std::cout << "\t \t \t NO MULTIPLE WEEK P12E DATA IN H ANALYSIS YET " << std::endl;
                    exit(EXIT_FAILURE);
                    Data->LoadMainData(Form("./Inputs/HInputs/P12E_%d.txt",NReactorPeriods));
                }
                break;
            case 2:// LBNL
                std::cout << "\t \t \t NO LBNL H ANALYSIS, OPTION NOT VALID " << std::endl;
                exit(EXIT_FAILURE);
                break;
            default:
                break;
        }
    }
    else//  Gd data
    {
        switch (DataSet)
        {
            case 0://  Simple reactor model used as input data
                std::cout << "\t Loading simple reactor model" << std::endl;
                break;
            case 1://   P12E
                if(1==Data->GetWeeks())
                {
                    std::cout << "\t Loading Gd P12E Data" << std::endl;
                    //                    Data->LoadMainData("./Inputs/GdInputs/P12E_Inclusive.txt");
                }
                else
                {
                    std::cout << "\t \t \t NO MULTIPLE WEEK P12E DATA IN Gd ANALYSIS YET " << std::endl;
                    exit(EXIT_FAILURE);
                    Data->LoadMainData(Form("./Inputs/GdInputs/P12E_%d.txt",NReactorPeriods));
                }
                break;
            case 2:// LBNL
                if(1==Data->GetWeeks())
                {
                    std::cout << "\t Loading LBNL Gd Data" << std::endl;
                    //                    Data->LoadMainData("./Inputs/GdInputs/Theta13-inputs_20week_inclusive.txt");
                }
                else
                {
                    std::cout << "\t Loading weekly LBNL Gd Data" << std::endl;
                    Data->LoadMainData(Form("./Inputs/GdInputs/Theta13-inputs_%dweek.txt",NReactorPeriods));
                }
                break;
            default:
                break;
        }
    }
    
    //To check non theta13 oscillation behaviour of the fitter. To avoid any oscillation at all then theta12 has to be set to 0 too.
    //    Data->SetSin22t12(0);
    //    Data->SetSin22t13(0);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    CovarianceMatrix3* Cov;
    
    Data->SetBCWModel(NL[0]);
    Data->SetLBNLModel(NL[1]);
    Data->SetUnifiedModel(NL[2]);
    //Order is important (Inititialize NL model before randomize it)
    VaryAccidentalMatrix=0;
    VaryLiHeMatrix=0;
    VaryFastNeutronsMatrix=0;
    VaryAmCMatrix=0;
    DistortLiHeMatrix=0;
    DistortFastNeutronsMatrix=0;
    DistortAmCMatrix=0;
    //Systematics
    IsotopeMatrix=0;
    ReactorPowerMatrix=0;
    EnergyScaleMatrix=1;
    
    Data->SetVaryAccidentalMatrix(VaryAccidentalMatrix);
    Data->SetVaryLiHeMatrix(VaryLiHeMatrix);
    Data->SetVaryFastNeutronsMatrix(VaryFastNeutronsMatrix);
    Data->SetVaryAmCMatrix(VaryAmCMatrix);
    Data->SetDistortLiHeMatrix(DistortLiHeMatrix);
    Data->SetDistortFastNeutronsMatrix(DistortFastNeutronsMatrix);
    Data->SetDistortAmCMatrix(DistortAmCMatrix);
    Data->SetIsotopeMatrix(IsotopeMatrix);
    Data->SetReactorPowerMatrix(ReactorPowerMatrix);
    Data->SetRelativeEnergyScaleMatrix(EnergyScaleMatrix);
    //        Data->SetAbsoluteEnergyScaleMatrix(AbsoluteScaleMatrix);
    //        Data->SetRelativeEnergyOffsetMatrix(EnergyOffsetMatrix);
    //        Data->SetAbsoluteEnergyOffsetMatrix(AbsoluteOffsetMatrix);
    Data->SetIAVMatrix(IAVMatrix);
    Data->SetNLMatrix(NLMatrix);
    Data->SetResolutionMatrix(ResolutionMatrix);
    Data->SetSin22t12Matrix(Sin22t12Matrix);
    Data->SetEfficiencyMatrix(EfficiencyMatrix);
   
    Int_t SampleStep = 20;
    Int_t MaxRepetitions = 15;
    if(Generate)//GenerateFiles
    {
        for(Int_t repetitions = 0; repetitions<MaxRepetitions; repetitions++)
        {
            std::cout << "REP: " << repetitions << std::endl;
            for(Int_t samples = NSamples-(SampleStep*4)+1; samples<NSamples+1-(SampleStep); samples+=SampleStep)
            {
                std::cout << "CALCULATING SAMPLE: " << samples << std::endl;
                
                Data->SetNSamples(samples);// 500 in the final version.
                
                Cov = new CovarianceMatrix3(Data);
                
                Cov->CovarianceMatrixMain(NULL);
                
                delete Cov;
            }
        }
    }
    else
    {
        TFile* File = new TFile("./Test/CovarianceMatrices/Gadolinium/Combine2/CovarianceMatricesRoot/RelativeEnergyScale.root");
        
        TH2D* CovMatrix500[MaxRepetitions];
        TH1D* DifferenceHisto[MaxRepetitions];
        
        bool CompareConsecutive = 1;//compare [cov sample - cov sample-1] or [cov sample - cov max sample]
        bool Maximum = 0;// maximum error or average
        Double_t Error = 0;
        Double_t MaxError;
        const Int_t MaxSamples = Int_t(NSamples/SampleStep);
        TH2D* CovMatrixNSample[MaxSamples][MaxRepetitions];

        for(Int_t repetitions = 1; repetitions<MaxRepetitions+1; repetitions++)
        {
            DifferenceHisto[repetitions] = new TH1D(Form("DifferenceH%d",repetitions),Form("Difference of Covariance Matrices vs Samples%d",repetitions), NSamples,0,NSamples);
           
            CovMatrix500[repetitions] = (TH2D*)File->FindObjectAny(Form("After Covariance Matrix%d %d;%d",0,NSamples+1,repetitions));
            for(Int_t samples = 1; samples<=(NSamples+1); samples+=SampleStep)
            {
                MaxError=0;
                
                CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions] = (TH2D*)File->FindObjectAny(Form("After Covariance Matrix%d %d;%d",0,samples,repetitions));
                
                std::cout << Int_t((samples-1)/SampleStep) << repetitions << Form("After Covariance Matrix%d %d;%d",0,samples,repetitions) << std::endl;
                if(Maximum)
                {
                    if(!CompareConsecutive)
                    {
                        if(samples==1)
                        {
                            continue;
                        }
                        else
                        {
                            for(Int_t i = 0; i<CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetXaxis()->GetNbins(); i++)
                            {
                                for(Int_t j = 0; j<CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetYaxis()->GetNbins(); j++)
                                {
                                    Error = TMath::Abs((CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetBinContent(i+1,j+1) - CovMatrix500[repetitions]->GetBinContent(i+1,j+1))/TMath::Abs(CovMatrix500[repetitions]->GetBinContent(i+1,j+1)));
                                    
                                    if(Error>MaxError)
                                    {
                                        MaxError =  Error;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        if(samples>SampleStep)
                        {
                            for(Int_t i = 0; i<CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetXaxis()->GetNbins(); i++)
                            {
                                for(Int_t j = 0; j<CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetYaxis()->GetNbins(); j++)
                                {
                                    Error = TMath::Abs((CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetBinContent(i+1,j+1) - CovMatrixNSample[Int_t((samples-SampleStep-1))/SampleStep][repetitions]->GetBinContent(i+1,j+1))/TMath::Abs(CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetBinContent(i+1,j+1)));
                                    
                                    if(Error>MaxError)
                                    {
                                        MaxError =  Error;
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    Error = 0;
                    
                    if(!CompareConsecutive)
                    {
                        if(samples==1)
                        {
                            continue;
                        }
                        else
                        {
                            for(Int_t i = 0; i<CovMatrixNSample[(samples-1)/SampleStep][repetitions]->GetXaxis()->GetNbins(); i++)
                            {
                                for(Int_t j = 0; j<CovMatrixNSample[(samples-1)/SampleStep][repetitions]->GetYaxis()->GetNbins(); j++)
                                {
                                    Error = Error + TMath::Abs((CovMatrixNSample[(samples-1)/SampleStep][repetitions]->GetBinContent(i+1,j+1) - CovMatrix500[repetitions]->GetBinContent(i+1,j+1))/TMath::Abs(CovMatrix500[repetitions]->GetBinContent(i+1,j+1)));
                                }
                            }
                        }
                    }
                    else
                    {
                        if(samples>SampleStep)
                        {
                            for(Int_t i = 0; i<CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetXaxis()->GetNbins(); i++)
                            {
                                for(Int_t j = 0; j<CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetYaxis()->GetNbins(); j++)
                                {
                                    Error = Error + TMath::Abs((CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetBinContent(i+1,j+1) - CovMatrixNSample[Int_t((samples-SampleStep-1))/SampleStep][repetitions]->GetBinContent(i+1,j+1))/CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetBinContent(i+1,j+1));
                                    
                                }
                            }
                    
                            Error = Error/(CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetYaxis()->GetNbins()*CovMatrixNSample[Int_t((samples-1)/SampleStep)][repetitions]->GetXaxis()->GetNbins());
                        }
                    }
                }
                
                
                std::cout << " SAMPLE " << samples << " error: " << Error << std::endl;
                
                DifferenceHisto[repetitions]->SetBinContent(samples,Error);
                
            }
            
        }
        
        
        for(Int_t repetitions = 1; repetitions<MaxRepetitions+1; repetitions++)
        {
            TCanvas* DifferenceC = new TCanvas("DifferenceC","CovMatrix difference vs NSamples");

            DifferenceHisto[repetitions]->SetStats(0);
            DifferenceHisto[repetitions]->Draw();
            DifferenceC->Print(Form("./Images/Test/AverageDifferenceCovMatrixVsNSamples%d.png",repetitions));
            delete DifferenceC;
        }
    }
}

void Test::LoganCrossCheck()
{
    
    //Load spectra
    TFile* LoganSpectrumF = new TFile("./CrossChecks/roofile_rawE.root");
    TH1D* LoganPromptH = (TH1D*)LoganSpectrumF->Get("hEp");//prompt spectrum
    TH1D* LoganDelayH = (TH1D*)LoganSpectrumF->Get("hEd");//delay spectrum
    delete LoganSpectrumF;
    
    TFile* MySpectrumF = new TFile("./ResponseMatrices/Hydrogen/NominalResponseMatrix37_39.root");
    
    TH2D* NominalResponseMatrixH = (TH2D*)MySpectrumF->Get("FineEvisEnu1");
    
    delete MySpectrumF;
    TH1D* TrueADH;
    
    TFile* MyTrueADF = TFile::Open("./RootOutputs/Hydrogen/NominalOutputs/Oscillation.root");
    MyTrueADF->cd("Total AD Spectra after oscillation");
    
    TrueADH = (TH1D*)gDirectory->Get("Total spectrum after oscillation at AD1");
    
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
    
    TH1D* ComparisonH = (TH1D*)LoganPromptH->Clone("Comparison");
    
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
    
    ComparisonC->Print("./Images/CrossChecks/LoganPrompt.eps");
    delete ComparisonC;

}

void Near_To_Far_Prediction_VS_Far_Prediction()
{
    //compare pure prediction and relative oscillation model outputs.
}
