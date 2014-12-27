#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2.h"
#include "TMath.h"
#include "NominalData.h"
#include "Prediction.h"
#include <math.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TCanvas.h"
#include "TF1.h"
#include "TGProgressBar.h"
#include "FitterGui.h"

#ifndef ActivateTest
const bool TestSamples = 0;
std::string TestString = "";
//To produce test plots
const bool TESTNEARFARCORR = 0;
const bool TESTFARCORR = 0;
const bool TESTFARANDNEARCORR = 0;
const bool TESTNOCORR = 0;
const bool TESTEH1CORR = 0;
#else
const bool TestSamples = 1;
std::string TestString = "/Test";
//To produce test plots
const bool TESTNEARFARCORR = 0;
const bool TESTFARCORR = 0;
const bool TESTFARANDNEARCORR = 0;
const bool TESTNOCORR = 0;
const bool TESTEH1CORR = 0;
#endif


class CovarianceMatrix3
{
private:
    
    Int_t NSamples;//Number of samples used in the generation of covariance matrices
    
    bool flagX;
    string CovString;
    string TitleString;
    
    NominalData* Nom;
    Prediction* Pred;
    
    Char_t optionN[10];
    std::string AnalysisString;

    Int_t Combine;
    
    Double_t Sin22t13;
    Double_t Dm2_31;
    
    bool Analysis;
    
    bool VaryAccidentalMatrix;
    bool VaryLiHeMatrix;
    bool VaryFastNeutronsMatrix;
    bool VaryAmCMatrix;
    bool DistortLiHeMatrix;
    bool DistortFastNeutronsMatrix;
    bool DistortAmCMatrix;
    
    bool IsotopeMatrix;
    bool ReactorPowerMatrix;
    bool RelativeEnergyOffsetMatrix;
    bool AbsoluteEnergyOffsetMatrix;
    bool AbsoluteEnergyScaleMatrix;
    bool RelativeEnergyScaleMatrix;
    bool IAVMatrixb;
    bool NLMatrix;
    bool ResolutionMatrix;
    bool Sin22t12Matrix;
    bool EfficiencyMatrix;
    
    enum Systematic{IsotopeE,PowerE, RelativeEnergyScaleE, AbsoluteEnergyScaleE, RelativeEnergyOffsetE, AbsoluteEnergyOffsetE,IAVE, NLE, ResolutionE,Sin22t12E,EfficiencyE};
    Systematic SystematicE;
    typedef Systematic SystematicType;
    enum Background{VaryAccidentalE, VaryLiHeE, VaryFastNeutronsE, VaryAmCE, DistortLiHeE, DistortFastNeutronsE, DistortAmCE};
    Background BackgroundE;
    typedef Background BackgroundType;
    
    //AD configuration parameters:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    Int_t Nweeks;
    Int_t NReactorPeriods;
    Int_t MaxNear;
    Int_t MaxFar;
    Int_t MaxBins;
    
    Int_t n_evis_bins;    
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    
    //Histograms
    TH1D* PredictionVisH[MaxNearDetectors][MaxFarDetectors];//  Prediction from near data with no alterations, nominal backgrounds may be added. Vis Binning (after detector response is applied)
    TH1D* AlteredPredictionVisH[MaxNearDetectors][MaxFarDetectors];
    
    TH1D* ReactorPredictionVisH[MaxFarDetectors];//Prediction from reactor with no alterations, nominal backgrounds may be added. Vis Binning (after detector response is applied)
    TH1D* AlteredReactorPredictionVisH[MaxFarDetectors];

    TH2D* Cov2H;
    TH2D* CovMatrix2H;
    TH2D* MatrixBeforeNormalizing;
    TH2D* MatrixAfterNormalizing;
    TH2D* TotalMatrixBeforeNormalizing;
    TH2D* TotalMatrixAfterNormalizing;
    TH2D* SqrtBeforeCovMatrix2H;
    TH2D* SqrtAfterCovMatrix2H;
    TH2D* SqrtCovMatrix2H;
    
    //Loading function
    void LoadNominalPredictions(Int_t);
    void LoadAlteredPredictions(Int_t);
    
    //Functions to generate the respective covariance matrices
    void GenerateCovarianceMatrix(Int_t,Int_t);
    void NormalizeCov(TH2D*);
    void SaveSpectrum(Int_t, Int_t);
    void SaveCovarianceMatrix(Int_t);
    
public:
    CovarianceMatrix3();
    CovarianceMatrix3(NominalData*);
    ~CovarianceMatrix3();
    void CovarianceMatrixMain(FitterGui*);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                         Nominal constructor, not used but coded in case it is needed
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CovarianceMatrix3 :: CovarianceMatrix3()
{
    std::cout << " the covariance matrix default constructor shouldn't be called" << std::endl;
    
    exit(EXIT_FAILURE);
    
    Nom = new NominalData(0,2);
    Nom->SetToyMC(1);//For covariance matrices use Toy MC.
    Pred = new Prediction(Nom);
    Analysis = Nom->GetAnalysis();
    if(Analysis)
    {
        AnalysisString = "Hydrogen";
    }
    else
    {
        AnalysisString = "Gadolinium";
    }
    
    Combine = Nom->GetCombineMode();
    NReactorPeriods=Nom->GetNReactorPeriods();
    
    Sin22t13 = Nom->GetSin22t13();
    Dm2_31 = Nom->GetDm231();
    
    VaryAccidentalMatrix = Nom->GetVaryAccidentalMatrix();
    VaryLiHeMatrix = Nom->GetVaryLiHeMatrix();
    VaryFastNeutronsMatrix = Nom->GetVaryFastNeutronsMatrix();
    VaryAmCMatrix = Nom->GetVaryAmCMatrix();
    DistortLiHeMatrix = Nom->GetDistortLiHeMatrix();
    DistortFastNeutronsMatrix = Nom->GetDistortFastNeutronsMatrix();
    DistortAmCMatrix = Nom->GetDistortAmCMatrix();
    
    IsotopeMatrix = Nom->GetIsotopeMatrix();
    ReactorPowerMatrix = Nom->GetReactorPowerMatrix();
    RelativeEnergyScaleMatrix = Nom->GetRelativeEnergyScaleMatrix();
    RelativeEnergyOffsetMatrix = Nom->GetRelativeEnergyOffsetMatrix();
    AbsoluteEnergyScaleMatrix = Nom->GetAbsoluteEnergyScaleMatrix();
    AbsoluteEnergyOffsetMatrix = Nom->GetAbsoluteEnergyOffsetMatrix();
    IAVMatrixb = Nom->GetIAVMatrix();
    NLMatrix = Nom->GetNLMatrix();
    ResolutionMatrix = Nom->GetResolutionMatrix();
    Sin22t12Matrix = Nom->GetSin22t12Matrix();
    EfficiencyMatrix = Nom->GetEfficiencyMatrix();
    
    if(Analysis&&!(VaryAccidentalMatrix||VaryLiHeMatrix||VaryFastNeutronsMatrix||VaryAmCMatrix||DistortLiHeMatrix||DistortFastNeutronsMatrix||DistortAmCMatrix||IsotopeMatrix||ReactorPowerMatrix||Sin22t12Matrix))
    {
        NSamples = 1;//nH Toy MC is made event by event, with a sufficiently dense number of events the random detector systematics can be applied to the whole spectrum at once. This is not the same for reactor or background systematics though.
    }
    else
    {
        NSamples = Nom->GetNSamples();//LBNL produces ~500 samples per covariance matrix with a different random variation for each sample
    }
    
    Nweeks = Nom->GetWeeks();
    NADs = Nom->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    if(NADs == 8)//    ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    n_evis_bins = Nom->GetVisibleBins();
    
    for (Int_t i = 0; i <= n_evis_bins; i++)
    {
        evis_bins[i] = Nom->GetVisibleBinningArray(i);
    }
    
    MaxBins = 9*n_evis_bins;
    
    BackgroundE = (CovarianceMatrix3::BackgroundType)(-1);//Reset BackgroundE. I put it here because VaryAccidentalMatrix is the first method to be called in the "automatic script"
    SystematicE=(CovarianceMatrix3::SystematicType)(-1);//    A number different to any of the different systematics included in the model
    
    delete Nom;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Constructor carrying all Data information
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CovarianceMatrix3 :: CovarianceMatrix3(NominalData* Data)
{
    Pred = new Prediction(Data);

    Data->SetToyMC(1);//For covariance matrices use Toy MC.

    Analysis = Data->GetAnalysis();
    
    if(Analysis)
    {
        AnalysisString = "Hydrogen";
    }
    else
    {
        AnalysisString = "Gadolinium";
    }
    
    Combine = Data->GetCombineMode();
    NReactorPeriods=Data->GetNReactorPeriods();

    if(Analysis&&!(VaryAccidentalMatrix||VaryLiHeMatrix||VaryFastNeutronsMatrix||VaryAmCMatrix||DistortLiHeMatrix||DistortFastNeutronsMatrix||DistortAmCMatrix||IsotopeMatrix||ReactorPowerMatrix||Sin22t12Matrix))
    {
        NSamples = 1;//nH Toy MC is made event by event, with a sufficiently dense number of events the random detector systematics can be applied to the whole spectrum at once. This is not the same for reactor or background systematics though.
    }
    else
    {
        NSamples = Data->GetNSamples();//LBNL produces ~500 samples per covariance matrix with a different random variation for each sample
    }
    
    Sin22t13 = Data->GetSin22t13();
    Dm2_31 = Data->GetDm231();
    
    VaryAccidentalMatrix = Data->GetVaryAccidentalMatrix();
    VaryLiHeMatrix = Data->GetVaryLiHeMatrix();
    VaryFastNeutronsMatrix = Data->GetVaryFastNeutronsMatrix();
    VaryAmCMatrix = Data->GetVaryAmCMatrix();
    DistortLiHeMatrix = Data->GetDistortLiHeMatrix();
    DistortFastNeutronsMatrix = Data->GetDistortFastNeutronsMatrix();
    DistortAmCMatrix = Data->GetDistortAmCMatrix();
    
    IsotopeMatrix = Data->GetIsotopeMatrix();
    ReactorPowerMatrix = Data->GetReactorPowerMatrix();
    RelativeEnergyScaleMatrix = Data->GetRelativeEnergyScaleMatrix();
    RelativeEnergyOffsetMatrix = Data->GetRelativeEnergyOffsetMatrix();
    AbsoluteEnergyScaleMatrix = Data->GetAbsoluteEnergyScaleMatrix();
    AbsoluteEnergyOffsetMatrix = Data->GetAbsoluteEnergyOffsetMatrix();
    IAVMatrixb = Data->GetIAVMatrix();
    NLMatrix = Data->GetNLMatrix();
    ResolutionMatrix = Data->GetResolutionMatrix();
    Sin22t12Matrix = Data->GetSin22t12Matrix();
    EfficiencyMatrix = Data->GetEfficiencyMatrix();
    
    Nweeks = Data->GetWeeks();
    NADs = Data->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    if(NADs == 8)//    ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    n_evis_bins = Data->GetVisibleBins();
    
    for (Int_t i = 0; i <= n_evis_bins; i++)
    {
        evis_bins[i] = Data->GetVisibleBinningArray(i);
    }
    
    MaxBins = 9*n_evis_bins;

    BackgroundE = (CovarianceMatrix3::BackgroundType)(-1);//Reset BackgroundE. I put it here because VaryAccidentalMatrix is the first method to be called in the "automatic script"
    SystematicE=(CovarianceMatrix3::SystematicType)(-1);//    A number different to any of the different systematics included in the model
}

CovarianceMatrix3 :: ~CovarianceMatrix3()
{
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                  Main function to calculate covariance matrices (Background + Systematics)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CovarianceMatrix3 :: CovarianceMatrixMain(FitterGui* FitterGui)
{
    
    std::cout << " Cov matrix n_evis_bins " << n_evis_bins << std::endl;
    
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
    }
    
    flagX=0;
    
    if(VaryAccidentalMatrix)
    {
        BackgroundE = VaryAccidentalE;
    }
    if(VaryLiHeMatrix)
    {
        BackgroundE = VaryLiHeE;
    }
    if(VaryFastNeutronsMatrix)
    {
        BackgroundE = VaryFastNeutronsE;
    }
    if(VaryAmCMatrix)
    {
        BackgroundE = VaryAmCE;
    }
    if(DistortLiHeMatrix)
    {
        BackgroundE = DistortLiHeE;
    }
    if(DistortFastNeutronsMatrix)
    {
        BackgroundE = DistortFastNeutronsE;
    }
    if(DistortAmCMatrix)
    {
        BackgroundE = DistortAmCE;
    }
    if(IsotopeMatrix)
    {
        SystematicE = IsotopeE;
    }
    if(ReactorPowerMatrix)
    {
        SystematicE = PowerE;
    }
    if(IAVMatrixb)
    {
        SystematicE = IAVE;
    }
    if(NLMatrix)
    {
        SystematicE = NLE;
    }
    if(RelativeEnergyScaleMatrix)
    {
        SystematicE = RelativeEnergyScaleE;
    }
    if(AbsoluteEnergyScaleMatrix)
    {
        SystematicE = AbsoluteEnergyScaleE;
    }
    if(RelativeEnergyOffsetMatrix)
    {
        SystematicE = RelativeEnergyOffsetE;
    }
    if(AbsoluteEnergyOffsetMatrix)
    {
        SystematicE = AbsoluteEnergyOffsetE;
    }
    if(ResolutionMatrix)
    {
        SystematicE = ResolutionE;
    }
    if(Sin22t12Matrix)
    {
        SystematicE = Sin22t12E;
    }
    if(EfficiencyMatrix)
    {
        SystematicE = EfficiencyE;
    }
    
    std::cout << "\t Randomize Background #" << BackgroundE << std::endl;
    std::cout << "\t Randomize Systematic #" << SystematicE << std::endl;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                      Initialize histograms, functions, etc.
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (Int_t week = 0; week<Nweeks; week++)
    {
        Cov2H = new TH2D(Form("Partial Covariance Matrix%d",week),Form("Partial Covariance Matrix%d",week),MaxBins,0,MaxBins,MaxBins,0,MaxBins);
        CovMatrix2H = new TH2D(Form("Covariance Matrix%d",week),Form("Covariance Matrix%d",week),MaxBins,0,MaxBins,MaxBins,0,MaxBins);
        
        if(SystematicE != (CovarianceMatrix3::SystematicType)(-1))
        {
            TotalMatrixBeforeNormalizing = new TH2D(Form("Total Before Covariance Matrix%d",week),Form("Total Before Covariance Matrix%d",week),MaxBins,0,MaxBins,MaxBins,0,MaxBins);
            TotalMatrixAfterNormalizing  = new TH2D(Form("Total After Covariance Matrix%d",week),Form("Total After Covariance Matrix%d",week),MaxBins,0,MaxBins,MaxBins,0,MaxBins);
            std::cout <<  "********************************************************************************************************" << std::endl;
            std::cout << "\t Nominal Prediction" << std::endl;
        }
        
        Pred->MakePrediction(Sin22t13,Dm2_31,0,week,1,1);//  Generate Nominal Prediction
        
        LoadNominalPredictions(week);
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                  Generate toy MC samples
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for (Int_t samples = 0; samples<NSamples; samples++)
        {
            std::cout <<  "********************************************************************************************************" << std::endl;
            std::cout << "\t Random Sample #" << samples << std::endl;
            if (flagX==0)
            {
                sprintf(optionN,"recreate");
                flagX=1;
            }
            else
            {
                sprintf(optionN,"update");
            }
            
            Pred->MakePrediction(Sin22t13,Dm2_31,1,week,1,1);// Prediction with random variations.
            LoadAlteredPredictions(week);
            
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                         Generate Covariance Matrix
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            Cov2H->Reset();
            GenerateCovarianceMatrix(week,samples);
            CovMatrix2H->Add(Cov2H);
            
            SaveSpectrum(samples, week);
  
            if(SystematicE != (CovarianceMatrix3::SystematicType)(-1))
            {
                TotalMatrixBeforeNormalizing->Add(MatrixBeforeNormalizing);
                TotalMatrixAfterNormalizing->Add(MatrixAfterNormalizing);
                delete MatrixBeforeNormalizing;
                delete MatrixAfterNormalizing;
            }
            
            for (Int_t far =0; far<MaxFar; far++)
            {
                delete AlteredReactorPredictionVisH[far];
                
                for (Int_t near = 0; near<MaxNear; near++)
                {
                    delete AlteredPredictionVisH[near][far];
                }
            }
            if(!TestSamples)
            {
               FitterGui->Update(samples);
            }
        }
        if(SystematicE != (CovarianceMatrix3::SystematicType)(-1))
        {
            TotalMatrixBeforeNormalizing->Scale(1./(NSamples));
            TotalMatrixAfterNormalizing->Scale(1./(NSamples));
        }
        CovMatrix2H->Scale(1./(NSamples));
        
        SaveCovarianceMatrix(week);
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                            Clean up the dust
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        std::cout << "\t Cleaning up " << std::endl;
        
        delete CovMatrix2H;
        delete Cov2H;
        
        if(SystematicE != (CovarianceMatrix3::SystematicType)(-1))
        {
            delete TotalMatrixAfterNormalizing;
            delete TotalMatrixBeforeNormalizing;
        }

        for (Int_t far =0; far<MaxFar; far++)
        {
            delete ReactorPredictionVisH[far];
            
            for (Int_t near = 0; near<MaxNear; near++)
            {
                delete PredictionVisH[near][far];
            }
        }
        
    }
    
    delete Pred;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Produces the Covariance Matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CovarianceMatrix3 :: GenerateCovarianceMatrix(Int_t week,Int_t samples)
{
    Double_t Cov[MaxBins][MaxBins];
    Double_t CopyCov[MaxBins][MaxBins];
    
    TH2D* CopyCov2H;
    if(SystematicE != (CovarianceMatrix3::SystematicType)(-1))
    {
        CopyCov2H = new TH2D(Form("Partial Before Covariance Matrix%d",week),Form("Partial Before Covariance Matrix%d",week),MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    }

    std::cout << "\t Generating Covariance Matrix" << std::endl;
    
    double factorFar,factorNear = 0;

        factorFar = 1;

        if(TESTNEARFARCORR)
        {
            for (Int_t fari=0; fari<MaxFar; fari++)
            {
                for (Int_t i = 0; i<n_evis_bins; i++)
                {
                    AlteredReactorPredictionVisH[fari]->SetBinContent(i+1,factorFar);
                    ReactorPredictionVisH[fari]->SetBinContent(i+1,1);
                }
                
                for (Int_t neari=0; neari<MaxNear; neari++)
                {
                    factorNear = (neari+2);

                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {
                        AlteredPredictionVisH[neari][fari]->SetBinContent(i+1,factorNear);
                        PredictionVisH[neari][fari]->SetBinContent(i+1,1);
                    }
                }
            }
        }
        if(TESTFARCORR)
        {
            factorNear = 1;

            for (Int_t fari=0; fari<MaxFar; fari++)
            {
                factorFar = fari+2;
                
                for (Int_t i = 0; i<n_evis_bins; i++)
                {
                    AlteredReactorPredictionVisH[fari]->SetBinContent(i+1,factorFar);
                    ReactorPredictionVisH[fari]->SetBinContent(i+1,1);
                }
                
                for (Int_t neari=0; neari<MaxNear; neari++)
                {
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {
                        AlteredPredictionVisH[neari][fari]->SetBinContent(i+1,factorNear);
                        PredictionVisH[neari][fari]->SetBinContent(i+1,1);
                    }
                }
            }
        }
        if(TESTFARANDNEARCORR)
        {
            for (Int_t fari=0; fari<MaxFar; fari++)
            {
                factorFar = (2);
                
                for (Int_t i = 0; i<n_evis_bins; i++)
                {
                    AlteredReactorPredictionVisH[fari]->SetBinContent(i+1,factorFar);
                    ReactorPredictionVisH[fari]->SetBinContent(i+1,1);
                }
                
                for (Int_t neari=0; neari<MaxNear; neari++)
                {
                    factorNear = (2);
                    
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {
                        AlteredPredictionVisH[neari][fari]->SetBinContent(i+1,factorNear);
                        PredictionVisH[neari][fari]->SetBinContent(i+1,1);
                    }
                }
            }
        }
        if(TESTNOCORR)
        {
            factorFar = 1;
            factorNear = 1;

            for (Int_t fari=0; fari<MaxFar; fari++)
            {
                for (Int_t i = 0; i<n_evis_bins; i++)
                {
                    AlteredReactorPredictionVisH[fari]->SetBinContent(i+1,factorFar);
                    ReactorPredictionVisH[fari]->SetBinContent(i+1,1);
                }
                
                for (Int_t neari=0; neari<MaxNear; neari++)
                {
                    
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {
                        AlteredPredictionVisH[neari][fari]->SetBinContent(i+1,factorNear);
                        PredictionVisH[neari][fari]->SetBinContent(i+1,1);
                    }
                }
            }
        }
        if(TESTEH1CORR)
        {
            factorFar = 1;

            for (Int_t fari=0; fari<MaxFar; fari++)
            {
                for (Int_t i = 0; i<n_evis_bins; i++)
                {
                    AlteredReactorPredictionVisH[fari]->SetBinContent(i+1,factorFar);
                    ReactorPredictionVisH[fari]->SetBinContent(i+1,1);
                }
                
                for (Int_t neari=0; neari<MaxNear; neari++)
                {
                    if(neari<ADsEH1)
                    {
                        factorNear = (2);
                    }
                    else
                    {
                        factorNear = 3;
                    }
                    
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {
                        AlteredPredictionVisH[neari][fari]->SetBinContent(i+1,factorNear);
                        PredictionVisH[neari][fari]->SetBinContent(i+1,1);
                    }
                }
            }
        }
    
    Int_t x =0;
    Int_t y =0;
    
    Int_t Ni1=0,Ni2=0,Ni3=0,Ni4 =0;
    Int_t Nj1=0,Nj2=0,Nj3=0,Nj4 =0;
    Int_t Fi1=0,Fi2=0,Fi3=0,Fi4 =0;
    Int_t Fj1=0,Fj2=0,Fj3=0,Fj4 =0;
    
    for (Int_t neari=0; neari<MaxNear; neari++)
    {
        //Logic for the 2D matrix index done up to 8 ADs
        if(neari==0){Ni1=1;Ni2=0;Ni3=0;Ni4 =0;}
        if(neari==1){Ni2++;}
        if(neari==2){Ni3++;}
        if(neari==3){Ni4++;}
        
        for (Int_t fari=0; fari<MaxFar; fari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(Ni1!=Ni2){Fi1=fari+1;Fi2=0;Fi3=0;Fi4 =0;}
            if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFar; Fi2=fari+1;}
            if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFar; Fi3=fari+1;}
            if(Ni3==Ni4&&Ni4==1){Fi3=MaxFar; Fi4=fari+1;}
            //                std::cout <<"\n"<< "Ni" << Ni1 << Ni2 << Ni3 << Ni4 << "\n";
            //                std::cout << "Fi"<< Fi1 << Fi2 << Fi3 << Fi4 << "\n";
            //                std::cout << "====================================" << icounter++ <<"\n";
            //                std::cout << Ni1*Fi1 << Ni2*Fi2 << Ni3*Fi3 << Ni4*Fi4;
            
            
            for (Int_t nearj=0; nearj<MaxNear; nearj++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4 =0;}
                if(nearj==1){Nj2++;}
                if(nearj==2){Nj3++;}
                if(nearj==3){Nj4++;}
                
                for (Int_t farj=0; farj<MaxFar; farj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(Nj1!=Nj2){Fj1=farj+1;Fj2=0;Fj3=0;Fj4 =0;}
                    if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFar; Fj2=farj+1;}
                    if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFar; Fj3=farj+1;}
                    if(Nj3==Nj4&&Nj4==1){Fj3=MaxFar; Fj4=farj+1;}
                    //                        std::cout <<"\n"<< "Nj"<< Nj1 << Nj2 << Nj3 << Nj4 << "\n";
                    //                        std::cout << "Fj"<< Fj1 << Fj2 << Fj3 << Fj4 << "\n";
                    //                        std::cout << "===================================="<< jcounter++ <<"\n";
                    //                        std::cout << Nj1*Fj1 << Nj2*Fj2 << Nj3*Fj3 << Nj4*Fj4 <<"\n";
                    //                        std::cout <<"\n";
                    
//                    #ifdef Print
//                        if(BackgroundE == (CovarianceMatrix3::BackgroundType)(-1))
//                        {
//                            TCanvas* PredC = new TCanvas("","");
//                            PredC->Divide(2,4);
//                            PredC->cd(1);
//                            PredictionVisH[neari][fari]->Draw();
//                            PredC->cd(2);
//                            PredictionVisH[nearj][farj]->Draw();
//                            PredC->cd(3);
//                            AlteredPredictionVisH[neari][fari]->Draw();
//                            PredC->cd(4);
//                            AlteredPredictionVisH[nearj][farj]->Draw();
//                            
//                            PredC->cd(5);
//                            ReactorPredictionVisH[fari]->Draw();
//                            PredC->cd(6);
//                            ReactorPredictionVisH[farj]->Draw();
//                            PredC->cd(7);
//                            AlteredReactorPredictionVisH[fari]->Draw();
//                            PredC->cd(8);
//                            AlteredReactorPredictionVisH[farj]->Draw();
//                            
//                            PredC->Print(Form("./Images/CovarianceMatrixPredictionsNeari%d_Fari%d_Nearj%d_Farj%d.eps",neari,fari,nearj,farj));
//                        }
//                    #endif
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {//columns
                        x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                        
                        for (Int_t j = 0; j<n_evis_bins; j++)
                        {//rows
                            y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                            
                            //  To calculate covariance matrix without background/systematic optimization:
                            
                            //                            Cov[x][y]=((FarRandomBackgroundSpectrumH[fari]->GetBinContent(i+1)+AlteredPredictionVisH[neari][fari]->GetBinContent(i+1))-(PredictionVisH[neari][fari]->GetBinContent(i+1)+FarBackgroundSpectrumH[fari]->GetBinContent(i+1)))*((FarRandomBackgroundSpectrumH[farj]->GetBinContent(j+1)+AlteredPredictionVisH[nearj][farj]->GetBinContent(j+1))-(PredictionVisH[nearj][farj]->GetBinContent(j+1)+FarBackgroundSpectrumH[farj]->GetBinContent(j+1)));
                            
//                            if(BackgroundE != (CovarianceMatrix3::BackgroundType)(-1))// Background has been varied
//                            {
//                                Cov[x][y]=(FarRandomBackgroundSpectrumH[fari]->GetBinContent(i+1)-FarBackgroundSpectrumH[fari]->GetBinContent(i+1))*(FarRandomBackgroundSpectrumH[farj]->GetBinContent(j+1)-FarBackgroundSpectrumH[farj]->GetBinContent(j+1));
//                            }
//                            else // Systematic has been varied, or no variations at all, in the latter case this should be 0. (Check)
//                            {
                                Cov[x][y]=(AlteredReactorPredictionVisH[fari]->GetBinContent(i+1) - ReactorPredictionVisH[fari]->GetBinContent(i+1) + AlteredPredictionVisH[neari][fari]->GetBinContent(i+1)-PredictionVisH[neari][fari]->GetBinContent(i+1))*(AlteredReactorPredictionVisH[farj]->GetBinContent(j+1) - ReactorPredictionVisH[farj]->GetBinContent(j+1)+ AlteredPredictionVisH[nearj][farj]->GetBinContent(j+1)-PredictionVisH[nearj][farj]->GetBinContent(j+1));
//                            }
                            
                            //                                std::cout << "COV VALUE" << Cov[x][y] << std::endl;
                            Cov2H->SetBinContent(x+1,y+1,Cov[x][y]);
                            
                            ////////////////////////////////Normilize systematic and then unnormalize when fitting for different sin2t13. Check also for NaNs.
                            if(SystematicE != (CovarianceMatrix3::SystematicType)(-1))
                            {
                                CopyCov[x][y]=Cov[x][y];
                                
                                Cov[x][y] = CopyCov[x][y]/(PredictionVisH[neari][fari]->GetBinContent(i+1)*PredictionVisH[nearj][farj]->GetBinContent(j+1));
                                
                                if((PredictionVisH[neari][fari]->GetBinContent(i+1)*PredictionVisH[nearj][farj]->GetBinContent(j+1))==0)
                                {
                                    
                                Cov[x][y] = 0;//Empty bins will produce NaNs.
                                    
                                    std::cout << "\t \t \t WARNING PREDi*PREDj IS 0 SO YOUR SYSTEMATIC MATRIX WILL BE NAN!!!!! FIX " << std::endl;
                                  //  exit(EXIT_FAILURE);
                                }
                            }
                            
                            Cov2H->SetBinContent(x+1,y+1,Cov[x][y]);
                            
                            if(SystematicE != (CovarianceMatrix3::SystematicType)(-1))
                            {
                                CopyCov2H->SetBinContent(x+1,y+1,CopyCov[x][y]);
                            }
                        }
                    }
                }
            }
        }
    }
    
    if(SystematicE != (CovarianceMatrix3::SystematicType)(-1))
    {
        MatrixBeforeNormalizing=(TH2D*)CopyCov2H->Clone(Form("cov matrix before normalizing, sys%d and bkgd%d", SystematicE,BackgroundE));
        MatrixAfterNormalizing=(TH2D*)Cov2H->Clone(Form("cov matrix after normalizing, sys%d and bkgd%d", SystematicE,BackgroundE));
        
        delete CopyCov2H;
    }

    //    NormalizeCov(Cov2H);
    //    Cov2H->Draw("colz");
}

void CovarianceMatrix3 :: SaveCovarianceMatrix(Int_t week)
{
    std::cout << "\t Saving Covariance Matrix" << std::endl;
    
    Char_t filenameCov[100];
    Char_t txtfilenameCov[100];
    
    //  In case there are no variations:

    sprintf(filenameCov,("."+TestString+"/CovarianceMatrices/"+ AnalysisString+ "/Combine%d/CovarianceMatricesRoot/NominalCovarianceMatrix.root").c_str(),Combine);
    sprintf(txtfilenameCov,("."+TestString+"/CovarianceMatrices/"+ AnalysisString+ "/Combine%d/CovarianceMatricesTxT/NominalCovarianceMatrix.txt").c_str(),Combine);
    
    //  Save Cov Matrix (Either background or systematics)
    switch (BackgroundE)
    {
        case 0://   Vary Accidentals
            CovString = "VaryAccidental";
            TitleString="Vary Accidental Covariance Matrix";
            break;
        case 1://   Vary LiHe
            CovString = "VaryLiHe";
            TitleString="Vary LiHe Covariance Matrix";
            break;
        case 2://   Vary Fast Neutrons
            CovString = "VaryFastNeutrons";
            TitleString="Vary FN Covariance Matrix";
            break;
        case 3://   Vary AmC
            CovString = "VaryAmC";
            TitleString="Vary AmC Covariance Matrix";
            break;
        case 4://   Distort LiHe
            CovString = "DistortLiHe";
            TitleString="Distort LiHe Covariance Matrix";
            break;
        case 5://   Distort Fast Neutrons
            CovString = "DistortFastNeutrons";
            TitleString="Distort FN Covariance Matrix";
            break;
        case 6://   Distort AmC
            CovString = "DistortAmC";
            TitleString="Distort AmC Covariance Matrix";
            break;
    }
    switch (SystematicE)
    {
        case 0://   Vary Reactor
            CovString = "Isotope";
            TitleString="Isotope Covariance Matrix";
            break;
        case 1://   Vary Reactor
            CovString = "ReactorPower";
            TitleString="Reactor Power Covariance Matrix";
            break;
        case 2://   Vary Energy Scale
            CovString = "RelativeEnergyScale";
            TitleString="Relative Energy Scale Covariance Matrix";
            break;
        case 3://   Vary Absolute Energy Scale
            CovString = "AbsoluteEnergyScale";
            TitleString="Absolute Energy Scale Covariance Matrix";
            break;
        case 4://   Vary Energy Offset
            CovString = "RelativeEnergyOffset";
            TitleString="Relative Energy Offset Covariance Matrix";
            break;
        case 5://   Vary Absolute Energy Offset
            CovString = "AbsoluteEnergyOffset";
            TitleString="Absolute Energy Offset Covariance Matrix";
            break;
        case 6://   Vary IAV
            CovString = "IAV";
            TitleString="IAV Covariance Matrix";
            break;
        case 7://   Vary NL
            CovString = "NL";
            TitleString="NL Covariance Matrix";
            break;
        case 8://   Vary Resolution
            CovString = "Resolution";
            TitleString="Resolution Covariance Matrix";
            break;
        case 9://   Vary Sin22t12
            CovString = "Sin22t12";
            TitleString = "Sin22t12 Covariance Matrix";
            break;
        case 10:
            CovString = "Efficiency";
            TitleString = "Efficiency Covariance Matrix";
            break;
    }
    
    sprintf(filenameCov,("."+TestString+"/CovarianceMatrices/"+AnalysisString+ "/Combine%d/CovarianceMatricesRoot/"+CovString+".root").c_str(),Combine);
    
    sprintf(txtfilenameCov,("."+TestString+"/CovarianceMatrices/"+AnalysisString+ "/Combine%d/CovarianceMatricesTxT/"+CovString+"Period%d.txt").c_str(),Combine,week);
    
    CovMatrix2H->SetTitle((TitleString).c_str());
    
    Char_t CharUpdate[20];
    
    if(TestSamples)
    {
      sprintf(CharUpdate,"update");
    }
    else
    {
        sprintf(CharUpdate,"recreate");
    }
    
    TFile* SaveCovarianceMatrixF = new TFile(filenameCov,CharUpdate);

    if(SystematicE != (CovarianceMatrix3::SystematicType)(-1))
    {
        if(TestSamples)
        {
            TotalMatrixBeforeNormalizing->Write(Form("Before Covariance Matrix%d %d",week,NSamples));
            TotalMatrixAfterNormalizing->Write(Form("After Covariance Matrix%d %d",week,NSamples));
        }
        else
        {
            TotalMatrixBeforeNormalizing->Write(Form("Before Covariance Matrix%d",week));
            TotalMatrixAfterNormalizing->Write(Form("After Covariance Matrix%d",week));
            CovMatrix2H->Write();
        }
    }
    else
    {
        CovMatrix2H->Write();
    }
    
    delete SaveCovarianceMatrixF;
    
    #ifdef PrintEps
    
    TH1D* DiagonalCov = new TH1D(("sigma "+CovString).c_str(),("sigma "+CovString).c_str(),Cov2H->GetXaxis()->GetNbins(),0,CovMatrix2H->GetXaxis()->GetNbins());
    
    for(Int_t i = 0; i < Cov2H->GetXaxis()->GetNbins(); i++)
    {
        DiagonalCov->SetBinContent(i+1,CovMatrix2H->GetBinContent(i+1,i+1));//Diagonal of the covariance matrix = sigma
    }
    
    TCanvas* sigmaC = new TCanvas("sigmaC","sigmaC");
    
    DiagonalCov->Draw();
    
    sigmaC->Print((Form("./Images/Sigma_Combine%d_",Combine)+CovString+".eps").c_str(),".eps");
    
    delete sigmaC;
    
    delete DiagonalCov;

    TCanvas* TestCovariance = new TCanvas("","");
    
    CovMatrix2H->Draw("colz");
    
    if(TESTFARCORR)
    {
        TestCovariance->Print("./Images/CovarianceMatrices/CovarianceMatrixTests/TestFar.eps");
    }
    if(TESTNEARFARCORR)
    {
        TestCovariance->Print("./Images/CovarianceMatrices/CovarianceMatrixTests/TestNeartoFar.eps");
    }
    if(TESTFARANDNEARCORR)
    {
        TestCovariance->Print("./Images/CovarianceMatrices/CovarianceMatrixTests/TestFarAndNeartoFar.eps");
    }
    if(TESTNOCORR)
    {
        TestCovariance->Print("./Images/CovarianceMatrices/CovarianceMatrixTests/TestNoCorr.eps");
    }
    if(TESTEH1CORR)
    {
        TestCovariance->Print("./Images/CovarianceMatrices/CovarianceMatrixTests/TestEH1Only.eps");
    }
    delete TestCovariance;
    #endif
    
    //Save in a txt file
    if(WriteOutput)//I will have to do the same here, select a covariance matrix name depending on what has been calculated
    {
        ofstream covf(txtfilenameCov);
        
        Int_t x =0;
        Int_t y =0;
        Int_t Ni1=0,Ni2=0,Ni3=0,Ni4 =0;
        Int_t Nj1=0,Nj2=0,Nj3=0,Nj4 =0;
        Int_t Fi1=0,Fi2=0,Fi3=0,Fi4 =0;
        Int_t Fj1=0,Fj2=0,Fj3=0,Fj4 =0;
        
        for (Int_t neari=0; neari<MaxNear; neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Ni1=1;Ni2=0;Ni3=0;Ni4 =0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<MaxFar; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Fi1=fari+1;Fi2=0;Fi3=0;Fi4 =0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFar; Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFar; Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=MaxFar; Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<MaxNear; nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<MaxFar; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Fj1=farj+1;Fj2=0;Fj3=0;Fj4 =0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFar; Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFar; Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=MaxFar; Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                            
                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                
                                covf << CovMatrix2H->GetBinContent(x+1,y+1)<< " ";
                            }
                        }
                    }
                }
            }
        }
        covf << std::endl;
        covf.close();
    }
}

void CovarianceMatrix3 :: SaveSpectrum(Int_t sample, Int_t week)
{

    std::cout << "\t Saving Spectrum" << std::endl;

    Char_t filenameSpec[200];

    switch (BackgroundE)
    {
        case 0://Vary Accidentals
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/VaryAccidentals.root").c_str(),Combine);
            break;
        case 1://Vary LiHe
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/VaryLiHe.root").c_str(),Combine);
            break;
        case 2://Vary Fast Neutrons
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/VaryFastNeutrons.root").c_str(),Combine);
            break;
        case 3://Vary AmC
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/VaryAmC.root").c_str(),Combine);
            break;
        case 4://Distort LiHe
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/DistortLiHe.root").c_str(),Combine);
            break;
        case 5://Distort Fast Neutrons
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/DistortFastNeutrons.root").c_str(),Combine);
            break;
        case 6://Distort AmC
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/DistortAmC.root").c_str(),Combine);
            break;
        default:
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/Nominal.root").c_str(),Combine);
            break;
    }
    switch (SystematicE)
    {
        case 0://Vary Reactor
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/Isotope.root").c_str(),Combine);
            break;
        case 1://Vary Reactor
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/Power.root").c_str(),Combine);
            break;
        case 2://Vary Energy Scale
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/RelativeEnergyScale.root").c_str(),Combine);
            break;
        case 3://Vary Energy Offset
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/RelativeEnergyOffset.root").c_str(),Combine);
            break;
        case 4://Vary Absolute Scale
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/AbsoluteEnergyScale.root").c_str(),Combine);
            break;
        case 5://Vary Absolute Offset
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/AbsoluteEnergyOffset.root").c_str(),Combine);
            break;
        case 6://Vary IAV
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/IAV.root").c_str(),Combine);
            break;
        case 7://Vary NL
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/NL.root").c_str(),Combine);
            break;
        case 8://Vary Resolution
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/Resolution.root").c_str(),Combine);
            break;
        case 9://Vary Sin22t12
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/Sin22t12.root").c_str(),Combine);
            break;
        case 10:
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/Efficiency.root").c_str(),Combine);
            break;
        default:
            sprintf(filenameSpec,("./CovarianceMatrices/"+AnalysisString+ "/Combine%d/Spectrum/Nominal.root").c_str(),Combine);
            break;
    }

    TFile* SpectrumAndPredictionsF = TFile::Open(filenameSpec,optionN);

        for (Int_t far =0; far<MaxFar; far++)
        {
            for (Int_t near = 0; near < MaxNear; near++)
            {
                TH1D* CopyPredictionVisH;

                PredictionVisH[near][far]->Write(Form("Nominal Far AD%d prediction from Near AD%d VisH Sample%i Period%d", far, near,sample, week));
                AlteredPredictionVisH[near][far]->Write(Form("Varied Far AD%d prediction from Near AD%d VisH Sample%i Period%d",far, near,sample, week));
                CopyPredictionVisH=(TH1D*)AlteredPredictionVisH[near][far]->Clone();
                CopyPredictionVisH->Add(PredictionVisH[near][far],-1);
                CopyPredictionVisH->Divide(PredictionVisH[near][far]);
                CopyPredictionVisH->Write(Form("Relative error Far AD%d from Near AD%d Spectra Sample%i Period%d",far,near,sample,week));
                
                delete CopyPredictionVisH;
            }
        }
    
    SpectrumAndPredictionsF->Close();
}

void CovarianceMatrix3 :: LoadNominalPredictions(Int_t week)
{
    
    //Far predictions from near prediction:
    TFile* PredictionF = new TFile(("./RootOutputs/"+AnalysisString+Form("/Spectra/Combine%d/NominalPredictedSpectrum.root",Combine)).c_str());

    for (Int_t near = 0; near<MaxNear; near++)
    {
        for (Int_t far =0; far<MaxFar; far++)
        {
            PredictionVisH[near][far]=(TH1D*)gDirectory->Get(Form("Combined Prediction AD%d from AD%d",far+1,near+1))->Clone(Form("Nominal Prediction AD%d from AD%d period%d",far+1,near+1,week));
        }
    }
    
    //Far predictions:
    for (Int_t far =0; far<MaxFar; far++)
    {
        ReactorPredictionVisH[far]=(TH1D*)gDirectory->Get(Form("Combined Reactor Prediction AD%d",far+1));
    }
    
    delete PredictionF;
}

void CovarianceMatrix3 :: LoadAlteredPredictions(Int_t week)//  This could be handled with only 1 load method and passing as parameter the TH1D that we want to use and return
{
    //Far predictions from near prediction:
    TFile* PredictionF1 = new TFile(("./RootOutputs/"+AnalysisString+ Form("/Spectra/Combine%d/PredictedSpectrum.root",Combine)).c_str());
    
    for (Int_t near = 0; near<MaxNear; near++)
    {
        for (Int_t far =0; far<MaxFar; far++)
        {
            AlteredPredictionVisH[near][far]=(TH1D*)gDirectory->Get(Form("Combined Prediction AD%d from AD%d",far+1,near+1));
            AlteredPredictionVisH[near][far]->SetName(Form("Altered Prediction AD%d from AD%d period%d",far+1,near+1,week));
        }
    }
    
    //Far predictions:
    for (Int_t far =0; far<MaxFar; far++)
    {
        AlteredReactorPredictionVisH[far]=(TH1D*)gDirectory->Get(Form("Combined Reactor Prediction AD%d",far+1));
    }
    
    delete PredictionF1;
    
}
