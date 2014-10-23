#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2.h"
#include "TMath.h"
#include "NominalData.h"
#include "Prediction.h"
#include <math.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TCanvas.h"
#include "TF1.h"

const bool WriteOutput=0;//To save the covariance matrices in a .txt file.

class CovarianceMatrix3
{
private:
    
    Int_t NSamples;//Number of samples used in the generation of covariance matrices

    bool flagX;
    TCanvas* c;
    NominalData* Nom;
    Prediction* Pred;
    TRandom3* rand;

    Int_t Combine;
    Double_t Sin22t13;
    
    bool VaryRate;
    bool Distort;
    
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
    
    enum Systematic{IsotopeE,PowerE, RelativeEnergyE, AbsoluteEnergyE, RelativeEnergyOffsetE, AbsoluteEnergyOffsetE, IAVE, NLE, ResolutionE,Sin22t12E};
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
    Int_t hall;
    
    bool LinearBinning;
    Int_t n_evis_bins;
    Int_t n_etrue_bins;
    Double_t InitialEnergy;
    Double_t FinalVisibleEnergy;
    Double_t BinWidth;

    Int_t MaxNear;
    Int_t MaxFar;
    Int_t MaxBins;

    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.
    
    //  Background rates:
    Double_t ScaleAcc[MaxDetectors][MaxPeriods];
    Double_t ScaleLiHe[MaxDetectors][MaxPeriods];
    Double_t ScaleFN[MaxDetectors][MaxPeriods];
    Double_t ScaleAmC[MaxDetectors][MaxPeriods];
    
    //Background errors:
    Double_t AccidentalError[MaxDetectors][MaxPeriods];
    Double_t LiHeError[MaxDetectors][MaxPeriods];
    Double_t FastNeutronsError[MaxDetectors][MaxPeriods];
    Double_t AmCError[MaxDetectors][MaxPeriods];
    
    //Distortion functions:
    Double_t DistortLiHe;
    Double_t DistortFN;
    Double_t DistortAmC;
    
    Double_t ScaleFactorAccidental[MaxDetectors][MaxPeriods];
    Double_t ScaleFactorLiHe[MaxDetectors][MaxPeriods];
    Double_t ScaleFactorFastNeutrons[MaxDetectors][MaxPeriods];
    Double_t ScaleFactorAmC[MaxDetectors][MaxPeriods];
    
    //Histograms
    TH1F* OriginalPredictionH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    TH1F* PredictionTrueH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];// Prediction with no alterations, nominal backgrounds may be added. True Binning.
    TH1F* PredictionVisH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];//  Prediction with no alterations, nominal backgrounds may be added. Vis Binning (after detector response is applied)
    TH1F* AlteredPredictionVisH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    TH1F* CopyPredictionVisH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    
    TH1F* OriginalADSpectrumH[MaxDetectors][MaxPeriods];//  Spectrum that may include reactor alterations if reactor covariance matrices are produced.
    TH1F* ADSpectrumTrueH[MaxDetectors][MaxPeriods];
    TH1F* ADSpectrumVisH[MaxDetectors][MaxPeriods];
    TH1F* AlteredADSpectrumVisH[MaxDetectors][MaxPeriods];
    
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
    TH1F* RandomAccidentalsH[MaxDetectors][MaxPeriods];
    TH1F* RandomLiHeH[MaxDetectors][MaxPeriods];
    TH1F* RandomFastNeutronsH[MaxDetectors][MaxPeriods];
    TH1F* RandomAmCH[MaxDetectors][MaxPeriods];
    
    TH2F* CovMatrix2H[MaxPeriods];
    TH2F* Cov2H[MaxPeriods];
    
    Double_t Cov[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t NormCov[9*MaxNbins][9*MaxNbins][MaxPeriods];
    
    //Loading function
    void LoadBackgrounds();
    void LoadNominalPredictions(Int_t);
    void LoadAlteredPredictions(Int_t);
    
    //Functions to vary background shapes. So far the same ones than LBNL.
    TF1* GetDistortionFunction(Double_t);
    TF1* GetFastNeutronsDistortionFunction(Double_t);
    void FluctuateBackgrounds(Int_t);

    //Functions to generate the respective covariance matrices
    void GenerateCovarianceMatrix(Int_t);
    void NormalizeCov(TH2F*);
    void SaveSpectrum(Int_t, Int_t);
    void SaveCovarianceMatrix(Int_t);
    
public:
    CovarianceMatrix3();
    CovarianceMatrix3(NominalData*);
    
    void CovarianceMatrixMain();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                             Nominal constructor, not used but coded in case it is needed
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CovarianceMatrix3 :: CovarianceMatrix3()
{
    Nom = new NominalData();
    c = new TCanvas("canvas");
    rand = new TRandom3();
    Pred = new Prediction(Nom);
    
    Combine = Nom->GetCombineMode();

    NSamples = Nom->GetNSamples();

    Sin22t13 = Nom->GetSin22t13();
    
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
    IAVMatrixb = Nom->GetIAVMatrix();
    NLMatrix = Nom->GetNLMatrix();
    ResolutionMatrix = Nom->GetResolutionMatrix();
    
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
    
    LinearBinning = Nom->GetBinning();
    
    //  Linear binning
    if(LinearBinning)
    {
        n_evis_bins = Nom->GetNbins();
        n_etrue_bins = Nom->GetNbins();
        
        for (Int_t i = 0; i <= n_evis_bins; i++)
        {
            evis_bins[i] = 0.2 * i + 0.7;
            enu_bins[i] = 0.2 * i + InitialEnergy;
        }
    }
    //  Non-linear binning
    else
    {
        n_evis_bins=37;
        n_etrue_bins=39;
        
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
    }
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            AccidentalError[AD][week]=Nom->GetAccidentalError(AD,week);
            LiHeError[AD][week]=Nom->GetLiHeError(AD,week);
            FastNeutronsError[AD][week]=Nom->GetFNError(AD,week);
            AmCError[AD][week]=Nom->GetAmCError(AD,week);
            
            if(AD<ADsEH1)
            {
                hall=1;
            }
            else if(AD >= ADsEH1 && AD < (ADsEH1+ADsEH2))
            {
                hall=2;
            }
            else
            {
                hall=3;
            }
            
            ScaleAcc[AD][week]=Nom->GetAccidentalEvents(AD,week);
            ScaleLiHe[AD][week]=Nom->GetLiHeEvents(hall,week);
            ScaleFN[AD][week]=Nom->GetFNEvents(hall,week);
            ScaleAmC[AD][week]=Nom->GetAmCEvents(hall,week);
        }
    }
    BackgroundE = (CovarianceMatrix3::BackgroundType)(-1);//Reset BackgroundE. I put it here because VaryAccidentalMatrix is the first method to be called in the "automatic script"
    SystematicE=(CovarianceMatrix3::SystematicType)(-1);//    A number different to any of the different systematics included in the model
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Constructor carrying all Data information
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CovarianceMatrix3 :: CovarianceMatrix3(NominalData* Data)
{
    Pred = new Prediction(Data);
    c = new TCanvas("canvas");
    rand = new TRandom3();

    Combine = Data->GetCombineMode();

    NSamples = Data->GetNSamples();

    Sin22t13 = Data->GetSin22t13();

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
    IAVMatrixb = Data->GetIAVMatrix();
    NLMatrix = Data->GetNLMatrix();
    ResolutionMatrix = Data->GetResolutionMatrix();

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
    
    LinearBinning = Data->GetBinning();
    
    //  Linear binning
    if(LinearBinning)
    {
        n_evis_bins = Data->GetNbins();
        n_etrue_bins = Data->GetNbins();
        
        for (Int_t i = 0; i <= n_evis_bins; i++)
        {
            evis_bins[i] = 0.2 * i + 0.7;
            enu_bins[i] = 0.2 * i + InitialEnergy;
        }
    }
    //  Non-linear binning
    else
    {
        n_evis_bins=37;
        n_etrue_bins=39;
        
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
    }
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            AccidentalError[AD][week]=Data->GetAccidentalError(AD,week);
            LiHeError[AD][week]=Data->GetLiHeError(AD,week);
            FastNeutronsError[AD][week]=Data->GetFNError(AD,week);
            AmCError[AD][week]=Data->GetAmCError(AD,week);
            
            if(AD<ADsEH1)
            {
                hall=1;
            }
            else if(AD >= ADsEH1 && AD < (ADsEH1+ADsEH2))
            {
                hall=2;
            }
            else
            {
                hall=3;
            }
            
            ScaleAcc[AD][week]=Data->GetAccidentalEvents(AD,week);
            ScaleLiHe[AD][week]=Data->GetLiHeEvents(hall,week);
            ScaleFN[AD][week]=Data->GetFNEvents(hall,week);
            ScaleAmC[AD][week]=Data->GetAmCEvents(hall,week);
        }
    }
    BackgroundE = (CovarianceMatrix3::BackgroundType)(-1);//Reset BackgroundE. I put it here because VaryAccidentalMatrix is the first method to be called in the "automatic script"
    SystematicE=(CovarianceMatrix3::SystematicType)(-1);//    A number different to any of the different systematics included in the model
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Main function to calculate covariance matrices (Background + Systematics)  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CovarianceMatrix3 :: CovarianceMatrixMain()
{
    if(IsotopeMatrix)
    {
        SystematicE = IsotopeE;
    }
    else if(ReactorPowerMatrix)
    {
        SystematicE = PowerE;
    }
    else if(IAVMatrixb)
    {
        SystematicE = IAVE;
    }
    else if(NLMatrix)
    {
        SystematicE = NLE;
    }
    else if(RelativeEnergyScaleMatrix)
    {
        SystematicE = RelativeEnergyE;
    }
    else if(ResolutionMatrix)
    {
        SystematicE = ResolutionE;
    }

    LoadBackgrounds();
    
    std::cout << "Randomize Background #" << BackgroundE << std::endl;
    std::cout << "Randomize Systematic #" << SystematicE << std::endl;
    
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
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                      Initialize histograms, functions, etc.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (Int_t week = 0; week<Nweeks; week++)
    {
        Cov2H[week] = new TH2F(Form("Covariance Matrix%d",week),Form("Covariance Matrix%d",week),MaxBins,0,MaxBins,MaxBins,0,MaxBins);
        CovMatrix2H[week] = new TH2F(Form("Covariance Matrix%d",week),Form("Covariance Matrix%d",week),MaxBins,0,MaxBins,MaxBins,0,MaxBins);
        
        Pred->MakeNominalPrediction(Sin22t13);//  Generate Nominal Prediction
        LoadNominalPredictions(week);
        
        TFile* SaveCanvasAmC1 = TFile::Open("./RootOutputs/AmCDistortionFunctions.root","recreate");//Reset
        SaveCanvasAmC1->Close();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                  Generate toy MC samples
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for (Int_t samples = 0; samples<NSamples; samples++)
        {
            Pred->MakePrediction(Sin22t13);
            LoadAlteredPredictions(week);
            
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                             Generate Covariance Matrix
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            Cov2H[week]->Reset();
            GenerateCovarianceMatrix(week);
            CovMatrix2H[week]->Add(Cov2H[week]);
            
            SaveSpectrum(samples, week);
        }
        
        CovMatrix2H[week]->Scale(1./(NSamples));
        // CovMatrix2H->Draw("colz");
        
        SaveCovarianceMatrix(week);
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                      Clean up the dust
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        delete CovMatrix2H[week];
        delete Cov2H[week];
    }
    delete c;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Produces the Covariance Matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CovarianceMatrix3 :: GenerateCovarianceMatrix(Int_t week)
{
    TFile* SaveVaryACCCanvasF1 = TFile::Open("./RootOutputs/AccidentalVariations.root","recreate");
    SaveVaryACCCanvasF1->Close();
    TFile* SaveVaryLiHeCanvasF1 = TFile::Open("./RootOutputs/LiHeVariations.root","recreate");
    SaveVaryLiHeCanvasF1->Close();
    TFile* SaveVaryFNCanvasF1 = TFile::Open("./RootOutputs/FNVariations.root","recreate");
    SaveVaryFNCanvasF1->Close();
    TFile* SaveVaryAmCCanvasF1 = TFile::Open("./RootOutputs/AmCVariations.root","recreate");
    SaveVaryAmCCanvasF1->Close();
//    TFile* SaveCanvasLiHeF1 = TFile::Open("./RootOutputs/LiHeDistortionFunctions.root","recreate");
//    SaveCanvasLiHeF1->Close();

    for(Int_t AD=0; AD<NADs; AD++)
    {
        BackgroundSpectrumH[AD][week]=(TH1F*)AccidentalsH[0][0]->Clone();
        BackgroundSpectrumH[AD][week]->Reset();
        RandomBackgroundSpectrumH[AD][week]=(TH1F*)BackgroundSpectrumH[AD][week]->Clone();
        
        //  BackgroundSpectrumH holds nominal backgrounds, without being distorted
        //  RandomBackgroundSpectrumH holds random backgrounds, either varying the rate or the shape
        
        FluctuateBackgrounds(week);//Vary backgrounds
        
        if(BackgroundE==VaryAccidentalE)
        {
            BackgroundSpectrumH[AD][week]->Add(AccidentalsH[AD][week]);
            RandomBackgroundSpectrumH[AD][week]->Add(RandomAccidentalsH[AD][week]);
        }
        if (BackgroundE == VaryLiHeE || BackgroundE == DistortLiHeE)
        {
            BackgroundSpectrumH[AD][week]->Add(LiHeH[AD][week]);
            RandomBackgroundSpectrumH[AD][week]->Add(RandomLiHeH[AD][week]);
            
            TFile* SaveVaryLiHeCanvasF = TFile::Open("./RootOutputs/LiHeVariations.root","update");
            TCanvas* VaryLiHeC = new TCanvas("LiHe rate variations");
            
            RandomLiHeH[AD][week]->Draw("same");
            
            VaryLiHeC->Write();
            delete VaryLiHeC;
            SaveVaryLiHeCanvasF->Close();
        }
        if (BackgroundE == VaryFastNeutronsE || BackgroundE == DistortFastNeutronsE)
        {
            BackgroundSpectrumH[AD][week]->Add(FastNeutronsH[AD][week]);
            RandomBackgroundSpectrumH[AD][week]->Add(RandomFastNeutronsH[AD][week]);
            
            TFile* SaveVaryFNCanvasF = TFile::Open("./RootOutputs/FNVariations.root","update");
            TCanvas* VaryFNC = new TCanvas("FN rate variations");
            
            RandomFastNeutronsH[AD][week]->Draw("same");
            
            VaryFNC->Write();
            delete VaryFNC;
            SaveVaryFNCanvasF->Close();

        }
        if (BackgroundE == VaryAmCE || BackgroundE == DistortAmCE)
        {
            BackgroundSpectrumH[AD][week]->Add(AmCH[AD][week]);
            RandomBackgroundSpectrumH[AD][week]->Add(RandomAmCH[AD][week]);
        }
    }
    
    
    //Combine matrices in 9x9, 2x2 or 1x1 prediction
    if (Combine == 1)
    {
        NearBackgroundSpectrumH[0][week]=(TH1F*)BackgroundSpectrumH[0][week]->Clone();
        NearBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[1][week]);//All near hall detectors together
        NearBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[2][week]);
        NearBackgroundSpectrumH[0][week]->Scale(1./3);

        FarBackgroundSpectrumH[0][week]=(TH1F*)BackgroundSpectrumH[3][week]->Clone();
        FarBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[4][week]);//All far hall detectors together
        FarBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[5][week]);
        FarBackgroundSpectrumH[0][week]->Scale(1./3);
        
        NearRandomBackgroundSpectrumH[0][week]=(TH1F*)RandomBackgroundSpectrumH[0][week]->Clone();
        NearRandomBackgroundSpectrumH[0][week]->Add(RandomBackgroundSpectrumH[1][week]);//All near hall detectors together
        NearRandomBackgroundSpectrumH[0][week]->Add(RandomBackgroundSpectrumH[2][week]);
        NearRandomBackgroundSpectrumH[0][week]->Scale(1./3);

        FarRandomBackgroundSpectrumH[0][week]=(TH1F*)RandomBackgroundSpectrumH[3][week]->Clone();
        FarRandomBackgroundSpectrumH[0][week]->Add(RandomBackgroundSpectrumH[4][week]);//All far hall detectors together
        FarRandomBackgroundSpectrumH[0][week]->Add(RandomBackgroundSpectrumH[5][week]);
        FarRandomBackgroundSpectrumH[0][week]->Scale(1./3);

        
    }
    else if(Combine == 2)
    {
        FarBackgroundSpectrumH[0][week]=(TH1F*)BackgroundSpectrumH[3][week]->Clone();//All far hall detectors together
        FarBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[4][week]);
        FarBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[5][week]);
        FarBackgroundSpectrumH[0][week]->Scale(1./3);
        
        NearBackgroundSpectrumH[0][week]=(TH1F*)BackgroundSpectrumH[0][week]->Clone();//DB
        NearBackgroundSpectrumH[0][week]->Add(BackgroundSpectrumH[1][week]);
        NearBackgroundSpectrumH[0][week]->Scale(0.5);
        NearBackgroundSpectrumH[1][week]=(TH1F*)BackgroundSpectrumH[2][week]->Clone();//LO
        
        FarRandomBackgroundSpectrumH[0][week]=(TH1F*)RandomBackgroundSpectrumH[3][week]->Clone();//All far hall detectors together
        FarRandomBackgroundSpectrumH[0][week]->Add(RandomBackgroundSpectrumH[4][week]);
        FarRandomBackgroundSpectrumH[0][week]->Add(RandomBackgroundSpectrumH[5][week]);
        FarRandomBackgroundSpectrumH[0][week]->Scale(1./3);

        NearRandomBackgroundSpectrumH[0][week]=(TH1F*)RandomBackgroundSpectrumH[0][week]->Clone();//DB
        NearRandomBackgroundSpectrumH[0][week]->Add(RandomBackgroundSpectrumH[1][week]);
        NearRandomBackgroundSpectrumH[0][week]->Scale(0.5);
        NearRandomBackgroundSpectrumH[1][week]=(TH1F*)RandomBackgroundSpectrumH[2][week]->Clone();//LO
    }
    else
    {
        NearBackgroundSpectrumH[0][week]=(TH1F*)BackgroundSpectrumH[0][week]->Clone();
        NearBackgroundSpectrumH[1][week]=(TH1F*)BackgroundSpectrumH[1][week]->Clone();
        NearBackgroundSpectrumH[2][week]=(TH1F*)BackgroundSpectrumH[2][week]->Clone();
        FarBackgroundSpectrumH[0][week]=(TH1F*)BackgroundSpectrumH[3][week]->Clone();
        FarBackgroundSpectrumH[1][week]=(TH1F*)BackgroundSpectrumH[4][week]->Clone();
        FarBackgroundSpectrumH[2][week]=(TH1F*)BackgroundSpectrumH[5][week]->Clone();
        
        NearRandomBackgroundSpectrumH[0][week]=(TH1F*)RandomBackgroundSpectrumH[0][week]->Clone();
        NearRandomBackgroundSpectrumH[1][week]=(TH1F*)RandomBackgroundSpectrumH[1][week]->Clone();
        NearRandomBackgroundSpectrumH[2][week]=(TH1F*)RandomBackgroundSpectrumH[2][week]->Clone();
        FarRandomBackgroundSpectrumH[0][week]=(TH1F*)RandomBackgroundSpectrumH[3][week]->Clone();
        FarRandomBackgroundSpectrumH[1][week]=(TH1F*)RandomBackgroundSpectrumH[4][week]->Clone();
        FarRandomBackgroundSpectrumH[2][week]=(TH1F*)RandomBackgroundSpectrumH[5][week]->Clone();
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
            if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3; Fi2=fari+1;}
            if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3; Fi3=fari+1;}
            if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3; Fi4=fari+1;}
            //                std::cout <<"\n"<< "Ni" << Ni1 << Ni2 << Ni3 << Ni4 << "\n";
            //                std::cout << "Fi"<< Fi1 << Fi2 << Fi3 << Fi4 << "\n";
            //                std::cout << "====================================" << icounter++ <<"\n";
            //                std::cout << Ni1*Fi1 << Ni2*Fi2 << Ni3*Fi3 << Ni4*Fi4;
            
            
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
                        
                        for (Int_t j = 0; j<n_evis_bins; j++)
                        {//rows
                            x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                            y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                            
                            Cov[x][y][week]=(AlteredPredictionVisH[fari][neari][week]->GetBinContent(i+1)+FarRandomBackgroundSpectrumH[fari][week]->GetBinContent(i+1)-PredictionVisH[fari][neari][week]->GetBinContent(i+1)-FarBackgroundSpectrumH[fari][week]->GetBinContent(i+1))*(AlteredPredictionVisH[farj][nearj][week]->GetBinContent(j+1)+FarRandomBackgroundSpectrumH[farj][week]->GetBinContent(j+1)-PredictionVisH[farj][nearj][week]->GetBinContent(j+1)-FarBackgroundSpectrumH[farj][week]->GetBinContent(j+1));
                            
                            //                                std::cout << "COV VALUE" << Cov[x][y][week] << std::endl;
                            ////////////////////////////////Normilize systematic and then unnormalize when fitting for different sin2t13. Check also NaNs.
                            if(SystematicE != (CovarianceMatrix3::SystematicType)(-1))
                            {
                                Cov[x][y][week]= Cov[x][y][week]/(PredictionVisH[fari][neari][week]->GetBinContent(i+1)*PredictionVisH[farj][nearj][week]->GetBinContent(j+1));
                            }
                            
                            Cov2H[week]->SetBinContent(x+1,y+1,Cov[x][y][week]);
                        }
                    }
                }
            }
        }
    }
    NormalizeCov(Cov2H[week]);
    //    Cov2H[week]->Draw("colz");
}

void CovarianceMatrix3 :: NormalizeCov(TH2F* Histo)
{
    Int_t x,y;
    
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
            if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3; Fi2=fari+1;}
            if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3; Fi3=fari+1;}
            if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3; Fi4=fari+1;}
            
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
                    if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3; Fj2=farj+1;}
                    if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3; Fj3=farj+1;}
                    if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3; Fj4=farj+1;}
                    
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {//columns
                        
                        for (Int_t j = 0; j<n_evis_bins; j++)
                        {//rows
                            x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
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

void CovarianceMatrix3 :: SaveCovarianceMatrix(Int_t week)
{
    Char_t filenameCov[100];
    
    sprintf(filenameCov,"./CovarianceMatrices/NominalCovarianceMatrix.root");

    //Save Cov Matrix (Either background or systematics)
    switch (BackgroundE)
    {
        case 0://Vary Accidentals
            sprintf(filenameCov,"./CovarianceMatrices/VaryAccidentalCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Vary Accidental Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Vary Accidental Covariance Matrix");
            break;
        case 1://Vary LiHe
            sprintf(filenameCov,"./CovarianceMatrices/VaryLiHeCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Vary LiHe Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Vary LiHe Covariance Matrix");
            break;
        case 2://Vary Fast Neutrons
            sprintf(filenameCov,"./CovarianceMatrices/VaryFastNeutronsCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Vary FN Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Vary FN Covariance Matrix");
            break;
        case 3://Vary AmC
            sprintf(filenameCov,"./CovarianceMatrices/VaryAmCCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Vary AmC Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Vary AmC Covariance Matrix");
            break;
        case 4://Distort LiHe
            sprintf(filenameCov,"./CovarianceMatrices/DistortLiHeCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Distort LiHe Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Distort LiHe Covariance Matrix");
            break;
        case 5://Distort Fast Neutrons
            sprintf(filenameCov,"./CovarianceMatrices/DistortFastNeutronsCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Distort FN Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Distort FN Covariance Matrix");
            break;
        case 6://Distort AmC
            sprintf(filenameCov,"./CovarianceMatrices/DistortAmCCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Distort AmC Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Distort AmC Covariance Matrix");
            break;
    }
    switch (SystematicE)
    {
        case 0://Vary Reactor
            sprintf(filenameCov,"./CovarianceMatrices/IsotopeCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Isotope Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Isotope Covariance Matrix");
            break;
        case 1://Vary Reactor
            sprintf(filenameCov,"./CovarianceMatrices/ReactorPowerCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Reactor Power Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Reactor Power Covariance Matrix");
            break;
        case 2://Vary Energy Scale
            sprintf(filenameCov,"./CovarianceMatrices/RelativeEnergyScaleCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Relative Energy Scale Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Relative Energy Scale Covariance Matrix");
            break;
        case 3:
            sprintf(filenameCov,"./CovarianceMatrices/AbsoluteEnergyScaleCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Absolute Energy Scale Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Absolute Energy Scale Covariance Matrix");
            break;
        case 4:
            sprintf(filenameCov,"./CovarianceMatrices/AbsoluteEnergyOffsetCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Absolute Energy Scale Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Absolute Energy Scale Covariance Matrix");
            break;
        case 5:
            sprintf(filenameCov,"./CovarianceMatrices/AbsoluteEnergyOffsetCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Absolute Energy Scale Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Absolute Energy Scale Covariance Matrix");
            break;
        case 6://Vary IAV
            sprintf(filenameCov,"./CovarianceMatrices/IAVCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("IAV Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("IAV Covariance Matrix");
            break;
        case 7://Vary NL
            sprintf(filenameCov,"./CovarianceMatrices/NLCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("NL Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("NL Covariance Matrix");
            break;
        case 8://Vary Resolution
            sprintf(filenameCov,"./CovarianceMatrices/ResolutionCovarianceMatrix.root");
            CovMatrix2H[week]->SetName(Form("Resolution Covariance Matrix%d",week));
            CovMatrix2H[week]->SetTitle("Resolution Covariance Matrix");
            break;
    }
    TFile* SaveCovarianceMatrixF = TFile::Open(filenameCov,"recreate");
    // CovMatrix2H[week]->Draw("colz");
    CovMatrix2H[week]->Write();
    SaveCovarianceMatrixF->Close();
    
    //Save in a txt file
    if(WriteOutput)//I will have to do the same here, select a covariance matrix name depending on what has been calculated
    {
        ofstream covf("CovarianceMatrices/CovarianceMatrix.txt");
        
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
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3; Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3; Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3; Fi4=fari+1;}
                
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
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3; Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3; Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3; Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            
                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                
                                covf << CovMatrix2H[week]->GetBinContent(x+1,y+1)<< " ";
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
    TFile* Prueba;
    Char_t filenameSpec[100];
    Char_t optionS[10];
    Char_t optionN[10];

    //Save Cov Matrix (Either background or systematics)
    
    switch (BackgroundE)
    {
        case 0://Vary Accidentals
            sprintf(filenameSpec,"./CovarianceMatrices/VaryAccidentalSpectrum.root");
            break;
        case 1://Vary LiHe
            sprintf(filenameSpec,"./CovarianceMatrices/VaryLiHeSpectrum.root");
            break;
        case 2://Vary Fast Neutrons
            sprintf(filenameSpec,"./CovarianceMatrices/VaryFastNeutronsSpectrum.root");
            break;
        case 3://Vary AmC
            sprintf(filenameSpec,"./CovarianceMatrices/VaryAmCSpectrum.root");
            break;
        case 4://Distort LiHe
            sprintf(filenameSpec,"./CovarianceMatrices/DistortLiHeSpectrum.root");
            break;
        case 5://Distort Fast Neutrons
            sprintf(filenameSpec,"./CovarianceMatrices/DistortFastNeutronsSpectrum.root");
            break;
        case 6://Distort AmC
            sprintf(filenameSpec,"./CovarianceMatrices/DistortAmCSpectrum.root");
            break;
        default:
            sprintf(filenameSpec,"./CovarianceMatrices/NominalSpectrum.root");
            break;
    }
    switch (SystematicE)
    {
        case 0://Vary Reactor
            sprintf(filenameSpec,"./CovarianceMatrices/IsotopeSpectrum.root");
            break;
        case 1://Vary Reactor
            sprintf(filenameSpec,"./CovarianceMatrices/PowerSpectrum.root");
            break;
        case 2://Vary Energy Scale
            sprintf(filenameSpec,"./CovarianceMatrices/EnergyScaleSpectrum.root");
            break;
        case 3:
            break;
        case 4:
            break;
        case 5:
            break;
        case 6://Vary IAV
            sprintf(filenameSpec,"./CovarianceMatrices/IAVSpectrum.root");
            break;
        case 7://Vary NL
            sprintf(filenameSpec,"./CovarianceMatrices/NLSpectrum.root");
            break;
        case 8://Vary Resolution
            sprintf(filenameSpec,"./CovarianceMatrices/ResolutionSpectrum.root");
            break;
        default:
//            sprintf(filenameSpec,"./CovarianceMatrices/NominalSpectrum.root");
            break;
    }
    if (sample==0)
    {
        sprintf(optionS,"recreate");
    }
    else
    {
        sprintf(optionS,"update");
    }
    
    TFile* SpectrumAndPredictionsF = TFile::Open(filenameSpec,optionS);

    for (Int_t AD = 0; AD<MaxFar; AD++)
    {
        ADSpectrumVisH[AD][week]->Write(Form("Nominal AD Vis Spectra Sample%i",sample));//
        AlteredADSpectrumVisH[AD][week]->Write(Form("Varied AD Vis Spectra Sample%i",sample));//Form("Near Vis H Varied Sample%i",sample)
        FarBackgroundSpectrumH[AD][week]->Write(Form("Nominal Background Sample%i",sample));
        FarRandomBackgroundSpectrumH[AD][week]->Write(Form("Varied Background Spectra Sample%i",sample));
    }
    
    for (Int_t near = 0; near < MaxNear; near++)
    {
        for (Int_t far =0; far<MaxFar; far++)
        {
            PredictionVisH[far][near][week]->Add(BackgroundSpectrumH[far+ADsEH1+ADsEH2][week],1);
            AlteredPredictionVisH[far][near][week]->Add(RandomBackgroundSpectrumH[far+ADsEH1+ADsEH2][week],1);
            
            PredictionVisH[far][near][week]->Write(Form("Nominal Far VisH Sample%i",sample));
            AlteredPredictionVisH[far][near][week]->Write(Form("Varied Far VisH Sample%i",sample));
            
            //Ratio of the variation
            CopyPredictionVisH[far][near][week]=(TH1F*)PredictionVisH[far][near][week]->Clone("Percentual variation");
            CopyPredictionVisH[far][near][week]->Add(AlteredPredictionVisH[far][near][week],-1);
            CopyPredictionVisH[far][near][week]->Divide(AlteredPredictionVisH[far][near][week]);//(Toy-Nominal)/Toy
            CopyPredictionVisH[far][near][week]->Write(Form("Percentual variation"));
            CopyPredictionVisH[far][near][week]->Draw("E1");
            CopyPredictionVisH[far][near][week]->Draw("C SAME");
        }
    }

    if(sample==5)
    {
        c->Write(); //here write canvas with 5 samples
    }

    SpectrumAndPredictionsF->Close();
    
    flagX=0;

    //Save a sample of the Backgrounds, Nominal Visible Spectrum and Random Visible Spectrum for each systematic.
    if (flagX==0)
    {
        sprintf(optionN,"recreate");
    }
    else
    {
        sprintf(optionN,"update");
    }
    Prueba = TFile::Open(Form("./Varied_And_Nominal_Spectrum_At_Systematic%d_and_Background%d.root",(Int_t)SystematicE,(Int_t)BackgroundE),optionN);

    for (Int_t near = 0; near < MaxNear; near++)
    {
        for (Int_t far =0; far<MaxFar; far++)
        {
            PredictionVisH[far][near][week]->Write();
            AlteredPredictionVisH[far][near][week]->Write();
            BackgroundSpectrumH[far+MaxNear][week]->Write(Form("Nominal Background Sample%i",sample));
            RandomBackgroundSpectrumH[far+MaxNear][week]->Write(Form("Varied Background Spectra Sample%i",sample));
        }
    }
    
    Prueba->Close();
    
    flagX=1;
}

void CovarianceMatrix3 :: LoadNominalPredictions(Int_t week)
{
TFile* PredictionF = TFile::Open("./RootOutputs/NominalPredictedSpectrum.root");
    
    for (Int_t near = 0; near<MaxNear; near++)
    {
        ADSpectrumVisH[near][week]=(TH1F*)gDirectory->Get(Form("Near Prediction AD%d",near));
    }
    for (Int_t near =0; near<MaxNear; near++)
    {
        for (Int_t far =0; far<MaxFar; far++)
        {
            PredictionVisH[far][near][week]=(TH1F*)gDirectory->Get(Form("Far Prediction AD%d from AD%d",far,near));
        }
    }

PredictionF->Close();
}

void CovarianceMatrix3 :: LoadAlteredPredictions(Int_t week)
{
    TFile* PredictionF = TFile::Open("./RootOutputs/PredictedSpectrum.root");

    for (Int_t near = 0; near<MaxNear; near++)
    {
        AlteredADSpectrumVisH[near][week]=(TH1F*)gDirectory->Get(Form("Near Prediction AD%d",near));
    }
    for (Int_t near =0; near<MaxNear; near++)
    {
        for (Int_t far =0; far<MaxFar; far++)
        {
            AlteredPredictionVisH[far][near][week]=(TH1F*)gDirectory->Get(Form("Far Prediction AD%d from AD%d",far,near));
        }
    }
    
    PredictionF->Close();
}

void CovarianceMatrix3 :: LoadBackgrounds()
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
      if(VaryAccidentalMatrix)
    {
        BackgroundE = VaryAccidentalE;
        VaryRate=1;
        Distort=0;
    }
    else if(VaryLiHeMatrix)
    {
        BackgroundE = VaryLiHeE;
        VaryRate=1;
        Distort=0;
    }
    else if(VaryFastNeutronsMatrix)
    {
        BackgroundE = VaryFastNeutronsE;
        VaryRate=1;
        Distort=0;
    }
    else if(VaryAmCMatrix)
    {
        BackgroundE = VaryAmCE;
        VaryRate=1;
        Distort=0;
    }
    else if(DistortLiHeMatrix)
    {
        BackgroundE = DistortLiHeE;
        DistortLiHe=0.2;
        VaryRate=0;
        Distort=1;
    }
    else if(DistortFastNeutronsMatrix)
    {
        BackgroundE = DistortFastNeutronsE;
        DistortFN=0.2;
        VaryRate=0;
        Distort=1;
    }
    else if(DistortAmCMatrix)
    {
        BackgroundE = DistortAmCE;
        DistortAmC=0.2;
        VaryRate=0;
        Distort=1;
    }
}

//Randomize the events (rate) inside errors taking into account the corresponding correlations between backgrounds and ADs. Also random shape variations are included using distortion functions.
void CovarianceMatrix3 :: FluctuateBackgrounds(Int_t week)
{
    for (Int_t AD =0; AD<NADs; AD++)
    {
        ScaleFactorAccidental[AD][week]=1.;
        ScaleFactorLiHe[AD][week]=1.;
        ScaleFactorFastNeutrons[AD][week]=1.;
        ScaleFactorAmC[AD][week]=1.;
        
        hall=2;
        if(AD<ADsEH1)
        {
            hall=1;
        }
        if(AD>=ADsEH1+ADsEH2)
        {
            hall=3;
        }
        
        if(VaryAccidentalMatrix)
        {
            //                rand->SetSeed(0);
            ScaleFactorAccidental[AD][week]=(1.+AccidentalError[AD][week]*rand->Gaus(0,1));
            //                std::cout << rand->Gaus(0,1) << std::endl;
        }
        if(VaryLiHeMatrix)
        {
            //                rand->SetSeed(0);
            ScaleFactorLiHe[AD][week]=(1.+LiHeError[hall][week]*rand->Gaus(0,1));
            //                std::cout << rand->Gaus(0,1) << std::endl;
        }
        if(VaryFastNeutronsMatrix)
        {
            //                rand->SetSeed(0);
            ScaleFactorFastNeutrons[AD][week]=(1.+FastNeutronsError[hall][week]*rand->Gaus(0,1));
            //                std::cout << rand->Gaus(0,1) << std::endl;
        }
        if(VaryAmCMatrix)
        {
            ScaleFactorAmC[AD][week]=(1.+AmCError[AD][week]*rand->Gaus(0,1));
        }
    }
    
    //Accidentals uncorrelated
    if(VaryLiHeMatrix)
    {   //  LiHe correlated by hall
        ScaleFactorLiHe[1][week]=ScaleFactorLiHe[0][week];// EH1
        ScaleFactorLiHe[4][week]=ScaleFactorLiHe[3][week];// EH2
        ScaleFactorLiHe[5][week]=ScaleFactorLiHe[3][week];// EH3
    }
    if(VaryFastNeutronsMatrix)
    {
        //  Fast neutrons by hall
        ScaleFactorFastNeutrons[1][week]=ScaleFactorFastNeutrons[0][week];// EH1
        ScaleFactorFastNeutrons[4][week]=ScaleFactorFastNeutrons[3][week];// EH2
        ScaleFactorFastNeutrons[5][week]=ScaleFactorFastNeutrons[3][week];// EH3
    }
    if(VaryAmCMatrix)
    {   //  AmC all correlated
        ScaleFactorAmC[1][week]=ScaleFactorAmC[0][week];
        ScaleFactorAmC[2][week]=ScaleFactorAmC[0][week];
        ScaleFactorAmC[3][week]=ScaleFactorAmC[0][week];
        ScaleFactorAmC[4][week]=ScaleFactorAmC[0][week];
        ScaleFactorAmC[5][week]=ScaleFactorAmC[0][week];
    }
    for (Int_t AD =0; AD<NADs; AD++)
    {
        if(VaryAccidentalMatrix)
        {
            RandomAccidentalsH[AD][week]=(TH1F*)AccidentalsH[AD][week]->Clone();
            RandomAccidentalsH[AD][week]->Scale(1.*ScaleFactorAccidental[AD][week]);
            
            TFile* SaveVaryACCCanvasF = TFile::Open("./RootOutputs/AccidentalVariations.root","update");
            TCanvas* VaryACCC = new TCanvas("Accidental rate variations");
            
            RandomAccidentalsH[AD][week]->Draw("same");
            
            VaryACCC->Write();
            delete VaryACCC;
            SaveVaryACCCanvasF->Close();
        }
        if(VaryLiHeMatrix)
        {
            RandomLiHeH[AD][week]=(TH1F*)LiHeH[AD][week]->Clone();
            RandomLiHeH[AD][week]->Scale(1.*ScaleFactorLiHe[AD][week]);
            
            TFile* SaveVaryLiHeCanvasF = TFile::Open("./RootOutputs/LiHeVariations.root","update");
            TCanvas* VaryLiHeC = new TCanvas("LiHe rate variations");
            
            RandomLiHeH[AD][week]->Draw("same");
            
            VaryLiHeC->Write();
            delete VaryLiHeC;
            SaveVaryLiHeCanvasF->Close();
        }
        if(VaryFastNeutronsMatrix)
        {
            RandomFastNeutronsH[AD][week]=(TH1F*)FastNeutronsH[AD][week]->Clone();
            RandomFastNeutronsH[AD][week]->Scale(1.*ScaleFactorFastNeutrons[AD][week]);
            
            TFile* SaveVaryFNCanvasF = TFile::Open("./RootOutputs/FNVariations.root","update");
            TCanvas* VaryFNC = new TCanvas("FN rate variations");
            
            RandomFastNeutronsH[AD][week]->Draw("same");
            
            VaryFNC->Write();
            delete VaryFNC;
            SaveVaryFNCanvasF->Close();
        }
        if(VaryAmCMatrix)
        {
            RandomAmCH[AD][week]=(TH1F*)AmCH[AD][week]->Clone();
            RandomAmCH[AD][week]->Scale(1.*ScaleFactorAmC[AD][week]);
            
            TFile* SaveVaryAmCCanvasF = TFile::Open("./RootOutputs/AmCVariations.root","update");
            TCanvas* VaryAmCC = new TCanvas("AmC rate variations");
            
            RandomAmCH[AD][week]->Draw("same");
            
            VaryAmCC->Write();
            delete VaryAmCC;
            SaveVaryAmCCanvasF->Close();
        }
        //                std::cout << ScaleFactorAccidental[AD][week] << " " << ScaleFactorLiHe[AD][week] <<" " <<  ScaleFactorFastNeutrons[AD][week] <<" " <<  ScaleFactorAmC[AD][week] << std::endl;
        //                RandomAmCH[AD][week]->Draw("same");
        //                RandomLiHeH[AD][week]->Draw("same");
        //                RandomFastNeutronsH[AD][week]->Draw("same");
        //                RandomAccidentalsH[AD][week]->Draw("same");
        
    }
    
    if(Distort)
    {
        if(DistortLiHeMatrix)//Distort LiHe shape for all ADs in the same way
        {
            TFile* SaveCanvasLiHe = TFile::Open("./RootOutputs/LiHeDistortions.root","recreate");
            TCanvas* LiHeC = new TCanvas("LiHe Distortions");
            TF1* func_LiHe = GetDistortionFunction(DistortLiHe);

            for(Int_t AD=0;AD<NADs;AD++)
            {
                RandomLiHeH[AD][week]=(TH1F*)LiHeH[AD][week]->Clone();
                RandomLiHeH[AD][week]->Multiply(func_LiHe);
                RandomLiHeH[AD][week]->Scale(LiHeH[AD][week]->Integral()/RandomLiHeH[AD][week]->Integral());
                RandomLiHeH[AD][week]->Draw("same");
            }
            LiHeC->Write();
            delete LiHeC;
            delete func_LiHe;
            
            SaveCanvasLiHe->Close();
        }
        if(DistortFastNeutronsMatrix)//Distort FN shape for all ADs in the same hall
        {
            TFile* SaveFNCanvas = TFile::Open("./RootOutputs/FNDistortions.root","recreate");
            TCanvas* FNC = new TCanvas("FN Distortions");
            
            //FN shape distortions are applied taking into account correlations in each EH
            TF1* func_FN = GetFastNeutronsDistortionFunction(DistortFN);
            
            for(Int_t AD=0;AD<ADsEH1;AD++)
            {
                RandomFastNeutronsH[AD][week]=(TH1F*)FastNeutronsH[AD][week]->Clone();
                RandomFastNeutronsH[AD][week]->Multiply(func_FN);
                RandomFastNeutronsH[AD][week]->Scale(RandomFastNeutronsH[AD][week]->Integral()/FastNeutronsH[AD][week]->Integral());
                RandomFastNeutronsH[AD][week]->Draw("same");
            }
            func_FN=GetFastNeutronsDistortionFunction(DistortFN);
            
            for(Int_t AD=ADsEH1;AD<(ADsEH1+ADsEH2);AD++)
            {
                RandomFastNeutronsH[AD][week]=(TH1F*)FastNeutronsH[AD][week]->Clone();
                RandomFastNeutronsH[AD][week]->Multiply(func_FN);
                RandomFastNeutronsH[AD][week]->Scale(RandomFastNeutronsH[AD][week]->Integral()/FastNeutronsH[AD][week]->Integral());
                RandomFastNeutronsH[AD][week]->Draw("same");
            }
            func_FN=GetFastNeutronsDistortionFunction(DistortFN);
            
            for(Int_t AD=(ADsEH1+ADsEH2);AD<(ADsEH1+ADsEH2+ADsEH3);AD++)
            {
                RandomFastNeutronsH[AD][week]=(TH1F*)FastNeutronsH[AD][week]->Clone();
                RandomFastNeutronsH[AD][week]->Multiply(func_FN);
                RandomFastNeutronsH[AD][week]->Scale(RandomFastNeutronsH[AD][week]->Integral()/FastNeutronsH[AD][week]->Integral());
                RandomFastNeutronsH[AD][week]->Draw("same");
            }
            
            FNC->Write();
            delete FNC;
            delete func_FN;
            SaveFNCanvas->Close();
        }
        if(DistortAmCMatrix)//Distort AmC shape for all ADs in the same way
        {
            TFile* SaveCanvasAmC = TFile::Open("./RootOutputs/AmCDistortions.root","recreate");
            TCanvas* AMCC = new TCanvas("AmC Distortions");
            
            TF1* func_AmC = GetDistortionFunction(DistortAmC);
            
            for(Int_t AD=0;AD<NADs;AD++)
            {
                RandomAmCH[AD][week]=(TH1F*)AmCH[AD][week]->Clone();
                RandomAmCH[AD][week]->Multiply(func_AmC);
                RandomAmCH[AD][week]->Scale(AmCH[AD][week]->Integral()/RandomAmCH[AD][week]->Integral());
                RandomAmCH[AD][week]->Draw("same");
            }
            AMCC->Write();
            delete AMCC;
            delete func_AmC;
            SaveCanvasAmC->Close();
        }
    }
//    for (Int_t AD =0; AD<NADs; AD++)
//    {
//        ScaleFactorAccidental[AD][week]=1.;
//        ScaleFactorLiHe[AD][week]=1.;
//        ScaleFactorFastNeutrons[AD][week]=1.;
//        ScaleFactorAmC[AD][week]=1.;
//
//        if(VaryAccidentalMatrix)
//        {
//            //                rand->SetSeed(0);
//            ScaleFactorAccidental[AD][week]=(1.+(AccidentalError[AD][week]*rand->Gaus(0,1)));
//            std::cout << "Accidental error " << AccidentalError[AD][week] << std::endl;
//        }
//        else if(VaryLiHeMatrix)
//        {
//            //                rand->SetSeed(0);
//            ScaleFactorLiHe[AD][week]=(1.+(LiHeError[AD][week]*rand->Gaus(0,1)));
//            std::cout <<"LiHe error " << LiHeError[AD][week] << std::endl;
//        }
//        else if(VaryFastNeutronsMatrix)
//        {
//            //                rand->SetSeed(0);
//            ScaleFactorFastNeutrons[AD][week]=(1.+(FastNeutronsError[AD][week]*rand->Gaus(0,1)));
//            std::cout <<"Fast Neutron error " << FastNeutronsError[AD][week] << std::endl;
//        }
//        if(VaryAmCMatrix)
//        {
//            //                rand->SetSeed(0);
//            ScaleFactorAmC[AD][week]=(1.+(AmCError[AD][week]*rand->Gaus(0,1)));
//            std::cout <<"AmC error " << AmCError[AD][week] << std::endl;
//        }
//    }
//    
//    //Accidentals uncorrelated
//    if(VaryLiHeMatrix)
//    {   //  LiHe correlated by hall
//        ScaleFactorLiHe[1][week]=ScaleFactorLiHe[0][week];// EH1
//        ScaleFactorLiHe[4][week]=ScaleFactorLiHe[3][week];// EH3
//        ScaleFactorLiHe[5][week]=ScaleFactorLiHe[3][week];// EH3
//    }
//    else if(VaryFastNeutronsMatrix)
//    {
//        //  Fast neutrons by hall
//        ScaleFactorFastNeutrons[1][week]=ScaleFactorFastNeutrons[0][week];// EH1
//        ScaleFactorFastNeutrons[4][week]=ScaleFactorFastNeutrons[3][week];// EH3
//        ScaleFactorFastNeutrons[5][week]=ScaleFactorFastNeutrons[3][week];// EH3
//    }
//    else if(VaryAmCMatrix)
//    {   //  AmC all correlated
//        ScaleFactorAmC[1][week]=ScaleFactorAmC[0][week];
//        ScaleFactorAmC[2][week]=ScaleFactorAmC[0][week];
//        ScaleFactorAmC[3][week]=ScaleFactorAmC[0][week];
//        ScaleFactorAmC[4][week]=ScaleFactorAmC[0][week];
//        ScaleFactorAmC[5][week]=ScaleFactorAmC[0][week];
//    }
//    for (Int_t AD =0; AD<NADs; AD++)
//    {
//        RandomAccidentalsH[AD][week]=(TH1F*)AccidentalsH[AD][week]->Clone();
//        RandomLiHeH[AD][week]=(TH1F*)LiHeH[AD][week]->Clone();
//        RandomFastNeutronsH[AD][week]=(TH1F*)FastNeutronsH[AD][week]->Clone();
//        RandomAmCH[AD][week]=(TH1F*)AmCH[AD][week]->Clone();
//        
//        if(VaryAccidentalMatrix)
//        {
//            RandomAccidentalsH[AD][week]->Scale(ScaleFactorAccidental[AD][week]);
//            
//            TFile* SaveVaryACCCanvasF = TFile::Open("./RootOutputs/AccidentalVariations.root","recreate");
//                TCanvas* VaryACCC = new TCanvas("Accidental rate variations");
//            
//                RandomAccidentalsH[AD][week]->Draw("same");
//            
//                VaryACCC->Write();
//                delete VaryACCC;
//            SaveVaryACCCanvasF->Close();
//        }
//        else if(VaryLiHeMatrix)
//        {
//            RandomLiHeH[AD][week]->Scale(ScaleFactorLiHe[AD][week]);
//            
//            TFile* SaveVaryLiHeCanvasF = TFile::Open("./RootOutputs/LiHeVariations.root","recreate");
//                TCanvas* VaryLiHeC = new TCanvas("LiHe rate variations");
//            
//                RandomLiHeH[AD][week]->Draw("same");
//            
//                VaryLiHeC->Write();
//                delete VaryLiHeC;
//            SaveVaryLiHeCanvasF->Close();
//        }
//        else if(VaryFastNeutronsMatrix)
//        {
//            RandomFastNeutronsH[AD][week]->Scale(ScaleFactorFastNeutrons[AD][week]);
//            
//            TFile* SaveVaryFNCanvasF = TFile::Open("./RootOutputs/FNVariations.root","recreate");
//                TCanvas* VaryFNC = new TCanvas("FN rate variations");
//            
//                RandomFastNeutronsH[AD][week]->Draw("same");
//            
//                VaryFNC->Write();
//                delete VaryFNC;
//            SaveVaryFNCanvasF->Close();
//        }
//        if(VaryAmCMatrix)
//        {
//            RandomAmCH[AD][week]->Scale(ScaleFactorAmC[AD][week]);
//            
//            TFile* SaveVaryAmCCanvasF = TFile::Open("./RootOutputs/AmCVariations.root","recreate");
//                TCanvas* VaryAmCC = new TCanvas("AmC rate variations");
//            
//                RandomAmCH[AD][week]->Draw("same");
//            
//                VaryAmCC->Write();
//                delete VaryAmCC;
//            SaveVaryAmCCanvasF->Close();
//        }
//        //                std::cout << ScaleFactorAccidental[AD][week] << " " << ScaleFactorLiHe[AD][week] <<" " <<  ScaleFactorFastNeutrons[AD][week] <<" " <<  ScaleFactorAmC[AD][week] << std::endl;
//    }
//    
//    if(DistortLiHeMatrix)// Distort LiHe shape for all ADs in the same way
//    {
//
//        TF1* func_LiHe = GetDistortionFunction(DistortLiHe);
//        TFile* SaveCanvasLiHeF = TFile::Open("./RootOutputs/LiHeDistortionFunctions.root","update");
//        TCanvas* LiHeFunctions = new TCanvas("LiHe Distortion Functions");
//            func_LiHe->Draw("same");
//            LiHeFunctions->Write();
//        delete LiHeFunctions;
//        SaveCanvasLiHeF->Close();
//        
//        TFile* SaveCanvasLiHe = TFile::Open("./RootOutputs/LiHeDistortions.root","recreate");
//        TCanvas* LiHeC = new TCanvas("LiHe Distortions");
//        
//        for(Int_t AD=0;AD<NADs;AD++)
//        {
//            RandomLiHeH[AD][week]->Multiply(func_LiHe);
//            RandomLiHeH[AD][week]->Scale(LiHeH[AD][week]->Integral()/RandomLiHeH[AD][week]->Integral());
//            RandomLiHeH[AD][week]->Draw("same");
//        }
//        LiHeC->Write();
//        delete LiHeC;
//        delete func_LiHe;
//        
//        SaveCanvasLiHe->Close();
//        
//    }
//    else if(DistortFastNeutronsMatrix)//Distort FN shape for all ADs in the same hall
//    {
//        TFile* SaveFNCanvas = TFile::Open("./RootOutputs/FNDistortions.root","recreate");
//        TCanvas* FNC = new TCanvas("FN Distortions");
//        
//        //  FN shape distortions are applied taking into account correlations in each EH
//        TF1* func_FN = GetFastNeutronsDistortionFunction(DistortFN);
//
//        for(Int_t AD=0;AD<ADsEH1;AD++)
//        {
//            RandomFastNeutronsH[AD][week]->Multiply(func_FN);
//            RandomFastNeutronsH[AD][week]->Scale(FastNeutronsH[AD][week]->Integral()/RandomFastNeutronsH[AD][week]->Integral());
//            RandomFastNeutronsH[AD][week]->Draw("same");
//        }
//        func_FN=GetFastNeutronsDistortionFunction(DistortFN);
//        
//        for(Int_t AD=ADsEH1;AD<(ADsEH1+ADsEH2);AD++)
//        {
//            RandomFastNeutronsH[AD][week]->Multiply(func_FN);
//            RandomFastNeutronsH[AD][week]->Scale(FastNeutronsH[AD][week]->Integral()/RandomFastNeutronsH[AD][week]->Integral());
//            RandomFastNeutronsH[AD][week]->Draw("same");
//        }
//        func_FN=GetFastNeutronsDistortionFunction(DistortFN);
//        
//        for(Int_t AD=(ADsEH1+ADsEH2);AD<(ADsEH1+ADsEH2+ADsEH3);AD++)
//        {
//            RandomFastNeutronsH[AD][week]->Multiply(func_FN);
//            RandomFastNeutronsH[AD][week]->Scale(FastNeutronsH[AD][week]->Integral()/RandomFastNeutronsH[AD][week]->Integral());
//            RandomFastNeutronsH[AD][week]->Draw("same");
//        }
//        FNC->Write();
//        delete FNC;
//        
//        delete func_FN;
//        SaveFNCanvas->Close();
//    }
//    else if(DistortAmCMatrix)//  Distort AmC shape for all ADs in the same way
//    {
//        TFile* SaveCanvasAmC = TFile::Open("./RootOutputs/AmCDistortions.root","recreate");
//        TCanvas* AMCC = new TCanvas("AmC Distortions");
//        
//        TF1* func_AmC = GetDistortionFunction(DistortAmC);
//        
//        for(Int_t AD=0;AD<NADs;AD++)
//        {
//            RandomAmCH[AD][week]->Multiply(func_AmC);
//            RandomAmCH[AD][week]->Scale(AmCH[AD][week]->Integral()/RandomAmCH[AD][week]->Integral());
//            RandomAmCH[AD][week]->Draw("same");
//        }
//        AMCC->Write();
//        delete AMCC;
//        delete func_AmC;
//        SaveCanvasAmC->Close();
//    }
//    
//    for(Int_t AD=0; AD<NADs; AD++)
//    {
//        //  RandomBackgroundSpectrumH holds random backgrounds, either varying the rate or the shape
//            RandomBackgroundSpectrumH[AD][week]=(TH1F*)BackgroundSpectrumH[AD][week]->Clone();
//            RandomBackgroundSpectrumH[AD][week]->Reset();
//            RandomBackgroundSpectrumH[AD][week]->Add(RandomAccidentalsH[AD][week]);
//            RandomBackgroundSpectrumH[AD][week]->Add(RandomLiHeH[AD][week]);
//            RandomBackgroundSpectrumH[AD][week]->Add(RandomFastNeutronsH[AD][week]);
//            RandomBackgroundSpectrumH[AD][week]->Add(RandomAmCH[AD][week]);
//    }
}

TF1* CovarianceMatrix3 :: GetDistortionFunction(Double_t amount)
{
    TF1 *func = new TF1("func","TMath::Abs([0]+[1]*x)",InitialEnergy,FinalVisibleEnergy);
    rand->SetSeed(0);
    Double_t slope=amount*rand->Gaus(0,1);
    Double_t anchor_point=3.5;
    //want offset to be set by requiring func at anchor point to be 1
    Double_t offset=(1-slope*anchor_point);
    func->SetParameter(0,offset);
    func->SetParameter(1,slope);

    return func;

}

TF1* CovarianceMatrix3 :: GetFastNeutronsDistortionFunction(Double_t amount)
{
    TF1 *func = new TF1("func","[0]/(0.2*pow(x,0.1))+[1]",InitialEnergy,FinalVisibleEnergy);
    rand->SetSeed(0);
    Double_t scaling =amount*rand->Gaus(0,1);
    func->SetParameter(0,scaling);
    // func->SetParameter(1,0); // see if it works
    //set offset so that func(FinalEnergy)=1;
    func->SetParameter(1,1-1*func->Eval(10));//I'm using 10 instead of 12 because when I fit the FN I use the range 10-100 MeV.
    //    func->Draw();
    return func;
}
