#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2.h"
#include "TMath.h"
#include "NominalData.h"
#include <math.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TCanvas.h"
#include "CrossSection.h"
#include "ReactorSpectrumMultiple.h"
#include "AntineutrinoSpectrum.h"
#include "Oscillation.h"
#include "OscillationReactor.h"
#include "FitBackgrounds2.h"
#include "TAxis.h"
#include "TArrayD.h"
#include "TTree.h"

#define UseChristineReactorModel
//#define Produce_Antineutrino_Spectrum_For_FirstTime
const bool WriteROOT = 1;
const bool ReadTxt = 0;//To use txt matrices or root files.
const bool WriteOutput=0;//To save the covariance matrices in a .txt file.

const bool PlotCorrelation = 1;

//  To test the LBNL Covariance Matrices:
const bool Rate =0;//Provisional, need to study how to do the rate analysis in the fitter.

const Int_t MaxSystematics =9;//(8)Systematics + Total Systematic

class Prediction
{
private:
    bool firstNominalPrediction;
    bool firstRandomPrediction;
    Double_t L[MaxSystematics][9*MaxNbins][9*MaxNbins];
    Double_t RenormToyMCSample[9*MaxNbins][9*MaxNbins];
    Double_t Correlation[9*MaxNbins][9*MaxNbins];
    
    bool flag_delete_ToyMCSample;
    TH1D* VariationHistoH[MaxSystematics];
    
    TTree *T;
    std::string AnalysisString;
    std::string RandomString;
    
    TH1D* TNominalHisto;
    TH1D* TDataHisto;
    TH1D* TVariationHisto;
    TH1D* TNominalHisto1;
    TH1D* TDataHisto1;
    TH1D* TVariationHisto1;
    
    NominalData* Data;
    TBenchmark* gBenchmark;
    
    bool ToyMCflag;
    
    void SetRandomS22t12(bool);
    void SetRandomIsotopeFraction(bool);
    void SetRandomReactorPower(bool);
    
    Int_t DataSet;
    
    bool BudgetTurnOff;
    bool BudgetTurnOn;
    bool RandomSin22t12;
    bool RandomIsotopeFraction;
    bool RandomReactorPower;
    
    std::string BkgCovDirectory;
    std::string SysCovDirectory;
    
    Double_t m_rel_escale[MaxDetectors];
    
    bool StatisticalFluctuation;
    
    bool IsotopeMatrix;
    bool ReactorPowerMatrix;
    bool Sin22t12Matrix;
    
    bool UseToyMCTree;
    
    Int_t NSteps;
    
    //Cell parameters:
    Int_t XCellLimit;
    Int_t YCellLimit;
    
    
    TRandom3* rand;
    OscillationReactor* OscRea;
    Oscillation* Osc;
    Char_t optionN[10];
    Char_t optionS[10];
    Char_t optionC[10];
    
    Double_t Sin22t13;//    Nominal
    Double_t DM2_ee;//      Nominal
    
    Double_t s22t13;
    Double_t Dm2_ee;
    
    Double_t s22t13start;
    Double_t s22t13end;
    Double_t dm2_eestart;
    Double_t dm2_eeend;
    
    Double_t SinWidth;
    Double_t DeltaWidth;
    
    Double_t Tsin22t13;
    Double_t Tdm2_ee;
    Int_t TExperiment;
    Int_t Experiment;
    Int_t NExperiments;
    
    Double_t s22t12;
    Double_t s22t12Nominal;
    Double_t s22t12Error;
    Double_t RandomSin22t12Error;
    
    Int_t MaxBins;
    
    //AD configuration parameters:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    
    Int_t Combine;
    Int_t MaxFarCombine;
    Int_t MaxNearCombine;
    Int_t MaxNearLoadOscModel;
    Int_t MaxFarLoadOscModel;
    
    Int_t hall;
    bool analysis;
    
    //Binning parameters:
    Int_t n_evis_bins;
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.

    Int_t Nweeks;
    Int_t TotalBins;
    
    TH1D* FluxHisto[NReactors][MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    
    //Histograms
    TH1D* PredictionVisH[MaxFarDetectors][MaxNearDetectors];//  Prediction with no alterations and no backgrounds.
    TH1D* ReactorPredictionVisH[MaxDetectors];//  Prediction with no alterations but nominal backgrounds have been added. Vis Binning (after detector response is applied)


    TH1D* CombinedReactorPredictionVisH[MaxFarDetectors];
    TH1D* CombinedPredictionVisH[MaxFarDetectors][MaxNearDetectors];
    TH1D* PredictionH[MaxFarDetectors][MaxNearDetectors];
    TH1D* PredictionDataH[MaxFarDetectors][MaxNearDetectors];
    
    TH1D* NearBackgroundSpectrumH[MaxNearDetectors];
    TH1D* FarBackgroundSpectrumH[MaxFarDetectors];
    TH1D* CombinedNearBackgroundSpectrumH[MaxNearDetectors];
    TH1D* CombinedFarBackgroundSpectrumH[MaxFarDetectors];
    
    // Data and prediction histograms
    
    TH1D* FarDataH[MaxFarDetectors][MaxNearDetectors];//  Prediction with no alterations, nominal backgrounds may be added. Vis Binning (after detector response is applied)
    TH1D* NearDataH[MaxNearDetectors];
    TH1D* CombinedNearDataH[MaxNearDetectors];
    TH1D* CopyFarDataH[MaxFarDetectors][MaxNearDetectors];
    TH1D* ToyFarDataH[MaxFarDetectors][MaxNearDetectors];
    
    //  Root Histograms
    TH2D* DAmCCovarianceMatrixH;
    TH2D* VAmCCovarianceMatrixH;
    TH2D* DFNCovarianceMatrixH;
    TH2D* VFNCovarianceMatrixH;
    TH2D* DLiHeCovarianceMatrixH;
    TH2D* VLiHeCovarianceMatrixH;
    TH2D* VAccCovarianceMatrixH;
    
    TH2D* IsotopeCovarianceMatrixH;
    TH2D* ReactorPowerCovarianceMatrixH;
    TH2D* RelativeEnergyScaleCovarianceMatrixH;
    //    TH2D* RelativeEnergyOffsetCovarianceMatrixH;
    //    TH2D* AbsoluteEnergyScaleCovarianceMatrixH;
    //    TH2D* AbsoluteEnergyOffsetCovarianceMatrixH;
    TH2D* IAVCovarianceMatrixH;
    TH2D* NLCovarianceMatrixH;
    TH2D* ResolutionCovarianceMatrixH;
    TH2D* Sin22t12CovarianceMatrixH;
    TH2D* EfficiencyCovarianceMatrixH;
    
    TH2D* RenormIsotopeCovarianceMatrixH;
    TH2D* RenormReactorPowerCovarianceMatrixH;
    TH2D* RenormRelativeEnergyScaleCovarianceMatrixH;
    //    TH2D* RenormRelativeEnergyOffsetCovarianceMatrixH;
    //    TH2D* RenormAbsoluteEnergyScaleCovarianceMatrixH;
    //    TH2D* RenormAbsoluteEnergyOffsetCovarianceMatrixH;
    TH2D* RenormIAVCovarianceMatrixH;
    TH2D* RenormNLCovarianceMatrixH;
    TH2D* RenormResolutionCovarianceMatrixH;
    TH2D* RenormSin22t12CovarianceMatrixH;
    TH2D* RenormEfficiencyCovarianceMatrixH;
    
    TH2D* StatisticalCovarianceMatrixH;
    TH2D* BackgroundsCovarianceMatrixH;
    TH2D* SystematicCovarianceMatrixH;
    
    TH2D* TotalCovarianceMatrixH;
    TH2D* InvTotalCovarianceMatrixH;
    
    TH2D* UnityH;
    //  Matrices
    Double_t DAmCCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t VAmCCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t DFNCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t VFNCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t DLiHeCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t VLiHeCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t VAccCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    
    Double_t IsotopeCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t ReactorPowerCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t RelativeEnergyScaleCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    //    Double_t AbsoluteEnergyScaleCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    //    Double_t RelativeEnergyOffsetCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    //    Double_t AbsoluteEnergyOffsetCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t IAVCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t NLCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t ResolutionCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t Sin22t12CovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t EfficiencyCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    
    Double_t RenormIsotopeCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t RenormReactorPowerCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t RenormRelativeEnergyScaleCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    //    Double_t RenormRelativeEnergyOffsetCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    //    Double_t RenormAbsoluteEnergyScaleCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    //    Double_t RenormAbsoluteEnergyOffsetCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t RenormIAVCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t RenormNLCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t RenormResolutionCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t RenormSin22t12CovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t RenormEfficiencyCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    
    Double_t BackgroundsCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t SystematicCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t InvTotalCovarianceMatrixM[9*MaxNbins][9*MaxNbins];
    Double_t CovStat[9*MaxNbins][9*MaxNbins];
    //
    //    Double_t UnityM[9*MaxNbins][9*MaxNbins];
    
    Double_t TotalCovarianceMatrixM[9*MaxNbins*9*MaxNbins];
    
    void GetOscillationPrediction(Int_t,Int_t);
    void GetReactorOscillationPrediction(Int_t,Int_t);
    
    void SetSystematic();
    
    void GenerateStatisticalCovarianceMatrix();
    void SaveStatisticalCovarianceMatrix(Int_t,Double_t,Double_t);
    void SaveChi2Components(Int_t,Double_t,Double_t);
    void InvertMatrix(Int_t);
    void CombineMatrices(Int_t);
    void LoadBackgrounds(Int_t,bool);
    
    TH2D* NormCov(TH2D*,TH2D*);
    void ApplyStatisticalFluctuation(TH1D*);
    void SaveCovarianceMatrices(Int_t);
    Double_t ChiSquare(Int_t);
    Double_t RateChiSquare(Int_t);
    void GenerateFluxCorrectedHistograms(NominalData*);
    void AddBackgroundsToReactorPrediction();
public:
    Prediction();
    Prediction(NominalData*);
    ~Prediction();
    
    NominalData* GetNominalData();
    void SetNominalData(NominalData*);
    
    void LoadData(Int_t,bool,Int_t,bool);
    void LoadRootCovarianceMatrices(Int_t);
    void LoadTxtCovarianceMatrices(Int_t);
    void DeleteData();
    void DeleteMatrices();
    void SetSin22t13(Double_t);
    void SetDM213(Double_t);
    void MakePrediction(Double_t,Double_t,bool,Int_t,bool,bool);
    
    TH1D* GetPrediction(Int_t,Int_t);
    TH1D* GetReactorPrediction(Int_t);
    
    TH1D* GetToyMCSample(Int_t);
    void DeletePrediction(Int_t,Int_t);
    
    Double_t CalculateChi2(Double_t,Double_t,Int_t,bool);
    void GenerateInverseMatrix(Double_t,Double_t,Int_t,Int_t,bool,Int_t);
    void GenerateToyMCTree(Int_t,Int_t,NominalData*);
    void ProduceCovToyMCSample(Int_t,TH1D**);
    
    void SetExperiment(Int_t);
};

Prediction :: ~Prediction()
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                               Clean up the dust
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(ToyMCflag&&!UseToyMCTree)
    {
        for(Int_t AD = 0; AD<NADs; AD++)
        {
            delete ReactorPredictionVisH[AD];
        }
    }
    
    
    for(Int_t week =0;week<Nweeks;week++)
    {
        for(Int_t AD =0;AD<NADs;AD++)
        {
            for(Int_t reactor =0;reactor<NReactors;reactor++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        delete FluxHisto[reactor][AD][week][idx][idy];
                    }
                }
            }
        }
    }
    if(flag_delete_ToyMCSample)
    {
        for(Int_t i = 0; i<MaxSystematics; i++)
        {
            delete VariationHistoH[i];
        }
    }
    delete rand;
    delete Data;
}

Prediction :: Prediction()
{
    firstNominalPrediction=0;
    firstRandomPrediction=0;
    std::cout << " the prediction default constructor shouldn't be called, except for Minuit?" << std::endl;
    
    exit(EXIT_FAILURE);
    
    Experiment = 0;
    
    Data = new NominalData(0,2);
    rand = new TRandom3(0);
    gBenchmark = new TBenchmark();
    
    NSteps = Data->GetNSteps();
    
    BudgetTurnOn = Data->GetTurnOnBudget();
    BudgetTurnOff = Data-> GetTurnOffBudget();
    
    Sin22t13 = Data->GetSin22t13();
    DM2_ee = Data->GetDm2ee();
    s22t12Nominal = Data->GetSin22t12();
    s22t12Error = Data->GetSin22t12Error()/s22t12Nominal;//    Relative error = AbsoluteError/NominalValue
    
    s22t13start= Data->GetSinStart();
    s22t13end=Data->GetSinEnd();
    dm2_eestart=Data->GetDmeeStart();
    dm2_eeend=Data->GetDmeeEnd();
    
    SinWidth = (s22t13end-s22t13start)/(NSteps-1);
    DeltaWidth = (dm2_eeend-dm2_eestart)/(NSteps-1);
    
    StatisticalFluctuation = Data->GetStatisticalFluctuation();
    
    analysis = Data->GetAnalysis();
    
    if(analysis)
    {
        AnalysisString = "Hydrogen";
        XCellLimit = VolumeX;
        YCellLimit = VolumeY;
    }
    else
    {
        AnalysisString = "Gadolinium";
        XCellLimit = 1;
        YCellLimit = 1;
    }
    
    Combine = Data->GetCombineMode();
    DataSet = Data->GetDataSet();
    
    IsotopeMatrix = Data->GetIsotopeMatrix();
    ReactorPowerMatrix = Data->GetReactorPowerMatrix();
    Sin22t12Matrix = Data->GetSin22t12Matrix();
    
    //Randomize parameters: 1 to active, 0 to desactivate
    RandomSin22t12 = 0;
    RandomIsotopeFraction = 0;
    RandomReactorPower = 0;
    
    Nweeks = Data->GetWeeks();
    
    TotalBins = MatrixBins;
    
    n_evis_bins = Data->GetVisibleBins();
    
    for (Int_t i = 0; i <= n_evis_bins; i++)
    {
        evis_bins[i] = Data->GetVisibleBinningArray(i);
    }
    
    NADs = Data->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    if(NADs == 8)// ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
#ifndef Produce_Antineutrino_Spectrum_For_FirstTime
    if(IsotopeMatrix||ReactorPowerMatrix)
    {
#endif
        
#ifdef UseChristineReactorModel
        Data->ReadChristineReactorSpectrum();
        
        if(IsotopeMatrix)
        {
            Data->ReadChristineCovMatrix();
        }
#else
        Data->ReadIHEPReactorSpectrum();
        
        if(IsotopeMatrix)
        {
        //Needed a covariance matrix or some way to vary the IHEP spectrum
        }
#endif

#ifndef Produce_Antineutrino_Spectrum_For_FirstTime
    }
#endif
    UseToyMCTree = Data->GetUseToyMCTree();
    
    BkgCovDirectory = Data->GetBkgCovDirectory();
    SysCovDirectory = Data->GetSysCovDirectory();
    
    GenerateFluxCorrectedHistograms(Data);

}

Prediction :: Prediction(NominalData* data)
{
    firstNominalPrediction=0;
    firstRandomPrediction=0;
    
    Data = new NominalData(data->GetAnalysis(),data->GetDataSet());
    
    Data->CopyData(data);
    
    rand = new TRandom3(0);
    gBenchmark = new TBenchmark();
    
    Experiment = 0;

    NSteps = Data->GetNSteps();

    BudgetTurnOn = Data->GetTurnOnBudget();
    BudgetTurnOff = Data-> GetTurnOffBudget();
    
    Sin22t13 = Data->GetSin22t13();
    DM2_ee = Data->GetDm2ee();
    s22t12Nominal = Data->GetSin22t12();
    s22t12Error = Data->GetSin22t12Error()/s22t12Nominal;//    Relative error = AbsoluteError/NominalValue
    
    s22t13start= Data->GetSinStart();
    s22t13end=Data->GetSinEnd();
    dm2_eestart=Data->GetDmeeStart();
    dm2_eeend=Data->GetDmeeEnd();
    
    SinWidth = (s22t13end-s22t13start)/(NSteps-1);
    DeltaWidth = (dm2_eeend-dm2_eestart)/(NSteps-1);
    
    StatisticalFluctuation = Data->GetStatisticalFluctuation();
    
    analysis = Data->GetAnalysis();
    
    if(analysis)
    {
        AnalysisString = "Hydrogen";
        XCellLimit = VolumeX;
        YCellLimit = VolumeY;
    }
    else
    {
        AnalysisString = "Gadolinium";
        XCellLimit = 1;
        YCellLimit = 1;
    }
    
    Combine = Data->GetCombineMode();
    DataSet = Data->GetDataSet();
    
    IsotopeMatrix = Data->GetIsotopeMatrix();
    ReactorPowerMatrix = Data->GetReactorPowerMatrix();
    Sin22t12Matrix = Data->GetSin22t12Matrix();
    
    //Randomize parameters: 1 to active, 0 to desactivate
    RandomSin22t12 = 0;
    RandomIsotopeFraction = 0;
    RandomReactorPower = 0;
    
    Nweeks = Data->GetWeeks();

    TotalBins = MatrixBins;
    
    n_evis_bins = Data->GetVisibleBins();
    
    for (Int_t i = 0; i <= n_evis_bins; i++)
    {
        evis_bins[i] = Data->GetVisibleBinningArray(i);
    }
    
    NADs = Data->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    if(NADs == 8)// ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
#ifndef Produce_Antineutrino_Spectrum_For_FirstTime
    if(IsotopeMatrix||ReactorPowerMatrix)
    {
#endif
        
#ifdef UseChristineReactorModel
        Data->ReadChristineReactorSpectrum();
        
        if(IsotopeMatrix)
        {
            Data->ReadChristineCovMatrix();
        }
#else
        Data->ReadIHEPReactorSpectrum();
        
        if(IsotopeMatrix)
        {
        //Needed a covariance matrix or some way to vary the IHEP spectrum
        }
#endif

#ifndef Produce_Antineutrino_Spectrum_For_FirstTime
    }
#endif
    UseToyMCTree = Data->GetUseToyMCTree();
    
    BkgCovDirectory = Data->GetBkgCovDirectory();
    SysCovDirectory = Data->GetSysCovDirectory();
    
    GenerateFluxCorrectedHistograms(Data);
}

void Prediction :: SetSystematic()
{
    this->SetRandomIsotopeFraction(0);
    this->SetRandomReactorPower(0);
    this->SetRandomS22t12(0);
    
    if(IsotopeMatrix)
    {
        this->SetRandomIsotopeFraction(1);
    }
    if(ReactorPowerMatrix)
    {
        this->SetRandomReactorPower(1);
    }
    if(Sin22t12Matrix)
    {
        this->SetRandomS22t12(1);
        Osc->SetSin22t12(s22t12);
        OscRea->SetSin22t12(s22t12);
    }
}

//Run one time with random parameters = 0 to generate the Expected Oscillation file.
void Prediction :: MakePrediction(Double_t sin22t13, Double_t dm2_ee, bool mode, Int_t week,bool ToyMC,bool CovMatrix)
{
    if(!mode)//Nominal
    {
        RandomString = "Nominal";
    }
    else
    {
        RandomString = "Random";
    }
    
    std::cout <<  "********************************************************************************************************" << std::endl;
    std::cout << "\t Making prediction" << std::endl;
    MaxNearLoadOscModel= ADsEH1+ADsEH2;
    MaxFarLoadOscModel= ADsEH3;
    
    if(Combine == 0)
    {
        MaxNearCombine = ADsEH1+ADsEH2;
        MaxFarCombine = ADsEH3;
    }
    else if (Combine == 1)
    {
        MaxNearCombine=1;
        MaxFarCombine = 1;
    }
    else if(Combine == 2)
    {
        MaxNearCombine=2;
        MaxFarCombine = 1;
    }
    
    std::cout << "\t \t \t  Local MAX FAR OSC MODEL : " << MaxFarCombine << std::endl;
    std::cout << "\t \t \t  Local MAX NEAR OSC MODEL : " << MaxNearCombine << std::endl;
    
    s22t13=sin22t13;
    Dm2_ee=dm2_ee;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                 Vary systematics
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    OscRea= new OscillationReactor(Data);
    
    Osc= new Oscillation(Data);

    if(mode)
    {
        SetSystematic();
    }
    
    TFile CheckF("./RootOutputs/Reactor/NominalOutputs/ReactorSpectrum.root");
    TFile CheckF1(("./RootOutputs/"+AnalysisString+"/NominalOutputs/AntineutrinoSpectrum.root").c_str());
    
#ifndef Produce_Antineutrino_Spectrum_For_FirstTime
    if((((IsotopeMatrix||ReactorPowerMatrix)&&mode)||((CheckF.IsZombie()||CheckF1.IsZombie())&&(!mode))))// No need to recalculate the spectrum if it's not varied     //This has to be run once if the nominal files have not been produced beforehand
    {
#endif
        ReactorSpectrumMultiple* Reactor = new ReactorSpectrumMultiple(Data);
        
        Reactor->MultipleReactorSpectrumMain(mode);
        
        delete Reactor;
        
        AntineutrinoSpectrum* Antineutrino = new AntineutrinoSpectrum(Data);
        Antineutrino->AntineutrinoSpectrumMain(mode);
        
        delete Antineutrino;
#ifndef Produce_Antineutrino_Spectrum_For_FirstTime
    }
#endif
    //    AntineutrinoSpectrum* Antineutrino = new AntineutrinoSpectrum(Data);//This has to be run once for the nominal case
    //
    //    Antineutrino->AntineutrinoSpectrumMain();
    //
    //    delete Antineutrino;
    
    //If we use external predictions we skip the nGd Toy MC prediction generation:
    
    Osc->SetOscillationParameters(s22t13,Dm2_ee);
    
    OscRea->SetOscillationParameters(s22t13,Dm2_ee);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                    Calculate Predictions
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout << "\t Calculating Prediction" << std::endl;
    
    if(ToyMC)//ToyMC uses reactor prediction
    {
        OscRea->OscillationFromReactorData(week,mode,CovMatrix);
        
        for (Int_t AD=0; AD<NADs; AD++)
        {
            GetReactorOscillationPrediction(AD,week);
        }
    }
    
    //Possible flaw in the relative oscillation analysis as used right now. Efficiencies used in the so called superhistograms are nominal efficiencies, when variations are applied we get a residual error due to this difference.
    

    Osc->SetDelayedEfficiency(OscRea->GetDelayedEfficiency());
    Osc->SetEfficiency(OscRea->GetEfficiency());
    
    for(Int_t AD =0;AD<NADs;AD++)
    {
        for(Int_t reactor =0;reactor<NReactors;reactor++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    Osc->SetFluxHistograms(FluxHisto[reactor][AD][week][idx][idy],reactor,AD,week,idx,idy);
                }
            }
        }
    }
    
    Osc->OscillationFromNearHallData(week,ToyMC,mode);//Calculates far hall predictions from data or toy mc near hall predictions
    
    for (Int_t far=0; far<MaxFarLoadOscModel; far++)
    {
        for (Int_t near = 0; near < MaxNearLoadOscModel; near++)
        {
            PredictionVisH[far][near] = new TH1D(Form("Far AD%d from AD%d sin_%f DM_%f",far+1,near+1,sin22t13,Dm2_ee),Form("Far AD%d from AD%d sin_%f DM_%f",far+1,near+1,sin22t13,Dm2_ee),n_evis_bins,evis_bins);
            
            GetOscillationPrediction(far,near);
        }
    }

    delete Osc;
    
    std::cout << "\t \t \t IN COMBINED PREDICTION:" << std::endl;
    std::cout << "\t \t \t MaxFarCombine" << MaxFarCombine<< std::endl;
    std::cout << "\t \t \t MaxNearCombine" << MaxNearCombine << std::endl;
    
    for (Int_t far=0; far<MaxFarLoadOscModel; far++)
    {
        for (Int_t near = 0; near < MaxNearLoadOscModel; near++)
        {
            CombinedPredictionVisH[far][near]=(TH1D*)PredictionVisH[far][near]->Clone();
            
            if(Combine!=0)
            {
                CombinedPredictionVisH[far][near]->Reset();
            }
        }
    }
    
    if(ToyMC)//ToyMC uses reactor prediction
    {
        for (Int_t far=ADsEH1+ADsEH2; far<NADs; far++)
        {
            CombinedReactorPredictionVisH[far-(ADsEH1+ADsEH2)] = (TH1D*)ReactorPredictionVisH[far]->Clone();
            
            CombinedReactorPredictionVisH[far-(ADsEH1+ADsEH2)]->SetTitle(Form("Combined Reactor Prediction AD%d", far+1));
            
            if(Combine!=0)
            {
                CombinedReactorPredictionVisH[far-(ADsEH1+ADsEH2)]->Reset();

                CombinedReactorPredictionVisH[0]->Add(ReactorPredictionVisH[far]);
            }
        }
    }
    //Statistical fluctuation has to be done in the far AD directly
    if (StatisticalFluctuation&&mode)
    {
        ApplyStatisticalFluctuation(CombinedReactorPredictionVisH[0]);
    }
    
    //Combine matrices in 9x9, 2x2 or 1x1 prediction
    if (Combine == 1)
    {
        for (Int_t far = 0; far < MaxFarLoadOscModel; far++)
        {
            for (Int_t near = 0; near < MaxNearLoadOscModel; near++)
            {
                //1x1
                
                CombinedPredictionVisH[0][0]->Add(PredictionVisH[far][near]);
            }
        }
        
        CombinedPredictionVisH[0][0]->Scale(1./(ADsEH1+ADsEH2));

    }
    else if(Combine == 2)
    {
        for (Int_t far=0; far<MaxFarLoadOscModel; far++)
        {
            for (Int_t near = 0; near < MaxNearLoadOscModel; near++)//If I decide to change the limits ADsEH1+ADsEH2 for MaxNearCombine etc I need to change the delete too!! Don't forget.
            {
                //2x2
                
                if(near<ADsEH1)
                {
                    CombinedPredictionVisH[0][0]->Add(PredictionVisH[far][near]);
                }
                if(near>=ADsEH1)
                {
                    CombinedPredictionVisH[0][1]->Add(PredictionVisH[far][near]);
                }
            }
        }
        CombinedPredictionVisH[0][0]->Scale(1./(ADsEH1));
        CombinedPredictionVisH[0][1]->Scale(1./(ADsEH2));

        //Statistical fluctuation has to be done in the near AD taken to the far

//        if (StatisticalFluctuation&&Mode)
//        {
//            std::cout << "STATISTICAL FLUCTUATION IN PREDICTION" << std::endl;
//            
//            ApplyStatisticalFluctuation(CombinedPredictionVisH[0][0]);
//            ApplyStatisticalFluctuation(CombinedPredictionVisH[0][1]);
//        }
    }
    
    Char_t FileName[200];
    
    std::cout << "\t Combination done" << std::endl;
    
    //test
    
    //    TFile* f2 = new TFile(Form("./ToyMCTrees/ToyMCTreeCombined%d.root",Combine));
    //    T = (TTree*)f2->Get("TVar");
    //
    //    T->SetBranchAddress("NominalHistoDayaBay",&TNominalHisto);
    //    T->SetBranchAddress("NominalHistoLingAo",&TNominalHisto1);
    //
    //    T->GetEntry(((sin22t13*100/0.2)*101+((dm2_ee-0.0015021)*100/(0.002))));//Calculate entry number for choosen sin22t13 and dm2_ee. The tree contains 101*101 grid points, with MaxFarCombine*MaxNearCombine predictions in each point.
    //    TH1D* CopyTHisto =(TH1D*)TNominalHisto->Clone();
    //    TH1D* CopyTHisto1 =(TH1D*)TNominalHisto1->Clone();
    //    TNominalHisto->Add(CombinedPredictionVisH[0][0],-1);
    //    TNominalHisto1->Add(CombinedPredictionVisH[0][1],-1);
    //
    //    std::cout << "THE DIFFERENCE SHOULD BE 0 IN THE NOMINAL CASE: " << TNominalHisto->Integral() << std::endl;
    //    std::cout << "Index: " << (sin22t13*100/0.2)*101+((dm2_ee-0.0015021)*100/(0.002)) << "sin: " << sin22t13 << "dm: " << dm2_ee << std::endl;
    //    TCanvas* treeC = new TCanvas("c","c");
    //
    //    treeC->Divide(3,2);
    //    treeC->cd(1);
    //    TNominalHisto->Draw();
    //    treeC->cd(2);
    //    CombinedPredictionVisH[0][0]->Draw();
    //    treeC->cd(3);
    //    CopyTHisto->Draw();
    //    treeC->cd(4);
    //    TNominalHisto1->Draw();
    //    treeC->cd(5);
    //    CombinedPredictionVisH[0][1]->Draw();
    //    treeC->cd(6);
    //    CopyTHisto1->Draw();
    //    treeC->Print("./Images/Test/TreeVsPrediction.eps");
    //
    //    delete treeC;
    //    delete f2;
    //    delete CopyTHisto;
    
    
    if(!mode)
    {
        sprintf(FileName,("./RootOutputs/"+AnalysisString+Form("/Spectra/Combine%d/NominalPredictedSpectrum.root",Combine)).c_str());
    }
    else
    {
        sprintf(FileName,("./RootOutputs/"+AnalysisString+Form("/Spectra/Combine%d/PredictedSpectrum.root",Combine)).c_str());
    }
    
    if(week==0)
    {
        sprintf(optionC,"recreate");
    }
    else
    {
        sprintf(optionC,"update");
    }
    
    TFile* PredictionF1 = new TFile(FileName,optionC);
    
    for (Int_t near = 0; near<MaxNearCombine; near++)
    {
        for (Int_t far =0; far<MaxFarCombine; far++)
        {
            CombinedPredictionVisH[far][near]->Write(Form("Combined Prediction AD%d from AD%d",far+1,near+1));
        }
    }
    for (Int_t near = 0; near<MaxNearLoadOscModel; near++)
    {
        for (Int_t far =0; far<MaxFarLoadOscModel; far++)
        {
            PredictionVisH[far][near]->Write(Form("Vis Prediction AD%d from AD%d",far+1,near+1));
        }
    }
    if(ToyMC)//ToyMC uses reactor prediction
    {
        for (Int_t far =0; far<MaxFarCombine; far++)
        {
            CombinedReactorPredictionVisH[far]->Write(Form("Combined Reactor Prediction AD%d",far+1));
        }
    }

    delete PredictionF1;

    std::cout << "\t Saved done" << std::endl;
#ifdef PrintEps
    if((firstNominalPrediction==0||firstRandomPrediction==0))//To print only once
    {
        if(!mode){firstNominalPrediction=1;}
        else{firstRandomPrediction=1;}
        TCanvas* CombC = new TCanvas("CombC","CombC",MaxNearCombine*400,MaxFarCombine*400);
        CombC->Divide(MaxNearCombine,MaxFarCombine);
        
        for (Int_t near = 0; near<MaxNearCombine; near++)
        {
            for (Int_t far =0; far<MaxFarCombine; far++)
            {
                CombC->cd(MaxFarCombine*near+far+1);
                CombinedPredictionVisH[far][near]->SetStats(1);
                CombinedPredictionVisH[far][near]->Draw();
                CombC->Modified();
            }
        }
        CombC->Update();
        CombC->Print(("./Images/"+AnalysisString+"/FitterInputs/"+RandomString+"CombinePredictions.eps").c_str());
        
        delete CombC;
        
        TCanvas* PredictionC = new TCanvas("PredictionC","Predictions",MaxNearLoadOscModel*400,MaxFarLoadOscModel*400);
        
        PredictionC->Divide(MaxNearLoadOscModel,MaxFarLoadOscModel);
        
        for (Int_t far=0; far<MaxFarLoadOscModel; far++)
        {
            for (Int_t near = 0; near < MaxNearLoadOscModel; near++)
            {
                PredictionC->cd(MaxNearLoadOscModel*far+near+1);
                //                PredictionVisH[far][near]->SetStats(0);
                PredictionVisH[far][near]->Draw("HIST");
            }
        }
        PredictionC->Print(("./Images/"+AnalysisString+"/FitterInputs/"+RandomString+"Predictions.eps").c_str(),".eps");
        delete PredictionC;
        
        TCanvas* AllPredictionC = new TCanvas("AllPredictionC","AllPredictionC",MaxFarLoadOscModel*400,400);
        
        TLegend *legend2=new TLegend(0.6,0.15,0.88,0.35);
        legend2->SetTextFont(72);
        legend2->SetTextSize(0.02);//Maybe I can remove the border of the legend to fusion it within the image
        legend2->SetFillColor(0);
        
        TLegend *legend3=new TLegend(0.6,0.15,0.88,0.35);
        legend3->SetTextFont(72);
        legend3->SetTextSize(0.02);//Maybe I can remove the border of the legend to fusion it within the image
        legend3->SetFillColor(0);
        
        TLegend *legend4=new TLegend(0.6,0.15,0.88,0.35);
        legend4->SetTextFont(72);
        legend4->SetTextSize(0.02);//Maybe I can remove the border of the legend to fusion it within the image
        legend4->SetFillColor(0);
        
        AllPredictionC->Divide(MaxFarLoadOscModel,1);
        
        for (Int_t far=0; far<MaxFarLoadOscModel; far++)
        {
            AllPredictionC->cd(far+1);
            
            for (Int_t near = 0; near < MaxNearLoadOscModel; near++)
            {
                PredictionVisH[far][near]->SetStats(1);
                PredictionVisH[far][near]->SetTitle(Form("Far AD%d",far+1));
                PredictionVisH[far][near]->SetLineColor(near+1);
                PredictionVisH[far][near]->SetLineWidth(1);
                
                PredictionVisH[far][near]->Draw("same");
                if(far==1)
                {
                    legend2->AddEntry(PredictionVisH[far][near],Form("Far AD%d from Near AD%d",far+1,near+1),"l");
                    legend2->Draw("same");
                }
                else if(far==2)
                {
                    legend3->AddEntry(PredictionVisH[far][near],Form("Far AD%d from Near AD%d",far+1,near+1),"l");
                    legend3->Draw("same");
                }
                else
                {
                    legend4->AddEntry(PredictionVisH[far][near],Form("Far AD%d from Near AD%d",far+1,near+1),"l");
                    legend4->Draw("same");
                }
            }
            
            AllPredictionC->Update();
        }
        AllPredictionC->Print(("./Images/"+AnalysisString+"/FitterInputs/"+RandomString+"AllPredictions.eps").c_str(),".eps");
        delete AllPredictionC;
        if(ToyMC)//ToyMC uses reactor prediction
        {
        TCanvas* AllReactorPredC = new TCanvas("AllReactorPredC","AllReactorPredC",MaxFarLoadOscModel*400,400);
        
        AllReactorPredC->Divide(ADsEH1+ADsEH2,1);
        
        for (Int_t far=ADsEH1+ADsEH2; far<NADs; far++)
        {
            AllReactorPredC->cd(far-ADsEH1-ADsEH2+1);

                ReactorPredictionVisH[far]->SetStats(1);
                ReactorPredictionVisH[far]->SetTitle(Form("Far AD%d",far+1));
                ReactorPredictionVisH[far]->SetLineColor(far-ADsEH1-ADsEH2+1);
                ReactorPredictionVisH[far]->SetLineWidth(1);
                ReactorPredictionVisH[far]->Draw("HIST");
        }
        AllReactorPredC->Print(("./Images/"+AnalysisString+"/FitterInputs/"+RandomString+"AllReactorPredictions.eps").c_str(),".eps");

        delete AllReactorPredC;
        }
    }
#endif
    
    for (Int_t far = 0; far<MaxFarLoadOscModel; far++)
    {
        for (Int_t near = 0; near<MaxNearLoadOscModel; near++)
        {
            delete PredictionVisH[far][near];
        }
    }
    
    if(ToyMC)
    {
        OscRea->FreeMemory();
    }
    delete OscRea;
}

void Prediction :: GetOscillationPrediction(Int_t far, Int_t near)
{
    TH1D* CopyOriginalPredictionH;
    
    for(Int_t j = 0; j < n_evis_bins; j++)
    {
        CopyOriginalPredictionH = Osc->GetOscillatedADSpectrum(far,near,j);
        
        //Integrate back to visible energy:
        
        PredictionVisH[far][near]->SetBinContent(j+1, CopyOriginalPredictionH->Integral());
    }
}

void Prediction :: GetReactorOscillationPrediction(Int_t AD, Int_t week)
{
    TH1D* SumCopyOriginalPredictionH;

    for(Int_t idx=0; idx<XCellLimit; idx++)
    {
        for(Int_t idy=0; idy<YCellLimit; idy++)
        {
            TH1D* CopyOriginalPredictionH;

            CopyOriginalPredictionH = OscRea->GetReactorOscillatedADSpectrum(AD,week,idx,idy);    //Given in Visible Energy already;
            
            if(idx==0&&idy==0)
            {
                SumCopyOriginalPredictionH=(TH1D*)CopyOriginalPredictionH->Clone();
            }
            else
            {
                SumCopyOriginalPredictionH->Add(CopyOriginalPredictionH);//Add all cells
            }
        }
    }
    ReactorPredictionVisH[AD]=(TH1D*)SumCopyOriginalPredictionH->Clone();
    
    delete SumCopyOriginalPredictionH;
}

void Prediction :: SetRandomS22t12(bool randomSin22t12)
{
    RandomSin22t12 = randomSin22t12;
    
    if(RandomSin22t12)
    {
        rand->SetSeed(0);
        RandomSin22t12Error = s22t12Error * rand->Gaus(0,1);
        s22t12 = ((1 + RandomSin22t12Error) * s22t12Nominal);
        std::cout << "\t \t \t \t Random Sin22theta12 is " << s22t12 << "\n";
        std::cout << "\t \t \t \t Random Sin22theta12Error is " <<  RandomSin22t12Error << "\n";
    }
}

void Prediction :: SetRandomReactorPower(bool random_reactor_power)
{
    RandomReactorPower=random_reactor_power;
}

void Prediction :: SetRandomIsotopeFraction(bool random_isotope_fraction)
{
    RandomIsotopeFraction=random_isotope_fraction;
}

TH1D* Prediction :: GetPrediction(Int_t far,Int_t near)
{
    return CombinedPredictionVisH[far][near];
}

TH1D* Prediction :: GetReactorPrediction(Int_t AD)
{
    return ReactorPredictionVisH[AD];
}

void Prediction :: DeletePrediction(Int_t far,Int_t near)
{
    delete CombinedPredictionVisH[far][near];
}

Double_t Prediction :: CalculateChi2(Double_t sen22t13,Double_t dm2_ee, Int_t week, bool ToyMC)
{
    if (Combine == 1)
    {
        MaxNearCombine=1;
        MaxFarCombine=1;
        MaxBins=n_evis_bins;
    }
    else if(Combine == 2)
    {
        MaxNearCombine=2;
        MaxFarCombine=1;
        MaxBins=2*n_evis_bins;
    }
    else
    {
        MaxNearCombine = ADsEH1+ADsEH2;
        MaxFarCombine = ADsEH3;
        MaxBins=9*n_evis_bins;
    }
    
    std::cout << "\t Calculating χ2" << std::endl;
    
    Double_t partial_chi2;
    
    if(!UseToyMCTree)
    {
        Prediction* ScanSinPred = new Prediction(Data);
        ScanSinPred->MakePrediction(sen22t13,dm2_ee,0,week,1,0);//Vary s22t13 in each step, nominal TOY MC prediction with NO VARIATIONS
        
        for (Int_t near=0; near<MaxNearCombine; near++)
        {
            for (Int_t far=0; far<MaxFarCombine; far++)
            {
                TH1D* CopyPredictionH;
                CopyPredictionH = ScanSinPred->GetPrediction(far,near);
                PredictionH[far][near]=(TH1D*)CopyPredictionH->Clone();
            }
        }
        
        for (Int_t near=0; near<MaxNearLoadOscModel; near++)
        {
            for (Int_t far=0; far<MaxFarLoadOscModel; far++)
            {
                ScanSinPred->DeletePrediction(far,near);
            }
        }
        delete ScanSinPred;
        
        if(ToyMC)
        {
            Prediction* ScanSinDataPred = new Prediction(Data);
            ScanSinDataPred->MakePrediction(sen22t13,dm2_ee,0,week,1,0);//Vary s22t13 in each step, nominal TOY MC prediction with NO VARIATIONS
            for (Int_t near=0; near<MaxNearCombine; near++)
            {
                for (Int_t far=0; far<MaxFarCombine; far++)
                {
                    TH1D* CopyPredictionH;
                    CopyPredictionH = ScanSinDataPred->GetPrediction(far,near);
                    PredictionDataH[far][near]=(TH1D*)CopyPredictionH->Clone();
                }
            }
            
            for (Int_t near=0; near<MaxNearLoadOscModel; near++)
            {
                for (Int_t far=0; far<MaxFarLoadOscModel; far++)
                {
                    ScanSinDataPred->DeletePrediction(far,near);
                }
            }
            
            delete ScanSinDataPred;
        }
        else
        {
            Prediction* ScanSinPred2 = new Prediction(Data);
            
            ScanSinPred2->MakePrediction(sen22t13,dm2_ee,0,week,0,0);//Vary s22t13 in each step, Data prediction with no variations
            
            for (Int_t near=0; near<MaxNearCombine; near++)
            {
                for (Int_t far=0; far<MaxFarCombine; far++)
                {
                    TH1D* CopyPredictionH;
                    CopyPredictionH = ScanSinPred2->GetPrediction(far,near);
                    PredictionDataH[far][near]=(TH1D*)CopyPredictionH->Clone();
                }
            }
            
            for (Int_t near=0; near<MaxNearLoadOscModel; near++)
            {
                for (Int_t far=0; far<MaxFarLoadOscModel; far++)
                {
                    ScanSinPred2->DeletePrediction(far,near);
                }
            }
            
            delete ScanSinPred2;
        }
        
        if(sen22t13==0&&(dm2_ee>0.0024&&dm2_ee<0.0025))
        {
            TCanvas* canvasC1 = new TCanvas("canvasC1","canvasC1");
            
            if(Combine==1)
            {
                canvasC1->Divide(1,2);
            }
            else if(Combine==2)
            {
                canvasC1->Divide(2,2);
            }
            
            canvasC1->cd(1);
            PredictionH[0][0]->Draw();
            canvasC1->cd(2);
            PredictionDataH[0][0]->Draw();
            
            if(Combine==2)
            {
                canvasC1->cd(3);
                PredictionH[0][1]->Draw();
                canvasC1->cd(4);
                PredictionDataH[0][1]->Draw();
            }
            
            canvasC1->Print(("./Images/"+AnalysisString+"/NoTreeCheck.eps").c_str());
            delete canvasC1;
        }
    }
    else
    {
        TFile* f;
        //        if(!StatisticalFluctuation)
        //        {
        f = new TFile(Form("./ToyMCTrees/ToyMCTreeCombined%d.root",Combine));
        //            f = new TFile("./ToyMCTrees/TreeOnesCombine2.root");//test
        
        T = (TTree*)f->Get("TNom");//change for TNom
        T->SetBranchAddress("sin22t13",&Tsin22t13);
        T->SetBranchAddress("dm2_ee",&Tdm2_ee);
        
        if(Combine == 1)
        {
            T->SetBranchAddress("NominalHisto",&TNominalHisto);
            T->SetBranchAddress("DataHisto",&TDataHisto);
        }
        else if(Combine == 2)
        {
            T->SetBranchAddress("NominalHistoDayaBay",&TNominalHisto);
            T->SetBranchAddress("DataHistoDayaBay",&TDataHisto);
            T->SetBranchAddress("NominalHistoLingAo",&TNominalHisto1);
            T->SetBranchAddress("DataHistoLingAo",&TDataHisto1);
        }
        else
        {
            std::cout << "cannot use fitter in combine 0 mode" << std::endl;
            exit(EXIT_FAILURE);
        }
        //        }
        //        else
        //        {
        //            f = new TFile(Form("./ToyMCTrees/VariationsToyMCTreeCombined%d.root",Combine));
        //            T = (TTree*)f->Get("TVar");
        //            T->SetBranchAddress("sin22t13",&Tsin22t13);
        //            T->SetBranchAddress("dm2_ee",&Tdm2_ee);
        //            if(Combine == 1)
        //            {
        //                T->SetBranchAddress("VariationHisto",&TNominalHisto);
        //            }
        //            else if(Combine == 2)
        //            {
        //                T->SetBranchAddress("VariationHistoDayaBay",&TNominalHisto);
        //                T->SetBranchAddress("VariationHistoLingAo",&TNominalHisto1);
        //            }
        //            else
        //            {
        //                std::cout << "cannot use fitter in combine 0 mode" << std::endl;
        //                exit(EXIT_FAILURE);
        //            }
        //        }
        
        std::cout << "Number of toys: " << T->GetEntries() << std::endl;
        
        T->GetEntry(((sen22t13-s22t13start)/SinWidth)*NSteps+((dm2_ee-dm2_eestart)/DeltaWidth));//Calculate entry number for choosen sin22t13 and dm2_ee. The tree contains 100*100 grid points, with MaxFarCombine*MaxNearCombine predictions in each point.
        std::cout << "TREE ENTRY: " << ((sen22t13-s22t13start)/SinWidth)*NSteps+((dm2_ee-dm2_eestart)/DeltaWidth) << " " << (dm2_ee-dm2_eestart)/DeltaWidth << " " << (sen22t13-s22t13start)/SinWidth << std::endl;
        
        std::cout << "Sin inside tree " << Tsin22t13 << " DM inside tree " << Tdm2_ee << std::endl;
        
        if(std::abs(Tsin22t13-sen22t13)>0.00000001||std::abs(Tdm2_ee-dm2_ee)>0.00000001)
        {
            std::cout << " TREE SIN / DELTAM IS DIFFERENT FROM EXPECTED VALUE" << std::endl << Tsin22t13 << " " << sen22t13 << std::endl << Tdm2_ee << " " << dm2_ee << std::endl;
            
//            exit(EXIT_FAILURE);
        }
        
        if(sen22t13==0&&(dm2_ee>0.0024&&dm2_ee<0.0025))
        {
            TCanvas* canvasC = new TCanvas("canvasC","canvasC");
            
            if(Combine==1)
            {
                canvasC->Divide(1,2);
            }
            else if(Combine==2)
            {
                canvasC->Divide(2,2);
            }
            
            canvasC->cd(1);
            TNominalHisto->Draw();
            canvasC->cd(2);
            TDataHisto->Draw();
            
            if(Combine==2)
            {
                canvasC->cd(3);
                TNominalHisto1->Draw();
                canvasC->cd(4);
                TDataHisto1->Draw();
            }
            
            canvasC->Print(("./Images/"+AnalysisString+"/TreeCheck.eps").c_str());
            delete canvasC;
        }
        
        for (Int_t near=0; near<MaxNearCombine; near++)
        {
            for (Int_t far=0; far<MaxFarCombine; far++)
            {
                if(near<1)//Already combined in DB and LO predictions
                {
                    PredictionH[far][near]=(TH1D*)TNominalHisto->Clone();
                }
                else
                {
                    PredictionH[far][near]=(TH1D*)TNominalHisto1->Clone();
                }
                if(ToyMC)//Here the prediction and the data are the same histogram
                {
                    PredictionDataH[far][near]=(TH1D*)PredictionH[far][near]->Clone();
                }
                else//here the real data is used.
                {
                    if(near<1)
                    {
                        PredictionDataH[far][near]=(TH1D*)TDataHisto->Clone();
                    }
                    else
                    {
                        PredictionDataH[far][near]=(TH1D*)TDataHisto1->Clone();
                    }
                }
            }
        }
        delete f;
    }
    
    //Prediction difference test:
    
    #ifdef PrintEps
      /*  TH1D* DifferenceH2;
        TCanvas* DifferenceC2 = new TCanvas("diffenreceC2","differenceC2");
        DifferenceC2->Divide(MaxNearCombine,1);
        for (Int_t near=0; near<MaxNearCombine; near++)
        {
            for (Int_t far=0; far<MaxFarCombine; far++)
            {
                DifferenceH2 = (TH1D*)PredictionH[far][near]->Clone();
                DifferenceH2->Add(PredictionDataH[far][near],-1);
                DifferenceC2->cd(far+near*MaxFarCombine+1);
                DifferenceH2->Draw();
            }
        }
        DifferenceC2->Print(("./Images/"+AnalysisString+"/NominalVsDataDifference.eps").c_str());
        delete DifferenceC2;
        delete DifferenceH2;*/ //Just to check that it was working.

        if(Combine==2)
        {
            TH1D* DifferenceH;
            
            DifferenceH = (TH1D*)PredictionH[0][0]->Clone();
            DifferenceH->Add(PredictionH[0][1],-1);
            DifferenceH->Divide(PredictionH[0][0]);
            
            TCanvas* DifferenceC = new TCanvas("diffenreceC","differenceC");
            DifferenceH->Draw();
            DifferenceC->Print(("./Images/"+AnalysisString+"/PredictionDifference.eps").c_str());
            delete DifferenceC;
            delete DifferenceH;
        }
    #endif
    
    if(!Rate)
    {
        partial_chi2 = this->ChiSquare(week);
    }
    else
    {
        partial_chi2 = this->RateChiSquare(week);
        std::cout << "Rate analysis is not properly done yet" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    std::cout <<  " "  << std::endl;
    std::cout << "\t Finished calculating χ2 for sin22t13: " << sen22t13 << " and Δm2ee: " << dm2_ee << std::endl;
    std::cout <<  " "  << std::endl;
    
    //    for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)
    //    {
    //        delete NearBackgroundSpectrumH[near];
    //    }
    //    for (Int_t far =0; far<ADsEH3; far++)
    //    {
    //        delete FarBackgroundSpectrumH[far];
    //    }
    for (Int_t near = 0; near<MaxNearCombine; near++)
    {
        //        delete CombinedNearBackgroundSpectrumH[near];
        
        for (Int_t far =0; far<MaxFarCombine; far++)
        {
            //            delete PredictionH[far][near];
            delete PredictionDataH[far][near];
        }
    }
    //
    //    for (Int_t far =0; far<MaxFarCombine; far++)
    //    {
    //        delete CombinedFarBackgroundSpectrumH[far];
    //    }
    //
    //    delete BackgroundsCovarianceMatrixH;
    //    delete SystematicCovarianceMatrixH;
    //    delete StatisticalCovarianceMatrixH;
    //    delete TotalCovarianceMatrixH;
    
    return partial_chi2;
}

void Prediction :: GenerateInverseMatrix(Double_t sen22t13,Double_t dm2_ee,Int_t week,Int_t step,bool ToyMC,Int_t GenSteps)
{
    if (step==0)
    {
        sprintf(optionS,"recreate");
    }
    else
    {
        sprintf(optionS,"update");
    }
    
    if (Combine == 1)
    {
        MaxNearCombine=1;
        MaxFarCombine=1;
        MaxBins=n_evis_bins;
    }
    else if(Combine == 2)
    {
        MaxNearCombine=2;
        MaxFarCombine=1;
        MaxBins=2*n_evis_bins;
    }
    else
    {
        MaxNearCombine = ADsEH1+ADsEH2;
        MaxFarCombine = ADsEH3;
        MaxBins=9*n_evis_bins;
    }
    
    if(!UseToyMCTree)
    {
        Prediction* ScanSinPred = new Prediction(Data);
        ScanSinPred->MakePrediction(sen22t13,dm2_ee,0,week,1,0);//Vary s22t13 in each step, nominal TOY MC prediction with NO VARIATIONS
        
        for (Int_t near=0; near<MaxNearCombine; near++)
        {
            for (Int_t far=0; far<MaxFarCombine; far++)
            {
                TH1D* CopyPredictionH;
                CopyPredictionH = ScanSinPred->GetPrediction(far,near);
                PredictionH[far][near]=(TH1D*)CopyPredictionH->Clone();
                if(ToyMC)
                {
                    PredictionDataH[far][near]=(TH1D*)PredictionH[far][near]->Clone();
                }
            }
        }
        
        for (Int_t near=0; near<MaxNearLoadOscModel; near++)
        {
            for (Int_t far=0; far<MaxFarLoadOscModel; far++)
            {
                ScanSinPred->DeletePrediction(far,near);
            }
        }
        
        delete ScanSinPred;
        
        if(!ToyMC)
        {
            Prediction* ScanSinPred2 = new Prediction(Data);
            
            ScanSinPred2->MakePrediction(sen22t13,dm2_ee,0,week,0,0);//Vary s22t13 in each step, Data prediction with no variations
            
            for (Int_t near=0; near<MaxNearCombine; near++)
            {
                for (Int_t far=0; far<MaxFarCombine; far++)
                {
                    TH1D* CopyPredictionH;
                    CopyPredictionH = ScanSinPred2->GetPrediction(far,near);
                    PredictionDataH[far][near]=(TH1D*)CopyPredictionH->Clone();
                }
            }
            
            for (Int_t near=0; near<MaxNearLoadOscModel; near++)
            {
                for (Int_t far=0; far<MaxFarLoadOscModel; far++)
                {
                    ScanSinPred2->DeletePrediction(far,near);
                }
            }
            
            delete ScanSinPred2;
        }
    }
    else
    {
        TFile *f;
        
//        if(!StatisticalFluctuation)
//        {
            f = new TFile(Form("./ToyMCTrees/ToyMCTreeCombined%d.root",Combine));
            T = (TTree*)f->Get("TNom");
        
            T->SetBranchAddress("sin22t13",&Tsin22t13);
            T->SetBranchAddress("dm2_ee",&Tdm2_ee);
        
            if(Combine == 1)
            {
                T->SetBranchAddress("NominalHisto",&TNominalHisto);
                T->SetBranchAddress("DataHisto",&TDataHisto);
            }
            else if(Combine == 2)
            {
                T->SetBranchAddress("NominalHistoDayaBay",&TNominalHisto);
                T->SetBranchAddress("DataHistoDayaBay",&TDataHisto);
                T->SetBranchAddress("NominalHistoLingAo",&TNominalHisto1);
                T->SetBranchAddress("DataHistoLingAo",&TDataHisto1);
            }
            else
            {
                std::cout << "cannot use fitter in combine 0 mode" << std::endl;
                exit(EXIT_FAILURE);
            }

//        }

        //            f = new TFile("./ToyMCTrees/TreeOnesCombine2.root");//To test flat data
        
        // Commented out because it doesn't make sense to use varied data to produce the inverse matrix, it only makes sense to use it for fake data, not for predictions.
        
        //        }
        //        else
        //        {
        //            f = new TFile(Form("./ToyMCTrees/FakeExperimentsCombined_%d.root",Combine));
        //            T = (TTree*)f->Get("TFake");
        //
        //            if(Combine == 1)
        //            {
        //                T->SetBranchAddress("VariationHisto",&TNominalHisto);
        //            }
        //            else if(Combine == 2)
        //            {
        //                T->SetBranchAddress("VariationHistoDayaBay",&TNominalHisto);
        //                T->SetBranchAddress("VariationHistoLingAo",&TNominalHisto1);
        //            }
        //            else
        //            {
        //                std::cout << "cannot use fitter in combine 0 mode" << std::endl;
        //                exit(EXIT_FAILURE);
        //            }
        //        }
        //
        T->GetEntry(((sen22t13-s22t13start)/SinWidth)*GenSteps+((dm2_ee-dm2_eestart)/(DeltaWidth)));//Calculate entry number for choosen sin22t13 and dm2_ee. The tree contains 100*100 grid points, with MaxFarCombine*MaxNearCombine predictions in each point.
        
        if(std::abs(Tsin22t13-sen22t13)>0.00000001||std::abs(Tdm2_ee-dm2_ee)>0.00000001)
        {
            std::cout << " TREE SIN / DELTAM IS DIFFERENT FROM EXPECTED VALUE" << std::endl << Tsin22t13 << sen22t13 << std::endl << Tdm2_ee << dm2_ee << std::endl;
            
//            exit(EXIT_FAILURE);
        }
        
        std::cout << "Inverse Matrix TREE ENTRY: " << ((sen22t13-s22t13start)/SinWidth)*GenSteps+(dm2_ee-dm2_eestart)/(DeltaWidth) << (sen22t13/0.2)*GenSteps <<(dm2_ee-dm2_eestart)/(DeltaWidth) << std::endl;
        
        for (Int_t near=0; near<MaxNearCombine; near++)
        {
            for (Int_t far=0; far<MaxFarCombine; far++)
            {
                if(near<1)
                {
                    PredictionH[far][near]=(TH1D*)TNominalHisto->Clone();
                }
                else
                {
                    PredictionH[far][near]=(TH1D*)TNominalHisto1->Clone();
                }
                
                if(ToyMC)
                {
                    PredictionDataH[far][near]=(TH1D*)PredictionH[far][near]->Clone();
                }
                else
                {
                    if(near<1)
                    {
                        PredictionDataH[far][near]=(TH1D*)TDataHisto->Clone();
                    }
                    else
                    {
                        PredictionDataH[far][near]=(TH1D*)TDataHisto1->Clone();
                    }
                }
            }
        }
        
        delete f;
    }
    
    LoadBackgrounds(week,0);//Always nominal background to substract from nominal prediction / statistical data
    
    GenerateStatisticalCovarianceMatrix();
    
    CombineMatrices(week);
    
    InvertMatrix(week);
    
    SaveStatisticalCovarianceMatrix(week,sen22t13,dm2_ee);
    
    if(WriteROOT)
    {
        SaveChi2Components(week,sen22t13,dm2_ee);
        SaveCovarianceMatrices(week);//Only saves one sample
    }
    
    for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)
    {
        delete NearBackgroundSpectrumH[near];
    }
    for (Int_t far =0; far<ADsEH3; far++)
    {
        delete FarBackgroundSpectrumH[far];
    }
    for (Int_t near = 0; near<MaxNearCombine; near++)
    {
        delete CombinedNearBackgroundSpectrumH[near];
        
        for (Int_t far =0; far<MaxFarCombine; far++)
        {
            delete PredictionH[far][near];
            delete PredictionDataH[far][near];
        }
    }
    
    for (Int_t far =0; far<MaxFarCombine; far++)
    {
        delete CombinedFarBackgroundSpectrumH[far];
    }
    
    delete BackgroundsCovarianceMatrixH;
    delete SystematicCovarianceMatrixH;
    delete StatisticalCovarianceMatrixH;
    delete TotalCovarianceMatrixH;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Produces the Statistical Covariance Matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Prediction :: GenerateStatisticalCovarianceMatrix()
{
    std::cout << "\t Generating StatisticalCovarianceMatrix" << std::endl;
    
    StatisticalCovarianceMatrixH = new TH2D("Statistical Covariance Matrix","Statistical Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    
    #ifdef PrintEps
        std::cout << " Printing" << std::endl;
        
        TCanvas* cstat = new TCanvas("TestStatisticalmatrix.eps","TestStatisticalmatrix.eps");
        cstat->Divide(3,1);
        cstat->cd(1);
        
        CombinedNearDataH[0]->Draw();
        
        cstat->cd(2);
        CombinedNearDataH[1]->Draw();
        
        cstat->cd(3);
        CombinedNearDataH[2]->Draw();
        
    cstat->Print(("./Images/"+AnalysisString+"/TestStatisticalmatrix.eps").c_str(), ".eps");
        
        delete cstat;
    #endif
    
    Double_t Sigma_Near[MaxFarDetectors][MaxNearDetectors][MaxNbins];
    Double_t Sigma_Far[MaxFarDetectors][MaxNearDetectors][MaxNbins];
    
    for (Int_t far=0; far<MaxFarCombine; far++)
    {
        for (Int_t near=0; near<MaxNearCombine; near++)
        {
            for (Int_t pts = 0; pts < n_evis_bins; pts++)
            {
                Sigma_Far[far][near][pts]=sqrt(PredictionDataH[far][near]->GetBinContent(pts+1)+CombinedFarBackgroundSpectrumH[far]->GetBinContent(pts+1));
                //                std::cout << "SIGMA FAR " << Sigma_Far[far][near][pts] << std::endl;
                //N OBSERVED, F PREDICTED FROM N OBSERVED, BACKGROUNDS PREDICTED (WITHOUT ANY FLUCTUATION)
                
                Sigma_Near[far][near][pts]=(PredictionDataH[far][near]->GetBinContent(pts+1)/CombinedNearDataH[near]->GetBinContent(pts+1))*sqrt(CombinedNearDataH[near]->GetBinContent(pts+1)+CombinedNearBackgroundSpectrumH[near]->GetBinContent(pts+1));
                //                std::cout << "SIGMA NEAR " << Sigma_Near[far][near][pts] << std::endl;
                
            }
        }
    }
    
    Int_t x =0;
    Int_t y =0;
    
    Int_t Ni1=0,Ni2=0,Ni3=0,Ni4=0;
    Int_t Nj1=0,Nj2=0,Nj3=0,Nj4=0;
    Int_t Fi1=0,Fi2=0,Fi3=0,Fi4=0;
    Int_t Fj1=0,Fj2=0,Fj3=0,Fj4=0;
    
    for (Int_t neari=0; neari<MaxNearCombine; neari++)
    {
        //Logic for the 2D matrix index done up to 8 ADs
        if(neari==0){Ni1=1;Ni2=0;Ni3=0;Ni4=0;}
        if(neari==1){Ni2++;}
        if(neari==2){Ni3++;}
        if(neari==3){Ni4++;}
        
        for (Int_t fari=0; fari<MaxFarCombine; fari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(Ni1!=Ni2){Fi1=fari+1;Fi2=0;Fi3=0;Fi4=0;}
            if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFarCombine; Fi2=fari+1;}
            if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFarCombine; Fi3=fari+1;}
            if(Ni3==Ni4&&Ni4==1){Fi3=MaxFarCombine; Fi4=fari+1;}
            
            for (Int_t nearj=0; nearj<MaxNearCombine; nearj++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                if(nearj==1){Nj2++;}
                if(nearj==2){Nj3++;}
                if(nearj==3){Nj4++;}
                
                for (Int_t farj=0; farj<MaxFarCombine; farj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(Nj1!=Nj2){Fj1=farj+1;Fj2=0;Fj3=0;Fj4=0;}
                    if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFarCombine; Fj2=farj+1;}
                    if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFarCombine; Fj3=farj+1;}
                    if(Nj3==Nj4&&Nj4==1){Fj3=MaxFarCombine; Fj4=farj+1;}
                    
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {//columns
                        x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                        
                        for (Int_t j = 0; j<n_evis_bins; j++)
                        {//rows
                            y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                            
                            if(i==j)
                            {
                                //Near component correlated
                                if(neari==nearj && fari!=farj)
                                {
                                    CovStat[x][y]=Sigma_Near[fari][neari][i]*Sigma_Near[farj][nearj][j];
                                }
                                //Far component correlated
                                if(fari==farj && neari!=nearj)
                                {
                                    CovStat[x][y]=Sigma_Far[fari][neari][i]*Sigma_Far[farj][nearj][j];
                                }
                                if(neari==nearj && fari==farj)
                                {
                                    //General covariance
                                    CovStat[x][y]=(Sigma_Near[fari][neari][i]*Sigma_Near[fari][neari][j])+(Sigma_Far[farj][nearj][i]*Sigma_Far[farj][nearj][j]);
                                }
                                //Uncorrelated terms
                                if(neari!=nearj && fari!=farj)
                                {
                                    CovStat[x][y]=0;
                                }
                            }
                            
                            StatisticalCovarianceMatrixH->SetBinContent(x+1,y+1,CovStat[x][y]);
                        }
                    }
                }
            }
        }
    }
    
    #ifdef PrintEps
        TCanvas* StatisticalCovarianceMatrixC = new TCanvas("StatisticalCovarianceMatrixC","Statistical Cov",400,400);
        StatisticalCovarianceMatrixH->SetStats(0);
        StatisticalCovarianceMatrixH->Draw("colz");
        
        StatisticalCovarianceMatrixC->Print(("./Images/"+AnalysisString+"/StatisticalCovarianceMatrix.eps").c_str(),".eps");
        
        delete StatisticalCovarianceMatrixC;
    #endif
    
    std::cout << "\t Finished Generating StatisticalCovarianceMatrix" << std::endl;
    std::cout <<  "\t ***********************************************************************************************" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                        Save Statistical Covariance Matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Prediction :: SaveStatisticalCovarianceMatrix(Int_t week,Double_t sen22t13,Double_t dm2_ee)
{
    std::cout <<  "\t ***********************************************************************************************" << std::endl;
    std::cout << "\t Saving StatisticalCovarianceMatrix" << optionS << std::endl;
    
    //Save statistical matrix
    TFile* SaveStatisticalCovMatrixF = TFile::Open(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/StatisticalCovarianceMatrix.root",Combine)).c_str(),optionS);
    
    StatisticalCovarianceMatrixH->Write(Form("Statistical Covariance Matrix for sin22t13 %f and Δm2ee %f period%d",sen22t13,dm2_ee,week));
    
    if(WriteROOT)
    {
        for (Int_t near = 0; near<MaxNearCombine; near++)
        {
            CombinedNearBackgroundSpectrumH[near]->Write(Form("Near Background Spectrum for sin22t13 %f and Δm2ee %f",sen22t13,dm2_ee));
            CombinedNearDataH[near]->Write(Form("Near Data Spectrum for sin22t13 %f and Δm2Dm2ee %f",sen22t13,dm2_ee));
            
            for (Int_t far =0; far<MaxFarCombine; far++)
            {
                FarDataH[far][near]->Write(Form("Far Data Spectrum for sin22t13 %f and Δm2Dm2ee %f",sen22t13,dm2_ee));
                PredictionH[far][near]->Write(Form("Far Prediction Spectrum for sin22t13 %f and Δm2Dm2ee %f",sen22t13,dm2_ee));
                PredictionDataH[far][near]->Write(Form("Far Prediction Spectrum from Data for sin22t13 %f and Δm2Dm2ee %f",sen22t13,dm2_ee));
            }
        }
        
        for (Int_t far =0; far<MaxFarCombine; far++)
        {
            CombinedFarBackgroundSpectrumH[far]->Write(Form("Far Background Spectrum for sin22t13 %f and Δm2Dm2ee %f",sen22t13,dm2_ee));
        }
    }
    SaveStatisticalCovMatrixF->Close();
    
    //Save in a txt file
    if(WriteOutput)
    {
        Char_t WriteC[50];
        sprintf(WriteC,"CovarianceMatrices/Combine%d/StatisticalCovarianceMatrix.txt",Combine);
        ofstream statf(WriteC);
        Int_t x =0;
        Int_t y =0;
        
        Int_t Ni1=0,Ni2=0,Ni3=0,Ni4=0;
        Int_t Nj1=0,Nj2=0,Nj3=0,Nj4=0;
        Int_t Fi1=0,Fi2=0,Fi3=0,Fi4=0;
        Int_t Fj1=0,Fj2=0,Fj3=0,Fj4=0;
        
        for (Int_t neari=0; neari<MaxNearCombine; neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Ni1=1;Ni2=0;Ni3=0;Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<MaxFarCombine; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Fi1=fari+1;Fi2=0;Fi3=0;Fi4=0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFarCombine;Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFarCombine;Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=MaxFarCombine;Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<MaxNearCombine; nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<MaxFarCombine; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Fj1=farj+1;Fj2=0;Fj3=0;Fj4=0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFarCombine; Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFarCombine; Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=MaxFarCombine; Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                            
                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                
                                statf << StatisticalCovarianceMatrixH->GetBinContent(x+1,y+1) << " ";
                            }
                        }
                    }
                }
            }
        }
        statf << std::endl;
        statf.close();
    }
    std::cout << "\t Finished Saving StatisticalCovarianceMatrix" << std::endl;
}

void Prediction :: SaveChi2Components(Int_t week,Double_t sen22t13, Double_t dm2_ee)
{
    TFile* SaveChi2Components1 = TFile::Open(("./ChiSquare/"+AnalysisString+Form("/Combine%d/Chi2ComponentsPeriod%d.root",Combine,week)).c_str(),optionS);
    
    for (Int_t far=0; far<MaxFarCombine; far++)
    {
        for (Int_t near=0; near<MaxNearCombine; near++)
        {
            PredictionDataH[far][near]->Write(Form("Prediction far%d, near%d for sin22t13 %f and Δm2Dm2ee %f",far,near,sen22t13,dm2_ee));
            FarDataH[far][near]->Write(Form("Far Data far%d, near%d for sin22t13 %f and Δm2Dm2ee %f",far,near,sen22t13,dm2_ee));
            CopyFarDataH[far][near]=(TH1D*)FarDataH[far][near]->Clone();
            CopyFarDataH[far][near]->Add(PredictionDataH[far][near],-1);
            CopyFarDataH[far][near]->Write(Form("Far Data - Prediction fari%d, neari%d for sin22t13 %f and Δm2Dm2ee %f",far,near,sen22t13,dm2_ee));
            delete CopyFarDataH[far][near];
        }
    }
    
    InvTotalCovarianceMatrixH->Write(Form("Inverse Matrix for sin22t13 %f and Δm2Dm2ee %f",sen22t13,dm2_ee));
    
    SaveChi2Components1->Close();
}

void Prediction :: LoadData(Int_t week,bool ToyMC,Int_t DataSteps,bool Mode)//Called from Fitter
{
    ToyMCflag = ToyMC;
    
    if(Combine == 0)
    {
        MaxNearCombine = ADsEH1+ADsEH2;
        MaxFarCombine = ADsEH3;
    }
    else if (Combine == 1)
    {
        MaxNearCombine=1;
        MaxFarCombine=1;
    }
    else if(Combine == 2)
    {
        MaxNearCombine=2;
        MaxFarCombine=1;
    }
    
    LoadBackgrounds(week,Mode);//To substract them from the data, the toy MC spectrum has background added, substract after to account for statistical fluctuations in the bkgd.

    if(ToyMC)//ToyMC
    {
        //create predictions
        Prediction* DataPred = new Prediction(Data);
        
        DataPred->MakePrediction(Sin22t13,DM2_ee,Mode,week,1,0);// Let's use the data for a fixed real sin22t13 and relative oscillation model. Set mode = 1 to apply variations if desired.
        
        //Need to scale data for the weighted expected number of events!
        
        for (Int_t near=0; near<ADsEH1+ADsEH2; near++)
        {
            NearDataH[near] = DataPred->GetReactorPrediction(near);
            
            if (StatisticalFluctuation)
            {
                ApplyStatisticalFluctuation(NearDataH[near]);
            }
        }
        
        //Substract nominal background from data after any statistical fluctuation has been applied.
        
        for (Int_t near=0; near<ADsEH1+ADsEH2; near++)
        {
            NearDataH[near]->Add(NearBackgroundSpectrumH[near],-1);
        }
        
        if(!UseToyMCTree)
        {
            for (Int_t far=0; far<MaxFarCombine; far++)
            {
                for (Int_t near=0; near<MaxNearCombine; near++)
                {
                    FarDataH[far][near] = DataPred->GetPrediction(far,near);
                }
            }
            
            for (Int_t near=0; near<MaxNearCombine; near++)
            {
                for (Int_t far=0; far<MaxFarCombine; far++)
                {

                    FarDataH[far][near]->Add(CombinedFarBackgroundSpectrumH[far]);//because the prediction is done without backgrounds, so add them here
                    
                    if (StatisticalFluctuation)
                    {
                        ApplyStatisticalFluctuation(FarDataH[far][near]);
                    }
                }
            }
        }
        
        //delete prediction
        
        for (Int_t near=0; near<MaxNearLoadOscModel; near++)
        {
            for (Int_t far=0; far<MaxFarLoadOscModel; far++)
            {
                DataPred->DeletePrediction(far,near);//To avoid leaks
            }
        }
        
        delete DataPred;
        
        if(!UseToyMCTree)//otherwise this will be overrriden by loading the predictions from the tree
        {
            for (Int_t near=0; near<MaxNearCombine; near++)
            {
                for (Int_t far=0; far<MaxFarCombine; far++)
                {
                    FarDataH[far][near]->Add(CombinedFarBackgroundSpectrumH[far],-1);
                }
            }
        }
        else
        {
            TFile* f[MaxSystematics];
          
            if(!StatisticalFluctuation)//Nominal mode, otherwise test mode
            {
                NExperiments = 1;//This will make the difference between running a big number of fake experiments for performance testing or just fitting a fake experiment for fitter fake data test.
                
                f[0] = new TFile(Form("./ToyMCTrees/ToyMCTreeCombined%d.root",Combine));
                
                T = (TTree*)f[0]->Get("TNom");
                T->SetBranchAddress("sin22t13",&Tsin22t13);
                T->SetBranchAddress("dm2_ee",&Tdm2_ee);
              
                if(Combine == 1)
                {
                    T->SetBranchAddress("NominalHisto",&TNominalHisto);
                    T->SetBranchAddress("DataHisto",&TDataHisto);
                }
                else if(Combine == 2)
                {
                    T->SetBranchAddress("NominalHistoDayaBay",&TNominalHisto);
                    T->SetBranchAddress("DataHistoDayaBay",&TDataHisto);
                    T->SetBranchAddress("NominalHistoLingAo",&TNominalHisto1);
                    T->SetBranchAddress("DataHistoLingAo",&TDataHisto1);
                }
                else
                {
                    std::cout << "cannot use fitter in combine 0 mode" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            else
            {
                for(Int_t SystematicI = 8; SystematicI < MaxSystematics; SystematicI++)
                {
                    //SystematicI = 8 to use all variations. Otherwise select accordingly
                    
                    if(Fake_Experiments)
                    {
                        std::cout << "USING 1000 FAKE EXPERIMENTS 21x21 VARIATIONS TREE" << std::endl;

                        f[SystematicI] = new TFile(Form("./ToyMCTrees/FakeExperiments%d_Combined_%d.root",SystematicI,Combine));
                        
                        T = (TTree*)f[SystematicI]->Get("TFake");
                        T->SetBranchAddress("sin22t13",&Tsin22t13);
                        T->SetBranchAddress("dm2_ee",&Tdm2_ee);
                        T->SetBranchAddress("Experiment",&TExperiment);
                        
                        if(Combine == 1)
                        {
                            T->SetBranchAddress(Form("VariationHisto"),&TNominalHisto);
                        }
                        else if(Combine == 2)
                        {
                            T->SetBranchAddress(Form("VariationHistoDayaBay"),&TNominalHisto);
                            T->SetBranchAddress(Form("VariationHistoLingAo"),&TNominalHisto1);
                        }
                        else
                        {
                            std::cout << "cannot use fitter in combine 0 mode" << std::endl;
                            exit(EXIT_FAILURE);
                        }
                    }
                    else
                    {
                        
                        std::cout << "USING 101x101 VARIATIONS TREE" << std::endl;

                        f[SystematicI] = new TFile(Form("./ToyMCTrees/Variations%d_ToyMCTreeCombined%d.root",SystematicI,Combine));
                        
                        T = (TTree*)f[SystematicI]->Get("TVar");
                        T->SetBranchAddress("sin22t13",&Tsin22t13);
                        T->SetBranchAddress("dm2_ee",&Tdm2_ee);
                        T->SetBranchAddress("Experiment",&TExperiment);
                        
                        if(Combine == 1)
                        {
                            T->SetBranchAddress(Form("VariationHisto"),&TNominalHisto);
                        }
                        else if(Combine == 2)
                        {
                            T->SetBranchAddress(Form("VariationHistoDayaBay"),&TNominalHisto);
                            T->SetBranchAddress(Form("VariationHistoLingAo"),&TNominalHisto1);
                        }
                        else
                        {
                            std::cout << "cannot use fitter in combine 0 mode" << std::endl;
                            exit(EXIT_FAILURE);
                        }
                    }
                }
            }
            
            std::cout << "Number of toys: " << T->GetEntries() << std::endl;
            
            T->GetEntry(Experiment+NExperiments*((Sin22t13-s22t13start)/SinWidth)*DataSteps+((DM2_ee-dm2_eestart)/(DeltaWidth)));//Calculate entry number for choosen sin22t13 and dm2_ee. The tree contains 100*100 grid points, with MaxFarCombine*MaxNearCombine predictions in each point.
            std::cout << "Data TREE ENTRY: " << Experiment << (Experiment+NExperiments*((Sin22t13-s22t13start)/SinWidth)*DataSteps+((DM2_ee-dm2_eestart)/(DeltaWidth))) << " " << (DM2_ee-dm2_eestart)/(DeltaWidth) << " " << (Sin22t13-s22t13start)/SinWidth << std::endl;
            
            std::cout << "   Data Sin inside tree " << Tsin22t13 << std::endl << "   Data DM inside tree " << Tdm2_ee << std::endl << std::endl;
            
            if(std::abs(Tsin22t13-Sin22t13)>0.00000001||std::abs(Tdm2_ee-DM2_ee)>0.00000001)
            {
                std::cout << " TREE SIN / DELTAM IS DIFFERENT FROM EXPECTED VALUE:" << std::endl << Tsin22t13 << " " << Sin22t13 << " " << std::endl << Tdm2_ee << " " << DM2_ee << std::endl;
                
                exit(EXIT_FAILURE);
            }
            
            for (Int_t near=0; near<MaxNearCombine; near++)
            {
                for (Int_t far=0; far<MaxFarCombine; far++)
                {
                    
                    //Dilema: we have two predictions (DB and LO) but data should be 1 unique histogram. The difference between them is less than 1% and probably this is negligible, but this implies that the total chi2 will never be 0. I am going to use an average of both of them.
                    
                    FarDataH[far][near]=(TH1D*)TNominalHisto->Clone();//No need to apply fluctuations if using tree, since they are already included
                    FarDataH[far][near]->Add(TNominalHisto1);
                    FarDataH[far][near]->Scale(1./2);
                    
                    // If we set the far data histograms to be different and exactly matching the predictions, then the chi2 should be 0.
                    if(near<1)
                    {
                        FarDataH[far][near]=(TH1D*)TNominalHisto->Clone();//No need to apply fluctuations if using tree, since they are already included

                    }
                    else
                    {
                        FarDataH[far][near]=(TH1D*)TNominalHisto1->Clone();//No need to apply fluctuations if using tree, since they are already included
                    }
                    
                    //Use an average of the predictions, this is arbitrary, the difference between predictions without fluctuations should be negligible.
                    //Before I expected that by using different FarDataH the chi2 should be minimal, but realistically the  real fardata must be the same since we don't have data based on different near hall information.
                }
            }
            if(!StatisticalFluctuation)
            {
                delete f[0];
            }
            else
            {
                for(Int_t SystematicI = 8; SystematicI < MaxSystematics; SystematicI++)
                {
                    delete f[SystematicI];
                }
            }
        }
        
        #ifdef PrintEps
            TCanvas* NearC = new TCanvas("NearC","Near Data", (ADsEH1+ADsEH2)*400,400);
            NearC->Divide(ADsEH1+ADsEH2,1);
            
            for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
            {
                NearC->cd(near+1);
                NearDataH[near]->SetStats(1);
                NearDataH[near]->Draw("HIST");
            }
            
            NearC->Update();
            
            TCanvas* FarC = new TCanvas("FarC","Far Data", MaxNearCombine*400,MaxFarCombine*400);
            FarC->Divide(MaxNearCombine,MaxFarCombine);
            
            for (Int_t near=0; near<MaxNearCombine; near++)
            {
                for(Int_t far = 0; far < MaxFarCombine; far++)
                {
                    FarC->cd(MaxFarCombine*near+far+1);
                    FarDataH[far][near]->Draw("HIST");
                }
            }
            FarC->Update();
            
            if (StatisticalFluctuation)
            {
                NearC->Print(("./Images/"+AnalysisString+"/FitterInputs/CombinedFluctuatedNearDataInFitter.eps").c_str(),".eps");
                FarC->Print(("./Images/"+AnalysisString+"/FitterInputs/CombinedFluctuatedFarDataInFitter.eps").c_str(),".eps");
                
            }
            else
            {
                NearC->Print(("./Images/"+AnalysisString+"/FitterInputs/CombinedNearToyMCDataInFitter.eps").c_str(),".eps");
                FarC->Print(("./Images/"+AnalysisString+"/FitterInputs/CombinedFarToyMCDataInFitter.eps").c_str(),".eps");
            }
            
            delete NearC;
            delete FarC;
        #endif
        
        for (Int_t far=0; far<ADsEH3; far++)
        {
            for (Int_t near=0; near<ADsEH1+ADsEH2; near++)
            {
                delete ToyFarDataH[far][near];
            }
        }
    }
    else
    {
        //  Real near data
        
        std::cout << " LOADING NEAR DATA " << std::endl;
        Char_t DataFile[100];
        Char_t FarDataSpec[100];
        Char_t NearDataSpec[100];
        
        if(analysis)
        {
            std::cout << "\t \t \t \t \t \t IBD DATA FILE NEEDED FOR HYDROGEN ANALYSIS" << std::endl;
            sprintf(DataFile,"./Inputs/HInputs/ibd_eprompt_shapes.root");
        }
        else
        {
            if(DataSet==2)
            {
                if(ADSimple) //P12E AD Simple
                {
                    sprintf(DataFile,"./Inputs/GdInputs/ibd_eprompt_shapes.root");//ADSimple
                }
                else //LBNL P14 AD Scaled
                {
                    sprintf(DataFile,"./Inputs/GdInputs/IHEP_data_lbnlbin_6AD.root");//ADScaled
                }
            }
            else if(DataSet==1)//P12E
            {
                sprintf(DataFile,"./Inputs/GdInputs/file_P12E_nGd_IBD_Acc_spectrum.root");
                //  This file has already data with accidentals substracted!
            }
            else if(DataSet==0)
            {
                std::cout << "ERROR, USING NEAR DATA WITH SIMPLE REACTOR MODEL" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
        TFile* DataF = TFile::Open(DataFile);
        
        for(Int_t far = 0; far < ADsEH3; far++)
        {
            if(Nweeks == 1)
            {
                if(DataSet==2)
                {
                    if(ADSimple)//P12E AD Simple
                    {
                        sprintf(FarDataSpec,"h_ibd_eprompt_inclusive_ad%i",far+ADsEH1+ADsEH2+1);//ADSimple
                    }
                    else//LBNL P14 AD Scaled
                    {
                        sprintf(FarDataSpec,"h_ibd_eprompt_inclusive_eh3_ad%i",far+1);
                    }
                }
                else if(DataSet==1)//P12E
                {
                    sprintf(FarDataSpec,"hist_AccSub_%d_EH3",far+1);
                }
            }
            else
            {
                if(DataSet==2)//LBNL
                {
                    sprintf(FarDataSpec,"h_ibd_eprompt_week%d_ad%d", week, far+1);
                }
            }
            
            FarDataH[far][0] = (TH1D*)gDirectory->Get(FarDataSpec);
            
            if(DataSet!=2)//need to rebin files to LBNL binning
            {
                FarDataH[far][0]=(TH1D*)FarDataH[far][0]->Rebin(n_evis_bins,Form("Rebinned Vis Far%d Data Spectrum",far),evis_bins);
            }

            FarDataH[far][0]->Add(FarBackgroundSpectrumH[far],-1);//Substract backgrounds from data
            
            for(Int_t i = 0; i<FarDataH[far][0]->GetXaxis()->GetNbins();i++)
            {
                if(FarDataH[far][0]->GetBinContent(i+1)<0)
                {
                    FarDataH[far][0]->SetBinContent(i+1,0);
                    std::cout << "Background larger than signal" << std::endl;
                    exit(EXIT_FAILURE);
                    
                }
            }
            if(Combine==2)
            {
                FarDataH[far][1]=(TH1D*)FarDataH[far][0]->Clone();
            }
        }
        
        #ifdef PrintEps
            TCanvas* FarDataC = new TCanvas("FarDataC","Far Data", ADsEH3*400,400);
            FarDataC->Divide(ADsEH3,1);
            
            for(Int_t far = 0; far < ADsEH3; far++)
            {
                FarDataC->cd(far+1);
                FarDataH[far][0]->Draw("HIST");
            }
            
            FarDataC->Update();
            
            FarDataC->Print(("./Images/"+AnalysisString+"/FitterInputs/RealFarDataInFitter.eps").c_str(),".eps");
            
            delete FarDataC;
        #endif
        
        for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            if(Nweeks == 1)
            {
                if(DataSet==2)//P12E given by LBNL
                {
                    sprintf(NearDataSpec,"h_ibd_eprompt_inclusive_ad%i",near+1);
                }
                else if(DataSet==1)//P12E given by Xiang Pan
                {
                    if(near<ADsEH1)
                    {
                        sprintf(NearDataSpec,"hist_AccSub_%d_EH1",near+1);
                    }
                    else
                    {
                        sprintf(NearDataSpec,"hist_AccSub_%d_EH2",near-ADsEH1+1);
                    }
                }
            }
            else
            {
                if(DataSet==2)//LBNL
                {
                    sprintf(NearDataSpec,"h_ibd_eprompt_week%d_ad%d", week, near+1);
                }
            }
            
            NearDataH[near] = (TH1D*)gDirectory->Get(NearDataSpec);
            
            
            //need to rebin files to LBNL binning
            
            if(DataSet!=2)//need to rebin files to LBNL binning
            {
                NearDataH[near]=(TH1D*)NearDataH[near]->Rebin(n_evis_bins,Form("Rebinned Vis Near%d Data Spectrum",near),evis_bins);
            }
            else if(DataSet==1)
            {

                NearDataH[near]->Add(NearBackgroundSpectrumH[near],-1);//Substract backgrounds from LBNL data
            }
            
            //Set 0 negative bins after subtraction.
            
            for(Int_t i = 0; i< NearDataH[near]->GetXaxis()->GetNbins();i++)
            {
                
                if( NearDataH[near]->GetBinContent(i+1)<0)
                {
                    std::cout << "Backgrounds larger than data by " <<  NearDataH[near]->GetBinContent(i+1)<< std::endl;
                    exit(EXIT_FAILURE);
                    
                    NearDataH[near]->SetBinContent(i+1,0);
                }
            }
        }
        
        DataF->Close();
    }
    
    #ifdef PrintEps
        TCanvas* NearDataC = new TCanvas("NearDataC","Near Data", (ADsEH1+ADsEH2)*400,400);
        
        NearDataC->Divide(ADsEH1+ADsEH2,1);
        
        for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            NearDataC->cd(near+1);
            NearDataH[near]->SetStats(1);
            NearDataH[near]->Draw("HIST");
        }
        
        NearDataC->Update();
        
        if(ToyMC&&StatisticalFluctuation)
        {
            NearDataC->Print(("./Images/"+AnalysisString+"/FitterInputs/FluctuatedToyMCDataInFitter.eps").c_str(),".eps");
        }
        else if(ToyMC)
        {
            NearDataC->Print(("./Images/"+AnalysisString+"/FitterInputs/ToyMCDataInFitter.eps").c_str(),".eps");
        }
        else
        {
            NearDataC->Print(("./Images/"+AnalysisString+"/FitterInputs/RealNearDataInFitter.eps").c_str(),".eps");
        }
        delete NearDataC;
    #endif
    
    for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
    {
        CombinedNearDataH[near]=(TH1D*)NearDataH[near]->Clone();
    }
    
    //Combine matrices in 9x9, 2x2 or 1x1 prediction
    if (Combine == 1)
    {
        for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            CombinedNearDataH[near]->Reset();
        }
        
        for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            //1x1
            
            CombinedNearDataH[0]->Add(NearDataH[near]);
        }
        
        //        CombinedNearDataH[0]->Scale(1./(ADsEH1+ADsEH2));
        if(!ToyMC)
        {
            for (Int_t far = 1; far < ADsEH3; far++)
            {
                FarDataH[0][0]->Add(FarDataH[far][0]);//  Far hall all together due to low statistics.
            }
        }
    }
    else if(Combine == 2)
    {
        for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            CombinedNearDataH[near]->Reset();
        }
        
        for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            //2x2
            if(near<ADsEH1)
            {
                CombinedNearDataH[0]->Add(NearDataH[near]);//DB
                //                CombinedNearDataH[0]->Scale(1./ADsEH1);//Since Predictions are scaled, I have to scale also near data
                
            }
            if(near>=ADsEH1)
            {
                CombinedNearDataH[1]->Add(NearDataH[near]);//LA
                //                CombinedNearDataH[1]->Scale(1./ADsEH2);//DB
            }
        }
        if(!ToyMC)
        {
            for (Int_t far = 1; far < ADsEH3; far++)
            {
                FarDataH[0][0]->Add(FarDataH[far][0]);//  Far hall all together due to low statistics.
                FarDataH[0][1]->Add(FarDataH[far][0]);//  Far hall all together due to low statistics.
            }
        }
    }
    for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
    {
        delete NearDataH[near];//Not used anymore;
    }
    for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)
    {
        delete NearBackgroundSpectrumH[near];
    }
    for (Int_t far =0; far<ADsEH3; far++)
    {
        delete FarBackgroundSpectrumH[far];
    }
    for (Int_t near = 0; near<MaxNearCombine; near++)
    {
        delete CombinedNearBackgroundSpectrumH[near];
    }
    
    for (Int_t far =0; far<MaxFarCombine; far++)
    {
        delete CombinedFarBackgroundSpectrumH[far];
    }
    
    #ifdef PrintEps
        TCanvas* CombinedNearC = new TCanvas("CombinedNearC","Combined Near Data", MaxNearCombine*400,400);
        
        CombinedNearC->Divide(MaxNearCombine,1);
        
        for(Int_t near = 0; near < MaxNearCombine; near++)
        {
            CombinedNearC->cd(near+1);
            
            CombinedNearDataH[near]->SetStats(1);
            CombinedNearDataH[near]->Draw("HIST");
        }
        
        CombinedNearC->Update();
        
        if(ToyMC&&StatisticalFluctuation)
        {
            CombinedNearC->Print(("./Images/"+AnalysisString+"/FitterInputs/CombinedFluctuatedNearDataInFitter.eps").c_str(),".eps");
        }
        else if(ToyMC)
        {
            CombinedNearC->Print(("./Images/"+AnalysisString+"/FitterInputs/CombinedToyMCNearDataInFitter.eps").c_str(),".eps");
        }
        else
        {
            CombinedNearC->Print(("./Images/"+AnalysisString+"/FitterInputs/CombinedRealNearDataInFitter.eps").c_str(),".eps");
        }
        delete CombinedNearC;
    
        TCanvas* CombinedFarC = new TCanvas("CombinedFarC","Combined Far Data", MaxNearCombine*400, MaxFarCombine*400);
        
        CombinedFarC->Divide(MaxNearCombine,MaxFarCombine);
        
        for(Int_t far = 0; far < MaxFarCombine; far++)
        {
            for(Int_t near = 0; near < MaxNearCombine; near++)
            {
                CombinedFarC->cd(MaxNearCombine*far+near+1);
                
                FarDataH[far][near]->Draw("HIST");
            }
        }
        
        CombinedFarC->Update();
        
        if(ToyMC&&StatisticalFluctuation)
        {
            CombinedFarC->Print(("./Images/"+AnalysisString+"/FitterInputs/CombinedFluctuatedFarDataInFitter.eps").c_str(),".eps");
        }
        else if(ToyMC)
        {
            CombinedFarC->Print(("./Images/"+AnalysisString+"/FitterInputs/CombinedToyMCFarDataInFitter.eps").c_str(),".eps");
        }
        else
        {
            CombinedFarC->Print(("./Images/"+AnalysisString+"/FitterInputs/CombinedRealFarDataInFitter.eps").c_str(),".eps");
        }
        
        delete CombinedFarC;
    #endif
}

void Prediction :: LoadBackgrounds(Int_t week,bool mode)
{
    Char_t backgroundC[100];
    
    if (mode)
    {
        RandomString = "Random";
    }
    else
    {
        RandomString = "Nominal";
    }
    
    std::cout <<  "\t ***********************************************************************************************" << std::endl;
    std::cout <<"\t \t IN LOAD" << RandomString << " BACKGROUNDS" << std::endl;
    
    std::cout << "\t \t \t MaxFarCombine" << MaxFarCombine<< std::endl;
    std::cout << "\t \t \t MaxNearCombine" << MaxNearCombine << std::endl;

    sprintf(backgroundC,("./RootOutputs/"+ AnalysisString+ "/Backgrounds/"+RandomString+"Backgrounds.root").c_str());

    TFile* BackgroundsF = TFile::Open(backgroundC);
    
    for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
    {
        NearBackgroundSpectrumH[near]=(TH1D*)gDirectory->Get((Form("Near AD%i ",near)+RandomString+ Form(" Background Period%d",week)).c_str());
    }
    for (Int_t far = 0; far < ADsEH3; far++)
    {
        FarBackgroundSpectrumH[far]=(TH1D*)gDirectory->Get((Form("Far AD%i ",far)+RandomString+ Form(" Background Period%d",week)).c_str());
    }
    
    BackgroundsF->Close();
    
    for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
    {
        CombinedNearBackgroundSpectrumH[near]=(TH1D*)NearBackgroundSpectrumH[near]->Clone();
        
        if(Combine!=0)
        {
            CombinedNearBackgroundSpectrumH[near]->Reset();
        }
    }
    for (Int_t far = 0; far < ADsEH3; far++)
    {
        CombinedFarBackgroundSpectrumH[far]=(TH1D*)FarBackgroundSpectrumH[far]->Clone();
        
        if(Combine!=0)
        {
            CombinedFarBackgroundSpectrumH[far]->Reset();
        }
    }
    
    //Combine matrices in 9x9, 2x2 or 1x1 prediction
    if (Combine == 1)
    {
        for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            CombinedNearBackgroundSpectrumH[0]->Add(NearBackgroundSpectrumH[near]);//All near hall detectors together
        }
        
        for (Int_t far = 0; far < ADsEH3; far++)
        {
            CombinedFarBackgroundSpectrumH[0]->Add(FarBackgroundSpectrumH[far]);//All far hall detectors together
        }
        
        //            CombinedNearBackgroundSpectrumH[0]->Scale(1./(ADsEH1+ADsEH2));//All near hall detectors together
    }
    else if(Combine == 2)
    {
        for (Int_t far = 0; far < ADsEH3; far++)
        {
            CombinedFarBackgroundSpectrumH[0]->Add(FarBackgroundSpectrumH[far]);//All far hall detectors together
        }
        
        //DB
        for(Int_t FirstHall = 0; FirstHall<ADsEH1;FirstHall++ )
        {
            CombinedNearBackgroundSpectrumH[0]->Add(NearBackgroundSpectrumH[FirstHall]);
            
        }
        //            CombinedNearBackgroundSpectrumH[0]->Scale(1./ADsEH1);
        
        //LO
        for(Int_t SecondHall = ADsEH1; SecondHall<ADsEH2;SecondHall++ )
        {
            CombinedNearBackgroundSpectrumH[1]->Add(NearBackgroundSpectrumH[SecondHall]);
            
        }
        //            CombinedNearBackgroundSpectrumH[1]->Scale(1./ADsEH2);
        
    }
    
    #ifdef PrintEps
    TCanvas* FarBackgroundsC = new TCanvas("FarBackgroundsC","FarBackgroundsC",MaxFarLoadOscModel*400,400);
    TCanvas* NearBackgroundsC = new TCanvas("NearBackgroundsC","NearBackgroundsC",(ADsEH1+ADsEH2)*400,400);

    FarBackgroundsC->Divide(ADsEH3);
    NearBackgroundsC->Divide(ADsEH1+ADsEH2);
        
    for (Int_t ad=0; ad<NADs/2; ad++)
    {
        FarBackgroundsC->cd(ad+1);
        FarBackgroundSpectrumH[ad]->Draw("HIST");

        NearBackgroundsC->cd(ad+1);
        NearBackgroundSpectrumH[ad]->Draw("HIST");
    }
    
    FarBackgroundsC->Print(("./Images/"+AnalysisString+"/FitterInputs/"+RandomString+"FarBackgrounds.eps").c_str(),".eps");
    NearBackgroundsC->Print(("./Images/"+AnalysisString+"/FitterInputs/"+RandomString+"NearBackgrounds.eps").c_str(),".eps");

    delete FarBackgroundsC;
    delete NearBackgroundsC;
    #endif
}

void Prediction :: LoadRootCovarianceMatrices(Int_t week)
{
    Char_t RootC[100];
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/VaryAccidental.root",Combine)).c_str());
    TFile* VAccCovarianceMatrixF = new TFile(RootC);
    VAccCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    VAccCovarianceMatrixH->SetName("Vacc Matrix");
    VAccCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/VaryFastNeutrons.root",Combine)).c_str());
    TFile* VFNCovarianceMatrixF = new TFile(RootC);
    VFNCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    VFNCovarianceMatrixH->SetName("VFN Matrix");
    VFNCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/VaryLiHe.root",Combine)).c_str());
    TFile* VLiHeCovarianceMatrixF = new TFile(RootC);
    VLiHeCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    VLiHeCovarianceMatrixH->SetName("VLiHe Matrix");
    VLiHeCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/VaryAmC.root",Combine)).c_str());
    TFile* VAmCCovarianceMatrixF = new TFile(RootC);
    VAmCCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    VAmCCovarianceMatrixH->SetName("VAmC Matrix");
    VAmCCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/DistortFastNeutrons.root",Combine)).c_str());
    TFile* DFNCovarianceMatrixF = new TFile(RootC);
    DFNCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    DFNCovarianceMatrixH->SetName("DFN Matrix");
    DFNCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/DistortLiHe.root",Combine)).c_str());
    TFile* DLiHeCovarianceMatrixF = new TFile(RootC);
    DLiHeCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    DLiHeCovarianceMatrixH->SetName("DLiHe Matrix");
    DLiHeCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/DistortAmC.root",Combine)).c_str());
    TFile* DAmCCovarianceMatrixF = new TFile(RootC);
    DAmCCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    DAmCCovarianceMatrixH->SetName("DAmC Matrix");
    DAmCCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/Isotope.root",Combine)).c_str());
    TFile* IsotopeCovarianceMatrixF = new TFile(RootC);
    IsotopeCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    IsotopeCovarianceMatrixH->SetName("Isotope Matrix");
    IsotopeCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/ReactorPower.root",Combine)).c_str());
    TFile* ReactorPowerCovarianceMatrixF = new TFile(RootC);
    ReactorPowerCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    ReactorPowerCovarianceMatrixH->SetName("Reactor Matrix");
    ReactorPowerCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/RelativeEnergyScale.root",Combine)).c_str());
    TFile* RelativeEnergyScaleCovarianceMatrixF = new TFile(RootC);
    RelativeEnergyScaleCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    RelativeEnergyScaleCovarianceMatrixH->SetName("Relative Energy Matrix");
    RelativeEnergyScaleCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/IAV.root",Combine)).c_str());
    TFile* IAVCovarianceMatrixF = new TFile(RootC);
    IAVCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    IAVCovarianceMatrixH->SetName("IAV Matrix");
    IAVCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/NL.root",Combine)).c_str());
    TFile* NLCovarianceMatrixF = new TFile(RootC);
    NLCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    NLCovarianceMatrixH->SetName("NL Matrix");
    NLCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/Resolution.root",Combine)).c_str());
    TFile* ResolutionCovarianceMatrixF = new TFile(RootC);
    ResolutionCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    ResolutionCovarianceMatrixH->SetName("Reso Matrix");
    ResolutionCovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/Sin22t12.root",Combine)).c_str());
    TFile* Sin22t12CovarianceMatrixF = new TFile(RootC);
    Sin22t12CovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    Sin22t12CovarianceMatrixH->SetName("Sin22t12 Matrix");
    Sin22t12CovarianceMatrixF->Close();
    
    sprintf(RootC,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/Efficiency.root",Combine)).c_str());
    TFile* EfficiencyCovarianceMatrixF = new TFile(RootC);
    EfficiencyCovarianceMatrixH = (TH2D*)gDirectory->Get(Form("Covariance Matrix%d",week));
    EfficiencyCovarianceMatrixH->SetName("Efficiency Matrix");
    EfficiencyCovarianceMatrixF->Close();
    
    // Just to check that they have been loaded properly
    TFile* SaveSystematicCovMatrixF = TFile::Open(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/SystematicCovarianceMatrixBeforeNormalize.root",Combine)).c_str(),"recreate");
    IsotopeCovarianceMatrixH->Write();
    ReactorPowerCovarianceMatrixH->Write();
    RelativeEnergyScaleCovarianceMatrixH->Write();
    IAVCovarianceMatrixH->Write();
    NLCovarianceMatrixH->Write();
    ResolutionCovarianceMatrixH->Write();
    Sin22t12CovarianceMatrixH->Write();
    EfficiencyCovarianceMatrixH->Write();
    SaveSystematicCovMatrixF->Close();
}

void Prediction :: LoadTxtCovarianceMatrices(Int_t week)
{
    Char_t filenameCov[100];
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/VaryAccidentalPeriod%d.txt",Combine,week)).c_str());
    
    ifstream covfile_vacc(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j <MaxBins; j++)
        {
            covfile_vacc >>  VAccCovarianceMatrixM[i][j];
        }
    }
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/VaryFastNeutronsPeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_vfn(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j <MaxBins; j++)
        {
            covfile_vfn >>  VFNCovarianceMatrixM[i][j];
        }
    }
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/VaryLiHePeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_vlihe(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j <MaxBins; j++)
        {
            covfile_vlihe >>  VLiHeCovarianceMatrixM[i][j];
        }
    }
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/VaryAmCPeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_vamc(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j <MaxBins; j++)
        {
            covfile_vamc >>  VAmCCovarianceMatrixM[i][j];
        }
    }
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/DistortFastNeutronsPeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_dfn(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j <MaxBins; j++)
        {
            covfile_dfn >>  DFNCovarianceMatrixM[i][j];
        }
    }
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/DistortLiHePeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_dlihe(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j < MaxBins; j++)
        {
            covfile_dlihe >>  DLiHeCovarianceMatrixM[i][j];
        }
    }
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/DistortAmCPeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_damc(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j < MaxBins; j++)
        {
            covfile_damc >>  DAmCCovarianceMatrixM[i][j];
        }
    }
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/IsotopePeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_iso(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j < MaxBins; j++)
        {
            covfile_iso >>  IsotopeCovarianceMatrixM[i][j];
        }
    }
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/ReactorPowerPeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_pow(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j <MaxBins; j++)
        {
            covfile_pow >> ReactorPowerCovarianceMatrixM[i][j];
        }
    }
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/IAVPeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_iav(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j <MaxBins; j++)
        {
            covfile_iav >> IAVCovarianceMatrixM[i][j];
        }
    }
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/NLPeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_nl(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j <MaxBins; j++)
        {
            covfile_nl >> NLCovarianceMatrixM[i][j];
        }
    }
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/ResolutionPeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_reso(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j <MaxBins; j++)
        {
            covfile_reso >> ResolutionCovarianceMatrixM[i][j];
        }
    }
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/RelativeEnergyScalePeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_relesc(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j <MaxBins; j++)
        {
            covfile_relesc >> RelativeEnergyScaleCovarianceMatrixM[i][j];
        }
    }
    
    //    sprintf(filenameCov,"./CovarianceMatrices/Combine%d/CovarianceMatricesTxT/RelativeEnergyOffsetPeriod%d.txt",Combine,week);
    //    ifstream covfile_reloff(filenameCov);
    //    for (Int_t i = 0; i < MaxBins; i++)
    //    {
    //        for (Int_t j = 0; j <MaxBins; j++)
    //        {
    //            covfile_reloff >> RelativeEnergyOffsetCovarianceMatrixM[i][j];
    //        }
    //    }
    //
    //    sprintf(filenameCov,"./CovarianceMatrices/Combine%d/CovarianceMatricesTxT/AbsoluteEnergyScalePeriod%d.txt",Combine,week);
    //    ifstream covfile_absrel(filenameCov);
    //    for (Int_t i = 0; i < MaxBins; i++)
    //    {
    //        for (Int_t j = 0; j <MaxBins; j++)
    //        {
    //            covfile_absrel >>  AbsoluteEnergyScaleCovarianceMatrixM[i][j];
    //        }
    //    }
    //    sprintf(filenameCov,"./CovarianceMatrices/Combine%d/CovarianceMatricesTxT/AbsoluteEnergyOffsetPeriod%d.txt",Combine,week);
    //    ifstream covfile_absoff(filenameCov);
    //    for (Int_t i = 0; i < MaxBins; i++)
    //    {
    //        for (Int_t j = 0; j <MaxBins; j++)
    //        {
    //            covfile_absoff >>  AbsoluteEnergyOffsetCovarianceMatrixM[i][j];
    //        }
    //    }
    //
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/Sin22t12Period%d.txt",Combine,week)).c_str());
    ifstream covfile_sin22t12(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j <MaxBins; j++)
        {
            covfile_sin22t12 >>  Sin22t12CovarianceMatrixM[i][j];
        }
    }
    
    sprintf(filenameCov,("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesTxT/EfficiencyPeriod%d.txt",Combine,week)).c_str());
    ifstream covfile_eff(filenameCov);
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j <MaxBins; j++)
        {
            covfile_eff >> EfficiencyCovarianceMatrixM[i][j];
        }
    }
    
    //  Initialize histograms
    
    VAccCovarianceMatrixH = new TH2D("Vary Accidentals Covariance Matrix","Vary Accidentals Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    VLiHeCovarianceMatrixH = new TH2D("Vary LiHe Covariance Matrix","Vary LiHe Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    VAmCCovarianceMatrixH = new TH2D("Vary AmC Covariance Matrix","Vary AmC Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    DFNCovarianceMatrixH = new TH2D("Distort FN Covariance Matrix","Distort FN Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    DLiHeCovarianceMatrixH = new TH2D("Distort LiHe Covariance Matrix","Distort LiHe Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    DAmCCovarianceMatrixH = new TH2D("Distort AmC Covariance Matrix","Distort AmC Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    
    IsotopeCovarianceMatrixH = new TH2D("Isotope Covariance Matrix","Isotope Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    ReactorPowerCovarianceMatrixH = new TH2D("Reactor Power Covariance Matrix","Reactor Power Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    RelativeEnergyScaleCovarianceMatrixH = new TH2D("Relative Energy Scale Covariance Matrix","Relative Energy Scale Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    //    RelativeEnergyOffsetCovarianceMatrixH = new TH2D("Relative Energy Offset Covariance Matrix","Relative Offset Scale Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    //    AbsoluteEnergyScaleCovarianceMatrixH = new TH2D("Absolute Energy Scale Covariance Matrix","Absolute Energy Scale Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    //    AbsoluteEnergyOffsetCovarianceMatrixH = new TH2D("Absolute Energy Offset Covariance Matrix","Absolute Offset Scale Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    IAVCovarianceMatrixH = new TH2D("IAV Covariance Matrix","IAV Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    NLCovarianceMatrixH = new TH2D("NL Covariance Matrix","NL Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    ResolutionCovarianceMatrixH = new TH2D("Resolution Covariance Matrix","Resolution Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    Sin22t12CovarianceMatrixH = new TH2D("Sin22t12 Covariance Matrix","Sin22t12 Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    EfficiencyCovarianceMatrixH = new TH2D("Efficiency Covariance Matrix","Efficiency Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    
}

void Prediction :: SaveCovarianceMatrices(Int_t week)
{
    std::cout << "SAVING COVARIANCE MATRICES " << std::endl;
    if(ReadTxt)
    {
        for (Int_t i = 0; i<MaxBins; i++)
        {
            for (Int_t j = 0; j<MaxBins; j++)
            {
                RenormIsotopeCovarianceMatrixH->SetBinContent(i+1,j+1,RenormIsotopeCovarianceMatrixM[i][j]);
                RenormReactorPowerCovarianceMatrixH->SetBinContent(i+1,j+1,RenormReactorPowerCovarianceMatrixM[i][j]);
                RenormRelativeEnergyScaleCovarianceMatrixH->SetBinContent(i+1,j+1,RenormRelativeEnergyScaleCovarianceMatrixM[i][j]);
                //                RelativeEnergyOffsetCovarianceMatrixH->SetBinContent(i+1,j+1,RelativeEnergyOffsetCovarianceMatrixM[i][j]);//Not used
                //                AbsoluteEnergyScaleCovarianceMatrixH->SetBinContent(i+1,j+1,AbsoluteEnergyScaleCovarianceMatrixM[i][j]);
                //                AbsoluteEnergyOffsetCovarianceMatrixH->SetBinContent(i+1,j+1,AbsoluteEnergyOffsetCovarianceMatrixM[i][j]);//Included in NL
                RenormIAVCovarianceMatrixH->SetBinContent(i+1,j+1,RenormIAVCovarianceMatrixM[i][j]);
                RenormNLCovarianceMatrixH->SetBinContent(i+1,j+1,RenormNLCovarianceMatrixM[i][j]);
                RenormResolutionCovarianceMatrixH->SetBinContent(i+1,j+1,RenormResolutionCovarianceMatrixM[i][j]);
                RenormSin22t12CovarianceMatrixH->SetBinContent(i+1,j+1,RenormSin22t12CovarianceMatrixM[i][j]);
                RenormEfficiencyCovarianceMatrixH->SetBinContent(i+1,j+1,RenormEfficiencyCovarianceMatrixM[i][j]);
                
                VAccCovarianceMatrixH->SetBinContent(i+1,j+1,VAccCovarianceMatrixM[i][j]);
                VLiHeCovarianceMatrixH->SetBinContent(i+1,j+1,VLiHeCovarianceMatrixM[i][j]);
                VFNCovarianceMatrixH->SetBinContent(i+1,j+1,VFNCovarianceMatrixM[i][j]);
                VAmCCovarianceMatrixH->SetBinContent(i+1,j+1,VAmCCovarianceMatrixM[i][j]);
                DLiHeCovarianceMatrixH->SetBinContent(i+1,j+1,DLiHeCovarianceMatrixM[i][j]);
                DFNCovarianceMatrixH->SetBinContent(i+1,j+1,DFNCovarianceMatrixM[i][j]);
                DAmCCovarianceMatrixH->SetBinContent(i+1,j+1,DAmCCovarianceMatrixM[i][j]);
                
                SystematicCovarianceMatrixH->SetBinContent(i+1,j+1,SystematicCovarianceMatrixM[i][j]);
                BackgroundsCovarianceMatrixH->SetBinContent(i+1,j+1,BackgroundsCovarianceMatrixM[i][j]);
                TotalCovarianceMatrixH->SetBinContent(i+1,j+1,TotalCovarianceMatrixM[j+i*MaxBins]);
            }
        }
    }
    
    TFile* SaveCovarianceMatricesF = new TFile(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesFitterPeriod%d.root",Combine,week)).c_str(),"recreate");
    VAccCovarianceMatrixH->Write("Vacc Matrix");
    VFNCovarianceMatrixH->Write("VFN Matrix");
    VLiHeCovarianceMatrixH->Write("VLiHe Matrix");
    VAmCCovarianceMatrixH->Write("VAmC Matrix");
    DFNCovarianceMatrixH->Write("DFN Matrix");
    DLiHeCovarianceMatrixH->Write("DLiHe Matrix");
    DAmCCovarianceMatrixH->Write("DAmC Matrix");
    
    RenormIsotopeCovarianceMatrixH->Write("Isotope Matrix");
    RenormReactorPowerCovarianceMatrixH->Write("Power Matrix");
    RenormRelativeEnergyScaleCovarianceMatrixH->Write("Relative Scale Matrix");
    //    RelativeEnergyOffsetCovarianceMatrixH->Write("Relative Offset Matrix");
    //    AbsoluteEnergyScaleCovarianceMatrixH->Write("Absolute Scale Matrix");
    //    AbsoluteEnergyOffsetCovarianceMatrixH->Write("Absolute Offset Matrix");
    RenormIAVCovarianceMatrixH->Write("IAV Matrix");
    RenormNLCovarianceMatrixH->Write("NL Matrix");
    RenormResolutionCovarianceMatrixH->Write("Reso Matrix");
    RenormSin22t12CovarianceMatrixH->Write("Sin22t12 Matrix");
    RenormEfficiencyCovarianceMatrixH->Write("Efficiency Matrix");
    
    BackgroundsCovarianceMatrixH->Write("Background Covariance Matrix");
    SystematicCovarianceMatrixH->Write("Systematic Covariance Matrix");
    StatisticalCovarianceMatrixH->Write("Statistical Covariance Matrix");
    TotalCovarianceMatrixH->Write("Total Covariance Matrix");
    InvTotalCovarianceMatrixH->Write("Inv Matrix");
    UnityH->Write("Unity Matrix");
    SaveCovarianceMatricesF->Close();
    
    if(PlotCorrelation)
    {
        TH2D* BackgroundsCorrelationMatrixH = NormCov(BackgroundsCovarianceMatrixH,BackgroundsCorrelationMatrixH);
        TH2D* SystematicCorrelationMatrixH = NormCov(SystematicCovarianceMatrixH,SystematicCorrelationMatrixH);
        TH2D* TotalCorrelationMatrixH = NormCov(TotalCovarianceMatrixH,TotalCorrelationMatrixH);
        TH2D* IsotopeCorrelationMatrixH = NormCov(RenormIsotopeCovarianceMatrixH,IsotopeCorrelationMatrixH);
        TH2D* ReactorPowerCorrelationMatrixH = NormCov(RenormReactorPowerCovarianceMatrixH,ReactorPowerCorrelationMatrixH);
        TH2D* RelativeEnergyScaleCorrelationMatrixH = NormCov(RenormRelativeEnergyScaleCovarianceMatrixH,RelativeEnergyScaleCorrelationMatrixH);
        TH2D* IAVCorrelationMatrixH = NormCov(RenormIAVCovarianceMatrixH,IAVCorrelationMatrixH);
        TH2D* NLCorrelationMatrixH = NormCov(RenormNLCovarianceMatrixH,NLCorrelationMatrixH);
        TH2D* ResolutionCorrelationMatrixH = NormCov(RenormResolutionCovarianceMatrixH,ResolutionCorrelationMatrixH);
        TH2D* Sin22t12CorrelationMatrixH = NormCov(RenormSin22t12CovarianceMatrixH,Sin22t12CorrelationMatrixH);
        TH2D* EfficiencyCorrelationMatrixH = NormCov(RenormEfficiencyCovarianceMatrixH,EfficiencyCorrelationMatrixH);
        TH2D* VAccCorrelationMatrixH = NormCov(VAccCovarianceMatrixH,VAccCorrelationMatrixH);
        TH2D* VFNCorrelationMatrixH = NormCov(VFNCovarianceMatrixH,VFNCorrelationMatrixH);
        TH2D* VLiHeCorrelationMatrixH = NormCov(VLiHeCovarianceMatrixH,VLiHeCorrelationMatrixH);
        TH2D* VAmCCorrelationMatrixH = NormCov(VAmCCovarianceMatrixH,VAmCCorrelationMatrixH);
        TH2D* DFNCorrelationMatrixH = NormCov(DFNCovarianceMatrixH,DFNCorrelationMatrixH);
        TH2D* DLiHeCorrelationMatrixH = NormCov(DLiHeCovarianceMatrixH,DLiHeCorrelationMatrixH);
        TH2D* DAmCCorrelationMatrixH = NormCov(DAmCCovarianceMatrixH,DAmCCorrelationMatrixH);
        TH2D* StatisticalCorrelationMatrixH = NormCov(StatisticalCovarianceMatrixH,StatisticalCorrelationMatrixH);
        //  To plot correlation matrices: (Not used in the code but used to show results)
        
        TFile* SaveNormF = TFile::Open(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CorrelationMatrices.root",Combine)).c_str(),"recreate");
        
        BackgroundsCorrelationMatrixH->Write();
        SystematicCorrelationMatrixH->Write();
        StatisticalCorrelationMatrixH->Write("Statistical Covariance Matrix");
        TotalCorrelationMatrixH->Write();
        
        VAccCorrelationMatrixH->Write("Vacc Matrix");
        VFNCorrelationMatrixH->Write("VFN Matrix");
        VLiHeCorrelationMatrixH->Write("VLiHe Matrix");
        VAmCCorrelationMatrixH->Write("VAmC Matrix");
        DFNCorrelationMatrixH->Write("DFN Matrix");
        DLiHeCorrelationMatrixH->Write("DLiHe Matrix");
        DAmCCorrelationMatrixH->Write("DAmC Matrix");
        
        IsotopeCorrelationMatrixH->Write("Isotope Matrix");
        ReactorPowerCorrelationMatrixH->Write("Power Matrix");
        RelativeEnergyScaleCorrelationMatrixH->Write("Relative Scale Matrix");
        IAVCorrelationMatrixH->Write("IAV Matrix");
        NLCorrelationMatrixH->Write("NL Matrix");
        ResolutionCorrelationMatrixH->Write("Reso Matrix");
        Sin22t12CorrelationMatrixH->Write("Sin22t12 Matrix");
        EfficiencyCorrelationMatrixH->Write("Efficiency Matrix");
        
        SaveNormF->Close();
        
        delete BackgroundsCorrelationMatrixH;
        delete SystematicCorrelationMatrixH;
        delete TotalCorrelationMatrixH;
        delete IsotopeCorrelationMatrixH;
        delete ReactorPowerCorrelationMatrixH;
        delete RelativeEnergyScaleCorrelationMatrixH;
        delete IAVCorrelationMatrixH;
        delete NLCorrelationMatrixH;
        delete ResolutionCorrelationMatrixH;
        delete Sin22t12CorrelationMatrixH;
        delete EfficiencyCorrelationMatrixH;
        delete VAccCorrelationMatrixH;
        delete VFNCorrelationMatrixH;
        delete VLiHeCorrelationMatrixH;
        delete VAmCCorrelationMatrixH;
        delete DFNCorrelationMatrixH;
        delete DLiHeCorrelationMatrixH;
        delete DAmCCorrelationMatrixH;
        delete StatisticalCorrelationMatrixH;
    }
    
    if(WriteOutput)
    {
        ofstream bkgf(("CovarianceMatrices/"+AnalysisString+Form("/Combine%d/BackgroundCovarianceMatrix.txt",Combine)).c_str());
        ofstream sysf(("CovarianceMatrices/"+AnalysisString+Form("/Combine%d/SystematicCovarianceMatrix.txt",Combine)).c_str());
        ofstream totalf(("CovarianceMatrices/"+AnalysisString+Form("/Combine%d/TotalCovarianceMatrix.txt",Combine)).c_str());
        
        Int_t x =0;
        Int_t y =0;
        
        Int_t Ni1=0,Ni2=0,Ni3=0,Ni4=0;
        Int_t Nj1=0,Nj2=0,Nj3=0,Nj4=0;
        Int_t Fi1=0,Fi2=0,Fi3=0,Fi4=0;
        Int_t Fj1=0,Fj2=0,Fj3=0,Fj4=0;
        
        for (Int_t neari=0; neari<MaxNearCombine; neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Ni1=1;Ni2=0;Ni3=0;Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<MaxFarCombine; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Fi1=fari+1;Fi2=0;Fi3=0;Fi4=0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFarCombine;Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFarCombine;Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=MaxFarCombine;Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<MaxNearCombine; nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<MaxFarCombine; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Fj1=farj+1;Fj2=0;Fj3=0;Fj4=0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFarCombine; Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFarCombine; Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=MaxFarCombine; Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                            
                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                
                                bkgf << BackgroundsCovarianceMatrixH->GetBinContent(x+1,y+1) << " ";
                                sysf << SystematicCovarianceMatrixH->GetBinContent(x+1,y+1) << " ";
                                totalf << TotalCovarianceMatrixH->GetBinContent(x+1,y+1) << " ";
                                
                            }
                        }
                    }
                }
            }
        }
        bkgf << std::endl;
        bkgf.close();
        sysf << std::endl;
        sysf.close();
        totalf << std::endl;
        totalf.close();
    }
    std::cout << "\t Finished Saving Covariance Matrices" << std::endl;
    
    for(Int_t i=0; i<MaxBins; i++)
    {
        for(Int_t j=0; j<MaxBins; j++)
        {
            if(i==j)
            {
                if(TMath::Abs(UnityH->GetBinContent(i+1,j+1)-1)>=0.000000001)
                {
                    std::cout << " INVERSE HAS NOT BEEN CALCULATED PROPERLY:" << UnityH->GetBinContent(i+1,j+1) << std::endl;
                    
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
}

void Prediction :: CombineMatrices(Int_t week)
{
    std::cout <<  "***********************************************************************************************" << std::endl;
    std::cout <<  "\t Combining Matrices " << std::endl;
    
    BackgroundsCovarianceMatrixH = new TH2D("Background Covariance Matrix","Background Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    SystematicCovarianceMatrixH = new TH2D("Systematic Covariance Matrix","Systematic Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    TotalCovarianceMatrixH = new TH2D("Total Covariance Matrix","Total Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    
    InvTotalCovarianceMatrixH = new TH2D("Inv Total Covariance Matrix","Inv Total Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    
    RenormIsotopeCovarianceMatrixH = new TH2D("Renorm Isotope Covariance Matrix","Renorm Isotope Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    RenormReactorPowerCovarianceMatrixH = new TH2D("Renorm Reactor Power Covariance Matrix","Renorm Reactor Power Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    RenormRelativeEnergyScaleCovarianceMatrixH = new TH2D("Renorm Relative Energy Scale Covariance Matrix","Renorm Relative Energy Scale Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    //    RenormRelativeEnergyOffsetCovarianceMatrixH = new TH2D("Renorm Relative Energy Offset Covariance Matrix","Renorm Relative Energy Offset Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    //    RenormAbsoluteEnergyScaleCovarianceMatrixH = new TH2D("Renorm Absolute Energy Scale Covariance Matrix","Renorm Absolute Energy Scale Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    //    RenormAbsoluteEnergyOffsetCovarianceMatrixH = new TH2D("Renorm Absolute Energy Offset Covariance Matrix","Renorm Absolute Energy Offset Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    RenormIAVCovarianceMatrixH = new TH2D("Renorm IAV Covariance Matrix","Renorm IAV Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    RenormNLCovarianceMatrixH = new TH2D("Renorm NL Covariance Matrix","Renorm NL Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    RenormResolutionCovarianceMatrixH = new TH2D("Renorm Resolution Covariance Matrix","Renorm Resolution Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    RenormSin22t12CovarianceMatrixH = new TH2D("Renorm Sin22t12 Covariance Matrix","Renorm Sin22t12 Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    RenormEfficiencyCovarianceMatrixH = new TH2D("Renorm Efficiency Covariance Matrix","Renorm Efficiency Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    
    std::cout <<  "\t Matrices Created " << std::endl;

    
    
    
    
    std::cout << SysCovDirectory << BkgCovDirectory << std::endl;
    //Test LBNL inputs, by loading external covariance matrices:
    if(strcmp((SysCovDirectory).c_str(),"")&&strcmp((BkgCovDirectory).c_str(),""))
    {
        TH2D* BigBackgroundsCovarianceMatrixH;
        TH2D* BigSystematicCovarianceMatrixH;
        
        std::cout << "\t USING LBNL COV MATRICES" << std::endl;
        TFile* LoadBkgdF = TFile::Open(BkgCovDirectory.c_str());
        BigBackgroundsCovarianceMatrixH = (TH2D*)gDirectory->Get("h_covmatrix");
        LoadBkgdF->Close();
        
        TFile* LoadSystematicF = TFile::Open(SysCovDirectory.c_str());
        BigSystematicCovarianceMatrixH= (TH2D*)gDirectory->Get("h_covmatrix");
        LoadSystematicF->Close();
        
        //need to rebin 9x9 matrices to 2x2 or 1x1:
        Double_t Xscale, Yscale;
        
        for(Int_t x = 1; x<=9;x++)
        {
            for(Int_t y = 1; y<=9;y++)
            {
                for(Int_t i = 0; i<MaxBins; i++)
                {
                    for(Int_t j = 0; j<MaxBins; j++)
                    {
                        if(Combine==2)
                        {
                            if(x<=6)
                            {
                                Xscale = 1./2;
                            }
                            else
                            {
                                Xscale = 1.;
                            }
                            if(y<=6)
                            {
                                Yscale = 1./2;
                            }
                            else
                            {
                                Yscale = 1.;
                            }
                            
                            BackgroundsCovarianceMatrixH->SetBinContent(i+1,j+1,BackgroundsCovarianceMatrixH->GetBinContent(i+1,j+1)+Xscale*Yscale*BigBackgroundsCovarianceMatrixH->GetBinContent(i*x+1,j*y+1));
                            SystematicCovarianceMatrixH->SetBinContent(i+1,j+1,SystematicCovarianceMatrixH->GetBinContent(i+1,j+1)+Xscale*Yscale*BigSystematicCovarianceMatrixH->GetBinContent(i*x+1,j*y+1));
                            
                        }
                        else if(Combine == 1)
                        {
                            BackgroundsCovarianceMatrixH->SetBinContent(i+1,j+1,BackgroundsCovarianceMatrixH->GetBinContent(i+1,j+1)+BigBackgroundsCovarianceMatrixH->GetBinContent(i*x+1,j*y+1));
                            SystematicCovarianceMatrixH->SetBinContent(i+1,j+1,SystematicCovarianceMatrixH->GetBinContent(i+1,j+1)+BigSystematicCovarianceMatrixH->GetBinContent(i*x+1,j*y+1));
                        }
                        else
                        {
                            BackgroundsCovarianceMatrixH->SetBinContent(i+1,j+1,BigBackgroundsCovarianceMatrixH->GetBinContent(i+1,j+1));
                            SystematicCovarianceMatrixH->SetBinContent(i+1,j+1,BigSystematicCovarianceMatrixH->GetBinContent(i+1,j+1));
                        }
                    }
                }
            }
        }
        
        delete BigBackgroundsCovarianceMatrixH;
        delete BigSystematicCovarianceMatrixH;
        
        //Rescale Systematic Covariance Matrix:
        
        Int_t x =0;
        Int_t y =0;
        Int_t Ni1=0,Ni2=0,Ni3=0,Ni4=0;
        Int_t Nj1=0,Nj2=0,Nj3=0,Nj4=0;
        Int_t Fi1=0,Fi2=0,Fi3=0,Fi4=0;
        Int_t Fj1=0,Fj2=0,Fj3=0,Fj4=0;
        
        for (Int_t neari=0; neari<MaxNearCombine; neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Ni1=1;Ni2=0;Ni3=0;Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<MaxFarCombine; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Fi1=fari+1;Fi2=0;Fi3=0;Fi4=0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFarCombine; Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFarCombine; Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=MaxFarCombine; Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<MaxNearCombine; nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<MaxFarCombine; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Fj1=farj+1;Fj2=0;Fj3=0;Fj4=0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFarCombine; Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFarCombine; Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=MaxFarCombine; Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                            
                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                
                                SystematicCovarianceMatrixH->SetBinContent(x+1,y+1,SystematicCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1)));
                            }
                        }
                    }
                }
            }
        }
        
        TotalCovarianceMatrixH->Add(StatisticalCovarianceMatrixH);
        TotalCovarianceMatrixH->Add(BackgroundsCovarianceMatrixH);
        TotalCovarianceMatrixH->Add(SystematicCovarianceMatrixH);
    }
    else
    {
        Int_t x =0;
        Int_t y =0;
        Int_t Ni1=0,Ni2=0,Ni3=0,Ni4=0;
        Int_t Nj1=0,Nj2=0,Nj3=0,Nj4=0;
        Int_t Fi1=0,Fi2=0,Fi3=0,Fi4=0;
        Int_t Fj1=0,Fj2=0,Fj3=0,Fj4=0;
        
        for (Int_t neari=0; neari<MaxNearCombine; neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Ni1=1;Ni2=0;Ni3=0;Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<MaxFarCombine; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Fi1=fari+1;Fi2=0;Fi3=0;Fi4=0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFarCombine; Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFarCombine; Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=MaxFarCombine; Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<MaxNearCombine; nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<MaxFarCombine; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Fj1=farj+1;Fj2=0;Fj3=0;Fj4=0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFarCombine; Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFarCombine; Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=MaxFarCombine; Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                            
                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                
                                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                //                                                              Denormalize Oscillation in Systematic Covariance Matrix
                                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                
                                if(ReadTxt)
                                {
                                    RenormIsotopeCovarianceMatrixM[x][y]=IsotopeCovarianceMatrixM[x][y]*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1));
                                    
                                    RenormReactorPowerCovarianceMatrixM[x][y]=ReactorPowerCovarianceMatrixM[x][y]*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1));
                                    
                                    RenormRelativeEnergyScaleCovarianceMatrixM[x][y]=RenormRelativeEnergyScaleCovarianceMatrixM[x][y]*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1));
                                    
                                    //                                RenormRelativeEnergyOffsetCovarianceMatrixM[x][y]=RenormRelativeEnergyOffsetCovarianceMatrixM[x][y]*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1));
                                    
                                    //                                RenormAbsoluteEnergyScaleCovarianceMatrixM[x][y]=RenormAbsoluteEnergyScaleCovarianceMatrixM[x][y]*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1));
                                    
                                    //                                RenormAbsoluteEnergyOffsetCovarianceMatrixM[x][y]=RenormAbsoluteEnergyOffsetCovarianceMatrixM[x][y]*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1));
                                    
                                    RenormIAVCovarianceMatrixM[x][y]=IAVCovarianceMatrixM[x][y]*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1));
                                    
                                    RenormNLCovarianceMatrixM[x][y]=NLCovarianceMatrixM[x][y]*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1));
                                    
                                    RenormResolutionCovarianceMatrixM[x][y]=ResolutionCovarianceMatrixM[x][y]*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1));
                                    
                                    RenormSin22t12CovarianceMatrixM[x][y]=Sin22t12CovarianceMatrixM[x][y]*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1));
                                    
                                    RenormEfficiencyCovarianceMatrixM[x][y]=EfficiencyCovarianceMatrixM[x][y]*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1));
                                    // Normal Fit
                                    
                                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                    //                                                             Add Systematic Covariance Matrix
                                    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                    SystematicCovarianceMatrixM[x][y]=
                                    RenormIsotopeCovarianceMatrixM[x][y]+
                                    RenormReactorPowerCovarianceMatrixM[x][y]+
                                    RenormRelativeEnergyScaleCovarianceMatrixM[x][y]+
                                    //                                RenormRelativeEnergyOffsetCovarianceMatrixM[x][y]+
                                    //                                RenormAbsoluteEnergyScaleCovarianceMatrixM[x][y]+
                                    //                                RenormAbsoluteEnergyOffsetCovarianceMatrixM[x][y]+
                                    RenormIAVCovarianceMatrixM[x][y]+
                                    RenormNLCovarianceMatrixM[x][y]+
                                    RenormResolutionCovarianceMatrixM[x][y]+
                                    RenormSin22t12CovarianceMatrixM[x][y]+
                                    RenormEfficiencyCovarianceMatrixM[x][y];
                                    
                                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                    //                                                             Add Background Covariance Matrix
                                    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                    BackgroundsCovarianceMatrixM[x][y]=
                                    VAccCovarianceMatrixM[x][y]+
                                    VLiHeCovarianceMatrixM[x][y]+
                                    VFNCovarianceMatrixM[x][y]+
                                    VAmCCovarianceMatrixM[x][y]+
                                    DLiHeCovarianceMatrixM[x][y]+
                                    DFNCovarianceMatrixM[x][y]+
                                    DAmCCovarianceMatrixM[x][y];
                                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                    //                                                      Add all matrices into a Total Covariance Matrix
                                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                    
                                    TotalCovarianceMatrixM[x*MaxBins+y]=CovStat[x][y]+BackgroundsCovarianceMatrixM[x][y]+SystematicCovarianceMatrixM[x][y];
                                    
                                    // Specific combinations to produce the error budget:
                                    if(BudgetTurnOff)
                                    {
                                        //  Substract the choosen systematic
                                        
                                        if(Data->GetVaryAccidentalBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-VAccCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetVaryLiHeBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-VLiHeCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetVaryFastNeutronsBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-VFNCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetVaryAmCBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-VAmCCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetDistortLiHeBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-DLiHeCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetDistortFastNeutronsBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-DFNCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetDistortAmCBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-DAmCCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetIsotopeBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-RenormIsotopeCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetReactorPowerBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-RenormReactorPowerCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetAbsoluteEnergyScaleBudget())
                                        {
                                            //                                        TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-RenormAbsoluteEnergyScaleCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetRelativeEnergyScaleBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-RenormRelativeEnergyScaleCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetAbsoluteEnergyOffsetBudget())
                                        {
                                            //                                        TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-RenormAbsoluteEnergyOffsetCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetRelativeEnergyOffsetBudget())
                                        {
                                            //                                        TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-RenormRelativeEnergyOffsetCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetIAVBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-RenormIAVCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetNLBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-RenormNLCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetResolutionBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-RenormResolutionCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetSin22t12Budget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-RenormSin22t12CovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetEfficiencyBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-RenormEfficiencyCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetBackgroundBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-BackgroundsCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetSystematicBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=TotalCovarianceMatrixM[x*MaxBins+y]-SystematicCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetTotalBudget())
                                        {
                                        }
                                    }
                                    else if(BudgetTurnOn)
                                    {
                                        //  Only one systematic at a time
                                        
                                        TotalCovarianceMatrixM[x*MaxBins+y]=CovStat[x][y];
                                        
                                        if(Data->GetVaryAccidentalBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=VAccCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetVaryLiHeBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=VLiHeCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetVaryFastNeutronsBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=VFNCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetVaryAmCBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=VAmCCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetDistortLiHeBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=DLiHeCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetDistortFastNeutronsBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=DFNCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetDistortAmCBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=DAmCCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetIsotopeBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=RenormIsotopeCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetReactorPowerBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=RenormReactorPowerCovarianceMatrixM[x][y];
                                        }
                                        //                                    else if(Data->GetAbsoluteEnergyScaleBudget())
                                        //                                    {
                                        //                                        //                                        TotalCovarianceMatrixM[x*MaxBins+y]+=RenormAbsoluteEnergyScaleCovarianceMatrixM[x][y];
                                        //                                    }
                                        else if(Data->GetRelativeEnergyScaleBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=RenormRelativeEnergyScaleCovarianceMatrixM[x][y];
                                        }
                                        //                                    else if(Data->GetAbsoluteEnergyOffsetBudget())
                                        //                                    {
                                        //                                        //                                        TotalCovarianceMatrixM[x*MaxBins+y]+=RenormAbsoluteEnergyOffsetCovarianceMatrixM[x][y];
                                        //                                    }
                                        //                                    else if(Data->GetRelativeEnergyOffsetBudget())
                                        //                                    {
                                        //                                        //                                        TotalCovarianceMatrixM[x*MaxBins+y]+=RenormRelativeEnergyOffsetCovarianceMatrixM[x][y];
                                        //                                    }
                                        else if(Data->GetIAVBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=RenormIAVCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetNLBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=RenormNLCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetResolutionBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=RenormResolutionCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetSin22t12Budget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=RenormSin22t12CovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetEfficiencyBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=RenormEfficiencyCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetBackgroundBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=BackgroundsCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetSystematicBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]+=SystematicCovarianceMatrixM[x][y];
                                        }
                                        else if(Data->GetTotalBudget())
                                        {
                                            TotalCovarianceMatrixM[x*MaxBins+y]=CovStat[x][y]+BackgroundsCovarianceMatrixM[x][y]+SystematicCovarianceMatrixM[x][y];
                                        }
                                    }
                                    
                                    //  Fill Histograms
                                    TotalCovarianceMatrixH->SetBinContent(x+1,y+1,TotalCovarianceMatrixM[x*MaxBins+y]);
                                    BackgroundsCovarianceMatrixH->SetBinContent(x+1,y+1,BackgroundsCovarianceMatrixM[x][y]);
                                    SystematicCovarianceMatrixH->SetBinContent(x+1,y+1,SystematicCovarianceMatrixM[x][y]);
                                    
                                    VAccCovarianceMatrixH->SetBinContent(x+1,y+1,VAccCovarianceMatrixM[x][y]);
                                    VLiHeCovarianceMatrixH->SetBinContent(x+1,y+1,VLiHeCovarianceMatrixM[x][y]);
                                    VFNCovarianceMatrixH->SetBinContent(x+1,y+1,VFNCovarianceMatrixM[x][y]);
                                    VAmCCovarianceMatrixH->SetBinContent(x+1,y+1,VAmCCovarianceMatrixM[x][y]);
                                    DLiHeCovarianceMatrixH->SetBinContent(x+1,y+1,DLiHeCovarianceMatrixM[x][y]);
                                    DFNCovarianceMatrixH->SetBinContent(x+1,y+1,DFNCovarianceMatrixM[x][y]);
                                    DAmCCovarianceMatrixH->SetBinContent(x+1,y+1,DAmCCovarianceMatrixM[x][y]);
                                    
                                    RenormIsotopeCovarianceMatrixH->SetBinContent(x+1,y+1,RenormIsotopeCovarianceMatrixM[x][y]);
                                    RenormReactorPowerCovarianceMatrixH->SetBinContent(x+1,y+1,RenormReactorPowerCovarianceMatrixM[x][y]);
                                    RenormRelativeEnergyScaleCovarianceMatrixH->SetBinContent(x+1,y+1,RenormRelativeEnergyScaleCovarianceMatrixM[x][y]);
                                    //                                RenormRelativeEnergyOffsetCovarianceMatrixH->SetBinContent(x+1,y+1,RenormRelativeEnergyOffsetCovarianceMatrixM[x][y]);
                                    //                                RenormAbsoluteEnergyScaleCovarianceMatrixH->SetBinContent(x+1,y+1,RenormAbsoluteEnergyScaleCovarianceMatrixM[x][y]);
                                    //                                RenormAbsoluteEnergyOffsetCovarianceMatrixH->SetBinContent(x+1,y+1,RenormAbsoluteEnergyOffsetCovarianceMatrixM[x][y]);
                                    RenormIAVCovarianceMatrixH->SetBinContent(x+1,y+1,RenormIAVCovarianceMatrixM[x][y]);
                                    RenormNLCovarianceMatrixH->SetBinContent(x+1,y+1,RenormNLCovarianceMatrixM[x][y]);
                                    RenormResolutionCovarianceMatrixH->SetBinContent(x+1,y+1,RenormResolutionCovarianceMatrixM[x][y]);
                                    RenormSin22t12CovarianceMatrixH->SetBinContent(x+1,y+1,RenormSin22t12CovarianceMatrixM[x][y]);
                                    RenormEfficiencyCovarianceMatrixH->SetBinContent(x+1,y+1,RenormEfficiencyCovarianceMatrixM[x][y]);
                                }
                                else
                                {
                                    RenormIsotopeCovarianceMatrixH->SetBinContent(x+1,y+1,IsotopeCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1)));
                                    
                                    RenormReactorPowerCovarianceMatrixH->SetBinContent(x+1,y+1,ReactorPowerCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1)));
                                    
                                    RenormRelativeEnergyScaleCovarianceMatrixH->SetBinContent(x+1,y+1,RelativeEnergyScaleCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1)));
                                    
                                    //                                RenormRelativeEnergyOffsetCovarianceMatrixH->SetBinContent(x+1,y+1,RelativeEnergyOffsetCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1)));
                                    
                                    //                                RenormAbsoluteEnergyScaleCovarianceMatrixH->SetBinContent(x+1,y+1,AbsoluteEnergyScaleCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1)));
                                    
                                    //                                RenormAbsoluteEnergyOffsetCovarianceMatrixH->SetBinContent(x+1,y+1,AbsoluteEnergyOffsetCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1)));
                                    
                                    RenormIAVCovarianceMatrixH->SetBinContent(x+1,y+1,IAVCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1)));
                                    
                                    RenormNLCovarianceMatrixH->SetBinContent(x+1,y+1,NLCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1)));
                                    
                                    RenormResolutionCovarianceMatrixH->SetBinContent(x+1,y+1,ResolutionCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1)));
                                    
                                    RenormSin22t12CovarianceMatrixH->SetBinContent(x+1,y+1,Sin22t12CovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1)));
                                    RenormEfficiencyCovarianceMatrixH->SetBinContent(x+1,y+1,EfficiencyCovarianceMatrixH->GetBinContent(x+1,y+1)*(PredictionH[fari][neari]->GetBinContent(i+1)*PredictionH[farj][nearj]->GetBinContent(j+1)));
                                    
                                }
                            }
                        }
                    }
                }
            }
        }
        
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                  Add all matrices into a single Total Covariance Matrix
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        if(!ReadTxt)
        {
            BackgroundsCovarianceMatrixH->Add(VAccCovarianceMatrixH);
            BackgroundsCovarianceMatrixH->Add(VLiHeCovarianceMatrixH);
            BackgroundsCovarianceMatrixH->Add(VFNCovarianceMatrixH);
            BackgroundsCovarianceMatrixH->Add(VAmCCovarianceMatrixH);
            BackgroundsCovarianceMatrixH->Add(DLiHeCovarianceMatrixH);
            BackgroundsCovarianceMatrixH->Add(DFNCovarianceMatrixH);
            BackgroundsCovarianceMatrixH->Add(DAmCCovarianceMatrixH);
            
            SystematicCovarianceMatrixH=(TH2D*)(RenormIAVCovarianceMatrixH->Clone("Systematic Covariance Matrix"));
            SystematicCovarianceMatrixH->Add(RenormIsotopeCovarianceMatrixH);
            SystematicCovarianceMatrixH->Add(RenormReactorPowerCovarianceMatrixH);
            SystematicCovarianceMatrixH->Add(RenormRelativeEnergyScaleCovarianceMatrixH);
            //        SystematicCovarianceMatrixH->Add(RenormRelativeEnergyOffsetCovarianceMatrixH);
            //        SystematicCovarianceMatrixH->Add(RenormAbsoluteEnergyScaleCovarianceMatrixH);
            //        SystematicCovarianceMatrixH->Add(RenormAbsoluteEnergyOffsetCovarianceMatrixH);
            SystematicCovarianceMatrixH->Add(RenormNLCovarianceMatrixH);
            SystematicCovarianceMatrixH->Add(RenormResolutionCovarianceMatrixH);
            SystematicCovarianceMatrixH->Add(RenormSin22t12CovarianceMatrixH);
            SystematicCovarianceMatrixH->Add(RenormEfficiencyCovarianceMatrixH);
            
            SystematicCovarianceMatrixH->SetTitle("Systematic Covariance Matrix");
            BackgroundsCovarianceMatrixH->SetTitle("Background Covariance Matrix");
            
            TotalCovarianceMatrixH->Add(StatisticalCovarianceMatrixH);
            TotalCovarianceMatrixH->Add(BackgroundsCovarianceMatrixH);
            TotalCovarianceMatrixH->Add(SystematicCovarianceMatrixH);
            
            // Specific combinations to produce the error budget:
            if(BudgetTurnOff)
            {
                //  Substract the choosen systematic
                
                if(Data->GetVaryAccidentalBudget())
                {
                    TotalCovarianceMatrixH->Add(VAccCovarianceMatrixH,-1);
                }
                else if(Data->GetVaryLiHeBudget())
                {
                    TotalCovarianceMatrixH->Add(VLiHeCovarianceMatrixH,-1);
                }
                else if(Data->GetVaryFastNeutronsBudget())
                {
                    TotalCovarianceMatrixH->Add(VFNCovarianceMatrixH,-1);
                }
                else if(Data->GetVaryAmCBudget())
                {
                    TotalCovarianceMatrixH->Add(VAmCCovarianceMatrixH,-1);
                }
                else if(Data->GetDistortLiHeBudget())
                {
                    TotalCovarianceMatrixH->Add(DLiHeCovarianceMatrixH,-1);
                }
                else if(Data->GetDistortFastNeutronsBudget())
                {
                    TotalCovarianceMatrixH->Add(DFNCovarianceMatrixH,-1);
                }
                else if(Data->GetDistortAmCBudget())
                {
                    TotalCovarianceMatrixH->Add(DAmCCovarianceMatrixH,-1);
                }
                else if(Data->GetIsotopeBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormIsotopeCovarianceMatrixH,-1);
                }
                else if(Data->GetReactorPowerBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormReactorPowerCovarianceMatrixH,-1);
                }
                else if(Data->GetAbsoluteEnergyScaleBudget())
                {
                    //                TotalCovarianceMatrixH->Add(RenormAbsoluteEnergyScaleCovarianceMatrixH,-1);
                }
                else if(Data->GetRelativeEnergyScaleBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormRelativeEnergyScaleCovarianceMatrixH,-1);
                }
                else if(Data->GetAbsoluteEnergyOffsetBudget())
                {
                    //                TotalCovarianceMatrixH->Add(RenormAbsoluteEnergyOffsetCovarianceMatrixH,-1);
                }
                else if(Data->GetRelativeEnergyOffsetBudget())
                {
                    //                TotalCovarianceMatrixH->Add(RenormRelativeEnergyOffsetCovarianceMatrixH,-1);
                }
                else if(Data->GetIAVBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormIAVCovarianceMatrixH,-1);
                }
                else if(Data->GetNLBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormNLCovarianceMatrixH,-1);
                }
                else if(Data->GetResolutionBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormResolutionCovarianceMatrixH,-1);
                }
                else if(Data->GetSin22t12Budget())
                {
                    TotalCovarianceMatrixH->Add(RenormSin22t12CovarianceMatrixH,-1);
                }
                else if(Data->GetEfficiencyBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormEfficiencyCovarianceMatrixH,-1);
                }
                else if(Data->GetBackgroundBudget())
                {
                    TotalCovarianceMatrixH->Add(BackgroundsCovarianceMatrixH,-1);
                }
                else if(Data->GetSystematicBudget())
                {
                    TotalCovarianceMatrixH->Add(SystematicCovarianceMatrixH,-1);
                }
                else if(Data->GetTotalBudget())
                {
                }
            }
            else if(BudgetTurnOn)
            {
                //  Only one systematic at a time
                TotalCovarianceMatrixH->Reset();
                TotalCovarianceMatrixH = (TH2D*)StatisticalCovarianceMatrixH->Clone();
                TotalCovarianceMatrixH->SetTitle("Error Turn-On Budget");
                
                if(Data->GetVaryAccidentalBudget())
                {
                    TotalCovarianceMatrixH->Add(VAccCovarianceMatrixH);
                }
                else if(Data->GetVaryLiHeBudget())
                {
                    TotalCovarianceMatrixH->Add(VLiHeCovarianceMatrixH);
                }
                else if(Data->GetVaryFastNeutronsBudget())
                {
                    TotalCovarianceMatrixH->Add(VFNCovarianceMatrixH);
                }
                else if(Data->GetVaryAmCBudget())
                {
                    TotalCovarianceMatrixH->Add(VAmCCovarianceMatrixH);
                }
                else if(Data->GetDistortLiHeBudget())
                {
                    TotalCovarianceMatrixH->Add(DLiHeCovarianceMatrixH);
                }
                else if(Data->GetDistortFastNeutronsBudget())
                {
                    TotalCovarianceMatrixH->Add(DFNCovarianceMatrixH);
                }
                else if(Data->GetDistortAmCBudget())
                {
                    TotalCovarianceMatrixH->Add(DAmCCovarianceMatrixH);
                }
                else if(Data->GetIsotopeBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormIsotopeCovarianceMatrixH);
                }
                else if(Data->GetReactorPowerBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormReactorPowerCovarianceMatrixH);
                }
                //            else if(Data->GetAbsoluteEnergyScaleBudget())
                //            {
                //                //                TotalCovarianceMatrixH->Add(RenormAbsoluteEnergyScaleCovarianceMatrixH);
                //            }
                else if(Data->GetRelativeEnergyScaleBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormRelativeEnergyScaleCovarianceMatrixH);
                }
                //            else if(Data->GetAbsoluteEnergyOffsetBudget())
                //            {
                //                //                TotalCovarianceMatrixH->Add(RenormAbsoluteEnergyOffsetCovarianceMatrixH);
                //            }
                //            else if(Data->GetRelativeEnergyOffsetBudget())
                //            {
                //                //                TotalCovarianceMatrixH->Add(RenormRelativeEnergyOffsetCovarianceMatrixH);
                //            }
                else if(Data->GetIAVBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormIAVCovarianceMatrixH);
                }
                else if(Data->GetNLBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormNLCovarianceMatrixH);
                }
                else if(Data->GetResolutionBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormResolutionCovarianceMatrixH);
                }
                else if(Data->GetSin22t12Budget())
                {
                    TotalCovarianceMatrixH->Add(RenormSin22t12CovarianceMatrixH);
                }
                else if(Data->GetEfficiencyBudget())
                {
                    TotalCovarianceMatrixH->Add(RenormEfficiencyCovarianceMatrixH);
                }
                else if(Data->GetBackgroundBudget())
                {
                    TotalCovarianceMatrixH->Add(BackgroundsCovarianceMatrixH);
                }
                else if(Data->GetSystematicBudget())
                {
                    TotalCovarianceMatrixH->Add(SystematicCovarianceMatrixH);
                }
                else if(Data->GetTotalBudget())
                {
                    TotalCovarianceMatrixH->Reset();
                    TotalCovarianceMatrixH->Add(StatisticalCovarianceMatrixH);
                    TotalCovarianceMatrixH->Add(BackgroundsCovarianceMatrixH);
                    TotalCovarianceMatrixH->Add(SystematicCovarianceMatrixH);
                }
            }
        }
        
        if(WriteROOT)
        {
            TFile* SaveSystematicF = TFile::Open(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/SystematicCovarianceMatrices.root",Combine)).c_str(),"recreate");//reno
            RenormIsotopeCovarianceMatrixH->Write();
            RenormReactorPowerCovarianceMatrixH->Write();
            RenormIAVCovarianceMatrixH->Write();
            RenormNLCovarianceMatrixH->Write();
            RenormResolutionCovarianceMatrixH->Write();
            RenormRelativeEnergyScaleCovarianceMatrixH->Write();
            //    RenormRelativeEnergyOffsetCovarianceMatrixH->Write();
            //    RenormAbsoluteEnergyScaleCovarianceMatrixH->Write();
            //    RenormAbsoluteEnergyOffsetCovarianceMatrixH->Write();
            RenormSin22t12CovarianceMatrixH->Write();
            RenormEfficiencyCovarianceMatrixH->Write();
            SystematicCovarianceMatrixH->Write();
            SaveSystematicF->Close();
            
            TFile* SaveBkg = TFile::Open(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/BackgroundCovarianceMatrices.root",Combine)).c_str(),"recreate");
            
            VAccCovarianceMatrixH->Write();
            VFNCovarianceMatrixH->Write();
            VLiHeCovarianceMatrixH->Write();
            VAmCCovarianceMatrixH->Write();
            DFNCovarianceMatrixH->Write();
            DLiHeCovarianceMatrixH->Write();
            DAmCCovarianceMatrixH->Write();
            BackgroundsCovarianceMatrixH->Write();
            
            SaveBkg->Close();
            
            TFile* SaveTotal= TFile::Open(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/TotalCovarianceMatrix.root",Combine)).c_str(),"recreate");
            
            TotalCovarianceMatrixH->Write("Total Covariance Matrix");
            
            SaveTotal->Close();
        }
    }

    std::cout <<  "\t Matrices Combined" << std::endl;
    
}

void Prediction :: ApplyStatisticalFluctuation(TH1D* Histo)
{
    for(Int_t VisibleEnergyIndex=1;VisibleEnergyIndex<=Histo->GetXaxis()->GetNbins();VisibleEnergyIndex++)
    {
        rand->SetSeed(0);
        Histo->SetBinContent(VisibleEnergyIndex,(Double_t)(rand->PoissonD(Histo->GetBinContent(VisibleEnergyIndex)*Histo->GetXaxis()->GetBinWidth(VisibleEnergyIndex))/Histo->GetXaxis()->GetBinWidth(VisibleEnergyIndex)));
    }
}

void Prediction :: InvertMatrix(Int_t week)
{
    TMatrixD* TotalCovarianceMatrix = new TMatrixD(MaxBins,MaxBins);
    
    if(!ReadTxt)
    {
        for(Int_t i = 0; i<MaxBins; i++)
        {
            for(Int_t j = 0; j<MaxBins; j++)
            {
                TotalCovarianceMatrixM[(j+MaxBins*i)]=TotalCovarianceMatrixH->GetBinContent(i+1,j+1);
            }
        }
    }
    
    TotalCovarianceMatrix->SetMatrixArray(&TotalCovarianceMatrixM[0]);
    
    TotalCovarianceMatrix->Invert();
    
    Double_t* InvTotalCovarianceMatrixArray = TotalCovarianceMatrix->GetMatrixArray();
    
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j < MaxBins; j++)
        {
            InvTotalCovarianceMatrixM[i][j]=InvTotalCovarianceMatrixArray[(i*MaxBins+j)];
            InvTotalCovarianceMatrixH->SetBinContent(i+1,j+1,InvTotalCovarianceMatrixM[i][j]);
            //            if(i==j)//test with diagonal matrix
            //            {
            //                InvTotalCovarianceMatrixH->SetBinContent(i+1,j+1,1);
            //            }
            //            else
            //            {
            //                InvTotalCovarianceMatrixH->SetBinContent(i+1,j+1,0);
            //            }
        }
    }
    
    UnityH=(TH2D*)TotalCovarianceMatrixH->Clone("Unity");
    UnityH->Reset();
    Double_t UnityM[MaxBins*MaxBins];
    
    for(Int_t i=0; i<MaxBins; i++)
    {
        for(Int_t j=0; j<MaxBins; j++)
        {
            UnityM[(j+MaxBins*i)] = 0;
            
            for(Int_t k=0; k<MaxBins; k++)
            {
                UnityM[(j+MaxBins*i)] += TotalCovarianceMatrixM[(k+MaxBins*i)] * InvTotalCovarianceMatrixM[k][j];
            }
        }
    }
    
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j < MaxBins; j++)
        {
            UnityH->SetBinContent(i+1,j+1,UnityM[(j+MaxBins*i)]);
        }
    }
    #ifdef PrintEps
        TCanvas* unityC = new TCanvas("unityC","unityC",1200,400);
        unityC->Divide(3,1);
        unityC->cd(1);
        TotalCovarianceMatrixH->SetStats(kFALSE);
        TotalCovarianceMatrixH->Draw("colz");
        unityC->cd(2);
        InvTotalCovarianceMatrixH->SetStats(kFALSE);
        InvTotalCovarianceMatrixH->Draw("colz");
        unityC->cd(3);
        UnityH->SetStats(kFALSE);
        UnityH->SetTitle("Unity Matrix");
        UnityH->Draw("colz");
        unityC->Update();
        unityC->Print("./Images/Test/TestInvert.eps", ".eps");
        delete unityC;
    #endif
    
    delete TotalCovarianceMatrix;
}

TH2D* Prediction :: NormCov(TH2D* Histo,TH2D* CopyHisto)
{
    CopyHisto = (TH2D*)Histo->Clone();
    
    for (Int_t i = 0; i<Histo->GetXaxis()->GetNbins(); i++)
    {
        for (Int_t j = 0; j<Histo->GetYaxis()->GetNbins(); j++)
        {
            CopyHisto->SetBinContent(i+1,j+1,Histo->GetBinContent(i+1,j+1)/(sqrt(Histo->GetBinContent(i+1,i+1)*Histo->GetBinContent(j+1,j+1))));
        }
    }
    
    return CopyHisto;
}

void Prediction :: ProduceCovToyMCSample(Int_t week,TH1D** NominalPredictionH)
{
    if (Combine == 1)
    {
        MaxNearCombine=1;
        MaxFarCombine=1;
        MaxBins=n_evis_bins;
    }
    else if(Combine == 2)
    {
        MaxNearCombine=2;
        MaxFarCombine=1;
        MaxBins=2*n_evis_bins;
    }
    else
    {
        MaxNearCombine = ADsEH1+ADsEH2;
        MaxFarCombine = ADsEH3;
        MaxBins=9*n_evis_bins;
    }
    
    TH2D* ToyMCSample[MaxSystematics];
    
    ToyMCSample[0] = (TH2D*)IsotopeCovarianceMatrixH->Clone();
    ToyMCSample[1] = (TH2D*)ReactorPowerCovarianceMatrixH->Clone();
    ToyMCSample[2] = (TH2D*)RelativeEnergyScaleCovarianceMatrixH->Clone();
    ToyMCSample[3] = (TH2D*)IAVCovarianceMatrixH->Clone();
    ToyMCSample[4] = (TH2D*)NLCovarianceMatrixH->Clone();
    ToyMCSample[5] = (TH2D*)ResolutionCovarianceMatrixH->Clone();
    ToyMCSample[6] = (TH2D*)Sin22t12CovarianceMatrixH->Clone();
    ToyMCSample[7] = (TH2D*)EfficiencyCovarianceMatrixH->Clone();
    
    ToyMCSample[8]  = (TH2D*)IsotopeCovarianceMatrixH->Clone();//Total systematics all together
    
    for(Int_t SystematicI = 1;SystematicI<MaxSystematics-1; SystematicI++)
    {
        ToyMCSample[8]->Add(ToyMCSample[SystematicI]);
    }
    
    Int_t x =0;
    Int_t y =0;
    Int_t Ni1=0,Ni2=0,Ni3=0,Ni4=0;
    Int_t Nj1=0,Nj2=0,Nj3=0,Nj4=0;
    Int_t Fi1=0,Fi2=0,Fi3=0,Fi4=0;
    Int_t Fj1=0,Fj2=0,Fj3=0,Fj4=0;
    
    for(Int_t SystematicI = 8;SystematicI<MaxSystematics; SystematicI++)
    {
        
        VariationHistoH[SystematicI] = new TH1D(Form("RandomPredictionFromSystematic%d",SystematicI),Form("RandomPredictionFromSystematic%d",SystematicI),MaxBins,0,MaxBins);
        
        for (Int_t i = 0; i<MaxBins; i++)
        {
            if(i<n_evis_bins)
            {
                VariationHistoH[SystematicI]->SetBinContent(i+1, NominalPredictionH[0+MaxFarCombine*0]->GetBinContent(i+1));
            }
            else
            {
                VariationHistoH[SystematicI]->SetBinContent(i+1, NominalPredictionH[0+MaxFarCombine*1]->GetBinContent(i+1));
            }
        }
        
        for (Int_t neari=0; neari<MaxNearCombine; neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Ni1=1;Ni2=0;Ni3=0;Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<MaxFarCombine; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Fi1=fari+1;Fi2=0;Fi3=0;Fi4=0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFarCombine; Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFarCombine; Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=MaxFarCombine; Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<MaxNearCombine; nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<MaxFarCombine; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Fj1=farj+1;Fj2=0;Fj3=0;Fj4=0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFarCombine; Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFarCombine; Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=MaxFarCombine; Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                            
                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                
                                RenormToyMCSample[x][y]=(ToyMCSample[SystematicI]->GetBinContent(x+1,y+1)*(NominalPredictionH[fari+MaxFarCombine*neari]->GetBinContent(i+1)*NominalPredictionH[farj+MaxFarCombine*nearj]->GetBinContent(j+1)));
                            }
                        }
                    }
                }
            }
        }
        
        #ifdef PrintEps//just a check
            TH2D* PrintCov2H = new TH2D("Cov2H","Cov2H",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
            TH2D* CorrelationH = new TH2D("Cov2H","Cov2H",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
            
            for (Int_t i = 0; i<MaxBins; i++)
            {
                for (Int_t j = 0; j<MaxBins; j++)
                {
                    PrintCov2H->SetBinContent(i+1,j+1,RenormToyMCSample[i][j]);
                    Correlation[i][j] = RenormToyMCSample[i][j]/(sqrt(RenormToyMCSample[i][i]*RenormToyMCSample[j][j]));
                    CorrelationH->SetBinContent(i+1,j+1,Correlation[i][j]);
                    
                }
            }
            
            TCanvas* CovarianceMatrixC = new TCanvas("","");
            
            CovarianceMatrixC->Divide(2,2);
            
            CovarianceMatrixC->cd(1);
            PrintCov2H->Draw("colz");
            CovarianceMatrixC->cd(2);
            CorrelationH->Draw("colz");
            
            CovarianceMatrixC->cd(3);
            NominalPredictionH[0]->Draw();
            
            if(Combine==2)
            {
                CovarianceMatrixC->cd(4);
                NominalPredictionH[1]->Draw();
            }
            
            CovarianceMatrixC->Print(("./Images/"+AnalysisString+Form("/RandomCovarianceMatrix%d.eps",SystematicI)).c_str());
            
            delete CovarianceMatrixC;
        #endif
        
        TMatrixD renormmatrix(9*MaxNbins,9*MaxNbins, &RenormToyMCSample[0][0]);
        renormmatrix.ResizeTo(MaxBins,MaxBins);
        
//        for(Int_t i = 0; i<MaxBins;i++)
//        {
//            for(Int_t j = 0; j<MaxBins;j++)
//            {
//                if(renormmatrix[i][j]!=renormmatrix[j][i])
//                {
//                    std::cout<<"It isn't a symmetric matrix" << std::endl;
//                }
//                else
//                {
//                    std::cout << renormmatrix[i][j] << " = " << renormmatrix[i][j]  << std::endl;
//                }
//            }
//        }
        TDecompChol cholrenorm(renormmatrix);//  M = L*U
        cholrenorm.Decompose();
        TMatrixD cmat(cholrenorm.GetU());//   U
        TMatrixD tcmat(cmat.Transpose(cmat));// L
        
        Double_t* tmp_matrix = tcmat.GetMatrixArray();
        
        #ifdef PrintEps
            TCanvas* MatrixC = new TCanvas("","");
            
            MatrixC->Divide(3,1);
            
            MatrixC->cd(1);
            renormmatrix.Draw("colz");
            MatrixC->cd(2);
            cmat.Draw("colz");
            MatrixC->cd(3);
            tcmat.Draw("colz");
            
            MatrixC->Print(("./Images/"+AnalysisString+Form("/MatrixRenormCheck%d.eps",SystematicI)).c_str());
            
            delete MatrixC;
        #endif
        
        for (Int_t i = 0; i < MaxBins; i++)
        {
            for (Int_t j = 0; j < MaxBins; j++)
            {
                L[SystematicI][i][j] = tmp_matrix[i*MaxBins+j];
            }
        }
        
        renormmatrix.~TMatrixD();
        cmat.~TMatrixD();
        tcmat.~TMatrixD();
        delete ToyMCSample[SystematicI];
    }
}

TH1D* Prediction :: GetToyMCSample(Int_t Systematic)
{
    flag_delete_ToyMCSample = 1;
    
    Double_t ranvec[MaxBins];
    
    for (Int_t i = 0; i < MaxBins; i++)
    {
        rand->SetSeed(0);
        ranvec[i] = rand->Gaus(0,1);
    }
    
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j < MaxBins; j++)
        {
            VariationHistoH[Systematic]->SetBinContent(i+1, VariationHistoH[Systematic]->GetBinContent(i+1)+L[Systematic][i][j] * ranvec[j]);
        }
    }
    
    return VariationHistoH[Systematic];
}

Double_t Prediction :: ChiSquare(Int_t week)
{
    Double_t chi2=0;
    
    Int_t x =0;
    Int_t y =0;
    
    Int_t Ni1=0,Ni2=0,Ni3=0,Ni4=0;
    Int_t Nj1=0,Nj2=0,Nj3=0,Nj4=0;
    Int_t Fi1=0,Fi2=0,Fi3=0,Fi4=0;
    Int_t Fj1=0,Fj2=0,Fj3=0,Fj4=0;
    
    for (Int_t neari=0; neari<MaxNearCombine; neari++)
    {
        //Logic for the 2D matrix index done up to 8 ADs
        if(neari==0){Ni1=1;Ni2=0;Ni3=0;Ni4=0;}
        if(neari==1){Ni2++;}
        if(neari==2){Ni3++;}
        if(neari==3){Ni4++;}
        
        for (Int_t fari=0; fari<MaxFarCombine; fari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(Ni1!=Ni2){Fi1=fari+1;Fi2=0;Fi3=0;Fi4=0;}
            if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFarCombine; Fi2=fari+1;}
            if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFarCombine; Fi3=fari+1;}
            if(Ni3==Ni4&&Ni4==1){Fi3=MaxFarCombine; Fi4=fari+1;}
            
            for (Int_t nearj=0; nearj<MaxNearCombine; nearj++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                if(nearj==1){Nj2++;}
                if(nearj==2){Nj3++;}
                if(nearj==3){Nj4++;}
                
                for (Int_t farj=0; farj<MaxFarCombine; farj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(Nj1!=Nj2){Fj1=farj+1;Fj2=0;Fj3=0;Fj4=0;}
                    if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFarCombine; Fj2=farj+1;}
                    if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFarCombine; Fj3=farj+1;}
                    if(Nj3==Nj4&&Nj4==1){Fj3=MaxFarCombine; Fj4=farj+1;}
                    
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {//columns
                        x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                        
                        for (Int_t j = 0; j<n_evis_bins; j++)
                        {//rows
                            y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                            
                            if(ReadTxt)
                            {
                                chi2+= ((FarDataH[fari][neari]->GetBinContent(i+1)-PredictionDataH[fari][neari]->GetBinContent(i+1))*InvTotalCovarianceMatrixM[x][y]*(FarDataH[farj][nearj]->GetBinContent(j+1)-PredictionDataH[farj][nearj]->GetBinContent(j+1)));
                            }
                            else
                            {
                                //
                                
                                chi2+= ((FarDataH[fari][neari]->GetBinContent(i+1)-PredictionDataH[fari][neari]->GetBinContent(i+1))*InvTotalCovarianceMatrixH->GetBinContent(x+1,y+1)*(FarDataH[farj][nearj]->GetBinContent(j+1)-PredictionDataH[farj][nearj]->GetBinContent(j+1)));
                                
                                //                                std::cout << "x : " << x << "; y : " << y << "; Difference i : " << (FarDataH[fari]->GetBinContent(i+1)-PredictionDataH[fari][neari]->GetBinContent(i+1)) << "; Inv x,y : " << InvTotalCovarianceMatrixH->GetBinContent(x+1,y+1) << "; Difference j : " << FarDataH[farj]->GetBinContent(j+1)-PredictionDataH[farj][nearj]->GetBinContent(j+1) << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
    
    std::cout << " Total χ2: " << chi2 << std::endl;
    
    return chi2;
}

Double_t Prediction :: RateChiSquare(Int_t week)
{
    Double_t chi2=0;
    
    TH1D* RateFarDataH[MaxFarCombine][MaxNearCombine];
    TH1D* RatePredictionDataH[MaxFarCombine][MaxNearCombine];
    
    for (Int_t fari=0; fari<MaxFarCombine; fari++)
    {
        for (Int_t neari=0; neari<MaxNearCombine; neari++)
        {
            RateFarDataH[fari][neari] = (TH1D*)FarDataH[fari][neari]->Rebin(n_evis_bins,Form("RateFarDataH%d",fari));
            
            RatePredictionDataH[fari][neari] = (TH1D*)PredictionDataH[fari][neari]->Rebin(n_evis_bins,Form("RatePredictionDataH%d%d",fari,neari));
        }
    }
    
    TMatrixD* RateTotalCovarianceMatrix = new TMatrixD(MaxBins,MaxBins);
    
    Double_t RateTotalCovarianceMatrixM[MaxBins*MaxBins];
    
    for(Int_t i = 0; i<MaxBins; i++)
    {
        for(Int_t j = 0; j<MaxBins; j++)
        {
            RateTotalCovarianceMatrixM[(j+MaxBins*i)]=StatisticalCovarianceMatrixH->GetBinContent(i+1,j+1);
        }
    }
    TH2D* RateInvTotalCovarianceMatrixH = (TH2D*)StatisticalCovarianceMatrixH->Clone();
    RateInvTotalCovarianceMatrixH->Reset();
    
    RateTotalCovarianceMatrix->SetMatrixArray(&RateTotalCovarianceMatrixM[0]);
    
    RateTotalCovarianceMatrix->Invert();
    
    Double_t* RateInvTotalCovarianceMatrixArray = RateTotalCovarianceMatrix->GetMatrixArray();
    
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j < MaxBins; j++)
        {
            RateInvTotalCovarianceMatrixH->SetBinContent(i+1,j+1,RateInvTotalCovarianceMatrixArray[(i*MaxBins+j)]);
        }
    }
    
    RateInvTotalCovarianceMatrixH->Rebin2D(n_evis_bins,n_evis_bins);
    
    #ifdef PrintEps
        TCanvas* TestRateC = new TCanvas("TestRateChi2","TestRateChi2");
        TestRateC->Divide(3,1);
        TestRateC->cd(1);
        
        RateFarDataH[0][0]->Draw();
        
        TestRateC->cd(2);
        RatePredictionDataH[0][0]->Draw();
        
        TestRateC->cd(3);
        RateInvTotalCovarianceMatrixH->Draw("colz");
        
        TestRateC->Print("./Images/TestRateChi2.eps", ".eps");
        delete TestRateC;
    #endif
    
    Int_t x =0;
    Int_t y =0;
    
    Int_t Ni1=0,Ni2=0,Ni3=0,Ni4=0;
    Int_t Nj1=0,Nj2=0,Nj3=0,Nj4=0;
    Int_t Fi1=0,Fi2=0,Fi3=0,Fi4=0;
    Int_t Fj1=0,Fj2=0,Fj3=0,Fj4=0;
    
    for (Int_t neari=0; neari<MaxNearCombine; neari++)
    {
        //Logic for the 2D matrix index done up to 8 ADs
        if(neari==0){Ni1=1;Ni2=0;Ni3=0;Ni4=0;}
        if(neari==1){Ni2++;}
        if(neari==2){Ni3++;}
        if(neari==3){Ni4++;}
        
        for (Int_t fari=0; fari<MaxFarCombine; fari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(Ni1!=Ni2){Fi1=fari+1;Fi2=0;Fi3=0;Fi4=0;}
            if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFarCombine; Fi2=fari+1;}
            if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFarCombine; Fi3=fari+1;}
            if(Ni3==Ni4&&Ni4==1){Fi3=MaxFarCombine; Fi4=fari+1;}
            
            for (Int_t nearj=0; nearj<MaxNearCombine; nearj++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                if(nearj==1){Nj2++;}
                if(nearj==2){Nj3++;}
                if(nearj==3){Nj4++;}
                
                for (Int_t farj=0; farj<MaxFarCombine; farj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(Nj1!=Nj2){Fj1=farj+1;Fj2=0;Fj3=0;Fj4=0;}
                    if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFarCombine; Fj2=farj+1;}
                    if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFarCombine; Fj3=farj+1;}
                    if(Nj3==Nj4&&Nj4==1){Fj3=MaxFarCombine; Fj4=farj+1;}
                    
                    for (Int_t i = 0; i<1; i++)
                    {//columns
                        x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*1;
                        
                        for (Int_t j = 0; j<1; j++)
                        {//rows
                            y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*1;
                            
                            chi2+= ((RateFarDataH[fari][neari]->GetBinContent(i+1)-RatePredictionDataH[fari][neari]->GetBinContent(i+1))*RateInvTotalCovarianceMatrixH->GetBinContent(x+1,y+1)*(RateFarDataH[farj][nearj]->GetBinContent(j+1)-RatePredictionDataH[farj][nearj]->GetBinContent(j+1)));
                            
                            std::cout << "x : " << x << "; y : " << y << "; Rate far data i : " << RateFarDataH[fari][neari]->GetBinContent(i+1) << "; Rate Prediction i : " << RatePredictionDataH[fari][neari]->GetBinContent(i+1) << "; Inv x,y : " << RateInvTotalCovarianceMatrixH->GetBinContent(x+1,y+1) << "; Rate far data j : " << RateFarDataH[farj][nearj]->GetBinContent(j+1) << "; Rate Prediction j : " << RatePredictionDataH[farj][nearj]->GetBinContent(j+1) << std::endl;
                        }
                    }
                }
            }
        }
    }
    
    std::cout << " Total χ2: " << chi2 << std::endl;
    
    return chi2;
}

void Prediction :: DeleteData()
{
    for (Int_t near=0; near<(ADsEH1+ADsEH2); near++)
    {
        delete CombinedNearDataH[near];
        
        for (Int_t far=0; far<ADsEH3; far++)
        {
            delete FarDataH[far][near];
        }
    }
}

void Prediction :: DeleteMatrices()
{
    delete VAccCovarianceMatrixH;
    delete VLiHeCovarianceMatrixH;
    delete VAmCCovarianceMatrixH;
    delete DFNCovarianceMatrixH;
    delete DLiHeCovarianceMatrixH;
    delete DAmCCovarianceMatrixH;
    delete IsotopeCovarianceMatrixH;
    delete ReactorPowerCovarianceMatrixH;
    delete RelativeEnergyScaleCovarianceMatrixH;
    delete IAVCovarianceMatrixH;
    delete NLCovarianceMatrixH;
    delete ResolutionCovarianceMatrixH;
    delete Sin22t12CovarianceMatrixH;
    delete EfficiencyCovarianceMatrixH;
    
    delete InvTotalCovarianceMatrixH;
}
void Prediction :: SetSin22t13(Double_t st13)
{
    Sin22t13 = st13;
}
void Prediction :: SetDM213(Double_t dm13)
{
    DM2_ee = dm13;
}

NominalData* Prediction :: GetNominalData()
{
    return Data;
}

void Prediction :: SetNominalData(NominalData* data)
{
    Data = data;
} 

void Prediction :: GenerateFluxCorrectedHistograms(NominalData* LocalData)
{
    Oscillation* LocalOsc = new Oscillation(LocalData);
  
    LocalOsc->LoadFluxHisto();

    for(Int_t week =0;week<Nweeks;week++)
    {
        for(Int_t AD =0;AD<NADs;AD++)
        {
            for(Int_t reactor =0;reactor<NReactors;reactor++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        TH1D* CopyFluxHisto = LocalOsc->GetFluxHisto(reactor,AD,week,idx,idy);//NReactorPeriods and periods to be fit respectively, first is used to correct the flux for weekly efficiencies, the second to fit one period of data or more.
                
                        FluxHisto[reactor][AD][week][idx][idy] = (TH1D*)CopyFluxHisto->Clone();
                    }
                }
            }
        }
    }
    delete LocalOsc;
}

void Prediction :: SetExperiment(Int_t experiment)
{
    Experiment = experiment;
}
