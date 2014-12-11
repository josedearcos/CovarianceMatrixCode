//In this class I store the nominal data, standard values for each parameter.
#pragma once
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include <sstream>
#include <string>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TH2.h"
#include "TBenchmark.h"

bool TestAllTheSame = 0;

const bool Print = 1;//To save results in .eps files

const bool DeltaMee = 0;//Use Δm^2ee instead of Δm32 and Δm31 values

const bool ADSimple = 1;

const bool LBNLBinning = 1;
const bool LoganBinning = 0;

const Int_t MaxExperiments = 1000;
const Int_t MaxPeriods = 32;
const Int_t MaxDetectors = 8;
const Int_t NIsotopes = 4;
const Int_t NHalls = 3;
const Int_t NReactors=6;
const Int_t MaxNearDetectors =4;
const Int_t MaxFarDetectors =4;
//const Int_t MaxNbins=51;
const Int_t MaxNbins=240;//temporally to match logan's binning

const Int_t MatrixBins = 240;
const Int_t MaxNearADs =4;
const Int_t MaxFarADs =4;
const Int_t Nbins=240-36;

// To control how to produce the predictions with variations
const bool Fake_Experiments = 0;
const bool CholeskyVariations = 0;

class NominalData
{
private:
    
    std::string ResponseDirectory;
    std::string PredictionsDirectory;
    std::string ToyMCSamplesDirectory;
    std::string SysCovDirectory;
    std::string BkgCovDirectory;
    
    //Analysis type:
    bool isH;
    Int_t DataSet;
    bool ToyMC;
    //Binning type:
    bool LinearBinning;//0 is non linear (LBNL like binning), 1 is linear (51 bins);
    
    //Detectors:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    
    Int_t hierarchy;
    Int_t Combine;
    
    Double_t ProtonsPerKton; //protons per kton
    Double_t DetectorMassGdLs[MaxDetectors];
    Double_t DetectorMassLs[MaxDetectors];
    
    Double_t m_detectorEfficiency_Dt;
    Double_t m_detectorEfficiency_Ep;
    Double_t m_detectorEfficiency_Ed_nominal;
    Double_t m_detectorEfficiency_flash;
    Double_t m_detectorEfficiency_nGd;
    Double_t m_detectorEfficiency_spill;
    
    //Toy samples
    Int_t NSamples;
    Int_t NSteps;
    //Oscillation parameters
    Double_t s22t12;
    Double_t s22t12Error;
    Double_t s22t13;
    Double_t s22t13Error;
    Double_t dm2_ee;
    Double_t dm2_31;
    Double_t dm2_32;
    Double_t dm2_21;
    Double_t sin_start;
    Double_t sin_end;
    Double_t dmee_start;
    Double_t dmee_end;
    Int_t NReactorPeriods;

    //Reactor parameters
    Double_t IsotopeFrac[NIsotopes];
    Double_t IsotopeFracError[NIsotopes];
    
    Double_t ReactorPower[NReactors];
    Double_t ReactorPowerError[NReactors];
    
    Double_t EnergyPerFission[NReactors];
    Double_t EnergyPerFissionError[NReactors];
    
    //Background relative errors, rates and events:
    Double_t AccidentalError[MaxDetectors*MaxPeriods];
    Double_t LiHeError[MaxDetectors*MaxPeriods];
    Double_t FastNeutronError[MaxDetectors*MaxPeriods];
    Double_t AmCError[MaxDetectors*MaxPeriods];
    
    Double_t AccidentalRate[MaxDetectors*MaxPeriods];
    Double_t FastNeutronRate[MaxDetectors*MaxPeriods];
    Double_t LiHeRate[MaxDetectors*MaxPeriods];
    Double_t AmCRate[MaxDetectors*MaxPeriods];
    
    Double_t AccidentalEvents[MaxDetectors*MaxPeriods];
    Double_t FastNeutronEvents[MaxDetectors*MaxPeriods];
    Double_t LiHeEvents[MaxDetectors*MaxPeriods];
    Double_t AmCEvents[MaxDetectors*MaxPeriods];
    
    Double_t ObservedEvents[MaxDetectors*MaxPeriods];
    
    // Days and efficiencies:
    Double_t FullTime[MaxDetectors*MaxPeriods];
    Double_t MuonEff[MaxDetectors*MaxPeriods];
    Double_t MultiEff[MaxDetectors*MaxPeriods];
    
    //IAV error:
    Double_t IAVError;
    
    //Resolution errors:
    Double_t ResolutionError;
    Double_t ResolutionErrorUncorrelated;
    
    //Energy scale:
    Double_t m_abs_escale;
    Double_t m_abs_escale_error;
    
    Double_t m_abs_eoffset;
    Double_t m_abs_eoffset_error;
    
    Double_t m_rel_eoffset_error;
    Double_t m_rel_escale[MaxDetectors];
    Double_t m_rel_escale_error[MaxDetectors];
    Double_t m_rel_escale_nominal[MaxDetectors];
    Double_t m_rel_eoffset[MaxDetectors];
    
    Double_t DetectorEfficiencyRelativeError;
    
    //Binning variables:
    Int_t Nweeks;
    Int_t NBins;
    Double_t Emin;
    Double_t Emax;
    Double_t EVisMin;
    Double_t EVisMax;
    
    //Reactor
    Double_t binWidth;
    Double_t m_dNdE_nom[NReactors*MatrixBins];
    Double_t m_dNdE_mcov[NReactors*MatrixBins][NReactors*MatrixBins];
    Int_t m_nSamples;
    Double_t m_eMin;
    Double_t m_eMax;
    
    //Reactor Covariance Matrix
    Double_t L[NReactors*MatrixBins*NReactors*MatrixBins]; // lower triangle of the reactor covariance matrix

public:
    
    //binning variables
    Int_t n_evis_bins;
    Int_t n_etrue_bins;
    
    Double_t evis_bins[MatrixBins+1];//+1 because these arrays contain binning limits
    Double_t enu_bins[MatrixBins+1];
    
    bool TurnOnBudget;
    bool TurnOffBudget;
    
    bool BCW;
    bool LBNL;
    bool Unified;
    
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
    bool IAVMatrix;
    bool NLMatrix;
    bool ResolutionMatrix;
    bool Sin22t12Matrix;
    bool EfficiencyMatrix;
    
    bool VaryAccidentalBudget;
    bool VaryLiHeBudget;
    bool VaryFastNeutronsBudget;
    bool VaryAmCBudget;
    bool DistortLiHeBudget;
    bool DistortFastNeutronsBudget;
    bool DistortAmCBudget;
    bool IsotopeBudget;
    bool ReactorPowerBudget;
    bool RelativeEnergyOffsetBudget;
    bool AbsoluteEnergyOffsetBudget;
    bool AbsoluteEnergyScaleBudget;
    bool RelativeEnergyScaleBudget;
    bool IAVBudget;
    bool NLBudget;
    bool ResolutionBudget;
    bool Sin22t12Budget;
    bool EfficiencyBudget;
    bool SystematicBudget;
    bool BackgroundBudget;
    bool TotalBudget;
    
    bool StatisticalFluctuation;
    bool UseToyMCTree;
    
    NominalData(bool,Int_t);
    void CopyData(NominalData*);

    void ReadChristineCovMatrix();
    void ReadChristineReactorSpectrum();
    
    void SetDataSet(Int_t);
    void SetCombineMode(Int_t);
    void SetWeeks(Int_t);
    void SetSin22t12(Double_t);
    void SetSin22t13(Double_t);
    void SetDm2ee(Double_t);
    void SetDm231(Double_t);
    void SetDm232(Double_t);
    void SetDm221(Double_t);
    void SetSinStart(Double_t);
    void SetSinEnd(Double_t);
    void SetDmeeStart(Double_t);
    void SetDmeeEnd(Double_t);
    
    void SetNSamples(Int_t);
    void SetNSteps(Int_t);
    void SetEmin(Double_t);
    void SetEmax(Double_t);
    void SetEVisMin(Double_t);
    void SetEVisMax(Double_t);
    void SetADs(Int_t);
    void SetAnalysis(bool);
    void SetToyMC(bool);
    void SetStatisticalFluctuation(bool);
    void SetUseToyMCTree(bool);
    void SetBinning(bool);
    
    void SetAllRandomSystematics(bool);
    
    void LoadMainData(const Char_t*);//To correct for efficiencies using Gd/H data info from a txt file.
    void SetResponseDirectory(std::string);
    void SetToyMCSamplesDirectory(std::string);
    void SetPredictionDirectory(std::string);
    void SetBkgCovDirectory(std::string);
    void SetSysCovDirectory(std::string);
    
    void SetHierarchy(Int_t);
    
    //setters:
    void SetBCWModel(bool);
    void SetLBNLModel(bool);
    void SetUnifiedModel(bool);
    
    void SetVaryAccidentalMatrix(bool);
    void SetVaryLiHeMatrix(bool);
    void SetVaryFastNeutronsMatrix(bool);
    void SetVaryAmCMatrix(bool);
    void SetDistortLiHeMatrix(bool);
    void SetDistortFastNeutronsMatrix(bool);
    void SetDistortAmCMatrix(bool);
    
    void SetIsotopeMatrix(bool);
    void SetReactorPowerMatrix(bool);
    void SetRelativeEnergyScaleMatrix(bool);
    void SetRelativeEnergyOffsetMatrix(bool);
    void SetAbsoluteEnergyScaleMatrix(bool);
    void SetAbsoluteEnergyOffsetMatrix(bool);
    void SetIAVMatrix(bool);
    void SetNLMatrix(bool);
    void SetResolutionMatrix(bool);
    void SetSin22t12Matrix(bool);
    void SetEfficiencyMatrix(bool);
    void SetVaryAccidentalBudget(bool);
    void SetVaryLiHeBudget(bool);
    void SetVaryFastNeutronsBudget(bool);
    void SetVaryAmCBudget(bool);
    void SetDistortLiHeBudget(bool);
    void SetDistortFastNeutronsBudget(bool);
    void SetDistortAmCBudget(bool);
    
    void SetIsotopeBudget(bool);
    void SetReactorPowerBudget(bool);
    void SetRelativeEnergyScaleBudget(bool);
    void SetRelativeEnergyOffsetBudget(bool);
    void SetAbsoluteEnergyScaleBudget(bool);
    void SetAbsoluteEnergyOffsetBudget(bool);
    void SetIAVBudget(bool);
    void SetNLBudget(bool);
    void SetResolutionBudget(bool);
    void SetSin22t12Budget(bool);
    void SetEfficiencyBudget(bool);
    void SetBackgroundBudget(bool);
    void SetSystematicBudget(bool);
    void SetTotalBudget(bool);
    
    void SetTurnOnBudget(bool);
    void SetTurnOffBudget(bool);
    
    void SetNReactorPeriods(Int_t);
    
    //getters
    Double_t GetTrueBinningArray(Int_t);
    Double_t GetVisibleBinningArray(Int_t);
    Int_t GetTrueBins();
    Int_t GetVisibleBins();

    Int_t GetDataSet();

    bool GetBCWModel();
    bool GetLBNLModel();
    bool GetUnifiedModel();
    
    bool GetTurnOnBudget();
    bool GetTurnOffBudget();
    
    bool GetVaryAccidentalMatrix();
    bool GetVaryLiHeMatrix();
    bool GetVaryFastNeutronsMatrix();
    bool GetVaryAmCMatrix();
    bool GetDistortLiHeMatrix();
    bool GetDistortFastNeutronsMatrix();
    bool GetDistortAmCMatrix();
    
    bool GetIsotopeMatrix();
    bool GetReactorPowerMatrix();
    bool GetRelativeEnergyScaleMatrix();
    bool GetAbsoluteEnergyOffsetMatrix();
    bool GetRelativeEnergyOffsetMatrix();
    bool GetAbsoluteEnergyScaleMatrix();
    bool GetIAVMatrix();
    bool GetNLMatrix();
    bool GetResolutionMatrix();
    bool GetSin22t12Matrix();
    bool GetEfficiencyMatrix();
    
    bool GetVaryAccidentalBudget();
    bool GetVaryLiHeBudget();
    bool GetVaryFastNeutronsBudget();
    bool GetVaryAmCBudget();
    bool GetDistortLiHeBudget();
    bool GetDistortFastNeutronsBudget();
    bool GetDistortAmCBudget();
    
    bool GetIsotopeBudget();
    bool GetReactorPowerBudget();
    bool GetAbsoluteEnergyScaleBudget();
    bool GetRelativeEnergyScaleBudget();
    bool GetAbsoluteEnergyOffsetBudget();
    bool GetRelativeEnergyOffsetBudget();
    bool GetIAVBudget();
    bool GetNLBudget();
    bool GetResolutionBudget();
    bool GetSin22t12Budget();
    bool GetEfficiencyBudget();
    bool GetBackgroundBudget();
    bool GetSystematicBudget();
    bool GetTotalBudget();

    Double_t GetReactorCovMatrix(Int_t,Int_t);
    Double_t GetNominalReactorSpectrum(Int_t,Int_t);
    Int_t GetReactorSamples();
    Double_t GetReactorBinWidth();
    Double_t GetReactorEmin();
    Double_t GetReactorEmax();
    
    std::string GetResponseDirectory();
    std::string GetToyMCSamplesDirectory();
    std::string GetPredictionDirectory();
    std::string GetBkgCovDirectory();
    std::string GetSysCovDirectory();
    
    Int_t GetHierarchy();
    Int_t GetNSamples();
    Int_t GetNSteps();
    bool GetToyMC();
    bool GetStatisticalFluctuation();
    bool GetUseToyMCTree();
    
    Double_t GetSin22t13();
    
    Double_t GetSinStart();
    Double_t GetSinEnd();
    Double_t GetDmeeStart();
    Double_t GetDmeeEnd();
    
    Double_t GetDm2ee();
    Double_t GetDm231();
    Double_t GetDm232();
    Double_t GetDm221();
    Double_t GetSin22t12();
    Double_t GetSin22t12Error();
    
    Double_t GetIsotopeFraction(Int_t);
    Double_t GetReactorPower(Int_t);
    Double_t GetEnergyPerFission(Int_t);
    
    Double_t GetIsotopeFractionError(Int_t);
    Double_t GetReactorPowerError(Int_t);
    Double_t GetEnergyPerFissionError(Int_t);
    Double_t GetAbsoluteEnergyScale();
    Double_t GetAbsoluteEnergyScaleError();
    Double_t GetAbsoluteEnergyOffset();
    Double_t GetAbsoluteEnergyOffsetError();
    
    Double_t GetRelativeEnergyScale(Int_t);
    Double_t GetRelativeEnergyError(Int_t);
    Double_t GetRelativeEnergyOffset(Int_t);
    Double_t GetRelativeEnergyOffsetError();
    
    Double_t GetResolutionError();
    Double_t GetResoUncorrelatedError();
    
    Double_t GetDetectorProtonsGdLs(Int_t);
    Double_t GetDetectorProtonsLs(Int_t);
    
    Double_t GetEnergyDelayedCutDetectorEfficiency();
    Double_t GetDetectorEfficiencyRelativeError();
    Double_t GetDetectorEfficiency(Int_t, Int_t);
    Double_t GetFullTime(Int_t,Int_t);
    Int_t GetCombineMode();
    Int_t GetWeeks();
    Double_t GetEmin();
    Double_t GetEmax();
    Double_t GetEVisMin();
    Double_t GetEVisMax();
    Int_t GetADs();
    bool GetAnalysis();
    bool GetBinning();
    void CalculateBinning();
    //IAV thickness error
    Double_t GetIAVError();
    // Background relative errors
    Double_t GetAccidentalError(Int_t,Int_t);
    Double_t GetLiHeError(Int_t,Int_t);
    Double_t GetFNError(Int_t,Int_t);
    Double_t GetAmCError(Int_t,Int_t);
    
    Double_t GetIBDEvents(Int_t, Int_t);
    Double_t GetObservedEvents(Int_t,Int_t);

    Double_t GetAccidentalRate(Int_t,Int_t);
    Double_t GetLiHeRate(Int_t,Int_t);
    Double_t GetFNRate(Int_t,Int_t);
    Double_t GetAmCRate(Int_t,Int_t);
    
    Double_t GetAccidentalEvents(Int_t,Int_t);
    Double_t GetLiHeEvents(Int_t,Int_t);
    Double_t GetFNEvents(Int_t,Int_t);
    Double_t GetAmCEvents(Int_t,Int_t);
    
    Int_t GetNReactorPeriods();
};

NominalData :: NominalData(bool ish,Int_t dataSet)
{
    DataSet = dataSet;
    isH = ish;
    
    BCW=0;
    LBNL=0;
    Unified=1;
    
    TurnOnBudget = 0;
    TurnOffBudget = 0;
    
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
    RelativeEnergyScaleMatrix=0;
    AbsoluteEnergyScaleMatrix=0;
    RelativeEnergyOffsetMatrix=0;
    AbsoluteEnergyOffsetMatrix=0;
    IAVMatrix=0;
    NLMatrix=0;
    ResolutionMatrix=0;
    Sin22t12Matrix=0;
    
    VaryAccidentalBudget=0;
    VaryLiHeBudget=0;
    VaryFastNeutronsBudget=0;
    VaryAmCBudget=0;
    DistortLiHeBudget=0;
    DistortFastNeutronsBudget=0;
    DistortAmCBudget=0;
    //Systematics
    IsotopeBudget=0;
    ReactorPowerBudget=0;
    RelativeEnergyScaleBudget=0;
    IAVBudget=0;
    NLBudget=0;
    ResolutionBudget=0;
    Sin22t12Budget=0;
    EfficiencyBudget=0;
    
    ToyMC=0;//Data is default;
    Nweeks = 1;
    NReactorPeriods=20;
    NADs = 6;
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    ProtonsPerKton = 7.1638e31; //protons per kton
    DetectorMassLs[0] = 0.0215735;// kton
    DetectorMassLs[1] = 0.0215196;
    DetectorMassLs[2] = 0.0215872;
    DetectorMassLs[3] = 0.0215662;
    DetectorMassLs[4] = 0.0214088;
    DetectorMassLs[5] = 0.0216526;
    
    DetectorMassGdLs[0] = 0.019941;// kton
    DetectorMassGdLs[1] = 0.019966;
    DetectorMassGdLs[2] = 0.019891;
    DetectorMassGdLs[3] = 0.019913;
    DetectorMassGdLs[4] = 0.019991;
    DetectorMassGdLs[5] = 0.019892;
    
    if(TestAllTheSame)
    {
        for(Int_t i = 0; i<NADs; i++)
        {
            DetectorMassGdLs[i]=0.02;
        }
    }
    
    m_detectorEfficiency_Dt = 0.986;
    m_detectorEfficiency_Ep = 0.9988;
    m_detectorEfficiency_Ed_nominal = 0.909;
    m_detectorEfficiency_flash = 0.9998;
    m_detectorEfficiency_nGd = 0.838;
    m_detectorEfficiency_spill = 1.050;
    
    DetectorEfficiencyRelativeError = 0.0014;

    NSamples = 500;
    NSteps = 101;
    
    Emin = 1.8;
    Emax = 12;
    EVisMax = 12;
    LinearBinning=0;
    CalculateBinning();//default
    
    //Nominal isotope fraction errors:
    for(Int_t i=0;i<NIsotopes;i++)
    {
        IsotopeFracError[i]=0.05;
    }
    
    //Reactor power (It's an average of first 140 days of data) This should be coded to get the data from a txt file instead.)
    
    //test:
    if(TestAllTheSame)
    {
        ReactorPower[0] = 2.178;// GW
        ReactorPower[1] = 2.178;
        ReactorPower[2] = 2.178;
        ReactorPower[3] = 2.178;
        ReactorPower[4] = 2.178;
        ReactorPower[5] = 2.178;
        
        IsotopeFrac[0] = 0.25; // NmlU235Frac
        IsotopeFrac[1] = 0.25; // NmlU238Frac
        IsotopeFrac[2] = 0.25; // NmlPu239Frac
        IsotopeFrac[3] = 0.25; // NmlPu241Frac
    }
    else
    {
        ReactorPower[0] = 2.178;// GW
        ReactorPower[1] = 2.876;
        ReactorPower[2] = 2.389;
        ReactorPower[3] = 2.444;
        ReactorPower[4] = 2.820;
        ReactorPower[5] = 2.787;
        
        //Nominal isotope fractions
        IsotopeFrac[0] = 0.64; // NmlU235Frac
        IsotopeFrac[1] = 0.08; // NmlU238Frac
        IsotopeFrac[2] = 0.25; // NmlPu239Frac
        IsotopeFrac[3] = 0.03; // NmlPu241Frac
    }
    //Reactor power errors:
    for(Int_t r=0;r<NReactors;r++)
    {
        ReactorPowerError[r]=0.005;
    }
    
    StatisticalFluctuation=0;
    
    //Energy per fission
    EnergyPerFission[0]=201.92;//MeV/fission
    EnergyPerFission[1]=205.52;
    EnergyPerFission[2]=209.99;
    EnergyPerFission[3]=213.60;
    
    //Energy per fission errors
    EnergyPerFissionError[0]=0.002278;
    EnergyPerFissionError[1]=0.004671;
    EnergyPerFissionError[2]=0.002857;
    EnergyPerFissionError[3]=0.003043;
    
    //Resolution errors
    ResolutionError=0.002;//Due to parameter uncertainty
    ResolutionErrorUncorrelated=0.002;//Due to energy scale difference between detectors
    
    //Energy Scale
    m_abs_escale =1;
    m_abs_escale_error = 0.01;
    
    m_abs_eoffset = 0.0;
    m_abs_eoffset_error = 0.08;  //(MeV)
    
    m_rel_eoffset_error = 0.013; //  (MeV)
    
    for(Int_t AD=0; AD<NADs; AD++)
    {
        m_rel_escale[AD] = 1.0;
        m_rel_escale_error[AD] = 0.0035; // 0.35%
        m_rel_escale_nominal[AD] = m_rel_escale[AD];
        m_rel_eoffset[AD] = 0.0;
    }
    
    //Oscillation parameters and errors:from PDG 2013
    s22t13 = 0.09;//0.09 ± 0.009
    s22t13Error = 0.01;
    s22t12 =0.857;//0.857 ± 0.024
    s22t12Error = 0.024;
    dm2_32 = 2.41e-3;//eV2
    dm2_21 = 7.50e-5;//eV2

    hierarchy=1;//-1 for inverted           // Set as external input to be set in the GUI   !!!!
    
    dm2_31=dm2_32+hierarchy*dm2_21;
    dm2_ee = dm2_32+hierarchy*5.21e-5;
    
    sin_start = 0;
    sin_end = 0.2;
    
    dmee_start = dm2_ee - 0.001;
    dmee_end = dm2_ee + 0.001;
    
    //For Hydrogen and 6 ADs: from http://dayabay.ihep.ac.cn/DocDB/0085/008556/020/Main-version3.pdf
    
    if(this->GetAnalysis())
    {
        if(DataSet==2)
        {            
            //p12b values
            if(Nweeks==1)
            {
                ObservedEvents[0]=74136;
                ObservedEvents[1]=74783;
                ObservedEvents[2]=69083;
                ObservedEvents[3]=20218;
                ObservedEvents[4]=20366;
                ObservedEvents[5]=21527;
                
                FullTime[0] = 190.9954;// Days;
                FullTime[1] = 190.9954;// Days;
                FullTime[2] = 189.6464 ;// Days;
                FullTime[3] = 189.7857;// Days;
                FullTime[4] = 189.7857;// Days;
                FullTime[5] = 189.7857;// Days;
                
                MultiEff[0]= 0.9836;
                MultiEff[1]= 0.9837;
                MultiEff[2]= 0.9842;
                MultiEff[3]= 0.9832;
                MultiEff[4]= 0.9827;
                MultiEff[5]= 0.9827;
                
                MuonEff[0]= 0.7946;
                MuonEff[1]= 0.7912;
                MuonEff[2]= 0.8338;
                MuonEff[3]= 0.9815;
                MuonEff[4]= 0.9815;
                MuonEff[5]= 0.9812;
                
                //Hydrogen data here, taken from the εμεm corrected data and dividing by that coefficient:

                AccidentalError[0]=0.10160843525865;
                AccidentalError[1]=0.101182017951218;
                AccidentalError[2]=0.0902700154097309;
                AccidentalError[3]=0.057902182109611;
                AccidentalError[4]=0.0578725809910863;
                AccidentalError[5]=0.0674940890483319;
                
                FastNeutronError[0]=0.211032903998735;
                FastNeutronError[1]=0.210147268052529;
                FastNeutronError[2]=0.2379845860802;
                FastNeutronError[3]=0.0386014547397407;
                FastNeutronError[4]=0.0385817206607242;
                FastNeutronError[5]=0.0385680508847611;
            
                LiHeError[0]=1.07861262043798;
                LiHeError[1]=1.07408603671293;
                LiHeError[2]=0.878081058985565;
                LiHeError[3]=0.125454727904157;
                LiHeError[4]=0.125390592147354;
                LiHeError[5]=0.125346165375473;
                
                AmCError[0]=0.0390801674071732;
                AmCError[1]=0.0311329286003747;
                AmCError[2]=0.0328254601489931;
                AmCError[3]=0.0289510910548055;
                AmCError[4]=0.0289362904955431;
                AmCError[5]=0.0289260381635708;
                
                //Hydrogen data here, taken from the εμεm corrected data and dividing by that coefficient:
            
                AccidentalRate[0]=50.778238638208;// per day
                AccidentalRate[1]=49.8577452650692;
                AccidentalRate[2]=47.2842616574847;
                AccidentalRate[3]=59.9340203187068;
                AccidentalRate[4]=61.7855296790011;
                AccidentalRate[5]=65.7609082243815;
                
                FastNeutronRate[0]=1.63355099761984;// per day
                FastNeutronRate[1]=1.62669551936958;
                FastNeutronRate[2]=1.12427201010301;
                FastNeutronRate[3]=0.0965036368493517;
                FastNeutronRate[4]=0.0964543016518104;
                FastNeutronRate[5]=0.0964201272119027;
                
                LiHeRate[0]=2.14940920739452;// per day
                LiHeRate[1]=2.14038884127576;
                LiHeRate[2]=1.75616211797113;
                LiHeRate[3]=0.250909455808314;
                LiHeRate[4]=0.250781184294707;
                LiHeRate[5]=0.250692330750947;
                
                AmCRate[0] = 0.0703443013329117;//per day
                AmCRate[1] = 0.070049089350843;
                AmCRate[2] = 0.0738572853352344;
                AmCRate[3] = 0.057902182109611;
                AmCRate[4] = 0.0578725809910863;
                AmCRate[5] = 0.0578520763271416;
            }
            else
            {
                std::cout << "need to do some averaging sum like yasu's, for reference I add his lines of code here " << std::endl;
                exit(EXIT_FAILURE);
//                if (FirstMakeSuperPrediction){
//                    sprintf(dummyname,"CombCorrEvtsSpec_%i",idet);
//                    CombCorrEvtsSpec[idet] = (TH1F*)tdper[0].CorrEvtsSpec[idet]->Clone(dummyname);
//                    sprintf(dummyname,"CombCorrBgEvtsSpec_%i",idet);
//                    CombCorrBgEvtsSpec[idet] = (TH1F*)tdper[0].CorrBgEvtsSpec[idet]->Clone(dummyname);
//                    //note: do not set FirstMakeSuperPrediction to false as that is done below
//                }
//                CombCorrEvtsSpec[idet]->Reset();
//                CombCorrBgEvtsSpec[idet]->Reset();
//                
//                for(int ii=0;ii<Nperiods;++ii){
//                    float factor=tdper[ii].MuonVetoEff[idet]
//                    *tdper[ii].DMCEff[idet]
//                    *tdper[ii].Livetime[idet]
//                    *tdper[ii].TargetMass[idet]/tdper[ii].TargetMass[0];
//                    
//                    weightedsum+=tdper[ii].CorrEvts[idet]*factor;
//                    livsum+=factor;
//                    weightederr+=pow(tdper[ii].ErrEvts[idet]*factor,2);
//                    weightedbg+=tdper[ii].CorrBgEvts[idet]*factor;
//                    
//                    CombCorrEvtsSpec[idet]->Add(tdper[ii].CorrEvtsSpec[idet],factor);
//                    CombCorrBgEvtsSpec[idet]->Add(tdper[ii].CorrBgEvtsSpec[idet],factor);
//                    
//                    // weightedsum+=tdper[ii].CorrEvts[idet]*tdper[ii].Livetime[idet];
//                    // livsum+=tdper[ii].Livetime[idet];
//                    // weightederr+=pow(tdper[ii].ErrEvts[idet]*tdper[ii].Livetime[idet],2);
//                    // weightedbg+=tdper[ii].CorrBgEvts[idet]*tdper[ii].Livetime[idet];
//                    
//                    // CombCorrEvtsSpec[idet]->Add(tdper[ii].CorrEvtsSpec[idet],tdper[ii].Livetime[idet]);
//                    
//                }
//                CombLivetime[idet]=livsum;
//                CombCorrEvts[idet]=weightedsum*1./livsum;
//                CombErrEvts[idet]=sqrt(weightederr)*1./livsum;
//                CombCorrBgEvts[idet]=weightedbg*1./livsum;
//                CombCorrEvtsSpec[idet]->Scale(1./livsum);
//                CombCorrBgEvtsSpec[idet]->Scale(1./livsum);
            }
        }
    }
    //For Gd and 6 ADs: http://dayabay.ihep.ac.cn/DocDB/0085/008556/016/Main-version2.pdf page 5
    else
    {
        if(DataSet==2)
        {
            // P12E Values in Theta13-inputs_32week_inclusive.txt
            if(Nweeks==1)
            {
                ObservedEvents[0]=101362;
                ObservedEvents[1]=102605;
                ObservedEvents[2]=92960;
                ObservedEvents[3]=13965;
                ObservedEvents[4]=13894;
                ObservedEvents[5]=13731;
                
                // Efficiencies, FullTime, etc
                FullTime[0] = 190.99;// Days;
                FullTime[1] = 190.99;// Days;
                FullTime[2] = 189.65;// Days;
                FullTime[3] = 189.79;// Days;
                FullTime[4] = 189.79;// Days;
                FullTime[5] = 189.79;// Days;
                
                MultiEff[0]= 0.97586;
                MultiEff[1]= 0.9762;
                MultiEff[2]= 0.9774;
                MultiEff[3]= 0.97642;
                MultiEff[4]= 0.97615;
                MultiEff[5]= 0.97598;
                
                MuonEff[0]= 0.81583;
                MuonEff[1]= 0.81245;
                MuonEff[2]= 0.84761;
                MuonEff[3]= 0.98076;
                MuonEff[4]= 0.9802;
                MuonEff[5]= 0.98018;
                
                AccidentalRate[0]=9.5508;// per day
                AccidentalRate[1]=9.3657;
                AccidentalRate[2]=7.441;
                AccidentalRate[3]=2.9645;
                AccidentalRate[4]=2.9238;
                AccidentalRate[5]=2.8741;
                
                FastNeutronRate[0]=0.92;// per day
                FastNeutronRate[1]=0.92;
                FastNeutronRate[2]=0.62;
                FastNeutronRate[3]=0.04;
                FastNeutronRate[4]=0.04;
                FastNeutronRate[5]=0.04;
                
                LiHeRate[0]=2.4;// per day
                LiHeRate[1]=2.4;
                LiHeRate[2]=1.2;
                LiHeRate[3]=0.22;
                LiHeRate[4]=0.22;
                LiHeRate[5]=0.22;
                
                AmCRate[0] = 0.26;//per day
                AmCRate[1] = 0.26;
                AmCRate[2] = 0.26;
                AmCRate[3] = 0.26;
                AmCRate[4] = 0.26;
                AmCRate[5] = 0.26;
                
                //test:
                if(TestAllTheSame)
                {
                    for(Int_t i = 0; i<NADs; i++)
                    {
                        for(Int_t j = 0; j<NReactorPeriods; j++)
                        {
                            ObservedEvents[i+j*MaxDetectors]=100000;
                            MuonEff[i+j*MaxDetectors]=1;
                            MultiEff[i+j*MaxDetectors]=1;
                            FullTime[i+j*MaxDetectors]=100;
                            AccidentalRate[i+j*MaxDetectors]=0;
                            FastNeutronRate[i+j*MaxDetectors]=0;
                            LiHeRate[i+j*MaxDetectors]=0;
                            AmCRate[i+j*MaxDetectors]=0;
                        }
                        
                    }
                }

                
                AccidentalError[0]=0.074594;//Absolute uncertainty AD1
                AccidentalError[1]=0.06356;//Absolute uncertainty AD2
                AccidentalError[2]=0.046238;//Absolute uncertainty AD3
                AccidentalError[3]=0.037342;//Absolute uncertainty AD4
                AccidentalError[4]=0.037342;//Absolute uncertainty AD5
                AccidentalError[5]=0.037314;//Absolute uncertainty AD5
                
                FastNeutronError[0]=0.46;//Absolute uncertainty EH1
                FastNeutronError[1]=0.46;//Absolute uncertainty EH2
                FastNeutronError[2]=0.31;//Absolute uncertainty EH2
                FastNeutronError[3]=0.02;//Absolute uncertainty EH3
                FastNeutronError[4]=0.02;//Absolute uncertainty EH3
                FastNeutronError[5]=0.02;//Absolute uncertainty EH3
                
                LiHeError[0]=0.86;//Absolute uncertainty EH1
                LiHeError[1]=0.86;//Absolute uncertainty EH1
                LiHeError[2]=0.63;//Absolute uncertainty EH2
                LiHeError[3]=0.06;//Absolute uncertainty EH3
                LiHeError[4]=0.06;//Absolute uncertainty EH3
                LiHeError[5]=0.06;//Absolute uncertainty EH3
                
                AmCError[0]=0.12;//Absolute uncertainty EH1
                AmCError[1]=0.12;//Absolute uncertainty EH2
                AmCError[2]=0.12;//Absolute uncertainty EH3
                AmCError[3]=0.12;//Absolute uncertainty EH3
                AmCError[4]=0.12;//Absolute uncertainty EH3
                AmCError[5]=0.12;//Absolute uncertainty EH3
            }
            else
            {
                
                LoadMainData(Form("./Inputs/Theta13-inputs_20week.txt"));// Need to do it Nweek dependent
                
//                string mainmatrixname = "./Inputs/Theta13-inputs_20week.txt";
            }
        }
    }
    
    //IAV thickness relative error
    
    IAVError=0.04;
}

void NominalData :: CopyData(NominalData * data)
{
    ResponseDirectory = data->ResponseDirectory;
    PredictionsDirectory = data->PredictionsDirectory;
    ToyMCSamplesDirectory = data->ToyMCSamplesDirectory;
    SysCovDirectory = data->SysCovDirectory;
    BkgCovDirectory = data->BkgCovDirectory;
    
    //Analysis type:
     isH = data->isH ;
     DataSet = data->DataSet ;
     ToyMC = data->ToyMC ;
    //Binning type:
     LinearBinning = data->LinearBinning ;//0 is non linear (LBNL like binning), 1 is linear (51 bins);
    
    //Detectors:
     NADs = data->NADs ;
     ADsEH1 = data->ADsEH1;
     ADsEH2 = data->ADsEH2;
     ADsEH3 = data->ADsEH3;
    
     hierarchy = data->hierarchy;
     Combine = data->Combine;
    
     ProtonsPerKton = data->ProtonsPerKton; //protons per kton
     std::copy(std::begin(data->DetectorMassGdLs), std::end(data->DetectorMassGdLs), std::begin(DetectorMassGdLs));
     std::copy(std::begin(data->DetectorMassLs), std::end(data->DetectorMassLs), std::begin(DetectorMassLs));

     m_detectorEfficiency_Dt = data->m_detectorEfficiency_Dt;
     m_detectorEfficiency_Ep = data->m_detectorEfficiency_Ep;
     m_detectorEfficiency_Ed_nominal = data->m_detectorEfficiency_Ed_nominal;
     m_detectorEfficiency_flash = data->m_detectorEfficiency_flash;
     m_detectorEfficiency_nGd = data->m_detectorEfficiency_nGd;
     m_detectorEfficiency_spill = data->m_detectorEfficiency_spill;
    
    //Toy samples
     NSamples = data->NSamples;
     NSteps = data->NSteps;
    //Oscillation parameters
     s22t12 = data->s22t12;
     s22t12Error = data->s22t12Error;
     s22t13 = data->s22t13;
     s22t13Error = data->s22t13Error;
     dm2_ee = data->dm2_ee;
     dm2_31 = data->dm2_31;
     dm2_32 = data->dm2_32;
     dm2_21 = data->dm2_21;
     sin_start = data->sin_start;
     sin_end = data->sin_end;
     dmee_start = data->dmee_start;
     dmee_end = data->dmee_end;
     NReactorPeriods = data->NReactorPeriods;
    
    //Reactor parameters
    std::copy(std::begin(data->IsotopeFrac), std::end(data->IsotopeFrac), std::begin(IsotopeFrac));
    std::copy(std::begin(data->IsotopeFracError), std::end(data->IsotopeFracError), std::begin(IsotopeFracError));
    std::copy(std::begin(data->ReactorPower), std::end(data->ReactorPower), std::begin(ReactorPower));
    std::copy(std::begin(data->ReactorPowerError), std::end(data->ReactorPowerError), std::begin(ReactorPowerError));
    std::copy(std::begin(data->EnergyPerFission), std::end(data->EnergyPerFission), std::begin(EnergyPerFission));
    std::copy(std::begin(data->EnergyPerFissionError), std::end(data->EnergyPerFissionError), std::begin(EnergyPerFissionError));
    //Background relative errors, rates and events:
    std::copy(std::begin(data->AccidentalError), std::end(data->AccidentalError), std::begin(AccidentalError));
    std::copy(std::begin(data->LiHeError), std::end(data->LiHeError), std::begin(LiHeError));
    std::copy(std::begin(data->FastNeutronError), std::end(data->FastNeutronError), std::begin(FastNeutronError));
    std::copy(std::begin(data->AmCError), std::end(data->AmCError), std::begin(AmCError));

    std::copy(std::begin(data->AccidentalRate), std::end(data->AccidentalRate), std::begin(AccidentalRate));
    std::copy(std::begin(data->FastNeutronRate), std::end(data->FastNeutronRate), std::begin(FastNeutronRate));
    std::copy(std::begin(data->LiHeRate), std::end(data->LiHeRate), std::begin(LiHeRate));
    std::copy(std::begin(data->AmCRate), std::end(data->AmCRate), std::begin(AmCRate));

    std::copy(std::begin(data->AccidentalEvents), std::end(data->AccidentalEvents), std::begin(AccidentalEvents));
    std::copy(std::begin(data->FastNeutronEvents), std::end(data->FastNeutronEvents), std::begin(FastNeutronEvents));
    std::copy(std::begin(data->LiHeEvents), std::end(data->LiHeEvents), std::begin(LiHeEvents));
    std::copy(std::begin(data->AmCEvents), std::end(data->AmCEvents), std::begin(AmCEvents));

    std::copy(std::begin(data->ObservedEvents), std::end(data->ObservedEvents), std::begin(ObservedEvents));
    
    // Days and efficiencies:
    std::copy(std::begin(data->FullTime), std::end(data->FullTime), std::begin(FullTime));
    std::copy(std::begin(data->MuonEff), std::end(data->MuonEff), std::begin(MuonEff));
    std::copy(std::begin(data->MultiEff), std::end(data->MultiEff), std::begin(MultiEff));
    
    //IAV error:
     IAVError = data->IAVError;
    
    //Resolution errors:
     ResolutionError = data->ResolutionError;
     ResolutionErrorUncorrelated = data->ResolutionErrorUncorrelated;
    
    //Energy scale:
     m_abs_escale = data->m_abs_escale;
     m_abs_escale_error = data->m_abs_escale_error;
    
     m_abs_eoffset = data->m_abs_eoffset;
     m_abs_eoffset_error = data->m_abs_eoffset_error;
    
     m_rel_eoffset_error = data->m_rel_eoffset_error;
     std::copy(std::begin(data->m_rel_escale), std::end(data->m_rel_escale), std::begin(m_rel_escale));
     std::copy(std::begin(data->m_rel_escale_error), std::end(data->m_rel_escale_error), std::begin(m_rel_escale_error));
     std::copy(std::begin(data->m_rel_escale_nominal), std::end(data->m_rel_escale_nominal), std::begin(m_rel_escale_nominal));
     std::copy(std::begin(data->m_rel_eoffset), std::end(data->m_rel_eoffset), std::begin(m_rel_eoffset));
    
     DetectorEfficiencyRelativeError = data->DetectorEfficiencyRelativeError;
    
    //Binning variables:
     Nweeks = data->Nweeks;
     NBins = data->NBins;
     Emin = data->Emin;
     Emax = data->Emax;
     EVisMin = data->EVisMin;
     EVisMax = data->EVisMax;
     n_etrue_bins = data->n_etrue_bins;
     n_evis_bins = data->n_evis_bins;
     std::copy(std::begin(data->enu_bins), std::end(data->enu_bins), std::begin(enu_bins));
     std::copy(std::begin(data->evis_bins), std::end(data->evis_bins), std::begin(evis_bins));

    //Reactor
     binWidth = data->binWidth;
     std::copy(std::begin(data->m_dNdE_nom), std::end(data->m_dNdE_nom), std::begin(m_dNdE_nom));

     m_nSamples = data->m_nSamples;
     m_eMin = data->m_eMin;
     m_eMax = data->m_eMax;
    
    //Reactor Covariance Matrix
     std::copy(std::begin(data->L), std::end(data->L), std::begin(L));
    
     TurnOnBudget = data->TurnOnBudget;
     TurnOffBudget = data->TurnOffBudget;
    
     BCW = data->BCW;
     LBNL = data->LBNL;
     Unified = data->Unified;
    
     VaryAccidentalMatrix = data->VaryAccidentalMatrix;
     VaryLiHeMatrix = data->VaryLiHeMatrix;
     VaryFastNeutronsMatrix = data->VaryFastNeutronsMatrix;
     VaryAmCMatrix = data->VaryAmCMatrix;
     DistortLiHeMatrix = data->DistortLiHeMatrix;
     DistortFastNeutronsMatrix = data->DistortFastNeutronsMatrix;
     DistortAmCMatrix = data->DistortAmCMatrix;
     IsotopeMatrix = data->IsotopeMatrix;
     ReactorPowerMatrix = data->ReactorPowerMatrix;
     RelativeEnergyOffsetMatrix = data->RelativeEnergyOffsetMatrix;
     AbsoluteEnergyOffsetMatrix = data->AbsoluteEnergyOffsetMatrix;
     AbsoluteEnergyScaleMatrix = data->AbsoluteEnergyScaleMatrix;
     RelativeEnergyScaleMatrix = data->RelativeEnergyScaleMatrix;
     IAVMatrix = data->IAVMatrix;
     NLMatrix = data->NLMatrix;
     ResolutionMatrix = data->ResolutionMatrix;
     Sin22t12Matrix = data->Sin22t12Matrix;
     EfficiencyMatrix = data->EfficiencyMatrix;
    
     VaryAccidentalBudget = data->VaryAccidentalBudget;
     VaryLiHeBudget = data->VaryLiHeBudget;
     VaryFastNeutronsBudget = data->VaryFastNeutronsBudget;
     VaryAmCBudget = data->VaryAmCBudget;
     DistortLiHeBudget = data->DistortLiHeBudget;
     DistortFastNeutronsBudget = data->DistortFastNeutronsBudget;
     DistortAmCBudget = data->DistortAmCBudget;
     IsotopeBudget = data->IsotopeBudget;
     ReactorPowerBudget = data->ReactorPowerBudget;
     RelativeEnergyOffsetBudget = data->RelativeEnergyOffsetBudget;
     AbsoluteEnergyOffsetBudget = data->AbsoluteEnergyOffsetBudget;
     AbsoluteEnergyScaleBudget = data->AbsoluteEnergyScaleBudget;
     RelativeEnergyScaleBudget = data->RelativeEnergyScaleBudget;
     IAVBudget = data->IAVBudget;
     NLBudget = data->NLBudget;
     ResolutionBudget = data->ResolutionBudget;
     Sin22t12Budget = data->Sin22t12Budget;
     EfficiencyBudget = data->EfficiencyBudget;
     SystematicBudget = data->SystematicBudget;
     BackgroundBudget = data->BackgroundBudget;
     TotalBudget = data->TotalBudget;
    
     StatisticalFluctuation = data->StatisticalFluctuation;
     UseToyMCTree = data->UseToyMCTree;
    
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                              Analysis, Binning and Model controllers
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NominalData :: SetDataSet(Int_t dataset)
{
    DataSet=dataset;
}

void NominalData :: SetVaryAccidentalMatrix(bool VaryAccMatrix)
{
    VaryAccidentalMatrix=VaryAccMatrix;
}

void NominalData :: SetVaryLiHeMatrix(bool VaryLiMatrix)
{
    VaryLiHeMatrix=VaryLiMatrix;
}

void NominalData :: SetVaryFastNeutronsMatrix(bool VaryFNmatrix)
{
    VaryFastNeutronsMatrix=VaryFNmatrix;
}

void NominalData :: SetVaryAmCMatrix(bool VaryAmCmatrix)
{
    VaryAmCMatrix=VaryAmCmatrix;
}

void NominalData :: SetDistortLiHeMatrix(bool DistortLiMatrix)
{
    DistortLiHeMatrix=DistortLiMatrix;
}

void NominalData :: SetDistortFastNeutronsMatrix(bool DistortFNmatrix)
{
    DistortFastNeutronsMatrix=DistortFNmatrix;
}

void NominalData :: SetDistortAmCMatrix(bool DistortAmCmatrix)
{
    DistortAmCMatrix=DistortAmCmatrix;
}

void NominalData :: SetIsotopeMatrix(bool Isotopematrix)
{
    IsotopeMatrix=Isotopematrix;
}

void NominalData :: SetReactorPowerMatrix(bool ReactorPowermatrix)
{
    ReactorPowerMatrix=ReactorPowermatrix;
}

void NominalData :: SetRelativeEnergyOffsetMatrix(bool RelativeEnergyOffsetmatrix)
{
    RelativeEnergyOffsetMatrix = RelativeEnergyOffsetmatrix;
}

void NominalData :: SetAbsoluteEnergyOffsetMatrix(bool AbsoluteEnergyOffsetmatrix)
{
    AbsoluteEnergyOffsetMatrix = AbsoluteEnergyOffsetmatrix;
}

void NominalData :: SetAbsoluteEnergyScaleMatrix(bool AbsoluteEnergyScalematrix)
{
    AbsoluteEnergyScaleMatrix=AbsoluteEnergyScalematrix;
}

void NominalData :: SetRelativeEnergyScaleMatrix(bool RelativeEnergyScalematrix)
{
    RelativeEnergyScaleMatrix=RelativeEnergyScalematrix;
}

void NominalData :: SetIAVMatrix(bool IAVmatrix)
{
    IAVMatrix=IAVmatrix;
}

void NominalData :: SetNLMatrix(bool NLmatrix)
{
    NLMatrix=NLmatrix;
}

void NominalData :: SetResolutionMatrix(bool Resolutionmatrix)
{
    ResolutionMatrix=Resolutionmatrix;
}

void NominalData :: SetSin22t12Matrix(bool Sin22t12matrix)
{
    Sin22t12Matrix=Sin22t12matrix;
}

void NominalData :: SetEfficiencyMatrix(bool Efficiencymatrix)
{
    EfficiencyMatrix = Efficiencymatrix;
}

void NominalData :: SetVaryAccidentalBudget(bool VaryAccBudget)
{
    VaryAccidentalBudget=VaryAccBudget;
}

void NominalData :: SetVaryLiHeBudget(bool VaryLiBudget)
{
    VaryLiHeBudget=VaryLiBudget;
}

void NominalData :: SetVaryFastNeutronsBudget(bool VaryFNBudget)
{
    VaryFastNeutronsBudget = VaryFNBudget;
}

void NominalData :: SetVaryAmCBudget(bool varyAmCBudget)
{
    VaryAmCBudget = varyAmCBudget;
}

void NominalData :: SetDistortLiHeBudget(bool DistortLiBudget)
{
    DistortLiHeBudget=DistortLiBudget;
}

void NominalData :: SetDistortFastNeutronsBudget(bool DistortFNBudget)
{
    DistortFastNeutronsBudget=DistortFNBudget;
}

void NominalData :: SetDistortAmCBudget(bool distortAmCBudget)
{
    DistortAmCBudget=distortAmCBudget;
}

void NominalData :: SetIsotopeBudget(bool isotopeBudget)
{
    IsotopeBudget=isotopeBudget;
}

void NominalData :: SetReactorPowerBudget(bool reactorPowerBudget)
{
    ReactorPowerBudget=reactorPowerBudget;
}

void NominalData :: SetRelativeEnergyOffsetBudget(bool relativeEnergyOffsetBudget)
{
    RelativeEnergyOffsetBudget = relativeEnergyOffsetBudget;
}

void NominalData :: SetAbsoluteEnergyOffsetBudget(bool absoluteEnergyOffsetBudget)
{
    AbsoluteEnergyOffsetBudget = absoluteEnergyOffsetBudget;
}

void NominalData :: SetAbsoluteEnergyScaleBudget(bool absoluteEnergyScaleBudget)
{
    AbsoluteEnergyScaleBudget=absoluteEnergyScaleBudget;
}

void NominalData :: SetRelativeEnergyScaleBudget(bool relativeEnergyScaleBudget)
{
    RelativeEnergyScaleBudget=relativeEnergyScaleBudget;
}

void NominalData :: SetIAVBudget(bool iAVBudget)
{
    IAVBudget=iAVBudget;
}

void NominalData :: SetNLBudget(bool nLBudget)
{
    NLBudget=nLBudget;
}

void NominalData :: SetResolutionBudget(bool resolutionBudget)
{
    ResolutionBudget=resolutionBudget;
}

void NominalData :: SetSin22t12Budget(bool sin22t12budget)
{
    Sin22t12Budget=sin22t12budget;
}

void NominalData :: SetEfficiencyBudget(bool Efficiencybudget)
{
    EfficiencyBudget=Efficiencybudget;
}

void NominalData :: SetSystematicBudget(bool systematicBudget)
{
    SystematicBudget=systematicBudget;
}

void NominalData :: SetBackgroundBudget(bool backgroundBudget)
{
    BackgroundBudget=backgroundBudget;
}

void NominalData :: SetTotalBudget(bool totalBudget)
{
    TotalBudget=totalBudget;
}

void NominalData :: SetBCWModel(bool bcw)
{
    BCW=bcw;
}

void NominalData :: SetLBNLModel(bool lbnl)
{
    LBNL=lbnl;
}

void NominalData :: SetUnifiedModel(bool unified)
{
    Unified=unified;
}

Int_t NominalData :: GetDataSet()
{
    return DataSet;
}

bool NominalData :: GetVaryAccidentalMatrix()
{
    return VaryAccidentalMatrix;
}

bool NominalData :: GetVaryLiHeMatrix()
{
    return VaryLiHeMatrix;
}

bool NominalData :: GetVaryFastNeutronsMatrix()
{
    return VaryFastNeutronsMatrix;
}

bool NominalData :: GetVaryAmCMatrix()
{
    return VaryAmCMatrix;
}

bool NominalData :: GetDistortLiHeMatrix()
{
    return DistortLiHeMatrix;
}

bool NominalData :: GetDistortFastNeutronsMatrix()
{
    return DistortFastNeutronsMatrix;
}

bool NominalData :: GetDistortAmCMatrix()
{
    return DistortAmCMatrix;
}

bool NominalData :: GetIsotopeMatrix()
{
    return IsotopeMatrix;
}

bool NominalData :: GetReactorPowerMatrix()
{
    return ReactorPowerMatrix;
}

bool NominalData :: GetRelativeEnergyOffsetMatrix()
{
    return RelativeEnergyOffsetMatrix;
}

bool NominalData :: GetAbsoluteEnergyOffsetMatrix()
{
    return AbsoluteEnergyOffsetMatrix;
}

bool NominalData :: GetAbsoluteEnergyScaleMatrix()
{
    return AbsoluteEnergyScaleMatrix;
}

bool NominalData :: GetRelativeEnergyScaleMatrix()
{
    return RelativeEnergyScaleMatrix;
}

bool NominalData :: GetIAVMatrix()
{
    return IAVMatrix;
}

bool NominalData :: GetNLMatrix()
{
    return NLMatrix;
}

bool NominalData :: GetResolutionMatrix()
{
    return ResolutionMatrix;
}

bool NominalData :: GetSin22t12Matrix()
{
    return Sin22t12Matrix;
}

bool NominalData :: GetEfficiencyMatrix()
{
    return EfficiencyMatrix;
}

bool NominalData :: GetVaryAccidentalBudget()
{
    return VaryAccidentalBudget;
}

bool NominalData :: GetVaryLiHeBudget()
{
    return VaryLiHeBudget;
}

bool NominalData :: GetVaryFastNeutronsBudget()
{
    return VaryFastNeutronsBudget;
}

bool NominalData :: GetVaryAmCBudget()
{
    return VaryAmCBudget;
}

bool NominalData :: GetDistortLiHeBudget()
{
    return DistortLiHeBudget;
}

bool NominalData :: GetDistortFastNeutronsBudget()
{
    return DistortFastNeutronsBudget;
}

bool NominalData :: GetDistortAmCBudget()
{
    return DistortAmCBudget;
}

bool NominalData :: GetIsotopeBudget()
{
    return IsotopeBudget;
}

bool NominalData :: GetReactorPowerBudget()
{
    return ReactorPowerBudget;
}

bool NominalData :: GetRelativeEnergyOffsetBudget()
{
    return RelativeEnergyOffsetBudget;
}

bool NominalData :: GetAbsoluteEnergyOffsetBudget()
{
    return AbsoluteEnergyOffsetBudget;
}

bool NominalData :: GetAbsoluteEnergyScaleBudget()
{
    return AbsoluteEnergyScaleBudget;
}

bool NominalData :: GetRelativeEnergyScaleBudget()
{
    return RelativeEnergyScaleBudget;
}

bool NominalData :: GetIAVBudget()
{
    return IAVBudget;
}

bool NominalData :: GetNLBudget()
{
    return NLBudget;
}

bool NominalData :: GetResolutionBudget()
{
    return ResolutionBudget;
}

bool NominalData :: GetSin22t12Budget()
{
    return Sin22t12Budget;
}

bool NominalData :: GetEfficiencyBudget()
{
    return EfficiencyBudget;
}

bool NominalData :: GetSystematicBudget()
{
    return SystematicBudget;
}

bool NominalData :: GetBackgroundBudget()
{
    return BackgroundBudget;
}

bool NominalData :: GetTotalBudget()
{
    return TotalBudget;
}

bool NominalData :: GetBCWModel()
{
    return BCW;
}

bool NominalData :: GetLBNLModel()
{
    return LBNL;
}

bool NominalData :: GetUnifiedModel()
{
    return Unified;
}

Int_t NominalData :: GetNSamples()
{
    return NSamples;
}

Int_t NominalData :: GetNSteps()
{
    return NSteps;
}

Int_t NominalData :: GetADs()
{
    return NADs;
}

bool NominalData :: GetAnalysis()
{
    return isH;
}

bool NominalData :: GetToyMC()
{
    return ToyMC;
}

bool NominalData :: GetStatisticalFluctuation()
{
    return StatisticalFluctuation;
}

void NominalData :: SetAnalysis(bool IsH)
{
    isH = IsH;
}

void NominalData :: SetADs(Int_t ads)
{
    NADs = ads;
}

void NominalData :: SetNSamples(Int_t nsamples)
{
    NSamples = nsamples;
}

void NominalData :: SetNSteps(Int_t nsteps)
{
    NSteps = nsteps;
}

void NominalData :: SetToyMC(bool toyMC)
{
    ToyMC = toyMC;
}


void NominalData :: SetStatisticalFluctuation(bool statisticalfluctuation)
{
    StatisticalFluctuation = statisticalfluctuation;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                              Binning inputs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NominalData :: SetBinning(bool linearBinning)
{
    LinearBinning=linearBinning;
}

void NominalData :: SetEmin(Double_t emin)
{
    Emin = emin;
}

void NominalData :: SetEmax(Double_t emax)
{
    Emax = emax;
}

void NominalData :: SetEVisMin(Double_t evismin)
{
    EVisMin = evismin;
}

void NominalData :: SetEVisMax(Double_t evismax)
{
    EVisMax = evismax;
}

bool  NominalData :: GetBinning()
{
    return LinearBinning;
}

Double_t NominalData :: GetEmin()
{
    return Emin;
}

Double_t NominalData :: GetEmax()
{
    return Emax;
}

Double_t NominalData :: GetEVisMin()
{
    return EVisMin;
}

Double_t NominalData :: GetEVisMax()
{
    return EVisMax;
}

void NominalData :: SetCombineMode(Int_t combine)
{
    Combine = combine;
}

Int_t NominalData :: GetCombineMode()
{
    return Combine;
}

void NominalData :: SetWeeks(Int_t nweeks)
{
    Nweeks = nweeks;
}

Int_t NominalData :: GetWeeks()
{
    return Nweeks;
    
//    return 1; // to check fluxH
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                              Reactor inputs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t NominalData :: GetIsotopeFraction(Int_t isotope)
{
    return IsotopeFrac[isotope];
}

Double_t NominalData :: GetIsotopeFractionError(Int_t isotope)
{
    return IsotopeFracError[isotope];
}

Double_t NominalData :: GetReactorPower(Int_t reactor)
{
    return ReactorPower[reactor];
}

Double_t NominalData :: GetReactorPowerError(Int_t reactor)
{
    return ReactorPowerError[reactor];
}

Double_t NominalData :: GetEnergyPerFission(Int_t isotope)
{
    return EnergyPerFission[isotope];
}

Double_t NominalData :: GetEnergyPerFissionError(Int_t isotope)
{
    return EnergyPerFissionError[isotope];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                              Oscillation inputs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NominalData :: SetSin22t13(Double_t S22t13)
{
    s22t13 = S22t13;
}

Double_t NominalData :: GetSin22t13()
{
    return s22t13;
}

void NominalData :: SetDm2ee(Double_t Dm2ee)
{
    dm2_ee = Dm2ee;
}

void NominalData :: SetDm231(Double_t Dm231)
{
    dm2_31 = Dm231;
}

void NominalData :: SetDm232(Double_t Dm232)
{
    dm2_32 = Dm232;
}

void NominalData :: SetDm221(Double_t Dm221)
{
    dm2_21 = Dm221;
}

Double_t NominalData :: GetDm2ee()
{
    return dm2_ee;
}

Double_t NominalData :: GetDm231()
{
    return dm2_31;
}

Double_t NominalData :: GetDm232()
{
    return dm2_32;
}

Double_t NominalData :: GetDm221()
{
    return dm2_21;
}

void NominalData :: SetSin22t12(Double_t S22t12)
{
    s22t12 = S22t12;
}

Double_t NominalData :: GetSin22t12()
{
    return s22t12;
}

Double_t NominalData :: GetSin22t12Error()
{
    return s22t12Error;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                              Detector inputs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t NominalData :: GetAbsoluteEnergyScale()
{
    return m_abs_escale;
}

Double_t NominalData :: GetAbsoluteEnergyScaleError()
{
    return m_abs_escale_error;
}

Double_t NominalData :: GetAbsoluteEnergyOffset()
{
    return m_abs_eoffset;
}

Double_t NominalData :: GetAbsoluteEnergyOffsetError()
{
    return m_abs_eoffset_error;
}

Double_t NominalData :: GetRelativeEnergyScale(Int_t AD)
{
    return m_rel_escale[AD];
}

Double_t NominalData :: GetRelativeEnergyError(Int_t AD)
{
    return m_rel_escale_error[AD];
}

Double_t NominalData :: GetRelativeEnergyOffset(Int_t AD)
{
    return m_rel_eoffset[AD];
}

Double_t NominalData :: GetRelativeEnergyOffsetError()
{
    return m_rel_eoffset_error;
}

Double_t NominalData :: GetResolutionError()
{
    return ResolutionError;
}

Double_t NominalData :: GetResoUncorrelatedError()
{
    return ResolutionErrorUncorrelated;
}

Double_t NominalData :: GetDetectorProtonsGdLs(Int_t detector)
{
    return ProtonsPerKton*DetectorMassGdLs[detector];
}

Double_t NominalData :: GetDetectorProtonsLs(Int_t detector)
{
    return ProtonsPerKton*DetectorMassLs[detector];
}


Double_t NominalData :: GetDetectorEfficiency(Int_t detector, Int_t week)
{
    //This is for nGd, will need to adapt for nH
    return MuonEff[detector+week*MaxDetectors]*MultiEff[detector+week*MaxDetectors]*m_detectorEfficiency_spill*m_detectorEfficiency_nGd*m_detectorEfficiency_Dt*m_detectorEfficiency_Ep*m_detectorEfficiency_Ed_nominal*m_detectorEfficiency_flash;
}

Double_t NominalData :: GetEnergyDelayedCutDetectorEfficiency()
{
    return m_detectorEfficiency_Ed_nominal;
}
Double_t NominalData :: GetDetectorEfficiencyRelativeError()
{
    return DetectorEfficiencyRelativeError;
}
Double_t NominalData :: GetFullTime(Int_t detector,Int_t week)
{
    return FullTime[detector+week*MaxDetectors];
}

Double_t NominalData :: GetIAVError()
{
    return IAVError;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                              Background inputs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Errors given as absolute error/(corrected rate) = relative error.
Double_t NominalData ::  GetAccidentalError(Int_t ad,Int_t week)
{
    return (AccidentalError[ad+week*MaxDetectors])/(AccidentalRate[ad+week*MaxDetectors]);
}

Double_t NominalData :: GetLiHeError(Int_t ad,Int_t week)
{
    return (LiHeError[ad+week*MaxDetectors])/(LiHeRate[ad+week*MaxDetectors]);
}

Double_t NominalData :: GetFNError(Int_t ad,Int_t week)
{
    return (FastNeutronError[ad+week*MaxDetectors])/(FastNeutronRate[ad+week*MaxDetectors]);
}

Double_t NominalData :: GetAmCError(Int_t ad,Int_t week)
{
    return (AmCError[ad+week*MaxDetectors])/(AmCRate[ad+week*MaxDetectors]);
}

Double_t NominalData :: GetAccidentalRate(Int_t ad,Int_t week)
{
    return AccidentalRate[ad+week*MaxDetectors];
}

Double_t NominalData :: GetLiHeRate(Int_t ad,Int_t week)
{
    return LiHeRate[ad+week*MaxDetectors];
}

Double_t NominalData :: GetFNRate(Int_t ad,Int_t week)
{
    return FastNeutronRate[ad+week*MaxDetectors];
}

Double_t NominalData :: GetAmCRate(Int_t ad,Int_t week)
{
    return AmCRate[ad+week*MaxDetectors];
}

Double_t NominalData :: GetIBDEvents(Int_t ad,Int_t week)
{
    return ((ObservedEvents[ad+week*MaxDetectors]/(FullTime[ad+week*MaxDetectors]*MuonEff[ad+week*MaxDetectors]*MultiEff[ad+week*MaxDetectors]))-(AccidentalRate[ad+week*MaxDetectors]+LiHeRate[ad+week*MaxDetectors]+FastNeutronRate[ad+week*MaxDetectors]+AmCRate[ad+week*MaxDetectors]))*FullTime[ad+week*MaxDetectors]*MuonEff[ad+week*MaxDetectors]*MultiEff[ad+week*MaxDetectors];//Apply background suppression for corrected events, then uncorrect it.
}

Double_t NominalData ::  GetAccidentalEvents(Int_t ad,Int_t week)
{
    return AccidentalRate[ad+week*MaxDetectors]*FullTime[ad+week*MaxDetectors]*MuonEff[ad+week*MaxDetectors]*MultiEff[ad+week*MaxDetectors];
}

Double_t NominalData :: GetLiHeEvents(Int_t ad,Int_t week)
{
    return LiHeRate[ad+week*MaxDetectors]*FullTime[ad+week*MaxDetectors]*MuonEff[ad+week*MaxDetectors]*MultiEff[ad+week*MaxDetectors];
}

Double_t NominalData :: GetFNEvents(Int_t ad,Int_t week)
{
    return FastNeutronRate[ad+week*MaxDetectors]*FullTime[ad+week*MaxDetectors]*MuonEff[ad+week*MaxDetectors]*MultiEff[ad+week*MaxDetectors];
}

Double_t NominalData :: GetAmCEvents(Int_t ad,Int_t week)
{
    return AmCRate[ad+week*MaxDetectors]*FullTime[ad+week*MaxDetectors]*MuonEff[ad+week*MaxDetectors]*MultiEff[ad+week*MaxDetectors];
}

bool NominalData :: GetTurnOnBudget()
{
    return TurnOnBudget;
}

bool NominalData :: GetTurnOffBudget()
{
    return TurnOffBudget;
}

void NominalData :: SetTurnOnBudget(bool turnOnBudget)
{
    TurnOnBudget = turnOnBudget;
}

void NominalData :: SetTurnOffBudget(bool turnOffBudget)
{
    TurnOffBudget = turnOffBudget;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                  This will be used to compare LBNL results. This is their method to input data from a txt file
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void NominalData :: LoadMainData(const Char_t* mainmatrixname)
{
    std::string line;
    Int_t linenum=0;//<---caution: only increments for lines that do not begin with #
    Int_t week=0;
    ifstream mainfile(mainmatrixname);
    while(!mainfile.eof())
    {
        getline(mainfile,line);
        std::string firstchar = line.substr(0,1);
        
        if(firstchar=="#") continue;//<-- ignore lines with comments
        
        //Special numbers
        if(linenum == 0)
        {
            Nweeks=atoi(line.c_str());
            std::cout << "The number of periods is " << Nweeks << std::endl;
        }
        
        if(linenum == 1){}
        //        {
        //            if(atoi(firstchar.c_str())==0) isMC=true;
        //            if(atoi(firstchar.c_str())==1) isMC=false;
        //            cout << "Simflag: " << isMC << " (1-->MC, 0-->Data)" << endl;
        //        }
        
        if(linenum > 2)
        {
            std::cout << "reading " << line << std::endl;
            std::istringstream iss(line);
            Int_t row=0;
            Int_t column=0;
            Double_t readvals[6]={0};
            while(iss)
            {
                std::string sub; iss >> sub;
                if(column==0) week=atoi(sub.c_str());
                
                if(column==1) row=atoi(sub.c_str());
                
                if(column>1 && sub!="") readvals[column-2]=atof(sub.c_str());
                column+=1;
            }//looping over columns
            
            //-->dates
            //-->obs events
            if(row==1)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    ObservedEvents[AD+(week-1)*MaxDetectors]=readvals[AD];
                }
            }
            //
            //            //Calculate Statistical Errors
            //
            //            for(Int_t i=0; i<NADs; i+
            //            {
            //                ErrEvts[AD+(week-1)*MaxDetectors]=sqrt(ObsEvts[AD]);
            //            }
            
            //-->FullTimes
            if(row==2)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    FullTime[AD+(week-1)*MaxDetectors]=readvals[AD];
                }
            }
            //-->muon efficiencies
            if(row==3)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    MuonEff[AD+(week-1)*MaxDetectors]=readvals[AD];
                }
            }
            //-->dmc efficiencies
            if(row==4)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    MultiEff[AD+(week-1)*MaxDetectors]=readvals[AD];
                }
            }
            //-->target masses
            if(row==8)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    DetectorMassGdLs[AD]=readvals[AD]/1000000;//ktons
                }
            }
            //-->bg events
            if(row==9)
            {
                //                for(Int_t AD=0;AD<NADs;AD++)
                //                {
                //                    BgEvts[AD+(week-1)*MaxDetectors]=readvals[AD];
                //                    BgEvts[AD+(week-1)*MaxDetectors]*=MultiEff[AD+(week-1)*MaxDetectors]*MuonEff[AD+(week-1)*MaxDetectors];
                //                }
            }
            //-->bg systematics
            if(row==10) continue;//<-- we don't care about background systematics
            
            //-->acc bg
            if(row==11)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    AccidentalRate[AD+(week-1)*MaxDetectors]=readvals[AD];
                }
            }
            if(row==12)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    AccidentalError[AD+(week-1)*MaxDetectors]=readvals[AD];
                }
            }
            
            //-->li9 bg
            if(row==13)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    LiHeRate[AD+(week-1)*MaxDetectors]=readvals[AD];
                }
            }
            if(row==14)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    LiHeError[AD+(week-1)*MaxDetectors]=readvals[AD];
                }
            }
            //-->fast-n bg
            if(row==15)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    FastNeutronRate[AD+(week-1)*MaxDetectors]=readvals[AD];
                }
            }
            if(row==16)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    FastNeutronError[AD+(week-1)*MaxDetectors]=readvals[AD];
                }
            }
            //-->amc bg
            if(row==17)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    AmCRate[AD+(week-1)*MaxDetectors]=readvals[AD];
                }
            }
            if(row==18)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    AmCError[AD+(week-1)*MaxDetectors]=readvals[AD];
                }
            }
            //-->aln bg
            //            if(row==19)
            //            {
            //                for(Int_t AD=0;AD<NADs;AD++)
            //                {
            //                    AlnEvts[AD+(week-1)*MaxDetectors]=readvals[AD];
            //                }
            //            }
            //            if(row==20)
            //            {
            //                for(Int_t AD=0;AD<NADs;AD++)
            //                {
            //                    AlnErr[AD+(week-1)*MaxDetectors]=readvals[AD];
            //                }
            //            }
        }
        linenum++;//only lines >2
    }
}

void NominalData :: ReadChristineReactorSpectrum()
{
    std::cout << " READING CHRISTINE SPECTRUM " << std::endl;
    
    std::ifstream fileData("./ReactorInputs/p12c_blinded/combined/nNu_Nom_combined.txt");
    
    std::string lines;
    Double_t eMin=0;
    Double_t eMax=0;
    Double_t eNu=0;
    Double_t dNdE[NReactors];
    
    for (Int_t reactor = 0; reactor < NReactors; reactor++)
    {
        dNdE[reactor]=0;
    }
    
    Int_t curSample = 0;
    binWidth = 0;
    
    while(!fileData.eof())//218
    {
        fileData >> eNu;
        
        for (Int_t i = 0; i < NReactors; i++)
        {
            fileData >> dNdE[i];
        }
        
        if(curSample==0)
        {
            eMin=eNu;
            eMax=eNu;
        }
        if(eNu>eMax)
        {
            eMax=eNu;
        }
        if(curSample==1)
        {
            binWidth = eMax-eMin;
        }
        for (Int_t reactor = 0; reactor < NReactors; reactor++)
        {
            m_dNdE_nom[reactor+NReactors*curSample] = dNdE[reactor];
        }
        
        curSample++;
    }
    
    curSample--;//219, from 1.8 to 12.75 in 0.05 steps.
    
    fileData.close();
    
    m_eMin = eMin;
    m_eMax = eMax;
    m_nSamples = curSample;
}

void NominalData :: ReadChristineCovMatrix()
{
    std::cout << " READING CHRISTINE COVARIANCE MATRIX " << std::endl;

    // Read covarianvematrix
    std::ifstream fileData_mcov("./ReactorInputs/p12c_blinded/combined/nNu_Mcov_combined.txt");
    
    for (Int_t i = 0; i < m_nSamples*NReactors; i++)
    {
        for (Int_t j = 0; j < m_nSamples*NReactors; j++)
        {
            fileData_mcov >> m_dNdE_mcov[i][j];
            
            //                    std::cout << "COV MATRIX" << m_dNdE_mcov[i][j] << std::endl;
        }
    }
    
    fileData_mcov.close();
    
    TMatrixD covmatrix(NReactors * MatrixBins, NReactors * MatrixBins, &m_dNdE_mcov[0][0]);//   Fix dimensions, then we resize it to avoid 0's in the empty spaces.
    covmatrix.ResizeTo(m_nSamples * NReactors, m_nSamples * NReactors);
//    if(Print)
//    {
//        TCanvas* c2 = new TCanvas("","");
//        covmatrix.Draw("colz");
//        c2->Print("./Images/Reactor/ReactorCovMatrix.eps");
//        delete c2;
//    }
    TDecompChol chol(covmatrix);//  M = L*U
    chol.Decompose();
    TMatrixD cmat(chol.GetU());//   U
    TMatrixD tcmat(cmat.Transpose(cmat));// L
    
    Double_t* tmp_matrix = tcmat.GetMatrixArray();
    
    for (Int_t i = 0; i < m_nSamples * NReactors; i++)
    {
        for (Int_t j = 0; j <  m_nSamples * NReactors; j++)
        {
            L[j+i*m_nSamples*NReactors] = tmp_matrix[i * m_nSamples * NReactors + j];
        }
    }
    covmatrix.~TMatrixD();
    cmat.~TMatrixD();
    tcmat.~TMatrixD();
}

Double_t NominalData :: GetReactorCovMatrix(Int_t i,Int_t j)
{
    return L[j+i*m_nSamples*NReactors];
}

Double_t NominalData :: GetNominalReactorSpectrum(Int_t reactor, Int_t sample)
{
    return m_dNdE_nom[reactor+NReactors*sample];
}

Int_t NominalData :: GetReactorSamples()
{
    return m_nSamples;
}

Double_t NominalData :: GetReactorBinWidth()
{
    return binWidth;
}

Double_t NominalData :: GetReactorEmin()
{
    return m_eMin;
}

Double_t NominalData :: GetReactorEmax()
{
    return m_eMax;
}

bool NominalData :: GetUseToyMCTree()
{
    return UseToyMCTree;
}

void NominalData :: SetUseToyMCTree(bool UseToyMCTreeb)
{
    UseToyMCTree = UseToyMCTreeb;
}


void NominalData :: SetAllRandomSystematics(bool set)
{
    if(set==1)
    {
        StatisticalFluctuation = 1;
        
        VaryAccidentalMatrix= 1;
        VaryLiHeMatrix= 1;
        VaryFastNeutronsMatrix= 1;
        VaryAmCMatrix= 1;
        DistortLiHeMatrix= 1;
        DistortFastNeutronsMatrix= 1;
        DistortAmCMatrix= 1;
        IsotopeMatrix= 1;
        ReactorPowerMatrix= 1;
        RelativeEnergyOffsetMatrix= 1;
        AbsoluteEnergyOffsetMatrix= 1;
        AbsoluteEnergyScaleMatrix= 1;
        RelativeEnergyScaleMatrix= 1;
        IAVMatrix= 1;
        NLMatrix= 1;
        ResolutionMatrix= 1;
        Sin22t12Matrix= 1;
        EfficiencyMatrix= 1;
    }
    else
    {
        StatisticalFluctuation = 0;
        
        VaryAccidentalMatrix= 0;
        VaryLiHeMatrix= 0;
        VaryFastNeutronsMatrix= 0;
        VaryAmCMatrix= 0;
        DistortLiHeMatrix= 0;
        DistortFastNeutronsMatrix= 0;
        DistortAmCMatrix= 0;
        IsotopeMatrix= 0;
        ReactorPowerMatrix= 0;
        RelativeEnergyOffsetMatrix= 0;
        AbsoluteEnergyOffsetMatrix= 0;
        AbsoluteEnergyScaleMatrix= 0;
        RelativeEnergyScaleMatrix= 0;
        IAVMatrix= 0;
        NLMatrix= 0;
        ResolutionMatrix= 0;
        Sin22t12Matrix= 0;
        EfficiencyMatrix= 0;
    }
}

void NominalData :: SetResponseDirectory(std::string directory)
{
    ResponseDirectory = directory;
}

std::string NominalData :: GetResponseDirectory()
{
    return ResponseDirectory;
}

void NominalData :: SetPredictionDirectory(std::string directory)
{
    PredictionsDirectory = directory;
}

std::string NominalData :: GetPredictionDirectory()
{
    return PredictionsDirectory;
}

void NominalData :: SetToyMCSamplesDirectory(std::string directory)
{
     ToyMCSamplesDirectory = directory;
}

std::string NominalData :: GetToyMCSamplesDirectory()
{
    return ToyMCSamplesDirectory;
}

void NominalData :: SetBkgCovDirectory(std::string directory)
{
    BkgCovDirectory = directory;
}

std::string NominalData :: GetBkgCovDirectory()
{
    return BkgCovDirectory;
}

void NominalData :: SetSysCovDirectory(std::string directory)
{
    SysCovDirectory = directory;
}

std::string NominalData :: GetSysCovDirectory()
{
    return SysCovDirectory;
}

void NominalData :: SetHierarchy(Int_t Hierarchy)
{
    hierarchy = Hierarchy;
}

Int_t NominalData :: GetHierarchy()
{
    return hierarchy;
}

Double_t NominalData :: GetSinStart()
{
    return sin_start;
}

Double_t NominalData :: GetSinEnd()
{
    return sin_end;
}

Double_t NominalData :: GetDmeeStart()
{
    return dmee_start;
}

Double_t NominalData :: GetDmeeEnd()
{
    return dmee_end;
}

Double_t NominalData :: GetObservedEvents(Int_t i,Int_t j)
{
    return ObservedEvents[i+j*MaxDetectors];
}

Int_t NominalData :: GetNReactorPeriods()
{
    return NReactorPeriods;
}

void NominalData :: SetSinStart(Double_t sinstart)
{
    sin_start = sinstart;
}
void NominalData :: SetSinEnd(Double_t sinend)
{
    sin_end = sinend;
}
void NominalData :: SetDmeeStart(Double_t dmeestart)
{
    dmee_start = dmeestart;
}
void NominalData :: SetDmeeEnd(Double_t dmeeend)
{
    dmee_end = dmeeend;
}
void NominalData :: SetNReactorPeriods(Int_t reactorperiods)
{
    NReactorPeriods = reactorperiods;
}

void NominalData :: CalculateBinning()
{
    
    //Gadolinium:

    if(LinearBinning)//  Linear binning
    {
        if(isH)//Hydrogen analysis linear binning
        {
            n_evis_bins=42;//match visible cuts
            n_etrue_bins=39;//match reactor data, this might be different when using IHEP reactor model for P14
           
            EVisMin = 1.5;//1.5 lower visible energy cut for Hydrogen Analysis

            if(LoganBinning)//To test our different spectra
            {
                EVisMin = 0;//1.5 for Hydrogen Analysis
                n_evis_bins = 240;
            }
        }
        else//Linear binning for Gadolinium Analysis
        {
            n_evis_bins = 60;//0.2 steps from 0 to 12
            n_etrue_bins=39;//match reactor flux data for P12E
            EVisMin = 0;//show all spectrum
        }
        
        Double_t TrueBinWidth = (Emax - Emin)/n_etrue_bins;//if n_etrue_bins = 51 and Emin = 1.8, 0.2
        Double_t VisBinWidth = (EVisMax - EVisMin)/n_evis_bins;//if 240 and Emin = 0 it is 0.05

        //Linear binning
        for (Int_t i = 0; i <= n_evis_bins; i++)
        {
            evis_bins[i] = VisBinWidth * i + EVisMin;
            enu_bins[i] = TrueBinWidth * i + Emin;
        }
    }
    else
    {
        if(!isH)// LBNL Non-linear binning for Gadolinium Analysis
        {
            n_evis_bins=37;
            n_etrue_bins=39;
            
            for (Int_t i = 0; i <= n_etrue_bins; i++)
            {
                enu_bins[i] = 0.2 * i + Emin;
            }
            
            evis_bins[0] = 0.7;
            
            for (Int_t i = 0; i < n_evis_bins-1; i++)
            {
                evis_bins[i+1] = 0.2 * i + 1.0;
            }
            evis_bins[n_evis_bins] = 12;
        }
    }
}

Double_t NominalData :: GetTrueBinningArray(Int_t true_index)
{
    return enu_bins[true_index];
}
Double_t NominalData :: GetVisibleBinningArray(Int_t vis_index)
{
    return evis_bins[vis_index];
}
Int_t NominalData :: GetTrueBins()
{
    return n_etrue_bins;
}
Int_t NominalData :: GetVisibleBins()
{
    return n_evis_bins;
}
