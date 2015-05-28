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

// Used to debug, if you don't want it to run, uncomment the following line: (#define NDEBUG)
//#define NDEBUG

#include <assert.h>

#define UseVolumes //To use 2 volumes, otherwise 100 cells.
//#define EREC_COMPARISON // To calculate the response matrices using ERec, and check the fit result of theta 13 using ERec data and prediction.

#ifdef UseVolumes//By default, this code could also work with cells, but they are statistically limited, and it is only useful for MC purposes.

//From Logan : These are relative to the number of protons in GdLS, so they do not sum to 1; we simply need to normalize them to 1.  Basically, LS efficiency is 0.539/0.146 times larger than GdLS.
    const Double_t LS_volume_efficiency = 0.539/(0.539+0.146);
    const Double_t GdLS_volume_efficiency = 0.146/(0.539+0.146);
#else
    const Double_t LS_volume_efficiency = 1;
    const Double_t GdLS_volume_efficiency = 1;
#endif

#define PrintEps//To save results in .eps files
//#define BlindedAnalysis // To use blinded reactor model and distances. Not all the files are operative (only those coming from Christine's reactor model)
const bool DeltaMee = 0;//Use Δm^2ee instead of Δm32 and Δm31 values

const bool ADSimple = 1;

const bool LoganBinning = 1;//240 visible bins from 0 to 12 MeV to check coherence between our predicted  spectra

const Int_t MaxExperiments = 1000;
const Int_t MaxPeriods = 101;//To be increased in further Data productions
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

/// For the nH Toy MC:
const Double_t IBDthreshold = 1.80607;  //( (Mn+Me)^2-Mp^2 ) / (2Mp) =
// = ( (939.565378+0.510998928)^2 - 938.272046^2 ) / (2*938.272046)
//1.80433;  // Mneutron+Mpositron-Mproton
const Double_t EnergyScale = 0.982;  // data_centercell/toy_centercell = 0.982 for nH gamma toy non-uniformity and AdSimple data
//const double EnergyScale = 1.023;  // data_centercell/toy_centercell = 1.02 when using same non-uniformity in both data and toy
const Int_t MaxCellNum = 401;

/// vetex region
const Int_t R2_binnum = 10;                      // ---> set option: divide the volume to sub-regions
const Double_t R2_lower  = 0;// m2                  // ---> set option
const Double_t R2_upper  = 4;// m2                  // ---> set option
const Double_t R2_binwidth = 4.0/R2_binnum;                  // ---> set option

const Int_t Z_binnum  = 10;                      // ---> set option
const Double_t Z_lower   = -2;// m                  // ---> set option
const Double_t Z_upper   = 2;// m                   // ---> set option

#ifdef UseVolumes
const Int_t VolumeX = 2;                      // ---> set option: divide the volume to sub-regions (2, GdLs-Ls) (10 cells = R2_binnum)

const Double_t VolumeX_lower  = 0;// m2                  // ---> set option
const Double_t VolumeX_upper  = 2;// m2                  // ---> set option

const Int_t VolumeY = 1;                      // ---> set option: divide the volume to sub-regions (1, Ls) (10 cells)

const Double_t VolumeY_lower  = 0;// m2                  // ---> set option
const Double_t VolumeY_upper  = 1;// m2                  // ---> set option
#else
const Int_t VolumeX = R2_binnum;                      // ---> set option: divide the volume to sub-regions (2, GdLs-Ls) (10 cells = R2_binnum)

const Double_t VolumeX_lower  = 0;// m2                  // ---> set option
const Double_t VolumeX_upper  = 4;// m2                  // ---> set option

const Int_t VolumeY = Z_binnum;                      // ---> set option: divide the volume to sub-regions (2, GdLs-Ls) (10 cells = R2_binnum)

const Double_t VolumeY_lower  = -2;// m2                  // ---> set option
const Double_t VolumeY_upper  = 2;// m2                  // ---> set option
#endif

class NominalData
{
private:
    
    std::string ResponseDirectory;
    std::string PredictionsDirectory;
    std::string ToyMCSamplesDirectory;
    std::string SysCovDirectory;
    std::string BkgCovDirectory;
    Char_t DistanceFileName[100];

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
    
    Double_t ProtonsPerKtonGdLs; //protons per kton GdLS
    Double_t ProtonsPerKtonLs; //protons per kton LS volume
    Double_t DetectorMass[MaxDetectors*VolumeX];
    Double_t nHCaptureFraction[VolumeX];
    
    Double_t m_detectorEfficiency_Dt;
    Double_t m_detectorEfficiency_Ep;
    Double_t m_detectorEfficiency_Ed_nominal;
    Double_t m_detectorEfficiency_flash;
    Double_t m_detectorEfficiency_nGd;
    Double_t m_detectorEfficiency_spill;
    
    Double_t m_detectorGlobalEfficiency[VolumeX];

    //nH Efficiency map
//    TH2D *h2d_Ep_ratio2center[MaxDetectors*MaxPeriods];
    
    Double_t ADIntegral[MaxDetectors];
    Double_t CellIntegral[MaxDetectors*VolumeX*VolumeY];
    Double_t PercentualEvents[MaxDetectors*VolumeX*VolumeY];
    
    bool delete_nH_maps_flag;
    
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
    Int_t DataPeriods;
    
    //Reactor parameters
    Double_t IsotopeFrac[NIsotopes];
    Double_t IsotopeFracError[NIsotopes];
    
    Double_t ReactorPower[NReactors];
    Double_t ReactorPowerError[NReactors];
    
    Double_t EnergyPerFission[NReactors];
    Double_t EnergyPerFissionError[NReactors];
    
    //Background relative errors, rates and events:
    Double_t AccidentalError[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    Double_t LiHeError[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    Double_t FastNeutronError[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    Double_t AmCError[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    
    Double_t AccidentalRate[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    Double_t FastNeutronRate[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    Double_t LiHeRate[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    Double_t AmCRate[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    
    Double_t AccidentalEvents[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    Double_t FastNeutronEvents[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    Double_t LiHeEvents[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    Double_t AmCEvents[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    
    Double_t ObservedEvents[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    Double_t IBDEvents[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    Double_t StatisticalError[MaxDetectors*MaxPeriods*VolumeX*VolumeY];
    // Days and efficiencies:
    Double_t FullTime[MaxDetectors*MaxPeriods];//Full time is the same for both volumes
    Double_t MuonEff[MaxDetectors*MaxPeriods];
    Double_t MultiEff[MaxDetectors*MaxPeriods];
    Double_t InclusiveFullTime[MaxDetectors*MaxPeriods];//Full time is the same for both volumes
    Double_t InclusiveMuonEff[MaxDetectors*MaxPeriods];
    Double_t InclusiveMultiEff[MaxDetectors*MaxPeriods];
    
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
    Double_t m_rel_escale;
    Double_t m_rel_escale_error;
    Double_t m_rel_escale_nominal;
    Double_t m_rel_eoffset;
    
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
    bool OAVMatrix;
    bool NLMatrix;
    bool ResolutionMatrix;
    bool Sin22t12Matrix;
    bool EfficiencyMatrix;
    bool AllDetectorSystematicsMatrix;
    
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
    bool OAVBudget;
    bool NLBudget;
    bool ResolutionBudget;
    bool Sin22t12Budget;
    bool EfficiencyBudget;
    bool SystematicBudget;
    bool BackgroundBudget;
    bool TotalBudget;
    
    bool StatisticalFluctuation;
    bool UseToyMCTree;
    
    Double_t ADdistances[MaxDetectors*NReactors];
    
    NominalData(bool,Int_t);
    void CopyData(NominalData*);
    ~NominalData();
    
    void ReadChristineCovMatrix();//P12C
    void ReadChristineReactorSpectrum();
    
    void ReadIHEPReactorSpectrum();//P14A
    bool IHEPReactorModel;
    bool UsingIHEPReactorModel();
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
    
    void LoadOriginalGDMainData(const Char_t*);//LBNL inputs
    void LoadHydrogenMainData();//To correct for efficiencies using Gd/H data info from a txt file.
    
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
    void SetOAVMatrix(bool);
    void SetIAVMatrix(bool);
    void SetNLMatrix(bool);
    void SetResolutionMatrix(bool);
    void SetSin22t12Matrix(bool);
    void SetEfficiencyMatrix(bool);
    void SetAllDetectorSystematicsMatrix(bool);
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
    void SetOAVBudget(bool);
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
    void SetDataPeriods(Int_t);
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
    bool GetOAVMatrix();
    bool GetNLMatrix();
    bool GetResolutionMatrix();
    bool GetSin22t12Matrix();
    bool GetEfficiencyMatrix();
    bool GetAllDetectorSystematicsMatrix();
    
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
    bool GetOAVBudget();
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
    
    Double_t GetRelativeEnergyScale();
    Double_t GetRelativeEnergyError();
    Double_t GetRelativeEnergyOffset();
    Double_t GetRelativeEnergyOffsetError();
    
    Double_t GetResolutionError();
    Double_t GetResoUncorrelatedError();
    
    Double_t GetDetectorProtons(Int_t,Int_t);
    
    Double_t GetEnergyDelayedCutDetectorEfficiency();
    Double_t GetDetectorEfficiencyRelativeError();
    Double_t GetDetectorEfficiency(Int_t, Int_t, Int_t, Int_t,bool);
    Double_t GetFullTime(Int_t,Int_t,bool);
    
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
    Double_t GetAccidentalError(Int_t,Int_t,Int_t,Int_t);
    Double_t GetLiHeError(Int_t,Int_t,Int_t,Int_t);
    Double_t GetFNError(Int_t,Int_t,Int_t,Int_t);
    Double_t GetAmCError(Int_t,Int_t,Int_t,Int_t);
    
    Double_t GetIBDEvents(Int_t, Int_t,Int_t, Int_t);
    Double_t GetObservedEvents(Int_t,Int_t,Int_t, Int_t);
    
    Double_t GetAccidentalRate(Int_t,Int_t,Int_t,Int_t);
    Double_t GetLiHeRate(Int_t,Int_t,Int_t,Int_t);
    Double_t GetFNRate(Int_t,Int_t,Int_t,Int_t);
    Double_t GetAmCRate(Int_t,Int_t,Int_t,Int_t);
    
    Double_t GetAccidentalEvents(Int_t,Int_t,Int_t,Int_t);
    Double_t GetLiHeEvents(Int_t,Int_t,Int_t,Int_t);
    Double_t GetFNEvents(Int_t,Int_t,Int_t,Int_t);
    Double_t GetAmCEvents(Int_t,Int_t,Int_t,Int_t);
    
    Int_t GetNReactorPeriods();
    Int_t GetDataPeriods();
    
    void GetnHInclusiveData(Double_t*);
    void GetnGdInclusiveData(Double_t*);
    void CorrectnHEvents();
    void CorrectnGdEvents();
    void CalculateInclusiveFullTime(Double_t*,Double_t*);
    void CalculateInclusiveEfficiencies(Double_t*,Double_t*);
    void ReadToyEventsByCell();
    Double_t GetEventsByCell(Int_t,Int_t,Int_t);
    void ReadDistances(Char_t*);
    Double_t GetDistances(Int_t,Int_t);
};

NominalData :: NominalData(bool ish,Int_t dataSet)
{
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //                              ANALYSIS, DATA SET AND NON-LINEARITY MODEL TO BE USED
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    DataSet = dataSet;
    isH = ish;
    
    delete_nH_maps_flag = 0;
    
    BCW=0;
    LBNL=0;
    Unified=1;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //                                          BINNING VARIABLES AND FITTER PARAMETERS
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
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
    OAVBudget=0;
    NLBudget=0;
    ResolutionBudget=0;
    Sin22t12Budget=0;
    EfficiencyBudget=0;
    
    StatisticalFluctuation=0;

    NSamples = 500;
    NSteps = 101;
    
    Emin = 1.8;
    Emax = 9.6;
    EVisMax = 12;
    LinearBinning=0;
    //CalculateBinning();//Calculated in setbinning
    
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
    
    ToyMC=0;//Data is default;
    Nweeks = 1;
    NReactorPeriods=20;
    DataPeriods = 101;
    NADs = 6;
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //                                          LOAD DISTANCES BETWEEN ADS AND CORES
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef BlindedAnalysis
    sprintf(DistanceFileName,"./Distances/blinded_baseline.txt");
    std::cout << " NEED TO ADD A BLINDED DISTANCE FILE" << std::endl;
    exit(EXIT_FAILURE);
#else
    sprintf(DistanceFileName,"./Distances/unblinded_baseline.txt");
#endif
    
    if(NADs == 8) //ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
#ifdef BlindedAnalysis
        sprintf(DistanceFileName,"./Distances/blinded_baseline.txt");
        std::cout << " NEED TO ADD A BLINDED 8AD DISTANCE FILE" << std::endl;
        exit(EXIT_FAILURE);
#else
        sprintf(DistanceFileName,"./Distances/unblinded_baseline8ADs.txt");//change file to calculate distances to a file that has the information for the 8 ADs
        std::cout << "NEED TO ADD A TXT FILE WITH THE 8AD DISTANCES" << std::endl;
        exit(EXIT_FAILURE);
#endif
    }
    
    ReadDistances(DistanceFileName);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //                                          LOAD DETECTOR MASS
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef TestAllTheSame
    for(Int_t i = 0; i<NADs; i++)
    {
        for(Int_t j = 0; j<VolumeX; j++)
        {
            DetectorMass[i+MaxDetectors*j]=0.02;
        }
    }
    
    ProtonsPerKtonGdLs = 7e31;
    ProtonsPerKtonLs = 7e31;
#else
    ProtonsPerKtonGdLs = 7.1638e31; //protons per kton GdLs
    ProtonsPerKtonLs = 7.116e31; //protons per kton Ls

    //Gd-LS Volume:
    DetectorMass[0] = 0.0215735;// kton
    DetectorMass[1] = 0.0215196;
    DetectorMass[2] = 0.0215872;
    DetectorMass[3] = 0.0215662;
    DetectorMass[4] = 0.0214088;
    DetectorMass[5] = 0.0216526;
    //LS volume:
    DetectorMass[MaxDetectors] = 0.019941;// kton
    DetectorMass[MaxDetectors+1] = 0.019966;
    DetectorMass[MaxDetectors+2] = 0.019891;
    DetectorMass[MaxDetectors+3] = 0.019913;
    DetectorMass[MaxDetectors+4] = 0.019991;
    DetectorMass[MaxDetectors+5] = 0.019892;
#endif
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //                                      REACTOR PARAMETERS
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    IHEPReactorModel=0;

    //Nominal isotope fraction errors:
    for(Int_t i=0;i<NIsotopes;i++)
    {
        IsotopeFracError[i]=0.05;
    }
    
    //test:
#ifdef TestAllTheSame
    
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
    
    //Energy per fission
    EnergyPerFission[0]=200;//MeV/fission
    EnergyPerFission[1]=200;
    EnergyPerFission[2]=200;
    EnergyPerFission[3]=200;
    
#else
    
    //Reactor power (It's an average of first 140 days of data) This should be coded to get the data from a txt file instead.)

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
    
    //Energy per fission
    EnergyPerFission[0]=201.92;//MeV/fission
    EnergyPerFission[1]=205.52;
    EnergyPerFission[2]=209.99;
    EnergyPerFission[3]=213.60;
    
#endif

    //Reactor power errors:
    for(Int_t r=0;r<NReactors;r++)
    {
        ReactorPowerError[r]=0.005;
    }
        
    //Energy per fission errors
    EnergyPerFissionError[0]=0.002278;
    EnergyPerFissionError[1]=0.004671;
    EnergyPerFissionError[2]=0.002857;
    EnergyPerFissionError[3]=0.003043;
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //                                          DETECTOR RESPONSE PARAMETERS
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    //For Hydrogen and 6 ADs: from http://dayabay.ihep.ac.cn/DocDB/0085/008556/020/Main-version3.pdf
    if(this->GetAnalysis())
    {
        //IAV error
        IAVError=0.001;//0.1% bin-to-bin uncorrelated error.
        
        //Resolution errors
        ResolutionError=0.02;//Due to parameter uncertainty
        ResolutionErrorUncorrelated=0.02;//Due to energy scale difference between detectors
        
        //Attenuation length -> Relative Energy Scale
        m_rel_escale_error = 0.0065; // 0.65%
        
        //Hard coded values, now we use the txt inputs so the following lines should be meaningless:
        if(DataSet==1)
        {
            std::cout << "USING P14E values right now" << std::endl;
            
            exit(EXIT_FAILURE);
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
                
                MultiEff[0]= 0.9917;
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
                std::cout << "If weekly data is provided you need to weight the weekly efficiencies, see GetDataInclusive() " << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        if(DataSet==2)//P12E
        {
//            //p12b values
//            if(Nweeks==1)
//            {
//                LoadHydrogenMainData(Form("./Inputs/HInputs/nH_GdLS_table_Inclusive.txt"));
//                LoadHydrogenMainData(Form("./Inputs/HInputs/nH_LS_table_Inclusive.txt"));
//            }
//            else
//            {
            
            
//                LoadHydrogenMainData();// Need to do it NReactorPeriods dependent instead of just 101. DataSet chooses the (101)NReactorPeriods so far.
//                CorrectnHEvents();
            
            
            
//            }
        }
    }
    //For Gd and 6 ADs: http://dayabay.ihep.ac.cn/DocDB/0085/008556/016/Main-version2.pdf page 5
    else
    {
        //IAV thickness relative error
        IAVError=0.04;
        
        //Energy Scale
        m_abs_escale =1;
        m_abs_escale_error = 0.01;
        
        m_abs_eoffset = 0.0;
        m_abs_eoffset_error = 0.08;  //(MeV)
        
        m_rel_eoffset_error = 0.013; //  (MeV)
        
        m_rel_escale_error = 0.0035; // 0.35%
        
        m_rel_escale = 1.0;
        m_rel_escale_nominal = m_rel_escale;
        m_rel_eoffset = 0.0;
        
        //Efficiency errors:
        m_detectorEfficiency_Dt = 0.986;
        m_detectorEfficiency_Ep = 0.9988;
        m_detectorEfficiency_Ed_nominal = 0.909;
        m_detectorEfficiency_flash = 0.9999;//No 0.9998 because Mineral Oil is not part of the fidutial volume.
        m_detectorEfficiency_nGd = 0.838;
        m_detectorEfficiency_spill = 1.050;
        
        m_detectorGlobalEfficiency[0] = m_detectorEfficiency_spill*m_detectorEfficiency_nGd*m_detectorEfficiency_Dt*m_detectorEfficiency_Ep*m_detectorEfficiency_Ed_nominal*m_detectorEfficiency_flash;
      
        DetectorEfficiencyRelativeError = 0.0014;
        
        //Resolution errors
        ResolutionError=0.002;//Due to parameter uncertainty
        ResolutionErrorUncorrelated=0.002;//Due to energy scale difference between detectors
        
        if(DataSet==2)
        {
            // P12E Values in Theta13-inputs_32week_inclusive.txt
            // These are hard coded values taken from the txt file, use LoadnGdMainData to read from the txt file:
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
#ifdef TestAllTheSame
                    for(Int_t i = 0; i<NADs; i++)
                    {
                        for(Int_t j = 0; j<MaxPeriods; j++)
                        {
                            ObservedEvents[i+j*NADs]=100000;
                            MuonEff[i+j*NADs]=1;
                            MultiEff[i+j*NADs]=1;
                            FullTime[i+j*NADs]=100;
                            AccidentalRate[i+j*NADs]=1;
                            FastNeutronRate[i+j*NADs]=1;
                            LiHeRate[i+j*NADs]=1;
                            AmCRate[i+j*NADs]=1;
                        }
                    }
#endif
                
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
//                LoadOriginalGDMainData(Form("./Inputs/GdInputs/Theta13-inputs_20week.txt"));// Need to do it Nweek dependent
            }
        }
    }
}

NominalData:: ~NominalData()
{
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
    
    ProtonsPerKtonGdLs = data->ProtonsPerKtonGdLs; //protons per kton
    ProtonsPerKtonLs = data->ProtonsPerKtonLs; //protons per kton
    std::copy(std::begin(data->DetectorMass), std::end(data->DetectorMass), std::begin(DetectorMass));
    
    m_detectorEfficiency_Dt = data->m_detectorEfficiency_Dt;
    m_detectorEfficiency_Ep = data->m_detectorEfficiency_Ep;
    m_detectorEfficiency_Ed_nominal = data->m_detectorEfficiency_Ed_nominal;
    m_detectorEfficiency_flash = data->m_detectorEfficiency_flash;
    m_detectorEfficiency_nGd = data->m_detectorEfficiency_nGd;
    m_detectorEfficiency_spill = data->m_detectorEfficiency_spill;
    std::copy(std::begin(data->m_detectorGlobalEfficiency), std::end(data->m_detectorGlobalEfficiency), std::begin(m_detectorGlobalEfficiency));
    std::copy(std::begin(data->nHCaptureFraction), std::end(data->nHCaptureFraction), std::begin(nHCaptureFraction));
    
    std::copy(std::begin(data->ADdistances), std::end(data->ADdistances), std::begin(ADdistances));
    std::copy(std::begin(data->DistanceFileName), std::end(data->DistanceFileName), std::begin(DistanceFileName));

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
    
    //Reactor parameters
    std::copy(std::begin(data->IsotopeFrac), std::end(data->IsotopeFrac), std::begin(IsotopeFrac));
    std::copy(std::begin(data->IsotopeFracError), std::end(data->IsotopeFracError), std::begin(IsotopeFracError));
    std::copy(std::begin(data->ReactorPower), std::end(data->ReactorPower), std::begin(ReactorPower));
    std::copy(std::begin(data->ReactorPowerError), std::end(data->ReactorPowerError), std::begin(ReactorPowerError));
    std::copy(std::begin(data->EnergyPerFission), std::end(data->EnergyPerFission), std::begin(EnergyPerFission));
    std::copy(std::begin(data->EnergyPerFissionError), std::end(data->EnergyPerFissionError), std::begin(EnergyPerFissionError));
    IHEPReactorModel = data->IHEPReactorModel;
    
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
    std::copy(std::begin(data->IBDEvents), std::end(data->IBDEvents), std::begin(IBDEvents));
    std::copy(std::begin(data->StatisticalError), std::end(data->StatisticalError), std::begin(StatisticalError));

    std::copy(std::begin(data->ADIntegral), std::end(data->ADIntegral), std::begin(ADIntegral));
    std::copy(std::begin(data->CellIntegral), std::end(data->CellIntegral), std::begin(CellIntegral));
    std::copy(std::begin(data->PercentualEvents), std::end(data->PercentualEvents), std::begin(PercentualEvents));

    // Days and efficiencies:
    std::copy(std::begin(data->FullTime), std::end(data->FullTime), std::begin(FullTime));
    std::copy(std::begin(data->MuonEff), std::end(data->MuonEff), std::begin(MuonEff));
    std::copy(std::begin(data->MultiEff), std::end(data->MultiEff), std::begin(MultiEff));
    std::copy(std::begin(data->InclusiveFullTime), std::end(data->InclusiveFullTime), std::begin(InclusiveFullTime));
    std::copy(std::begin(data->InclusiveMuonEff), std::end(data->InclusiveMuonEff), std::begin(InclusiveMuonEff));
    std::copy(std::begin(data->InclusiveMultiEff), std::end(data->InclusiveMultiEff), std::begin(InclusiveMultiEff));
    
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
    m_rel_escale = data->m_rel_escale;
    m_rel_escale_error = data->m_rel_escale_error;
    m_rel_escale_nominal = data->m_rel_escale_nominal;
    m_rel_eoffset = data->m_rel_eoffset;
    
    DetectorEfficiencyRelativeError = data->DetectorEfficiencyRelativeError;
    
    //Binning variables:
    Nweeks = data->Nweeks;
    NReactorPeriods = data->NReactorPeriods;
    DataPeriods = data->DataPeriods;
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
    OAVMatrix = data->OAVMatrix;
    NLMatrix = data->NLMatrix;
    ResolutionMatrix = data->ResolutionMatrix;
    Sin22t12Matrix = data->Sin22t12Matrix;
    EfficiencyMatrix = data->EfficiencyMatrix;
    AllDetectorSystematicsMatrix = data->AllDetectorSystematicsMatrix;
    
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
    OAVBudget = data->OAVMatrix;
    
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

void NominalData :: SetVaryAccidentalMatrix(bool varyAccMatrix)
{
    VaryAccidentalMatrix=varyAccMatrix;
}

void NominalData :: SetVaryLiHeMatrix(bool varyLiMatrix)
{
    VaryLiHeMatrix=varyLiMatrix;
}

void NominalData :: SetVaryFastNeutronsMatrix(bool varyFNmatrix)
{
    VaryFastNeutronsMatrix=varyFNmatrix;
}

void NominalData :: SetVaryAmCMatrix(bool varyAmCmatrix)
{
    VaryAmCMatrix=varyAmCmatrix;
}

void NominalData :: SetDistortLiHeMatrix(bool distortLiMatrix)
{
    DistortLiHeMatrix=distortLiMatrix;
}

void NominalData :: SetDistortFastNeutronsMatrix(bool distortFNmatrix)
{
    DistortFastNeutronsMatrix=distortFNmatrix;
}

void NominalData :: SetDistortAmCMatrix(bool distortAmCmatrix)
{
    DistortAmCMatrix=distortAmCmatrix;
}

void NominalData :: SetIsotopeMatrix(bool isotopematrix)
{
    IsotopeMatrix=isotopematrix;
}

void NominalData :: SetReactorPowerMatrix(bool reactorPowermatrix)
{
    ReactorPowerMatrix=reactorPowermatrix;
}

void NominalData :: SetRelativeEnergyOffsetMatrix(bool relativeEnergyOffsetmatrix)
{
    RelativeEnergyOffsetMatrix = relativeEnergyOffsetmatrix;
}

void NominalData :: SetAbsoluteEnergyOffsetMatrix(bool absoluteEnergyOffsetmatrix)
{
    AbsoluteEnergyOffsetMatrix = absoluteEnergyOffsetmatrix;
}

void NominalData :: SetAbsoluteEnergyScaleMatrix(bool absoluteEnergyScalematrix)
{
    AbsoluteEnergyScaleMatrix=absoluteEnergyScalematrix;
}

void NominalData :: SetRelativeEnergyScaleMatrix(bool relativeEnergyScalematrix)
{
    RelativeEnergyScaleMatrix=relativeEnergyScalematrix;
}

void NominalData :: SetIAVMatrix(bool iAVmatrix)
{
    IAVMatrix=iAVmatrix;
}

void NominalData :: SetOAVMatrix(bool oAVmatrix)
{
    OAVMatrix=oAVmatrix;
}

void NominalData :: SetNLMatrix(bool nLmatrix)
{
    NLMatrix=nLmatrix;
}

void NominalData :: SetResolutionMatrix(bool resolutionmatrix)
{
    ResolutionMatrix=resolutionmatrix;
}

void NominalData :: SetSin22t12Matrix(bool sin22t12matrix)
{
    Sin22t12Matrix=sin22t12matrix;
}

void NominalData :: SetEfficiencyMatrix(bool efficiencymatrix)
{
    EfficiencyMatrix = efficiencymatrix;
}

void NominalData :: SetAllDetectorSystematicsMatrix(bool allDetectorSystematicsMatrix)
{
    AllDetectorSystematicsMatrix = allDetectorSystematicsMatrix;
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

void NominalData :: SetOAVBudget(bool oAVBudget)
{
    OAVBudget=oAVBudget;
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

bool NominalData :: GetOAVMatrix()
{
    return OAVMatrix;
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

bool NominalData :: GetAllDetectorSystematicsMatrix()
{
    return AllDetectorSystematicsMatrix;
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

bool NominalData :: GetOAVBudget()
{
    return OAVBudget;
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
    CalculateBinning();
}

void NominalData :: SetEmin(Double_t emin)
{
    Emin = emin;
    CalculateBinning();
}

void NominalData :: SetEmax(Double_t emax)
{
    Emax = emax;
    CalculateBinning();
}

void NominalData :: SetEVisMin(Double_t evismin)
{
    EVisMin = evismin;
    CalculateBinning();
}

void NominalData :: SetEVisMax(Double_t evismax)
{
    EVisMax = evismax;
    CalculateBinning();
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

Double_t NominalData :: GetRelativeEnergyScale()
{
    return m_rel_escale;
}

Double_t NominalData :: GetRelativeEnergyError()
{
    return m_rel_escale_error;
}

Double_t NominalData :: GetRelativeEnergyOffset()
{
    return m_rel_eoffset;
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

Double_t NominalData :: GetDetectorProtons(Int_t detector,Int_t idx)
{
    if(idx == 0)//GdLs
    {
        return ProtonsPerKtonGdLs*DetectorMass[detector+MaxDetectors*idx];
    }
    else
    {
        return ProtonsPerKtonLs*DetectorMass[detector+MaxDetectors*idx];
    }
}

Double_t NominalData :: GetDetectorEfficiency(Int_t detector, Int_t week, Int_t idx, Int_t idy,bool FluxMode)
{
    if(!isH)//Gd
    {
        if(FluxMode)//Weekly even if Nweeks == 1
        {
            return MuonEff[detector+week*MaxDetectors]*MultiEff[detector+week*MaxDetectors]*m_detectorGlobalEfficiency[idx];
        }
        
        if(Nweeks==1)//Inclusive
        {
            return InclusiveMuonEff[detector+week*MaxDetectors]*InclusiveMultiEff[detector+week*MaxDetectors]*m_detectorGlobalEfficiency[idx];
        }
        else
        {
            return MuonEff[detector+week*MaxDetectors]*MultiEff[detector+week*MaxDetectors]*m_detectorGlobalEfficiency[idx];
        }
    }
    else//nH
    {

        Double_t VolumeEfficiency = 1;
        
#ifdef UseVolumes//VolumeEfficiency should be equal to m_detectorGlobalEfficiency[idx]*nHCaptureFraction[idx];

        //There might be a small difference due to spill-in (accounted in Logan's calculation) which might be included in the delayed efficiency, needs to be checked, in the meantime use Xiangpan's values where
        
        //Absolute detector efficiency: MuonEff[detector+week*MaxDetectors]*MultiEff[detector+week*MaxDetectors]*m_detectorGlobalEfficiency[idx]*nHCaptureFraction[idx]
        if(idx==0)
        {
            VolumeEfficiency = GdLS_volume_efficiency;
        }//GdLs
        else
        {
            VolumeEfficiency = LS_volume_efficiency;
        }//Ls
#endif
        //return VolumeEfficiency * MuonEff[detector+week*MaxDetectors]*MultiEff[detector+week*MaxDetectors]; //using Logan's numbers.
        
//        std::cout <<  MuonEff[detector+week*MaxDetectors]*MultiEff[detector+week*MaxDetectors]*m_detectorGlobalEfficiency[idx]*nHCaptureFraction[idx] << std::endl;
        
        if(FluxMode)//Weekly even if Nweeks == 1
        {
            return MuonEff[detector+week*MaxDetectors]*MultiEff[detector+week*MaxDetectors]*m_detectorGlobalEfficiency[idx]*nHCaptureFraction[idx];
        }
        
        if(Nweeks==1)
        {
            return MuonEff[detector+week*MaxDetectors]*MultiEff[detector+week*MaxDetectors]*m_detectorGlobalEfficiency[idx]*nHCaptureFraction[idx];
        }
        else
        {
            return InclusiveMuonEff[detector+week*MaxDetectors]*InclusiveMultiEff[detector+week*MaxDetectors]*m_detectorGlobalEfficiency[idx]*nHCaptureFraction[idx];
        }
        //Efficiency of Prompt energy, Delayed energy, time, distance included in global efficiency. Flasher is not (Very small) nH capture included in delayed energy cut in the new version, right now working with out it included, that's why we multiply by the nH capture fraction.
    }
}

Double_t NominalData :: GetEnergyDelayedCutDetectorEfficiency()
{
    return m_detectorEfficiency_Ed_nominal;
}

Double_t NominalData :: GetDetectorEfficiencyRelativeError()
{
    return DetectorEfficiencyRelativeError;
}

Double_t NominalData :: GetFullTime(Int_t detector,Int_t week, bool FluxMode)
{
    if(FluxMode)
    {
        return FullTime[detector+week*MaxDetectors];
    }
    if(Nweeks==1)
    {
        return InclusiveFullTime[detector+week*MaxDetectors];
    }
    else
    {
        return FullTime[detector+week*MaxDetectors];
    }
}


Double_t NominalData :: GetIAVError()
{
    return IAVError;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                              Background inputs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Errors given as absolute error/(corrected rate) = relative error.
Double_t NominalData ::  GetAccidentalError(Int_t ad,Int_t week,Int_t idx, Int_t idy)
{
    return (AccidentalError[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods])/(AccidentalEvents[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods]);
}

Double_t NominalData :: GetLiHeError(Int_t ad,Int_t week,Int_t idx, Int_t idy)
{
    return (LiHeError[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods])/(LiHeEvents[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods]);
}

Double_t NominalData :: GetFNError(Int_t ad,Int_t week,Int_t idx, Int_t idy)
{
    return (FastNeutronError[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods])/(FastNeutronEvents[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods]);
}

Double_t NominalData :: GetAmCError(Int_t ad,Int_t week,Int_t idx, Int_t idy)
{
    return (AmCError[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods])/(AmCEvents[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods]);
}

Double_t NominalData :: GetAccidentalRate(Int_t ad,Int_t week,Int_t idx, Int_t idy)
{
    return AccidentalRate[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
}

Double_t NominalData :: GetLiHeRate(Int_t ad,Int_t week,Int_t idx, Int_t idy)
{
    return LiHeRate[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
}

Double_t NominalData :: GetFNRate(Int_t ad,Int_t week,Int_t idx, Int_t idy)
{
    return FastNeutronRate[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
}

Double_t NominalData :: GetAmCRate(Int_t ad,Int_t week,Int_t idx, Int_t idy)
{
    return AmCRate[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
}

Double_t NominalData :: GetIBDEvents(Int_t ad,Int_t week,Int_t idx, Int_t idy)//IBD events corrected for efficiencies
{
    return IBDEvents[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
}

Double_t NominalData ::  GetAccidentalEvents(Int_t ad,Int_t week,Int_t idx, Int_t idy)
{
     //Hydrogen data is given with the backgrounds rate with AD effects already included, Gd data has not efficiency corrections included anymore, but this was fixed in LoadOriginalGDMainData
    
    std::cout << " ACCIDENTAL EVENTS IN NOMINAL DATA : " << AccidentalEvents[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods] << std::endl;
    
    return AccidentalEvents[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
    
}

Double_t NominalData :: GetLiHeEvents(Int_t ad,Int_t week,Int_t idx, Int_t idy)
{
    //Hydrogen data is given with the backgrounds rate with AD effects already included, Gd data has not efficiency corrections included anymore, but this was fixed in LoadOriginalGDMainData

        return LiHeEvents[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
}

Double_t NominalData :: GetFNEvents(Int_t ad,Int_t week,Int_t idx, Int_t idy)
{
    //Hydrogen data is given with the backgrounds rate with AD effects already included, Gd data has not efficiency corrections included anymore, but this was fixed in LoadOriginalGDMainData

        return FastNeutronEvents[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
}

Double_t NominalData :: GetAmCEvents(Int_t ad,Int_t week,Int_t idx, Int_t idy)
{
    //Hydrogen data is given with the backgrounds rate with AD effects already included, Gd data has not efficiency corrections included anymore, but this was fixed in LoadOriginalGDMainData

        return AmCEvents[ad+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
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
void NominalData :: LoadOriginalGDMainData(const Char_t* mainmatrixname)
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
            DataPeriods=atoi(line.c_str());
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
                    ObservedEvents[AD+(week-1)*MaxDetectors]=readvals[AD];//GD VolumeX = 0, no need to add further indexes here
                    //Calculate Statistical Errors
                    StatisticalError[AD+(week-1)*MaxDetectors]=sqrt(ObservedEvents[AD+(week-1)*MaxDetectors]);
                }
            }
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
                    DetectorMass[AD]=readvals[AD]/1000000.;//ktons
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
                    AccidentalEvents[AD+(week-1)*MaxDetectors]=readvals[AD]*FullTime[AD+(week-1)*MaxDetectors]*MuonEff[AD+(week-1)*MaxDetectors]*MultiEff[AD+(week-1)*MaxDetectors];//Backgrounds are not corrected anymore in Theta13-inputs txt file, add here AD efficiencies
                }
            }
            if(row==12)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    AccidentalError[AD+(week-1)*MaxDetectors]=readvals[AD]*FullTime[AD+(week-1)*MaxDetectors]*MuonEff[AD+(week-1)*MaxDetectors]*MultiEff[AD+(week-1)*MaxDetectors];
                }
            }
            
            //-->li9 bg
            if(row==13)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    LiHeEvents[AD+(week-1)*MaxDetectors]=readvals[AD]*FullTime[AD+(week-1)*MaxDetectors]*MuonEff[AD+(week-1)*MaxDetectors]*MultiEff[AD+(week-1)*MaxDetectors];
                }
            }
            if(row==14)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    LiHeError[AD+(week-1)*MaxDetectors]=readvals[AD]*FullTime[AD+(week-1)*MaxDetectors]*MuonEff[AD+(week-1)*MaxDetectors]*MultiEff[AD+(week-1)*MaxDetectors];
                }
            }
            //-->fast-n bg
            if(row==15)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    FastNeutronEvents[AD+(week-1)*MaxDetectors]=readvals[AD]*FullTime[AD+(week-1)*MaxDetectors]*MuonEff[AD+(week-1)*MaxDetectors]*MultiEff[AD+(week-1)*MaxDetectors];
                }
            }
            if(row==16)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    FastNeutronError[AD+(week-1)*MaxDetectors]=readvals[AD]*FullTime[AD+(week-1)*MaxDetectors]*MuonEff[AD+(week-1)*MaxDetectors]*MultiEff[AD+(week-1)*MaxDetectors];
                }
            }
            //-->amc bg
            if(row==17)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    AmCEvents[AD+(week-1)*MaxDetectors]=readvals[AD]*FullTime[AD+(week-1)*MaxDetectors]*MuonEff[AD+(week-1)*MaxDetectors]*MultiEff[AD+(week-1)*MaxDetectors];
                }
            }
            if(row==18)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    AmCError[AD+(week-1)*MaxDetectors]=readvals[AD]*FullTime[AD+(week-1)*MaxDetectors]*MuonEff[AD+(week-1)*MaxDetectors]*MultiEff[AD+(week-1)*MaxDetectors];
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

    std::cout << " DataPeriods from file: " << DataPeriods << std::endl;
    
    for(Int_t AD=0;AD<NADs;AD++)
    {
        for(Int_t week = 0; week<DataPeriods; week++)
        {
            IBDEvents[AD+week*MaxDetectors] = (ObservedEvents[AD+week*MaxDetectors]-(AccidentalEvents[AD+week*MaxDetectors]+LiHeEvents[AD+week*MaxDetectors]+FastNeutronEvents[AD+week*MaxDetectors]+AmCEvents[AD+week*MaxDetectors]));//Apply background  suppression for events with AD effects, then correct AD effects and scale to days
        }
    }
    
    //Copy values
    for(Int_t AD=0;AD<NADs;AD++)
    {
        for(Int_t week = 0; week<DataPeriods; week++)
        {
            Double_t correction_factor = (FullTime[AD+week*MaxDetectors]*MuonEff[AD+week*MaxDetectors]*MultiEff[AD+week*MaxDetectors]*GetDetectorProtons(AD,0))/GetDetectorProtons(0, 0);
            
            if(correction_factor==0)//Avoid NanS
            {
                StatisticalError[AD+week*MaxDetectors] = 0;
                IBDEvents[AD+week*MaxDetectors] = 0;
                ObservedEvents[AD+week*MaxDetectors] = 0;
                
                AccidentalEvents[AD+week*MaxDetectors] = 0;
                AccidentalError[AD+week*MaxDetectors] = 0;
                LiHeEvents[AD+week*MaxDetectors] = 0;
                LiHeError[AD+week*MaxDetectors] = 0;
                FastNeutronEvents[AD+week*MaxDetectors] = 0;
                FastNeutronError[AD+week*MaxDetectors] = 0;
                AmCEvents[AD+week*MaxDetectors] = 0;
                AmCError[AD+week*MaxDetectors] = 0;
            }
            else
            {
                //Correct all events and rates:
                
                StatisticalError[AD+week*MaxDetectors] = StatisticalError[AD+week*MaxDetectors]/correction_factor;
                IBDEvents[AD+week*MaxDetectors] = IBDEvents[AD+week*MaxDetectors]/correction_factor;
                ObservedEvents[AD+week*MaxDetectors] = ObservedEvents[AD+week*MaxDetectors]/correction_factor;
                
                AccidentalEvents[AD+week*MaxDetectors] = AccidentalEvents[AD+week*MaxDetectors]/correction_factor;
                AccidentalError[AD+week*MaxDetectors] = AccidentalError[AD+week*MaxDetectors]/correction_factor;
                LiHeEvents[AD+week*MaxDetectors] = LiHeEvents[AD+week*MaxDetectors]/correction_factor;
                LiHeError[AD+week*MaxDetectors] = LiHeError[AD+week*MaxDetectors]/correction_factor;
                FastNeutronEvents[AD+week*MaxDetectors] = FastNeutronEvents[AD+week*MaxDetectors]/correction_factor;
                FastNeutronError[AD+week*MaxDetectors] = FastNeutronError[AD+week*MaxDetectors]/correction_factor;
                AmCEvents[AD+week*MaxDetectors] = AmCEvents[AD+week*MaxDetectors]/correction_factor;
                AmCError[AD+week*MaxDetectors] = AmCError[AD+week*MaxDetectors]/correction_factor;
            }
        }
    }
    
    if(Nweeks==1)//If fit is not weekly generate inclusive data from weekly data.
    {
        GetnGdInclusiveData(IBDEvents);
        GetnGdInclusiveData(StatisticalError);
        GetnGdInclusiveData(ObservedEvents);
        GetnGdInclusiveData(AccidentalEvents);
        GetnGdInclusiveData(AccidentalError);
        GetnGdInclusiveData(LiHeEvents);
        GetnGdInclusiveData(LiHeError);
        GetnGdInclusiveData(FastNeutronEvents);
        GetnGdInclusiveData(FastNeutronError);
        GetnGdInclusiveData(AmCEvents);
        GetnGdInclusiveData(AmCError);

        for(Int_t AD=0;AD<NADs;AD++)
        {
            //To make sure the inclusive data has changed
            
            std::cout << "Inclusive observed events AD: " << AD << " is : " << ObservedEvents[AD+week*MaxDetectors] << "+/-" << StatisticalError[AD+week*MaxDetectors] << std::endl;
            
            std::cout << "Inclusive ibd events AD: " << AD << " is : " << IBDEvents[AD+week*MaxDetectors] << std::endl;
            
            std::cout << "Inclusive Accidental events AD: " << AD << " is : " << AccidentalEvents[AD+week*MaxDetectors]  << "+/-" << AccidentalError[AD+week*MaxDetectors] << std::endl;
            
            std::cout << "Inclusive LiHe events AD: " << AD << " is : " << LiHeEvents[AD+week*MaxDetectors] << "+/-" << LiHeError[AD+week*MaxDetectors] << std::endl;
            std::cout << "Inclusive FN events AD: " << AD << " is : " << FastNeutronEvents[AD+week*MaxDetectors] << "+/-" << FastNeutronError[AD+week*MaxDetectors] << std::endl;
            
            std::cout << "Inclusive AmC events AD: " << AD << " is : " << AmCEvents[AD+week*MaxDetectors]
            << "+/-" << AmCError[AD+week*MaxDetectors] << std::endl;
        }
    }
    
    //The following are needed in the superflux calculation:
    CalculateInclusiveEfficiencies(InclusiveMultiEff,MultiEff);
    CalculateInclusiveEfficiencies(InclusiveMuonEff,MuonEff);
    CalculateInclusiveFullTime(InclusiveFullTime,FullTime);
    
}


void NominalData :: LoadHydrogenMainData()
{
    Char_t* mainmatrixname;
    
    for(Int_t Volumes = 0; Volumes<VolumeX*VolumeY; Volumes++)
    {
#ifdef UseVolumes
        if(Volumes==0)
        {
            mainmatrixname = Form("./Inputs/HInputs/nH_GdLS_table.txt");
        }
        else
        {
            mainmatrixname = Form("./Inputs/HInputs/nH_LS_table.txt");
        }
#else
        std::cout << " NEED INPUT FILES GIVEN PER CELL " << std::endl;
        exit(EXIT_FAILURE);
#endif
      
        Int_t VolumeIndex;
        std::string line;
        Int_t linenum=0;//<---caution: only increments for lines that do not begin with #
        Int_t week=0;
        ifstream mainfile(mainmatrixname);
        Int_t AD6Index = 0;
        while(!mainfile.eof())
        {
            getline(mainfile,line);
            std::string firstchar = line.substr(0,1);
            
            if(firstchar=="#") continue;//<-- ignore lines with comments
            
            //Special numbers
            if(linenum == 0)
            {
                DataPeriods=atoi(line.c_str());
            }
            
            if(linenum == 1)
            {
                if(atoi(firstchar.c_str())==0) VolumeIndex=0;
                if(atoi(firstchar.c_str())==1) VolumeIndex=1;
                std::cout << "Volume (Gd = 0, Ls = 1): " << VolumeIndex << std::endl;
            }
            
            if(linenum > 2)
            {
                std::cout << "reading " << line << std::endl;
                std::istringstream iss(line);
                Int_t row=0;
                Int_t column=0;
                Double_t readvals[8]={0};
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
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        ObservedEvents[AD6Index+(week-1)*MaxDetectors+VolumeIndex*MaxDetectors*DataPeriods]=readvals[AD];
                        StatisticalError[AD6Index+(week-1)*MaxDetectors+VolumeIndex*MaxDetectors*DataPeriods]=sqrt(ObservedEvents[AD6Index+(week-1)*MaxDetectors+VolumeIndex*MaxDetectors*DataPeriods]);
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
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        FullTime[AD6Index+(week-1)*MaxDetectors]=readvals[AD];
                    }
                }
                //-->muon efficiencies
                if(row==3)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        MuonEff[AD6Index+(week-1)*MaxDetectors]=readvals[AD];
                    }
                }
                //-->dmc efficiencies
                if(row==4)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        MultiEff[AD6Index+(week-1)*MaxDetectors]=readvals[AD];
                    }
                }
                //-->target masses
                if(row==5)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        DetectorMass[AD6Index+MaxDetectors*VolumeIndex]=readvals[AD]/1000000.;
                    }
                }
                //row 6 -->Hydrogen Capture Fraction
                if(row==6)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        nHCaptureFraction[VolumeIndex]=readvals[AD]/100;
                    }
                }
                //row 7 -->Detector Global Efficiency (Efficiency of Prompt energy, Delayed energy, time, distance
                if(row==7)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        m_detectorGlobalEfficiency[VolumeIndex]=readvals[AD];
                    }
                }
                //row 8 -->Detector Global Efficiency error
                if(row==8)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        DetectorEfficiencyRelativeError=readvals[AD]/100;//given in %
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
                //-->acc bg
                if(row==10)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        AccidentalEvents[AD6Index+(week-1)*MaxDetectors+VolumeIndex*MaxDetectors*DataPeriods]=readvals[AD]*FullTime[AD6Index+(week-1)*MaxDetectors];
                    }
                }
                if(row==11)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        AccidentalError[AD6Index+(week-1)*MaxDetectors+VolumeIndex*MaxDetectors*DataPeriods]=readvals[AD]*FullTime[AD6Index+(week-1)*MaxDetectors];
                    }
                }
                //-->fast-n bg
                if(row==12)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        FastNeutronEvents[AD6Index+(week-1)*MaxDetectors+VolumeIndex*MaxDetectors*DataPeriods]=readvals[AD]*FullTime[AD6Index+(week-1)*MaxDetectors];
                    }
                }
                if(row==13)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        FastNeutronError[AD6Index+(week-1)*MaxDetectors+VolumeIndex*MaxDetectors*DataPeriods]=readvals[AD]*FullTime[AD6Index+(week-1)*MaxDetectors];
                    }
                }
                //-->li9 bg
                if(row==14)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        LiHeEvents[AD6Index+(week-1)*MaxDetectors+VolumeIndex*MaxDetectors*DataPeriods]=readvals[AD]*FullTime[AD6Index+(week-1)*MaxDetectors];
                    }
                }
                if(row==15)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        LiHeError[AD6Index+(week-1)*MaxDetectors+VolumeIndex*MaxDetectors*DataPeriods]=readvals[AD]*FullTime[AD6Index+(week-1)*MaxDetectors];
                    }
                }
                //-->amc bg
                if(row==16)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        AmCEvents[AD6Index+(week-1)*MaxDetectors+VolumeIndex*MaxDetectors*DataPeriods]=readvals[AD]*FullTime[AD6Index+(week-1)*MaxDetectors];
                    }
                }
                if(row==17)
                {
                    for(Int_t AD=0;AD<MaxDetectors;AD++)
                    {
                        AD6Index = AD;
                        
                        if(NADs != MaxDetectors)//6AD analysis
                        {
                            if(AD==3||AD==(MaxDetectors-1))
                            {
                                continue;//6AD analysis doesn't use AD4, AD8 info
                            }
                            if(AD>3)
                            {
                                AD6Index = AD-1;//Correct index
                            }
                        }
                        
                        AmCError[AD6Index+(week-1)*MaxDetectors+VolumeIndex*MaxDetectors*DataPeriods]=readvals[AD]*FullTime[AD6Index+(week-1)*MaxDetectors];
                    }
                }
            }
            linenum++;//only lines >2
        }
    }
    
    for(Int_t AD=0;AD<NADs;AD++)
    {
        for(Int_t week = 0; week<DataPeriods; week++)
        {
            for(Int_t idx=0; idx<VolumeX; idx++)
            {
                for(Int_t idy=0; idy<VolumeY; idy++)
                {
                    Int_t IndexVariant = AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods;
                    
                    IBDEvents[IndexVariant] = (ObservedEvents[IndexVariant]-(AccidentalEvents[IndexVariant]+LiHeEvents[IndexVariant]+FastNeutronEvents[IndexVariant]+AmCEvents[IndexVariant]));//Apply background  suppression for events with AD effects, then correct AD effects and scale to days
                }
            }
        }
    }
    
    CorrectnHEvents();
}

void NominalData :: CorrectnHEvents()
{
    for(Int_t AD=0;AD<NADs;AD++)
    {
        for(Int_t week = 0; week<DataPeriods; week++)
        {
            for(Int_t idx=0; idx<VolumeX; idx++)
            {
                for(Int_t idy=0; idy<VolumeY; idy++)
                {
                    Int_t IndexVariant = AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods;
                    Int_t IndexConstant = AD+week*MaxDetectors;
                    Double_t correction_factor = (FullTime[IndexConstant]*MuonEff[IndexConstant]*MultiEff[IndexConstant]*GetDetectorProtons(AD,idx))/GetDetectorProtons(0,idx);
                    
                    if(correction_factor==0)//Avoid NaNs
                    {
                        IBDEvents[IndexVariant] = 0;
                        
                        ObservedEvents[IndexVariant] = 0;
                        
                        StatisticalError[IndexVariant] = 0;
                        
                        AccidentalEvents[IndexVariant] = 0;
                        AccidentalError[IndexVariant] = 0;
                        LiHeEvents[IndexVariant] = 0;
                        LiHeError[IndexVariant] = 0;
                        FastNeutronEvents[IndexVariant] =0;
                        FastNeutronError[IndexVariant] = 0;
                        AmCEvents[IndexVariant] = 0;
                        AmCError[IndexVariant] = 0;
                    }
                    else
                    {
                    IBDEvents[IndexVariant] = IBDEvents[IndexVariant]/correction_factor;
                    
                    ObservedEvents[IndexVariant] = ObservedEvents[IndexVariant]/correction_factor;
                    
                    StatisticalError[IndexVariant] = StatisticalError[IndexVariant]/correction_factor;
                    
                    AccidentalEvents[IndexVariant] = AccidentalEvents[IndexVariant]/correction_factor;
                    AccidentalError[IndexVariant] = AccidentalError[IndexVariant]/correction_factor;
                    LiHeEvents[IndexVariant] = LiHeEvents[IndexVariant]/correction_factor;
                    LiHeError[IndexVariant] = LiHeError[IndexVariant]/correction_factor;
                    FastNeutronEvents[IndexVariant] = FastNeutronEvents[IndexVariant]/correction_factor;
                    FastNeutronError[IndexVariant] = FastNeutronError[IndexVariant]/correction_factor;
                    AmCEvents[IndexVariant] = AmCEvents[IndexVariant]/correction_factor;
                    AmCError[IndexVariant] = AmCError[IndexVariant]/correction_factor;
                    }
                    
                    //                        std::cout << "? " << ObservedEvents[IndexVariant] << " = " << ObservedEvents[IndexVariant]/correction_factor << std::endl;
                    
                }
            }
        }
    }

    if(Nweeks==1)//If the fit is not weekly, generate inclusive data from weekly data.
    {
        GetnHInclusiveData(IBDEvents);
        GetnHInclusiveData(ObservedEvents);
        GetnHInclusiveData(StatisticalError);
        GetnHInclusiveData(AccidentalEvents);
        GetnHInclusiveData(AccidentalError);
        GetnHInclusiveData(LiHeEvents);
        GetnHInclusiveData(LiHeError);
        GetnHInclusiveData(FastNeutronEvents);
        GetnHInclusiveData(FastNeutronError);
        GetnHInclusiveData(AmCEvents);
        GetnHInclusiveData(AmCError);
        
        for(Int_t AD=0;AD<NADs;AD++)
        {
            for(Int_t idx=0; idx<VolumeX; idx++)
            {
                for(Int_t idy=0; idy<VolumeY; idy++)
                {
                    Int_t IndexInclusive = AD+idx*MaxDetectors*DataPeriods;
                    //To make sure the inclusive data has changed
                    
                    std::cout << "Inclusive observed events AD: " << AD << " is : " << ObservedEvents[IndexInclusive] << " in " << idx << std::endl;
                    
                    std::cout << "Inclusive IBD events AD: " << AD << " is : " << IBDEvents[IndexInclusive] << " in " << idx << std::endl;

                    std::cout << "Inclusive BKGD events AD: " << AD << " is : " << AccidentalEvents[IndexInclusive] << AccidentalEvents[IndexInclusive] << " + " << LiHeEvents[IndexInclusive] << " + " << FastNeutronEvents[IndexInclusive] << " + " << AmCEvents[IndexInclusive] << " = " << (AccidentalEvents[IndexInclusive]+LiHeEvents[IndexInclusive]+FastNeutronEvents[IndexInclusive]+AmCEvents[IndexInclusive]) << " in " << idx << std::endl;

                    std::cout <<  ObservedEvents[IndexInclusive] <<  " = " << IBDEvents[IndexInclusive]  + (AccidentalEvents[IndexInclusive]+LiHeEvents[IndexInclusive]+FastNeutronEvents[IndexInclusive]+AmCEvents[IndexInclusive]) << std::endl;
                    
                        //To make sure the inclusive data has changed
                        
//                        std::cout << "Inclusive observed events AD: " << AD << " is : " << ObservedEvents[IndexInclusive] << "+/-" << StatisticalError[IndexInclusive] << std::endl;
                    
//                        std::cout << "Inclusive ibd events AD: " << AD << " is : \t " << IBDEvents[IndexInclusive] << std::endl;
                    
                        std::cout << "Inclusive Accidental events AD: " << AD << " is : \t" << AccidentalEvents[IndexInclusive]  << "+/-" << AccidentalError[IndexInclusive] << std::endl;
                        
                        std::cout << "Inclusive LiHe events AD: " << AD << " is : \t" << LiHeEvents[IndexInclusive] << "+/-" << LiHeError[IndexInclusive] << std::endl;
                        std::cout << "Inclusive FN events AD: " << AD << " is : \t" << FastNeutronEvents[IndexInclusive] << "+/-" << FastNeutronError[IndexInclusive] << std::endl;
                        
                        std::cout << "Inclusive AmC events AD: " << AD << " is : " << AmCEvents[IndexInclusive]
                        << "+/-" << AmCError[IndexInclusive] << std::endl;
                        
                        
                    
                }
            }
        }
        
        
//        for(Int_t AD=0;AD<NADs;AD++)
//        {
//            for(Int_t week=0;week<DataPeriods;week++)
//            {
//                InclusiveMultiEff[AD+week*MaxDetectors] = MultiEff[AD+week*MaxDetectors];
//                InclusiveMuonEff[AD+week*MaxDetectors] = MuonEff[AD+week*MaxDetectors];
//                InclusiveFullTime[AD+week*MaxDetectors] = FullTime[AD+week*MaxDetectors];
//                
//            }
//        }
    }
    //The following are needed in the superflux calculation:
    CalculateInclusiveEfficiencies(InclusiveMultiEff,MultiEff);
    CalculateInclusiveEfficiencies(InclusiveMuonEff,MuonEff);
    CalculateInclusiveFullTime(InclusiveFullTime,FullTime);
}

void NominalData :: ReadChristineReactorSpectrum()
{
    
    std::ifstream fileData;
    
#ifdef BlindedAnalysis
    fileData.open("./ReactorInputs/p12c_blinded/combined/nNu_Nom_combined.txt",std::fstream::in);
    std::cout << " READING BLINDED CHRISTINE SPECTRUM " << std::endl;
#else
    fileData.open("./ReactorInputs/p12c_unblinded/combined/nNu_Nom_combined_huber-french.txt",std::fstream::in);
    std::cout << " READING UNBLINDED CHRISTINE SPECTRUM " << std::endl;
#endif
    
    if(!fileData.is_open())
    {
        std::cout << "Christine Reactor file not open" << std::endl;
        
        exit(EXIT_FAILURE);
    }
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
            m_dNdE_nom[reactor+NReactors*curSample] = dNdE[reactor]* 1.0e18;
        }
        
        curSample++;
    }
    
    curSample--;//219, from 1.8 to 12.75 in 0.05 steps.
    
    fileData.close();
    
    m_eMin = eMin;
    m_eMax = eMax;
    m_nSamples = curSample;
    
    std::cout << " P12C REACTOR FROM CHRISTINE'S TXT FILE: " << std::endl;
    std::cout << " MINIMUM ENERGY: " << m_eMin << std::endl;
    std::cout << " MAX ENERGY: " << m_eMax << std::endl;
    std::cout << " SAMPLES: " << m_nSamples << std::endl;
}

void NominalData :: ReadIHEPReactorSpectrum()
{
    std::cout << " READING IHEP SPECTRUM " << std::endl;
    
    const Int_t NIHEPReactorBins = 32+1; // 1.5, 1.75, ..., 9.25, 9.5 in 0.25 MeV steps
    Double_t eMin=1.5;
    Double_t eMax=9.5;
    binWidth = 0.25;
    
    Int_t curPeriod = 0;
    std::fstream IHEPfileData;
    Double_t trash;
    
    for (Int_t reactor = 0; reactor < NReactors; reactor++)
    {
        curPeriod = 0;
        
        Double_t activePeriods = 0;//Account for the number of weeks it's been on.
        Double_t m_dNdE_temporal = 0;
        
        switch(reactor)
        {
            case 0://Daya Bay A
                IHEPfileData.open("./ReactorInputs/P14A/DayaBayA_2011-12-24_2013-11-27.txt",std::fstream::in);
                break;
            case 1://Daya Bay B
                IHEPfileData.open("./ReactorInputs/P14A/DayaBayB_2011-12-24_2013-11-27.txt",std::fstream::in);
                break;
            case 2://LingAo IA
                IHEPfileData.open("./ReactorInputs/P14A/LingAoIA_2011-12-24_2013-11-27.txt",std::fstream::in);
                break;
            case 3://LingAo IB
                IHEPfileData.open("./ReactorInputs/P14A/LingAoIB_2011-12-24_2013-11-27.txt",std::fstream::in);
                break;
            case 4://LingAo IIA
                IHEPfileData.open("./ReactorInputs/P14A/LingAoIIA_2011-12-24_2013-11-27.txt",std::fstream::in);
                break;
            case 5://LingAo IIB
                IHEPfileData.open("./ReactorInputs/P14A/LingAoIIB_2011-12-24_2013-11-27.txt",std::fstream::in);
                break;
            default:
                std::cout << "Reactor data file not found" << std::endl;
                
                exit(EXIT_FAILURE);
        }
        
        if(!IHEPfileData.is_open())
        {
            std::cout << "P14A Reactor file not open" << std::endl;
            
            exit(EXIT_FAILURE);
        }
        
        while(!IHEPfileData.eof())
        {
            for(Int_t idr = 0; idr<5; idr++)//first 5 columns are the reactor ID
            {
                IHEPfileData >> trash;
                
                //std::cout << trash << std::endl;
                
            }
            for (Int_t bin = 0; bin <NIHEPReactorBins; bin++)
            {
                IHEPfileData >> m_dNdE_temporal;
                
                m_dNdE_nom[reactor+NReactors*bin] = m_dNdE_temporal+m_dNdE_nom[reactor+NReactors*bin];
                // std::cout << "reactor : " << reactor << " bin: " << bin << ", reactor spectrum: "<< m_dNdE_nom[reactor+NReactors*bin] << std::endl;
            }
            
            //Last value seems to be way off.
            IHEPfileData >> trash;
            
            
            if(m_dNdE_temporal!=0)
            {
                activePeriods++;
            }
            
            curPeriod++;
        }
        
        IHEPfileData.close();
        
        curPeriod--;
        
        //Average p14A nominal reactor:
        for (Int_t bin = 0; bin <=NIHEPReactorBins; bin++)
        {
            m_dNdE_nom[reactor+NReactors*bin] = m_dNdE_nom[reactor+NReactors*bin]/activePeriods;
        }
        
        std::cout << "Reactor " << reactor << " has been active for: " << activePeriods << " periods" << std::endl;
        
    }
    
    m_eMin = eMin;
    m_eMax = eMax;
    m_nSamples = NIHEPReactorBins;
    
    NReactorPeriods = curPeriod;
    
    std::cout << " P14A REACTOR FROM IHEP TXT FILE: " << std::endl;
    std::cout << " PERIODS: " << curPeriod << std::endl;
    std::cout << " MINIMUM ENERGY: " << m_eMin << std::endl;
    std::cout << " MAX ENERGY: " << m_eMax << std::endl;
    std::cout << " SAMPLES: " << m_nSamples << std::endl;
    
    IHEPReactorModel = 1;
}

bool NominalData :: UsingIHEPReactorModel()
{
    return IHEPReactorModel;
}
void NominalData :: ReadChristineCovMatrix()
{
    std::cout << " READING CHRISTINE COVARIANCE MATRIX " << std::endl;
    
    // Read covarianvematrix
#ifdef BlindedAnalysis
    std::ifstream fileData_mcov("./ReactorInputs/p12c_blinded/combined/nNu_Mcov_combined.txt");
#else
    std::ifstream fileData_mcov("./ReactorInputs/p12c_unblinded/combined/nNu_Mcov_combined_huber-french_u238cor.txt");
#endif
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
    //    #ifdef PrintEps
    //        TCanvas* c2 = new TCanvas("","");
    //        covmatrix.Draw("colz");
    //        c2->Print("./Images/Reactor/ReactorCovMatrix.eps");
    //        delete c2;
    //    #endif
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
        OAVMatrix= 1;
        NLMatrix= 1;
        ResolutionMatrix= 1;
        Sin22t12Matrix= 1;
        EfficiencyMatrix= 1;
        AllDetectorSystematicsMatrix= 1;
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
        OAVMatrix = 0;
        NLMatrix= 0;
        ResolutionMatrix= 0;
        Sin22t12Matrix= 0;
        EfficiencyMatrix= 0;
        AllDetectorSystematicsMatrix= 0;
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

Double_t NominalData :: GetObservedEvents(Int_t i,Int_t j,Int_t idx, Int_t idy)
{
    return ObservedEvents[i+j*MaxDetectors+idx*MaxPeriods*MaxDetectors+idy*MaxPeriods*MaxDetectors*VolumeX];
}

Int_t NominalData :: GetNReactorPeriods()
{
    if(!isH)
    {
        NReactorPeriods = 20;//Hard coded
    }
    return NReactorPeriods;
}

Int_t NominalData :: GetDataPeriods()
{
    return DataPeriods;
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

void NominalData :: SetDataPeriods(Int_t dataperiods)
{
    DataPeriods = dataperiods;
}

void NominalData :: CalculateBinning()
{
    std::cout << " Recalculating binning " << std::endl;
    
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
            n_etrue_bins = 39;//match reactor flux data for P12E
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
        
        std::cout << " True bin width " << TrueBinWidth << std::endl;
        std::cout << " Vis bin width " << VisBinWidth << std::endl;
        
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
        else
        {
            n_evis_bins=34;
            n_etrue_bins=39;
            
            for (Int_t i = 0; i <= n_etrue_bins; i++)
            {
                enu_bins[i] = 0.2 * i + Emin;
            }
            
            evis_bins[0] = 1.5;
            
            for (Int_t i = 0; i < n_evis_bins-1; i++)
            {
                evis_bins[i+1] = 0.2 * i + 1.6;
            }
            evis_bins[n_evis_bins] = 12;
        }
    }
    
    std::cout << " n_evis_bins " << n_evis_bins << std::endl;
    std::cout << "  n_etrue_bins " << n_etrue_bins << std::endl;
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

void NominalData :: ReadToyEventsByCell()//Included in nHToyMC.h right now
{
    std::string line;
    
    string SystematicS;
    
    //Detector systematics have different matrices:
    if(IAVMatrix)
    {
        SystematicS = "IAV";
    }
    else if(OAVMatrix)
    {
        SystematicS = "OAV";
    }
    else if(NLMatrix)
    {
        SystematicS = "NL";
    }
    else if(RelativeEnergyScaleMatrix)
    {
        SystematicS = "RelativeEnergyScale";
    }
    else if(ResolutionMatrix)
    {
        SystematicS = "Resolution";
    }
//    else if(EfficiencyMatrix)
//    {
//        SystematicS = "Efficiency";
//    }
    else if(AllDetectorSystematicsMatrix)
    {
        SystematicS = "AllDetectorSystematics";
    }
    else
    {
        SystematicS = "Nominal";
    }
    
#ifndef EREC_COMPARISON
    
    //Save txt file in case we want to use it externally:
    ifstream mainfile(("./Inputs/HInputs/"+SystematicS+"ToyMCEventRatio.txt").c_str());
#else
    
    ifstream mainfile(("./Inputs/HInputs/E_REC_"+SystematicS+"ToyMCEventRatio.txt").c_str());
    
#endif
    
    
    Int_t linenum=0;//<---caution: only increments for lines that do not begin with #
    
    while(!mainfile.eof())
    {
        std::getline(mainfile,line);
        std::string firstchar = line.substr(0,1);
        
        if(firstchar=="#") continue;//<-- ignore lines with comments
        
        std::istringstream iss(line);
        
        if(linenum == 0)
        {
            Int_t AD = 0;
            
            iss >> AD >> ADIntegral[AD];
            
            std::cout << "line " << linenum << " the integral in AD " << AD << " is " << ADIntegral[AD] << std::endl;
            
            linenum++;
        }
        else if(linenum <=VolumeX*VolumeY)
        {
            Int_t AD = 0;
            Int_t idx = 0;
            Int_t idy = 0;
            
            iss >> AD >> idx >> idy >> CellIntegral[AD+NADs*idx+NADs*VolumeX*idy];
            //                if(column==0) AD=atoi(firstchar.c_str());
            //
            //                if(column==1) idx=atoi(sub.c_str());
            //
            //                if(column==2) idy=atoi(sub.c_str());
            //
            //                if(column==3)
            
            std::cout << "line " << linenum << " the integral in AD " << AD << " cell " << idx << " , " << idy << " is " << CellIntegral[AD+NADs*idx+NADs*VolumeX*idy] << std::endl;
            
            linenum++;
            
            if(linenum>VolumeX*VolumeY)
            {
                linenum = 0;//reset counter
            }
        }
    }

    //Draw it:
#ifndef UseVolume
    TH2D* MapEvents_ratio2center = new TH2D("MapEventsRatio2center","MapEventsRatio2center",VolumeX,VolumeX_lower,VolumeX_upper,VolumeY,VolumeY_lower,VolumeY_upper);
    
    for(Int_t AD = 0; AD<NADs; AD++)
    {
        for(Int_t idx = 0; idx<VolumeX; idx++)
        {
            for(Int_t idy = 0; idy<VolumeY; idy++)
            {
                MapEvents_ratio2center->SetBinContent(idx+1,idy+1,CellIntegral[AD+NADs*idx+NADs*VolumeX*idy]/CellIntegral[AD+NADs*0+NADs*VolumeX*Int_t(VolumeY/2)]);
            }
        }
    }
    
    TCanvas* MapEventC = new TCanvas("MapEventRatio2Center","MapEventRatio2Center");
    
    MapEventC->cd(1);
    
    MapEvents_ratio2center->Draw("colz");
    
    MapEventC->Print("./Images/Hydrogen/Detector/NominalDataMapEventRatio2Center.eps");
    
    delete MapEvents_ratio2center;
    delete MapEventC;
#endif
    
    TH2D* MapEvents_ratio2total = new TH2D("MapEventsRatio2total","MapEventsRatio2total",VolumeX,VolumeX_lower,VolumeX_upper,VolumeY,VolumeY_lower,VolumeY_upper);
    
    for(Int_t AD = 0; AD<NADs; AD++)
    {
        for(Int_t idx = 0; idx<VolumeX; idx++)
        {
            for(Int_t idy = 0; idy<VolumeY; idy++)
            {
                PercentualEvents[AD+NADs*idx+NADs*VolumeX*idy] = (CellIntegral[AD+NADs*idx+NADs*VolumeX*idy]/ADIntegral[AD]);

                MapEvents_ratio2total->SetBinContent(idx+1,idy+1,PercentualEvents[AD+NADs*idx+NADs*VolumeX*idy]);
            }
        }
    }
    
    TCanvas* MapEventTotalC = new TCanvas("MapEventRatio2Total","MapEventRatio2Total");
    
    MapEventTotalC->cd(1);
    
    MapEvents_ratio2total->Draw("colz");
    
    MapEventTotalC->Print("./Images/Hydrogen/Detector/NominalMapEventRatio2Total.eps");
    
    delete MapEventTotalC;
    
    TFile* SaveFile = new TFile("./Inputs/HInputs/MapEventRatio2Total.root","recreate");
    {
        MapEvents_ratio2total->Write();
    }
    
    delete SaveFile;
    
    delete MapEvents_ratio2total;
}

Double_t NominalData :: GetEventsByCell(Int_t AD,Int_t idx, Int_t idy)
{
    return PercentualEvents[AD+NADs*idx+NADs*VolumeX*idy];
}

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

//void NominalData :: GetnHInclusiveData(Double_t* InclusiveData)
//{
//    //Weight weekly rates to get a inclusive rate to use in the fitter
//    Double_t WeeklyData[DataPeriods*MaxDetectors*VolumeX*VolumeY];
//    
//    for(Int_t AD=0;AD<NADs;AD++)
//    {
//        for(Int_t week=0;week<DataPeriods;week++)
//        {
//            Int_t IndexConstant = AD+week*MaxDetectors;
//
//            for(Int_t idx=0; idx<VolumeX; idx++)
//            {
//                Int_t IndexWeekly = AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods;
//
//                for(Int_t idy=0; idy<VolumeY; idy++)
//                {
//                    WeeklyData[IndexWeekly] = InclusiveData[IndexWeekly];
//         
//                    
////                        std::cout << " Weekly observed events AD: " << AD << "Volume " << idx << " is : " << WeeklyData[IndexWeekly] << " for week : " << week << " scaled " << (FullTime[IndexConstant]*MuonEff[IndexConstant]*MultiEff[IndexConstant]*GetDetectorProtons(AD,idx))/GetDetectorProtons(0, idx) << std::endl;
//                
//                    
//                    InclusiveData[IndexWeekly] = 0;
//                }
//            }
//        }
//    }
//
//    for(Int_t AD=0;AD<NADs;AD++)
//    {
//        for(Int_t idx=0; idx<VolumeX; idx++)
//        {
//            for(Int_t idy=0; idy<VolumeY; idy++)
//            {
//                Double_t scalesum=0;
//
//                for(Int_t week=0;week<DataPeriods;week++)
//                {
//                    Double_t scale = FullTime[AD+week*MaxDetectors]*MuonEff[AD+week*MaxDetectors]*MultiEff[AD+week*MaxDetectors]*GetDetectorProtons(AD,idx)/GetDetectorProtons(0,idx);
//                 
////                    if(week==0)
////                    {
////                        std:: cout << FullTime[AD+week*MaxDetectors] << " " << MuonEff[AD+week*MaxDetectors] << " " << MultiEff[AD+week*MaxDetectors]  << " Volume " << idx << " This should be different for volumes 0,1 : " << GetDetectorProtons(AD,idx) << " " << GetDetectorProtons(0,idx) << std::endl;
////                        std::cout<< " WEEK " << week << " SCALE: " << scale << " IDX, AD " << idx << AD << std::endl;
////                    }
//                    
//                    WeeklyData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods]*=scale;//At this point The background and observed data rates/errors already have been corected for the muon and multiplicity efficiencies, also for livetime and protons.
//                    
//                    if(DataPeriods==Nweeks)//The data and the fit periods are the same, this can be use to fit inclusive data or fit the data weekly (studies of theta13/deltam fit depending on the weekly variations)
//                    {
//                        InclusiveData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods] = WeeklyData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
//                    }
//                    else
//                    {
//                        scalesum+=scale;
//                        
//                        if(ErrorMode)
//                        {
//                            InclusiveData[AD+idx*MaxDetectors*DataPeriods]+=pow(WeeklyData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods],2);
//                        }
//                        else
//                        {
//                            InclusiveData[AD+idx*MaxDetectors*DataPeriods]+=WeeklyData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
//                        }
//                    }
//                    
//                    if(scale!=0)
//                    {
//                    WeeklyData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods]*=1./scale;//Return the vector to its original state
//                    }
//                }
//                
//                if(scalesum==0)
//                {
//                    InclusiveData[AD+idx*MaxDetectors*DataPeriods]*=0;
//                    std::cout << " WEEKLY DATA IS EMPTY, FILE NOT READ? " << std::endl;
//                    
//                    exit(EXIT_FAILURE);
//                }
//                else
//                {
//                    if(ErrorMode)
//                    {
//                        InclusiveData[AD+idx*MaxDetectors*DataPeriods]=sqrt(InclusiveData[AD+idx*MaxDetectors*DataPeriods])*1./scalesum;
//                    }
//                    else
//                    {
//                        InclusiveData[AD+idx*MaxDetectors*DataPeriods]*=1./scalesum;
//                    }
//                }
//            }
//        }
//    }
//}

void NominalData :: GetnHInclusiveData(Double_t* InclusiveData)
{
    //Weight weekly rates to get a inclusive rate to use in the fitter
    Double_t WeeklyData[DataPeriods*MaxDetectors*VolumeX*VolumeY];
    
    for(Int_t AD=0;AD<NADs;AD++)
    {
        for(Int_t week=0;week<DataPeriods;week++)
        {
            for(Int_t idx=0; idx<VolumeX; idx++)
            {
                Int_t IndexWeekly = AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods;
                
                for(Int_t idy=0; idy<VolumeY; idy++)
                {
                    WeeklyData[IndexWeekly] = InclusiveData[IndexWeekly];
                    
                    
                    //                        std::cout << " Weekly observed events AD: " << AD << "Volume " << idx << " is : " << WeeklyData[IndexWeekly] << " for week : " << week << " scaled " << (FullTime[IndexConstant]*MuonEff[IndexConstant]*MultiEff[IndexConstant]*GetDetectorProtons(AD,idx))/GetDetectorProtons(0, idx) << std::endl;
                    
                    
                    InclusiveData[IndexWeekly] = 0;
                }
            }
        }
    }
    
    for(Int_t AD=0;AD<NADs;AD++)
    {
        for(Int_t idx=0; idx<VolumeX; idx++)
        {
            for(Int_t idy=0; idy<VolumeY; idy++)
            {
                for(Int_t week=0;week<DataPeriods;week++)
                {
                    InclusiveData[AD+idx*MaxDetectors*DataPeriods]+=WeeklyData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
                }
            }
        }
    }
}

//void NominalData :: GetnGdInclusiveData(Double_t* InclusiveData, bool ErrorMode)
//{
//    
//    //Weight weekly rates to get a inclusive rate to use in the fitter
//    Double_t WeeklyData[DataPeriods*MaxDetectors*VolumeX*VolumeY];
//    
//    for(Int_t AD=0;AD<NADs;AD++)
//    {
//        for(Int_t week=0;week<DataPeriods;week++)
//        {
//            for(Int_t idx=0; idx<VolumeX; idx++)
//            {
//                for(Int_t idy=0; idy<VolumeY; idy++)
//                {
//                    WeeklyData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods] = InclusiveData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
//                    
//                    InclusiveData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods] = 0;
//                }
//            }
//        }
//    }
//
//    //Weight weekly rates to get a inclusive rate to use in the fitter
//    
//    for(Int_t AD=0;AD<NADs;AD++)
//    {
//        Double_t scalesum=0;
//        
//        for(Int_t week=0;week<DataPeriods;week++)
//        {
//            Double_t scale = FullTime[AD+week*MaxDetectors]*MuonEff[AD+week*MaxDetectors]*MultiEff[AD+week*MaxDetectors]*GetDetectorProtons(AD,0)/GetDetectorProtons(0,0);
//
//            if(DataPeriods==Nweeks)//The data and the fit periods are the same, this can be use to fit inclusive data or fit the data weekly (studies of theta13/deltam fit depending on the weekly variations)
//            {
//                InclusiveData[AD+week*MaxDetectors] = WeeklyData[AD+week*MaxDetectors]*scale;
//            }
//            else
//            {
//                scalesum+=scale;
//                
//                WeeklyData[AD+week*MaxDetectors]*=FullTime[AD+week*MaxDetectors]*GetDetectorProtons(AD,0)/GetDetectorProtons(0,0);//The background and observed data rates already include the muon and multiplicity efficiencies, here include livetime and protons.
//                
//                if(ErrorMode)
//                {
//                    InclusiveData[AD]+=pow(WeeklyData[AD+week*MaxDetectors],2);
//                }
//                else
//                {
//                    InclusiveData[AD]+=WeeklyData[AD+week*MaxDetectors];
//                }
//            }
//        }
//        
//        if(scalesum==0)
//        {
//            InclusiveData[AD]*=0;
//            std::cout << " WEEKLY DATA IS EMPTY, FILE NOT READ? " << std::endl;
//            
//            exit(EXIT_FAILURE);
//        }
//        else
//        {
//            if(ErrorMode)
//            {
//                InclusiveData[AD]=sqrt(InclusiveData[AD])*1./scalesum;
//            }
//            else
//            {
//                InclusiveData[AD]*=1./scalesum;
//            }
//        }
//    }
//}

void NominalData :: GetnGdInclusiveData(Double_t* InclusiveData)
{
    //Just add events and absolute errors
    
    //Sum weekly rates to get a inclusive rate to use in the fitter
    Double_t WeeklyData[DataPeriods*MaxDetectors*VolumeX*VolumeY];
    
    for(Int_t AD=0;AD<NADs;AD++)
    {
        for(Int_t week=0;week<DataPeriods;week++)
        {
            for(Int_t idx=0; idx<VolumeX; idx++)
            {
                for(Int_t idy=0; idy<VolumeY; idy++)
                {
                    WeeklyData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods] = InclusiveData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods];
                    
                    InclusiveData[AD+week*MaxDetectors+idx*MaxDetectors*DataPeriods] = 0;
                }
            }
        }
    }
    
    //Sum weekly rates to get a inclusive rate to use in the fitter
    
    for(Int_t AD=0;AD<NADs;AD++)
    {
        for(Int_t week=0;week<DataPeriods;week++)
        {
            InclusiveData[AD]+=WeeklyData[AD+week*MaxDetectors];
        }
    }
}

void NominalData:: CalculateInclusiveEfficiencies(Double_t* InclusiveData, Double_t* WeeklyData)
{
    for(Int_t AD=0;AD<NADs;AD++)
    {
        Double_t ScaleSum = 0;

        InclusiveData[AD] = 0;

        for(Int_t week=0;week<DataPeriods;week++)
        {
            ScaleSum += FullTime[AD+week*MaxDetectors];

            if(WeeklyData[AD+week*MaxDetectors]!=0)
            {
                InclusiveData[AD] += WeeklyData[AD+week*MaxDetectors]*FullTime[AD+week*MaxDetectors];
            }
        }
        
        InclusiveData[AD]*=1./ScaleSum;
        
        std::cout << " AVERAGE EFFICIENCY AD IN THE WHOLE PERIOD" << AD << " IS : " << InclusiveData[AD] << std::endl;
    }
}

void NominalData:: CalculateInclusiveFullTime(Double_t* InclusiveData, Double_t* WeeklyData)
{
    
    for(Int_t AD=0;AD<NADs;AD++)
    {
        InclusiveData[AD] = 0;
        
        for(Int_t week=0;week<DataPeriods;week++)
        {
            InclusiveData[AD] += WeeklyData[AD+week*MaxDetectors];//Total fulltimes are additive
        }
        
        std::cout << " TOTAL FULL TIME AD " << AD << " IS : " << InclusiveData[AD] << std::endl;
    }
}

//Reads baseline distances from file
void NominalData :: ReadDistances(Char_t* distanceFileName)
{
    std::ifstream infile(distanceFileName);
    std::string line;
    
    getline(infile,line); //To throw away the first line
    
    //Get baselines from text file into a Matrix
    for(Int_t i=0;i<NADs;i++)
    {
        infile >> ADdistances[i*NReactors+0] >> ADdistances[i*NReactors+1] >> ADdistances[i*NReactors+2] >> ADdistances[i*NReactors+3] >> ADdistances[i*NReactors+4] >> ADdistances[i*NReactors+5];
        
        std::cout << " Distances in Nominal Data : " << std::endl;
        
        printf("Baseline AD%d to D1 is: %f \n", i+1, ADdistances[i*NReactors+0]);
        printf("Baseline AD%d to D2 is: %f \n", i+1, ADdistances[i*NReactors+1]);
        printf("Baseline AD%d to L1 is: %f \n", i+1, ADdistances[i*NReactors+2]);
        printf("Baseline AD%d to L2 is: %f \n", i+1, ADdistances[i*NReactors+3]);
        printf("Baseline AD%d to L3 is: %f \n", i+1, ADdistances[i*NReactors+4]);
        printf("Baseline AD%d to L4 is: %f \n", i+1, ADdistances[i*NReactors+5]);
        //        for(int j=0;j<NReactors;j++)
        //        {
        //            ADdistances[i*NReactors+j]=ADdistances[i*NReactors+j]*km; //Km instead of meters
        //        }
    }
    
    //Test:
#ifdef TestAllTheSame
    
    for(Int_t i=0;i<NADs;i++)
    {
        for(Int_t j=0;j<NReactors;j++)
        {
            ADdistances[i*NReactors+j]=1000;
        }
    }
}
#endif
infile.close();
}

Double_t NominalData :: GetDistances(Int_t AD, Int_t Reactor)
{
    return  ADdistances[Reactor+AD*NReactors];
}