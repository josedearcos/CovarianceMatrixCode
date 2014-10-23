#pragma once
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
#include <vector>
#include "NominalData.h"
#include "Prediction.h"
#include <math.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TCanvas.h"

const bool StatisticalFluctuation=0;//To fluctuate each bin as a statistical fluctuation
const bool WriteOutput=0;//To save the covariance matrices in a .txt file.

//NL
const Int_t n_bcw_positron_nl = 1000;
const Int_t n_unified_nl_points = 500;
const Int_t m_num_unified_nl_pars = 4; // 4 marginal curves in the final model

//Particle masses
const Double_t Me = 0.510999; // MeV
const Double_t Mn = 939.565; // MeV
const Double_t Mp = 938.272; // MeV
const Int_t Combine = 0;

class CovarianceMatrix3
{
private:
    NominalData* Nom;
    TRandom3* rand;
    Prediction* Pred;
    TCanvas* c;

    enum Systematic{IsotopeE,PowerE, RelativeEnergyE, AbsoluteEnergyE, RelativeEnergyOffsetE, AbsoluteEnergyOffsetE, IAVE, NLE, ResolutionE,Sin22t12E};
    Systematic SystematicE;
    typedef Systematic SystematicType;
    enum Background{VaryAccidentalE, VaryLiHeE, VaryFastNeutronsE, VaryAmCE, DistortLiHeE, DistortFastNeutronsE, DistortAmCE};
    Background BackgroundE;
    typedef Background BackgroundType;
    enum NLModel{BCWE, LBNLE, IHEPE, UnifiedE};//Add here new NL models
    NLModel NLModelE;
    
    //AD configuration parameters:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    Int_t MaxNear;
    Int_t MaxFar;
    Int_t MaxBins;
    Int_t hall;
    
    //Binning parameters:
    bool LinearBinning;
    bool SimpleReactorModel;
    bool flag;
    bool flagX;
    Int_t n_evis_bins;
    Int_t n_etrue_bins;

    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    Int_t Nweeks;
    Double_t BinWidth;
    Int_t TotalBins;
    Int_t Sum;

    Int_t NSamples;//Number of samples used in the generation of covariance matrices

    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.
    
    Double_t LiveTime[MaxDetectors][MaxPeriods];
    Double_t DetectorEfficiency[MaxDetectors][MaxPeriods];
    Double_t IBDRate[MaxDetectors][MaxPeriods];
    
    //  Background rates:
    Double_t ScaleAcc[MaxDetectors][MaxPeriods];
    Double_t ScaleLiHe[MaxDetectors][MaxPeriods];
    Double_t ScaleFN[MaxDetectors][MaxPeriods];
    Double_t ScaleAmC[MaxDetectors][MaxPeriods];
    
    //Response Matrix:
    Double_t PositronEnergy;
    Int_t PositronEnergyIndex;
    Double_t Norma[MaxNbins+1];
    TH1F* PositronTrueSpectrumH[MaxDetectors][MatrixBins];
    TH1F* PositronIAVSpectrumH[MaxDetectors][MatrixBins];
    TH1F* PositronNLSpectrumH[MaxDetectors][MatrixBins];
    TH1F* PositronVisibleSpectrumH[MaxDetectors][MatrixBins];
    
    //IAV:
    Double_t IAVNominalError; // relative uncertainty of the IAV thickness
    Double_t IAVError[MaxDetectors]; // relative uncertainty of the IAV thickness
    Double_t NominalIAVMatrix[MatrixBins][MatrixBins];
    Double_t IAVMatrix[MaxDetectors][MatrixBins][MatrixBins];
    Double_t NominalIAVMatrixFrac[MatrixBins];
    Double_t IAVMatrixFrac[MaxDetectors][MatrixBins];
    
    //Non linearity:
    // Detector response non-linearlity function
    TF1 * NLF;
    
    //Interpolation vectors
    std::vector<Int_t> sign;
    std::vector<Int_t> EnergyIdx;
    std::vector<Double_t> Energy;
    std::vector<Double_t> dEtrue;
    std::vector<Double_t> binScaling;
    std::vector<Double_t> dNdE;
    
    //BCW Model
    Double_t m_bcw_elec_nl_par[5];
    Double_t m_bcw_elec_nl_par_nominal[5];
    Double_t m_bcw_elec_nl_par_error[5];
    
    Double_t m_bcw_positron_nl_e[n_bcw_positron_nl];
    Double_t m_bcw_positron_nl_fac[n_bcw_positron_nl];
    
    Double_t m_rel_escale[MaxDetectors];
    Double_t m_rel_escale_nominal[MaxDetectors];
    Double_t m_rel_escale_error[MaxDetectors];
    Double_t m_rel_eoffset[MaxDetectors];
    TGraph* g_bcw_elec_nl_error[2];
    
    //LBNL Model
    Double_t total_err;
    Double_t m_lbnl_nl_par[3];
    Double_t m_lbnl_nl_par_nominal[3];
    Double_t m_lbnl_nl_par_error[3];
    Double_t m_lbnl_positron_nl_e[300];
    Double_t m_lbnl_positron_nl_fac[300];
    Double_t m_lbnl_positron_nl_err[3][300];
    
    //IHEP Model
    Double_t m_ihep_nl_par[5];
    Double_t m_ihep_nl_par_nominal[5];
    Double_t m_ihep_nl_par_error[5];

    //Unified Model
    Double_t m_unified_nl_par[m_num_unified_nl_pars];
    Double_t m_unified_nl_par_nominal[m_num_unified_nl_pars];
    Double_t m_unified_nl_par_error[m_num_unified_nl_pars];
    TGraph* g_unified_positron_nl;
    TGraph* g_unified_positron_nl_pulls[10];
    
    Double_t m_unified_positron_nl_e[n_unified_nl_points];
    Double_t m_unified_positron_nl_fac[n_unified_nl_points];
    Double_t m_unified_positron_nl_err[10][n_unified_nl_points];
    
    Double_t m_unified_nl_par_covmatrix[4][4];
    Double_t m_unified_nl_par_covmatrix_l[4][4];
    
    // Resolution
    // Detector resolution function
    TF1 * ResoF;
    
    Double_t ResolutionRange; //Range of resolution
    Double_t ResolutionError; // Detector energy resolution parameter
    Double_t ResolutionErrorUncorrelated; // Detector energy resolution parameter
    Double_t ResolutionBias[MaxDetectors]; // Random biases of the detector resolution.
    
    Double_t m_abs_escale;
    Double_t m_abs_escale_nominal;
    Double_t m_abs_escale_error;
    
    Double_t m_abs_eoffset;
    Double_t m_abs_eoffset_error;
    
    Double_t m_rel_eoffset_error;
    
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

    //Choose matrix
    bool RandomSin22t13;
    bool RandomIsotopeFraction;
    bool RandomReactorPower;
    bool VaryAccidentalMatrix;
    bool VaryLiHeMatrix;
    bool VaryFastNeutronsMatrix;
    bool VaryAmCMatrix;
    bool DistortLiHeMatrix;
    bool DistortFastNeutronsMatrix;
    bool DistortAmCMatrix;
    
    //Histograms
    TH1F* FDBH[MaxPeriods];
    TH1F* FLOH[MaxPeriods];
    TH1F* OriginalPredictionH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];//  Spectrum with no reactor alterations. 240 bins.
    TH1F* PredictionTrueH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];// Prediction with no alterations, nominal backgrounds may be added. True Binning.
    TH1F* PredictionVisH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];//  Prediction with no alterations, nominal backgrounds may be added. Vis Binning (after detector response is applied)
    TH1F* CopyPredictionVisH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    TH1F* AlteredPredictionVisH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];// Prediction with background alterations. Vis Binning (after detector response is applied)
    TH1F* CopyAlteredPredictionVisH[MaxFarDetectors][MaxNearDetectors][MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    TH1F* OriginalADSpectrumH[MaxDetectors][MaxPeriods];//  Spectrum that may include reactor alterations if reactor covariance matrices are produced.
    TH1F* ADSpectrumTrueH[MaxDetectors][MaxPeriods];
    TH1F* CopyADSpectrumVisH[MaxDetectors][MaxPeriods];
    TH1F* ADSpectrumVisH[MaxDetectors][MaxPeriods];
    TH1F* ADSpectrumVisWithBkgdH[MaxDetectors][MaxPeriods];
    
    TH1F* BackgroundSpectrumH[MaxDetectors][MaxPeriods];
    TH1F* RandomBackgroundSpectrumH[MaxDetectors][MaxPeriods];
    TH1F* AccidentalsH[MaxDetectors][MaxPeriods];
    TH1F* LiHeH[MaxDetectors][MaxPeriods];
    TH1F* FastNeutronsH[MaxDetectors][MaxPeriods];
    TH1F* AmCH[MaxDetectors][MaxPeriods];
    TH1F* RandomAccidentalsH[MaxDetectors][MaxPeriods];
    TH1F* RandomLiHeH[MaxDetectors][MaxPeriods];
    TH1F* RandomFastNeutronsH[MaxDetectors][MaxPeriods];
    TH1F* RandomAmCH[MaxDetectors][MaxPeriods];
    
    TH2F* EnergyMatrixH[MaxDetectors];
    TH2F* NominalResponseMatrixH[MaxDetectors];
    TH2F* CovMatrix2H;
    TH2F* Cov2H;
    
    Double_t Sigma_Near[MaxFarDetectors][MaxNearDetectors][MaxPeriods][MaxNbins];
    Double_t Sigma_Far[MaxFarDetectors][MaxNearDetectors][MaxPeriods][MaxNbins];
    
    Double_t CovStat[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t NormCovStat[9*MaxNbins][9*MaxNbins][MaxPeriods];
    TH2F* StatCov2H;
    
    Double_t Cov[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t NormCov[9*MaxNbins][9*MaxNbins][MaxPeriods];

    //Main Data loaded from Theta13-inputs_20week.txt if 20 weeks is used. Otherwise change input file.
    Int_t ObsEvts[MaxDetectors][MaxPeriods];
    Double_t Livetime[MaxDetectors][MaxPeriods];
    Double_t MuonVetoEff[MaxDetectors][MaxPeriods];
    Double_t DMCEff[MaxDetectors][MaxPeriods];
    Double_t TargetMass[MaxDetectors][MaxPeriods];
    Double_t BgEvts[MaxDetectors][MaxPeriods];
    Double_t AccEvts[MaxDetectors][MaxPeriods];
    Double_t AccErr[MaxDetectors][MaxPeriods];
    Double_t Li9Evts[MaxDetectors][MaxPeriods];
    Double_t Li9Err[MaxDetectors][MaxPeriods];
    Double_t FnEvts[MaxDetectors][MaxPeriods];
    Double_t FnErr[MaxDetectors][MaxPeriods];
    Double_t AmcEvts[MaxDetectors][MaxPeriods];
    Double_t AmcErr[MaxDetectors][MaxPeriods];
    Double_t AlnEvts[MaxDetectors][MaxPeriods];
    Double_t AlnErr[MaxDetectors][MaxPeriods];
    Double_t ErrEvts[MaxDetectors][MaxPeriods];
    
    //Load functions
    void LoadPredictions();
    void LoadBackgrounds();
    void LoadIavCorrection();
    void LoadNearHall();
    void SetMode();
    
    //Functions to vary background shapes. So far the same ones than LBNL.
    void FluctuateBackgrounds();
    TF1* GetDistortionFunction(Double_t);
    TF1* GetFastNeutronsDistortionFunction(Double_t);
    
    //Functions to calculate the energy matrix
    TF1* VisibleF;
    TF1* GetNeutrinoToVisibleFunction(Int_t);
    void GetEnergyShift(Int_t);
    void GetIAVShift(Int_t);
    void GetNLShift(Int_t);
    void GetResolutionShift(Int_t);
    void GetStatisticalFluctuation(TH1F*);
    void CreateEnergyMatrix();
    void Interpolation(TF1*);

    //Functions to generate the respective covariance matrices
    void GenerateCovarianceMatrix();
    void SaveSpectrum(Int_t);
    void SaveCovarianceMatrix();
    void GenerateStatisticalCovarianceMatrix();
    void SaveStatisticalCovarianceMatrix();
    
    Double_t VisibleEnergy0F(Double_t *, Double_t *);// Zeroth order true to visible function.
    Double_t VisibleEnergy1F(Double_t *, Double_t *);// First order true to visible function.
    
    Double_t NLBCWF(Double_t*, Double_t*);
    Double_t NLLBNLF(Double_t*, Double_t*);
    Double_t NLUnifiedF(Double_t*, Double_t*);

    Double_t ResolutionF(Double_t *, Double_t *);
    
public:
    CovarianceMatrix3();
    CovarianceMatrix3(NominalData*);

    void CovarianceMatrixMain();
    
    void SetVaryAccidentalMatrix(bool);
    void SetVaryLiHeMatrix(bool);
    void SetVaryFastNeutronsMatrix(bool);
    void SetVaryAmCMatrix(bool);
    
    void SetDistortLiHeMatrix(bool);
    void SetDistortFastNeutronsMatrix(bool);
    void SetDistortAmCMatrix(bool);
    
    void SetIsotopeMatrix(bool);
    void SetReactorPowerMatrix(bool);
    
    void SetIAVMatrix(bool);
    void SetNLMatrix(bool);
    void SetResolutionMatrix(bool);

    void SetAbsoluteEnergyScaleMatrix(bool);
    void SetRelativeEnergyScaleMatrix(bool);
    void SetRelativeEnergyOffset(bool);
    void SetAbsoluteEnergyOffset(bool);

    void RandomAbsoluteEnergyScaleMatrix();
    void RandomRelativeEnergyScaleMatrix();
    void RandomRelativeEnergyOffsetMatrix();
    void RandomAbsoluteEnergyOffsetMatrix();
    
    void RandomIAVMatrix();
    void RandomNLMatrix();
    void RandomResolutionMatrix();
    
    void SetBCWModel(bool);
    void SetLBNLModel(bool);
    void SetIHEPModel(bool);
    void SetUnifiedModel(bool);

    void LoadNLParameters();
    void SetNLParameters();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                             Nominal constructor, not used but coded in case it is needed
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CovarianceMatrix3 :: CovarianceMatrix3()
{
    Nom = new NominalData();
    rand = new TRandom3();
    Pred = new Prediction();
    
    SimpleReactorModel = Nom->GetSimpleReactorModel();
    
    NSamples = Nom->GetNSamples();
    Nweeks = Nom->GetWeeks();
    
    InitialEnergy = Nom->GetEmin();
    FinalEnergy = Nom->GetEmax();
    InitialVisibleEnergy = Nom->GetEVisMin();
    FinalVisibleEnergy =  Nom->GetEVisMax();
    
    TotalBins = MatrixBins;
    BinWidth=(FinalVisibleEnergy-InitialVisibleEnergy)/TotalBins;
    
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
        
    NADs = Nom->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    if(NADs == 8)//    ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    //  IAV Error from Bryce
    IAVNominalError=Nom->GetIAVError();
    //  Non uniformity
    m_abs_escale = Nom->GetAbsoluteEnergyScale();
    m_abs_escale_nominal = m_abs_escale;
    m_abs_escale_error = Nom->GetAbsoluteEnergyScaleError();
    
    ResolutionError = Nom->GetResolutionError();
    ResolutionErrorUncorrelated = Nom->GetResoUncorrelatedError();
    ResolutionRange = 8;// Why 8σ? Seems chosen trivially but in my opinion it's the way to limit the range of the convolution, for a Normal distribution this covers up to 99.99999... of the area

    SystematicE=(CovarianceMatrix3::SystematicType)-1;//    A number different to any of the different systematics included in the model
    DistortLiHe= 0;
    DistortFN  = 0;
    DistortAmC = 0;
    
    //  Scaling Backgrounds to the data rates:

    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            LiveTime[AD][week] = Nom->GetLiveTime(AD,week);
            DetectorEfficiency[AD][week] = Nom->GetDetectorEfficiency(AD,week);
            IBDRate[AD][week] = Nom->GetIBDRate(AD,week);
            
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
        //  Relative energy scale
        m_rel_escale[AD] = Nom->GetRelativeEnergyScale(AD);
        m_rel_escale_error[AD] = Nom->GetRelativeEnergyError(AD);
        m_rel_escale_nominal[AD] = m_rel_escale[AD];
        m_rel_eoffset[AD] = Nom->GetRelativeEnergyOffset(AD);
    }    
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Constructor carrying all Data information
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CovarianceMatrix3 :: CovarianceMatrix3(NominalData* Data)
{
    rand = new TRandom3();
    Pred = new Prediction(Data);

    SimpleReactorModel = Data->GetSimpleReactorModel();

    NSamples = Data->GetNSamples();
    Nweeks = Data->GetWeeks();
    
    InitialEnergy = Data->GetEmin();
    FinalEnergy = Data->GetEmax();
    InitialVisibleEnergy = Data->GetEVisMin();
    FinalVisibleEnergy = Data->GetEVisMax();
    TotalBins = MatrixBins;
    BinWidth=(FinalVisibleEnergy-InitialVisibleEnergy)/TotalBins;
    
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
        
    NADs = Data->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;

    if(NADs == 8)// ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    //  IAV Error from Bryce
    IAVNominalError=Data->GetIAVError();
    //  Non uniformity
    m_abs_escale = Data->GetAbsoluteEnergyScale();
    m_abs_escale_nominal = m_abs_escale;
    m_abs_escale_error = Data->GetAbsoluteEnergyScaleError();
    
    ResolutionError = Data->GetResolutionError();
    ResolutionErrorUncorrelated = Data->GetResoUncorrelatedError();
    
    ResolutionRange = 8; // Why 8σ? Seems chosen trivially but in my opinion it's the way to limit the range of the convolution, for a Normal distribution this covers up to 99.99999... of the area
    
    SystematicE=(CovarianceMatrix3::SystematicType)-1;//    A number different to any of the systematics included in the model
    DistortLiHe= 0;
    DistortFN  = 0;
    DistortAmC = 0;
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            LiveTime[AD][week] = Data->GetLiveTime(AD,week);
            DetectorEfficiency[AD][week] = Data->GetDetectorEfficiency(AD,week);
            IBDRate[AD][week] = Data->GetIBDRate(AD,week);

            //Background uncertainties
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
            //  Scaling Backgrounds to the data rates:
            ScaleAcc[AD][week]=Data->GetAccidentalEvents(AD,week);
            ScaleLiHe[AD][week]=Data->GetLiHeEvents(hall,week);
            ScaleFN[AD][week]=Data->GetFNEvents(hall,week);
            ScaleAmC[AD][week]=Data->GetAmCEvents(hall,week);
        }
        //  Relative energy scale
        m_rel_escale[AD] = Data->GetRelativeEnergyScale(AD);
        m_rel_escale_error[AD] = Data->GetRelativeEnergyError(AD);
        m_rel_escale_nominal[AD] = m_rel_escale[AD];
        m_rel_eoffset[AD] = Data->GetRelativeEnergyOffset(AD);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Main function to calculate covariance matrices (Background + Systematics)  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CovarianceMatrix3 :: CovarianceMatrixMain()
{
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
    
    if(SystematicE==NLE)
    {
        TFile* NLMatrix1F = TFile::Open("./CrossChecks/EnergyMatrixComparison/NLalteredResponseMatrix.root","recreate");
        NLMatrix1F->Close();
    }
    if(SystematicE==IAVE)
    {
        TFile* IAVMatrixF = TFile::Open("./CrossChecks/EnergyMatrixComparison/IAValteredResponseMatrix.root","recreate");
        IAVMatrixF->Close();
    }
    if(SystematicE==ResolutionE)
    {
        TFile* ResoMatrix1F = TFile::Open("./CrossChecks/EnergyMatrixComparison/ResoalteredResponseMatrix.root","recreate");
        ResoMatrix1F->Close();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                      Initialize histograms, functions, etc.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Cov2H = new TH2F("Covariance Matrix","Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);

    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t AD = 0; AD < NADs;  AD++)
        {
            ADSpectrumVisH[AD][week] = new TH1F(Form("Nominal Visible AD%d Spectrum",AD),Form("Nominal Visible AD%d Spectrum",AD),n_evis_bins,evis_bins);
            ADSpectrumVisWithBkgdH[AD][week] = new TH1F(Form("Varied Visible AD%d Spectrum",AD),Form("Varied Visible AD%d Spectrum",AD),n_evis_bins,evis_bins);
        }
        
        for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)
        {
            for (Int_t far = 0; far<ADsEH3; far++)
            {
                PredictionVisH[far][near][week] = new TH1F(Form("Far Nominal Visible Spectrum at AD%d from AD%d",far+1,near+1),Form("Far Nominal Visible Spectrum at AD%d from AD%d",far+1,near+1),n_evis_bins,evis_bins);
                AlteredPredictionVisH[far][near][week] = new TH1F(Form("Far Visible Spectrum at AD%d from AD%d with variations due to Systematic%d or Background%d",far+1,near+1,(Int_t)SystematicE,(Int_t)BackgroundE),Form("Far Visible Spectrum at AD%d from AD%d with variations due to Systematic%d or Background%d",far+1,near+1,(Int_t)SystematicE,(Int_t)BackgroundE),n_evis_bins,evis_bins);
                
//                Tree->Branch(Form("Far Nominal Visible Spectrum at AD%d from AD%d",far+1,near+1),"TH1F",&PredictionVisH[far][near][week],32000,0);
//                Tree->Branch(Form("Far Visible Spectrum at AD%d from AD%d with variations due to Systematic%d or Background%d",far+1,near+1,(Int_t)SystematicE,(Int_t)BackgroundE),"TH1F",& AlteredPredictionVisH[far][near][week],32000,0);
            }
        }
    }
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        EnergyMatrixH[AD] = new TH2F(Form("EvisEnu%d",AD),Form("EvisEnu%d",AD),n_etrue_bins,enu_bins,n_evis_bins,evis_bins);
    }
    CovMatrix2H = new TH2F("Covariance Matrix","Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    c = new TCanvas("canvas");
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                             Set kinematic energy shift function
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    VisibleF = GetNeutrinoToVisibleFunction(0);//0 for 0th order, 1 for 1st order

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                              Set nominal resolution parameters
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ResoF = new TF1("ResoF",this,&CovarianceMatrix3::ResolutionF,0,20,3,"CovarianceMatrix3","ResolutionF");
//    ResoF->SetParameters(0.022,0.077,0.018); // based on Bryce's TN
    ResoF ->SetParameters(0.0148,0.0869,0.0271); // based on Doc-8982

    //BCW values http://dayabay.ihep.ac.cn/DocDB/0087/008768/013/6AdAnalysis-BCW.pdf are different! 0.13226034931, 0.32604048828, 0.26196488314

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                   Set nominal NL parameters
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    LoadNLParameters();//nominal parameters
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                         Load spectra
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Pred->MakePrediction();//Generate Nominal Prediction
    LoadPredictions();
    LoadNearHall();//Order matters, this have to be after LoadPredictions()
    LoadBackgrounds();
    LoadIavCorrection();
    
    //Create Nominal Response Matrix:
    SetNLParameters();
    CreateEnergyMatrix();
    TFile* EnergyMatrixDataF1 = TFile::Open("./CrossChecks/EnergyMatrixComparison/NominalResponseMatrix.root","recreate");
   
    for (Int_t AD =0; AD<NADs; AD++)
    {
        NominalResponseMatrixH[AD] = (TH2F*)EnergyMatrixH[AD]->Clone(Form("NominalEvisEnu%i",AD));
        NominalResponseMatrixH[AD]->Write();
    }
     EnergyMatrixDataF1->Close();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                     Generate statistical matrix
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    StatCov2H = new TH2F("Statistical Covariance Matrix","Statistical Covariance Matrix",MaxBins,0,MaxBins,MaxBins,0,MaxBins);
    GenerateStatisticalCovarianceMatrix();
    SaveStatisticalCovarianceMatrix();
    delete StatCov2H;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                     Generate toy MC samples
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (Int_t samples = 0; samples<NSamples; samples++)
    {
        flag=0;//first sample flag

        std::cout << samples << " SAMPLE " << std::endl;
        std::cout << "Randomize Background #" << BackgroundE << std::endl;
        std::cout << "Randomize Systematic #" << SystematicE << std::endl;
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                             Reset to nominal parameters after systematic has been varied
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //When I use the automatic script the SystematicE increases in steps of 1, after finishing with the current systematic variations I load the nominal values again.
        if(SystematicE>IAVE && SystematicE<=(IAVE+1))
        {
            for (Int_t AD = 0; AD < NADs; AD++)
            {
                for(Int_t i=0;i<MatrixBins;i++)
                {
                    for(Int_t j=0;j<MatrixBins;j++)
                    {
                        IAVMatrix[AD][i][j] = NominalIAVMatrix[i][j];//reset
                    }
                }
            }
        }
        
        if(SystematicE>ResolutionE && SystematicE<=(ResolutionE+1))
        {
            for(Int_t i=0;i<NADs;i++)
            {
                ResolutionBias[i] = 0;//reset
            }
        }
        
        SetNLParameters();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                  Vary systematics
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        switch (SystematicE)
        {
            case 0://Vary Reactor
                Pred->SetRandomIsotopeFraction(1);
                Pred->MakePrediction();
                Pred->SetRandomIsotopeFraction(0);//reset after calculation
                break;
            case 1://Vary Reactor
                Pred->SetRandomReactorPower(1);
                Pred->MakePrediction();
                Pred->SetRandomReactorPower(0);//reset after calculation
                break;
            case 2://Vary Relative Energy Scale
                RandomRelativeEnergyScaleMatrix();
                break;
            case 3:
                break;
            case 4:
                break;
            case 5:
                break;
            case 6://Vary IAV
                RandomIAVMatrix();
                break;
            case 7://   Vary NL
                //this case is handle by SetNLParameters();
                break;
            case 8://   Vary Resolution
                RandomResolutionMatrix();
                break;
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                     Create Response Matrix with variations included
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        CreateEnergyMatrix();
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                             Generate Covariance Matrix
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Cov2H->Reset();
        GenerateCovarianceMatrix();
        CovMatrix2H->Add(Cov2H);
        
        SaveSpectrum(samples);
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                 Clean up the dust
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for (Int_t TrueEnergyIndex= 0; TrueEnergyIndex<TotalBins; TrueEnergyIndex++)
        {
            for (Int_t AD = 0; AD<NADs; AD++)
            {
                delete PositronTrueSpectrumH[AD][TrueEnergyIndex];
                delete PositronIAVSpectrumH[AD][TrueEnergyIndex];
                delete PositronNLSpectrumH[AD][TrueEnergyIndex];
                delete PositronVisibleSpectrumH[AD][TrueEnergyIndex];
            }
//            if(SystematicE==NLE || SystematicE==NLE+1)
//            {
//                //If I need to interpolate again I need to empty the vectors first.
//                dEtrue.pop_back();
//                Energy.pop_back();
//                EnergyIdx.pop_back();
//                sign.pop_back();
//            }
        }

    }

    CovMatrix2H->Scale(1./(NSamples));
        // CovMatrix2H->Draw("colz");

    SaveCovarianceMatrix();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                      Clean up the dust
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    delete c;
    delete NLF;
    delete ResoF;
    delete VisibleF;
    delete CovMatrix2H;
    delete Cov2H;
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            delete ADSpectrumTrueH[AD][week];
            delete ADSpectrumVisH[AD][week];
            delete ADSpectrumVisWithBkgdH[AD][week];
        }
        delete EnergyMatrixH[AD];
    }
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)
        {
            for (Int_t far = 0; far<ADsEH3; far++)
            {
                delete PredictionVisH[far][near][week];
                delete AlteredPredictionVisH[far][near][week];
                delete PredictionTrueH[far][near][week];
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Produces the Statistical Covariance Matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CovarianceMatrix3 :: GenerateStatisticalCovarianceMatrix()
{
    //Nominal in this one:
    for (Int_t week = 0; week<Nweeks; week++)
    {
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
            ADSpectrumVisH[AD][week]->Reset();

            for(Int_t i=1; i<=n_evis_bins; i++)
            {
                for(Int_t j=1; j<=n_etrue_bins; j++)
                {
                    //  ADSpectrumVisH should be multiplied by the nominal response matrix, produced with all systematics = 0 if is not data, 
                    ADSpectrumVisH[AD][week]->SetBinContent(i,ADSpectrumVisH[AD][week]->GetBinContent(i) + NominalResponseMatrixH[AD]->GetBinContent(i,j) * ADSpectrumTrueH[AD][week]->GetBinContent(j));
                }
            }
        }
        for (Int_t far=0; far<MaxFar; far++)
        {
            for (Int_t near = 0; near < MaxNear; near++)
            {
                //  Add detector effects //Multiply Matrix and Prediction (Far Spectrum)
                PredictionVisH[far][near][week]->Reset();                
                for(Int_t i=1; i<=n_evis_bins; i++)
                {
                    for(Int_t j=1; j<=n_etrue_bins; j++)
                    {
                        //  Nominal in this one:
                        PredictionVisH[far][near][week]->SetBinContent(i,PredictionVisH[far][near][week]->GetBinContent(i) + NominalResponseMatrixH[far+ADsEH1+ADsEH2]->GetBinContent(i,j) * PredictionTrueH[far][near][week]->GetBinContent(j));
                    }
                }
            }
        }
        
        if (Combine == 1)
        {
            ADSpectrumVisH[0][week]->Add(ADSpectrumVisH[1][week]);//All near hall detectors together
            ADSpectrumVisH[0][week]->Add(ADSpectrumVisH[2][week]);
            
            for (Int_t far=0; far<ADsEH3; far++)
            {
                for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
                {
                    //1x1
                    if(far>0||near>0)
                    {
                        PredictionVisH[0][0][week]->Add(PredictionVisH[far][near][week]);
                    }
                }
            }
        }
        else if(Combine == 2)
        {
            for (Int_t far=0; far<ADsEH3; far++)
            {
                for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
                {
                    //1x1
                    if(far>0||near>0)
                    {
                        if(near<ADsEH1)
                        {
                            std::cout << "DB" << std::endl;
                            PredictionVisH[0][0][week]->Add(PredictionVisH[far][near][week]);
                        }
                        if(near>=ADsEH1)
                        {
                            std::cout << "LO" << std::endl;
                            PredictionVisH[0][1][week]->Add(PredictionVisH[far][near][week]);
                        }
                    }
                }
            }
            PredictionVisH[0][0][week]->Scale(0.5);
            
            ADSpectrumVisH[0][week]->Add(ADSpectrumVisH[1][week]);//DB
            ADSpectrumVisH[0][week]->Scale(0.5);
            ADSpectrumVisH[1][week]=(TH1F*)ADSpectrumVisH[2][week]->Clone();//LO
            
        }

        for (Int_t far=0; far<MaxFar; far++)
        {
            for (Int_t near=0; near<MaxNear; near++)
            {
                for (Int_t pts = 0; pts < n_evis_bins; pts++)
                {
                    Sigma_Far[far][near][week][pts]=sqrt(PredictionVisH[far][near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[far+(ADsEH1+ADsEH2)][week]->GetBinContent(pts+1));
                    
                    //N OBSERVED, F PREDICTED FROM N OBSERVED, BACKGROUNDS PREDICTED (WITHOUT ANY FLUCTUATION)
                    Sigma_Near[far][near][week][pts]=(PredictionVisH[far][near][week]->GetBinContent(pts+1)/ADSpectrumVisH[near][week]->GetBinContent(pts+1))*sqrt(ADSpectrumVisH[near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[near][week]->GetBinContent(pts+1));
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
                            }
                        }
                    }
                }
            }
        }
    x = 0;
    y = 0;
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
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Fi4=fari+1;}
                
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
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            
                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                
                                NormCovStat[x][y][week]=(CovStat[x][y][week])/(sqrt(CovStat[x][x][week]*CovStat[y][y][week]));
                                
                                StatCov2H->SetBinContent(x+1,y+1,NormCovStat[x][y][week]);
                            }
                        }
                    }
                }
            }
        }
    }
    // StatCov2H->Draw("colz");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Save Statistical Covariance Matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CovarianceMatrix3 :: SaveStatisticalCovarianceMatrix()
{
    //Save statistical matrix
    TFile* SaveStatisticalCovMatrixF = TFile::Open("./CovarianceMatrices/StatisticalCovarianceMatrix.root","recreate");
    StatCov2H->Write();
    SaveStatisticalCovMatrixF->Close();
    
    //Save in a txt file
    if(WriteOutput)
    {
        ofstream statf("CovarianceMatrices/StatisticalCovarianceMatrix.txt");
        Int_t x =0;
        Int_t y =0;
        for (Int_t week = 0; week<Nweeks; week++)
        {
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
                    if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Fi2=fari+1;}
                    if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Fi3=fari+1;}
                    if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Fi4=fari+1;}
                    
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
                                    
                                    statf << NormCovStat[x][y][week] << " ";
                                }
                            }
                        }
                    }
                }
            }
            statf << std::endl;
        }
        statf.close();
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                 Produces the Covariance Matrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CovarianceMatrix3 :: GenerateCovarianceMatrix()
{
    TFile* Prueba;
    Char_t optionN[10];

    for (Int_t week = 0; week<Nweeks; week++)
    {
        for(Int_t AD=0; AD<NADs; AD++)
        {
            BackgroundSpectrumH[AD][week]=(TH1F*)AccidentalsH[0][0]->Clone();
            BackgroundSpectrumH[AD][week]->Reset();

            //  BackgroundSpectrumH holds nominal backgrounds, without being distorted
            //  RandomBackgroundSpectrumH holds random backgrounds, either varying the rate or the shape

            FluctuateBackgrounds();//Vary backgrounds
            
            BackgroundSpectrumH[AD][week]->Add(AccidentalsH[AD][week]);
            BackgroundSpectrumH[AD][week]->Add(LiHeH[AD][week]);
            BackgroundSpectrumH[AD][week]->Add(FastNeutronsH[AD][week]);
            BackgroundSpectrumH[AD][week]->Add(AmCH[AD][week]);
            
            RandomBackgroundSpectrumH[AD][week]=(TH1F*)BackgroundSpectrumH[AD][week]->Clone();

            if (VaryAccidentalMatrix||VaryAmCMatrix||VaryFastNeutronsMatrix||VaryLiHeMatrix||DistortFastNeutronsMatrix||DistortLiHeMatrix||DistortAmCMatrix)
            {
                RandomBackgroundSpectrumH[AD][week]->Reset();
                RandomBackgroundSpectrumH[AD][week]->Add(RandomAccidentalsH[AD][week]);
                RandomBackgroundSpectrumH[AD][week]->Add(RandomLiHeH[AD][week]);
                RandomBackgroundSpectrumH[AD][week]->Add(RandomFastNeutronsH[AD][week]);
                RandomBackgroundSpectrumH[AD][week]->Add(RandomAmCH[AD][week]);
            }
            
                      
            //                                              Add detector effects to AD spectra by multiplying the response matrix and AD True Spectrum
            ADSpectrumVisWithBkgdH[AD][week]->Reset();
            //  ADSpectrumVisWithBkgdH is multiplied by the one fluctuated.
        
            for(Int_t i=1; i<=n_evis_bins; i++)
            {
                for(Int_t j=1; j<=n_etrue_bins; j++)
                {
                    //  Fluctuations in this one:
                    ADSpectrumVisWithBkgdH[AD][week]->SetBinContent(i,ADSpectrumVisWithBkgdH[AD][week]->GetBinContent(i) + EnergyMatrixH[AD]->GetBinContent(i,j) * ADSpectrumTrueH[AD][week]->GetBinContent(j));
                }
            }
            
            if (StatisticalFluctuation)
            {
                GetStatisticalFluctuation(ADSpectrumVisH[AD][week]);
            }
        }
        flagX=0;

        for (Int_t far=0; far<MaxFar; far++)
        {
            for (Int_t near = 0; near < MaxNear; near++)
            {
                //  Add detector effects //Multiply Matrix and Prediction (Far Spectrum)
                AlteredPredictionVisH[far][near][week]->Reset();

                for(Int_t i=1; i<=n_evis_bins; i++)
                {
                    for(Int_t j=1; j<=n_etrue_bins; j++)
                    {
                        //  Fluctuations in this one:
                        AlteredPredictionVisH[far][near][week]->SetBinContent(i,AlteredPredictionVisH[far][near][week]->GetBinContent(i) + EnergyMatrixH[far+ADsEH1+ADsEH2]->GetBinContent(i,j) *  PredictionTrueH[far][near][week]->GetBinContent(j));
                    }
                }
                
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
                    PredictionVisH[far][near][week]->Write();
                    AlteredPredictionVisH[far][near][week]->Write();
                Prueba->Close();

                if (StatisticalFluctuation)
                {
                    GetStatisticalFluctuation(AlteredPredictionVisH[far][near][week]);
                }
                flagX=1;
            }
        }
        if (Combine == 1)
        {
            
            ADSpectrumVisWithBkgdH[0][week]->Add(ADSpectrumVisWithBkgdH[1][week]);//All near hall detectors together
            ADSpectrumVisWithBkgdH[0][week]->Add(ADSpectrumVisWithBkgdH[2][week]);
            
            for (Int_t far=0; far<ADsEH3; far++)
            {
                for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
                {
                    //1x1
                    if(far>0||near>0)
                    {
                        AlteredPredictionVisH[0][0][week]->Add(AlteredPredictionVisH[far][near][week]);
                    }
                }
            }
        }
        else if(Combine == 2)
        {
            for (Int_t far=0; far<ADsEH3; far++)
            {
                for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
                {
                    //1x1
                    if(far>0||near>0)
                    {
                        if(near<ADsEH1)
                        {
                            AlteredPredictionVisH[0][0][week]->Add(AlteredPredictionVisH[far][near][week]);
                        }
                        if(near>=ADsEH1)
                        {
                            AlteredPredictionVisH[0][1][week]->Add(AlteredPredictionVisH[far][near][week]);
                        }
                    }
                }
            }
            AlteredPredictionVisH[0][0][week]->Scale(0.5);
            
            ADSpectrumVisWithBkgdH[0][week]->Add(ADSpectrumVisWithBkgdH[1][week]);//DB
            ADSpectrumVisWithBkgdH[0][week]->Scale(0.5);
            ADSpectrumVisWithBkgdH[1][week]=(TH1F*)ADSpectrumVisWithBkgdH[2][week]->Clone();//LO
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
                        
                        Prueba = TFile::Open(Form("./Varied_And_Nominal_Spectrum_At_Systematic%d_and_Background%d.root",(Int_t)SystematicE,(Int_t)BackgroundE),"update");
                        CopyAlteredPredictionVisH[fari][neari][farj][nearj][week]=(TH1F*)AlteredPredictionVisH[fari][neari][week]->Clone(Form("Difference between Fari %d Neari %d and Farj %d Nearj %d",fari,neari,farj,nearj));
                        CopyAlteredPredictionVisH[fari][neari][farj][nearj][week]->Add(PredictionVisH[farj][nearj][week],-1);
                        CopyAlteredPredictionVisH[fari][neari][farj][nearj][week]->Write();
                        Prueba->Close();
                                                                                                                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            //std::cout << "x = " << i << "+ ("<< Ni1*Fi1<< "+"<< Ni2*Fi2 << "+"<< Ni3*Fi3 << "+"<< Ni4*Fi4 << "-1)*41"<< " that is "<< i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*TotalBins << "\n" ;

                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;

                                Cov[x][y][week]=(AlteredPredictionVisH[fari][neari][week]->GetBinContent(i+1)+RandomBackgroundSpectrumH[fari+(ADsEH1+ADsEH2)][week]->GetBinContent(i+1)-PredictionVisH[fari][neari][week]->GetBinContent(i+1)-BackgroundSpectrumH[fari+(ADsEH1+ADsEH2)][week]->GetBinContent(i+1))*(AlteredPredictionVisH[farj][nearj][week]->GetBinContent(j+1)+RandomBackgroundSpectrumH[farj+(ADsEH1+ADsEH2)][week]->GetBinContent(j+1)-PredictionVisH[farj][nearj][week]->GetBinContent(j+1)-BackgroundSpectrumH[farj+(ADsEH1+ADsEH2)][week]->GetBinContent(j+1));
                                
//                                std::cout << "COV VALUE" << Cov[x][y][week] << std::endl;
////////////////////////////////Normilize systematic and then unnormalize when fitting for different sin2t13. Check also NaNs.
//                                if(SystematicE != (CovarianceMatrix3::SystematicType)(-1))
//                                {
//                                    Cov[x][y][week]= Cov[x][y][week]/(PredictionVisH[fari][neari][week]->GetBinContent(i+1)*PredictionVisH[farj][nearj][week]->GetBinContent(j+1));
//                                }
                            }
                        }
                    }
                }
            }
        }

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

                                if((Cov[x][x][week]*Cov[y][y][week])==0)
                                {
                                    NormCov[x][y][week]=1;//To avoid (0/0) nans when bin contents are practically the same; (this happens when bins are empty, in this case the correlation should be 1)
                                    //                                    std::cout << "Norm cov tried to be inf" << std::endl;
                                }
                                else
                                {
                                    NormCov[x][y][week]=(Cov[x][y][week])/(sqrt(Cov[x][x][week]*Cov[y][y][week]));
                                }
                                Cov2H->SetBinContent(x+1,y+1,NormCov[x][y][week]);
                                //                            std::cout << NormCov[x][y][week];

                            }
                            //  std::cout <<NormCov[x][y][week];
                            // std::cout << " " << Cov[x][y][week] << " " <<  Cov[x][x][week] << " " << Cov[y][y][week] << " " << sqrt(Cov[x][x][week]*Cov[y][y][week]) << "\n";
                        }
                    }
                }
            }
        }
    }
        //    Cov2H->Draw("colz");
}

void CovarianceMatrix3 :: RandomIAVMatrix()
{
    for (Int_t AD=0; AD<NADs; AD++)
    {
        IAVError[AD] = (1+IAVNominalError*rand->Gaus(0,1)); //Each AD is varied individually.
        for(Int_t i=0; i<MatrixBins; i++)
        {
            IAVMatrixFrac[AD][i]=IAVError[AD]*NominalIAVMatrixFrac[i];
            IAVMatrix[AD][i][i]= 1-IAVMatrixFrac[AD][i];//Diagonal terms
            for(Int_t j=0; j<MatrixBins; j++)
            {
                if(i!=j)
                {
                    IAVMatrix[AD][i][j]=IAVError[AD]*NominalIAVMatrix[i][j];
                }
                //            std::cout<<"IAV ORIGINAL" << IAVMatrix[i][j] << std::endl;
                
                
                //            std::cout<< "IAV ALTERED" << IAVMatrix[i][j] << std::endl;
            }
        }
    }
}

void CovarianceMatrix3 :: RandomAbsoluteEnergyScaleMatrix()
{
    m_abs_escale = m_abs_escale_nominal + m_abs_escale_error * rand->Gaus(0,1);
}

void CovarianceMatrix3 :: RandomRelativeEnergyScaleMatrix()
{
    for(Int_t AD=0;AD<NADs;AD++)
    {
        Double_t rel_escale_shift = m_rel_escale_error[AD] * rand->Gaus(0,1);
        m_rel_escale[AD] = m_rel_escale_nominal[AD] + rel_escale_shift;
        
        //        // recalculate detection efficiency
        //        m_detectorEfficiency_Ed[AD] = m_detectorEfficiency_Ed_nominal + rel_escale_shift * 0.24; // 0.24 is a magic factor that convert energy scale shift to delayed energy cut efficiency
        //        m_detectorEfficiency[AD]
        //        = pred->tdper[0].DMCEff[AD]
        //        * pred->tdper[0].MuonVetoEff[AD]
        //        * m_detectorEfficiency_Dt
        //        * m_detectorEfficiency_Ep
        //        * m_detectorEfficiency_Ed[AD]
        //        * m_detectorEfficiency_flash
        //        * m_detectorEfficiency_nGd
        //        * m_detectorEfficiency_spill;
    }
}

void CovarianceMatrix3 :: RandomAbsoluteEnergyOffsetMatrix()
{
    m_abs_eoffset =  m_abs_eoffset_error * rand->Gaus(0,1);
}

void CovarianceMatrix3 :: RandomRelativeEnergyOffsetMatrix()
{
    for(Int_t AD=0;AD<NADs;AD++)
    {
        m_rel_eoffset[AD] =  m_rel_eoffset_error * rand->Gaus(0,1);
    }
}

void CovarianceMatrix3 :: RandomResolutionMatrix()
{
    Double_t corr_bias = ResolutionError  * rand->Gaus(0,1);
    
    for(Int_t AD=0;AD<NADs;AD++)
    {
        ResolutionBias[AD] = corr_bias + ResolutionErrorUncorrelated * rand->Gaus(0,1);
    }
}

void CovarianceMatrix3 :: Interpolation(TF1* func)
{
    dNdE.clear();
    dEtrue.clear();
    sign.clear();
    EnergyIdx.clear();
    Energy.clear();
    binScaling.clear();
    
    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
    {
        Energy.push_back(func->GetX(VisibleEnergyIndex*BinWidth));

        EnergyIdx.push_back((Int_t)(Energy[VisibleEnergyIndex]/BinWidth));
//        std::cout << "IAV Energy " << IAVEnergy[VisibleEnergyIndex] << std::endl;
//        std::cout << "IAV Energy Idx " << IAVEnergyIdx[VisibleEnergyIndex] << std::endl;
//        std::cout << " VisibleEnergyIndex" << VisibleEnergyIndex << std::endl;
        
        binScaling.push_back(func->Derivative(Energy[VisibleEnergyIndex]));
        
//        std::cout << "BIN SCALING" << binScaling[VisibleEnergyIndex] << std::endl;
        
        sign.push_back(1);
        dEtrue.push_back(Energy[VisibleEnergyIndex] - EnergyIdx[VisibleEnergyIndex]*BinWidth);
//        std::cout << "dEtrue" << dEtrue[VisibleEnergyIndex] << std::endl;
        if (dEtrue[VisibleEnergyIndex] < 0)
        {
            sign.push_back(-1);
        }
//        std::cout << "sign" << sign[VisibleEnergyIndex] << std::endl;
    }
}

void CovarianceMatrix3 :: CreateEnergyMatrix()
{
    for(Int_t TrueEnergyIndex=0; TrueEnergyIndex<TotalBins; TrueEnergyIndex++)
    {
        GetEnergyShift(TrueEnergyIndex);
        GetIAVShift(TrueEnergyIndex);
        GetNLShift(TrueEnergyIndex);
        GetResolutionShift(TrueEnergyIndex);
    }
    
    TFile* EnergyMatrixDataF = TFile::Open("./CrossChecks/EnergyMatrixComparison/ResponseMatrix.root","recreate");

    // This should be used with the relative oscillation method
    for(Int_t AD=0; AD<NADs; AD++)
    {
        //Fill the Matrix
            EnergyMatrixH[AD]->Reset();
            
            for (Int_t i = 0; i < n_etrue_bins; i++)
            {
                for (Int_t j = 0; j < n_evis_bins; j++)
                {
                    for(Int_t TrueEnergyIndex = Int_t((enu_bins[i])*(MatrixBins/FinalVisibleEnergy)); TrueEnergyIndex<=Int_t((enu_bins[i+1])*(MatrixBins/FinalVisibleEnergy))&&TrueEnergyIndex<MatrixBins; TrueEnergyIndex++)
                    {
                        for(Int_t VisibleBin = Int_t(evis_bins[j]*(MatrixBins/FinalVisibleEnergy)); VisibleBin<=Int_t(evis_bins[j+1]*(MatrixBins/FinalVisibleEnergy)); VisibleBin++)
                        {
                            EnergyMatrixH[AD]->SetBinContent(i+1,j+1,EnergyMatrixH[AD]->GetBinContent(i+1,j+1)+PositronVisibleSpectrumH[AD][TrueEnergyIndex]->GetBinContent(VisibleBin+1));
                        }
                    }
                }
            }
        // UNCOMMENT FOLLOWING LINES TO SAVE ENERGY MATRIX SLICES IN EACH PRODUCTION STEP
//        for(Int_t TrueEnergyIndex = 0; TrueEnergyIndex<MatrixBins; TrueEnergyIndex++)
//        {
//            PositronTrueSpectrumH[AD][TrueEnergyIndex]->Write();
//            PositronIAVSpectrumH[AD][TrueEnergyIndex]->Write();
//            PositronNLSpectrumH[AD][TrueEnergyIndex]->Write();
//            PositronVisibleSpectrumH[AD][TrueEnergyIndex]->Write();
//        }

        EnergyMatrixH[AD]->Write();//Save Matrix before rebinning

        //Normalize the Matrix
        for (Int_t j = 0; j < n_evis_bins; j++)
        {
            Norma[j] = EnergyMatrixH[AD]->ProjectionX("VisibleEnergySlide",j,j+1)->Integral();//Both methods give the same solution, I assume the first is faster.

            for (Int_t i = 0; i < n_etrue_bins; i++)
            {
                if(Norma[j]!=0)
                {
                    EnergyMatrixH[AD]->SetBinContent(i+1,j+1,EnergyMatrixH[AD]->GetBinContent(i+1,j+1)/Norma[j]);//Normalization so Σi E(i,j) = 1; (Σ(x axis) =1)
                    //After checking the normalization the sum doesn't match exactly 1, but 0.999... There's some rounding error of a few 0.X% in the worst case <0.4%
                }
            }
        }
        //Save Matrix after rebinning
        EnergyMatrixH[AD]->Write();
        
        //Code to check correct normalization:
//        Double_t Normx[MaxNbins]={};
//        Double_t Normy[MaxNbins]={};
//        
//        for(Int_t i=0;i<n_etrue_bins;i++)
//        {
//            for(Int_t j=0;j<n_evis_bins;j++)
//            {
//                //Check normalization
//                Normy[i]=Normy[i]+EnergyMatrixH[AD]->GetBinContent(i+1,j+1);
//            }
//
//        }
//        for(Int_t j=0;j<n_evis_bins;j++)
//        {
//            for(Int_t i=0;i<n_etrue_bins;i++)
//            {
//                Normx[j]=Normx[j]+EnergyMatrixH[AD]->GetBinContent(i+1,j+1);
//            }
//        }
//        for(Int_t j=0;j<n_evis_bins;j++)
//        {
//            TH1D *px =EnergyMatrixH[AD]->ProjectionX("x",j,j+1);
//            std::cout << "Norma in X" << Normx[j]<<std::endl;
//            std::cout << "Norma in X " << px->Integral() <<std::endl;;
//        }
//        for(Int_t i=0;i<n_etrue_bins;i++)
//        {
//            TH1D *py =EnergyMatrixH[AD]->ProjectionY("y",i,i+1);
//            std::cout << "Norma in Y" << Normy[i]<<std::endl;
//            std::cout << "Norma in Y " << py->Integral() <<std::endl;;
//        }
    }
    EnergyMatrixDataF->Close();
    for(Int_t AD=0; AD<NADs; AD++)
    {
        if(SystematicE==NLE)
        {
            TFile* NLMatrixF = TFile::Open("./CrossChecks/EnergyMatrixComparison/NLalteredResponseMatrix.root","update");
            EnergyMatrixH[AD]->Write();
            NLMatrixF->Close();
        }
        if(SystematicE==IAVE)
        {
            TFile* IAVMatrix1F = TFile::Open("./CrossChecks/EnergyMatrixComparison/IAValteredResponseMatrix.root","update");
            EnergyMatrixH[AD]->Write();
            IAVMatrix1F->Close();
        }
        if(SystematicE==ResolutionE)
        {
            TFile* ResoMatrix1F = TFile::Open("./CrossChecks/EnergyMatrixComparison/ResoalteredResponseMatrix.root","update");
            EnergyMatrixH[AD]->Write();
            ResoMatrix1F->Close();
        }
    }
}

void CovarianceMatrix3 :: GetEnergyShift(Int_t TrueEnergyIndex)
{
    Double_t NeutrinoEnergy = (TrueEnergyIndex)*BinWidth;
    PositronEnergy = VisibleF->Eval(NeutrinoEnergy);
    PositronEnergyIndex=(Int_t)(PositronEnergy/BinWidth);

    for(Int_t AD=0; AD<NADs; AD++)
    {
        //Calculate Energy Shift
        PositronTrueSpectrumH[AD][TrueEnergyIndex]= new TH1F(Form("PositronTrueSpectrum%d,%d",AD,TrueEnergyIndex),Form("PositronTrueSpectrum%d,%d",AD,TrueEnergyIndex), TotalBins,InitialVisibleEnergy,FinalVisibleEnergy);
        
        //Reset values to 0
        if (PositronEnergy < 1.015)//is about the minimum visible energy (1.8-(939.565-938.272-0.51099) = 1.017999 for 0th order, 1.01574 for 1st order calculation)
        {
            for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
            {
                PositronTrueSpectrumH[AD][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, 0);
            }
        }
        else
        {
            PositronTrueSpectrumH[AD][TrueEnergyIndex]->SetBinContent(PositronEnergyIndex+1,OriginalADSpectrumH[AD][0]->GetBinContent(PositronEnergyIndex-Int_t(1.015*(MatrixBins/FinalVisibleEnergy))+1));
        }
    }
}

void CovarianceMatrix3 :: GetIAVShift(Int_t TrueEnergyIndex)
{
    for(int AD=0; AD<NADs; AD++)
    {
        //IAV
        PositronIAVSpectrumH[AD][TrueEnergyIndex]=(TH1F*)PositronTrueSpectrumH[AD][TrueEnergyIndex]->Clone(Form("IAVSpectrum%d,%d",AD,TrueEnergyIndex));
        PositronIAVSpectrumH[AD][TrueEnergyIndex]->SetTitle("IAV Spectrum");
        //Calculate IAV Shift
        if (PositronEnergy >= 1.015)//1.8-0.78 is about the minimum visible energy (1.8-(939.565-938.272-0.51099) = 1.017999 for 0th order, 1.01574 for 1st order calculation)
        {
            //Reset values to 0
            for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
            {
               PositronIAVSpectrumH[AD][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, 0);
            }
            
            for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TrueEnergyIndex+1; VisibleEnergyIndex++)
            {
                PositronIAVSpectrumH[AD][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, PositronIAVSpectrumH[AD][TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1) + IAVMatrix[AD][PositronEnergyIndex][VisibleEnergyIndex] * PositronTrueSpectrumH[AD][TrueEnergyIndex]->GetBinContent(PositronEnergyIndex+1));
            }
        }
    }
}

void CovarianceMatrix3 :: GetNLShift(Int_t TrueEnergyIndex)
{
    for(Int_t AD=0; AD<NADs; AD++)
    {
        //NL
        PositronNLSpectrumH[AD][TrueEnergyIndex] = (TH1F*)PositronIAVSpectrumH[AD][TrueEnergyIndex]->Clone(Form("NLSpectrum%d,%d",AD,TrueEnergyIndex));
        PositronNLSpectrumH[AD][TrueEnergyIndex]->SetTitle("NL Spectrum");
        
        //Calculate NL Shift
        for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
        {
            if(EnergyIdx[VisibleEnergyIndex]==(TotalBins-1))
            {
                dNdE.push_back(PositronIAVSpectrumH[AD][TrueEnergyIndex]->GetBinContent(EnergyIdx[VisibleEnergyIndex]+1));
            }
            else
            {
                dNdE.push_back((TMath::Abs(dEtrue[VisibleEnergyIndex])/BinWidth)*PositronIAVSpectrumH[AD][TrueEnergyIndex]->GetBinContent(EnergyIdx[VisibleEnergyIndex]+sign[VisibleEnergyIndex]+1)+(1 - TMath::Abs(dEtrue[VisibleEnergyIndex])/BinWidth)* PositronIAVSpectrumH[AD][TrueEnergyIndex]->GetBinContent(EnergyIdx[VisibleEnergyIndex]+1));
            }
        }
        for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
        {
            PositronNLSpectrumH[AD][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, dNdE[VisibleEnergyIndex]/binScaling[VisibleEnergyIndex]);
            dNdE.pop_back();//Without this line I keep pushing them down the stack            
        }
    }
}

void CovarianceMatrix3 :: GetResolutionShift(Int_t TrueEnergyIndex)
{
    for(Int_t AD=0; AD<NADs; AD++)
    {
        //Resolution
        PositronVisibleSpectrumH[AD][TrueEnergyIndex] = (TH1F*)PositronNLSpectrumH[AD][TrueEnergyIndex]->Clone(Form("VisibleSpectrum%d,%d",AD,TrueEnergyIndex));
        PositronVisibleSpectrumH[AD][TrueEnergyIndex]->SetTitle("Visible Spectrum");

        //Calculate Resolution effect
        for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
        {
            Double_t sigma = (ResoF->Eval(VisibleEnergyIndex*BinWidth) + ResolutionBias[AD]) * VisibleEnergyIndex*BinWidth;//Sigma = FReso(Energy) * Energy
            Double_t minDetE = VisibleEnergyIndex*BinWidth - ResolutionRange*sigma;
            Double_t maxDetE = VisibleEnergyIndex*BinWidth + ResolutionRange*sigma;
            Int_t minDetEIdx = (Int_t)((minDetE)/BinWidth);
            Int_t maxDetEIdx = (Int_t)((maxDetE)/BinWidth);
            
            if(minDetEIdx < 0){minDetEIdx = 0;}
            if(maxDetEIdx >= TotalBins){maxDetEIdx = TotalBins-1;}
            
            for(Int_t detIdx=minDetEIdx; detIdx<=maxDetEIdx; detIdx++)
            {
                if(detIdx==0)
                {
                    continue;
                }
                Double_t gausFactor = TMath::Gaus((VisibleEnergyIndex-detIdx)*BinWidth,0,sigma,true);
                
                PositronVisibleSpectrumH[AD][TrueEnergyIndex]->SetBinContent(detIdx+1, (PositronVisibleSpectrumH[AD][TrueEnergyIndex]->GetBinContent(detIdx+1)+PositronNLSpectrumH[AD][TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1)*gausFactor));
            }
        }
    }
}

void CovarianceMatrix3 :: GetStatisticalFluctuation(TH1F* Histo)
{
    for(int VisibleEnergyIndex=1;VisibleEnergyIndex<=Histo->GetXaxis()->GetNbins();VisibleEnergyIndex++)
    {
        Histo->SetBinContent(VisibleEnergyIndex,(Double_t)(rand->Poisson(Histo->GetBinContent(VisibleEnergyIndex))));
    }
}

void CovarianceMatrix3 :: LoadNLParameters()
{
    switch (NLModelE)
    {
        case 0://BCW NL Model
        {
            std::cout << "LOAD BCW NL PARAMETERS" << std::endl;
            ifstream bcw_positron_data("bcw_nl_data/positron.dat");
            for (Int_t i = 0; i < n_bcw_positron_nl; i++)
            {
                bcw_positron_data >> m_bcw_positron_nl_e[i] >> m_bcw_positron_nl_fac[i];
            }
            bcw_positron_data.close();
            
            ifstream bcw_elec_data("bcw_nl_data/par.dat");
            for (Int_t i = 0; i < 3; i++)
            {
                bcw_elec_data >> m_bcw_elec_nl_par_nominal[i]  >> m_bcw_elec_nl_par_error[i];
                m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];
                
//                std::cout <<  m_bcw_elec_nl_par_error[i] << "IS ELEC NL PAR ERROR" << std::endl;
//                std::cout <<  m_bcw_elec_nl_par[i] << "IS ELEC NL PAR" << std::endl;

            }
            bcw_elec_data.close();
            
            for (Int_t i = 3; i < 5; i++)
            { // those describe additional uncertainty for bcw model that are described in Doc-XXXX
                m_bcw_elec_nl_par[i] = 0;
                m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];
                m_bcw_elec_nl_par_error[i] = 1;
            }
            
            TFile *bcw_ele_err_file = TFile::Open("bcw_nl_data/ele_err.root");
            g_bcw_elec_nl_error[0] = (TGraph*)bcw_ele_err_file->Get("g_up")->Clone();
            g_bcw_elec_nl_error[1] = (TGraph*)bcw_ele_err_file->Get("g_down")->Clone();
            bcw_ele_err_file->Close();
        }
            break;
        case 1://LBNL NL Model
        {
            std::cout << "LOAD LBNL NL PARAMETERS" << std::endl;
            ifstream lbnl_positron_data("lbnl_nl_data/lbnl_positron_nl.txt");
            if (!lbnl_positron_data.is_open())
            {
                std::cout << "Error: cannot find LBNL non-linearity curve!!!" << std::endl;
                exit(0);
            }
            for (Int_t i = 0; i < 300; i++)
            {
                lbnl_positron_data >> m_lbnl_positron_nl_e[i] >> m_lbnl_positron_nl_fac[i] >> total_err >> m_lbnl_positron_nl_err[0][i] >> m_lbnl_positron_nl_err[1][i] >> m_lbnl_positron_nl_err[2][i];
            }
            for (Int_t i = 0; i < 3; i++)
            {
                m_lbnl_nl_par[i] = 0;
                m_lbnl_nl_par_nominal[i] = 0;
                m_lbnl_nl_par_error[i] = 1.0;
            }
        }
            break;
        case 2://IHEP NL Model
        {
            
        }
            break;
        case 3://Unified NL Model
        {
            for (Int_t i = 0; i < m_num_unified_nl_pars; i++)
            {
                m_unified_nl_par[i] = 0.0;
                m_unified_nl_par_nominal[i] = 0.0;
                m_unified_nl_par_error[i] = 1.0;
            }
            
            TFile *unified_nl_file = TFile::Open("unified_nl_data/nl_models_final.root");
            
            // Use IHEP I as the nominal model
            g_unified_positron_nl = (TGraph*)unified_nl_file->Get("positron_0")->Clone();
            g_unified_positron_nl_pulls[0] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",1))->Clone();
            g_unified_positron_nl_pulls[1] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",2))->Clone();
            g_unified_positron_nl_pulls[2] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",3))->Clone();
            g_unified_positron_nl_pulls[3] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",4))->Clone();
            
            unified_nl_file->Close();
            
            // Copy into arrays to speed up
            for (Int_t ie = 0; ie < n_unified_nl_points; ie++)
            {
                Double_t e = 1.022 + 0.02 * ie;
                m_unified_positron_nl_e[ie] = e;
                m_unified_positron_nl_fac[ie] = g_unified_positron_nl->Eval(e);
                for (Int_t i = 0; i < m_num_unified_nl_pars; i++)
                {
                    m_unified_positron_nl_err[i][ie] = g_unified_positron_nl_pulls[i]->Eval(e) - m_unified_positron_nl_fac[ie];
                }
            }
        }
            break;
    }
}

void CovarianceMatrix3 :: SetNLParameters()
{
    switch (NLModelE)
    {
        case 0://BCW NL Model
            std::cout << "Using BCW NL Model"<< std::endl;
            if(SystematicE==NLE)
            {   //Randomize BCW nonlinear model parameters
                std::cout << "Randomizing BCW NL parameters" << std::endl;
                
                for (Int_t i = 0; i < 5; i++)
                {
                    m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i] + m_bcw_elec_nl_par_error[i] * rand->Gaus(0,1);
                    //                    std::cout <<  m_bcw_elec_nl_par_nominal[i] << " AND " <<  m_bcw_elec_nl_par_error[i] << "AND" << rand->Gaus(0,1) << "IS" <<  m_bcw_elec_nl_par[i]<< std::endl;
                }
                for (Int_t AD = 0; AD < NADs; AD++)
                {
                    NLF->SetParameters(m_bcw_elec_nl_par[0],m_bcw_elec_nl_par[1],m_bcw_elec_nl_par[2], m_bcw_elec_nl_par[3], m_bcw_elec_nl_par[4], m_abs_escale * m_rel_escale[AD], m_abs_eoffset+m_rel_eoffset[AD]);
                }
            }
            else if (flag==0 && SystematicE!=NLE)
            {
                for (Int_t i = 0; i < 5; i++)
                {
                    m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];//reset
                }
                for (Int_t AD = 0; AD < NADs; AD++)
                {
                    NLF->SetParameters(m_bcw_elec_nl_par[0],m_bcw_elec_nl_par[1],m_bcw_elec_nl_par[2], m_bcw_elec_nl_par[3], m_bcw_elec_nl_par[4], m_abs_escale * m_rel_escale[AD], m_abs_eoffset+m_rel_eoffset[AD]);
                }
                flag=1;
            }
            break;
        case 1://LBNL NL Model
            std::cout << "Using LBNL NL Model"<< std::endl;
            if(SystematicE==NLE)
            {   //Randomize LBNL nonlinear model parameters
                std::cout << "Randomizing LBNL NL parameters" << std::endl;
                for (Int_t i = 0; i < 3; i++)
                {
                    m_lbnl_nl_par[i] = m_lbnl_nl_par_nominal[i] + m_lbnl_nl_par_error[i] * rand->Gaus(0,1);//electronics
                    //                std::cout << m_lbnl_nl_par_nominal[i] << " AND " <<  m_lbnl_nl_par_error[i] << "AND" << rand->Gaus(0,1) << "IS" << m_lbnl_nl_par[i] << std::endl;
                }
                for (Int_t AD = 0; AD < NADs; AD++)
                {
                    NLF->SetParameters(m_lbnl_nl_par[0],m_lbnl_nl_par[1],m_lbnl_nl_par[2],m_abs_escale * m_rel_escale[AD],m_abs_eoffset + m_rel_eoffset[AD]);
                }
            }
            else if (flag==0 && SystematicE!=NLE)
            {    
                for (Int_t i = 0; i < 3; i++)
                {
                    m_lbnl_nl_par[i] = m_lbnl_nl_par_nominal[i];//reset
                }
                for (Int_t AD = 0; AD < NADs; AD++)
                {
                    NLF->SetParameters(m_lbnl_nl_par[0],m_lbnl_nl_par[1],m_lbnl_nl_par[2],m_abs_escale * m_rel_escale[AD],m_abs_eoffset + m_rel_eoffset[AD]);
                }
                flag=1;
            }
            break;
        case 2:
            if(SystematicE==NLE)
            {   //Randomize IHEP nonlinear model parameters
                std::cout << "Randomizing IHEP NL parameters" << std::endl;

            }
            else if (flag==0 && SystematicE!=NLE)
            {
            }
            break;
        case 3://Unified NL Model
            std::cout << "Using Unified NL Model"<< std::endl;
            if(SystematicE==NLE)
            {   //Randomize Unified nonlinear model parameters
                std::cout << "Randomizing Unified NL parameters" << std::endl;

            }
            else if (flag==0 && SystematicE!=NLE)
            {
                for (Int_t i = 0; i < m_num_unified_nl_pars; i++)
                {
                    NLF->SetParameter(i,m_unified_nl_par[i]);
                }
                
                for (Int_t AD = 0; AD < NADs; AD++)
                {
                    NLF->SetParameter(m_num_unified_nl_pars, m_abs_escale * m_rel_escale[AD]);
                    NLF->SetParameter(m_num_unified_nl_pars+1, m_abs_eoffset + m_rel_eoffset[AD]);
                }
            }
            //

            
            //                //Read covariance matrix
            //                TMatrixD covmatrix_unified(4,4);
            //
            //                TFile* f_covmatrix = new TFile(data->getString("NonlinearCovmatrixFilename"));
            //                f_covmatrix->ls();
            //                Double_t* mat_tmp = ((TMatrixD*)f_covmatrix->Get("coeffmatrix"))->GetMatrixArray();
            //                for (Int_t i = 0; i < 4; i++)
            //                {
            //                    for (Int_t j = 0; j < 4; j++)
            //                    {
            //                        m_unified_nl_par_covmatrix[i][j] =  mat_tmp[i*4+j];
            //                    }
            //                }
            //                else
            //                {
            //                    for (Int_t i = 0; i < 4; i++)
            //                    {
            //                        for (Int_t j = 0; j < 4; j++)
            //                        {
            //                            if (i == j)
            //                                m_unified_nl_par_covmatrix[i][j] = 1;
            //                            else
            //                                m_unified_nl_par_covmatrix[i][j] = 0;
            //                        }
            //                    }
            //                }
            //
            //                covmatrix_unified.SetMatrixArray(&m_unified_nl_par_covmatrix[0][0]);
            //                std::cout << "Covariance matrix for UNIFIED non-liner parameters:" << std::endl;
            //                covmatrix_unified.Print();
            //
            //                TDecompChol chol_unified(covmatrix_unified);
            //                chol_unified.Decompose();
            //
            //                TMatrixD cmat_unified(chol_unified.GetU());
            //                TMatrixD tcmat_unified(cmat_unified.Transpose(cmat_unified));
            //
            //                Double_t * tmp_matrix_unified = tcmat_unified.GetMatrixArray();
            //
            //                for (Int_t i = 0; i < 4; i++)
            //                {
            //                    for (Int_t j = 0; j < 4; j++)
            //                    {
            //                        m_unified_nl_par_covmatrix_l[i][j] = tmp_matrix_unified[i*4 + j];
            //                    }
            //                }
            
            //            //Randomize Unified NL model
            //            Double_t ranvec_uniform[4];
            //            Double_t ranvec[4];
            //
            //            for (Int_t i = 0; i < 4; i++)
            //            {
            //                ranvec_uniform[i] = rand->Gaus(0,1);
            //            }
            //            for (Int_t i = 0; i < 4; i++)
            //            {
            //                ranvec[i] = 0;
            //                for (Int_t j = 0; j < 4; j++)
            //                {
            //                    ranvec[i] += m_unified_nl_par_covmatrix_l[i][j]*ranvec_uniform[j];
            //                }
            //            }
            //            for (Int_t i = 0; i < 4; i++)
            //            {
            //                m_unified_nl_par[i] = m_unified_nl_par_nominal[i] + m_unified_nl_par_error[i] * ranvec[i];
            //                
            //            }
    }
    Interpolation(NLF);//Interpolate the new function
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
}

//Randomize the events (rate) inside errors taking into account the corresponding correlations between backgrounds and ADs. Also random shape variations are included using distortion functions.
void CovarianceMatrix3 :: FluctuateBackgrounds()
{
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t AD =0; AD<NADs; AD++)
        {
            ScaleFactorAccidental[AD][week]=1.;
            ScaleFactorLiHe[AD][week]=1.;
            ScaleFactorFastNeutrons[AD][week]=1.;
            ScaleFactorAmC[AD][week]=1.;
            
            if(VaryAccidentalMatrix)
            {
//                rand->SetSeed(0);
                ScaleFactorAccidental[AD][week]=(1.+(AccidentalError[AD][week]*rand->Gaus(0,1)));
                std::cout << "Accidental error " << AccidentalError[AD][week] << std::endl;
            }
            if(VaryLiHeMatrix)
            {
//                rand->SetSeed(0);
                ScaleFactorLiHe[AD][week]=(1.+(LiHeError[AD][week]*rand->Gaus(0,1)));
                std::cout <<"LiHe error " << LiHeError[AD][week] << std::endl;
            }
            if(VaryFastNeutronsMatrix)
            {
//                rand->SetSeed(0);
                ScaleFactorFastNeutrons[AD][week]=(1.+(FastNeutronsError[AD][week]*rand->Gaus(0,1)));
                std::cout <<"Fast Neutron error " << FastNeutronsError[AD][week] << std::endl;
            }
            if(VaryAmCMatrix)
            {
//                rand->SetSeed(0);
                ScaleFactorAmC[AD][week]=(1.+(AmCError[AD][week]*rand->Gaus(0,1)));
                std::cout <<"AmC error " << AmCError[AD][week] << std::endl;
            }
        }
        
        //Accidentals uncorrelated
        if(VaryLiHeMatrix)
        {   //  LiHe correlated by hall
            ScaleFactorLiHe[1][week]=ScaleFactorLiHe[0][week];// EH1
            ScaleFactorLiHe[4][week]=ScaleFactorLiHe[3][week];// EH3
            ScaleFactorLiHe[5][week]=ScaleFactorLiHe[3][week];// EH3
        }
        if(VaryFastNeutronsMatrix)
        {
            //  Fast neutrons by hall
            ScaleFactorFastNeutrons[1][week]=ScaleFactorFastNeutrons[0][week];// EH1
            ScaleFactorFastNeutrons[4][week]=ScaleFactorFastNeutrons[3][week];// EH3
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
            RandomAccidentalsH[AD][week]=(TH1F*)AccidentalsH[AD][week]->Clone();
            RandomLiHeH[AD][week]=(TH1F*)LiHeH[AD][week]->Clone();
            RandomFastNeutronsH[AD][week]=(TH1F*)FastNeutronsH[AD][week]->Clone();
            RandomAmCH[AD][week]=(TH1F*)AmCH[AD][week]->Clone();

            if(VaryAccidentalMatrix)
            {
                RandomAccidentalsH[AD][week]->Scale(ScaleFactorAccidental[AD][week]);
            }
            if(VaryLiHeMatrix)
            {
                RandomLiHeH[AD][week]->Scale(ScaleFactorLiHe[AD][week]);
            }
            if(VaryFastNeutronsMatrix)
            {
                RandomFastNeutronsH[AD][week]->Scale(ScaleFactorFastNeutrons[AD][week]);
            }
            if(VaryAmCMatrix)
            {
                RandomAmCH[AD][week]->Scale(ScaleFactorAmC[AD][week]);
            }
//                std::cout << ScaleFactorAccidental[AD][week] << " " << ScaleFactorLiHe[AD][week] <<" " <<  ScaleFactorFastNeutrons[AD][week] <<" " <<  ScaleFactorAmC[AD][week] << std::endl;
//                RandomAmCH[AD][week]->Draw("same");
//                RandomLiHeH[AD][week]->Draw("same");
//                RandomFastNeutronsH[AD][week]->Draw("same");
//                RandomAccidentalsH[AD][week]->Draw("same");
        }
        if(DistortLiHeMatrix)// Distort LiHe shape for all ADs in the same way
        {
            TF1* func_LiHe = GetDistortionFunction(DistortLiHe);
            
            for(Int_t AD=0;AD<NADs;AD++)
            {
                RandomLiHeH[AD][week]=(TH1F*)LiHeH[AD][week]->Clone();
                RandomLiHeH[AD][week]->Multiply(func_LiHe);
                //                    if(LiHeH[AD][week]->Integral()!=0)//REMOVE THIS LINE AFTER LIHE IS INCLUDED IN H ANALYSIS
                //                    {
                RandomLiHeH[AD][week]->Scale(LiHeH[AD][week]->Integral()/RandomLiHeH[AD][week]->Integral());
                //                    }
                //                      std::cout<< "Integral division " << (LiHeH[AD][week]->Integral()/RandomLiHeH[AD][week]->Integral())<<"\n";
                //  func_LiHe->Draw("same");
            }
            delete func_LiHe;
        }
        if(DistortFastNeutronsMatrix)//Distort FN shape for all ADs in the same hall
        {
            //  FN shape distortions are applied taking into account correlations in each EH
            TF1* func_FN = GetFastNeutronsDistortionFunction(DistortFN);
            
            for(Int_t AD=0;AD<ADsEH1;AD++)
            {
                RandomFastNeutronsH[AD][week]=(TH1F*)FastNeutronsH[AD][week]->Clone();
                RandomFastNeutronsH[AD][week]->Multiply(func_FN);
                RandomFastNeutronsH[AD][week]->Scale(FastNeutronsH[AD][week]->Integral()/RandomFastNeutronsH[AD][week]->Integral());
            }
            func_FN=GetFastNeutronsDistortionFunction(DistortFN);
            
            for(Int_t AD=ADsEH1;AD<(ADsEH1+ADsEH2);AD++)
            {
                RandomFastNeutronsH[AD][week]=(TH1F*)FastNeutronsH[AD][week]->Clone();
                RandomFastNeutronsH[AD][week]->Multiply(func_FN);
                RandomFastNeutronsH[AD][week]->Scale(FastNeutronsH[AD][week]->Integral()/RandomFastNeutronsH[AD][week]->Integral());
            }
            func_FN=GetFastNeutronsDistortionFunction(DistortFN);
            
            for(Int_t AD=(ADsEH1+ADsEH2);AD<(ADsEH1+ADsEH2+ADsEH3);AD++)
            {
                RandomFastNeutronsH[AD][week]=(TH1F*)FastNeutronsH[AD][week]->Clone();
                RandomFastNeutronsH[AD][week]->Multiply(func_FN);
                RandomFastNeutronsH[AD][week]->Scale(FastNeutronsH[AD][week]->Integral()/RandomFastNeutronsH[AD][week]->Integral());
            }
            delete func_FN;
        }
        if(DistortAmCMatrix)//  Distort AmC shape for all ADs in the same way
        {
            TF1* func_AmC = GetDistortionFunction(DistortAmC);
            
            for(Int_t AD=0;AD<NADs;AD++)
            {
                RandomAmCH[AD][week]=(TH1F*)AmCH[AD][week]->Clone();
                RandomAmCH[AD][week]->Multiply(func_AmC);
                RandomAmCH[AD][week]->Scale(AmCH[AD][week]->Integral()/RandomAmCH[AD][week]->Integral());
                //  RandomAmCH[AD][week]->Draw("same");
                //  func_AmC->Draw("same");
            }
            delete func_AmC;
        }
    }
}

void CovarianceMatrix3 :: LoadPredictions()
{
    std::cout << "Loading predictions" << std::endl;

    TFile* FarHallPredictionsF = TFile::Open("./RootOutputs/NominalFarSpectrumFraction.root");// If reactor model is upgraded we need to produce this file again.
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        FDBH[week]=(TH1F*)gDirectory->Get("Far Prediction From DayaBay");
        FLOH[week]=(TH1F*)gDirectory->Get("Far Prediction From Ling Ao");
        
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            for (Int_t far =0; far<ADsEH3; far++)
            {
                OriginalPredictionH[far][near][week] = (TH1F*)gDirectory->Get(Form("AD%i Far Spectrum prediction from near AD%i",far+1,near+1));
                

                PredictionTrueH[far][near][week] = new TH1F(Form("AD%i Far Spectrum from near AD%i",far+1,near+1),Form("AD%i Far Spectrum from near AD%i",far+1,near+1),n_etrue_bins,enu_bins);

                if(SimpleReactorModel)
                {
                    PredictionTrueH[far][near][week]=(TH1F*)OriginalPredictionH[far][near][week]->Rebin(n_etrue_bins,Form("AD%i Far Spectrum from near AD%i",far+1,near+1),enu_bins);
                }
                else
                {
                    PredictionTrueH[far][near][week]=(TH1F*)OriginalPredictionH[far][near][week]->Clone();
                }
            }
        }
    }
    FarHallPredictionsF->Close();
    
    std::cout << "Finished Loading predictions" << std::endl;

}

void CovarianceMatrix3 :: LoadNearHall()
{
    std::cout << "Loading Near Hall" << std::endl;

    Char_t filenameAD[100];
    if(SystematicE==Sin22t12E)
    {
        sprintf(filenameAD,"./RootOutputs/RandomOutputs/RandomOscillation.root");
    }
    else if(SystematicE == IsotopeE)
    {
    sprintf(filenameAD,"./RootOutputs/RandomOutputs/RandomIsotopeOscillation.root");
    }
    else if(SystematicE == PowerE)
    {
    sprintf(filenameAD,"./RootOutputs/RandomOutputs/RandomPowerOscillation.root");
    }
    else
    {
        sprintf(filenameAD,"./RootOutputs/NominalOutputs/Oscillation.root");
    }
    
    TFile* ADSpectrumF = TFile::Open(filenameAD);
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t AD = 0; AD<NADs; AD++)
        {
            sprintf(filenameAD,"Total spectrum after oscillation at AD%i",AD+1);
            ADSpectrumF->cd("Total AD Spectra after oscillation");
            OriginalADSpectrumH[AD][week] = (TH1F*)gDirectory->Get(filenameAD);
            ADSpectrumTrueH[AD][week] = new TH1F(Form("Spectrum at AD%i",AD+1),Form("Spectrum at AD%i",AD+1),n_etrue_bins,enu_bins);
            if(SimpleReactorModel)
            {
                ADSpectrumTrueH[AD][week]=(TH1F*)OriginalADSpectrumH[AD][week]->Rebin(n_etrue_bins,Form("True Spectrum at AD%i",AD+1),enu_bins);
                ADSpectrumTrueH[AD][week]->Scale(DetectorEfficiency[AD][week]*LiveTime[AD][week]*IBDRate[AD][week]/ADSpectrumTrueH[AD][week]->Integral());//Normalize to expected real data rate. This shouldn't be used for real data inputs, or for oscillationreactor if the reactor model is accuratelly constructed. Check if IBDRate is already corrected for AD efficiencies, if it is I have to remove the following contribution: DetectorEfficiency[AD][week]
            }
            else
            {
                ADSpectrumTrueH[AD][week]=(TH1F*)OriginalADSpectrumH[AD][week]->Clone();
            }
        }
    }
    ADSpectrumF->Close();
    
    TFile* SpectrumF = TFile::Open("./RootOutputs/TrueSpectrum.root","recreate");
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t AD = 0; AD<NADs; AD++)
        {
            ADSpectrumTrueH[AD][week]->Write();
        }
        for (Int_t near =0; near<MaxNear; near++)
        {            
            for (Int_t far =0; far<MaxFar; far++)
            {
                PredictionTrueH[far][near][week]->Write();
            }
        }
    }
    SpectrumF->Close();
    
    std::cout << "Finished Loading Near Hall" << std::endl;

}

void CovarianceMatrix3 :: LoadIavCorrection()// From Bryce Littlejohn's results. //I don't use different IAV matrix for each AD since the analysis it's been done for just 1 of them, assume identical ADs.
{
    TFile* f = TFile::Open("./IavDistortion/IAVDistortion.root");
    TH2F* Correction = (TH2F*)f->Get("Correction");
    
    std::cout << "Reading IAV correction file" << std::endl;
    
    for(Int_t i=0;i<TotalBins;i++)
    { // i: true positron energy bin; j: distorted energy bin
        //Total events in input spectrum bin i
        
        NominalIAVMatrixFrac[i]=0;
        if (Correction->Integral(i+1,i+1,0,-1) > 0)
        {
            for(Int_t j=0;j<i+1;j++)
            {
                NominalIAVMatrix[i][j]= Correction->GetBinContent(i+1,j+1)/Correction->Integral(i+1,i+1,0,-1);
                if(i!=j)
                {
                     NominalIAVMatrixFrac[i] += NominalIAVMatrix[i][j];
                }
            }
        }
        else
        {
            for(Int_t j=0;j<i+1;j++)
            {
                if (i==j)
                {
                    NominalIAVMatrix[i][j] = 1;
                }
                else
                {
                    NominalIAVMatrix[i][j] = 0;
                }
            }
            NominalIAVMatrixFrac[i]=0;
        }
    }
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for(Int_t i=0;i<TotalBins;i++)
        {
            IAVMatrixFrac[AD][i] = NominalIAVMatrixFrac[i];//   Copy that will be varied when calculated random IAV matrix.
            
            for(Int_t j=0;j<TotalBins;j++)
            {
                IAVMatrix[AD][i][j] = NominalIAVMatrix[i][j];// Copy that will be varied when calculated random IAV matrix.
            }
        }
    }

    f->Close();
}

void CovarianceMatrix3 :: SaveCovarianceMatrix()
{
    Char_t filenameCov[100];
    
    sprintf(filenameCov,"./CovarianceMatrices/NominalCovarianceMatrix.root");

    //Save Cov Matrix (Either background or systematics)
    switch (BackgroundE)
    {
        case 0://Vary Accidentals
            sprintf(filenameCov,"./CovarianceMatrices/VaryAccidentalCovarianceMatrix.root");
            CovMatrix2H->SetName("Vary Accidental Covariance Matrix");
            CovMatrix2H->SetTitle("Vary Accidental Covariance Matrix");
            break;
        case 1://Vary LiHe
            sprintf(filenameCov,"./CovarianceMatrices/VaryLiHeCovarianceMatrix.root");
            CovMatrix2H->SetName("Vary LiHe Covariance Matrix");
            CovMatrix2H->SetTitle("Vary LiHe Covariance Matrix");
            break;
        case 2://Vary Fast Neutrons
            sprintf(filenameCov,"./CovarianceMatrices/VaryFastNeutronsCovarianceMatrix.root");
            CovMatrix2H->SetName("Vary FN Covariance Matrix");
            CovMatrix2H->SetTitle("Vary FN Covariance Matrix");
            break;
        case 3://Vary AmC
            sprintf(filenameCov,"./CovarianceMatrices/VaryAmCCovarianceMatrix.root");
            CovMatrix2H->SetName("Vary AmC Covariance Matrix");
            CovMatrix2H->SetTitle("Vary AmC Covariance Matrix");
            break;
        case 4://Distort LiHe
            sprintf(filenameCov,"./CovarianceMatrices/DistortLiHeCovarianceMatrix.root");
            CovMatrix2H->SetName("Distort LiHe Covariance Matrix");
            CovMatrix2H->SetTitle("Distort LiHe Covariance Matrix");
            break;
        case 5://Distort Fast Neutrons
            sprintf(filenameCov,"./CovarianceMatrices/DistortFastNeutronsCovarianceMatrix.root");
            CovMatrix2H->SetName("Distort FN Covariance Matrix");
            CovMatrix2H->SetTitle("Distort FN Covariance Matrix");
            break;
        case 6://Distort AmC
            sprintf(filenameCov,"./CovarianceMatrices/DistortAmCCovarianceMatrix.root");
            CovMatrix2H->SetName("Distort AmC Covariance Matrix");
            CovMatrix2H->SetTitle("Distort AmC Covariance Matrix");
            break;
    }
    switch (SystematicE)
    {
        case 0://Vary Reactor
            sprintf(filenameCov,"./CovarianceMatrices/IsotopeCovarianceMatrix.root");
            CovMatrix2H->SetName("Isotope Covariance Matrix");
            CovMatrix2H->SetTitle("Isotope Covariance Matrix");
            break;
        case 1://Vary Reactor
            sprintf(filenameCov,"./CovarianceMatrices/ReactorPowerCovarianceMatrix.root");
            CovMatrix2H->SetName("Reactor Power Covariance Matrix");
            CovMatrix2H->SetTitle("Reactor Power Covariance Matrix");
            break;
        case 2://Vary Energy Scale
            sprintf(filenameCov,"./CovarianceMatrices/RelativeEnergyScaleCovarianceMatrix.root");
            CovMatrix2H->SetName("Relative Energy Scale Covariance Matrix");
            CovMatrix2H->SetTitle("Relative Energy Scale Covariance Matrix");
            break;
        case 3:
            sprintf(filenameCov,"./CovarianceMatrices/AbsoluteEnergyScaleCovarianceMatrix.root");
            CovMatrix2H->SetName("Absolute Energy Scale Covariance Matrix");
            CovMatrix2H->SetTitle("Absolute Energy Scale Covariance Matrix");
            break;
        case 4:
            sprintf(filenameCov,"./CovarianceMatrices/AbsoluteEnergyOffsetCovarianceMatrix.root");
            CovMatrix2H->SetName("Absolute Energy Scale Covariance Matrix");
            CovMatrix2H->SetTitle("Absolute Energy Scale Covariance Matrix");
            break;
        case 5:
            sprintf(filenameCov,"./CovarianceMatrices/AbsoluteEnergyOffsetCovarianceMatrix.root");
            CovMatrix2H->SetName("Absolute Energy Scale Covariance Matrix");
            CovMatrix2H->SetTitle("Absolute Energy Scale Covariance Matrix");
            break;
        case 6://Vary IAV
            sprintf(filenameCov,"./CovarianceMatrices/IAVCovarianceMatrix.root");
            CovMatrix2H->SetName("IAV Covariance Matrix");
            CovMatrix2H->SetTitle("IAV Covariance Matrix");
            break;
        case 7://Vary NL
            sprintf(filenameCov,"./CovarianceMatrices/NLCovarianceMatrix.root");
            CovMatrix2H->SetName("NL Covariance Matrix");
            CovMatrix2H->SetTitle("NL Covariance Matrix");
            break;
        case 8://Vary Resolution
            sprintf(filenameCov,"./CovarianceMatrices/ResolutionCovarianceMatrix.root");
            CovMatrix2H->SetName("Resolution Covariance Matrix");
            CovMatrix2H->SetTitle("Resolution Covariance Matrix");
            break;
    }
    TFile* SaveCovarianceMatrixF = TFile::Open(filenameCov,"recreate");
    // CovMatrix2H->Draw("colz");
    CovMatrix2H->Write();
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

void CovarianceMatrix3 :: SaveSpectrum(Int_t sample)
{
    Char_t filenameSpec[100];
    Char_t optionS[10];

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

    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t AD = 0; AD<NADs; AD++)
        {
            //            ADSpectrumVisH[AD][week]->Add(BackgroundSpectrumH[AD][week],-1);
            //            ADSpectrumVisWithBkgdH[AD][week]->Add(RandomBackgroundSpectrumH[AD][week],-1);
            ADSpectrumVisH[AD][week]->Write(Form("Nominal AD Vis Spectra Sample%i",sample));//
            ADSpectrumVisWithBkgdH[AD][week]->Write(Form("Varied AD Vis Spectra Sample%i",sample));//Form("Near Vis H Varied Sample%i",sample)
        }
        
        for (Int_t near = 0; near < MaxNear; near++)
        {
            for (Int_t far =0; far<MaxFar; far++)
            {
                //                PredictionVisH[far][near][week]->Add(BackgroundSpectrumH[far+ADsEH1+ADsEH2][week],-1);
                //                AlteredPredictionVisH[far][near][week]->Add(RandomBackgroundSpectrumH[far+ADsEH1+ADsEH2][week],-1);
                
                PredictionVisH[far][near][week]->Write(Form("Nominal Far VisH Sample%i",sample));
                AlteredPredictionVisH[far][near][week]->Write(Form("Varied Far VisH Sample%i",sample));
                
                //Ratio of the variation
                CopyPredictionVisH[far][near][week]=(TH1F*)PredictionVisH[far][near][week]->Clone("Percentual variation");
                CopyPredictionVisH[far][near][week]->Add(AlteredPredictionVisH[far][near][week],-1);
                CopyPredictionVisH[far][near][week]->Divide(AlteredPredictionVisH[far][near][week]);//(Toy-Nominal)/Toy
                CopyPredictionVisH[far][near][week]->Draw("E1");
                CopyPredictionVisH[far][near][week]->Draw("C SAME");
            }
        }
    }
    
    if(sample==5)
    {
        c->Write(); //here write canvas with 20 samples
    }

    SpectrumAndPredictionsF->Close();
}

void CovarianceMatrix3 :: SetVaryAccidentalMatrix(bool VaryAccMatrix)
{
    BackgroundE = (CovarianceMatrix3::BackgroundType)(-1);//Reset BackgroundE. I put it here because VaryAccidentalMatrix is the first method to be called in the "automatic script"
    VaryAccidentalMatrix=VaryAccMatrix;
    if(VaryAccidentalMatrix)
    {
        BackgroundE = VaryAccidentalE;
    }
}

void CovarianceMatrix3 :: SetVaryLiHeMatrix(bool VaryLiMatrix)
{
    VaryLiHeMatrix=VaryLiMatrix;
    if(VaryLiHeMatrix)
    {
        BackgroundE = VaryLiHeE;
    }
}

void CovarianceMatrix3 :: SetVaryFastNeutronsMatrix(bool VaryFNmatrix)
{
  
    VaryFastNeutronsMatrix=VaryFNmatrix;
    if(VaryFastNeutronsMatrix)
    {
        BackgroundE = VaryFastNeutronsE;
    }
}

void CovarianceMatrix3 :: SetVaryAmCMatrix(bool VaryAmCmatrix)
{
    VaryAmCMatrix=VaryAmCmatrix;
    if(VaryAmCMatrix)
    {
        BackgroundE = VaryAmCE;
    }
}

void CovarianceMatrix3 :: SetDistortLiHeMatrix(bool DistortLiMatrix)
{
    DistortLiHeMatrix=DistortLiMatrix;
    if(DistortLiHeMatrix)
    {
        BackgroundE = DistortLiHeE;
        DistortLiHe=0.2;
    }
}

void CovarianceMatrix3 :: SetDistortFastNeutronsMatrix(bool DistortFNmatrix)
{
    DistortFastNeutronsMatrix=DistortFNmatrix;
    if(DistortFastNeutronsMatrix)
    {
        BackgroundE = DistortFastNeutronsE;
        DistortFN=0.2;
    }
}

void CovarianceMatrix3 :: SetDistortAmCMatrix(bool DistortAmCmatrix)
{
    DistortAmCMatrix=DistortAmCmatrix;
    if(DistortAmCMatrix)
    {
        BackgroundE = DistortAmCE;
        DistortAmC=0.2;
    }
}

void CovarianceMatrix3 :: SetIsotopeMatrix(bool IsotopeMatrix)
{    
    if(IsotopeMatrix)
    {
        SystematicE = IsotopeE;
    }
}

void CovarianceMatrix3 :: SetReactorPowerMatrix(bool ReactorPowerMatrix)
{
    if(ReactorPowerMatrix)
    {
        SystematicE = PowerE;
    }
}

void CovarianceMatrix3 :: SetRelativeEnergyOffset(bool RelativeEnergyOffsetMatrix)
{
    if(RelativeEnergyOffsetMatrix)
    {
        SystematicE = RelativeEnergyOffsetE;
    }
}

void CovarianceMatrix3 :: SetAbsoluteEnergyOffset(bool AbsoluteEnergyOffsetMatrix)
{
    if(AbsoluteEnergyOffsetMatrix)
    {
        SystematicE = AbsoluteEnergyOffsetE;
    }
}

void CovarianceMatrix3 :: SetAbsoluteEnergyScaleMatrix(bool AbsoluteEnergyScaleMatrix)
{
    if(AbsoluteEnergyScaleMatrix)
    {
        SystematicE = AbsoluteEnergyE;
    }
}

void CovarianceMatrix3 :: SetRelativeEnergyScaleMatrix(bool RelativeEnergyScaleMatrix)
{
    if(RelativeEnergyScaleMatrix)
    {
        SystematicE = RelativeEnergyE;
    }
}

void CovarianceMatrix3 :: SetIAVMatrix(bool IAVMatrixb)
{
    if(IAVMatrixb)
    {
        SystematicE = IAVE;
    }
}

void CovarianceMatrix3 :: SetNLMatrix(bool NLMatrix)
{
    if(NLMatrix)
    {
        SystematicE = NLE;
    }
}

void CovarianceMatrix3 :: SetResolutionMatrix(bool ResolutionMatrix)
{
    if(ResolutionMatrix)
    {
        SystematicE = ResolutionE;
    }
}

void CovarianceMatrix3 :: SetBCWModel(bool bcw)
{
    if(bcw)
    {
        NLModelE=BCWE;
        NLF = new TF1("NLF",this,&CovarianceMatrix3::NLBCWF,InitialVisibleEnergy,FinalVisibleEnergy,7,"CovarianceMatrix3","NLBCWF");
    }
}

void CovarianceMatrix3 :: SetLBNLModel(bool lbnl)
{
    if(lbnl)
    {
        NLModelE=LBNLE;
        NLF = new TF1("NLF",this,&CovarianceMatrix3::NLLBNLF,InitialVisibleEnergy,FinalVisibleEnergy,5,"CovarianceMatrix3","NLLBNLF");
    }
}

void CovarianceMatrix3 :: SetIHEPModel(bool ihep)
{
    if(ihep)
    {
        NLModelE=IHEPE;
    }
}

void CovarianceMatrix3 :: SetUnifiedModel(bool unified)
{
    if(unified)
    {
        NLModelE=UnifiedE;
        NLF = new TF1("NLF",this,&CovarianceMatrix3::NLUnifiedF,InitialVisibleEnergy,FinalVisibleEnergy,m_num_unified_nl_pars+2,"CovarianceMatrix3","NLUnifiedF");
    }
}

TF1* CovarianceMatrix3 :: GetNeutrinoToVisibleFunction(Int_t order)
{
    if(order==0)
    {
        Double_t Correction = Mn-Mp-Me;
        VisibleF = new TF1("VisibleF",this,&CovarianceMatrix3::VisibleEnergy0F,InitialEnergy,FinalVisibleEnergy,1,"CovarianceMatrix3","VisibleEnergy0F");
        VisibleF->SetParameter(0,Correction);
        std::cout<<"Zeroth order"<<std::endl;
    }
    if(order==1)
    {
        VisibleF = new TF1("VisibleF",this,&CovarianceMatrix3::VisibleEnergy1F,InitialEnergy,FinalVisibleEnergy,1,"CovarianceMatrix3","VisibleEnergy1F");
        VisibleF->SetParameter(0,0);//Maybe I can improve this using a MC simulation of the angular distribution calculated in http://authors.library.caltech.edu/2796/1/VOGprd99.pdf
        std::cout<<"First order"<<std::endl;
    }
    return VisibleF;
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

// First order true to visible function.
//(x-(Mn-Mp))*(1 - x/Mn*(1.0 - [0]*sqrt(1 - Me*Me/(x-(Mn-Mp))/(x-(Mn-Mp)))))- ((Mn-Mp)*(Mn-Mp) - Me*Me)/(2*Mn)-Me]

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Formula is different from previous version (Mp instead of Mn in 1st order, check which is correct)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t CovarianceMatrix3 :: VisibleEnergy1F(Double_t* x, Double_t* par)
{
    Double_t Enu = x[0];

    Double_t Delta = Mn-Mp;
    
    Double_t Ee0 = Enu - Delta; 
    
    Double_t gamma0 = Ee0/Me;
    Double_t v0 = sqrt(1 - 1/gamma0/gamma0);
    
    Double_t costheta = -0.034*v0+2.4*Enu/Mp;
    
    Double_t y2 = (Delta*Delta - Me*Me)/2.;
    
    Double_t Ee1 = Ee0 * (1 - Enu/Mp*(1.0 - v0*costheta)) - y2/Mp; 
    
    return Ee1 + Me;
}

Double_t CovarianceMatrix3 :: VisibleEnergy0F(Double_t* x, Double_t* par)
{
    return x[0]-par[0];
}

Double_t CovarianceMatrix3 :: ResolutionF(Double_t* x, Double_t* par)
{
    Double_t e_orig = x[0];
    Double_t e_sigma = 1.0;
    
    if (e_orig > 0)//To avoid NaNs.
    {
        e_sigma = TMath::Sqrt(par[0]*par[0] + par[1]*par[1]/e_orig + par[2]*par[2]/e_orig/e_orig);
    }
    
    return e_sigma;
}

Double_t CovarianceMatrix3 :: NLBCWF(Double_t* x, Double_t* par)
{
    //Input error
    Double_t e_positron_true = x[0];
    Double_t escale_par = par[5]; // Add flat energy scale parameter
    Double_t escale_offset = par[6]; // Add fixed energy offset
    Double_t scinti_nl_fac = 1;

    if (e_positron_true < m_bcw_positron_nl_e[0])
    {
        scinti_nl_fac = m_bcw_positron_nl_fac[0];
    }
    else if (e_positron_true > m_bcw_positron_nl_e[n_bcw_positron_nl-1])
    {
        scinti_nl_fac = m_bcw_positron_nl_fac[n_bcw_positron_nl-1];
    }
    else
    {
        for (Int_t i = 0; i < n_bcw_positron_nl-1; i++)
        {
            if (e_positron_true >= m_bcw_positron_nl_e[i] && e_positron_true < m_bcw_positron_nl_e[i+1])
            {
                scinti_nl_fac = ((m_bcw_positron_nl_e[i+1] - e_positron_true)*m_bcw_positron_nl_fac[i] + (e_positron_true - m_bcw_positron_nl_e[i])*m_bcw_positron_nl_fac[i+1]) / (m_bcw_positron_nl_e[i+1] - m_bcw_positron_nl_e[i]);
                break;
            }
        }
    }
    Double_t visibleE = scinti_nl_fac * e_positron_true;
    //  double electronicsCorrection = exp(par[0] + par[1] * visibleE) + exp(par[2] + par[3] * visibleE);
    
    Double_t err_offset = 0;
    Double_t err_shift = 0;
    
    Double_t elec_err_min_x = 3.45; // MeV
    Double_t par3_up = g_bcw_elec_nl_error[0]->Eval(elec_err_min_x) - 1;
    Double_t par3_down = g_bcw_elec_nl_error[1]->Eval(elec_err_min_x) - 1;
    
    if (par[3] > 0)
    {
        err_offset = par[3]*par3_up;
    }
    else
    {
        err_offset = par[3]*par3_down;
    }
    
    Double_t par4_up = 0;
    Double_t par4_down = 0;
    
    if (visibleE > elec_err_min_x)
    {
        par4_up = g_bcw_elec_nl_error[0]->Eval(visibleE)- g_bcw_elec_nl_error[0]->Eval(elec_err_min_x);
        par4_down = g_bcw_elec_nl_error[1]->Eval(visibleE)- g_bcw_elec_nl_error[1]->Eval(elec_err_min_x);
    }
    else
    {
        par4_up = g_bcw_elec_nl_error[1]->Eval(visibleE)- g_bcw_elec_nl_error[1]->Eval(elec_err_min_x);
        par4_down = g_bcw_elec_nl_error[0]->Eval(visibleE)- g_bcw_elec_nl_error[0]->Eval(elec_err_min_x);
    }
    if (par[4] > 0) err_shift = par[4]*par4_up;
    else err_shift = par[4]*par4_down;
    
    Double_t electronicsCorrection = exp(par[0] + par[1] * visibleE) + par[2] + err_offset + err_shift;
    
    Double_t final_energy =  visibleE * electronicsCorrection * escale_par + escale_offset;
    
    return final_energy;
}

// LBNL non-linearityr function, based on beta-gamma spectra
Double_t CovarianceMatrix3 :: NLLBNLF(Double_t * x, Double_t * par)
{
    //par[0]: flat energy scale shift for electron
    //par[1]: size of energy scale shift propotional to exp(-1.5*eVis) for electron
    //par[2]: size of energy scale shift due to Ge68 calibration point
    
    Double_t e_positron_true = x[0];
    Double_t escale_par = par[3]; // Add flat energy scale parameter
    Double_t escale_offset = par[4]; // Add fixed energy offset
    
    Double_t scinti_nl_fac = 1;
    Double_t err[3];
    
    if (e_positron_true < m_lbnl_positron_nl_e[0])
    {
        scinti_nl_fac = m_lbnl_positron_nl_fac[0];
        for (Int_t ierr = 0; ierr < 3; ierr++)
        {
            err[ierr] = m_lbnl_positron_nl_err[ierr][0];
        }
    }
    else if (e_positron_true > m_lbnl_positron_nl_e[299])
    {
        scinti_nl_fac = m_lbnl_positron_nl_fac[299];
        for (Int_t ierr = 0; ierr < 3; ierr++)
        {
            err[ierr] = m_lbnl_positron_nl_err[ierr][299];
        }
    }
    else
    {
        for (Int_t i = 0; i < 299; i++)
        {
            if (e_positron_true >= m_lbnl_positron_nl_e[i] && e_positron_true < m_lbnl_positron_nl_e[i+1])
            {
                scinti_nl_fac = ((m_lbnl_positron_nl_e[i+1] - e_positron_true)*m_lbnl_positron_nl_fac[i] + (e_positron_true - m_lbnl_positron_nl_e[i])*m_lbnl_positron_nl_fac[i+1]) / (m_lbnl_positron_nl_e[i+1] - m_lbnl_positron_nl_e[i]);
                for (Int_t ierr = 0; ierr < 3; ierr++)
                {
                    err[ierr] = ((m_lbnl_positron_nl_e[i+1] - e_positron_true)*m_lbnl_positron_nl_err[ierr][i] + (e_positron_true - m_lbnl_positron_nl_e[i])*m_lbnl_positron_nl_err[ierr][i+1]) / (m_lbnl_positron_nl_e[i+1] - m_lbnl_positron_nl_e[i]);
                }
                break;
            }
        }
    }
    double random_nl_fac = scinti_nl_fac;

    for (Int_t ierr = 0; ierr < 3; ierr++)
    {
        random_nl_fac += par[ierr]*err[ierr];
    }
    
    double visibleE = random_nl_fac * e_positron_true;
    //  std::cout << "NL " << e_positron_true << " " << visibleE/e_positron_true << std::endl;
    
    return visibleE  * escale_par + escale_offset;
}

Double_t CovarianceMatrix3 :: NLUnifiedF(Double_t * x, Double_t * par)
{
    Double_t e_positron_true = x[0];
    Double_t escale_par = par[m_num_unified_nl_pars]; // Add flat energy scale parameter
    Double_t escale_offset = par[m_num_unified_nl_pars+1]; // Add fixed energy offset
    
    Double_t scinti_nl_fac = 1;
    Double_t err[10];
    
    if (e_positron_true < m_unified_positron_nl_e[0])
    {
        scinti_nl_fac = m_unified_positron_nl_fac[0];
        
        for (Int_t ierr = 0; ierr < m_num_unified_nl_pars; ierr++)
        {
            err[ierr] = m_unified_positron_nl_err[ierr][0];
        }
    }
    else if (e_positron_true > m_unified_positron_nl_e[n_unified_nl_points-1])
    {
        scinti_nl_fac = m_unified_positron_nl_fac[n_unified_nl_points-1];
        
        for (Int_t ierr = 0; ierr < m_num_unified_nl_pars; ierr++)
        {
            err[ierr] = m_unified_positron_nl_err[ierr][n_unified_nl_points-1];
        }
    }
    else
    {
        for (Int_t i = 0; i < n_unified_nl_points-1; i++)
        {
            if (e_positron_true >= m_unified_positron_nl_e[i] && e_positron_true < m_unified_positron_nl_e[i+1])
            {
                scinti_nl_fac = ((m_unified_positron_nl_e[i+1] - e_positron_true)*m_unified_positron_nl_fac[i] + (e_positron_true - m_unified_positron_nl_e[i])*m_unified_positron_nl_fac[i+1]) / (m_unified_positron_nl_e[i+1] - m_unified_positron_nl_e[i]);
                for (Int_t ierr = 0; ierr < m_num_unified_nl_pars; ierr++)
                {
                    err[ierr] = ((m_unified_positron_nl_e[i+1] - e_positron_true)*m_unified_positron_nl_err[ierr][i] + (e_positron_true - m_unified_positron_nl_e[i])*m_unified_positron_nl_err[ierr][i+1]) / (m_unified_positron_nl_e[i+1] - m_unified_positron_nl_e[i]);
                }
                break;
            }
        }
    }
    Double_t random_nl_fac = scinti_nl_fac;
    
    for (Int_t ierr = 0; ierr < m_num_unified_nl_pars; ierr++)
    {
        random_nl_fac += par[ierr]*err[ierr];
    }
    double visibleE = random_nl_fac * e_positron_true;
    return visibleE  * escale_par + escale_offset;
}