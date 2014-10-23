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
#include "Prediction.h"
#include "ReactorSpectrumMultiple.h"
#include "NominalData.h"
#include <math.h>
#include <TMatrixD.h>
#include <TDecompChol.h>

const Int_t NSamples = 1;//Number of samples used in the generation of covariance matrices
const Int_t MaxPeriods = 1;
const Int_t MaxDetectors = 8;
const Int_t MaxNearDetectors =4;
const Int_t MaxFarDetectors =4;
const Int_t Halls=3;
const Int_t MaxNbins=51;
const Int_t MatrixBins=240;//Same than IAV matrix for Gd. When the IAV H matrix is produced this can be selected through NominalData.h

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

class CovarianceMatrix3
{
private:
    TRandom3* rand;
    Prediction Pred;

    enum Systematic{ReactorE, EnergyE, IAVE, NLE, ResolutionE,StatisticalE};
    Systematic SystematicE;
    enum Background{VaryAccidentalE, VaryLiHeE, VaryFastNeutronsE, VaryAmCE, DistortLiHeE, DistortFastNeutronsE, DistortAmCE };
    Background BackgroundE;
    enum NLModel{BCWE, LBNLE, IHEPE, UnifiedE};//Add here new NL models
    NLModel NLModelE;
    
    //AD configuration parameters:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    Int_t hall;
    
    //Binning parameters:
    bool LinearBinning;
    Int_t n_evis_bins;
    Int_t n_etrue_bins;
    Int_t n_cov_matrix_bins;
    Double_t Difference[MaxNbins+1];

    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    Int_t Nweeks;
    Double_t BinWidth;
    Int_t TotalBins;
    char* OutputFileName;
    
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.
    Double_t cov_matrix_bins[9*(MaxNbins)+1];
    //Response Matrix:
    Double_t PositronEnergy;
    Int_t PositronEnergyIndex;
    Double_t Norma[MaxNbins+1];
    TH1F* PositronTrueSpectrumH[MaxDetectors][MatrixBins];
    TH1F* PositronIAVSpectrumH[MaxDetectors][MatrixBins];
    TH1F* PositronNLSpectrumH[MaxDetectors][MatrixBins];
    TH1F* PositronVisibleSpectrumH[MaxDetectors][MatrixBins];
    
    //IAV:
    Double_t IAVError; // relative uncertainty of the IAV thickness
    Double_t IAVError[MaxDetectors]; // relative uncertainty of the IAV thickness
    Double_t IAVMatrix[MatrixBins][MatrixBins];
    
    //Non linearity:
    // Detector response non-linearlity function
    TF1 * nl_func;
    
    //Interpolation vectors
    vector<Int_t> sign;
    vector<Int_t> EnergyIdx;
    vector<Double_t> Energy;
    vector<Double_t> dEtrue;
    vector<Double_t> binScaling;
    vector<Double_t> dNdE;
    
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
    
    //Resolution
    // Detector resolution function
    TF1 * reso_func;
    
    Double_t resolutionRange; //Range of resolution
    Double_t m_detectorResolution; // Detector energy resolution parameter
    Double_t m_detectorResolution_nominal; // Detector energy resolution parameter
    Double_t m_detectorResolution_error; // Detector energy resolution parameter
    Double_t m_detectorResolution_error_uncorr; // Detector energy resolution parameter
    Double_t m_detectorResolution_bias[MaxDetectors]; // Random biases of the detector resolution.
    
    Double_t m_abs_escale;
    Double_t m_abs_escale_nominal;
    Double_t m_abs_escale_error;
    
    Double_t m_abs_eoffset;
    Double_t m_abs_eoffset_error;
    
    Double_t m_rel_eoffset_error;
    
    //Background errors:
    Double_t HAccidentalError[MaxDetectors];
    Double_t HLiHeError[Halls];
    Double_t HFastNeutronsError[Halls];
    Double_t HAmCError[Halls];
    
    //Distortion functions:
    Double_t DistortAcc;
    Double_t DistortLiHe;
    Double_t DistortFN;
    Double_t DistortAmC;

    Double_t ScaleFactorAccidental[MaxDetectors][MaxPeriods];
    Double_t ScaleFactorLiHe[MaxDetectors][MaxPeriods];
    Double_t ScaleFactorFastNeutrons[MaxDetectors][MaxPeriods];
    Double_t ScaleFactorAmC[MaxDetectors][MaxPeriods];

    //Choose matrix
    bool AccidentalMatrix;
    bool LiHeMatrix;
    bool FastNeutronsMatrix;
    bool AmCMatrix;
    
    bool DistortBackgrounds;
    bool VariateRate;
    bool AddBackgrounds;

    char* OutputFileName;
    
    //Histograms
    TH1F* OriginalPredictionH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    TH1F* PredictionTrueH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    TH1F* PredictionVisH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    TH1F* PredictionVisWithBkgdH[MaxNearDetectors][MaxNearDetectors][MaxPeriods];
    TH1F* PredictionLA[MaxPeriods];
    TH1F* PredictionDB[MaxPeriods];
    TH1F* PredictionLAWithBkgd[MaxPeriods];
    TH1F* PredictionDBWithBkgd[MaxPeriods];

    TH1F* OriginalNearHallSpectrumH[MaxNearDetectors][MaxPeriods];
    TH1F* NearHallSpectrumTrueH[MaxNearDetectors][MaxPeriods];
    TH1F* NearHallSpectrumVisH[MaxNearDetectors][MaxPeriods];
    TH1F* NearHallVisWithBkgdH[MaxNearDetectors][MaxPeriods];
    
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
    TH2F* CovMatrixDB2H;
    TH2F* CovMatrixLA2H;
    TH2F* StatMatrixDB2H;
    TH2F* StatMatrixLA2H;
    TH2F* StatCovDB2H;
    TH2F* StatCovLA2H;
    Double_t CovDB[MaxNbins][MaxNbins][MaxPeriods];
    Double_t CovLA[MaxNbins][MaxNbins][MaxPeriods];
    TH2F* CovDB2H;
    TH2F* CovLA2H;

    Double_t Sigma_Near[MaxFarDetectors][MaxNearDetectors][MaxPeriods][MaxNbins];
    Double_t Sigma_Far[MaxFarDetectors][MaxNearDetectors][MaxPeriods][MaxNbins];
    Double_t CovStat[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t NormCovStat[9*MaxNbins][9*MaxNbins][MaxPeriods];
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
    void LoadResponseMatrix();
    void CombinePredictions();
    void ApplyDetectorResponse(TH1F* Prediction, TH2F* Response);

    //Functions to vary background shapes. So far the same ones than LBNL.
    void FluctuateBackgrounds();
    TF1* GetDistortionFunction(Double_t amount);
    TF1* GetFastNeutronsDistortionFunction(Double_t amount);

    //Functions to calculate the energy matrix
    TF1* GetNeutrinoToVisibleFunction();
    void GetEnergyShift(Int_t TrueEnergyIndex);
    void GetIAVShift(Int_t TrueEnergyIndex);
    void GetNLShift(Int_t TrueEnergyIndex);
    void GetResolutionShift(Int_t TrueEnergyIndex);
    void GetStatisticalFluctuation(TH1F* Histo);
    void CreateEnergyMatrix();
    void Interpolation(TF1* func);

    //Functions to generate the respective covariance matrices
    void GenerateCovarianceMatrix();
    TH2F* NormalizeCovariance(Double_t Cov[][MaxNbins][MaxPeriods]);
    void GenerateStatisticalCovarianceMatrix();
    
    Double_t VisibleEnergy0F(Double_t * x, Double_t * par);// Zeroth order true to visible function.
    Double_t VisibleEnergy1F(Double_t * x, Double_t * par);// First order true to visible function.
    
    Double_t nl_func(Double_t * x, Double_t * par);
    Double_t nl_func_bcw(Double_t* x, Double_t* par);
    Double_t nl_func_lbnl(Double_t* x, Double_t* par);
    Double_t nl_func_unified(Double_t* x, Double_t* par);

    Double_t reso_func_bcw(Double_t * x, Double_t * par);
    
public:
    CovarianceMatrix3();
    CovarianceMatrix3(NominalData* Data);

    void CovarianceMatrixMain();

    void SetVaryAccidentalMatrix(bool AccMatrix);
    void SetVaryLiHeMatrix(bool LiMatrix);
    void SetVaryFastNeutronsMatrix(bool FNmatrix);
    void SetVaryAmCMatrix(bool AmCmatrix);
    
    void SetDistortLiHeMatrix(bool LiMatrix);
    void SetDistortFastNeutronsMatrix(bool FNmatrix);
    void SetDistortAmCMatrix(bool AmCmatrix);
    
    void SetIsotopeMatrix(bool IsotopeMatrix);
    void SetReactorPowerMatrix(bool ReactorPowerMatrix);
    void SetIAVMatrix(bool IAVMatrix);
    void SetNLMatrix(bool NLMatrix);
    void SetResolutionMatrix(bool ResolutionMatrix);
    void SetStatisticalMatrix(bool StatisticalMatrix);

    void RandomEnergyScaleMatrix();
    void RandomIAVMatrix();
    void RandomNLMatrix();
    void RandomResolutionMatrix();
    
    void SetBCWModel(bool BCW);
    void SetLBNLModel(bool LBNL);
    void SetIHEPModel(bool IHEP);
    void SetUnifiedModel(bool Unified);

    void LoadNLParameters();

};
CovarianceMatrix3 :: CovarianceMatrix3()
{
    Nom = new NominalData();
    rand = new TRandom3();
    Pred = new Prediction();

    n_evis_bins = Nom->GetNbins();
    n_etrue_bins = Nom->GetNbins();
    InitialEnergy = Nom->GetEmin();
    FinalEnergy = Nom->GetEmax();
    InitialVisibleEnergy = Nom->GetEVisMin();
    FinalVisibleEnergy =  Nom->GetEVisMax();
    TotalBins = MatrixBins;
    BinWidth=(FinalVisibleEnergy-InitialVisibleEnergy)/TotalBins;
    
    LinearBinning = Nom->GetBinning();
    
    //Linear binning
    if(LinearBinning)
    {
        for (Int_t i = 0; i <= n_evis_bins; i++)
        {
            evis_bins[i] = 0.2 * i + InitialEnergy;
            enu_bins[i] = 0.2 * i + InitialEnergy;
        }
        for (Int_t j = 0; j < 9; j++)
        {
            Double_t Sum=0;
            for (Int_t i = 0; i <= n_evis_bins; i++)
            {
                cov_matrix_bins[i+j*n_evis_bins]=j*evis_bins[n_evis_bins]+Sum;
                Sum+=(evis_bins[i+1]-evis_bins[i]);
            }
        }
        n_cov_matrix_bins=9*(n_evis_bins);
    }
    //Non-linear binning
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
        
        n_cov_matrix_bins=9*(n_evis_bins);
        
        for (Int_t i = 0; i <= n_evis_bins; i++)
        {
            Difference[i]+=(evis_bins[i+1]-evis_bins[i]);
            //            cout << Difference[i] << " IS SUM" << i << "AT I " << endl;
        }
//
//        for (Int_t j = 0; j < 9; j++)
//        {
//            cov_matrix_bins[0]=0;
//            Double_t Sum=0;
//            
//            for (Int_t i = 0; i < n_evis_bins; i++)
//            {
//                Sum=Sum+Difference[i];
//                cov_matrix_bins[1+i+j*n_evis_bins]=j*(evis_bins[n_evis_bins]-0.7)+Sum;
//                cout << cov_matrix_bins[1+i+j*n_evis_bins] << "cov_matrix_bins" << endl;
//
//            }
//            //            cov_matrix_bins[n_cov_matrix_bins]= cov_matrix_bins[n_cov_matrix_bins-1]+ Difference[n_evis_bins-1];
//            Difference[i] =0;
//        }
        for (Int_t j = 0; j < 9; j++)
        {
            Double_t Sum =0;
            
            for (Int_t i = 0; i < n_evis_bins; i++)
            {
                cov_matrix_bins[i+j*n_evis_bins]=j*evis_bins[n_evis_bins]+Sum;
                
                Sum+=(evis_bins[i+1]-evis_bins[i]);
            }
        }
        n_cov_matrix_bins=9*(n_evis_bins);
        
    }
    
    Nweeks = Nom->GetWeeks();
    
    NADs = Nom->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    if(NADs == 8) //ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    for(Int_t i=0;i<NADs;i++)
    {
        HAccidentalError[i]=Nom->GetHAccidentalError(i);
    }
    
    for(Int_t i=0;i<3;i++)
    {
        HLiHeError[i]=Nom->GetHLiHeError(i);
    }
    
    for(Int_t i=0;i<3;i++)
    {
        HFastNeutronsError[i]=Nom->GetHFNError(i);
    }
    
    for(Int_t i=0;i<3;i++)
    {
        HAmCError[i]=Nom->GetHAmCError(i);
    }
    
    IAVError=Nom->GetIAVError();
    
    m_abs_escale =1;
    m_abs_escale_nominal =1;
    m_abs_escale_error = 0.01;
    
    m_abs_eoffset = 0.0;
    m_abs_eoffset_error = 0.08;  //(MeV)
    
    m_rel_eoffset_error = 0.013; //(MeV)
    
    for(int idet=0; idet<3; idet++)
    {
        m_rel_escale[idet] = 1.0;
        m_rel_escale_error[idet] = 0.0035; // 0.35%
        m_rel_escale_nominal[idet] = m_rel_escale[idet];
        m_rel_eoffset[idet] = 0.0;
    }
       
    SystematicE=-1;//A number different to any of the different systematics included in the model
    DistortAcc = 0;
    DistortLiHe= 0;
    DistortFN  = 0;
    DistortAmC = 0;
}

CovarianceMatrix3 :: CovarianceMatrix3(NominalData* Data)
{
    rand = new TRandom3();
    Pred = new Prediction(Data);
    
    n_evis_bins = Data->GetNbins();
    n_etrue_bins = Data->GetNbins();
    InitialEnergy = Data->GetEmin();
    FinalEnergy = Data->GetEmax();
    InitialVisibleEnergy =Data->GetEVisMin();
    FinalVisibleEnergy = Data->GetEVisMax();
    LinearBinning = Data->GetBinning();
    TotalBins = MatrixBins;
    BinWidth=(FinalVisibleEnergy-InitialVisibleEnergy)/TotalBins;
    
    //Linear binning
    if(LinearBinning)
    {
        for (Int_t i = 0; i <= n_evis_bins; i++)
        {
            evis_bins[i] = 0.2 * i + InitialEnergy;
            enu_bins[i] = 0.2 * i + InitialEnergy;
        }
        for (Int_t j = 0; j < 9; j++)
        {
            Double_t Sum=0;
            for (Int_t i = 0; i <= n_evis_bins; i++)
            {
                cov_matrix_bins[i+j*n_evis_bins]=j*evis_bins[n_evis_bins]+Sum;                
                Sum+=(evis_bins[i+1]-evis_bins[i]);
            }
        }
        n_cov_matrix_bins=9*(n_evis_bins);
    }
    //Non-linear binning
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
        
        for (Int_t i = 0; i < n_evis_bins; i++)
        {
            Difference[i]+=(evis_bins[i+1]-evis_bins[i]);
            //            cout << Difference[i] << " IS SUM" << i << "AT I " << endl;
        }
        //
        //        for (Int_t j = 0; j < 9; j++)
        //        {
        //            cov_matrix_bins[0]=0;
        //            Double_t Sum=0;
        //
        //            for (Int_t i = 0; i < n_evis_bins; i++)
        //            {
        //                Sum=Sum+Difference[i];
        //                cov_matrix_bins[1+i+j*n_evis_bins]=j*(evis_bins[n_evis_bins]-0.7)+Sum;
        //                cout << cov_matrix_bins[1+i+j*n_evis_bins] << "cov_matrix_bins" << endl;
        //
        //            }
        //            //            cov_matrix_bins[n_cov_matrix_bins]= cov_matrix_bins[n_cov_matrix_bins-1]+ Difference[n_evis_bins-1];
        //            Difference[i] =0;
        //        }
        for (Int_t j = 0; j < 9; j++)
        {
            Double_t Sum =0;
            
            for (Int_t i = 0; i <= n_evis_bins; i++)
            {
                cov_matrix_bins[i+j*n_evis_bins]=j*evis_bins[n_evis_bins]+Sum;
                
                Sum+=(evis_bins[i+1]-evis_bins[i]);
            }
        }
        n_cov_matrix_bins=9*(n_evis_bins);
    }
    
    Nweeks = Data->GetWeeks();

    NADs = Data->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;

    if(NADs == 8) //ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    for(Int_t i=0;i<NADs;i++)
    {
        HAccidentalError[i]=Data->GetHAccidentalError(i);
    }
    
    for(Int_t i=0;i<3;i++)
    {
        HLiHeError[i]=Data->GetHLiHeError(i);
    }
    
    for(Int_t i=0;i<3;i++)
    {
        HFastNeutronsError[i]=Data->GetHFNError(i);
    }

    for(Int_t i=0;i<3;i++)
    {
        HAmCError[i]=Data->GetHAmCError(i);
    }
    
    IAVError=Data->GetIAVError();
  
    m_abs_escale =1;
    m_abs_escale_nominal =1;
    m_abs_escale_error = 0.01;
    
    m_abs_eoffset = 0.0;
    m_abs_eoffset_error = 0.08;  //(MeV)
    
    m_rel_eoffset_error = 0.013; //  (MeV)
    
    for(int idet=0; idet<3; idet++)
    {
        m_rel_escale[idet] = 1.0;
        m_rel_escale_error[idet] = 0.0035; // 0.35%
        m_rel_escale_nominal[idet] = m_rel_escale[idet];
        m_rel_eoffset[idet] = 0.0;
    }
    
    SystematicE=-1;//A number different to any of the different systematics included in the model
    DistortAcc = 0;
    DistortLiHe= 0;
    DistortFN  = 0;
    DistortAmC = 0;
}

void CovarianceMatrix3 :: CovarianceMatrixMain()
{
    LoadBackgrounds();
    LoadPredictions();
    LoadNearHall();
    LoadIavCorrection();
    LoadResponseMatrix();

    cout << "Randomize Background #" << BackgroundE << endl;
    cout << "Randomize Systematic #" << SystematicE << endl;

//    CovMatrix2H = new TH2F("Covariance Matrix","Covariance Matrix",n_evis_bins*9,0,n_evis_bins*9,n_evis_bins*9,0,n_evis_bins*9);
//    StatMatrix2H = new TH2F("Statistical Cov Matrix","Statistical Cov Matrix",n_evis_bins*9,0,n_evis_bins*9,n_evis_bins*9,0,n_evis_bins*9);
    CovMatrix2H = new TH2F("Covariance Matrix","Covariance Matrix",n_cov_matrix_bins,cov_matrix_bins,n_cov_matrix_bins,cov_matrix_bins);
    StatMatrix2H = new TH2F("Statistical Cov Matrix","Statistical Cov Matrix",n_cov_matrix_bins,cov_matrix_bins,n_cov_matrix_bins,cov_matrix_bins);

    CovarianceMatrix3* NLObject = new CovarianceMatrix3();//Finally the standard constructor is useful
    switch (NLModelE)
    {
        case 0://BCW NL Model
            NLObject->SetBCWModel(1);
            break;
        case 1://LBNL NL Model
            NLObject->SetLBNLModel(1);
            break;
        case 2://IHEP NL Model
            NLObject->SetIHEPModel(1);
            break;
        case 3://Unified NL Model
            NLObject->SetUnifiedModel(1);
            break;
        default:
    }
    NLObject->LoadNLParameters();
    for (Int_t samples = 0; samples<NSamples; samples++)
    {
        switch (SystematicE)
        {
            case 0://Vary Reactor
                Pred.MakePrediction();
                break;
            case 1://Vary Energy Scale
                RandomEnergyScaleMatrix();
                break;
            case 2://Vary IAV
                RandomIAVMatrix();
                break;
            case 3://Vary NL
                NLObject->LoadNLParameters();//Reset
                RandomNLMatrix();
                break;
            case 4://Vary Resolution
                RandomResolutionMatrix();
                break;
            case 5://Statistical fluctuation
//                GenerateStatisticalCovarianceMatrix();
                break;
            default://Add nominal systematics
        }
        
        switch (NLModelE)
        {
            case 0://BCW NL Model
                cout << "Using BCW NL Model"<< endl;
                nl_func = new TF1("nl_func",NLObject,&CovarianceMatrix3::nl_func_bcw,InitialVisibleEnergy,FinalVisibleEnergy,7,"CovarianceMatrix3","nl_func_bcw");
                for (Int_t idet = 0; idet < ADsEH3; idet++)
                {
                    //NL function set up
                    nl_func->SetParameters(m_bcw_elec_nl_par[0],m_bcw_elec_nl_par[1],m_bcw_elec_nl_par[2], m_bcw_elec_nl_par[3], m_bcw_elec_nl_par[4], m_abs_escale * m_rel_escale[idet], m_abs_eoffset+m_rel_eoffset[idet]);
                }
                break;
            case 1://LBNL NL Model
                cout << "Using LBNL NL Model"<< endl;
                nl_func = new TF1("nl_func",NLObject,&CovarianceMatrix3::nl_func_lbnl,InitialVisibleEnergy,FinalVisibleEnergy,5,"CovarianceMatrix3","nl_func_lbnl");
                for (Int_t idet = 0; idet < ADsEH3; idet++)
                {
                    nl_func->SetParameters(m_lbnl_nl_par[0],m_lbnl_nl_par[1],m_lbnl_nl_par[2],m_abs_escale * m_rel_escale[idet],m_abs_eoffset + m_rel_eoffset[idet]);
                }
                
                break;
            case 2:
                break;
            case 3://Unified NL Model
                cout << "Using Unified NL Model"<< endl;
                nl_func = new TF1("nl_func",NLObject,&CovarianceMatrix3::nl_func_unified,InitialVisibleEnergy,FinalVisibleEnergy,m_num_unified_nl_pars+2,"CovarianceMatrix3","nl_func_unified");
  
                for (Int_t i = 0; i < m_num_unified_nl_pars; i++)
                {
                    nl_func->SetParameter(i,m_unified_nl_par[i]);
                }
                
                for (Int_t idet = 0; idet < ADsEH3; idet++)
                {
                    nl_func->SetParameter(m_num_unified_nl_pars, m_abs_escale * m_rel_escale[idet]);
                    nl_func->SetParameter(m_num_unified_nl_pars+1, m_abs_eoffset + m_rel_eoffset[idet]);
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
//                cout << "Covariance matrix for UNIFIED non-liner parameters:" << endl;
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

            default:
        }

        //Resolution function set up
        CovarianceMatrix3* ResoObject = new CovarianceMatrix3();//Finally the standard constructor is useful
        
        reso_func = new TF1("reso_func",ResoObject,&CovarianceMatrix3::reso_func_bcw,0,20,3,"CovarianceMatrix3","reso_func_bcw");
        reso_func->SetParameters(0.022,0.077,0.018); // based on Bryce's TN
        //BCW values http://dayabay.ihep.ac.cn/DocDB/0087/008768/013/6AdAnalysis-BCW.pdf are different! 0.13226034931, 0.32604048828, 0.26196488314
        
        //Resolution
        resolutionRange = 8; // Why 8σ? Seems chosen trivially but in my opinion it's the way to limit the range of the convolution, for a Normal distribution this covers up to 99.99999... of the area
        for(Int_t i=0;i<NADs;i++)
        {
            m_detectorResolution_bias[i] = 0;//this will be used to distort it
        }
        
        //Energy shift function set up
        VisibleF=GetNeutrinoToVisibleFunction(0);//0 for 0th order, 1 for 1st order
        Interpolation(nl_func);
        for(Int_t TrueEnergyIndex=0; TrueEnergyIndex<TotalBins; TrueEnergyIndex++)
        {
            GetEnergyShift(TrueEnergyIndex);
            GetIAVShift(TrueEnergyIndex);
            GetNLShift(TrueEnergyIndex);
            GetResolutionShift(TrueEnergyIndex);
        }
        
        CreateEnergyMatrix();

        if(SystematicE!=StatisticalE)
        {
            CovMatrixDB2H=(TH2F*)CovDB2H->Clone();
            CovMatrixLA2H->(TH2F*)CovLA2H->Clone();

            CovMatrixDB2H->Reset();
            CovMatrixDB2H->Reset();

            GenerateCovarianceMatrix();
            CovMatrixDB2H->Add(CovDB2H);
            CovMatrixLA2H->Add(CovLA2H);
        }
        else
        {
            StatMatrixDB2H=(TH2F*)StatCovDB2H->Clone();
            StatMatrixLA2H->(TH2F*)StatCovLA2H->Clone();
            
             StatMatrixDB2H->Reset();
            StatMatrixLA2H->Reset();
            
            StatMatrixLA2H->Add(StatCovDB2H);
            StatMatrixLA2H->Add(StatCovLA2H);
        }
        //Clean up
        delete nl_func;
        delete reso_func;
        delete VisibleF;
        ResoObject->~CovarianceMatrix3();
        NLObject->~CovarianceMatrix3();
        
        for (Int_t idet = 0; idet<ADsEH3; idet++)
        {
            delete EnergyMatrixH[idet];
            for (Int_t TrueEnergyIndex= 0; TrueEnergyIndex<TotalBins; TrueEnergyIndex++)
            {
                dEtrue.pop_back();
                Energy.pop_back();
                EnergyIdx.pop_back();
                sign.pop_back();
                delete PositronTrueSpectrumH[idet][TrueEnergyIndex];
                delete PositronIAVSpectrumH[idet][TrueEnergyIndex];
                delete PositronNLSpectrumH[idet][TrueEnergyIndex];
                delete PositronVisibleSpectrumH[idet][TrueEnergyIndex];
            }
        }
    }
    //Save results
    if(SystematicE!=StatisticalE)
    {
        CovMatrix2H->Scale(1./(NSamples));
        // CovMatrix2H->Draw("colz");
    }
    else
    {
        StatMatrix2H->Scale(1./(NSamples));
    }
    
    SaveTotalSpectrum();
    SaveCovarianceMatrix();
}

void CovarianceMatrix3 :: ApplyDetectorResponse(TH1F* Prediction, TH2F* Response)
{
    //One should be multiplied by the nominal Energy Matrix, and the other one by the one fluctuated.
    //This matrix will be the nominal one:
    for(Int_t j=1; j<=n_evis_bins; j++)
    {
        Int_t Sum = 0;
        
        for(Int_t i=1;i<=n_etrue_bins;i++)
        {
            Sum = Sum+(Response->GetBinContent(i,j)*Prediction->GetBinContent(i));
        }
        if (Sum<0)
        {
            cout << "SUM NEGATIVE!!!!" << endl;//It should never be negative, just in case.
            Sum=0;//Fix temporary negative sum in the last bins
        }
        Prediction->SetBinContent(j,Sum);
    }
}

void CovarianceMatrix3 :: GenerateCovarianceMatrix()
{
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for(Int_t AD=0; AD<NADs; AD++)
        {
            if(AddBackgrounds)
            {
                BackgroundSpectrumH[AD][week]=(TH1F*)AccidentalsH[0][0]->Clone();
                BackgroundSpectrumH[AD][week]->Reset();
                RandomBackgroundSpectrumH[AD][week]=(TH1F*)BackgroundSpectrumH[AD][week]->Clone();

                //Add nominal backgrounds, without being distorted
                
                BackgroundSpectrumH[AD][week]->Add(AccidentalsH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(LiHeH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(FastNeutronsH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(AmCH[AD][week]);
                
                FluctuateBackgrounds();//Vary backgrounds
                
                RandomBackgroundSpectrumH[AD][week]->Add(RandomAccidentalsH[AD][week]);
                RandomBackgroundSpectrumH[AD][week]->Add(RandomLiHeH[AD][week]);
                RandomBackgroundSpectrumH[AD][week]->Add(RandomFastNeutronsH[AD][week]);
                RandomBackgroundSpectrumH[AD][week]->Add(RandomAmCH[AD][week]);
            }
        }
        
        for (Int_t near=0; near<(ADsEH1+ADsEH2); near++)
        {
            //Add detector effects //Multiply Matrix and NearHallSpectrum

            NearHallSpectrumVisH[near][week] = new TH1F(Form("Corrected Random Near AD%d",near),Form("Corrected Random Near AD%d",near),n_evis_bins,evis_bins);
            NearHallVisWithBkgdH[near][week] = (TH1F*)NearHallSpectrumVisH[near][week]->Clone();

            //One should be multiplied by the nominal Energy Matrix, and the other one by the one fluctuated.
            //Nominal in this one:
            ApplyDetectorResponse(NearHallSpectrumVisH[near][week], NominalResponseMatrixH[0]);
            //Fluctuations in this one:
            ApplyDetectorResponse(NearHallVisWithBkgdH[near][week], EnergyMatrixH[0]);
            
            for (Int_t far=0; far<(ADsEH3); far++)
            {
                //Add detector effects //Multiply Matrix and Prediction (Far Spectrum)
                
                PredictionVisH[far][near][week] = new TH1F(Form("Visible Far Spectrum at AD%d from AD%d",far+1,near+1),Form("Nominal Visible Far Spectrum at AD%d from AD%d",far+1,near+1),n_evis_bins,evis_bins);
                PredictionVisWithBkgdH[far][near][week] = (TH1F*)PredictionVisH[far][near][week]->Clone(Form("Visible Far Spectrum at AD%d from AD%d with variations due to Systematic%d or Background%d",far+1,near+1,SystematicE,BackgroundE));
                
                //Nominal in this one:
                ApplyDetectorResponse(PredictionVisH[far][near][week],NominalResponseMatrixH[far]);
                //Fluctuations in this one:
                ApplyDetectorResponse(PredictionVisWithBkgdH[far][near][week],EnergyMatrixH[far]);
               
                if(AddBackgrounds)
                {
                    //(ADD BACKGROUNDS AFTER RESPONSE MATRIX SINCE THEY ALREADY INCLUDE MC EFFECTS OF THE DETECTOR)
                    //Add nominal backgrounds 
                    PredictionVisH[far][near][week]->Add(BackgroundSpectrumH[far+ADsEH1+ADsEH2][week]);
                    //Add random backgrounds 
                    PredictionVisWithBkgdH[far][near][week]->Add(RandomBackgroundSpectrumH[far+ADsEH1+ADsEH2][week]);
                }
            }
        }
        
        CombinePredictions();

        for (Int_t i = 0; i<n_evis_bins; i++)
        {//columns
            for (Int_t j = 0; j<n_evis_bins; j++)
            {//rows
                CovDB[i][j][week]=(PredictionDBWithBkgd[week]->GetBinContent(i+1)-PredictionDB[week]->GetBinContent(i+1))*(PredictionDBWithBkgd[week]->GetBinContent(j+1)-PredictionDB[week]->GetBinContent(j+1));
                CovLA[i][j][week]=(PredictionLAWithBkgd[week]->GetBinContent(i+1)-PredictionLA[week]->GetBinContent(i+1))*(PredictionLAWithBkgd[week]->GetBinContent(j+1)-PredictionLA[week]->GetBinContent(j+1));
            }
        }
        
        TH2F* CovDB2H = NormalizeCovariance(CovDB);
        TH2F* CovLA2H = NormalizeCovariance(CovLA);

        //Save spectrum without backgrounds
        
        TFile* CombinedPredictionsF = TFile::Open("./RootOutputs/CombinedNominalSpectrum.root","recreate");
        for (Int_t week = 0; week<Nweeks; week++)
        {
            PredictionDB[week]->Write();
            PredictionLA[week]->Write();
        }
        CombinedPredictionsF->Close();
        
        //Save spectrum with backgrounds
        
        TFile* CombinedPredictionsF = TFile::Open("./RootOutputs/CombinedVariedSpectrum.root","recreate");
        for (Int_t week = 0; week<Nweeks; week++)
        {
            PredictionDBWithBkgd[week]->Write();
            PredictionLAWithBkgd[week]->Write();
        }
        CombinedPredictionsF->Close();
        
        //Save spectrum without backgrounds
        
        TFile* PredictionsF = TFile::Open("./RootOutputs/DetectorEvisSpectrum.root","recreate");
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
            {
                for (Int_t far =0; far<ADsEH3; far++)
                {
                    PredictionVisH[far][near][week]->Write();
                }
                NearHallSpectrumVisH[near][week]->Write();
            }
        }
        PredictionsF->Close();
        
        //Save spectrum with backgrounds
        TFile* CorrectedPredictionsF = TFile::Open("./RootOutputs/DetectorCorrectedEvisSpectrum.root","recreate");
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
            {
                for (Int_t far =0; far<ADsEH3; far++)
                {
                    PredictionVisWithBkgdH[far][near][week]->Write();
                }
                NearHallVisWithBkgdH[near][week]->Write();
            }
        }
        CorrectedPredictionsF->Close();
    }
}

////////////////////////////////////////////////////////////////////////////////////
//// I need to correct backgrounds for efficiencies, and maybe events too (Check)
////////////////////////////////////////////////////////////////////////////////////
void CovarianceMatrix3 :: GenerateStatisticalCovarianceMatrix()//FIX THIS TO BE ABLE TO RUN IT ONLY ONCE, SAME THAN COV MATRIX AND APPLY STATISTICAL FLUCTUATION TO THE WHOLE SPECTRA+BACKGROUND
{
//    for (Int_t week = 0; week<Nweeks; week++)
//    {
//        FluctuateBackgrounds();//Just to se
//        GetStatisticalFluctuation(PredictionVisWithBkgdH[week]);
//
//        //Add nominal backgrounds
//        for(Int_t AD=0; AD<NADs; AD++)
//        {
//            BackgroundSpectrumH[AD][week]=(TH1F*)AccidentalsH[0][0]->Clone();
//            BackgroundSpectrumH[AD][week]->Reset();
//
//            BackgroundSpectrumH[AD][week]->Add(AccidentalsH[AD][week]);
//            BackgroundSpectrumH[AD][week]->Add(LiHeH[AD][week]);
//            BackgroundSpectrumH[AD][week]->Add(FastNeutronsH[AD][week]);
//            BackgroundSpectrumH[AD][week]->Add(AmCH[AD][week]);
//
//        }
//        for (Int_t far=0; far<ADsEH3; far++)
//        {
//            for (Int_t near=0; near<(ADsEH1+ADsEH2); near++)
//            {
//                for (Int_t pts = 0; pts < n_evis_bins; pts++)
//                {
//                    Sigma_Far[far][near][week][pts]=sqrt(PredictionVisH[far][near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[ADsEH1+ADsEH2+far][week]->GetBinContent(pts+1));
//                    
//                    if ((NearHallSpectrumVisH[near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[near][week]->GetBinContent(pts+1))!=0)
//                    {
//                        if((NearHallSpectrumVisH[near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[near][week]->GetBinContent(pts+1))<0)
//                        {
//                            Sigma_Near[far][near][week][pts]=(PredictionVisH[far][near][week]->GetBinContent(pts+1)/NearHallSpectrumVisH[near][week]->GetBinContent(pts+1))*sqrt(-1*(NearHallSpectrumVisH[near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[near][week]->GetBinContent(pts+1)));
//                        }
//                        else
//                        {
//                            Sigma_Near[far][near][week][pts]=(PredictionVisH[far][near][week]->GetBinContent(pts+1)/NearHallSpectrumVisH[near][week]->GetBinContent(pts+1))*sqrt(NearHallSpectrumVisH[near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[near][week]->GetBinContent(pts+1));
//                        }
//                    }
//                    else
//                    {
//                        Sigma_Near[far][near][week][pts]=0;
//                    }
//                }
//            }
//        }
//    }
//    
////    StatCov2H = new TH2F("Statistical Covariance Matrix","Statistical Covariance Matrix",n_evis_bins*9,0,n_evis_bins*9,n_evis_bins*9,0,n_evis_bins*9);
//    StatCov2H = new TH2F("Statistical Covariance Matrix","Statistical Covariance Matrix",n_cov_matrix_bins,cov_matrix_bins,n_cov_matrix_bins,cov_matrix_bins);
//
//
//    for (Int_t week = 0; week<Nweeks; week++)
//    {
//        for (Int_t i = 0; i<n_evis_bins; i++)
//        {//columns
//            for (Int_t j = 0; j<n_evis_bins; j++)
//            {//rows
//                //Near component correlated
//                if(neari==nearj && fari!=farj)
//                {
//                    CovStat[x][y][week]=Sigma_Near[fari][neari][week][i]*Sigma_Near[farj][nearj][week][j];
//                }
//                //Far component correlated
//                if(fari==farj && neari!=nearj)
//                {
//                    CovStat[x][y][week]=Sigma_Far[fari][neari][week][i]*Sigma_Far[farj][nearj][week][j];
//                }
//                if(neari==nearj && fari==farj)
//                {
//                    //General covariance
//                    CovStat[x][y][week]=(Sigma_Near[fari][neari][week][i]*Sigma_Near[fari][neari][week][j])+(Sigma_Far[farj][nearj][week][i]*Sigma_Far[farj][nearj][week][j]);
//                }
//                //Uncorrelated terms
//                if(neari!=nearj && fari!=farj)
//                {
//                    CovStat[x][y][week]=0;
//                }
//            }
//        }
//    }
//    Int_t x =0;
//    Int_t y =0;
//    for (Int_t week = 0; week<Nweeks; week++)
//    {
//        for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
//        {
//            //Logic for the 2D matrix index done up to 8 ADs
//            if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
//            if(neari==1){Ni2++;}
//            if(neari==2){Ni3++;}
//            if(neari==3){Ni4++;}
//            
//            for (Int_t fari=0; fari<ADsEH3; fari++)
//            {
//                //Logic for the 2D matrix index done up to 8 ADs
//                if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
//                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
//                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
//                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
//                
//                for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
//                {
//                    //Logic for the 2D matrix index done up to 8 ADs
//                    if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
//                    if(nearj==1){Nj2++;}
//                    if(nearj==2){Nj3++;}
//                    if(nearj==3){Nj4++;}
//                    
//                    for (Int_t farj=0; farj<ADsEH3; farj++)
//                    {
//                        //Logic for the 2D matrix index done up to 8 ADs
//                        if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
//                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
//                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
//                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
//                        for (Int_t i = 0; i<n_evis_bins; i++)
//                        {//columns
//                            for (Int_t j = 0; j<n_evis_bins; j++)
//                            {//rows
//                                x= i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
//                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
//                                
//                                if (CovStat[x][x][week]*CovStat[y][y][week]<0))
//                                {
//                                    NormCovStat[x][y][week]=-1*(CovStat[x][y][week])/(sqrt(-1*CovStat[x][x][week]*CovStat[y][y][week]));
//                                }
//                                else if ((CovStat[x][x][week]||CovStat[y][y][week])==0)
//                                {
//                                    NormCovStat[x][y][week]=0;
//                                }
//                                else
//                                {
//                                    NormCovStat[x][y][week]=(CovStat[x][y][week])/(sqrt(CovStat[x][x][week]*CovStat[y][y][week]));
//                                }
//                                
//                                StatCov2H->SetBinContent(x+1,y+1,NormCovStat[x][y][week]);
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//  // StatCov2H->Draw("colz");
}

void CovarianceMatrix3 :: RandomEnergyScaleMatrix()
{
    
}
void CovarianceMatrix3 :: RandomIAVMatrix()
{
    for (Int_t AD=0; AD<NADs; AD++)
    {
        IAVError[AD]= (1+IAVError*rand->Gaus(0,1)); //Each AD is varied individually.
    }
    
    for(Int_t i=0; i<MatrixBins; i++)
    {
        for(Int_t j=0; i<MatrixBins; i++)
        {
//            cout<<"IAV ORIGINAL" << IAVMatrix[i][j] << endl;
            
            IAVMatrix[i][j]=IAVError[AD]*IAVMatrix[i][j];
            
//            cout<< "IAV ALTERED" << IAVMatrix[i][j] << endl;
        }
    }
}
void CovarianceMatrix3 :: RandomNLMatrix()
{
    switch (NLModelE)
    {
        case 0://BCW NL Model
            //Randomize BCW nonlinear model parameters
            for (Int_t i = 0; i < 5; i++)
            {
                m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i] + m_bcw_elec_nl_par_error[i] * rand->Gaus(0,1);
                cout << " IF ANY OF THIS NUMBERS IS NEGATIVE I HAVE TO CHANGE THIS PART OF THE CODE, CHECK THAT NO PARAMETER IS NEGATIVE" << m_bcw_elec_nl_par[i] << endl;
            }
            break;
        case 1:
            //Randomize LBNL nonlinear model parameters
            for (Int_t i = 0; i < 3; i++)
            {
                m_lbnl_nl_par[i] = m_lbnl_nl_par_nominal[i] + m_lbnl_nl_par_error[i] * rand->Gaus(0,1);
                cout << " IF ANY OF THIS NUMBERS IS NEGATIVE I HAVE TO CHANGE THIS PART OF THE CODE, CHECK THAT NO PARAMETER IS NEGATIVE" <<  m_lbnl_nl_par[i]  << endl;
            }
            break;

        case 2:
            //Randomize IHEP nonlinear model parameters
            break;
        case 3:
//            //Randomize Unified NL model
//            Double_t ranvec_uniform[4];
//            Double_t ranvec[4];
//            
//            for (Int_t i = 0; i < 4; i++)
//            {
//                ranvec_uniform[i] = ran->Gaus(0,1);
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
        
            break;

        default:
    }
}
void CovarianceMatrix3 :: RandomResolutionMatrix()
{
    
}

void CovarianceMatrix3 :: RandomStatisticalMatrix()
{
                    
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
//        cout << "IAV Energy " << IAVEnergy[VisibleEnergyIndex] << endl;
//        cout << "IAV Energy Idx " << IAVEnergyIdx[VisibleEnergyIndex] << endl;
//        cout << " VisibleEnergyIndex" << VisibleEnergyIndex << endl;
        
        binScaling.push_back(func->Derivative(Energy[VisibleEnergyIndex]));
        
//        cout << "BIN SCALING" << binScaling[VisibleEnergyIndex] << endl;
        
        sign.push_back(1);
        dEtrue.push_back(Energy[VisibleEnergyIndex] - EnergyIdx[VisibleEnergyIndex]*BinWidth);
//        cout << "dEtrue" << dEtrue[VisibleEnergyIndex] << endl;
        if (dEtrue[VisibleEnergyIndex] < 0)
        {
            sign.push_back(-1);
        }
//        cout << "sign" << sign[VisibleEnergyIndex] << endl;
    }
}

void CovarianceMatrix3 :: CreateEnergyMatrix()
{
    TFile* SaveSpectrumDataF = TFile::Open("./CrossChecks/EnergyMatrixComparison/Energy240.root","recreate");

    // This should be used with the relative oscillation method
    for(int idet=0; idet<ADsEH3; idet++)
    {
        //Fill the Matrix
            EnergyMatrixH[idet] = new TH2F(Form("EvisEnu%d",idet),Form("EvisEnu%d",idet),n_etrue_bins,enu_bins,n_evis_bins,evis_bins);
            EnergyMatrixH[idet]->Reset();
            
            for (Int_t i = 0; i < n_etrue_bins; i++)
            {
                for (Int_t j = 0; j < n_evis_bins; j++)
                {
                    for(Int_t TrueEnergyIndex = Int_t((enu_bins[i])*(MatrixBins/FinalVisibleEnergy)); TrueEnergyIndex<=Int_t((enu_bins[i+1])*(MatrixBins/FinalVisibleEnergy))&&TrueEnergyIndex<MatrixBins; TrueEnergyIndex++)
                    {
                        for(Int_t VisibleBin = Int_t(evis_bins[j]*(MatrixBins/FinalVisibleEnergy)); VisibleBin<=Int_t(evis_bins[j+1]*(MatrixBins/FinalVisibleEnergy)); VisibleBin++)
                        {
                            EnergyMatrixH[idet]->SetBinContent(i+1,j+1,EnergyMatrixH[idet]->GetBinContent(i+1,j+1)+PositronVisibleSpectrumH[idet][TrueEnergyIndex]->GetBinContent(VisibleBin+1));
                        }
                    }
                }
            }
        // UNCOMMENT FOLLOWING LINES TO SAVE ENERGY MATRIX SLICES IN EACH PRODUCTION STEP
        for(Int_t TrueEnergyIndex = 0; TrueEnergyIndex<MatrixBins; TrueEnergyIndex++)
        {
            PositronTrueSpectrumH[idet][TrueEnergyIndex]->Write();
            PositronIAVSpectrumH[idet][TrueEnergyIndex]->Write();
            PositronNLSpectrumH[idet][TrueEnergyIndex]->Write();
            PositronVisibleSpectrumH[idet][TrueEnergyIndex]->Write();
        }

        EnergyMatrixH[idet]->Write();//Save Matrix before rebinning

        //Normalize the Matrix
        for (Int_t j = 0; j < n_evis_bins; j++)
        {
            Norma[j] = EnergyMatrixH[idet]->Integral(0,n_etrue_bins,j,j+1);
           // Norma[j] = EnergyMatrixH[idet]->ProjectionX("VisibleEnergySlide",j,j+1)->Integral();Both methods give the same solution, I assume the first is faster.

        }
        
        for (Int_t j = 0; j < n_evis_bins; j++)
        {
            for (Int_t i = 0; i < n_etrue_bins; i++)
            {
                if(Norma[j]!=0)
                {
                    EnergyMatrixH[idet]->SetBinContent(i+1,j+1,EnergyMatrixH[idet]->GetBinContent(i+1,j+1)/Norma[j]);//Normalization so Σi E(i,j) = 1; (Σ(x axis) =1)
                    //After checking the normalization the sum doesn't match exactly 1, but 0.999... There's some rounding error of a few 0.X% in the worst case <0.4%
                }
            }
        }
        //Save Matrix after rebinning
        EnergyMatrixH[idet]->Write();

        //Code to check correct normalization:
//        Double_t Normx[MaxNbins]={};
//        Double_t Normy[MaxNbins]={};
//        
//        for(Int_t i=0;i<n_etrue_bins;i++)
//        {
//            for(Int_t j=0;j<n_evis_bins;j++)
//            {
//                //Check normalization
//                Normy[i]=Normy[i]+EnergyMatrixH[idet]->GetBinContent(i+1,j+1);
//            }
//
//        }
//        for(Int_t j=0;j<n_evis_bins;j++)
//        {
//            for(Int_t i=0;i<n_etrue_bins;i++)
//            {
//                Normx[j]=Normx[j]+EnergyMatrixH[idet]->GetBinContent(i+1,j+1);
//            }
//        }
//        for(Int_t j=0;j<n_evis_bins;j++)
//        {
//            TH1D *px =EnergyMatrixH[idet]->ProjectionX("x",j,j+1);
//            cout << "Norma in X" << Normx[j]<<endl;
//            cout << "Norma in X " << px->Integral() <<endl;;
//        }
//        for(Int_t i=0;i<n_etrue_bins;i++)
//        {
//            TH1D *py =EnergyMatrixH[idet]->ProjectionY("y",i,i+1);
//            cout << "Norma in Y" << Normy[i]<<endl;
//            cout << "Norma in Y " << py->Integral() <<endl;;
//        }
    }
    SaveSpectrumDataF->Close();
}

void CovarianceMatrix3 :: GetEnergyShift(Int_t TrueEnergyIndex)
{
    Double_t NeutrinoEnergy = (TrueEnergyIndex)*BinWidth;
    PositronEnergy = VisibleF->Eval(NeutrinoEnergy);
    
    PositronEnergyIndex=(Int_t)(PositronEnergy/BinWidth);

    for(int idet=0; idet<ADsEH3; idet++)
    {
        //Calculate Energy Shift        
        PositronTrueSpectrumH[idet][TrueEnergyIndex]= new TH1F(Form("PositronTrueSpectrum%d,%d",idet,TrueEnergyIndex),Form("PositronTrueSpectrum%d,%d",idet,TrueEnergyIndex), TotalBins,InitialVisibleEnergy,FinalVisibleEnergy);
        
        //Reset values to 0
        if (PositronEnergy < 1.015)//is about the minimum visible energy (1.8-(939.565-938.272-0.51099) = 1.017999 for 0th order, 1.01574 for 1st order calculation)
        {
            for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
            {
                PositronTrueSpectrumH[idet][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, 0);
            }
        }
        else
        {
            PositronTrueSpectrumH[idet][TrueEnergyIndex]->SetBinContent(PositronEnergyIndex+1,OriginalPredictionH[idet][0][0]->GetBinContent(PositronEnergyIndex-Int_t(1.015*(MatrixBins/FinalVisibleEnergy))+1));
        }
    }
}

void CovarianceMatrix3 :: GetIAVShift(Int_t TrueEnergyIndex)
{
    for(int idet=0; idet<ADsEH3; idet++)
    {
        //IAV
        PositronIAVSpectrumH[idet][TrueEnergyIndex]=(TH1F*)PositronTrueSpectrumH[idet][TrueEnergyIndex]->Clone(Form("IAVSpectrum%d,%d",idet,TrueEnergyIndex));
        PositronIAVSpectrumH[idet][TrueEnergyIndex]->SetTitle("IAV Spectrum");
        //Calculate IAV Shift
        if (PositronEnergy >= 1.015)//1.8-0.78 is about the minimum visible energy (1.8-(939.565-938.272-0.51099) = 1.017999 for 0th order, 1.01574 for 1st order calculation)
        {
            //Reset values to 0
            for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
            {
               PositronIAVSpectrumH[idet][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, 0);
            }
            
            for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TrueEnergyIndex+1; VisibleEnergyIndex++)
            {
                PositronIAVSpectrumH[idet][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, PositronIAVSpectrumH[idet][TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1) + IAVMatrix[PositronEnergyIndex][VisibleEnergyIndex] * PositronTrueSpectrumH[idet][TrueEnergyIndex]->GetBinContent(PositronEnergyIndex+1));
            }
        }
    }
}

void CovarianceMatrix3 :: GetNLShift(Int_t TrueEnergyIndex)
{
    for(Int_t idet=0; idet<ADsEH3; idet++)
    {
        //NL
        PositronNLSpectrumH[idet][TrueEnergyIndex] = (TH1F*)PositronIAVSpectrumH[idet][TrueEnergyIndex]->Clone(Form("NLSpectrum%d,%d",idet,TrueEnergyIndex));
        PositronNLSpectrumH[idet][TrueEnergyIndex]->SetTitle("NL Spectrum");
        
        //Calculate NL Shift
        for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
        {
            if(EnergyIdx[VisibleEnergyIndex]==TotalBins-1)
            {
                dNdE.push_back(PositronIAVSpectrumH[idet][TrueEnergyIndex]->GetBinContent(EnergyIdx[VisibleEnergyIndex]+1));
            }
            else
            {
                dNdE.push_back((TMath::Abs(dEtrue[VisibleEnergyIndex])/BinWidth)*PositronIAVSpectrumH[idet][TrueEnergyIndex]->GetBinContent(EnergyIdx[VisibleEnergyIndex]+sign[VisibleEnergyIndex]+1)+(1 - TMath::Abs(dEtrue[VisibleEnergyIndex])/BinWidth)* PositronIAVSpectrumH[idet][TrueEnergyIndex]->GetBinContent(EnergyIdx[VisibleEnergyIndex]+1));
            }
        }
        for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
        {
            PositronNLSpectrumH[idet][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, dNdE[VisibleEnergyIndex]/binScaling[VisibleEnergyIndex]);
            dNdE.pop_back();//Without this line I keep pushing them down the stack            
        }
    }
}

void CovarianceMatrix3 :: GetResolutionShift(Int_t TrueEnergyIndex)
{
    for(Int_t idet=0; idet<ADsEH3; idet++)
    {
        //Resolution
        PositronVisibleSpectrumH[idet][TrueEnergyIndex] = (TH1F*)PositronNLSpectrumH[idet][TrueEnergyIndex]->Clone(Form("VisibleSpectrum%d,%d",idet,TrueEnergyIndex));
        PositronVisibleSpectrumH[idet][TrueEnergyIndex]->SetTitle("Visible Spectrum");

        //Calculate Resolution effect
        for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
        {
            Double_t sigma = (reso_func->Eval(VisibleEnergyIndex*BinWidth) + m_detectorResolution_bias[idet]) * VisibleEnergyIndex*BinWidth;//Sigma = FReso(Energy) * Energy
            Double_t minDetE = VisibleEnergyIndex*BinWidth - resolutionRange*sigma;
            Double_t maxDetE = VisibleEnergyIndex*BinWidth + resolutionRange*sigma;
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
                
                PositronVisibleSpectrumH[idet][TrueEnergyIndex]->SetBinContent(detIdx+1, (PositronVisibleSpectrumH[idet][TrueEnergyIndex]->GetBinContent(detIdx+1)+PositronNLSpectrumH[idet][TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1)*gausFactor));
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

TF1* CovarianceMatrix3 :: GetNeutrinoToVisibleFunction(Int_t order)
{
    if(order==0)
    {
        Double_t Correction = Mn-Mp-Me;
        TF1 *VisibleF = new TF1("VisibleF",VisibleEnergy0F,InitialEnergy,FinalVisibleEnergy,1);
        VisibleF->SetParameter(0,Correction);
        cout<<"Zeroth order"<<endl;
    }
    if(order==1)
    {
        TF1 *VisibleF = new TF1("VisibleF",VisibleEnergy1F,InitialEnergy,FinalVisibleEnergy,1);
        VisibleF->SetParameter(0,0);//Maybe I can improve this using a MC simulation of the angular distribution calculated in http://authors.library.caltech.edu/2796/1/VOGprd99.pdf
        cout<<"First order"<<endl;
    }
    return VisibleF;
}

void CovarianceMatrix3 :: LoadNLParameters()
{
    switch (NLModelE)
    {
        case 0://BCW NL Model
            cout << "LOAD BCW NL PARAMETERS" << endl;
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
            }
            bcw_elec_data.close();
            
            for (Int_t i = 3; i < 5; i++)
            { // those describe additional uncertainty for bcw model that are described in Doc-XXXX
                m_bcw_elec_nl_par[i] = 0;
                m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];
                m_bcw_elec_nl_par_error[i] = 1;
            }
            
            TFile *bcw_ele_err_file = new TFile("bcw_nl_data/ele_err.root");
            g_bcw_elec_nl_error[0] = (TGraph*)bcw_ele_err_file->Get("g_up")->Clone();
            g_bcw_elec_nl_error[1] = (TGraph*)bcw_ele_err_file->Get("g_down")->Clone();
            bcw_ele_err_file->Close();
            break;
        case 1://LBNL NL Model
            cout << "LOAD LBNL NL PARAMETERS" << endl;
            ifstream lbnl_positron_data("lbnl_nl_data/lbnl_positron_nl.txt");
            if (!lbnl_positron_data.is_open())
            {
                cout << "Error: cannot find LBNL non-linearity curve!!!" << endl;
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
            break;
        case 2://IHEP NL Model
            break;
        case 3://Unified NL Model
            
            for (Int_t i = 0; i < m_num_unified_nl_pars; i++)
            {
                m_unified_nl_par[i] = 0.0;
                m_unified_nl_par_nominal[i] = 0.0;
                m_unified_nl_par_error[i] = 1.0;
            }
            
            TFile *unified_nl_file = new TFile("unified_nl_data/nl_models_final.root");
            if (!unified_nl_file->IsOpen())
            {
                cout << "Error: cannot find the unified non-linearity curve!!!" << endl;
                exit(0);
            }
            
            // Use IHEP I as the nominal model
            g_unified_positron_nl = (TGraph*)unified_nl_file->Get("positron_0")->Clone();
            g_unified_positron_nl_pulls[0] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",1))->Clone();
            g_unified_positron_nl_pulls[1] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",2))->Clone();
            g_unified_positron_nl_pulls[2] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",3))->Clone();
            g_unified_positron_nl_pulls[3] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",4))->Clone();
            
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
            break;
        default:
    }
}

void CovarianceMatrix3 :: LoadResponseMatrix()
{
    switch (NLModelE)
    {
        case 0://BCW NL Model
            TFile* ResponseF = TFile::Open("./NominalResponseMatrices/NominalResponseBCWModel.root");
            break;
        case 1://LBNL NL Model
            TFile* ResponseF = TFile::Open("./NominalResponseMatrices/NominalResponseLBNLModel.root");
            break;
        case 2://IHEP NL Model
            TFile* ResponseF = TFile::Open("./NominalResponseMatrices/NominalResponseIHEPModel.root");
            break;
        case 3://Unified NL Model
            TFile* ResponseF = TFile::Open("./NominalResponseMatrices/NominalResponseUnifiedModel.root");
            break;
        default:
    }

    
    for (Int_t AD =0; AD<ADsEH3; AD++)
    {
        NominalResponseMatrixH[AD] = (TH2F*)gDirectory->Get(Form("EvisEnu%i;2",AD));
    }
    
    ResponseF->Close();
}

void CovarianceMatrix3 :: LoadBackgrounds()
{
    TFile* BackgroundsF = TFile::Open("./BackgroundSpectrum/Backgrounds.root");
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t AD =0; AD<NADs; AD++)
        {
            AccidentalsH[AD][week]= (TH1F*)gDirectory->Get(Form("Accidentals_AD%i",AD+1));
            LiHeH[AD][week]=(TH1F*)gDirectory->Get("hist_Bkg_StrongAmC");//Missing LiHe inputs so far
            FastNeutronsH[AD][week]=(TH1F*)gDirectory->Get("FN");
            AmCH[AD][week]=(TH1F*)gDirectory->Get("hist_Bkg_StrongAmC");
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
            ScaleFactorAccidental[AD][week]=1;
            ScaleFactorLiHe[AD][week]=1;
            ScaleFactorFastNeutrons[AD][week]=1;
            ScaleFactorAmC[AD][week]=1;
            
            hall=2;
            if(AD<ADsEH1)
            {
                hall=1;
            }
            if(AD>=ADsEH1+ADsEH2)
            {
                hall=3;
            }
            
            if(AccidentalMatrix)
            {
                rand->SetSeed(0);
                ScaleFactorAccidental[AD][week]=(1+HAccidentalError[AD]*rand->Gaus(0,1));
            }
            if(LiHeMatrix)
            {
                rand->SetSeed(0);
                ScaleFactorLiHe[AD][week]=(1+HLiHeError[hall]*rand->Gaus(0,1));
            }
            if(FastNeutronsMatrix)
            {
                rand->SetSeed(0);
                ScaleFactorFastNeutrons[AD][week]=(1+HFastNeutronsError[hall]*rand->Gaus(0,1));
            }
            if(AmCMatrix)
            {
                rand->SetSeed(0);
                ScaleFactorAmC[AD][week]=(1+HAmCError[hall]*rand->Gaus(0,1));
            }
        }
        
        //Accidentals uncorrelated
        //LiHe and Fast neutrons correlated by hall
        
        ScaleFactorLiHe[1][week]=ScaleFactorLiHe[0][week];
        ScaleFactorLiHe[4][week]=ScaleFactorLiHe[3][week];
        ScaleFactorLiHe[5][week]=ScaleFactorLiHe[3][week];
        ScaleFactorFastNeutrons[1][week]=ScaleFactorFastNeutrons[0][week];
        ScaleFactorFastNeutrons[4][week]=ScaleFactorFastNeutrons[3][week];
        ScaleFactorFastNeutrons[5][week]=ScaleFactorFastNeutrons[3][week];
        
        //AmC all correlated
        
        ScaleFactorAmC[1][week]=ScaleFactorAmC[0][week];
        ScaleFactorAmC[2][week]=ScaleFactorAmC[0][week];
        ScaleFactorAmC[3][week]=ScaleFactorAmC[0][week];
        ScaleFactorAmC[4][week]=ScaleFactorAmC[0][week];
        ScaleFactorAmC[5][week]=ScaleFactorAmC[0][week];
        
        for (Int_t AD =0; AD<NADs; AD++)
        {
            RandomAccidentalsH[AD][week]=(TH1F*)AccidentalsH[AD][week]->Clone();
            RandomLiHeH[AD][week]=(TH1F*)LiHeH[AD][week]->Clone();
            RandomFastNeutronsH[AD][week]=(TH1F*)FastNeutronsH[AD][week]->Clone();
            RandomAmCH[AD][week]=(TH1F*)AmCH[AD][week]->Clone();
            
            if(VariateRate)
            {
                RandomAccidentalsH[AD][week]->Scale(ScaleFactorAccidental[AD][week]);
                RandomLiHeH[AD][week]->Scale(ScaleFactorLiHe[AD][week]);
                RandomFastNeutronsH[AD][week]->Scale(ScaleFactorFastNeutrons[AD][week]);
                RandomAmCH[AD][week]->Scale(ScaleFactorAmC[AD][week]);
                RandomAmCH[AD][week]->Draw("same");
                //            RandomLiHeH[AD][week]->Draw("same");
                //           RandomFastNeutronsH[AD][week]->Draw("same");
                // RandomAccidentalsH[AD][week]->Draw("same");
            }
        }
        
        if(DistortBackgrounds)
        {
            if(AccidentalMatrix)
            {
                TF1* func_acc=GetDistortionFunction(DistortAcc);
                
                for(Int_t iAD=0;iAD<NADs;iAD++)
                {
                    RandomAccidentalsH[iAD][week]->Multiply(func_acc);
                    RandomAccidentalsH[iAD][week]->Scale(RandomAccidentalsH[iAD][week]->Integral()/AccidentalsH[iAD][week]->Integral());
                    
                    //      cout<<RandomAccidentalsH[AD][week]->Integral()/AccidentalsH[AD][week]->Integral()<<"\n";
                    //      func_acc->Draw();
                }
                //                    RandomAccidentalsH[0][week]->Draw("same");
                
                delete func_acc;
            }
            if(LiHeMatrix)
            {
                TF1* func_LiHe=GetDistortionFunction(DistortLiHe);
                for(Int_t iAD=0;iAD<NADs;iAD++)
                {
                    RandomLiHeH[iAD][week]->Multiply(func_LiHe);
                    RandomLiHeH[iAD][week]->Scale(RandomLiHeH[iAD][week]->Integral()/LiHeH[iAD][week]->Integral());
                    //                    cout<<RandomLiHeH[AD][week]->Integral()/LiHeH[AD][week]->Integral()<<"\n";
                    //      func_LiHe->Draw("same");
                }
                delete func_LiHe;
            }
            if(FastNeutronsMatrix)
            {
                //FN shape distortions are applied taking into account correlations in each EH
                //EH1
                
                TF1* func_FN=GetFastNeutronsDistortionFunction(DistortFN);
                for(Int_t iAD=0;iAD<ADsEH1;iAD++)
                {
                    RandomFastNeutronsH[iAD][week]->Multiply(func_FN);
                    RandomFastNeutronsH[iAD][week]->Scale(RandomFastNeutronsH[iAD][week]->Integral()/FastNeutronsH[iAD][week]->Integral());
                }
                //EH2
                
                func_FN=GetFastNeutronsDistortionFunction(DistortFN);
                for(Int_t iAD=ADsEH1;iAD<ADsEH2+ADsEH1;iAD++)
                {
                    RandomFastNeutronsH[iAD][week]->Multiply(func_FN);
                    RandomFastNeutronsH[iAD][week]->Scale(RandomFastNeutronsH[iAD][week]->Integral()/FastNeutronsH[iAD][week]->Integral());
                }
                //EH3
                
                func_FN=GetFastNeutronsDistortionFunction(DistortFN);
                for(Int_t iAD=ADsEH1+ADsEH2;iAD<ADsEH1+ADsEH2+ADsEH3;iAD++)
                {
                    RandomFastNeutronsH[iAD][week]->Multiply(func_FN);
                    RandomFastNeutronsH[iAD][week]->Scale(RandomFastNeutronsH[iAD][week]->Integral()/FastNeutronsH[iAD][week]->Integral());
                }
                delete func_FN;
            }
            if(AmCMatrix)
            {
                TF1* func_AmC=GetDistortionFunction(DistortAmC);
                for(Int_t iAD=0;iAD<NADs;iAD++)
                {
                    RandomAmCH[iAD][week]->Multiply(func_AmC);
                    RandomAmCH[iAD][week]->Scale(RandomAmCH[iAD][week]->Integral()/AmCH[iAD][week]->Integral());
                    //RandomAmCH[iAD][week]->Draw("same");
                    //func_AmC->Draw("same");
                }
                delete func_AmC;
            }
        }
    }
}

void CovarianceMatrix3 :: AddBackgrounds(TH1F* Histogram, TH1F* Backgrounds)
{
    for (Int_t i = 0; i<4; i++)
    {
        Histogram->Add(Backgrounds);
    }
}

void CovarianceMatrix3 :: LoadPredictions()
{
    TFile* FarHallPredictionsF = TFile::Open("./RootOutputs/FarSpectrumFraction.root");

    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            for (Int_t far =0; far<ADsEH3; far++)
            {
                OriginalPredictionH[far][near][week] = (TH1F*)gDirectory->Get(Form("AD%i Far Spectrum prediction from near AD%i",far+1,near+1));
                PredictionTrueH[far][near][week]= new TH1F(OriginalPredictionH[far][near][week]->GetName(), OriginalPredictionH[far][near][week]->GetTitle(),n_etrue_bins,enu_bins);

                for (Int_t i = 0; i < n_etrue_bins; i++)
                {
                    if(enu_bins[i]>=InitialEnergy)
                    {
                        for(Int_t TrueBin=Int_t((enu_bins[i]-InitialEnergy)*OriginalPredictionH[far][near][week]->GetXaxis()->GetNbins()/(FinalEnergy-InitialEnergy)); TrueBin<=Int_t((enu_bins[i+1]-InitialEnergy)*OriginalPredictionH[far][near][week]->GetXaxis()->GetNbins()/(FinalEnergy-InitialEnergy)); TrueBin++)
                        {
                            PredictionTrueH[far][near][week]->SetBinContent(i+1, PredictionTrueH[far][near][week]->GetBinContent(i+1) + OriginalPredictionH[far][near][week]->GetBinContent(TrueBin+1));
                        }
                    }
                }
            }
        }
    }
    FarHallPredictionsF->Close();
}

void CovarianceMatrix3 :: CombinePredictions()
{
    for (Int_t week = 0; week<Nweeks; week++)
    {
        PredictionLA[week] = new TH1F("Combined Spectrum in LA site","Combined Spectrum in LA site",n_evis_bins,evis_bins);
        PredictionDB[week] = new TH1F("Combined Spectrum in DB site","Combined Spectrum in LA site",n_evis_bins,evis_bins);
        PredictionLAWithBkgd[week] = new TH1F(Form("Combined Far Spectrum in LA site with variations due to Systematic%d or Background%d",SystematicE,BackgroundE),Form("Combined Far Spectrum from LA site with variations due to Systematic%d or Background%d",SystematicE,BackgroundE),n_evis_bins,evis_bins);
        PredictionDBWithBkgd[week] = new TH1F(Form("Combined Far Spectrum from DB site with variations due to Systematic%d or Background%d",SystematicE,BackgroundE),Form("Combined Far Spectrum from DB site with variations due to Systematic%d or Background%d",SystematicE,BackgroundE),n_evis_bins,evis_bins);
        
        for (Int_t near =0; near<(ADsEH1); near++)//EH1
        {
            for (Int_t far =0; far<ADsEH3; far++)
            {                
                for (Int_t i = 0; i < n_etrue_bins; i++)
                {
                     PredictionDB[week]->Add(PredictionVisH[far][near][week]);
                     PredictionDBWithBkgd[week]->Add(PredictionVisWithBkgdH[far][near][week]);
                }
            }
        }
        for (Int_t near =0; near<(ADsEH2); near++)//EH2
        {
            for (Int_t far =0; far<ADsEH3; far++)
            {
                for (Int_t i = 0; i < n_etrue_bins; i++)
                {
                    PredictionLA[week]->Add(PredictionVisH[far][near][week]);
                    PredictionLAWithBkgd[week]->Add(PredictionVisWithBkgdH[far][near][week]);
                }
            }
        }

    }
}

TH2F* CovarianceMatrix3 :: NormalizeCovariance(Double_t Cov[][MaxNbins][MaxPeriods])
{
    TH2F* Cov2H = new TH2F("Covariance Matrix","Covariance Matrix",n_evis_bins,evis_bins,n_evis_bins,evis_bins);

    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t i = 0; i<n_evis_bins; i++)
        {//columns
            for (Int_t j = 0; j<n_evis_bins; j++)
            {//rows
                if((Cov[i][i][week]||Cov[j][j][week])==0)
                {
                    NormCov[i][j][week]=0;//To avoid (0/0) nans when bin contents are practically the same; (this happens when backgrounds are not varied)
                    cout << "Norm cov tried to be inf" << endl;
                }
                else if((Cov[i][i][week]*Cov[j][j][week])<0))
                {
                    NormCov[i][j][week]=-1*(Cov[i][j][week])/(sqrt(-1*Cov[i][i][week]*Cov[j][j][week]));
                    cout << "Negative square root" << endl;
                }
                else
                {
                    NormCov[i][j][week]=(Cov[i][j][week])/(sqrt(Cov[i][i][week]*Cov[j][j][week]));
                }
                Cov2H->SetBinContent(i+1,j+1,NormCov[i][j][week]);
            }
        }
    }
    
    //  cout <<NormCov[x][y][week];
    // cout << " " << Cov[x][y][week] << " " <<  Cov[x][x][week] << " " << Cov[y][y][week] << " " << sqrt(Cov[x][x][week]*Cov[y][y][week]) << "\n";
    
    //    Cov2H->Draw("colz");
    return Cov2H;
}

void CovarianceMatrix3 :: LoadNearHall()
{
    Char_t filenameNear[1024];

    TFile* NearHallDataF = TFile::Open("./RootOutputs/NominalOutputs/Oscillation.root"); //This should be real data. It has to be fixed.
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            sprintf(filenameNear,"Total spectrum after oscillation at AD%i",near+1);
            NearHallDataF->cd("Total AD Spectra after oscillation");
            OriginalNearHallSpectrumH[near][week] = (TH1F*)gDirectory->Get(filenameNear);
            NearHallSpectrumTrueH[near][week]= new TH1F(OriginalNearHallSpectrumH[near][week]->GetName(),OriginalNearHallSpectrumH[near][week]->GetTitle(),n_etrue_bins,enu_bins);

            for (Int_t i = 0; i < n_etrue_bins; i++)
            {
                if(enu_bins[i]>=InitialEnergy)
                {
                    for(Int_t TrueBin=Int_t((enu_bins[i]-InitialEnergy)*OriginalNearHallSpectrumH[near][week]->GetXaxis()->GetNbins()/(FinalEnergy-InitialEnergy)); TrueBin<=Int_t((enu_bins[i+1]-InitialEnergy)*OriginalNearHallSpectrumH[near][week]->GetXaxis()->GetNbins()/(FinalEnergy-InitialEnergy)); TrueBin++)
                    {
                            NearHallSpectrumTrueH[near][week]->SetBinContent(i+1, NearHallSpectrumTrueH[near][week]->GetBinContent(i+1)+ OriginalNearHallSpectrumH[near][week]->GetBinContent(TrueBin+1));
                    }
                }
            }
        }
    }
    NearHallDataF->Close();
    
    TFile* SpectrumF = TFile::Open("./RootOutputs/TrueSpectrum.root","recreate");
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            NearHallSpectrumTrueH[near][week]->Write();

            for (Int_t far =0; far<ADsEH3; far++)
            {
                PredictionTrueH[far][near][week]->Write();
            }
        }
    }
    SpectrumF->Close();
}

void CovarianceMatrix3 :: LoadIavCorrection()//From Bryce Littlejohn's results. //I don't use different IAV matrix for each AD since the analysis it's been done for just 1 of them, assume identical ADs.
{
    TFile * f = new TFile("./IavDistortion/IAVDistortion.root");
    TH2F * Correction = (TH2F*)f->Get("Correction");
    
    cout << "Reading IAV correction file" << endl;
    
    for(Int_t i=0;i<TotalBins;i++)
    { // i: true positron energy bin; j: distorted energy bin
        //Total events in input spectrum bin i
        if (Correction->Integral(i+1,i+1,0,-1) > 0)
        {
            for(Int_t j=0;j<i+1;j++)
            {
                IAVMatrix[i][j]= Correction->GetBinContent(i+1,j+1)/Correction->Integral(i+1,i+1,0,-1);
            }
        }
        else
        {
            for(Int_t j=0;j<i+1;j++)
            {
                if (i==j)
                {
                    IAVMatrix[i][j] = 1;
                }
                else
                {
                    IAVMatrix[i][j] = 0;
                }
            }
        }
    }
    f->Close();
}

void CovarianceMatrix3 :: SaveTotalSpectrum()
{
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            NearHallVisWithBkgdH[near][week]->Add(BackgroundSpectrumH[near][week]);

            for (Int_t far =0; far<ADsEH3; far++)
            {
                PredictionVisWithBkgdH[far][near][week]->Add(BackgroundSpectrumH[ADsEH1+ADsEH2+far][week]);
            }
        }
     }

    TFile* SaveSpectrumDataF = TFile::Open("./RootOutputs/SpectrumWithBackgrounds.root","recreate");
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            NearHallVisWithBkgdH[near][week]->Write();
            
            for (Int_t far =0; far<ADsEH3; far++)
            {
                PredictionVisWithBkgdH[far][near][week]->Write();
            }
        }
    }
   SaveSpectrumDataF->Close();
}

void CovarianceMatrix3 :: SaveCovarianceMatrix()
{
    Char_t filenameCov[1024];
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
        default:
    }
    switch (SystematicE)
    {
        case 0://Vary Reactor
            sprintf(filenameCov,"./CovarianceMatrices/ReactorCovarianceMatrix.root");
            CovMatrix2H->SetName("Reactor Covariance Matrix");
            CovMatrix2H->SetTitle("Reactor Covariance Matrix");
            break;
        case 1://Vary Energy Scale
            sprintf(filenameCov,"./CovarianceMatrices/EnergyScaleCovarianceMatrix.root");
            CovMatrix2H->SetName("EnergyScale Covariance Matrix");
            CovMatrix2H->SetTitle("EnergyScale Covariance Matrix");
            break;
        case 2://Vary IAV
            sprintf(filenameCov,"./CovarianceMatrices/IAVCovarianceMatrix.root");
            CovMatrix2H->SetName("IAV Covariance Matrix");
            CovMatrix2H->SetTitle("IAV Covariance Matrix");
            break;
        case 3://Vary NL
            sprintf(filenameCov,"./CovarianceMatrices/NLCovarianceMatrix.root");
            CovMatrix2H->SetName("NL Covariance Matrix");
            CovMatrix2H->SetTitle("NL Covariance Matrix");
            break;
        case 4://Vary Resolution
            sprintf(filenameCov,"./CovarianceMatrices/ResolutionCovarianceMatrix.root");
            CovMatrix2H->SetName("Resolution Covariance Matrix");
            CovMatrix2H->SetTitle("Resolution Covariance Matrix");
            break;
        case 5://Statistical
            //Save statistical matrix
            TFile* SaveStatisticalCovMatrixF = TFile::Open("./CovarianceMatrices/StatisticalCovarianceMatrix.root","recreate");
            StatCov2H->Write();
            SaveStatisticalCovMatrixF->Close();
            if(WriteOutput)
            {
                ofstream statf("CovarianceMatrices/StatisticalCovarianceMatrix.txt");
                for (Int_t week = 0; week<Nweeks; week++)
                {
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {//columns
                        for (Int_t j = 0; j<n_evis_bins; j++)
                        {//rows
                            statf << NormCovStat[i][j][week] << " ";
                        }
                    }
                    statf << endl;
                }
                statf.close();
            }
            break;
        default:
    }
    
    if(SystematicE!=5)
    {
    TFile* SaveCovarianceMatrixF = TFile::Open(filenameCov,"recreate");
    // CovMatrix2H->Draw("colz");
        CovMatrix2H->Write();
    SaveCovarianceMatrixF->Close();
    }
    //Save in a txt file
    if(WriteOutput)
    {
        ofstream covf("CovarianceMatrices/CovarianceMatrix.txt");
        for (Int_t i = 0; i<n_evis_bins; i++)
        {//columns
            for (Int_t j = 0; j<n_evis_bins; j++)
            {//rows
                covf << CovMatrix2H->GetBinContent(i+1,j+1) << " ";
            }
        }
     
        covf << endl;
        covf.close();
    }
}

void CovarianceMatrix3 :: SetVaryAccidentalMatrix(bool AccMatrix)
{
    BackgroundE = -1;//Reset BackgroundE
    AccidentalMatrix=AccMatrix;
    if(AccMatrix)
    {
        AddBackgrounds=1;//To add backgrounds
        VariateRate = 1;
        DistortBackgrounds = 0;
        BackgroundE = VaryAccidentalE;
        DistortAcc=0.2;
    }
}

void CovarianceMatrix3 :: SetVaryLiHeMatrix(bool LiMatrix)
{
    LiHeMatrix=LiMatrix;
    if(LiMatrix)
    {
        AddBackgrounds=1;//To add backgrounds
        VariateRate = 1;
        DistortBackgrounds = 0;

        BackgroundE = VaryLiHeE;
        DistortLiHe=0.2;
    }
}

void CovarianceMatrix3 :: SetVaryFastNeutronsMatrix(bool FNmatrix)
{
  
    FastNeutronsMatrix=FNmatrix;
    if(FNmatrix)
    {
        AddBackgrounds=1;//To add backgrounds
        VariateRate = 1;
        DistortBackgrounds = 0;
        BackgroundE = VaryFastNeutronsE;
        DistortFN=0.2;
    }
}

void CovarianceMatrix3 :: SetVaryAmCMatrix(bool AmCmatrix)
{
   
    AmCMatrix=AmCmatrix;
    if(AmCMatrix)
    {
        AddBackgrounds=1;//To add backgrounds
        VariateRate = 1;
        DistortBackgrounds = 0;
        BackgroundE = VaryAmCE;
        DistortAmC=0.2;
    }
}

void CovarianceMatrix3 :: SetDistortLiHeMatrix(bool LiMatrix)
{
  
    LiHeMatrix=LiMatrix;
    if(LiMatrix)
    {
        AddBackgrounds=1;//To add backgrounds
        VariateRate = 0;
        DistortBackgrounds=1;
        BackgroundE = DistortLiHeE;
        DistortLiHe=0.2;
    }
}

void CovarianceMatrix3 :: SetDistortFastNeutronsMatrix(bool FNmatrix)
{
    FastNeutronsMatrix=FNmatrix;
    if(FNmatrix)
    {
        AddBackgrounds=1;//To add backgrounds
        VariateRate = 0;
        DistortBackgrounds=1;
        BackgroundE = DistortFastNeutronsE;
        DistortFN=0.2;
    }
}

void CovarianceMatrix3 :: SetDistortAmCMatrix(bool AmCmatrix)
{
    
    AmCMatrix=AmCmatrix;
    if(AmCMatrix)
    {
        AddBackgrounds=1;//To add backgrounds
        VariateRate = 0;
        DistortBackgrounds=1;
        BackgroundE = DistortAmCE;
        DistortAmC=0.2;
    }
}

void CovarianceMatrix3 :: SetIsotopeMatrix(bool IsotopeMatrix)
{
    if(IsotopeMatrix)
    {
        AddBackgrounds=0;//To add backgrounds
        VariateRate = 0;
        DistortBackgrounds=0;
        Pred.SetRandomIsotopeFraction(IsotopeMatrix);
        SystematicE = ReactorE;
    }
}

void CovarianceMatrix3 :: SetReactorPowerMatrix(bool ReactorPowerMatrix)
{
    if(ReactorPowerMatrix)
    {
        AddBackgrounds=0;//To add backgrounds
        VariateRate = 0;
        DistortBackgrounds=0;
        Pred.SetRandomReactorPower(ReactorPowerMatrix);
        SystematicE = ReactorE;
    }
}

void CovarianceMatrix3 :: SetEnergyScaleMatrix(bool EnergyScaleMatrix)
{
    if(EnergyScaleMatrix)
    {
        AddBackgrounds=0;//To add backgrounds
        VariateRate = 0;
        DistortBackgrounds=0;
        RandomEnergyScaleMatrix();
        SystematicE = EnergyE;
    }
}

void CovarianceMatrix3 :: SetIAVMatrix(bool IAVMatrix)
{
    if(IAVMatrix)
    {
        AddBackgrounds=0;//To add backgrounds
        VariateRate = 0;
        DistortBackgrounds=0;
        RandomIAVMatrix();
        SystematicE = IAVE;
    }
}

void CovarianceMatrix3 :: SetNLMatrix(bool NLMatrix)
{
    if(NLMatrix)
    {
        AddBackgrounds=0;//To add backgrounds
        VariateRate = 0;
        DistortBackgrounds=0;
        RandomNLMatrix();
        SystematicE = NLE;
    }
}

void CovarianceMatrix3 :: SetResolutionMatrix(bool ResolutionMatrix)
{
    if(ResolutionMatrix)
    {
        AddBackgrounds=0;//To add backgrounds
        VariateRate = 0;
        DistortBackgrounds=0;
        RandomResolutionMatrix();
        SystematicE = ResolutionE;
    }
}
                                              
void CovarianceMatrix3 :: SetStatisticalMatrix(bool StatisticalMatrix)
{
   if(StatisticalMatrix)
   {
       AddBackgrounds=1;//To add backgrounds
       VariateRate = 0;
       DistortBackgrounds=0;
       
       SystematicE = StatisticalE;
   }
}

void CovarianceMatrix3 :: SetBCWModel(bool bcw)
{
    if(bcw)
    {
        NLModelE=BCWE;
    }
}
    
void CovarianceMatrix3 :: SetLBNLModel(bool lbnl)
{
    if(lbnl)
    {
        NLModelE=LBNLE;
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
    }
}

TF1* CovarianceMatrix3 :: GetNeutrinoToVisibleFunction(Int_t order)
{
    if(order==0)
    {
        Double_t Correction = Mn-Mp-Me;
        TF1 *VisibleF = new TF1("VisibleF",VisibleEnergy0F,InitialEnergy,FinalVisibleEnergy,1);
        VisibleF->SetParameter(0,Correction);
        cout<<"Zeroth order"<<endl;
    }
    if(order==1)
    {
        TF1 *VisibleF = new TF1("VisibleF",VisibleEnergy1F,InitialEnergy,FinalVisibleEnergy,1);
        VisibleF->SetParameter(0,0);//Maybe I can improve this using a MC simulation of the angular distribution calculated in http://authors.library.caltech.edu/2796/1/VOGprd99.pdf
        cout<<"First order"<<endl;
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
    func->Draw();
    return func;
}

// First order true to visible function.
//(x-(Mn-Mp))*(1 - x/Mn*(1.0 - [0]*sqrt(1 - Me*Me/(x-(Mn-Mp))/(x-(Mn-Mp)))))- ((Mn-Mp)*(Mn-Mp) - Me*Me)/(2*Mn)-Me]

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Formula is different from previous version (Mp instead of Mn in 1st order, check which is correct)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t VisibleEnergy1F(Double_t* x, Double_t* par)
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

Double_t VisibleEnergy0F(Double_t* x, Double_t* par)
{
    return x[0]-par[0];
}

Double_t CovarianceMatrix3 :: reso_func_bcw(Double_t* x, Double_t* par)
{
    Double_t e_orig = x[0];
    Double_t e_sigma = 1.0;
    
    if (e_orig > 0)//To avoid NaNs.
    {
        e_sigma = TMath::Sqrt(par[0]*par[0] + par[1]*par[1]/e_orig + par[2]*par[2]/e_orig/e_orig);
    }
    
    return e_sigma;
}

Double_t CovarianceMatrix3 :: nl_func_bcw(Double_t* x, Double_t* par)
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
Double_t CovarianceMatrix3 :: nl_func_lbnl(Double_t * x, Double_t * par)
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
    //  cout << "NL " << e_positron_true << " " << visibleE/e_positron_true << endl;
    
    return visibleE  * escale_par + escale_offset;
}

Double_t CovarianceMatrix3 :: nl_func_unified(Double_t * x, Double_t * par)
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