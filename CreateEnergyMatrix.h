#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include <vector>
#include <math.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TCanvas.h"
#include "OscillationReactor.h"
#include "NominalData.h"

const Double_t m_nSamples = 240;
const Double_t m_binWidth = 12/m_nSamples;
const Double_t m_eMin = 0;

//const bool NormalizeRow = 0; // Column normalization is the correct one
const bool RebinEnergyMatrix = 1; // To check that the fine matrix multiplication works in the same way than the coarse one. There's difference, and it seems best to apply the rebin version. 
const bool FlatEnergyMatrix = 1;// This one shouldn't make any change since the matrix is normalized. DayaBay asked to use a flat spectrum to produce the EvisEnu matrix. Yasu says that for the kind of apporach we use it is more accurate to use the best predictiondetector shape.

////NL
//const Int_t n_bcw_positron_nl = 1000;
//const Int_t n_unified_nl_points = 500;
//const Int_t m_num_unified_nl_pars = 4; // 4 marginal curves in the final model
//
////Particle masses
//const Double_t Me = 0.510999; // MeV
//const Double_t Mn = 939.565; // MeV
//const Double_t Mp = 938.272; // MeV

class CreateEnergyMatrix
{
private:
    NominalData* Data;
    std::string AnalysisString;

    bool BCW;
    bool LBNL;
    bool Unified;
    
    bool IsotopeMatrix;
    bool ReactorPowerMatrix;
    bool RelativeEnergyOffsetMatrix;
    bool AbsoluteEnergyOffsetMatrix;
    bool AbsoluteEnergyScaleMatrix;
    bool RelativeEnergyScaleMatrix;
    bool IAVMatrixb;
    bool NLMatrix;
    bool ResolutionMatrix;
    
    enum NLModel{BCWE, LBNLE, UnifiedE};//Add here new NL models
    NLModel NLModelE;
    
    //AD configuration parameters:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    
    bool analysis;
    //Binning parameters:
    bool LinearBinning;
    
    Int_t LimitTrue;
    Int_t LimitVis;
    
    Int_t n_evis_bins;
    Int_t n_etrue_bins;
    
    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    Int_t Nweeks;
    Double_t BinWidth;
    
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.
    
    Double_t DetectorEfficiency[MaxDetectors][MaxPeriods];
    
    //Response Matrix:
    Double_t PositronEnergy;
    Int_t PositronEnergyIndex;
    
    Double_t Norma[MatrixBins];
    Double_t NormaPos[MatrixBins];
    Double_t NormaIAV[MatrixBins];
    Double_t NormaNL[MatrixBins];
    Double_t NormaReso[MatrixBins];
    Double_t NormaTrans[MatrixBins];
    Double_t NormaInv[MatrixBins];
    
    Double_t FineNorma[MatrixBins];
    Double_t FineNormaPos[MatrixBins];
    Double_t FineNormaIAV[MatrixBins];
    Double_t FineNormaNL[MatrixBins];
    Double_t FineNormaReso[MatrixBins];
    
    TH1D* OscDeltaPositronSpectrumH[MatrixBins];
    TH1D* OscDeltaIAVSpectrumH[MatrixBins];
    TH1D* OscDeltaNLSpectrumH[MatrixBins];
    TH1D* OscDeltaVisibleSpectrumH[MatrixBins];
    
    TH1D* OscDeltaPositronSpectrumSumH[MatrixBins];
    TH1D* OscDeltaIAVSpectrumSumH[MatrixBins];
    TH1D* OscDeltaNLSpectrumSumH[MatrixBins];
    TH1D* OscDeltaVisibleSpectrumSumH[MatrixBins];
    
    TH1D* RebinnedOscDeltaPositronSpectrumSumH[MatrixBins];
    TH1D* RebinnedOscDeltaIAVSpectrumSumH[MatrixBins];
    TH1D* RebinnedOscDeltaNLSpectrumSumH[MatrixBins];
    TH1D* RebinnedOscDeltaVisibleSpectrumSumH[MatrixBins];
    TH1D* TotalOscillatedSpectrumAD[MaxDetectors];
    TH2D* EnergyMatrixH;
    TH2D* RowEnergyMatrixH;
    
    TH2D* NoNormalizedEnergyMatrixH;
    TH2D* NoNormalizedEnergyMatrixPosH;
    TH2D* NoNormalizedEnergyMatrixIAVH;
    TH2D* NoNormalizedEnergyMatrixNLH;
    TH2D* NoNormalizedEnergyMatrixResoH;
    
    //IAV:
    Double_t IAVNominalError; // relative uncertainty of the IAV thickness
    Double_t IAVError[MaxDetectors]; // relative uncertainty of the IAV thickness
    Double_t NominalIAVMatrix[MatrixBins][MatrixBins];
    Double_t IAVMatrix[MaxDetectors][MatrixBins][MatrixBins];
    Double_t NominalIAVMatrixFrac[MatrixBins];
    Double_t IAVMatrixFrac[MaxDetectors][MatrixBins];
    TH2F* IAVMatrixH[MaxDetectors];
    //Non linearity:
    // Detector response non-linearlity function
    TF1 * NLF;
    
    //Interpolation vectors
    Int_t sign[MatrixBins];
    Int_t EnergyIdx[MatrixBins];
    Double_t Energy[MatrixBins];
    Double_t dEtrue[MatrixBins];
    Double_t binScaling[MatrixBins];
    Double_t dNdE[MatrixBins];
    
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
    
    //Load functions
    void LoadIavCorrection();
    
    //Functions to calculate the energy matrix
    TF1* VisibleF;
    TF1* GetNeutrinoToPositronFunction(bool);
    
    void Interpolation(TF1*);
    
    Double_t VisibleEnergy0F(Double_t *, Double_t *);// Zeroth order true to visible function.
    Double_t VisibleEnergy1F(Double_t *, Double_t *);// First order true to visible function.
    
    Double_t NLBCWF(Double_t*, Double_t*);
    Double_t NLLBNLF(Double_t*, Double_t*);
    Double_t NLUnifiedF(Double_t*, Double_t*);
    
    Double_t ResolutionF(Double_t *, Double_t *);
    
    void LoadNLParameters();
    void SetNLParameters();
    
    void SetSystematic();
    
    Double_t EnergyVector[MatrixBins];
    Double_t AntineutrinoSpectrum[MatrixBins];
    
    Double_t m_positronTrueSpectrum[MatrixBins];
    Double_t m_positronIavDistortedSpectrum[MatrixBins];
    Double_t m_positronNLSpectrum[MatrixBins];
    Double_t m_positronDetectedSpectrum[MatrixBins];
    
public:
    CreateEnergyMatrix();
    CreateEnergyMatrix(NominalData*);
    
    void SetBCWModel(bool);
    void SetLBNLModel(bool);
    void SetUnifiedModel(bool);
    
    void  GenerateEnergyMatrix(Double_t, Double_t, Int_t);
    void  GenerateLBNLEnergyMatrix(Double_t, Double_t, Int_t);
    
    void  GetOscEnergyShift(Int_t,Int_t);
    void  GetOscIAVShift(Int_t,Int_t);
    void  GetOscNLShift(Int_t,Int_t);
    void  GetOscResolutionShift(Int_t,Int_t);
    void updatePositronTrue(Double_t,Double_t);
    void LoadADSpectrum();
};

//CreateEnergyMatrix :: CreateEnergyMatrix()
//{
//    Data = new NominalData();
//
//    analysis = Data->GetAnalysis();
//    BCW = Data->GetBCWModel();
//    LBNL = Data->GetLBNLModel();
//    Unified = Data->GetUnifiedModel();
//
//    IsotopeMatrix = Data->GetIsotopeMatrix();
//    ReactorPowerMatrix = Data->GetReactorPowerMatrix();
//    RelativeEnergyScaleMatrix = Data->GetRelativeEnergyScaleMatrix();
//    IAVMatrixb = Data->GetIAVMatrix();
//    NLMatrix = Data->GetNLMatrix();
//    ResolutionMatrix = Data->GetResolutionMatrix();
//
//    Nweeks = Data->GetWeeks();
//
//    InitialEnergy = Data->GetEmin();
//    FinalEnergy = Data->GetEmax();
//    InitialVisibleEnergy = Data->GetEVisMin();
//    FinalVisibleEnergy =  Data->GetEVisMax();
//
//    BinWidth=(FinalVisibleEnergy-InitialVisibleEnergy)/MatrixBins;
//
//    LinearBinning = Data->GetBinning();
//
//    //  Linear binning
//    if(LinearBinning)
//    {
//        n_evis_bins = Data->GetNbins();
//        n_etrue_bins = Data->GetNbins();
//
//        for (Int_t i = 0; i <= n_evis_bins; i++)
//        {
//            evis_bins[i] = 0.2 * i + 0.7;
//            enu_bins[i] = 0.2 * i + InitialEnergy;
//        }
//    }
//    //  Non-linear binning
//    else
//    {
//        n_evis_bins=37;
//        n_etrue_bins=39;
//
//        for (Int_t i = 0; i <= n_etrue_bins; i++)
//        {
//            enu_bins[i] = 0.2 * i + InitialEnergy;
//        }
//
//        evis_bins[0] = 0.7;
//        for (Int_t i = 0; i < n_evis_bins-1; i++)
//        {
//            evis_bins[i+1] = 0.2 * i + 1.0;
//        }
//        evis_bins[n_evis_bins] = FinalVisibleEnergy;
//    }
//
//    NADs = Data->GetADs();
//    ADsEH1 = 2;
//    ADsEH2 = 1;
//    ADsEH3 = 3;
//
//    if(NADs == 8)//    ADs can only be 6 or 8
//    {
//        ADsEH2 = 2;
//        ADsEH3 = 4;
//    }
//
//    //  IAV Error from Bryce
//    IAVNominalError=Data->GetIAVError();
//    //  Non uniformity
//    m_abs_escale = Data->GetAbsoluteEnergyScale();
//    m_abs_escale_nominal = m_abs_escale;
//    m_abs_escale_error = Data->GetAbsoluteEnergyScaleError();
//
//    ResolutionError = Data->GetResolutionError();
//    ResolutionErrorUncorrelated = Data->GetResoUncorrelatedError();
//    ResolutionRange = 8;// Why 8σ? Seems chosen trivially but in my opinion it's the way to limit the range of the convolution, for a Normal distribution this covers up to 99.99999... of the area
//
//    //  Scaling Backgrounds to the data rates:
//
//    for (Int_t AD = 0; AD<NADs; AD++)
//    {
//        for (Int_t week = 0; week<Nweeks; week++)
//        {
//            DetectorEfficiency[AD][week] = Data->GetDetectorEfficiency(AD,week);
//        }
//        //  Relative energy scale
//        m_rel_escale[AD] = Data->GetRelativeEnergyScale(AD);
//        m_rel_escale_error[AD] = Data->GetRelativeEnergyError(AD);
//        m_rel_escale_nominal[AD] = m_rel_escale[AD];
//        m_rel_eoffset[AD] = Data->GetRelativeEnergyOffset(AD);
//    }
//}

CreateEnergyMatrix :: CreateEnergyMatrix(NominalData* data)
{
    Data = new NominalData(0,2);
    Data->CopyData(data);

    analysis = Data->GetAnalysis();
   
    if(analysis)
    {
        AnalysisString = "Hydrogen";
    }
    else
    {
        AnalysisString = "Gadolinium";
    }
    BCW = Data->GetBCWModel();
    LBNL = Data->GetLBNLModel();
    Unified = Data->GetUnifiedModel();
    
    IsotopeMatrix = Data->GetIsotopeMatrix();
    ReactorPowerMatrix = Data->GetReactorPowerMatrix();
    RelativeEnergyScaleMatrix = Data->GetRelativeEnergyScaleMatrix();
    IAVMatrixb = Data->GetIAVMatrix();
    NLMatrix = Data->GetNLMatrix();
    ResolutionMatrix = Data->GetResolutionMatrix();
    
    Nweeks = Data->GetWeeks();
    
    InitialEnergy = Data->GetEmin();
    FinalEnergy = Data->GetEmax();
    InitialVisibleEnergy = Data->GetEVisMin();
    FinalVisibleEnergy = Data->GetEVisMax();
    
    BinWidth=(FinalVisibleEnergy-InitialVisibleEnergy)/MatrixBins;
    
    n_evis_bins = Data->GetVisibleBins();
    
    for (Int_t i = 0; i <= n_evis_bins; i++)
    {
        evis_bins[i] = Data->GetVisibleBinningArray(i);
    }
    
    n_etrue_bins = Data->GetTrueBins();
    
    for (Int_t i = 0; i <= n_etrue_bins; i++)
    {
        enu_bins[i] = Data->GetTrueBinningArray(i);
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
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            DetectorEfficiency[AD][week] = Data->GetDetectorEfficiency(AD,week);
            
        }
        //  Relative energy scale
        m_rel_escale[AD] = Data->GetRelativeEnergyScale(AD);
        m_rel_escale_error[AD] = Data->GetRelativeEnergyError(AD);
        m_rel_escale_nominal[AD] = m_rel_escale[AD];
        m_rel_eoffset[AD] = Data->GetRelativeEnergyOffset(AD);
    }
}

void CreateEnergyMatrix :: GenerateEnergyMatrix(Double_t sin22t13, Double_t dm2_31, Int_t week)
{
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                           Set kinematic energy shift function
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    GetNeutrinoToPositronFunction(1);//0 for 0th order, 1 for 1st order
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                           Set nominal IAV Matrix
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    LoadIavCorrection();
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                            Set nominal NL parameters
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    this->SetBCWModel(BCW);
    this->SetLBNLModel(LBNL);
    this->SetUnifiedModel(Unified);
    
    LoadNLParameters();//nominal parameters
    SetNLParameters();
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                              Set nominal resolution parameters
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ResoF = new TF1("ResoF",this,&CreateEnergyMatrix::ResolutionF,0,20,3,"CreateEnergyMatrix","ResolutionF");
    //    ResoF->SetParameters(0.022,0.077,0.018); // based on Bryce's TN
    ResoF->SetParameters(0.0148,0.0869,0.0271); // based on Doc-8982
    
    //BCW values http://dayabay.ihep.ac.cn/DocDB/0087/008768/013/6AdAnalysis-BCW.pdf are different! 0.13226034931, 0.32604048828, 0.26196488314
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //                                          Calculate distortions
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    LoadADSpectrum();
    
    for(Int_t TrueEnergyIndex=0; TrueEnergyIndex<MatrixBins; TrueEnergyIndex++)
    {
        GetOscEnergyShift(TrueEnergyIndex,week);
        GetOscIAVShift(TrueEnergyIndex,week);
        GetOscNLShift(TrueEnergyIndex,week);
        GetOscResolutionShift(TrueEnergyIndex,week);
    }
    
    
    TH2D* EnergyMatrixPosH;
    TH2D* EnergyMatrixIAVH;
    TH2D* EnergyMatrixNLH;
    TH2D* EnergyMatrixResoH;
    
    std::cout << "\t Generating Energy Matrix" << std::endl;
    
    Char_t EnergyMatrixC[50];
    
    sprintf(EnergyMatrixC,("./ResponseMatrices/"+AnalysisString+"/NominalResponseMatrix.root").c_str());
    
    TFile* EnergyMatrixDataF = TFile::Open(EnergyMatrixC,"recreate");
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                          Matrices with the original 240x240 binning
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    if(RebinEnergyMatrix)
    {
        std::cout << "\t Rebin Matrix from Fine to Coarse binning" << std::endl;
        
        LimitVis = n_evis_bins;
        LimitTrue = n_etrue_bins;
        
        EnergyMatrixH = new TH2D("EvisEnu","EvisEnu",LimitTrue,enu_bins,LimitVis,evis_bins);
        EnergyMatrixPosH = new TH2D("EvisEnuPos","EvisEnuPos",LimitTrue,enu_bins,LimitVis,evis_bins);
        EnergyMatrixIAVH = new TH2D("EvisEnuIAV","EvisEnuIAV",LimitTrue,enu_bins,LimitVis,evis_bins);
        EnergyMatrixNLH = new TH2D("EvisEnuNL","EvisEnuNL",LimitTrue,enu_bins,LimitVis,evis_bins);
        EnergyMatrixResoH = new TH2D("EvisEnuReso","EvisEnuReso",LimitTrue,enu_bins,LimitVis,evis_bins);
        RowEnergyMatrixH = new TH2D("EnuEvis","EnuEvis", LimitTrue,enu_bins,LimitVis,evis_bins);
        
        for (Int_t i = 0; i < n_etrue_bins; i++)
        {
            OscDeltaPositronSpectrumSumH[i] = new TH1D(Form("Fine Positron Spectrum Vis for True Energy Index %d",i),Form("Fine Positron Spectrum Vis for True Energy Index %d",i),MatrixBins,0,FinalVisibleEnergy);
            OscDeltaIAVSpectrumSumH[i]= new TH1D(Form("Fine IAV Spectrum Vis for True Energy Index %d",i),Form("Fine Spectrum IAV Vis for True Energy Index %d",i),MatrixBins,0,FinalVisibleEnergy);
            OscDeltaNLSpectrumSumH[i]= new TH1D(Form("Fine NL Spectrum Vis for True Energy Index %d",i),Form("Fine Spectrum NL Vis for True Energy Index %d",i),MatrixBins,0,FinalVisibleEnergy);
            
            OscDeltaVisibleSpectrumSumH[i]= new TH1D(Form("Fine Visible Spectrum Vis for True Energy Index %d",i),Form("Fine Visible Spectrum Vis for True Energy Index %d",i),MatrixBins,0,FinalVisibleEnergy);
            
            for (Int_t TrueEnergyIndex = Int_t(enu_bins[i]*MatrixBins/FinalVisibleEnergy); TrueEnergyIndex < Int_t(enu_bins[i+1]*MatrixBins/FinalVisibleEnergy); TrueEnergyIndex++)//Although the results are integers I need to add Int_t() if I don't want the precision errors to mess up the rebinning.
            {
                OscDeltaPositronSpectrumSumH[i]->Add(OscDeltaPositronSpectrumH[TrueEnergyIndex]);
                OscDeltaIAVSpectrumSumH[i]->Add(OscDeltaIAVSpectrumH[TrueEnergyIndex]);
                OscDeltaNLSpectrumSumH[i]->Add(OscDeltaNLSpectrumH[TrueEnergyIndex]);
                OscDeltaVisibleSpectrumSumH[i]->Add(OscDeltaVisibleSpectrumH[TrueEnergyIndex]);
            }
            
            RebinnedOscDeltaPositronSpectrumSumH[i]=(TH1D*)OscDeltaPositronSpectrumSumH[i]->Rebin(n_evis_bins,Form("Rebinned Positron Spectrum Vis for True Energy Index %d",i),evis_bins);
            RebinnedOscDeltaIAVSpectrumSumH[i]=(TH1D*)OscDeltaIAVSpectrumSumH[i]->Rebin(n_evis_bins,Form("Rebinned IAV Spectrum Vis for True Energy Index %d",i),evis_bins);
            RebinnedOscDeltaNLSpectrumSumH[i]=(TH1D*)OscDeltaNLSpectrumSumH[i]->Rebin(n_evis_bins,Form("Rebinned NL Spectrum Vis for True Energy Index %d",i),evis_bins);
            RebinnedOscDeltaVisibleSpectrumSumH[i]=(TH1D*)OscDeltaVisibleSpectrumSumH[i]->Rebin(n_evis_bins,Form("Rebinned Visible Spectrum Vis for True Energy Index %d",i),evis_bins);
            
            for (Int_t j = 0; j < n_evis_bins; j++)
            {
                EnergyMatrixH->SetBinContent(i+1,j+1, RebinnedOscDeltaVisibleSpectrumSumH[i]->GetBinContent(j+1));
                EnergyMatrixPosH->SetBinContent(i+1,j+1,RebinnedOscDeltaPositronSpectrumSumH[i]->GetBinContent(j+1));
                EnergyMatrixIAVH->SetBinContent(i+1,j+1,RebinnedOscDeltaIAVSpectrumSumH[i]->GetBinContent(j+1));
                EnergyMatrixNLH->SetBinContent(i+1,j+1,RebinnedOscDeltaNLSpectrumSumH[i]->GetBinContent(j+1));
                EnergyMatrixResoH->SetBinContent(i+1,j+1,RebinnedOscDeltaVisibleSpectrumSumH[i]->GetBinContent(j+1));

            }
        }
        for (Int_t i = 0; i < n_etrue_bins; i++)
        {
            for (Int_t j = 0; j < n_evis_bins; j++)
            {
                RowEnergyMatrixH->SetBinContent(i+1,j+1,EnergyMatrixH->GetBinContent(i+1,j+1));
                
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //                          Maybe instead of transposing the correct option is to invert it. Study this!!
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            }
        }
    }
    else
    {
        LimitVis = MatrixBins;
        LimitTrue = MatrixBins;
        
        EnergyMatrixH= new TH2D("EvisEnu","EvisEnu",LimitTrue,0,LimitTrue,LimitVis,0,LimitVis);
        EnergyMatrixPosH= new TH2D("EvisEnuPos","EvisEnuPos",LimitTrue,0,LimitTrue,LimitVis,0,LimitVis);
        EnergyMatrixIAVH= new TH2D("EvisEnuIAV","EvisEnuIAV",LimitTrue,0,LimitTrue,LimitVis,0,LimitVis);
        EnergyMatrixNLH= new TH2D("EvisEnuNL","EvisEnuNL",LimitTrue,0,LimitTrue,LimitVis,0,LimitVis);
        EnergyMatrixResoH= new TH2D("EvisEnuReso","EvisEnuReso",LimitTrue,0,LimitTrue,LimitVis,0,LimitVis);
        RowEnergyMatrixH= new TH2D("EnuEvis","EnuEvis",LimitVis,0,LimitVis, LimitTrue,0,LimitTrue);
        
        for (Int_t i = 0; i < MatrixBins; i++)
        {
            for (Int_t j = 0; j < MatrixBins; j++)
            {
                EnergyMatrixH->SetBinContent(i+1,j+1,EnergyMatrixH->GetBinContent(i+1,j+1)+OscDeltaVisibleSpectrumH[i]->GetBinContent(j+1));
                
                EnergyMatrixPosH->SetBinContent(i+1,j+1,EnergyMatrixPosH->GetBinContent(i+1,j+1)+OscDeltaPositronSpectrumH[i]->GetBinContent(j+1));
                
                EnergyMatrixIAVH->SetBinContent(i+1,j+1,EnergyMatrixIAVH->GetBinContent(i+1,j+1)+OscDeltaIAVSpectrumH[i]->GetBinContent(j+1));
                
                EnergyMatrixNLH->SetBinContent(i+1,j+1,EnergyMatrixNLH->GetBinContent(i+1,j+1)+OscDeltaNLSpectrumH[i]->GetBinContent(j+1));
                
                EnergyMatrixResoH->SetBinContent(i+1,j+1,EnergyMatrixResoH->GetBinContent(i+1,j+1)+OscDeltaVisibleSpectrumH[i]->GetBinContent(j+1));
            }
        }
        for (Int_t i = 0; i < MatrixBins; i++)
        {
            for (Int_t j = 0; j < MatrixBins; j++)
            {
                RowEnergyMatrixH->SetBinContent(j+1,i+1,EnergyMatrixH->GetBinContent(i+1,j+1));
            }
        }
        
    }
    
    // UNCOMMENT FOLLOWING LINES TO SAVE ENERGY MATRIX SLICES IN EACH PRODUCTION STEP
    //    for(Int_t TrueEnergyIndex = 0; TrueEnergyIndex<MatrixBins; TrueEnergyIndex++)
    //    {
    //        OscDeltaPositronSpectrumH[TrueEnergyIndex]->Write();
    //        OscDeltaIAVSpectrumH[TrueEnergyIndex]->Write();
    //        OscDeltaNLSpectrumH[TrueEnergyIndex]->Write();
    //        OscDeltaVisibleSpectrumH[TrueEnergyIndex]->Write();
    //    }
    if(RebinEnergyMatrix)
    {
        for(Int_t TrueEnergyIndex = 0; TrueEnergyIndex<n_etrue_bins; TrueEnergyIndex++)
        {
            //        OscDeltaPositronSpectrumSumH[TrueEnergyIndex]->Write();
            //        OscDeltaIAVSpectrumSumH[TrueEnergyIndex]->Write();
            //        OscDeltaNLSpectrumSumH[TrueEnergyIndex]->Write();
            //        OscDeltaVisibleSpectrumSumH[TrueEnergyIndex]->Write();
            
            delete OscDeltaPositronSpectrumSumH[TrueEnergyIndex];
            delete OscDeltaIAVSpectrumSumH[TrueEnergyIndex];
            delete OscDeltaNLSpectrumSumH[TrueEnergyIndex];
            delete OscDeltaVisibleSpectrumSumH[TrueEnergyIndex];
            
            delete RebinnedOscDeltaPositronSpectrumSumH[TrueEnergyIndex];
            delete RebinnedOscDeltaIAVSpectrumSumH[TrueEnergyIndex];
            delete RebinnedOscDeltaNLSpectrumSumH[TrueEnergyIndex];
            delete RebinnedOscDeltaVisibleSpectrumSumH[TrueEnergyIndex];
        }
    }
    //
    NoNormalizedEnergyMatrixH = (TH2D*)EnergyMatrixH->Clone("NoNormalizedEvisEnu");
    NoNormalizedEnergyMatrixPosH = (TH2D*)EnergyMatrixPosH->Clone("PosNoNormalizedEvisEnu");
    NoNormalizedEnergyMatrixIAVH = (TH2D*)EnergyMatrixIAVH->Clone("IAVNoNormalizedEvisEnu");
    NoNormalizedEnergyMatrixNLH = (TH2D*)EnergyMatrixNLH->Clone("NLNoNormalizedEvisEnu");
    NoNormalizedEnergyMatrixResoH = (TH2D*)EnergyMatrixResoH->Clone("ResoNoNormalizedEvisEnu");
    
    NoNormalizedEnergyMatrixPosH->Write();
    NoNormalizedEnergyMatrixIAVH->Write();
    NoNormalizedEnergyMatrixNLH->Write();
    NoNormalizedEnergyMatrixResoH->Write();
    
    NoNormalizedEnergyMatrixH->Write();//Save Matrix before rebinning
    RowEnergyMatrixH->Write("NoNormalizedEnuEvis");//Save Matrix before rebinning
    
    //Normalize the Matrix
    for(Int_t j=0;j<LimitTrue;j++)
    {
        Norma[j]=0;
        NormaPos[j]=0;
        NormaIAV[j]=0;
        NormaNL[j]=0;
        NormaReso[j]=0;
        for(Int_t i=0;i<LimitVis;i++)
        {
            Norma[j] = Norma[j]+EnergyMatrixH->GetBinContent(j+1,i+1);
            NormaPos[j] = NormaPos[j]+EnergyMatrixPosH->GetBinContent(j+1,i+1);
            NormaIAV[j] = NormaIAV[j]+EnergyMatrixIAVH->GetBinContent(j+1,i+1);
            NormaNL[j] = NormaNL[j]+EnergyMatrixNLH->GetBinContent(j+1,i+1);
            NormaReso[j] = NormaReso[j]+EnergyMatrixResoH->GetBinContent(j+1,i+1);
        }
    }
    
    for(Int_t i=0;i<LimitVis;i++)
    {
        NormaTrans[i]=0;
        for(Int_t j=0;j<LimitTrue;j++)
        {
            NormaTrans[i] = NormaTrans[i]+RowEnergyMatrixH->GetBinContent(j+1,i+1);
        }
    }
    
    for (Int_t i = 0; i < LimitVis; i++)
    {
        for (Int_t j = 0; j < LimitTrue; j++)
        {
            if(Norma[j]!=0)
            {
                EnergyMatrixH->SetBinContent(j+1,i+1,EnergyMatrixH->GetBinContent(j+1,i+1)/Norma[j]);//Normalization so Σj E(i,j) = 1; (Σ(y axis) =1)
            }
            if(NormaPos[j]!=0)
            {
                EnergyMatrixPosH->SetBinContent(j+1,i+1,EnergyMatrixPosH->GetBinContent(j+1,i+1)/NormaPos[j]);
            }
            if(NormaIAV[j]!=0)
            {
                EnergyMatrixIAVH->SetBinContent(j+1,i+1,EnergyMatrixIAVH->GetBinContent(j+1,i+1)/NormaIAV[j]);
            }
            if(NormaNL[j]!=0)
            {
                EnergyMatrixNLH->SetBinContent(j+1,i+1,EnergyMatrixNLH->GetBinContent(j+1,i+1)/NormaNL[j]);
            }
            if(NormaReso[j]!=0)
            {
                EnergyMatrixResoH->SetBinContent(j+1,i+1,EnergyMatrixResoH->GetBinContent(j+1,i+1)/NormaReso[j]);
            }
            if(NormaTrans[i]!=0)
            {
                RowEnergyMatrixH->SetBinContent(j+1,i+1,RowEnergyMatrixH->GetBinContent(j+1,i+1)/NormaTrans[i]);
            }
        }
    }
    
    //Save Matrix after rebinning
    EnergyMatrixH->Write();
    EnergyMatrixPosH->Write();
    EnergyMatrixIAVH->Write();
    EnergyMatrixNLH->Write();
    EnergyMatrixResoH->Write();
    RowEnergyMatrixH->Write();
    
    EnergyMatrixDataF->Close();
    if(Print)
    {
        TCanvas* EnergyC = new TCanvas("EnergyC","EnergyC");
        EnergyC->SetLogz();
        EnergyMatrixH->SetStats(0);
        EnergyMatrixH->SetTitle("E_{vis} - E_{#nu}");
        EnergyMatrixH->Draw("colz");
        
        EnergyC->Print(("./Images/"+AnalysisString+"/Detector/ResponseMatrix.eps").c_str(),".eps");
        
        delete EnergyC;
    }
    delete NoNormalizedEnergyMatrixPosH;
    delete NoNormalizedEnergyMatrixIAVH;
    delete NoNormalizedEnergyMatrixNLH;
    delete NoNormalizedEnergyMatrixResoH;
    delete NoNormalizedEnergyMatrixH;
    
    delete EnergyMatrixIAVH;
    delete EnergyMatrixNLH;
    delete EnergyMatrixPosH;
    delete EnergyMatrixResoH;
    delete ResoF;
    std::cout << "\t Finished Generating Energy Matrix" << std::endl;
}

void CreateEnergyMatrix :: SetBCWModel(bool bcw)
{
    if(bcw)
    {
        NLModelE=BCWE;
        NLF = new TF1("NLF",this,&CreateEnergyMatrix::NLBCWF,InitialVisibleEnergy,FinalVisibleEnergy,7,"CreateEnergyMatrix","NLBCWF");
    }
}

void CreateEnergyMatrix :: SetLBNLModel(bool lbnl)
{
    if(lbnl)
    {
        NLModelE=LBNLE;
        NLF = new TF1("NLF",this,&CreateEnergyMatrix::NLLBNLF,InitialVisibleEnergy,FinalVisibleEnergy,5,"CreateEnergyMatrix","NLLBNLF");
    }
}

void CreateEnergyMatrix :: SetUnifiedModel(bool unified)
{
    if(unified)
    {
        NLModelE=UnifiedE;
        NLF = new TF1("NLF",this,&CreateEnergyMatrix::NLUnifiedF,InitialVisibleEnergy,FinalVisibleEnergy,m_num_unified_nl_pars+2,"CreateEnergyMatrix","NLUnifiedF");
    }
}

TF1* CreateEnergyMatrix :: GetNeutrinoToPositronFunction(bool order)
{
    if(order==0)
    {
        Double_t Correction = Mn-Mp-Me;
        VisibleF = new TF1("VisibleF",this,&CreateEnergyMatrix::VisibleEnergy0F,InitialEnergy,FinalVisibleEnergy,1,"CreateEnergyMatrix","VisibleEnergy0F");
        VisibleF->SetParameter(0,Correction);
        std::cout<<"Zeroth order"<<std::endl;
    }
    if(order==1)
    {
        VisibleF = new TF1("VisibleF",this,&CreateEnergyMatrix::VisibleEnergy1F,InitialEnergy,FinalVisibleEnergy,1,"CreateEnergyMatrix","VisibleEnergy1F");
        VisibleF->SetParameter(0,-2);//Maybe I can improve this using a MC simulation of the angular distribution calculated in http://authors.library.caltech.edu/2796/1/VOGprd99.pdf
        std::cout<<"First order"<<std::endl;
        
    }
    
    for(Int_t i = 0; i < MatrixBins; i++)
    {
        EnergyVector[i] = (0.5+i)*BinWidth;//set energy vector for later interpolation purposes
    }
    
    return VisibleF;
}


void CreateEnergyMatrix :: Interpolation(TF1* func)
{
    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<MatrixBins; VisibleEnergyIndex++)
    {
        Energy[VisibleEnergyIndex]=(func->GetX(EnergyVector[VisibleEnergyIndex]));//GetX
        EnergyIdx[VisibleEnergyIndex]=((Int_t)(Energy[VisibleEnergyIndex]/BinWidth));
        
        //        std::cout << "IAV Energy " << IAVEnergy[VisibleEnergyIndex] << std::endl;
        //        std::cout << "IAV Energy Idx " << IAVEnergyIdx[VisibleEnergyIndex] << std::endl;
        //        std::cout << " VisibleEnergyIndex" << VisibleEnergyIndex << std::endl;
        
        binScaling[VisibleEnergyIndex]=(func->Derivative(Energy[VisibleEnergyIndex]));
        
        //        std::cout << "BIN SCALING" << binScaling[VisibleEnergyIndex] << std::endl;
        
        sign[VisibleEnergyIndex]=(1);
        
        dEtrue[VisibleEnergyIndex]=(Energy[VisibleEnergyIndex] - (0.5+EnergyIdx[VisibleEnergyIndex])*BinWidth);
        
        //        std::cout << "dEtrue" << dEtrue[VisibleEnergyIndex] << std::endl;
        if (dEtrue[VisibleEnergyIndex] < 0)
        {
            sign[VisibleEnergyIndex]=(-1);
            
        }
        //        std::cout << "sign" << sign[VisibleEnergyIndex] << std::endl;
    }
}

void CreateEnergyMatrix :: GetOscEnergyShift(Int_t TrueEnergyIndex, Int_t week)
{
    
    OscDeltaPositronSpectrumH[TrueEnergyIndex] = new TH1D(Form("PositronTrueSpectrum,W%d%d",week,TrueEnergyIndex),Form("PositronTrueSpectrum,W%d%d",week,TrueEnergyIndex), MatrixBins,InitialVisibleEnergy,FinalVisibleEnergy);

//    PositronEnergy = TrueEnergyIndex*BinWidth;
//    
//    Double_t e_nu = VisibleF->GetX(PositronEnergy);
//    
//    Int_t e_nu_Idx = (Int_t)((e_nu)/BinWidth);
    
    
//    if (PositronEnergy < 1.022)
//    {
//        OscDeltaPositronSpectrumH[TrueEnergyIndex]->SetBinContent(TrueEnergyIndex+1, 0);
//        return;
//    }
//    
//    if(e_nu_Idx>=MatrixBins)
//    {
//        OscDeltaPositronSpectrumH[TrueEnergyIndex]->SetBinContent(TrueEnergyIndex+1, 0);
//        return;
//    }
//    
//    Double_t binScaling = VisibleF->Derivative(e_nu,0,0.0001);
//    Double_t dNdE = 0;
//    Int_t sign = 1;
//    Double_t dE = e_nu - (0.5 + e_nu_Idx)*BinWidth;
//    
//    if (dE < 0)
//    {
//        sign = -1;
//    }
    
//    if(e_nu_Idx==MatrixBins-1)
//    {
//        dNdE = 1;//flat
//    }
//    else
//    {
//        dNdE =(TMath::Abs(dE)/BinWidth)+(1 - TMath::Abs(dE)/BinWidth);
//    }
//    
    //Flat:
//    if(FlatEnergyMatrix)
//    {
//        OscDeltaPositronSpectrumH[TrueEnergyIndex]->SetBinContent(TrueEnergyIndex-Int_t(1.022/BinWidth)+1,dNdE/binScaling);
//    }
    //    else
    //    {
    //    Limit has to be 1.8 if oscillation is used
            //Interpolate Visible Function
//            Double_t binScaling = VisibleF->Derivative(e_nu,0,0.0001);
//            Double_t dNdE = 0;
//            Int_t sign = 1;
//            Double_t dE = e_nu - (0.5 + e_nu_Idx)*BinWidth;
//    
//            if (dE < 0)
//            {
//                sign = -1;
//            }
//    
//            if(e_nu_Idx==MatrixBins-1)
//            {
//                dNdE = OriginalPredictionH[far][near]->GetBinContent(e_nu_Idx-Int_t(Limit/BinWidth)+1);
//            }
//            else
//            {
//                dNdE =
//                (TMath::Abs(dE)/BinWidth)
//                * OriginalPredictionH[far][near]->GetBinContent(e_nu_Idx-Int_t(Limit/BinWidth)+sign+1)
//                +(1 - TMath::Abs(dE)/BinWidth)
//                * OriginalPredictionH[far][near]->GetBinContent(e_nu_Idx-Int_t(Limit/BinWidth)+1);
//            }
//    
//           OscDeltaPositronSpectrumH[TrueEnergyIndex]->SetBinContent(TrueEnergyIndex-Int_t(1.022/BinWidth)+1,dNdE/binScaling);
//        }
    
    
    Double_t Limit = 1.8;
    
    PositronEnergy = EnergyVector[TrueEnergyIndex];

    Double_t e_nu = VisibleF->GetX(PositronEnergy);
    
    Int_t e_nu_Idx = (Int_t)((e_nu)/BinWidth);//TRUE BINWIDTH
    
    if(e_nu_Idx>=MatrixBins)
    {
        OscDeltaPositronSpectrumH[TrueEnergyIndex]->SetBinContent(TrueEnergyIndex+1, 0);
        return;
    }

    if (PositronEnergy < 1.022)
    {
        OscDeltaPositronSpectrumH[TrueEnergyIndex]->SetBinContent(TrueEnergyIndex+1, 0);
        return;
    }
    
    Double_t binScaling = VisibleF->Derivative(e_nu,0,0.0001);
    Double_t dNdE = 0;
    Int_t sign = 1;
    Double_t dE = e_nu - (0.5 + e_nu_Idx)*BinWidth;
    
    if (dE < 0)
    {
        sign = -1;
    }
    
    if(e_nu_Idx==MatrixBins-1)
    {
        dNdE = TotalOscillatedSpectrumAD[0]->GetBinContent(e_nu_Idx-Int_t(Limit/BinWidth)+1);//substract 1.8 in the index because the histogram starts in 1.8 MeV instead of 0 MeV.
    }
    else
    {
        dNdE = (TMath::Abs(dE)/BinWidth) * TotalOscillatedSpectrumAD[0]->GetBinContent(e_nu_Idx-Int_t(Limit/BinWidth)+sign+1)
        +(1 - TMath::Abs(dE)/BinWidth) * TotalOscillatedSpectrumAD[0]->GetBinContent(e_nu_Idx-Int_t(Limit/BinWidth)+1);
    }
    
    OscDeltaPositronSpectrumH[TrueEnergyIndex]->SetBinContent(TrueEnergyIndex+1,dNdE/binScaling);
}

void CreateEnergyMatrix :: GetOscIAVShift(Int_t TrueEnergyIndex, Int_t week)
{
    //IAV
    OscDeltaIAVSpectrumH[TrueEnergyIndex] = new TH1D(Form("IAVSpectrumW%d%d",week,TrueEnergyIndex),Form("IAVSpectrum,W%d%d",week,TrueEnergyIndex), MatrixBins,InitialVisibleEnergy,FinalVisibleEnergy);
    
    //Calculate IAV Shift
    if (PositronEnergy >= 1.022)//2*0.511
    {
        for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TrueEnergyIndex+1; VisibleEnergyIndex++)
        {
            OscDeltaIAVSpectrumH[TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, OscDeltaIAVSpectrumH[TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1) + IAVMatrix[0][TrueEnergyIndex-Int_t(1.022/BinWidth)][VisibleEnergyIndex] * OscDeltaPositronSpectrumH[TrueEnergyIndex]->GetBinContent(TrueEnergyIndex+1));
        }
    }
}

void CreateEnergyMatrix :: GetOscNLShift(Int_t TrueEnergyIndex, Int_t week)
{
    //NL
    OscDeltaNLSpectrumH[TrueEnergyIndex] = new TH1D(Form("NLSpectrum,W%d%d",week,TrueEnergyIndex),Form("NLSpectrum,W%d%d",week,TrueEnergyIndex), MatrixBins,InitialVisibleEnergy,FinalVisibleEnergy);
    
    //Calculate NL Shift
    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<MatrixBins; VisibleEnergyIndex++)
    {
        if(EnergyIdx[VisibleEnergyIndex]==(MatrixBins-1))
        {
            dNdE[VisibleEnergyIndex]=(OscDeltaIAVSpectrumH[TrueEnergyIndex]->GetBinContent(EnergyIdx[VisibleEnergyIndex]+1));
        }
        else
        {
            dNdE[VisibleEnergyIndex] = ((TMath::Abs(dEtrue[VisibleEnergyIndex])/BinWidth)*OscDeltaIAVSpectrumH[TrueEnergyIndex]->GetBinContent(EnergyIdx[VisibleEnergyIndex]+sign[VisibleEnergyIndex]+1)+(1 - TMath::Abs(dEtrue[VisibleEnergyIndex])/BinWidth)* OscDeltaIAVSpectrumH[TrueEnergyIndex]->GetBinContent(EnergyIdx[VisibleEnergyIndex]+1));
        }
    }
    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<MatrixBins; VisibleEnergyIndex++)
    {
        OscDeltaNLSpectrumH[TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, dNdE[VisibleEnergyIndex]/binScaling[VisibleEnergyIndex]);
    }
}

void CreateEnergyMatrix :: GetOscResolutionShift(Int_t TrueEnergyIndex, Int_t week)
{
    //Resolution
    
    
    OscDeltaVisibleSpectrumH[TrueEnergyIndex] = new TH1D(Form("VisibleSpectrum,W%d%d",week,TrueEnergyIndex),Form("VisibleSpectrum,W%d%d",week,TrueEnergyIndex), MatrixBins,InitialVisibleEnergy,FinalVisibleEnergy);
    
    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<MatrixBins; VisibleEnergyIndex++)
    {
        OscDeltaVisibleSpectrumH[TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, OscDeltaNLSpectrumH[TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1));//Copy histogram, I don't use clone() because of memory leaks.
        
        if(OscDeltaNLSpectrumH[TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1)==0)
        {
            continue;
        }
    }
    
    //Calculate Resolution effect
    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<MatrixBins; VisibleEnergyIndex++)
    {
        Double_t sigma;
        
        sigma = (ResoF->Eval(VisibleEnergyIndex*BinWidth)) * VisibleEnergyIndex*BinWidth;
        
        Double_t minDetE = VisibleEnergyIndex*BinWidth - ResolutionRange*sigma;
        Double_t maxDetE = VisibleEnergyIndex*BinWidth + ResolutionRange*sigma;
        Int_t minDetEIdx = (Int_t)((minDetE)/BinWidth);
        Int_t maxDetEIdx = (Int_t)((maxDetE)/BinWidth);
        
        if(minDetEIdx < 0){minDetEIdx = 0;}
        if(maxDetEIdx >= MatrixBins){maxDetEIdx = MatrixBins-1;}
        
        for(Int_t detIdx=minDetEIdx; detIdx<=maxDetEIdx; detIdx++)
        {
            if(detIdx==0)
            {
                continue;
            }
            Double_t gausFactor = TMath::Gaus((VisibleEnergyIndex-detIdx)*BinWidth,0,sigma,true);
            
            OscDeltaVisibleSpectrumH[TrueEnergyIndex]->SetBinContent(detIdx+1, (OscDeltaVisibleSpectrumH[TrueEnergyIndex]->GetBinContent(detIdx+1)+OscDeltaNLSpectrumH[TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1)*gausFactor));
        }
    }
    // Scale reso to original NL energy.
    if ((OscDeltaVisibleSpectrumH[TrueEnergyIndex]->Integral())!=0)
    {
        OscDeltaVisibleSpectrumH[TrueEnergyIndex]->Scale(OscDeltaNLSpectrumH[TrueEnergyIndex]->Integral()/OscDeltaVisibleSpectrumH[TrueEnergyIndex]->Integral());
    }
}
void CreateEnergyMatrix :: LoadNLParameters()
{
    switch (NLModelE)
    {
        case 0://BCW NL Model
        {
            std::cout << "LOAD BCW NL PARAMETERS" << std::endl;
            ifstream bcw_positron_data("./Inputs/bcw_nl_data/positron.dat");
            for(Int_t i = 0; i < n_bcw_positron_nl; i++)
            {
                bcw_positron_data >> m_bcw_positron_nl_e[i] >> m_bcw_positron_nl_fac[i];
            }
            bcw_positron_data.close();
            
            ifstream bcw_elec_data("./Inputs/bcw_nl_data/par.dat");
            for(Int_t i = 0; i < 3; i++)
            {
                bcw_elec_data >> m_bcw_elec_nl_par_nominal[i]  >> m_bcw_elec_nl_par_error[i];
                m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];
                
                //                std::cout <<  m_bcw_elec_nl_par_error[i] << "IS ELEC NL PAR ERROR" << std::endl;
                //                std::cout <<  m_bcw_elec_nl_par[i] << "IS ELEC NL PAR" << std::endl;
                
            }
            bcw_elec_data.close();
            
            for(Int_t i = 3; i < 5; i++)
            { // those describe additional uncertainty for bcw model that are described in Doc-XXXX
                m_bcw_elec_nl_par[i] = 0;
                m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];
                m_bcw_elec_nl_par_error[i] = 1;
            }
            
            TFile *bcw_ele_err_file = TFile::Open("./Inputs/bcw_nl_data/ele_err.root");
            g_bcw_elec_nl_error[0] = (TGraph*)bcw_ele_err_file->Get("g_up")->Clone();
            g_bcw_elec_nl_error[1] = (TGraph*)bcw_ele_err_file->Get("g_down")->Clone();
            bcw_ele_err_file->Close();
        }
            break;
        case 1://LBNL NL Model
        {
            std::cout << "LOAD LBNL NL PARAMETERS" << std::endl;
            ifstream lbnl_positron_data("./Inputs/lbnl_nl_data/lbnl_positron_nl.txt");
            if(!lbnl_positron_data.is_open())
            {
                std::cout << "Error: cannot find LBNL non-linearity curve!!!" << std::endl;
                exit(0);
            }
            for(Int_t i = 0; i < 300; i++)
            {
                lbnl_positron_data >> m_lbnl_positron_nl_e[i] >> m_lbnl_positron_nl_fac[i] >> total_err >> m_lbnl_positron_nl_err[0][i] >> m_lbnl_positron_nl_err[1][i] >> m_lbnl_positron_nl_err[2][i];
            }
            for(Int_t i = 0; i < 3; i++)
            {
                m_lbnl_nl_par[i] = 0;
                m_lbnl_nl_par_nominal[i] = 0;
                m_lbnl_nl_par_error[i] = 1.0;
            }
        }
            break;
        case 2://Unified NL Model
        {
            for(Int_t i = 0; i < m_num_unified_nl_pars; i++)
            {
                m_unified_nl_par[i] = 0.0;
                m_unified_nl_par_nominal[i] = 0.0;
                m_unified_nl_par_error[i] = 1.0;
            }
            
            TFile *unified_nl_file = TFile::Open("./Inputs/unified_nl_data/nl_models_final.root");
            
            // Use IHEP I as the nominal model
            g_unified_positron_nl = (TGraph*)unified_nl_file->Get("positron_0")->Clone();
            g_unified_positron_nl_pulls[0] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",1))->Clone();
            g_unified_positron_nl_pulls[1] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",2))->Clone();
            g_unified_positron_nl_pulls[2] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",3))->Clone();
            g_unified_positron_nl_pulls[3] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",4))->Clone();
            
            unified_nl_file->Close();
            
            // Copy into arrays to speed up
            for(Int_t ie = 0; ie < n_unified_nl_points; ie++)
            {
                Double_t e = 1.022 + 0.02 * ie;
                m_unified_positron_nl_e[ie] = e;
                m_unified_positron_nl_fac[ie] = g_unified_positron_nl->Eval(e);
                for (Int_t i = 0; i < m_num_unified_nl_pars; i++)
                {
                    m_unified_positron_nl_err[i][ie] = g_unified_positron_nl_pulls[i]->Eval(e) - m_unified_positron_nl_fac[ie];
                }
            }
            
            TMatrixD covmatrix_unified(4,4);
            
            //            //Read covariance matrix
            bool m_correlateUnifiedNonlinearPars = 0; //Fix this when I have the covariance matrix.
            if (m_correlateUnifiedNonlinearPars)
            {
                TFile* f_covmatrix = new TFile("NonlinearCovmatrixFilename");
                
                f_covmatrix->ls();
                Double_t * mat_tmp = ((TMatrixD*)f_covmatrix->Get("coeffmatrix"))->GetMatrixArray();
                for (Int_t i = 0; i < 4; i++)
                {
                    for (Int_t j = 0; j < 4; j++)
                    {
                        m_unified_nl_par_covmatrix[i][j] =  mat_tmp[i*4+j];
                    }
                }
            }
            else
            {
                for(Int_t i = 0; i < 4; i++)
                {
                    for(Int_t j = 0; j < 4; j++)
                    {
                        if (i == j)
                        {
                            m_unified_nl_par_covmatrix[i][j] = 1;
                        }
                        else
                        {
                            m_unified_nl_par_covmatrix[i][j] = 0;
                        }
                    }
                }
            }
            
            covmatrix_unified.SetMatrixArray(&m_unified_nl_par_covmatrix[0][0]);
            //            std::cout << "Covariance matrix for UNIFIED non-liner parameters:" << std::endl;
            //            covmatrix_unified.Print();
            
            TDecompChol chol_unified(covmatrix_unified);
            chol_unified.Decompose();
            
            TMatrixD cmat_unified(chol_unified.GetU());
            //  cmat.Print();
            TMatrixD tcmat_unified(cmat_unified.Transpose(cmat_unified));
            // std::cout << "Cholesky matrix---------" << std::endl;
            // tcmat.Print();
            
            Double_t* tmp_matrix_unified = tcmat_unified.GetMatrixArray();
            
            for (Int_t i = 0; i < 4; i++)
            {
                for (Int_t j = 0; j < 4; j++)
                {
                    m_unified_nl_par_covmatrix_l[i][j] = tmp_matrix_unified[i*4 + j];
                    //      cout << "\t" << m_unified_nl_par_covmatrix_l[i][j];
                }
            }
        }
            break;
    }
}

void CreateEnergyMatrix :: SetNLParameters()
{
    switch (NLModelE)
    {
        case 0://BCW NL Model
            std::cout << "Using BCW NL Model"<< std::endl;
            
            for (Int_t i = 0; i < 5; i++)
            {
                m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];//reset
            }
            for (Int_t AD = 0; AD < NADs; AD++)
            {
                NLF->SetParameters(m_bcw_elec_nl_par[0],m_bcw_elec_nl_par[1],m_bcw_elec_nl_par[2], m_bcw_elec_nl_par[3], m_bcw_elec_nl_par[4], m_abs_escale * m_rel_escale[AD], m_abs_eoffset+m_rel_eoffset[AD]);
            }
            
            break;
        case 1://LBNL NL Model
            std::cout << "Using LBNL NL Model"<< std::endl;
            
            for (Int_t i = 0; i < 3; i++)
            {
                m_lbnl_nl_par[i] = m_lbnl_nl_par_nominal[i];//reset
            }
            for (Int_t AD = 0; AD < NADs; AD++)
            {
                NLF->SetParameters(m_lbnl_nl_par[0],m_lbnl_nl_par[1],m_lbnl_nl_par[2],m_abs_escale * m_rel_escale[AD],m_abs_eoffset + m_rel_eoffset[AD]);
            }
            
            break;
        case 2://Unified NL Model
            std::cout << "Using Unified NL Model"<< std::endl;
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
    Interpolation(NLF);//Interpolate the new function
}

void CreateEnergyMatrix :: LoadIavCorrection()// From Bryce Littlejohn's results. //I don't use different IAV matrix for each AD since the analysis it's been done for just 1 of them, assume identical ADs.
{
    TFile* f = TFile::Open("./IavDistortion/IAVDistortion.root");
    TH2F* Correction = (TH2F*)f->Get("Correction");
    
    std::cout << "Reading IAV correction file" << std::endl;
    
    for(Int_t i=0;i<MatrixBins;i++)
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
        for(Int_t i=0;i<MatrixBins;i++)
        {
            IAVMatrixFrac[AD][i] = NominalIAVMatrixFrac[i];//   Copy that will be varied when calculated random IAV matrix.
            
            for(Int_t j=0;j<MatrixBins;j++)
            {
                IAVMatrix[AD][i][j] = NominalIAVMatrix[i][j];// Copy that will be varied when calculated random IAV matrix.
            }
        }
    }
    f->Close();
}

// First order true to visible function.
//(x-(Mn-Mp))*(1 - x/Mn*(1.0 - [0]*sqrt(1 - Me*Me/(x-(Mn-Mp))/(x-(Mn-Mp)))))- ((Mn-Mp)*(Mn-Mp) - Me*Me)/(2*Mn)-Me]

Double_t CreateEnergyMatrix :: VisibleEnergy1F(Double_t* energ, Double_t* par)
{
    Double_t Enu = energ[0];
    Double_t Delta = Mn-Mp;
    Double_t gamma0 = (Enu-Delta)/Me;
    
    if((1-1/gamma0/gamma0) <=0)//To avoid NaNs
    {
        return (Enu-Delta)*(1-Enu/Mp*(1-2.4*Enu/Mp))-(Delta*Delta-Me*Me)/(2*Mp)+Me;
    }
    
    Double_t ve0 = sqrt(1 - 1/gamma0/gamma0);
    Double_t costheta = -0.034*ve0+2.4*Enu/Mp;
    
    return (Enu-Delta)*(1-Enu/Mp*(1-costheta))-(Delta*Delta-Me*Me)/(2*Mp)+Me;
    
    //    Double_t Enu = energ[0];
    //    Double_t costheta = par[0];
    //
    //    Double_t Me = 0.510999; // MeV
    //    Double_t Mn = 939.565; // MeV
    //    Double_t Mp = 938.272; // MeV
    //    Double_t Delta = Mn-Mp;
    //
    //    Double_t Ee0 = Enu - Delta;
    //    // if (Ee0 < Me){    // not allowed!
    //    //   return -1;
    //    // }
    //
    //    Double_t gamma0 = Ee0/Me;
    //    Double_t beta0 = sqrt(1 - 1/gamma0/gamma0);
    //    Double_t y2 = (Delta*Delta - Me*Me)/2.;
    //
    //    if (costheta < -1 || costheta > 1) { // then use average cos theta
    //        costheta = - 0.034*beta0 + 2.4*Enu/Mp;
    //    }
    //
    //    Double_t Ee1 = Ee0 * (1 - Enu/Mp*(1.0 - costheta*beta0)) - y2/Mp;
    //
    //    //   if (Ee1 < Me){    // not allowed!
    //    //     return -1;
    //    //   }
    //
    //    return Ee1 + Me; // Add Me to include annihiration gamma energy
}

Double_t CreateEnergyMatrix :: VisibleEnergy0F(Double_t* energ, Double_t* par)
{
    return energ[0]-par[0];
}

Double_t CreateEnergyMatrix :: ResolutionF(Double_t* x, Double_t* par)
{
    Double_t e_orig = x[0];
    Double_t e_sigma = 1.0;
    
    if (e_orig > 0)//To avoid NaNs.
    {
        e_sigma = TMath::Sqrt(par[0]*par[0] + par[1]*par[1]/e_orig + par[2]*par[2]/e_orig/e_orig);
    }
    
    return e_sigma;
}

Double_t CreateEnergyMatrix :: NLBCWF(Double_t* x, Double_t* par)
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
    if (par[4] > 0)
    {
        err_shift = par[4]*par4_up;
    }
    else err_shift = par[4]*par4_down;
    
    Double_t electronicsCorrection = exp(par[0] + par[1] * visibleE) + par[2] + err_offset + err_shift;
    
    Double_t final_energy =  visibleE * electronicsCorrection * escale_par + escale_offset;
    
    return final_energy;
}

// LBNL non-linearityr function, based on beta-gamma spectra
Double_t CreateEnergyMatrix :: NLLBNLF(Double_t * x, Double_t * par)
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
    Double_t random_nl_fac = scinti_nl_fac;
    
    for (Int_t ierr = 0; ierr < 3; ierr++)
    {
        random_nl_fac += par[ierr]*err[ierr];
    }
    
    Double_t visibleE = random_nl_fac * e_positron_true;
    //  std::cout << "NL " << e_positron_true << " " << visibleE/e_positron_true << std::endl;
    
    return visibleE  * escale_par + escale_offset;
}

Double_t CreateEnergyMatrix :: NLUnifiedF(Double_t * x, Double_t * par)
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
    Double_t visibleE = random_nl_fac * e_positron_true;
    return visibleE  * escale_par + escale_offset;
}

void CreateEnergyMatrix :: GenerateLBNLEnergyMatrix(Double_t sin22t13, Double_t dm2_31, Int_t week)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                           Set kinematic energy shift function
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    GetNeutrinoToPositronFunction(1);//0 for 0th order, 1 for 1st order
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                           Set nominal IAV Matrix
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    LoadIavCorrection();
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                            Set nominal NL parameters
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    this->SetBCWModel(BCW);
    this->SetLBNLModel(LBNL);
    this->SetUnifiedModel(Unified);
    
    LoadNLParameters();//nominal parameters
    SetNLParameters();
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                              Set nominal resolution parameters
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ResoF = new TF1("ResoF",this,&CreateEnergyMatrix::ResolutionF,0,20,3,"CreateEnergyMatrix","ResolutionF");
    //    ResoF->SetParameters(0.022,0.077,0.018); // based on Bryce's TN
    ResoF->SetParameters(0.0148,0.0869,0.0271); // based on Doc-8982
    
    //BCW values http://dayabay.ihep.ac.cn/DocDB/0087/008768/013/6AdAnalysis-BCW.pdf are different! 0.13226034931, 0.32604048828, 0.26196488314
    
    TH2D* h_evis_vs_enu;
    TH2D* h_evis_vs_enu_pos;
    TH2D* h_evis_vs_enu_iav;
    TH2D* h_evis_vs_enu_nl;
    TH2D* h_evis_vs_enu_reso;

    h_evis_vs_enu = new TH2D("EvisEnu","EvisEnu",240,0,12,n_evis_bins,evis_bins);
    h_evis_vs_enu_pos = new TH2D("EvisEnuPos","EvisEnuPos",240,0,12,n_evis_bins,evis_bins);
    h_evis_vs_enu_iav = new TH2D("EvisEnuIAV","EvisEnuIAV",240,0,12,n_evis_bins,evis_bins);
    h_evis_vs_enu_nl = new TH2D("EvisEnuNL","EvisEnuNL",240,0,12,n_evis_bins,evis_bins);
    h_evis_vs_enu_reso = new TH2D("EvisEnuReso","EvisEnuReso",240,0,12,n_evis_bins,evis_bins);
    
    TAxis * xa = h_evis_vs_enu->GetXaxis();
    
    for(int ibin_enu=0;ibin_enu< xa->GetNbins();++ibin_enu)
    {
        std::cout << ibin_enu << " / " << xa->GetNbins() << std::endl;
        
        this->updatePositronTrue(xa->GetBinLowEdge(ibin_enu+1),xa->GetBinUpEdge(ibin_enu+1));
        
        for(int ibin=0;ibin<MatrixBins;++ibin)
        {
            h_evis_vs_enu->Fill(xa->GetBinCenter(ibin_enu+1), EnergyVector[ibin], m_positronDetectedSpectrum[ibin]*m_binWidth); // don't include background spectra for this true enu_vs_evis convertion.
            
            h_evis_vs_enu_pos->Fill(xa->GetBinCenter(ibin_enu+1), EnergyVector[ibin], m_positronTrueSpectrum[ibin]*m_binWidth); // don't include background spectra for this true enu_vs_evis convertion.
            
            h_evis_vs_enu_iav->Fill(xa->GetBinCenter(ibin_enu+1), EnergyVector[ibin], m_positronIavDistortedSpectrum[ibin]*m_binWidth); // don't include background spectra for this true enu_vs_evis convertion.
            
            h_evis_vs_enu_nl->Fill(xa->GetBinCenter(ibin_enu+1), EnergyVector[ibin], m_positronNLSpectrum[ibin]*m_binWidth); // don't include background spectra for this true enu_vs_evis convertion.
            h_evis_vs_enu_reso->Fill(xa->GetBinCenter(ibin_enu+1), EnergyVector[ibin], m_positronDetectedSpectrum[ibin]*m_binWidth); // don't include background spectra for this true enu_vs_evis convertion.
            
        }
        
    } //ibin_enu loop
    
    
    std::ofstream fout("matrix_evis_to_enu_unified.txt");
    
    TH2F* RebinHist = new TH2F("Energy","Energy",n_etrue_bins,enu_bins,n_evis_bins,evis_bins);
    TH2F* RebinHistpos = new TH2F("Energypos","Energypos",n_etrue_bins,enu_bins,n_evis_bins,evis_bins);
    TH2F* RebinHistiav = new TH2F("Energyiav","Energyiav",n_etrue_bins,enu_bins,n_evis_bins,evis_bins);
    TH2F* RebinHistnl = new TH2F("Energynl","Energynl",n_etrue_bins,enu_bins,n_evis_bins,evis_bins);
    TH2F* RebinHistreso = new TH2F("Energyreso","Energyreso",n_etrue_bins,enu_bins,n_evis_bins,evis_bins);
    
    for (Int_t iEvisBin = 0; iEvisBin < n_evis_bins; iEvisBin++)
    {
        //    h_enu_tmp->Reset();
        TH1D * htmp = h_evis_vs_enu->ProjectionX(Form("h_enu_%d",iEvisBin),iEvisBin+1,iEvisBin+1);
        TH1D * htmppos = h_evis_vs_enu_pos->ProjectionX(Form("h_enu_pos_%d",iEvisBin),iEvisBin+1,iEvisBin+1);
        TH1D * htmpiav = h_evis_vs_enu_iav->ProjectionX(Form("h_enu_iav_%d",iEvisBin),iEvisBin+1,iEvisBin+1);
        TH1D * htmpnl = h_evis_vs_enu_nl->ProjectionX(Form("h_enu_nl_%d",iEvisBin),iEvisBin+1,iEvisBin+1);
        TH1D * htmpreso = h_evis_vs_enu_reso->ProjectionX(Form("h_enu_reso_%d",iEvisBin),iEvisBin+1,iEvisBin+1);
        
        TH1D * htmp_rebin = (TH1D*)htmp->Rebin(n_etrue_bins,Form("h_enu_rebin_%d",iEvisBin),enu_bins);
        TH1D * htmp_rebinpos = (TH1D*)htmppos->Rebin(n_etrue_bins,Form("h_enu_rebin_pos_%d",iEvisBin),enu_bins);
        TH1D * htmp_rebiniav = (TH1D*)htmpiav->Rebin(n_etrue_bins,Form("h_enu_rebin_iav_%d",iEvisBin),enu_bins);
        TH1D * htmp_rebinnl = (TH1D*)htmpnl->Rebin(n_etrue_bins,Form("h_enu_rebin_nl_%d",iEvisBin),enu_bins);
        TH1D * htmp_rebinreso = (TH1D*)htmpreso->Rebin(n_etrue_bins,Form("h_enu_rebin_reso_%d",iEvisBin),enu_bins);
        
        Double_t norm = htmp_rebin->Integral();
        Double_t normpos = htmp_rebinpos->Integral();
        Double_t normiav = htmp_rebiniav->Integral();
        Double_t normnl = htmp_rebinnl->Integral();
        Double_t normreso = htmp_rebinreso->Integral();
        
        htmp_rebin->Scale(1./norm);
        htmp_rebinpos->Scale(1./normpos);
        htmp_rebiniav->Scale(1./normiav);
        htmp_rebinnl->Scale(1./normnl);
        htmp_rebinreso->Scale(1./normreso);
        
        for (Int_t iEnuBin = 0; iEnuBin < n_etrue_bins; iEnuBin++)
        {
            RebinHist->SetBinContent(iEvisBin+1,iEnuBin+1,htmp_rebin->GetBinContent(iEnuBin+1));
            RebinHistpos->SetBinContent(iEvisBin+1,iEnuBin+1,htmp_rebinpos->GetBinContent(iEnuBin+1));
            RebinHistiav->SetBinContent(iEvisBin+1,iEnuBin+1,htmp_rebiniav->GetBinContent(iEnuBin+1));
            RebinHistnl->SetBinContent(iEvisBin+1,iEnuBin+1,htmp_rebinnl->GetBinContent(iEnuBin+1));
            RebinHistreso->SetBinContent(iEvisBin+1,iEnuBin+1,htmp_rebinreso->GetBinContent(iEnuBin+1));
        }
    }
    
    fout.close();
    
    std::cout << "\t Generating Energy Matrix" << std::endl;
    
    //
    //    std::ifstream input("matrix_evis_to_enu_unified.txt");
    //    std::string line;
    //
    //    Double_t Energy[n_etrue_bins][n_evis_bins];
    //
    //    for(Int_t i=0;i<n_evis_bins;i++)
    //    {
    //        for(Int_t j=0;j<n_etrue_bins;j++)
    //        {
    //            input >> Energy[j][i];
    //        }
    //    }
    //    for(Int_t i=0;i<n_etrue_bins;i++)
    //    {
    //        for(Int_t j=0;j<n_evis_bins;j++)
    //        {
    //
    //            Hist->SetBinContent(i+1,j+1,Energy[i][j]);
    //        }
    //    }
    
    TFile* EnergyMatrixDataF = new TFile(("./ResponseMatrices/"+AnalysisString+"/NominalResponseMatrix.root").c_str(),"recreate");
    
    RebinHist->Write("EnuEvis");
    RebinHistpos->Write("EvisEnuPos");
    RebinHistiav->Write("EvisEnuIAV");
    RebinHistnl->Write("EvisEnuNL");
    RebinHistreso->Write("EvisEnuReso");
    
    h_evis_vs_enu->Write();
    
    
    delete EnergyMatrixDataF;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                          Matrices with the original 240x240 binning
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //    if(RebinEnergyMatrix)
    //    {
    //        std::cout << "\t Rebin Matrix from Fine to Coarse binning" << std::endl;
    //
    //        LimitVis = n_evis_bins;
    //        LimitTrue = n_etrue_bins;
    //
    //        EnergyMatrixH = new TH2D("EvisEnu","EvisEnu",LimitTrue,enu_bins,LimitVis,evis_bins);
    //        EnergyMatrixPosH = new TH2D("EvisEnuPos","EvisEnuPos",LimitTrue,enu_bins,LimitVis,evis_bins);
    //        EnergyMatrixIAVH = new TH2D("EvisEnuIAV","EvisEnuIAV",LimitTrue,enu_bins,LimitVis,evis_bins);
    //        EnergyMatrixNLH = new TH2D("EvisEnuNL","EvisEnuNL",LimitTrue,enu_bins,LimitVis,evis_bins);
    //        EnergyMatrixResoH = new TH2D("EvisEnuReso","EvisEnuReso",LimitTrue,enu_bins,LimitVis,evis_bins);
    //        RowEnergyMatrixH = new TH2D("EnuEvis","EnuEvis", LimitTrue,enu_bins,LimitVis,evis_bins);
    //
    //        for (Int_t i = 0; i < n_etrue_bins; i++)
    //        {
    //            OscDeltaPositronSpectrumSumH[i] = new TH1D(Form("Fine Positron Spectrum Vis for True Energy Index %d",i),Form("Fine Positron Spectrum Vis for True Energy Index %d",i),MatrixBins,0,FinalVisibleEnergy);
    //            OscDeltaIAVSpectrumSumH[i]= new TH1D(Form("Fine IAV Spectrum Vis for True Energy Index %d",i),Form("Fine Spectrum IAV Vis for True Energy Index %d",i),MatrixBins,0,FinalVisibleEnergy);
    //            OscDeltaNLSpectrumSumH[i]= new TH1D(Form("Fine NL Spectrum Vis for True Energy Index %d",i),Form("Fine Spectrum NL Vis for True Energy Index %d",i),MatrixBins,0,FinalVisibleEnergy);
    //
    //            OscDeltaVisibleSpectrumSumH[i]= new TH1D(Form("Fine Visible Spectrum Vis for True Energy Index %d",i),Form("Fine Visible Spectrum Vis for True Energy Index %d",i),MatrixBins,0,FinalVisibleEnergy);
    //
    //            for (Int_t TrueEnergyIndex = Int_t(enu_bins[i]*MatrixBins/FinalVisibleEnergy); TrueEnergyIndex < Int_t(enu_bins[i+1]*MatrixBins/FinalVisibleEnergy); TrueEnergyIndex++)//Although the results are integers I need to add Int_t() if I don't want the precision errors to mess up the rebinning.
    //            {
    //                OscDeltaPositronSpectrumSumH[i]->Add(OscDeltaPositronSpectrumH[TrueEnergyIndex]);
    //                OscDeltaIAVSpectrumSumH[i]->Add(OscDeltaIAVSpectrumH[TrueEnergyIndex]);
    //                OscDeltaNLSpectrumSumH[i]->Add(OscDeltaNLSpectrumH[TrueEnergyIndex]);
    //                OscDeltaVisibleSpectrumSumH[i]->Add(OscDeltaVisibleSpectrumH[TrueEnergyIndex]);
    //            }
    //
    //            RebinnedOscDeltaPositronSpectrumSumH[i]=(TH1D*)OscDeltaPositronSpectrumSumH[i]->Rebin(n_evis_bins,Form("Rebinned Positron Spectrum Vis for True Energy Index %d",i),evis_bins);
    //            RebinnedOscDeltaIAVSpectrumSumH[i]=(TH1D*)OscDeltaIAVSpectrumSumH[i]->Rebin(n_evis_bins,Form("Rebinned IAV Spectrum Vis for True Energy Index %d",i),evis_bins);
    //            RebinnedOscDeltaNLSpectrumSumH[i]=(TH1D*)OscDeltaNLSpectrumSumH[i]->Rebin(n_evis_bins,Form("Rebinned NL Spectrum Vis for True Energy Index %d",i),evis_bins);
    //            RebinnedOscDeltaVisibleSpectrumSumH[i]=(TH1D*)OscDeltaVisibleSpectrumSumH[i]->Rebin(n_evis_bins,Form("Rebinned Visible Spectrum Vis for True Energy Index %d",i),evis_bins);
    //
    //            for (Int_t j = 0; j < n_evis_bins; j++)
    //            {
    //                EnergyMatrixH->SetBinContent(i+1,j+1, RebinnedOscDeltaVisibleSpectrumSumH[i]->GetBinContent(j+1));
    //                EnergyMatrixPosH->SetBinContent(i+1,j+1,RebinnedOscDeltaPositronSpectrumSumH[i]->GetBinContent(j+1));
    //                EnergyMatrixIAVH->SetBinContent(i+1,j+1,RebinnedOscDeltaIAVSpectrumSumH[i]->GetBinContent(j+1));
    //                EnergyMatrixNLH->SetBinContent(i+1,j+1,RebinnedOscDeltaNLSpectrumSumH[i]->GetBinContent(j+1));
    //                EnergyMatrixResoH->SetBinContent(i+1,j+1,RebinnedOscDeltaVisibleSpectrumSumH[i]->GetBinContent(j+1));
    //            }
    //        }
    //        for (Int_t i = 0; i < n_etrue_bins; i++)
    //        {
    //            for (Int_t j = 0; j < n_evis_bins; j++)
    //            {
    //                RowEnergyMatrixH->SetBinContent(i+1,j+1,EnergyMatrixH->GetBinContent(i+1,j+1));
    //
    //                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                //                          Maybe instead of transposing the correct option is to invert it. Study this!!
    //                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //            }
    //        }
    //    }
    //    else
    //    {
    //        LimitVis = MatrixBins;
    //        LimitTrue = MatrixBins;
    //
    //        EnergyMatrixH= new TH2D("EvisEnu","EvisEnu",LimitTrue,0,LimitTrue,LimitVis,0,LimitVis);
    //        EnergyMatrixPosH= new TH2D("EvisEnuPos","EvisEnuPos",LimitTrue,0,LimitTrue,LimitVis,0,LimitVis);
    //        EnergyMatrixIAVH= new TH2D("EvisEnuIAV","EvisEnuIAV",LimitTrue,0,LimitTrue,LimitVis,0,LimitVis);
    //        EnergyMatrixNLH= new TH2D("EvisEnuNL","EvisEnuNL",LimitTrue,0,LimitTrue,LimitVis,0,LimitVis);
    //        EnergyMatrixResoH= new TH2D("EvisEnuReso","EvisEnuReso",LimitTrue,0,LimitTrue,LimitVis,0,LimitVis);
    //        RowEnergyMatrixH= new TH2D("EnuEvis","EnuEvis",LimitVis,0,LimitVis, LimitTrue,0,LimitTrue);
    //
    //        for (Int_t i = 0; i < MatrixBins; i++)
    //        {
    //            for (Int_t j = 0; j < MatrixBins; j++)
    //            {
    //                EnergyMatrixH->SetBinContent(i+1,j+1,EnergyMatrixH->GetBinContent(i+1,j+1)+OscDeltaVisibleSpectrumH[i]->GetBinContent(j+1));
    //
    //                EnergyMatrixPosH->SetBinContent(i+1,j+1,EnergyMatrixPosH->GetBinContent(i+1,j+1)+OscDeltaPositronSpectrumH[i]->GetBinContent(j+1));
    //
    //                EnergyMatrixIAVH->SetBinContent(i+1,j+1,EnergyMatrixIAVH->GetBinContent(i+1,j+1)+OscDeltaIAVSpectrumH[i]->GetBinContent(j+1));
    //
    //                EnergyMatrixNLH->SetBinContent(i+1,j+1,EnergyMatrixNLH->GetBinContent(i+1,j+1)+OscDeltaNLSpectrumH[i]->GetBinContent(j+1));
    //
    //                EnergyMatrixResoH->SetBinContent(i+1,j+1,EnergyMatrixResoH->GetBinContent(i+1,j+1)+OscDeltaVisibleSpectrumH[i]->GetBinContent(j+1));
    //            }
    //        }
    //        for (Int_t i = 0; i < MatrixBins; i++)
    //        {
    //            for (Int_t j = 0; j < MatrixBins; j++)
    //            {
    //                RowEnergyMatrixH->SetBinContent(j+1,i+1,EnergyMatrixH->GetBinContent(i+1,j+1));
    //
    //                ////////////////////////////////////////////////////////////////////////////////////////
    //                //        Maybe instead of transposing the correct option is to invert it. Study this!!
    //                ////////////////////////////////////////////////////////////////////////////////////////
    //
    //            }
    //        }
    //
    //    }
    //
    //    // UNCOMMENT FOLLOWING LINES TO SAVE ENERGY MATRIX SLICES IN EACH PRODUCTION STEP
    //    //    for(Int_t TrueEnergyIndex = 0; TrueEnergyIndex<MatrixBins; TrueEnergyIndex++)
    //    //    {
    //    //        OscDeltaPositronSpectrumH[TrueEnergyIndex]->Write();
    //    //        OscDeltaIAVSpectrumH[TrueEnergyIndex]->Write();
    //    //        OscDeltaNLSpectrumH[TrueEnergyIndex]->Write();
    //    //        OscDeltaVisibleSpectrumH[TrueEnergyIndex]->Write();
    //    //    }
    //    if(RebinEnergyMatrix)
    //    {
    //        for(Int_t TrueEnergyIndex = 0; TrueEnergyIndex<n_etrue_bins; TrueEnergyIndex++)
    //        {
    //            //        OscDeltaPositronSpectrumSumH[TrueEnergyIndex]->Write();
    //            //        OscDeltaIAVSpectrumSumH[TrueEnergyIndex]->Write();
    //            //        OscDeltaNLSpectrumSumH[TrueEnergyIndex]->Write();
    //            //        OscDeltaVisibleSpectrumSumH[TrueEnergyIndex]->Write();
    //
    //            delete OscDeltaPositronSpectrumSumH[TrueEnergyIndex];
    //            delete OscDeltaIAVSpectrumSumH[TrueEnergyIndex];
    //            delete OscDeltaNLSpectrumSumH[TrueEnergyIndex];
    //            delete OscDeltaVisibleSpectrumSumH[TrueEnergyIndex];
    //
    //            delete RebinnedOscDeltaPositronSpectrumSumH[TrueEnergyIndex];
    //            delete RebinnedOscDeltaIAVSpectrumSumH[TrueEnergyIndex];
    //            delete RebinnedOscDeltaNLSpectrumSumH[TrueEnergyIndex];
    //            delete RebinnedOscDeltaVisibleSpectrumSumH[TrueEnergyIndex];
    //        }
    //    }
    //    //
    //    NoNormalizedEnergyMatrixH = (TH2D*)EnergyMatrixH->Clone("NoNormalizedEvisEnu");
    //    NoNormalizedEnergyMatrixPosH = (TH2D*)EnergyMatrixPosH->Clone("PosNoNormalizedEvisEnu");
    //    NoNormalizedEnergyMatrixIAVH = (TH2D*)EnergyMatrixIAVH->Clone("IAVNoNormalizedEvisEnu");
    //    NoNormalizedEnergyMatrixNLH = (TH2D*)EnergyMatrixNLH->Clone("NLNoNormalizedEvisEnu");
    //    NoNormalizedEnergyMatrixResoH = (TH2D*)EnergyMatrixResoH->Clone("ResoNoNormalizedEvisEnu");
    //
    //    NoNormalizedEnergyMatrixPosH->Write();
    //    NoNormalizedEnergyMatrixIAVH->Write();
    //    NoNormalizedEnergyMatrixNLH->Write();
    //    NoNormalizedEnergyMatrixResoH->Write();
    //
    //    NoNormalizedEnergyMatrixH->Write();//Save Matrix before rebinning
    //    RowEnergyMatrixH->Write("NoNormalizedEnuEvis");//Save Matrix before rebinning
    //
    //    //Normalize the Matrix
    //    for(Int_t j=0;j<LimitTrue;j++)
    //    {
    //        Norma[j]=0;
    //        NormaPos[j]=0;
    //        NormaIAV[j]=0;
    //        NormaNL[j]=0;
    //        NormaReso[j]=0;
    //        for(Int_t i=0;i<LimitVis;i++)
    //        {
    //            Norma[j] = Norma[j]+EnergyMatrixH->GetBinContent(j+1,i+1);
    //            NormaPos[j] = NormaPos[j]+EnergyMatrixPosH->GetBinContent(j+1,i+1);
    //            NormaIAV[j] = NormaIAV[j]+EnergyMatrixIAVH->GetBinContent(j+1,i+1);
    //            NormaNL[j] = NormaNL[j]+EnergyMatrixNLH->GetBinContent(j+1,i+1);
    //            NormaReso[j] = NormaReso[j]+EnergyMatrixResoH->GetBinContent(j+1,i+1);
    //        }
    //    }
    //
    //    for(Int_t i=0;i<LimitVis;i++)
    //    {
    //        NormaTrans[i]=0;
    //        for(Int_t j=0;j<LimitTrue;j++)
    //        {
    //            NormaTrans[i] = NormaTrans[i]+RowEnergyMatrixH->GetBinContent(j+1,i+1);
    //        }
    //    }
    //
    //    for (Int_t i = 0; i < LimitVis; i++)
    //    {
    //        for (Int_t j = 0; j < LimitTrue; j++)
    //        {
    //            if(Norma[j]!=0)
    //            {
    //                EnergyMatrixH->SetBinContent(j+1,i+1,EnergyMatrixH->GetBinContent(j+1,i+1)/Norma[j]);//Normalization so Σj E(i,j) = 1; (Σ(y axis) =1)
    //            }
    //            if(NormaPos[j]!=0)
    //            {
    //                EnergyMatrixPosH->SetBinContent(j+1,i+1,EnergyMatrixPosH->GetBinContent(j+1,i+1)/NormaPos[j]);
    //            }
    //            if(NormaIAV[j]!=0)
    //            {
    //                EnergyMatrixIAVH->SetBinContent(j+1,i+1,EnergyMatrixIAVH->GetBinContent(j+1,i+1)/NormaIAV[j]);
    //            }
    //            if(NormaNL[j]!=0)
    //            {
    //                EnergyMatrixNLH->SetBinContent(j+1,i+1,EnergyMatrixNLH->GetBinContent(j+1,i+1)/NormaNL[j]);
    //            }
    //            if(NormaReso[j]!=0)
    //            {
    //                EnergyMatrixResoH->SetBinContent(j+1,i+1,EnergyMatrixResoH->GetBinContent(j+1,i+1)/NormaReso[j]);
    //            }
    //            if(NormaTrans[i]!=0)
    //            {
    //                RowEnergyMatrixH->SetBinContent(j+1,i+1,RowEnergyMatrixH->GetBinContent(j+1,i+1)/NormaTrans[i]);
    //            }
    //        }
    //    }
    //    //    }
    //    //Save Matrix after rebinning
    //    EnergyMatrixH->Write();
    //    EnergyMatrixPosH->Write();
    //    EnergyMatrixIAVH->Write();
    //    EnergyMatrixNLH->Write();
    //    EnergyMatrixResoH->Write();
    //    RowEnergyMatrixH->Write();
    //
    //    EnergyMatrixDataF->Close();
    //
    //    delete NoNormalizedEnergyMatrixPosH;
    //    delete NoNormalizedEnergyMatrixIAVH;
    //    delete NoNormalizedEnergyMatrixNLH;
    //    delete NoNormalizedEnergyMatrixResoH;
    //    delete NoNormalizedEnergyMatrixH;
    //
    //    delete EnergyMatrixIAVH;
    //    delete EnergyMatrixNLH;
    //    delete EnergyMatrixPosH;
    //    delete EnergyMatrixResoH;
    
    delete ResoF;
    std::cout << "\t Finished Generating Energy Matrix" << std::endl;
}

void CreateEnergyMatrix::updatePositronTrue(Double_t eNu_min,Double_t eNu_max)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //                                          Calculate distortions
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    for(int idx=0; idx<m_nSamples; idx++)
    {
        AntineutrinoSpectrum[idx]=1;//Flat spectrum
    }
    
    for(int idx=0; idx<m_nSamples; idx++)
    {
        EnergyVector[idx] = (0.5+idx)*m_binWidth+m_eMin;
        
        double e_positron = EnergyVector[idx];
        if (e_positron < 1.022){
            m_positronTrueSpectrum[idx] = 0;
            continue;
        }
        double e_nu = VisibleF->GetX(e_positron);
        int e_nu_Idx = (int)((e_nu-m_eMin)/m_binWidth);
        if(e_nu_Idx<0 || e_nu_Idx>=m_nSamples){
            m_positronTrueSpectrum[idx] = 0;
            continue;
        }
        
        if (eNu_min >= 0 && (e_nu < eNu_min || e_nu > eNu_max)){ // select single true Enu bin if specified.
            m_positronTrueSpectrum[idx] = 0;
            continue;
        }
        
        double binScaling = VisibleF->Derivative(e_nu,0,0.0001);
        double dNdE = 0;
        int sign = 1;
        double dE = e_nu - (0.5 + e_nu_Idx)*m_binWidth;
        if (dE < 0) sign = -1;
        
        if(e_nu_Idx==m_nSamples-1){
            dNdE = AntineutrinoSpectrum[e_nu_Idx];
        }else{
            dNdE =
            (TMath::Abs(dE)/m_binWidth)
            *AntineutrinoSpectrum[e_nu_Idx+sign]
            +(1 - TMath::Abs(dE)/m_binWidth)
            * AntineutrinoSpectrum[e_nu_Idx];
        }
        //        cout << energy_nl << "\t" << energy_true << "\t" << dNdE/binScaling << endl;
        
        m_positronTrueSpectrum[idx] = dNdE / binScaling;
    }
    
    // ///////////////
    
    //    cout << "Updated total positron true counts: " << totalCounts << endl;
    
    // Apply IAV correction
    for(int idx=0; idx<m_nSamples; idx++)
    {
        m_positronIavDistortedSpectrum[idx] = 0;
    }
    for(int idx=0; idx<m_nSamples; idx++)
    {
        for(int jdx=0; jdx<idx+1; jdx++)
        {
            m_positronIavDistortedSpectrum[jdx]
            += IAVMatrix[0][idx][jdx] * m_positronTrueSpectrum[idx];
        }
        
    }
    
    
    
    for(int idx=0; idx<m_nSamples; idx++)
    {
        double energy_nl = EnergyVector[idx];
        double energy_true = NLF->GetX(energy_nl);
        
        int eTrueIdx = (int)((energy_true-m_eMin)/m_binWidth);
        if(eTrueIdx<0 || eTrueIdx>=m_nSamples)
        {
            m_positronNLSpectrum[idx] = 0;
            continue;
        }
        
        double binScaling = NLF->Derivative(energy_true);
        double dNdE = 0;
        int sign = 1;
        double dEtrue = energy_true - (0.5 + eTrueIdx)*m_binWidth;
        if (dEtrue < 0) sign = -1;
        
        if(eTrueIdx==m_nSamples-1)
        {
            dNdE = m_positronIavDistortedSpectrum[eTrueIdx];
            // }if(TMath::Abs(m_positronIavDistortedSpectrum[eTrueIdx+sign] -m_positronIavDistortedSpectrum[eTrueIdx])/m_positronIavDistortedSpectrum[eTrueIdx] > 0.5){
            //   dNdE = m_positronIavDistortedSpectrum[eTrueIdx];
        }
        else{
            // dNdE = (m_positronIavDistortedSpectrum[eTrueIdx]
            //         +((dEtrue/m_binWidth)
            //           *(m_positronIavDistortedSpectrum[eTrueIdx+1]
            //             -m_positronIavDistortedSpectrum[eTrueIdx])));
            
            // change the interpolation method better
            dNdE =
            (TMath::Abs(dEtrue)/m_binWidth)
            *m_positronIavDistortedSpectrum[eTrueIdx+sign]
            +(1 - TMath::Abs(dEtrue)/m_binWidth)
            * m_positronIavDistortedSpectrum[eTrueIdx];
        }
        //        cout << energy_nl << "\t" << energy_true << "\t" << dNdE/binScaling << endl;
        
        m_positronNLSpectrum[idx] = dNdE / binScaling;
        
        double resolutionRange = 8; // sigma
        for(int idx=0; idx<m_nSamples; idx++)
        {
            m_positronDetectedSpectrum[idx] = 0;
        }
        
        for(int idx=0; idx<m_nSamples; idx++)
        {
            
            if(m_positronNLSpectrum[idx]==0) continue;
            
            double energy_nl = EnergyVector[idx];
            //      double sigma = m_detectorResolution * TMath::Sqrt(energy_nl);
            double sigma = (ResoF->Eval(energy_nl) + ResolutionBias[0]) * energy_nl;
            double minDetE = energy_nl - resolutionRange*sigma;
            double maxDetE = energy_nl + resolutionRange*sigma;
            int minDetEIdx = (int)((minDetE-m_eMin)/m_binWidth);
            int maxDetEIdx = (int)((maxDetE-m_eMin)/m_binWidth);
            if(minDetEIdx < 0) minDetEIdx = 0;
            if(maxDetEIdx >= m_nSamples) maxDetEIdx = m_nSamples-1;
            for(int detIdx=minDetEIdx; detIdx<=maxDetEIdx; detIdx++)
            {
                if(detIdx==0) continue;
                double gausFactor = TMath::Gaus((energy_nl-EnergyVector[detIdx]),
                                                0,sigma,true);
                m_positronDetectedSpectrum[detIdx] += (m_positronNLSpectrum[idx]
                                                       * gausFactor
                                                       * m_binWidth);
            }
        }
    }
}

void CreateEnergyMatrix::LoadADSpectrum()
{
    TFile* OutputFile = new TFile(("./RootOutputs/"+AnalysisString+"/NominalOutputs/Oscillation.root").c_str());

    OutputFile->cd("Total AD Spectra after oscillation");

    for(Int_t AD = 0; AD<NADs; AD++)
    {
        TotalOscillatedSpectrumAD[AD] = (TH1D*)gDirectory->Get(Form("Total spectrum after oscillation at AD%i",AD+1));
    }
    delete OutputFile;
}