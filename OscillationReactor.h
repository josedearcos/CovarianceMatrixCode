#pragma once
#include "NominalData.h"
#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom3.h"
#include "nHToyMC.h"
// Used to debug, if you don't want it to run, uncomment the following line: (#define NDEBUG)
//#define NDEBUG
#include <assert.h>

//#define TestDiagonalMatrix
//
// Originally: Given a reactor antineutrino spectrum I output the respective spectrua for each AD after oscillation.
// Updated: Class that is able to calculate the spectra in each AD either from Reactor Data or Near Hall Data. 5/3/13.
//
//  Created by Jose de Arcos	 on 3/10/13.
//
//

//const Int_t NReactors=6;

//NL
const Int_t n_bcw_positron_nl = 1000;
const Int_t n_unified_nl_points = 500;
const Int_t m_num_unified_nl_pars = 4; // 4 marginal curves in the final model

//Particle masses
const Double_t Me = 0.510999; // MeV
const Double_t Mn = 939.565; // MeV
const Double_t Mp = 938.272; // MeV

class OscillationReactor
{
private:
    
    std::string NominalPredictionDirectory;
    std::string RandomPredictionDirectory;
    
    //External classes used
    NominalData* Nom;
    TRandom3* rand;
    bool isH;
    Int_t DataSet;
    Int_t week;
    bool CovMatrix;
    std::string AnalysisString;
#ifdef UseVolumes
    std::string VolumeString[2];
#endif
    //Survival probability calculation parameters:
    Int_t hierarchy;
    Double_t deltam2_32;
    Double_t deltam2_21;
    Double_t deltam2_31;
    Double_t dm2_ee;
    Double_t s22t13;
    Double_t s22t12;
    
    //Detector properties
    Double_t ADdistances[MaxDetectors][NReactors];
    
    bool Mode;
    
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
    bool Sin22t12Matrix;
    
    bool RelativeEnergyScaleMatrix;
    bool IAVMatrixb;
    bool OAVMatrixb;
    bool NLMatrix;
    bool ResolutionMatrix;
    bool EfficiencyMatrix;
    bool AllDetectorSystematicsMatrix;
    
    //AD configuration parameters:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    
    //Cell parameters:
    Int_t XCellLimit;
    Int_t YCellLimit;
    
    //Binning parameters:
    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    
    Double_t DetectorProtons[MaxDetectors][VolumeX][VolumeY];

    Double_t NominalDetectorEfficiency[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    Double_t MultiMuonEff[MaxDetectors][MaxPeriods];
    Double_t FullTime[MaxDetectors][MaxPeriods];

    Double_t DetectorEfficiencyRelativeError;
    
    Double_t BinWidth;
    Int_t TotalBins;
    
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Int_t n_evis_bins;
    Int_t n_etrue_bins;
    
    Int_t Nweeks;
    
    TH1D* OscillatedSpectrumAD[MaxDetectors];
    TH1D* TotalOscillatedSpectrumAD[MaxDetectors];
    Double_t IntegralNominalTotalOscillatedSpectrumAD[MaxDetectors];
    TH1D* NominalTotalOscillatedSpectrumAD[MaxDetectors];
    TH1D* ScaledOscillatedSpectrumAD[MaxDetectors][VolumeX][VolumeY];
    TH1D* OscillatedSpectrumEH1;
    TH1D* OscillatedSpectrumEH2;
    TH1D* OscillatedSpectrumEH3;
    
    TH1D* RatioEH1EH2;
    TH1D* RatioEH1EH3;
    TH1D* RatioEH2EH3;
    
    TH1D* DifferenceEH1EH2;
    TH1D* DifferenceEH1EH3;
    TH1D* DifferenceEH2EH3;
    
    TH1D* ReactorSpectrumH[NReactors];
    TH1D* NearHallSpectrumH[MaxDetectors];
    
    //Files
    Char_t DistanceFileName[200];
    Char_t OutputFileName[200];
    Char_t ReactorData[200];
    
    void GenerateVisibleSpectrum();
    void LoadNominalTotalOscillatedSpectrum();
    
    void PlotEH();
    void PlotRatioEH();
    void PlotDifferenceEH();
    
    void LoadReactorHistograms();
    
    Double_t OscProb(Double_t, Double_t, Double_t, Double_t);
    
    //Detector effects:
    
    TH1D* OscDeltaPositronSpectrumH[MatrixBins];
    TH1D* OscDeltaIAVSpectrumH[MatrixBins];
    TH1D* OscDeltaNLSpectrumH[MatrixBins];
    TH1D* OscDeltaVisibleSpectrumH[MatrixBins];
    
    void RandomAbsoluteEnergyScaleMatrix();
    void RandomAbsoluteEnergyOffsetMatrix();
    
    void RandomRelativeEnergyScaleMatrix();
    void RandomRelativeEnergyOffsetMatrix();
    
    void RandomIAVMatrix();
    void RandomResolutionMatrix();
    
    Double_t PositronEnergy;
    Double_t EnergyVector[MatrixBins];
    
    enum NLModel{BCWE, LBNLE, UnifiedE};//Add here new NL models
    NLModel NLModelE;
    
    void  GetOscEnergyShift(Int_t,Int_t,Int_t);
    void  GetOscIAVShift(Int_t,Int_t,Int_t);
    void  GetOscNLShift(Int_t,Int_t,Int_t);
    void  GetOscResolutionShift(Int_t,Int_t,Int_t);
    
    void Interpolation(TF1*);
    
    Double_t VisibleEnergy0F(Double_t *, Double_t *);// Zeroth order true to visible function.
    Double_t VisibleEnergy1F(Double_t *, Double_t *);// First order true to visible function.
    
    Double_t NLBCWF(Double_t*, Double_t*);
    Double_t NLLBNLF(Double_t*, Double_t*);
    Double_t NLUnifiedF(Double_t*, Double_t*);
    
    Double_t ResolutionF(Double_t *, Double_t *);
    
    void LoadNLParameters();
    void SetNLParameters(bool);
    void SetUpDetectorResponse();
    
    //Load functions
    void LoadIavCorrection();
    
    //Functions to calculate the energy matrix
    TF1* VisibleF;
    TF1* GetNeutrinoToPositronFunction(bool);
    
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
    
    // Resolution
    // Detector resolution function
    TF1 * ResoF;
    Double_t IBDEvents[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    Double_t ObservedEvents[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
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
    Double_t DetectorEfficiencyDelayed[MaxDetectors];
    Double_t DetectorEfficiency[MaxDetectors*MaxPeriods*VolumeX*VolumeY];//necessary in vector format to return it in the getter function
    
    Double_t NominalDetectorEfficiencyDelayed;
    
    TH1D* BackgroundSpectrumH[MaxDetectors][VolumeX][VolumeY];
    
    void SetBCWModel(bool);
    void SetLBNLModel(bool);
    void SetUnifiedModel(bool);
    void SetSystematic();
    void RandomEfficiency();
    void RandomPower();
    
    TH1D* VisibleHisto[MaxDetectors][VolumeX][VolumeY];
    
    void LoadExternalInputs();
    void LoadNominalBackgrounds();
    
public:
    OscillationReactor();
    OscillationReactor(NominalData*);
    ~OscillationReactor();
    void OscillationFromReactorData(Int_t,bool,bool);
    TH1D* GetReactorOscillatedADSpectrum(Int_t,Int_t,Int_t,Int_t);
    void SetOscillationParameters(Double_t,Double_t);
    void SetRelativeEnergy(Double_t,Int_t);
    void SetSin22t12(Double_t);
    void FreeMemory();
    Double_t* GetDelayedEfficiency();
    Double_t* GetEfficiency();
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//DEFAULT CONSTRUCTOR

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Default values
OscillationReactor :: OscillationReactor()
{
    std::cout << " the oscillation reactor default constructor shouldn't be called" << std::endl;
    
    exit(EXIT_FAILURE);
    
    Nom = new NominalData(0,2);
    
    DataSet = Nom->GetDataSet();

    rand = new TRandom3(0);
    
    RandomPredictionDirectory = Nom->GetToyMCSamplesDirectory();
    NominalPredictionDirectory = Nom->GetPredictionDirectory();
    
    hierarchy=Nom->GetHierarchy();
    deltam2_32=2.30e-3;//eV2
    deltam2_21=7.50e-5;//eV2
    deltam2_31=deltam2_32+hierarchy*deltam2_21;
    
    s22t13 = Nom->GetSin22t13();
    dm2_ee = Nom->GetDm2ee();
    s22t12 = Nom->GetSin22t12();
    
    //Hard code these variables since the studies are provided in the following energies:
    InitialEnergy = 1.8;
    FinalEnergy = 12;
    InitialVisibleEnergy = 0;
    FinalVisibleEnergy = 12;
    
    TotalBins = MatrixBins;
    BinWidth=(FinalVisibleEnergy-InitialVisibleEnergy)/TotalBins;
    
    n_evis_bins = Nom->GetVisibleBins();
    n_etrue_bins = Nom->GetTrueBins();

    for (Int_t i = 0; i <= n_evis_bins; i++)
    {
        evis_bins[i] = Nom->GetVisibleBinningArray(i);
    }
    
    Nweeks = Nom->GetWeeks();
    isH = Nom->GetAnalysis();
    
    if(isH)
    {
        AnalysisString = "Hydrogen";
        XCellLimit = VolumeX;
        YCellLimit = VolumeY;
        Nom->ReadToyEventsByCell();
    }
    else
    {
        AnalysisString = "Gadolinium";
        XCellLimit = 1;
        YCellLimit = 1;
    }
#ifdef UseVolumes
    VolumeString[0] = "Gd-LS";
    VolumeString[1] = "LS";
#endif
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    NADs = Nom->GetADs();
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for(Int_t reactor = 0; reactor < NReactors; reactor++)
        {
            ADdistances[AD][reactor] = Nom->GetDistances(AD,reactor)/1000;//IN KILOMETERS
        }
        
        for (Int_t week = 0; week <Nweeks; week++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    IBDEvents[AD][week][idx][idy] = Nom->GetIBDEvents(AD,week,idx,idy);
                    ObservedEvents[AD][week][idx][idy] = Nom->GetObservedEvents(AD,week,idx,idy);
                }
            }
        }
    }
    
    VaryAccidentalMatrix = Nom->GetVaryAccidentalMatrix();
    VaryLiHeMatrix = Nom->GetVaryLiHeMatrix();
    VaryFastNeutronsMatrix = Nom->GetVaryFastNeutronsMatrix();
    VaryAmCMatrix = Nom->GetVaryAmCMatrix();
    DistortLiHeMatrix = Nom->GetDistortLiHeMatrix();
    DistortFastNeutronsMatrix = Nom->GetDistortFastNeutronsMatrix();
    DistortAmCMatrix = Nom->GetDistortAmCMatrix();
    
    IsotopeMatrix = Nom->GetIsotopeMatrix();
    ReactorPowerMatrix = Nom->GetReactorPowerMatrix();
    Sin22t12Matrix = Nom->GetSin22t12Matrix();
    RelativeEnergyScaleMatrix = Nom->GetRelativeEnergyScaleMatrix();
    IAVMatrixb = Nom->GetIAVMatrix();
    OAVMatrixb = Nom->GetOAVMatrix();
    NLMatrix = Nom->GetNLMatrix();
    ResolutionMatrix = Nom->GetResolutionMatrix();
    EfficiencyMatrix = Nom->GetEfficiencyMatrix();
    AllDetectorSystematicsMatrix = Nom->GetAllDetectorSystematicsMatrix();
    
    //  IAV Error from Bryce
    IAVNominalError=Nom->GetIAVError();
    //  Non uniformity
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        ResolutionBias[AD]=0;
        
        //  Relative energy scale
        m_rel_escale[AD] = Nom->GetRelativeEnergyScale();
        m_rel_escale_error[AD] = Nom->GetRelativeEnergyError();
        m_rel_escale_nominal[AD] = m_rel_escale[AD];
        m_rel_eoffset[AD] = Nom->GetRelativeEnergyOffset();
    }
    m_abs_escale = Nom->GetAbsoluteEnergyScale();
    m_abs_eoffset = Nom->GetAbsoluteEnergyOffset();
    
    m_abs_escale_nominal = m_abs_escale;
    
    m_abs_escale_error = Nom->GetAbsoluteEnergyScaleError();
    m_abs_eoffset_error = Nom->GetAbsoluteEnergyOffsetError();
    
    m_rel_eoffset_error = Nom->GetRelativeEnergyOffsetError();
    
    ResolutionError = Nom->GetResolutionError();
    ResolutionErrorUncorrelated = Nom->GetResoUncorrelatedError();
    ResolutionRange = 8;// Why 8σ? Seems chosen trivially but in my opinion it's the way to limit the range of the convolution, for a Normal distribution this covers up to 99.99999... of the area
    
    BCW = Nom->GetBCWModel();
    LBNL = Nom->GetLBNLModel();
    Unified = Nom->GetUnifiedModel();
    NominalDetectorEfficiencyDelayed = Nom->GetEnergyDelayedCutDetectorEfficiency();
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        DetectorEfficiencyDelayed[AD] = NominalDetectorEfficiencyDelayed;
        
        for (Int_t week = 0; week <Nweeks; week++)
        {
            MultiMuonEff[AD][week] = Nom->GetMultiMuonEff(AD,week,0);
            FullTime[AD][week] = Nom->GetFullTime(AD,week,0);
        }
        
    }
    
    DetectorEfficiencyRelativeError = Nom->GetDetectorEfficiencyRelativeError();

    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
#ifdef UseVolumes //To use 2 volumes, otherwise 100 cells.
                
                DetectorProtons[AD][idx][idy] = Nom->GetDetectorProtons(AD,idx);
                
#else
                if((idy==0||idy==(YCellLimit-1)||idx>=(XCellLimit-4)))//nH LS
                {
                    if(Analysis)//Hydrogen LS
                    {
                        DetectorProtons[AD][idx][idy] = Nom->GetDetectorProtons(AD,1)/(XCellLimit+YCellLimit+(YCellLimit-2)*4);//52 cells (10+10+(10-2)+(10-2)+(10-2)+(10-2)) Share uniformly the mass
                        
                    }
                    else//nGd only 1 volume
                    {
                        DetectorProtons[AD][idx][idy] = Nom->GetDetectorProtons(AD,0);
                    }
                }
                else//GdLS
                {
                    DetectorProtons[AD][idx][idy] = Nom->GetDetectorProtons(AD,0)/((XCellLimit-4)*(YCellLimit-2));//48 cells, Share uniformly the mass
                }
#endif
                
                for (Int_t week = 0; week <Nweeks; week++)
                {
                    NominalDetectorEfficiency[AD][week][idx][idy] = Nom->GetDetectorEfficiency(AD,week,idx,idy,0);
                    
                    DetectorEfficiency[AD+NADs*week+NADs*Nweeks*idx+NADs*Nweeks*XCellLimit*idy]=NominalDetectorEfficiency[AD][week][idx][idy];
                    
                       std::cout << "\t \t \t Nominal Detector Efficiency in AD" << AD << " ,week: " << week << " Cell: " << idx << "," << idy << " inside OscillationReactor.h: " << DetectorEfficiency[AD+NADs*week+NADs*Nweeks*idx+NADs*Nweeks*XCellLimit*idy] << std::endl;
                }
            }//idy
        }//idx
    }//AD
    

    delete Nom;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                  USEFUL CONSTRUCTOR
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
OscillationReactor :: OscillationReactor(NominalData* Data)
{    
    DataSet = Data->GetDataSet();
    
    rand = new TRandom3(0);
    
    RandomPredictionDirectory = Data->GetToyMCSamplesDirectory();
    NominalPredictionDirectory = Data->GetPredictionDirectory();
    
    hierarchy=Data->GetHierarchy();
    deltam2_32=2.30e-3;//eV2
    deltam2_21=7.50e-5;//eV2
    deltam2_31=deltam2_32+hierarchy*deltam2_21;
    
    s22t13 = Data->GetSin22t13();
    dm2_ee = Data->GetDm2ee();
    s22t12 = Data->GetSin22t12();
    
    //Hard code these variables since the studies are provided in the following energies:
    InitialEnergy = 1.8;
    FinalEnergy = 12;
    InitialVisibleEnergy = 0;
    FinalVisibleEnergy = 12;
    
    TotalBins = MatrixBins;
    BinWidth=(FinalVisibleEnergy-InitialVisibleEnergy)/TotalBins;
    
    n_evis_bins = Data->GetVisibleBins();
    n_etrue_bins = Data->GetTrueBins();

    for (Int_t i = 0; i <= n_evis_bins; i++)
    {
        evis_bins[i] = Data->GetVisibleBinningArray(i);
    }
    
    Nweeks = Data->GetWeeks();
    isH = Data->GetAnalysis();
    
    if(isH)
    {
        AnalysisString = "Hydrogen";
        XCellLimit = VolumeX;
        YCellLimit = VolumeY;
        Data->ReadToyEventsByCell();
    }
    else
    {
        AnalysisString = "Gadolinium";
        XCellLimit = 1;
        YCellLimit = 1;
    }
    
#ifdef UseVolumes
    VolumeString[0] = "Gd-LS";
    VolumeString[1] = "LS";
#endif
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    NADs = Data->GetADs();
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for(Int_t reactor = 0; reactor < NReactors; reactor++)
        {
            ADdistances[AD][reactor] = Data->GetDistances(AD,reactor)/1000;//IN KILOMETERS
        }
        for (Int_t week = 0; week <Nweeks; week++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    IBDEvents[AD][week][idx][idy] = Data->GetIBDEvents(AD,week,idx,idy);
                    ObservedEvents[AD][week][idx][idy] = Data->GetObservedEvents(AD,week,idx,idy);
                }
            }
        }
    }
    
    VaryAccidentalMatrix = Data->GetVaryAccidentalMatrix();
    VaryLiHeMatrix = Data->GetVaryLiHeMatrix();
    VaryFastNeutronsMatrix = Data->GetVaryFastNeutronsMatrix();
    VaryAmCMatrix = Data->GetVaryAmCMatrix();
    DistortLiHeMatrix = Data->GetDistortLiHeMatrix();
    DistortFastNeutronsMatrix = Data->GetDistortFastNeutronsMatrix();
    DistortAmCMatrix = Data->GetDistortAmCMatrix();
    
    IsotopeMatrix = Data->GetIsotopeMatrix();
    ReactorPowerMatrix = Data->GetReactorPowerMatrix();
    Sin22t12Matrix = Data->GetSin22t12Matrix();
    RelativeEnergyScaleMatrix = Data->GetRelativeEnergyScaleMatrix();
    IAVMatrixb = Data->GetIAVMatrix();
    OAVMatrixb = Data->GetOAVMatrix();
    NLMatrix = Data->GetNLMatrix();
    ResolutionMatrix = Data->GetResolutionMatrix();
    EfficiencyMatrix = Data->GetEfficiencyMatrix();
    AllDetectorSystematicsMatrix = Data->GetAllDetectorSystematicsMatrix();
    
    //  IAV Error from Bryce
    IAVNominalError=Data->GetIAVError();
    //  Non uniformity
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        ResolutionBias[AD]=0;
        
        //  Relative energy scale
        m_rel_escale[AD] = Data->GetRelativeEnergyScale();
        m_rel_escale_error[AD] = Data->GetRelativeEnergyError();
        m_rel_escale_nominal[AD] = m_rel_escale[AD];
        m_rel_eoffset[AD] = Data->GetRelativeEnergyOffset();
    }
    m_abs_escale = Data->GetAbsoluteEnergyScale();
    m_abs_eoffset = Data->GetAbsoluteEnergyOffset();
    
    m_abs_escale_nominal = m_abs_escale;
    
    m_abs_escale_error = Data->GetAbsoluteEnergyScaleError();
    m_abs_eoffset_error = Data->GetAbsoluteEnergyOffsetError();
    
    m_rel_eoffset_error = Data->GetRelativeEnergyOffsetError();
    
    ResolutionError = Data->GetResolutionError();
    ResolutionErrorUncorrelated = Data->GetResoUncorrelatedError();
    ResolutionRange = 8;// Why 8σ? Seems chosen trivially but in my opinion it's the way to limit the range of the convolution, for a Normal distribution this covers up to 99.99999... of the area
    
    BCW = Data->GetBCWModel();
    LBNL = Data->GetLBNLModel();
    Unified = Data->GetUnifiedModel();

    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
#ifdef UseVolumes //To use 2 volumes, otherwise 100 cells.
                
                DetectorProtons[AD][idx][idy] = Data->GetDetectorProtons(AD,idx);
#else
                if((idy==0||idy==(YCellLimit-1)||idx>=(XCellLimit-4)))//nH LS
                {
                    if(Analysis)//Hydrogen LS
                    {
                        DetectorProtons[AD][idx][idy] = Data->GetDetectorProtons(AD,1)/(XCellLimit+YCellLimit+(YCellLimit-2)*4);//52 cells (10+10+(10-2)+(10-2)+(10-2)+(10-2)) Share uniformly the mass
                        
                    }
                    else//nGd only 1 volume
                    {
                        DetectorProtons[AD][idx][idy] = Data->GetDetectorProtons(AD,0);
                    }
                }
                else//GdLS
                {
                    DetectorProtons[AD][idx][idy] = Data->GetDetectorProtons(AD,0)/((XCellLimit-4)*(YCellLimit-2));//48 cells, Share uniformly the mass
                }
#endif
                
                for (Int_t week = 0; week <Nweeks; week++)
                {
                    NominalDetectorEfficiency[AD][week][idx][idy] = Data->GetDetectorEfficiency(AD,week,idx,idy,0);
                    
                    DetectorEfficiency[AD+NADs*week+NADs*Nweeks*idx+NADs*Nweeks*XCellLimit*idy]=NominalDetectorEfficiency[AD][week][idx][idy];
                    std::cout << "\t \t \t Nominal Detector Efficiency in AD" << AD << " ,week: " << week << " Cell: " << idx << "," << idy << " inside OscillationReactor.h: " << DetectorEfficiency[AD+NADs*week+NADs*Nweeks*idx+NADs*Nweeks*XCellLimit*idy] << std::endl;

                }
            }//idy
        }//idx
    }//AD
    
    NominalDetectorEfficiencyDelayed = Data->GetEnergyDelayedCutDetectorEfficiency();

    for (Int_t AD = 0; AD<NADs; AD++)
    {
        DetectorEfficiencyDelayed[AD] = NominalDetectorEfficiencyDelayed;
        
        for (Int_t week = 0; week <Nweeks; week++)
        {
            MultiMuonEff[AD][week] = Data->GetMultiMuonEff(AD,week,0);
            FullTime[AD][week] = Data->GetFullTime(AD,week,0);
        }
    }
    
    DetectorEfficiencyRelativeError = Data->GetDetectorEfficiencyRelativeError();
    
}
OscillationReactor :: ~OscillationReactor()
{
    delete rand;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                  GETTERS AND SETTERS
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OscillationReactor :: SetOscillationParameters(Double_t sin22t13,Double_t Dm2_ee)
{
    s22t13=sin22t13;
    dm2_ee=Dm2_ee;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//                      MAINSPECTRUM SHOULD CALL THIS IF OscillationReactor FROM REACTOR MODEL IS USED

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Main script to produce the histograms in each AD after oscillation
void OscillationReactor :: OscillationFromReactorData(Int_t Week,bool mode,bool covMatrix)
{
    Mode = mode;
    
    std::cout <<  "\t ***********************************************************************************************" << std::endl;
    std::cout << "\t \t \t Sin22t12 in OscillationReactor.h is: " << s22t12 << std::endl;
    std::cout << "\t \t \t Sin22t13 in OscillationReactor.h is: " << s22t13 << std::endl;
    std::cout << "\t \t \t Dm2_ee in OscillationReactor.h is: " << dm2_ee << std::endl;
    
    std::cout <<  "\t ***********************************************************************************************" << std::endl;
    
    week = Week;
   
    CovMatrix = covMatrix;
 
    if((Sin22t12Matrix||IsotopeMatrix||ReactorPowerMatrix)&&Mode)
    {
        sprintf(OutputFileName,("./RootOutputs/"+ AnalysisString+ "/RandomOutputs/RandomOscillation_Isotope_%d_Power_%d_Sin22t12_%d_.root").c_str(),IsotopeMatrix,ReactorPowerMatrix,Sin22t12Matrix);
    }
    else
    {
        sprintf(OutputFileName,("./RootOutputs/"+ AnalysisString+ "/NominalOutputs/Oscillation.root").c_str());
    }
    
    LoadReactorHistograms();//Reactor model is the same for all the analyses
    
    LoadNominalBackgrounds();
    
    //Don't load external predictions if the user doesn't choose a file in the nGd case, if he does:
    //Use external near hall predictions:
    if((strcmp((RandomPredictionDirectory).c_str(),"")&&strcmp((NominalPredictionDirectory).c_str(),"")))
    {
        std::cout << " Loading external inputs: " << std::endl;
        
        LoadExternalInputs();
    }
    else//Produce the predictions inside my code:
    {
        GenerateVisibleSpectrum();
    }
    
    TFile* OutputFile = new TFile(OutputFileName, "recreate");
    
    TDirectory* TotalDirectory = OutputFile->mkdir("Total AD Spectra after oscillation");
    TotalDirectory->cd();
    
    for(Int_t AD = 0; AD<NADs; AD++)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                VisibleHisto[AD][idx][idy]->Write(Form("Oscillation Prediction AD%d, week%d, cell %d,%d",AD+1,week,idx,idy));
            }
        }
        if((!strcmp((RandomPredictionDirectory).c_str(),"")&&!strcmp((NominalPredictionDirectory).c_str(),"")))
        {
            //external inputs don't have this histogram:
            TotalOscillatedSpectrumAD[AD]->Write();
        }
    }
    
    if((!strcmp((RandomPredictionDirectory).c_str(),"")&&!strcmp((NominalPredictionDirectory).c_str(),"")))
    {
        TDirectory* EHDirectory = OutputFile->mkdir("EH spectra after oscillation");
        EHDirectory->cd();
        PlotEH();
        
        TDirectory* RatioDirectory = OutputFile->mkdir("EH spectra ratio");
        RatioDirectory->cd();
        PlotRatioEH();
        
        TDirectory* DifferenceDirectory = OutputFile->mkdir("EH spectra difference");
        DifferenceDirectory->cd();
        PlotDifferenceEH();
    }
    delete OutputFile;
}

void OscillationReactor:: LoadNominalBackgrounds()
{
    TH1D* AccidentalsH[MaxDetectors][VolumeX][VolumeY];
    TH1D* LiHeH[MaxDetectors][VolumeX][VolumeY];
    TH1D* FastNeutronsH[MaxDetectors][VolumeX][VolumeY];
    TH1D* AmCH[MaxDetectors][VolumeX][VolumeY];
    std::string BackgroundS;
    
    if(isH)
    {
        BackgroundS = "HBackground";
    }
    else
    {
        BackgroundS = "GDBackground";
    }
    
    TFile* BackgroundsF = new TFile(("./BackgroundSpectrum/"+BackgroundS+"/Backgrounds.root").c_str());
    
    for (Int_t AD =0; AD<NADs; AD++)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                AccidentalsH[AD][idx][idy]=(TH1D*)gDirectory->Get(Form("Accidentals_AD%i_Volume%i",AD,idx));
                LiHeH[AD][idx][idy]=(TH1D*)gDirectory->Get(Form("LiHe_AD%i_Volume%i",AD,idx));
                FastNeutronsH[AD][idx][idy]=(TH1D*)gDirectory->Get(Form("FN_AD%i_Volume%i",AD,idx));
                AmCH[AD][idx][idy]=(TH1D*)gDirectory->Get(Form("AmC_AD%i_Volume%i",AD,idx));
            }
        }
    }
    
    delete BackgroundsF;
    
    //Nominal background spectrum
    for (Int_t AD =0; AD<NADs; AD++)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                BackgroundSpectrumH[AD][idx][idy] = new TH1D(Form("Background Spectrum_AD%i_Volume%i",AD,idx),Form("Background Spectrum_AD%i_Volume%i",AD,idx),n_evis_bins,evis_bins);
                BackgroundSpectrumH[AD][idx][idy]->Add(AccidentalsH[AD][idx][idy]);
                BackgroundSpectrumH[AD][idx][idy]->Add(LiHeH[AD][idx][idy]);
                BackgroundSpectrumH[AD][idx][idy]->Add(FastNeutronsH[AD][idx][idy]);
                BackgroundSpectrumH[AD][idx][idy]->Add(AmCH[AD][idx][idy]);//scaled in livetime and with ad efficiencies included.
                
                delete AccidentalsH[AD][idx][idy];
                delete LiHeH[AD][idx][idy];
                delete FastNeutronsH[AD][idx][idy];
                delete AmCH[AD][idx][idy];

            }
        }
    }
}

void OscillationReactor :: LoadExternalInputs()//This only works for week = 1 (inclusive)
{
    std::cout << " Loading external inputs " << std::endl;
    
    
    for (Int_t AD =0; AD<NADs; AD++)
    {
        //Load predictions: (need to format root file so there are 6/8 nominal histograms named: "h_nominal_ad%i",AD+1")
        
        TFile* LoadPredictionF;
        
        if(Mode)//return random toy sample
        {
            std::cout << " Need to read random tree here, not nominal " << std::endl;
            exit(EXIT_FAILURE);

            LoadPredictionF = new TFile((RandomPredictionDirectory).c_str());
            
            //Get random sample, need to see how the tree is built, this is nominal sample so far.
            if(NADs==8)
            {
                std::cout << " NEED RANDOM 8 AD SAMPLES " << std::endl;
                
                exit(EXIT_FAILURE);
            }
            else
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        VisibleHisto[AD][idx][idy]=(TH1D*)gDirectory->Get(Form("h_nominal_ad%i",AD+1));//For 6ADs (Yasu's inputs)
                    }
                }
            }
            
        }
        else//Nominal prediction
        {
            LoadPredictionF = new TFile((NominalPredictionDirectory).c_str());
            
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    if(NADs==8)
                    {
                        VisibleHisto[AD][idx][idy]=(TH1D*)gDirectory->Get(Form("h_nominal_stage_%i_ad_%i",week+1,AD+1));//For 8ADs (Henoch's inputs)
                    }
                    else
                    {
                        VisibleHisto[AD][idx][idy]=(TH1D*)gDirectory->Get(Form("h_nominal_ad%i",AD+1));//For 6ADs (Yasu's inputs)
                    }
                }
            }
            
        }
        
        delete LoadPredictionF;
        
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                //Already binned, shouldn't be necessary:

                VisibleHisto[AD][idx][idy]=(TH1D*)VisibleHisto[AD][idx][idy]->Rebin(n_evis_bins,Form("Visible Oscillation Prediction AD%d, week%d, cell %d,%d",AD+1,week, idx, idy),evis_bins);

                VisibleHisto[AD][idx][idy]->Scale(ObservedEvents[AD][week][idx][idy]/VisibleHisto[AD][idx][idy]->Integral());//Here events are in days and AD effect corrected, but backgrounds are included!
                
                //Substract Backgrounds in oscillation.h
                //  VisibleHisto[AD][idx][idy]->Add(BackgroundSpectrumH[AD][idx][idy],-1.);//backgrounds are in livetime and with AD efficiencies!
            }
        }
        
        //CORRECT FOR EFFICIENCIES IN OSCILLATION
        
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//                                          ANTINEUTRINO SURVIVAL PROBABILITY

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Calculates oscillation probability
Double_t OscillationReactor :: OscProb(Double_t L, Double_t E, Double_t S22t13,Double_t dm2ee)
{
    L = L*1000; //Distances are given in km, the 1.267 factor works for km/gev or m/mev
    
    Double_t theta13=TMath::ASin(sqrt(S22t13))*0.5;
    Double_t theta12=TMath::ASin(sqrt(s22t12))*0.5;
    Double_t dm231 = dm2ee-hierarchy*(5.21e-5+deltam2_21);
    //    std::cout << "dm231 is: " << dm231 << "calculated from dmee: " << dm2ee << std::endl;
    
    if(!DeltaMee)
    {
        //This calculation follows equation 9 in the TDR
        Double_t term1=s22t12*pow(cos(theta13),4)*pow(sin(1.267*deltam2_21*L/E),2);
        Double_t term2=pow(cos(theta12)*sin(2*theta13)*sin(1.267*dm231*L/E),2);
        Double_t term3=pow(sin(theta12)*sin(2*theta13)*sin(1.267*deltam2_32*L/E),2);
        return 1-(term1+term2+term3);
    }
    else
    {
        //Δm_ee aproximation
        Double_t term1=s22t12*pow(cos(theta13),4)*pow(sin(1.267*deltam2_21*L/E),2);
        Double_t term2=pow(sin(2*theta13)*sin(1.267*dm2ee*L/E),2);
        return 1-(term1+term2);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//                                      SPECTRA AFTER OSCILLATION IN EACH AD IS CALCULATED

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Calculates the total spectrum in each AD after oscillation from the contributions of all reactors.
void OscillationReactor :: SetUpDetectorResponse()
{
    // Apply detector distortions to get visible energy.
    
    #ifdef PrintEps
        TLegend *legend=new TLegend(0.6,0.15,0.88,0.35);
        legend->SetTextFont(72);
        legend->SetTextSize(0.02);
        legend->SetFillColor(0);
        
        TCanvas* VisibleC = new TCanvas("VisibleC","VisibleC");
        GetNeutrinoToPositronFunction(0);
        VisibleF->SetLineWidth(2);
        VisibleF->GetXaxis()->SetTitle("E_{#nu} [MeV]");
        VisibleF->GetXaxis()->SetTitleSize(0.05);
        VisibleF->GetYaxis()->SetTitleSize(0.05);
        VisibleF->GetYaxis()->SetTitle("E_{e^{+}}[MeV]");
        VisibleF->SetLineColor(46);
        legend->AddEntry(VisibleF,"0th order","l");
        VisibleF->Draw();
        GetNeutrinoToPositronFunction(1);
        VisibleF->SetLineWidth(2);
        VisibleF->SetLineColor(9);
        VisibleF->Draw("same");
        legend->AddEntry(VisibleF,"1st order","l");
        legend->Draw("same");
        
        VisibleC->Print(("./Images/"+AnalysisString+"/Detector/VisibleFunctions.eps").c_str(),".eps");
        delete VisibleC;
    #endif
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                           Set kinematic energy shift function
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    GetNeutrinoToPositronFunction(0);//0 for 0th order, 1 for 1st order
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                           Set nominal IAV Matrix
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    LoadIavCorrection();
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                            Set nominal NL parameters
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    SetBCWModel(BCW);
    SetLBNLModel(LBNL);
    SetUnifiedModel(Unified);
    
    LoadNLParameters();//nominal parameters
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                              Set nominal resolution parameters
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ResoF = new TF1("ResoF",this,&OscillationReactor::ResolutionF,0,20,3,"OscillationReactor","ResolutionF");
    //    ResoF->SetParameters(0.022,0.077,0.018); // based on Bryce's TN
    ResoF->SetParameters(0.0148,0.0869,0.0271); // based on Doc-8982
    
    //BCW values http://dayabay.ihep.ac.cn/DocDB/0087/008768/013/6AdAnalysis-BCW.pdf are different! 0.13226034931, 0.32604048828, 0.26196488314
    
    // GET ENERGY SHIFT< IAV< NL < RESOLUTION
    
    if(Mode)
    {
        SetSystematic();
    }
    
    SetNLParameters(Mode);
}

void OscillationReactor :: GenerateVisibleSpectrum()
{
    Char_t histnameTotal[100];//47
    Char_t histnameAD[100];//55
    
    std::cout << " Generating " << AnalysisString << " ToyMC predictions " << std::endl;
    if(!isH)//Gadolinium
    {
        SetUpDetectorResponse();//Need to load it before the random efficiency is applied
    }
    else
    {
        if(EfficiencyMatrix)
        {
            this->RandomEfficiency();//Uncertainty not included in the Toy MC, but here.
        }
    }
    
    LoadNominalTotalOscillatedSpectrum();
    
    //Oscilation probabilities calculated for all baselines
    for (Int_t AD = 0; AD <NADs; AD++)
    {
        if (Nweeks == 1)
        {
            sprintf(histnameTotal,"Total spectrum after oscillation at AD%i", AD+1);
        }
        else
        {
            sprintf(histnameTotal,"Total spectrum after oscillation at AD%i, week%i", AD+1, week+1);
        }
        
        TotalOscillatedSpectrumAD[AD] = new TH1D(histnameTotal, histnameTotal, Nbins, InitialEnergy, FinalVisibleEnergy);
        
        for (Int_t i = 0; i <NReactors; i++)
        {
            if (Nweeks == 1)
            {
                sprintf(histnameAD,"Spectrum after oscillation at AD%i from Reactor%i", AD+1, i+1);
            }
            else
            {
                sprintf(histnameAD,"Spectrum after oscillation at AD%i from Reactor%i, week%i", AD+1, i+1, week+1);
            }
            
            OscillatedSpectrumAD[i] = new TH1D(histnameAD, histnameAD, Nbins, InitialEnergy, FinalVisibleEnergy);
            
            //Oscillation:
            
            Int_t TotalPoints = ReactorSpectrumH[i]->GetXaxis()->GetNbins();
            
//            std::cout << "EVENTS IN REACTOR SPECTRUM BEFORE OSCILLATION : " << ReactorSpectrumH[i]->Integral() << std::endl;

            for(Int_t pts=1;pts<=TotalPoints;++pts)
            {
                Double_t Energy = ReactorSpectrumH[i]->GetXaxis()->GetBinCenter(pts);
                //std::cout << Energy <<"\n";
                OscillatedSpectrumAD[i]->SetBinContent(pts,ReactorSpectrumH[i]->GetBinContent(pts)*OscProb(ADdistances[AD][i],Energy,s22t13,dm2_ee)/(4*TMath::Pi()*ADdistances[AD][i]*ADdistances[AD][i]));
                
                //std::cout << " Survival Probability at energy " << Energy << " is " << OscProb(ADdistances[AD][i],Energy,s22t13,dm2_ee) << std::endl;
                
                //std::cout << " Reactor at energy " << Energy << " is " << ReactorSpectrumH[i]->GetBinContent(pts) << std::endl;

                //                    Double_t Bin = OscillatedSpectrumAD[i]->GetBinContent(pts);
                //                      std::cout << Bin <<"\n";
            }
            
//            std::cout << "EVENTS IN REACTOR OSCILLATED SPECTRUM AFTER OSCILLATION : " << OscillatedSpectrumAD[i]->Integral() << std::endl;

            TotalOscillatedSpectrumAD[AD]->Add(OscillatedSpectrumAD[i]);
            delete OscillatedSpectrumAD[i];
        }
        
        std::cout << "TOTAL REACTOR SPECTRUM IN AD : " << AD << " IS : " << TotalOscillatedSpectrumAD[AD]->Integral() << std::endl;

        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                ScaledOscillatedSpectrumAD[AD][idx][idy] = (TH1D*)TotalOscillatedSpectrumAD[AD]->Clone();
                
                std::cout << "Events in ScaledOscillatedSpectrumAD before : " << ScaledOscillatedSpectrumAD[AD][idx][idy]->Integral() << std::endl;

                std::cout << "EVENTS IN TOTAL OSCILLATED SPECTRUM : " << TotalOscillatedSpectrumAD[AD]->Integral() << std::endl;
                std::cout << "EVENTS IN TOTAL NOMINAL OSCILLATED SPECTRUM : " << IntegralNominalTotalOscillatedSpectrumAD[AD] << std::endl;

                ScaledOscillatedSpectrumAD[AD][idx][idy]->Scale(1./IntegralNominalTotalOscillatedSpectrumAD[AD]);//Divide over the nominal reactor spectrum area. This way we keep the reactor variations information.
                
//                std::cout << "Events in ScaledOscillatedSpectrumAD after : " << ScaledOscillatedSpectrumAD[AD][idx][idy]->Integral() << std::endl;

                Double_t EffLost;
                
                if(isH)
                {
                    EffLost = 1.05;
                }
                else
                {
                    EffLost = 1.0012;//0.9988
                }
                //the 1.05 hard coded is the efficiency loss because of the 1.5 MeV cut. Technically this is true energy, and the IBD events are visible energy with selection cuts included. This scaling would decrease the number of events and make the covariance matrices smaller, so we need to add a small percentage due to that loss. It has been measured to be 1.5 for the nH analysis.
                ScaledOscillatedSpectrumAD[AD][idx][idy]->Scale(EffLost*IBDEvents[AD][week][idx][idy]);//No backgrounds, given full time and AD effects included
//                std::cout << "Events in ScaledOscillatedSpectrumAD ~ IBD Events : " << ScaledOscillatedSpectrumAD[AD][idx][idy]->Integral() << std::endl;
                
//                std::cout << "Events in ScaledOscillatedSpectrumAD IBD EVENTS + EFFICIENCY : " << ScaledOscillatedSpectrumAD[AD][idx][idy]->Integral() << std::endl;

                ScaledOscillatedSpectrumAD[AD][idx][idy]->Scale(DetectorProtons[0][idx][idy]/(DetectorProtons[AD][idx][idy]*MultiMuonEff[AD][week]*FullTime[AD][week]));//in days and km^2 already calculated in CrossSection.//Correct events for efficiency
                
                ScaledOscillatedSpectrumAD[AD][idx][idy]->Scale(DetectorEfficiency[AD+NADs*week+NADs*Nweeks*idx+NADs*Nweeks*XCellLimit*idy]/NominalDetectorEfficiency[AD][week][idx][idy]);//If we vary the efficiency, it changes here, after the previous scaling, which is efficiency corrected.
//
//                std::cout << "SCALE : " << DetectorProtons[AD][idx][idy]*DetectorEfficiency[AD+NADs*week+NADs*Nweeks*idx+NADs*Nweeks*XCellLimit*idy]/MultiMuonEff[AD][week] << std::endl;//Per day, since the events are corrected I will generate the covariance matrices with corrected events (no weekly efficiencies)
//                
//                //Used to scale the predictions in the covariance matrix and get the proper scaling of the matrix according to the expected events in the detector. It would be better not to use it, to allow the covariance matrix to be free of data inputs, to do so we need a good stimation of the events in ScaledOscillatedSpectrum, which right now is off by a few of orders of magnitude (2-3)
//                
                if(EfficiencyMatrix)
                {
                    std::cout << " EFFICIENCY VARIATION : " << DetectorEfficiency[AD+NADs*week+NADs*Nweeks*idx+NADs*Nweeks*XCellLimit*idy]/NominalDetectorEfficiency[AD][week][idx][idy]  << std::endl;
                }
                std::cout << " EVENTS AFTER CROSS SECTION, REACTOR MODEL AND TARGET MASSES for s22t13,dm2_ee: " << s22t13 << dm2_ee << " is" << ScaledOscillatedSpectrumAD[AD][idx][idy]->Integral() << "amd for the nominal oscillation values should be equal to: " << IBDEvents[AD][week][idx][idy]*DetectorProtons[0][idx][idy]/(DetectorProtons[AD][idx][idy]*MultiMuonEff[AD][week]*FullTime[AD][week]) << std::endl;
                
//                std::cout << "BACKGROUND EVENTS : " << BackgroundSpectrumH[AD][idx][idy]->Integral() << ", COMPARE TO EVENTS IN ANTINEUTRINO SPECTRUM TO SEE IF MAKES SENSE " << std::endl;
                
                
                //ScaledOscillatedSpectrumAD[AD][idx][idy]->Scale(IBDEvents[AD][week][idx][idy]/VisibleHisto[AD][idx][idy]->Integral()/ScaledOscillatedSpectrumAD[AD][idx][idy]->Integral());//Scale the true spectrum to the expected corrected value of events.
                
                VisibleHisto[AD][idx][idy] = new TH1D(Form("Visible Oscillation Prediction AD%d, week%d, cell %d,%d",AD+1,week, idx, idy),Form("Visible Oscillation Prediction AD%d, week%d, cell %d,%d",AD+1,week, idx, idy), MatrixBins,0,FinalVisibleEnergy);
            }
        }
    }
    
    #ifdef PrintEps
        TCanvas* AntineutrinoSpectrumC = new TCanvas("AntineutrinoSpectrumC","AntineutrinoSpectrumC");
        
        AntineutrinoSpectrumC->Divide(NADs/2,2);
        
        for (Int_t AD = 0; AD <NADs; AD++)
        {
            AntineutrinoSpectrumC->cd(AD+1);
            TotalOscillatedSpectrumAD[AD]->SetStats(ShowStatBoxInPlots);
            TotalOscillatedSpectrumAD[AD]->GetXaxis()->SetTitle("E_{true} (MeV)");
            TotalOscillatedSpectrumAD[AD]->GetXaxis()->SetTitleSize(0.04);
            TotalOscillatedSpectrumAD[AD]->GetYaxis()->SetTitleSize(0.04);
            TotalOscillatedSpectrumAD[AD]->GetYaxis()->SetTitle("Events");
            TotalOscillatedSpectrumAD[AD]->GetYaxis()->SetTitleOffset(1.5);
            TotalOscillatedSpectrumAD[AD]->SetTitle(Form("AD%d",AD+1));
            
            TotalOscillatedSpectrumAD[AD]->Draw("HIST");
        }
        AntineutrinoSpectrumC->Print(("./Images/"+AnalysisString+"/ReactorSpectrumOscillatedInADs_NotScaled.eps").c_str());
        delete AntineutrinoSpectrumC;
    
    TCanvas* ScaledOscillatedC = new TCanvas("ScaledOscillatedC","ScaledOscillatedC");
    
    ScaledOscillatedC->Divide(NADs/2,2*XCellLimit);
    
    for(Int_t idy=0; idy<YCellLimit; idy++)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for (Int_t AD = 0; AD <NADs; AD++)
            {
                ScaledOscillatedC->cd(AD+idx*NADs+1);
                ScaledOscillatedSpectrumAD[AD][idx][idy]->SetStats(ShowStatBoxInPlots);
                ScaledOscillatedSpectrumAD[AD][idx][idy]->GetXaxis()->SetTitle("E_{true} (MeV)");
                ScaledOscillatedSpectrumAD[AD][idx][idy]->GetXaxis()->SetTitleSize(0.04);
                ScaledOscillatedSpectrumAD[AD][idx][idy]->GetYaxis()->SetTitleSize(0.04);
                ScaledOscillatedSpectrumAD[AD][idx][idy]->GetYaxis()->SetTitle("Events/day");
                ScaledOscillatedSpectrumAD[AD][idx][idy]->SetTitle((VolumeString[idx] + Form(" Volume, AD%d",AD+1)).c_str());
                
                ScaledOscillatedSpectrumAD[AD][idx][idy]->Draw();
            }
        }
    }
    ScaledOscillatedC->Print(("./Images/"+AnalysisString+"/ReactorSpectrumOscillatedInADs_Scaled.eps").c_str());
    delete ScaledOscillatedC;
    
    #endif

    if(!isH)//Gadolinium
    {
        
        //nGd Toy MC based on Yasu's work (but not identical, check if they agree)
        
        for (Int_t AD = 0; AD <NADs; AD++)
        {
            TH1D* ShiftHisto = new TH1D(Form("Shift Oscillation Prediction AD%d, week%d", AD+1,week),Form("Shift  Oscillation Prediction AD%d, week%d", AD+1,week),MatrixBins,0,FinalVisibleEnergy);
            
            TH1D* IAVHisto = new TH1D(Form("IAV Oscillation Prediction AD%d, week%d", AD+1,week),Form("IAV Oscillation Prediction AD%d, week%d", AD+1,week),MatrixBins,0,FinalVisibleEnergy);
            
            TH1D* NLHisto = new TH1D(Form("NL Oscillation Prediction AD%d, week%d", AD+1,week),Form("NL  Oscillation Prediction AD%d, week%d", AD+1,week),MatrixBins,0,FinalVisibleEnergy);
            
            for(Int_t TrueEnergyIndex=0; TrueEnergyIndex<TotalBins; TrueEnergyIndex++)
            {
                GetOscEnergyShift(AD,TrueEnergyIndex,week);
                GetOscIAVShift(AD,TrueEnergyIndex,week);
                GetOscNLShift(AD,TrueEnergyIndex,week);
                GetOscResolutionShift(AD,TrueEnergyIndex,week);
                
                for(Int_t VisibleEnergyIndex=1; VisibleEnergyIndex<=TotalBins; VisibleEnergyIndex++)
                {
                    ShiftHisto->SetBinContent(VisibleEnergyIndex, ShiftHisto->GetBinContent(VisibleEnergyIndex)+OscDeltaPositronSpectrumH[TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex));
                    IAVHisto->SetBinContent(VisibleEnergyIndex, IAVHisto->GetBinContent(VisibleEnergyIndex)+OscDeltaIAVSpectrumH[TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex));
                    NLHisto->SetBinContent(VisibleEnergyIndex, NLHisto->GetBinContent(VisibleEnergyIndex)+OscDeltaNLSpectrumH[TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex));
                    VisibleHisto[AD][0][0]->SetBinContent(VisibleEnergyIndex, VisibleHisto[AD][0][0]->GetBinContent(VisibleEnergyIndex)+OscDeltaVisibleSpectrumH[TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex));
                }
                
                delete OscDeltaPositronSpectrumH[TrueEnergyIndex];
                delete OscDeltaIAVSpectrumH[TrueEnergyIndex];
                delete OscDeltaNLSpectrumH[TrueEnergyIndex];
                delete OscDeltaVisibleSpectrumH[TrueEnergyIndex];
            }
            #ifdef PrintEps
                TCanvas* CheckC = new TCanvas("CheckC","CheckC",1600,1600);
                CheckC->Divide(2,2);
                CheckC->cd(1);
                ShiftHisto->SetStats(ShowStatBoxInPlots);
                ShiftHisto->SetTitle("IBD Kinematics");
                //            ShiftHisto->Scale(IBDEvents[AD][week]/VisibleHisto[AD]->Integral());//Here AD effects are included! CORRECT FOR EFFICIENCIES IN OSCILLATION
                
                ShiftHisto->Draw();
                CheckC->cd(2);
                IAVHisto->SetStats(ShowStatBoxInPlots);
                IAVHisto->SetTitle("IAV correction");
                //            IAVHisto->Scale(IBDEvents[AD][week]/VisibleHisto[AD]->Integral());//Here AD effects are included! CORRECT FOR EFFICIENCIES IN OSCILLATION
                
                IAVHisto->Draw();
                
                CheckC->cd(3);
                NLHisto->SetStats(ShowStatBoxInPlots);
                NLHisto->SetTitle("Non-linearity correction");
                //            NLHisto->Scale(IBDEvents[AD][week]/VisibleHisto[AD]->Integral());//Here AD effects are included! CORRECT FOR EFFICIENCIES IN OSCILLATION
                
                NLHisto->Draw();
                CheckC->cd(4);
                VisibleHisto[AD][0][0]->SetStats(ShowStatBoxInPlots);
                VisibleHisto[AD][0][0]->SetTitle("Visible energy");
                //            VisibleHisto[AD]->Scale(IBDEvents[AD][week]/VisibleHisto[AD]->Integral());//Here AD effects are included! CORRECT FOR EFFICIENCIES IN OSCILLATION
                VisibleHisto[AD][0][0]->Draw();
                CheckC->Update();
                CheckC->Print(("./Images/"+ AnalysisString+ Form("/Detector/SpectrumVisibleStepsAD%d.eps",AD+1)).c_str());
                delete CheckC;
                
                TCanvas* CheckC2 = new TCanvas("CheckC","CheckC");
                
                TLegend *legend2=new TLegend(0.6,0.15,0.88,0.35);
                legend2->SetTextFont(72);
                legend2->SetTextSize(0.02);
                legend2->SetFillColor(0);
                ShiftHisto->SetTitle("All effects");
                ShiftHisto->SetLineColor(1);
                ShiftHisto->GetXaxis()->SetTitle("E_{MeV}");
                ShiftHisto->Draw();
                legend2->AddEntry(ShiftHisto,"Kinematics ","l");
                
                IAVHisto->SetLineColor(2);
                IAVHisto->Draw("same");
                legend2->AddEntry(ShiftHisto,"IAV ","l");
                
                NLHisto->SetLineColor(3);
                NLHisto->Draw("same");
                legend2->AddEntry(ShiftHisto,"NL ","l");
                
                VisibleHisto[AD][0][0]->SetLineColor(4);
                legend2->AddEntry(VisibleHisto[AD][0][0],"Visible ","l");
                VisibleHisto[AD][0][0]->Draw("same");
                legend2->Draw("same");
                CheckC2->Update();
                
                CheckC2->Print(("./Images/"+ AnalysisString+ Form("/Detector/AllVisibleStepsAD%d.eps",AD+1)).c_str());
                //            VisibleHisto[AD]->Scale(VisibleHisto[AD]->Integral()/IBDEvents[AD][week]);//undo scaling
                delete CheckC2;
        #endif
            //
            delete ShiftHisto;
            delete IAVHisto;
            delete NLHisto;
        }
    }
    else//Hydrogen
    {
        //  Here generate Xiang Pan/Logan Toy MC matrix and multiply by the reactor spectrum (his Toy MC is produced from a flat antineutrino spectrum, with no oscillation/reactor information inside the prediction)
        
        TH2D* nHPredictionMatrix[MaxDetectors][VolumeX][VolumeY];
        
        Int_t ShiftBin = Int_t(InitialEnergy*MatrixBins/(FinalEnergy));
        
        Int_t Ybins;
        Int_t Xbins;
        
        for(Int_t AD = 0; AD<NADs; AD++)
        {
#ifdef TestDiagonalMatrix//Test diagonal matrix
                nHPredictionMatrix[AD][0][0] = new TH2D("DiagonalResponseMatrix","DiagonalResponseMatrix",n_evis_bins,evis_bins[0],evis_bins[n_evis_bins],TotalOscillatedSpectrumAD[0]->GetXaxis()->GetNbins(),InitialEnergy,FinalEnergy);
                
                Ybins = nHPredictionMatrix[AD][0][0]->GetYaxis()->GetNbins();
                Xbins = nHPredictionMatrix[AD][0][0]->GetXaxis()->GetNbins();

                for(Int_t i = 1; i<=Ybins; i++)//visible, 240 bins
                {
                    for(Int_t j = 1; j<=Xbins; j++)//true bins are not 240 //1.8 to 12MeV in 0.05 steps
                    {
                        if(i==j)
                        {
                            nHPredictionMatrix[AD][0][0]->SetBinContent(i+1,j+1,1);
                        }
                        else
                        {
                            nHPredictionMatrix[AD][0][0]->SetBinContent(i+1,j+1,0);
                        }
                    }
                }
#else//This is the normal code:
            
            string SystematicS;
            
            if(VaryAccidentalMatrix||VaryLiHeMatrix||VaryFastNeutronsMatrix||VaryAmCMatrix||DistortLiHeMatrix||DistortFastNeutronsMatrix||DistortAmCMatrix||ReactorPowerMatrix||IsotopeMatrix||Sin22t12Matrix||EfficiencyMatrix)
            {
                SystematicS = "NominalResponseMatrix";
            }
            else if(Mode)
            {
                if(RelativeEnergyScaleMatrix)
                {
                    SystematicS = "RelativeEnergyScaleResponseMatrix";
                }
                if(IAVMatrixb)
                {
                    SystematicS = "IAVResponseMatrix";
                }
                if(ResolutionMatrix)
                {
                    SystematicS = "ResolutionResponseMatrix";
                }
//                if(EfficiencyMatrix)
//                {
//                    SystematicS = "EfficiencyResponseMatrix";
//                }
                if(NLMatrix)
                {
                    SystematicS = "NLResponseMatrix";
                }
                if(OAVMatrixb)
                {
                    SystematicS = "OAVResponseMatrix";
                }
                if(AllDetectorSystematicsMatrix)
                {
                    SystematicS = "AllDetectorSystematicsResponseMatrix";
                }
            }
            else
            {
                SystematicS = "NominalResponseMatrix";
            }
            
#ifndef EREC_COMPARISON
            
            TFile* nHMatrixF = new TFile(("./ResponseMatrices/Hydrogen/"+SystematicS+Form("%i_%i.root",n_evis_bins,n_etrue_bins)).c_str(),"read");
#else
            TFile* nHMatrixF = new TFile(("./ResponseMatrices/Hydrogen/E_REC_"+SystematicS+Form("%i_%i.root",n_evis_bins,n_etrue_bins)).c_str(),"read");
#endif
            
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    nHPredictionMatrix[AD][idx][idy] = (TH2D*)nHMatrixF->Get(Form("FineEvisEnu%i,Cell%i,%i",AD+1,idx,idy));
                    
                    Ybins = nHPredictionMatrix[AD][idx][idy]->GetYaxis()->GetNbins();
                    Xbins = nHPredictionMatrix[AD][idx][idy]->GetXaxis()->GetNbins();
                    
                    if(Ybins != VisibleHisto[AD][idx][idy]->GetXaxis()->GetNbins())
                    {
                        std::cout << " RESPONSE MATRIX NEEDS REBINNING " << std::endl;
                        
                        exit(EXIT_FAILURE);
                    }
                    if(Xbins != VisibleHisto[AD][idx][idy]->GetXaxis()->GetNbins())
                    {
                        std::cout << " RESPONSE MATRIX NEEDS REBINNING " << std::endl;
                        
                        exit(EXIT_FAILURE);
                    }
                }
            }
            
            delete nHMatrixF;
            
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    for(Int_t i = 1; i<=Ybins; i++)//visible, 240 bins
                    {
                        for(Int_t j = 1; j<=Xbins; j++)//true bins are not 240 //1.8 to 12MeV in 0.05 steps
                        {
                            
                            VisibleHisto[AD][idx][idy]->SetBinContent(i,VisibleHisto[AD][idx][idy]->GetBinContent(i)+nHPredictionMatrix[AD][idx][idy]->GetBinContent(j+ShiftBin,i)*ScaledOscillatedSpectrumAD[AD][idx][idy]->GetBinContent(j));
                            //ScaledOscillatedSpectrum already has information about the ratio of events between GdLs Ls, but the covariance matrices
                            
                            //Try to minimize the statistical noise, but this will come at a cost of additional uncertainty in this process, just a trial to see if the covariance matrices look smoother:
                            //VisibleHisto[AD][idx][idy]->Smooth(0);
                            
                        }//idy
                    }//idx
                    //std::cout << "Visible bin: " << i << " , true bin: " << j << " - Matrix: " << nHPredictionMatrix[AD]->GetBinContent(j+ShiftBin,i) << " - True Spectrum Bin: " << TotalOscillatedSpectrumAD[AD]->GetBinContent(j) << std::endl;
                    
                    delete  nHPredictionMatrix[AD][idx][idy];
                }//j
                
                //std::cout << "Visible bin: " << i  << " - Visible Spectrum Bin: " << VisibleHisto[AD]->GetBinContent(i) << std::endl;
                
            }//i
#endif
        }//AD
    }
    
    for(Int_t idx=0; idx<XCellLimit; idx++)
    {
        for(Int_t idy=0; idy<YCellLimit; idy++)
        {
            for(Int_t AD = 0; AD<NADs; AD++)
            {
                VisibleHisto[AD][idx][idy]=(TH1D*)VisibleHisto[AD][idx][idy]->Rebin(n_evis_bins,Form("Oscillation Prediction AD%d, week%d, cell %d,%d",AD+1,week,idx,idy),evis_bins);
                
                if(!CovMatrix)//Scale data to expected number of events//Actually doing this here in the covariance matrix cancels the reactor effects.
                {
                    //Need to scale this after the rebin, otherwise some events are lost with the prompt energy cut!
                    VisibleHisto[AD][idx][idy]->Scale(IBDEvents[AD][week][idx][idy]/VisibleHisto[AD][idx][idy]->Integral());
                    //Here AD effects are included! CORRECT FOR EFFICIENCIES IN OSCILLATION
                }
                else
                {
                    VisibleHisto[AD][idx][idy]->Scale(FullTime[AD][week]*MultiMuonEff[AD][week]*DetectorProtons[AD][idx][idy]/DetectorProtons[0][idx][idy]);
                }
            }
        }
    }
    
    #ifdef PrintEps
    TCanvas* PredictionC = new TCanvas("PredictionC","PredictionC",700*(NADs/2),500*(XCellLimit*YCellLimit*2));

    PredictionC->Divide(NADs/2,XCellLimit*YCellLimit*2);
    
    for(Int_t idx=0; idx<XCellLimit; idx++)
    {
        for(Int_t idy=0; idy<YCellLimit; idy++)
        {
            for(Int_t AD = 0; AD<NADs; AD++)
            {
                VisibleHisto[AD][idx][idy]->SetStats(ShowStatBoxInPlots);
                
                PredictionC->cd(AD+idy*NADs+NADs*YCellLimit*idx+1);
                VisibleHisto[AD][idx][idy]->SetStats(ShowStatBoxInPlots);
                VisibleHisto[AD][idx][idy]->SetTitle((Form("Spectrum in AD%d, ",AD+1) + VolumeString[idx]).c_str());

                VisibleHisto[AD][idx][idy]->GetXaxis()->SetTitle("E_{vis} (MeV)");
                VisibleHisto[AD][idx][idy]->GetXaxis()->SetTitleSize(0.04);
                VisibleHisto[AD][idx][idy]->GetYaxis()->SetTitleSize(0.04);
                VisibleHisto[AD][idx][idy]->GetYaxis()->SetTitleOffset(1.5);
                VisibleHisto[AD][idx][idy]->GetYaxis()->SetTitle("Events");
                VisibleHisto[AD][idx][idy]->Draw("HIST");
            }
        }
    }
    
        PredictionC->Print(("./Images/"+AnalysisString+"/OscillationReactorPredictionWithoutBackgrounds.eps").c_str());
        delete PredictionC;
    #endif
    
    //TEST FOR NO OSCILLATION ACTUALLY PRODUCES THE SAME NEUTRINO EVENTS IN ALL DETECTORS BEFORE SCALING:
#ifdef NoOscillation
#ifdef TestAllTheSame
//    for (Int_t AD = 0; AD <NADs; AD++)
//    {
//        for(Int_t idx=0; idx<XCellLimit; idx++)
//        {
//            for(Int_t idy=0; idy<YCellLimit; idy++)
//            {
//                assert(Int_t(1000000*VisibleHisto[AD][idx][idy]-VisibleHisto[AD][idx][idy])<=1);
//            }
//        }
//    }
#endif
#endif
    
    for (Int_t AD = 0; AD <NADs; AD++)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                
                //Add backgrounds:
                
                //Add nominal backgrounds here, this way the predictions in the far hall will carry the background variations when they are substracted there
                VisibleHisto[AD][idx][idy]->Add(BackgroundSpectrumH[AD][idx][idy],1.);//full time and detector effects
                
                if(!Mode)//If there are not background variations the observed events should be equal to IBD events + background events
                {
                    if(!CovMatrix)
                    {
                    std::cout << "ASSERTING: " << VisibleHisto[AD][idx][idy]->Integral() << " - " << BackgroundSpectrumH[AD][idx][idy]->Integral() << " = " << IBDEvents[AD][week][idx][idy] << std::endl;
                    }
                    std::cout << "ASSERTING: " << BackgroundSpectrumH[AD][idx][idy]->Integral() << " + " << IBDEvents[AD][week][idx][idy] << " = " << BackgroundSpectrumH[AD][idx][idy]->Integral() + IBDEvents[AD][week][idx][idy] << "or = " << ObservedEvents[AD][week][idx][idy] << std::endl;
                    
                    std::cout << "DIFFERENCE IS : " << Int_t(100000*(BackgroundSpectrumH[AD][idx][idy]->Integral()+IBDEvents[AD][week][idx][idy] - ObservedEvents[AD][week][idx][idy])) << std::endl;
                    
//                    assert(Int_t(100000*(BackgroundSpectrumH[AD][idx][idy]->Integral()+IBDEvents[AD][week][idx][idy] - ObservedEvents[AD][week][idx][idy]))<=1);
                }
                
//                //Then scale to expected number of events.
//                if(!CovMatrix)
//                {
//                    VisibleHisto[AD][idx][idy]->Scale(ObservedEvents[AD][week][idx][idy]/VisibleHisto[AD][idx][idy]->Integral());
//                }
                
                delete BackgroundSpectrumH[AD][idx][idy];
            }
        }
    }
    
#ifdef PrintEps
    TCanvas* PredictionC2 = new TCanvas("PredictionC","PredictionC");
    
    PredictionC2->Divide(NADs/2,XCellLimit*YCellLimit*2);
    
    for(Int_t idx=0; idx<XCellLimit; idx++)
    {
        for(Int_t idy=0; idy<YCellLimit; idy++)
        {
            for(Int_t AD = 0; AD<NADs; AD++)
            {
                VisibleHisto[AD][idx][idy]->SetStats(ShowStatBoxInPlots);
                
                PredictionC2->cd(AD+idy*NADs+NADs*YCellLimit*idx+1);
                
                VisibleHisto[AD][idx][idy]->SetTitle((Form("Spectrum in AD%d, ",AD+1) + VolumeString[idx]).c_str());
                VisibleHisto[AD][idx][idy]->SetStats(ShowStatBoxInPlots);
                VisibleHisto[AD][idx][idy]->GetXaxis()->SetTitle("E_{vis} (MeV)");
                VisibleHisto[AD][idx][idy]->GetXaxis()->SetTitleSize(0.04);
                VisibleHisto[AD][idx][idy]->GetYaxis()->SetTitleSize(0.04);
                VisibleHisto[AD][idx][idy]->GetYaxis()->SetTitleOffset(1.5);
                VisibleHisto[AD][idx][idy]->GetYaxis()->SetTitle("Events");
                VisibleHisto[AD][idx][idy]->Draw("HIST");
            }
        }
    }
    
    PredictionC2->Print(("./Images/"+AnalysisString+"/OscillationReactorPredictionWithBackgrounds.eps").c_str());
    delete PredictionC2;
#endif
    
    TFile* PredictionsFromReactorF = new TFile(("./RootOutputs/"+AnalysisString+"/NominalOutputs/PredictionsFromReactor.root").c_str(),"recreate");
    
    for(Int_t AD = 0; AD<NADs; AD++)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                VisibleHisto[AD][idx][idy]->Write();
            }
        }
    }
    
    delete PredictionsFromReactorF;
}


TH1D* OscillationReactor :: GetReactorOscillatedADSpectrum(Int_t AD,Int_t week,Int_t idx, Int_t idy)
{
    return VisibleHisto[AD][idx][idy];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                          LOAD REACTOR DATA
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OscillationReactor :: LoadReactorHistograms()
{
    if((IsotopeMatrix||ReactorPowerMatrix)&&Mode)
    {
        sprintf(ReactorData,("./RootOutputs/"+AnalysisString+"/RandomOutputs/AntineutrinoSpectrum_Isotope_%d_Power_%d.root").c_str(),IsotopeMatrix,ReactorPowerMatrix);
    }
    else
    {
        sprintf(ReactorData,("./RootOutputs/"+AnalysisString+"/NominalOutputs/AntineutrinoSpectrum.root").c_str());
    }
    
    TFile* ReactorDataF = new TFile(ReactorData);
    
    Char_t filenameReactor[100];
    
    for(Int_t reactor = 0; reactor < NReactors; reactor++)
    {
        if(Nweeks==1)
        {
            sprintf(filenameReactor,"AntineutrinoSpectrumFromReactor%i",reactor+1);
        }
        else
        {
            sprintf(filenameReactor,"Week%i/%i", week, reactor);//Need to change this to accept weekly data properly
        }
        ReactorSpectrumH[reactor] = (TH1D*)gDirectory->Get(filenameReactor);
    }
    
    delete ReactorDataF;
    
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//                                                  PLOT DATA (EH, RATIO, DIFFERENCE)

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Saves all information in the root file
void OscillationReactor :: PlotEH()
{
    //Produce EH1, EH2 and EH3 plots:
    
    Char_t histnameEH1[100];//41
    Char_t histnameEH2[100];
    Char_t histnameEH3[100];
    
    if (Nweeks == 1)
    {
        sprintf(histnameEH1,"Spectrum after oscillation at EH1");
        sprintf(histnameEH2,"Spectrum after oscillation at EH2");
        sprintf(histnameEH3,"Spectrum after oscillation at EH3");
    }
    else
    {
        sprintf(histnameEH1,"Spectrum after oscillation at EH1, week%i", week+1);
        sprintf(histnameEH2,"Spectrum after oscillation at EH2, week%i", week+1);
        sprintf(histnameEH3,"Spectrum after oscillation at EH3, week%i", week+1);
        
    }
    
    OscillatedSpectrumEH1 = (TH1D*)TotalOscillatedSpectrumAD[0]->Clone(histnameEH1);
    OscillatedSpectrumEH1->SetTitle(histnameEH1);
    OscillatedSpectrumEH1->Add(TotalOscillatedSpectrumAD[1]);
    //        OscillatedSpectrumEH1->Scale(1./ADsEH1);
    
    OscillatedSpectrumEH2 = (TH1D*)TotalOscillatedSpectrumAD[2]->Clone(histnameEH2);
    OscillatedSpectrumEH2->SetTitle(histnameEH2);
    //        OscillatedSpectrumEH2->Scale(1./ADsEH2);
    
    OscillatedSpectrumEH3 = (TH1D*)TotalOscillatedSpectrumAD[3]->Clone(histnameEH3);
    OscillatedSpectrumEH3->SetTitle(histnameEH3);
    OscillatedSpectrumEH3->Add(TotalOscillatedSpectrumAD[4]);
    OscillatedSpectrumEH3->Add(TotalOscillatedSpectrumAD[5]);
    //        OscillatedSpectrumEH3->Scale(1./ADsEH3);
    
    OscillatedSpectrumEH1->Write();
    OscillatedSpectrumEH2->Write();
    OscillatedSpectrumEH3->Write();
}

void OscillationReactor :: PlotRatioEH()
{
    
    Char_t histnameRatioEH13[100];//21
    Char_t histnameRatioEH23[100];
    Char_t histnameRatioEH12[100];
    
    if (Nweeks == 1)
    {
        sprintf(histnameRatioEH13,"Ratio EH1 EH3");
        sprintf(histnameRatioEH23,"Ratio EH2 EH3");
        sprintf(histnameRatioEH12,"Ratio EH1 EH2");
    }
    else
    {
        sprintf(histnameRatioEH13,"Ratio EH1 EH3, week%i", week+1);
        sprintf(histnameRatioEH23,"Ratio EH2 EH3, week%i", week+1);
        sprintf(histnameRatioEH12,"Ratio EH1 EH2, week%i", week+1);
        
    }
    //Plot ratio of halls:
    RatioEH1EH3 = (TH1D*)OscillatedSpectrumEH1->Clone(histnameRatioEH13);
    RatioEH1EH3->SetTitle(histnameRatioEH13);
    RatioEH1EH3->GetYaxis()->SetTitle("EH1/EH3");
    RatioEH1EH3->Divide(OscillatedSpectrumEH3);
    
    RatioEH2EH3 = (TH1D*)OscillatedSpectrumEH2->Clone(histnameRatioEH23);
    RatioEH2EH3->SetTitle(histnameRatioEH23);
    RatioEH2EH3->GetYaxis()->SetTitle("EH2/EH3");
    RatioEH2EH3->Divide(OscillatedSpectrumEH3);
    
    RatioEH1EH2 = (TH1D*)OscillatedSpectrumEH1->Clone(histnameRatioEH12);
    RatioEH1EH2->SetTitle(histnameRatioEH12);
    RatioEH1EH2->GetYaxis()->SetTitle("EH1/EH2");
    RatioEH1EH2->Divide(OscillatedSpectrumEH2);
    
    RatioEH1EH3->Write();
    RatioEH2EH3->Write();
    RatioEH1EH2->Write();
}

void OscillationReactor :: PlotDifferenceEH()
{
    
    Char_t histnameDiffEH13[100];//28
    Char_t histnameDiffEH23[100];
    Char_t histnameDiffEH12[100];
    
    if (Nweeks == 1)
    {
        sprintf(histnameDiffEH13,"Difference EH1 EH3");
        sprintf(histnameDiffEH23,"Difference EH2 EH3");
        sprintf(histnameDiffEH12,"Difference EH1 EH2");
    }
    else
    {
        sprintf(histnameDiffEH13,"Difference EH1 EH3, week%i", week+1);
        sprintf(histnameDiffEH23,"Difference EH2 EH3, week%i", week+1);
        sprintf(histnameDiffEH12,"Difference EH1 EH2, week%i", week+1);
    }
    //Plot difference of halls:
    DifferenceEH1EH3 = (TH1D*)OscillatedSpectrumEH1->Clone(histnameDiffEH13);
    DifferenceEH1EH3->SetTitle(histnameDiffEH13);
    DifferenceEH1EH3->GetYaxis()->SetTitle("EH1 - EH3");
    DifferenceEH1EH3->Add(OscillatedSpectrumEH3,-1);
    
    DifferenceEH2EH3 = (TH1D*)OscillatedSpectrumEH2->Clone(histnameDiffEH23);
    DifferenceEH2EH3->SetTitle(histnameDiffEH23);
    DifferenceEH2EH3->GetYaxis()->SetTitle("EH2 - EH3");
    DifferenceEH2EH3->Add(OscillatedSpectrumEH3,-1);
    
    DifferenceEH1EH2 = (TH1D*)OscillatedSpectrumEH1->Clone(histnameDiffEH12);
    DifferenceEH1EH2->SetTitle(histnameDiffEH12);
    DifferenceEH1EH2->GetYaxis()->SetTitle("EH1 - EH2");
    DifferenceEH1EH2->Add(OscillatedSpectrumEH2,-1);
    
    DifferenceEH1EH3->Write();
    DifferenceEH2EH3->Write();
    DifferenceEH1EH2->Write();
}

void OscillationReactor :: FreeMemory()
{
    delete DifferenceEH1EH3;
    delete DifferenceEH2EH3;
    delete DifferenceEH1EH2;
    
    delete RatioEH1EH3;
    delete RatioEH2EH3;
    delete RatioEH1EH2;
    
    delete OscillatedSpectrumEH1;
    delete OscillatedSpectrumEH2;
    delete OscillatedSpectrumEH3;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                           Clean up the dust
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(!isH)
    {
        delete ResoF;
        delete NLF;
        delete VisibleF;
    }
    //Oscilation probabilities calculated for all baselines
    for (Int_t AD = 0; AD <NADs; AD++)
    {
        delete TotalOscillatedSpectrumAD[AD];
        delete NominalTotalOscillatedSpectrumAD[AD];
        
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                delete ScaledOscillatedSpectrumAD[AD][idx][idy];
                delete VisibleHisto[AD][idx][idy];
            }
        }
    }
    
    for (Int_t i = 0; i <NReactors; i++)
    {
        delete ReactorSpectrumH[i];
    }
    
}
void OscillationReactor :: SetSin22t12(Double_t S22t12)
{
    s22t12=S22t12;
}

void OscillationReactor :: SetBCWModel(bool bcw)
{
    if(bcw)
    {
        NLModelE=BCWE;
        NLF = new TF1("NLF",this,&OscillationReactor::NLBCWF,InitialVisibleEnergy,FinalVisibleEnergy,7,"OscillationReactor","NLBCWF");
    }
}

void OscillationReactor :: SetLBNLModel(bool lbnl)
{
    if(lbnl)
    {
        NLModelE=LBNLE;
        NLF = new TF1("NLF",this,&OscillationReactor::NLLBNLF,InitialVisibleEnergy,FinalVisibleEnergy,5,"OscillationReactor","NLLBNLF");
    }
}

void OscillationReactor :: SetUnifiedModel(bool unified)
{
    if(unified)
    {
        NLModelE=UnifiedE;
        NLF = new TF1("NLF",this,&OscillationReactor::NLUnifiedF,InitialVisibleEnergy,FinalVisibleEnergy,m_num_unified_nl_pars+2,"OscillationReactor","NLUnifiedF");
    }
}

TF1* OscillationReactor :: GetNeutrinoToPositronFunction(bool order)
{
    if(order==0)
    {
        Double_t Correction = Mn-Mp-Me;
        VisibleF = new TF1("VisibleF",this,&OscillationReactor::VisibleEnergy0F,InitialVisibleEnergy,FinalVisibleEnergy,1,"OscillationReactor","VisibleEnergy0F");
        VisibleF->SetParameter(0,Correction);
        std::cout << "\t \t \t Zeroth order" << std::endl;
    }
    if(order==1)//review formula, it seems to produce a large peak around 4MeV
    {
        VisibleF = new TF1("VisibleF",this,&OscillationReactor::VisibleEnergy1F,InitialVisibleEnergy,FinalVisibleEnergy,1,"OscillationReactor","VisibleEnergy1F");
        VisibleF->SetParameter(0,-2);//Maybe I can improve this using a MC simulation of the angular distribution calculated in http://authors.library.caltech.edu/2796/1/VOGprd99.pdf
        std::cout << "\t \t \t First order" << std::endl;
    }
    
    TFile* VisibleFile = new TFile("VisibleFunction.root","recreate");
    VisibleF->Write();
    delete VisibleFile;
    
    for(Int_t i = 0; i < MatrixBins; i++)
    {
        EnergyVector[i] = (0.5+i)*BinWidth;//set energy vector for later interpolation purposes
    }
    
    return VisibleF;
}

void OscillationReactor :: RandomIAVMatrix()
{
    for (Int_t AD=0; AD<NADs; AD++)
    {
        rand->SetSeed(0);
        IAVError[AD] = (1+IAVNominalError*rand->Gaus(0,1)); //Each AD is varied individually.
        
        for(Int_t i=0; i < MatrixBins; i++)
        {
            IAVMatrixFrac[AD][i]=IAVError[AD]*NominalIAVMatrixFrac[i];
            IAVMatrix[AD][i][i]= 1-IAVMatrixFrac[AD][i];//Diagonal terms
            
            for(Int_t j=0; j < MatrixBins; j++)
            {
                if(i!=j)
                {
                    IAVMatrix[AD][i][j]=IAVError[AD]*NominalIAVMatrix[i][j];
                }
            }
        }
        
        std::cout << "\t \t \t IAV error" << IAVError[AD] << std::endl;
    }
    #ifdef PrintEps
        TH2D*  NominalIAVMatrixH = new TH2D("Nominal IAV Matrix","Nominal IAV Matrix", MatrixBins,0,MatrixBins,MatrixBins,0,MatrixBins);
        NominalIAVMatrixH->Reset();
        
        for(Int_t i=0; i < MatrixBins; i++)
        {
            for(Int_t j=0; j < MatrixBins; j++)
            {
                NominalIAVMatrixH->SetBinContent(i+1,j+1,NominalIAVMatrix[i][j]);
            }
        }
        
        //        TFile* SaveIAVMatrixF = new TFile("IAVNominalMatrix.root","recreate");
        //
        //        NominalIAVMatrixH->Write();
        //
        //        delete SaveIAVMatrixF;
        
        TCanvas* IAVC = new TCanvas("IAVC","IAVC");
        IAVC->SetLogz();
        
        NominalIAVMatrixH->SetStats(ShowStatBoxInPlots);
        NominalIAVMatrixH->SetTitle("");
        NominalIAVMatrixH->Draw("colz");
        
        IAVC->Print(("./Images/"+ AnalysisString+ "/Detector/NominalIAVMatrix.eps").c_str(),".eps");
        delete IAVC;
        delete NominalIAVMatrixH;
    #endif
}

void OscillationReactor :: RandomAbsoluteEnergyScaleMatrix()
{
    rand->SetSeed(0);
    m_abs_escale = m_abs_escale_nominal + m_abs_escale_error * rand->Gaus(0,1);
    
    std::cout << "\t \t \t Absolute scale" << m_abs_escale << std::endl;
}

void OscillationReactor :: RandomAbsoluteEnergyOffsetMatrix()
{
    rand->SetSeed(0);
    m_abs_eoffset =  m_abs_eoffset_error * rand->Gaus(0,1);
    
    std::cout << "\t \t \t Absolute offset " << m_abs_eoffset << std::endl;
}

void OscillationReactor :: RandomRelativeEnergyOffsetMatrix()
{
    for(Int_t AD=0;AD<NADs;AD++)
    {
        rand->SetSeed(0);
        m_rel_eoffset[AD] =  m_rel_eoffset_error * rand->Gaus(0,1);
        std::cout << "\t \t \t Relative Offset " << m_rel_eoffset[AD] << std::endl;
    }
}

Double_t* OscillationReactor :: GetDelayedEfficiency()
{
    return DetectorEfficiencyDelayed;
}

void OscillationReactor:: RandomEfficiency()
{
    TRandom3* rand = new TRandom3(0);
    
    if(isH)
    {
        for (Int_t AD = 0; AD<NADs; AD++)
        {
            for (Int_t week = 0; week < Nweeks; week++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        rand->SetSeed(0);
                        DetectorEfficiency[AD+NADs*week+NADs*Nweeks*idx+NADs*Nweeks*XCellLimit*idy] = NominalDetectorEfficiency[AD][week][idx][idy]*(1 + DetectorEfficiencyRelativeError * rand->Gaus(0,1));
                        
                        std::cout << "\t \t \t Random Detector Efficiency in AD" << AD << " ,week: " << week << " Cell: " << idx << "," << idy << " inside OscillationReactor.h: " << DetectorEfficiency[AD+NADs*week+NADs*Nweeks*idx+NADs*Nweeks*XCellLimit*idy] << std::endl;
                    }
                }
            }
        }
    }
    else
    {
        for (Int_t AD = 0; AD<NADs; AD++)
        {
            for (Int_t week = 0; week < Nweeks; week++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        rand->SetSeed(0);
                        DetectorEfficiency[AD+NADs*week+NADs*Nweeks*idx+NADs*Nweeks*XCellLimit*idy]= (1 + DetectorEfficiencyRelativeError * rand->Gaus(0,1))
                        * DetectorEfficiencyDelayed[AD]
                        * NominalDetectorEfficiency[AD][week][idx][idy]/NominalDetectorEfficiencyDelayed;//I divide over the nominal detector efficiency delayed because it is included in the total NominalDetectorEfficiency magnitude.
                        
                        std::cout << "\t \t \t Random Detector Efficiency in AD" << AD << " ,week: " << week << " Cell: " << idx << "," << idy << " inside OscillationReactor.h: " << DetectorEfficiency[AD+NADs*week+NADs*Nweeks*idx+NADs*Nweeks*XCellLimit*idy] << std::endl;
                    }
                }
            }
        }
    }
    
    delete rand;
}

Double_t* OscillationReactor :: GetEfficiency()
{
    return DetectorEfficiency;
}

void OscillationReactor :: RandomRelativeEnergyScaleMatrix()
{
    for(Int_t AD=0;AD<NADs;AD++)
    {
        rand->SetSeed(0);
        Double_t rel_escale_shift = m_rel_escale_error[AD] * rand->Gaus(0,1);
        m_rel_escale[AD] = m_rel_escale_nominal[AD] + rel_escale_shift;
        std::cout << "\t \t \t Relative Scale" << m_rel_escale[AD] << std::endl;
        
        DetectorEfficiencyDelayed[AD] = NominalDetectorEfficiencyDelayed + rel_escale_shift * 0.24; // 0.24 is a magic factor that convert energy scale shift to delayed energy cut efficiency
    }
}

void OscillationReactor :: RandomResolutionMatrix()
{
    rand->SetSeed(0);
    Double_t corr_bias = ResolutionError * rand->Gaus(0,1);
    
    for(Int_t AD=0;AD<NADs;AD++)
    {
        rand->SetSeed(0);
        ResolutionBias[AD] = corr_bias + ResolutionErrorUncorrelated * rand->Gaus(0,1);
        std::cout << "\t \t \t Resolution Bias " << ResolutionBias[AD] << std::endl;
    }
}

void OscillationReactor :: Interpolation(TF1* func)
{
    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
    {
        Energy[VisibleEnergyIndex]=(func->GetX(EnergyVector[VisibleEnergyIndex]));//GetX
        
        EnergyIdx[VisibleEnergyIndex]=((Int_t)(Energy[VisibleEnergyIndex]/BinWidth));
        
        //        std::cout << "IAV Energy " << IAVEnergy[VisibleEnergyIndex] << std::endl;
        //        std::cout << "IAV Energy Idx " << IAVEnergyIdx[VisibleEnergyIndex] << std::endl;
        //        std::cout << " VisibleEnergyIndex" << VisibleEnergyIndex << std::endl;
        
        binScaling[VisibleEnergyIndex]=(func->Derivative(Energy[VisibleEnergyIndex]));
        
        //        std::cout << "BIN SCALING" << binScaling[VisibleEnergyIndex] << std::endl;
        
        sign[VisibleEnergyIndex]=(1);
        
        dEtrue[VisibleEnergyIndex]=(Energy[VisibleEnergyIndex] - EnergyIdx[VisibleEnergyIndex]*BinWidth);
        
        //        std::cout << "dEtrue" << dEtrue[VisibleEnergyIndex] << std::endl;
        if (dEtrue[VisibleEnergyIndex] < 0)
        {
            sign[VisibleEnergyIndex]=(-1);
            
        }
        //        std::cout << "sign" << sign[VisibleEnergyIndex] << std::endl;
    }
}

void OscillationReactor :: GetOscEnergyShift(Int_t AD, Int_t TrueEnergyIndex, Int_t week)
{
    Double_t Limit = 1.8;
    
    PositronEnergy = EnergyVector[TrueEnergyIndex];
    
    Double_t e_nu = VisibleF->GetX(PositronEnergy);
    
    Int_t e_nu_Idx = (Int_t)((e_nu)/BinWidth);//TRUE BINWIDTH
    
    OscDeltaPositronSpectrumH[TrueEnergyIndex] = new TH1D(Form("PositronTrueSpectrumAD%d,W%d,%d",AD,week,e_nu_Idx),Form("PositronTrueSpectrumAD%d,W%d,%d",AD,week,e_nu_Idx), TotalBins,InitialVisibleEnergy,FinalVisibleEnergy);
    
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
    
    if(e_nu_Idx>=TotalBins)
    {
        OscDeltaPositronSpectrumH[TrueEnergyIndex]->SetBinContent(TrueEnergyIndex+1, 0);
        return;
    }
    
    Double_t binScaling = VisibleF->Derivative(e_nu,0,0.0001);
    Double_t dNdE = 0;
    Int_t sign = 1;
    Double_t dE = e_nu - (0.5 + e_nu_Idx)*BinWidth;//use energy in the center of the bin for the visible and nl functions
    
    if (dE < 0)
    {
        sign = -1;
    }
    
    if(e_nu_Idx==MatrixBins-1)
    {
        dNdE = TotalOscillatedSpectrumAD[AD]->GetBinContent(e_nu_Idx-Int_t(Limit/BinWidth)+1);
    }
    else
    {
        dNdE = (TMath::Abs(dE)/BinWidth) * ScaledOscillatedSpectrumAD[AD][0][0]->GetBinContent(e_nu_Idx-Int_t(Limit/BinWidth)+sign+1)
        +(1 - TMath::Abs(dE)/BinWidth) * ScaledOscillatedSpectrumAD[AD][0][0]->GetBinContent(e_nu_Idx-Int_t(Limit/BinWidth)+1);
    }
    
    OscDeltaPositronSpectrumH[TrueEnergyIndex]->SetBinContent(TrueEnergyIndex+1,dNdE/binScaling);
    
}

void OscillationReactor :: GetOscIAVShift(Int_t AD, Int_t TrueEnergyIndex, Int_t week)
{
    //IAV
    OscDeltaIAVSpectrumH[TrueEnergyIndex] = new TH1D(Form("IAVSpectrumAD%d,W%d,%d",AD,week,TrueEnergyIndex),Form("IAVSpectrumAD%d,W%d,%d",AD,week,TrueEnergyIndex), TotalBins,InitialVisibleEnergy,FinalVisibleEnergy);
    
    //Calculate IAV Shift
    if (PositronEnergy >= 1.022)//2*0.511
    {
        Int_t IAVLimit;
        
        IAVLimit = TrueEnergyIndex;
        
        for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<IAVLimit+1; VisibleEnergyIndex++)
        {
            OscDeltaIAVSpectrumH[TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, OscDeltaIAVSpectrumH[TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1) + IAVMatrix[AD][TrueEnergyIndex][VisibleEnergyIndex] * OscDeltaPositronSpectrumH[TrueEnergyIndex]->GetBinContent(TrueEnergyIndex+1));
        }
    }
}

void OscillationReactor :: GetOscNLShift(Int_t AD, Int_t TrueEnergyIndex, Int_t week)
{
    //NL
    OscDeltaNLSpectrumH[TrueEnergyIndex] = new TH1D(Form("NLSpectrumAD%d,W%d,%d",AD,week,TrueEnergyIndex),Form("NLSpectrumAD%d,W%d,%d",AD,week,TrueEnergyIndex), TotalBins,InitialVisibleEnergy,FinalVisibleEnergy);
    
    //Calculate NL Shift
    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
    {
        if(EnergyIdx[VisibleEnergyIndex]==(TotalBins-1))
        {
            dNdE[VisibleEnergyIndex]=(OscDeltaIAVSpectrumH[TrueEnergyIndex]->GetBinContent(EnergyIdx[VisibleEnergyIndex]+1));
        }
        else
        {
            dNdE[VisibleEnergyIndex] = ((TMath::Abs(dEtrue[VisibleEnergyIndex])/BinWidth)*OscDeltaIAVSpectrumH[TrueEnergyIndex]->GetBinContent(EnergyIdx[VisibleEnergyIndex]+sign[VisibleEnergyIndex]+1)+(1 - TMath::Abs(dEtrue[VisibleEnergyIndex])/BinWidth)* OscDeltaIAVSpectrumH[TrueEnergyIndex]->GetBinContent(EnergyIdx[VisibleEnergyIndex]+1));
        }
    }
    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
    {
        OscDeltaNLSpectrumH[TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, dNdE[VisibleEnergyIndex]/binScaling[VisibleEnergyIndex]);
    }
}

void OscillationReactor :: GetOscResolutionShift(Int_t AD, Int_t TrueEnergyIndex, Int_t week)
{
    //Resolution
    OscDeltaVisibleSpectrumH[TrueEnergyIndex] = new TH1D(Form("VisibleSpectrumAD%d,W%d,%d",AD,week,TrueEnergyIndex),Form("VisibleSpectrumAD%d,W%d,%d",AD,week,TrueEnergyIndex), TotalBins,InitialVisibleEnergy,FinalVisibleEnergy);
    
    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
    {
        OscDeltaVisibleSpectrumH[TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, OscDeltaNLSpectrumH[TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1));//Copy histogram, I don't use clone() because of memory leaks.
    }
    
    //Calculate Resolution effect
    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
    {
        Double_t sigma;
        
        sigma = (ResoF->Eval(VisibleEnergyIndex*BinWidth) + ResolutionBias[AD]) * VisibleEnergyIndex*BinWidth;//Sigma = FReso(Energy) * Energy
        
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
            
            OscDeltaVisibleSpectrumH[TrueEnergyIndex]->SetBinContent(detIdx+1, (OscDeltaVisibleSpectrumH[TrueEnergyIndex]->GetBinContent(detIdx+1)+OscDeltaNLSpectrumH[TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1)*gausFactor));
        }
    }
    // Scale reso to original NL energy.
    if ((OscDeltaVisibleSpectrumH[TrueEnergyIndex]->Integral())!=0)
    {
        OscDeltaVisibleSpectrumH[TrueEnergyIndex]->Scale(OscDeltaNLSpectrumH[TrueEnergyIndex]->Integral()/OscDeltaVisibleSpectrumH[TrueEnergyIndex]->Integral());
    }
}

void OscillationReactor :: LoadNLParameters()
{
    switch (NLModelE)
    {
        case 0:
        {
            std::cout << "\t \t \t LOAD BCW NL PARAMETERS" << std::endl;
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
            
            TFile *bcw_ele_err_file = new TFile("./Inputs/bcw_nl_data/ele_err.root");
            g_bcw_elec_nl_error[0] = (TGraph*)bcw_ele_err_file->Get("g_up")->Clone();
            g_bcw_elec_nl_error[1] = (TGraph*)bcw_ele_err_file->Get("g_down")->Clone();
            delete bcw_ele_err_file;
        }
            break;
            
        case 1://LBNL NL Model
        {
            std::cout << "\t \t \t LOAD LBNL NL PARAMETERS" << std::endl;
            ifstream lbnl_positron_data("./Inputs/lbnl_nl_data/lbnl_positron_nl.txt");
            
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
            std::cout << "\t \t \t LOAD UNIFIED NL PARAMETERS" << std::endl;
            
            for(Int_t i = 0; i < m_num_unified_nl_pars; i++)
            {
                m_unified_nl_par[i] = 0.0;
                m_unified_nl_par_nominal[i] = 0.0;
                m_unified_nl_par_error[i] = 1.0;
            }
            
            TFile *unified_nl_file = new TFile("./Inputs/unified_nl_data/nl_models_final.root");
            
            // Use IHEP I as the nominal model
            g_unified_positron_nl = (TGraph*)unified_nl_file->Get("positron_0")->Clone();
            g_unified_positron_nl_pulls[0] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",1))->Clone();
            g_unified_positron_nl_pulls[1] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",2))->Clone();
            g_unified_positron_nl_pulls[2] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",3))->Clone();
            g_unified_positron_nl_pulls[3] =  (TGraph*)unified_nl_file->Get(Form("positron_%d",4))->Clone();
            
            delete unified_nl_file;
            
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
            
            #ifdef PrintEps
                TCanvas* NLCurvesC = new TCanvas("NLCurvesC","NLCurvesC");
                
                TLegend *legend3=new TLegend(0.6,0.15,0.88,0.35);
                legend3->SetTextFont(72);
                legend3->SetTextSize(0.02);
                legend3->SetFillColor(0);
                
                g_unified_positron_nl->Draw();
                legend3->AddEntry(g_unified_positron_nl,"Nominal non-linearity function","l");
                
                for(Int_t i = 0; i<4; i++)
                {
                    g_unified_positron_nl_pulls[i]->SetLineColor(i+2);
                    g_unified_positron_nl_pulls[i]->Draw("same");
                    legend3->AddEntry(g_unified_positron_nl_pulls[i],Form("Error function %d",i+1),"l");
                }
                legend3->Draw("same");
                NLCurvesC->Print(("./Images/"+ AnalysisString+ "/Detector/NonlinearityCurves.eps").c_str(),".eps");
                
                delete NLCurvesC;
            #endif
            delete g_unified_positron_nl;
            
            for (Int_t i = 0; i < m_num_unified_nl_pars; i++)
            {
                delete g_unified_positron_nl_pulls[i];
            }
        }
            break;
    }
}

void OscillationReactor :: SetNLParameters(bool mode)
{
    switch (NLModelE)
    {
        case 0://BCW NL Model
            std::cout << "\t \t \t Using BCW NL Model"<< std::endl;
            if(NLMatrix && Mode)
            {   //Randomize BCW nonlinear model parameters
                std::cout << "\t \t \t Randomizing BCW NL parameters" << std::endl;
                
                for (Int_t i = 0; i < 5; i++)
                {
                    rand->SetSeed(0);
                    m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i] + m_bcw_elec_nl_par_error[i] * rand->Gaus(0,1);
                    //                    std::cout <<  m_bcw_elec_nl_par_nominal[i] << " AND " <<  m_bcw_elec_nl_par_error[i] << "AND" << rand->Gaus(0,1) << "IS" <<  m_bcw_elec_nl_par[i]<< std::endl;
                }
                for (Int_t AD = 0; AD < NADs; AD++)
                {
                    rand->SetSeed(0);
                    NLF->SetParameters(m_bcw_elec_nl_par[0],m_bcw_elec_nl_par[1],m_bcw_elec_nl_par[2], m_bcw_elec_nl_par[3], m_bcw_elec_nl_par[4], m_abs_escale * m_rel_escale[AD], m_abs_eoffset+m_rel_eoffset[AD]);
                }
            }
            else
            {
                for (Int_t i = 0; i < 5; i++)
                {
                    m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];//reset
                }
                for (Int_t AD = 0; AD < NADs; AD++)
                {
                    NLF->SetParameters(m_bcw_elec_nl_par[0],m_bcw_elec_nl_par[1],m_bcw_elec_nl_par[2], m_bcw_elec_nl_par[3], m_bcw_elec_nl_par[4], m_abs_escale * m_rel_escale[AD], m_abs_eoffset+m_rel_eoffset[AD]);
                }
            }
            break;
        case 1://LBNL NL Model
            std::cout << "\t \t \t Using LBNL NL Model"<< std::endl;
            if(NLMatrix && Mode)
            {   //Randomize LBNL nonlinear model parameters
                std::cout << "\t \t \t Randomizing LBNL NL parameters" << std::endl;
                for (Int_t i = 0; i < 3; i++)
                {
                    rand->SetSeed(0);
                    m_lbnl_nl_par[i] = m_lbnl_nl_par_nominal[i] + m_lbnl_nl_par_error[i] * rand->Gaus(0,1);//electronics
                    //                std::cout << m_lbnl_nl_par_nominal[i] << " AND " <<  m_lbnl_nl_par_error[i] << "AND" << rand->Gaus(0,1) << "IS" << m_lbnl_nl_par[i] << std::endl;
                }
                for (Int_t AD = 0; AD < NADs; AD++)
                {
                    NLF->SetParameters(m_lbnl_nl_par[0],m_lbnl_nl_par[1],m_lbnl_nl_par[2],m_abs_escale * m_rel_escale[AD],m_abs_eoffset + m_rel_eoffset[AD]);
                }
            }
            else
            {
                for (Int_t i = 0; i < 3; i++)
                {
                    m_lbnl_nl_par[i] = m_lbnl_nl_par_nominal[i];//reset
                }
                for (Int_t AD = 0; AD < NADs; AD++)
                {
                    NLF->SetParameters(m_lbnl_nl_par[0],m_lbnl_nl_par[1],m_lbnl_nl_par[2],m_abs_escale * m_rel_escale[AD],m_abs_eoffset + m_rel_eoffset[AD]);
                }
            }
            break;
        case 2://Unified NL Model
            std::cout << "\t \t \t Using Unified NL Model"<< std::endl;
            if(NLMatrix  && Mode)
            {   //Randomize Unified nonlinear model parameters
                std::cout << "\t \t \t Randomizing Unified NL parameters" << std::endl;
                
                Double_t ranvec[m_num_unified_nl_pars];
                
                for (Int_t i = 0; i < m_num_unified_nl_pars; i++)
                {
                    rand->SetSeed(0);
                    ranvec[i] = rand->Gaus(0,1);
                }
                
                for (Int_t i = 0; i < m_num_unified_nl_pars; i++)
                {
                    m_unified_nl_par[i] = m_unified_nl_par_nominal[i] + m_unified_nl_par_error[i] * ranvec[i];
                }
                
                for (Int_t i = 0; i < m_num_unified_nl_pars; i++)
                {
                    NLF->SetParameter(i,m_unified_nl_par[i]);
                    
                    std::cout << "\t \t \t Random NL Unified Par: " << m_unified_nl_par[i] << std::endl;
                }
                
                for (Int_t AD = 0; AD < NADs; AD++)
                {
                    NLF->SetParameter(m_num_unified_nl_pars, m_abs_escale * m_rel_escale[AD]);
                    NLF->SetParameter(m_num_unified_nl_pars+1, m_abs_eoffset + m_rel_eoffset[AD]);
                    
                    std::cout << "\t \t \t Random Scale: " << m_abs_escale * m_rel_escale[AD] << "; Random Offset: " << m_abs_eoffset + m_rel_eoffset[AD] << std::endl;
                }
            }
            else
            {
                for (Int_t i = 0; i < m_num_unified_nl_pars; i++)
                {
                    NLF->SetParameter(i,m_unified_nl_par[i]);
                    
                    //                    std::cout << "\t \t \t Nominal NL Unified Par: " << m_unified_nl_par[i] << std::endl;
                }
                
                for (Int_t AD = 0; AD < NADs; AD++)
                {
                    NLF->SetParameter(m_num_unified_nl_pars, m_abs_escale * m_rel_escale[AD]);
                    NLF->SetParameter(m_num_unified_nl_pars+1, m_abs_eoffset + m_rel_eoffset[AD]);
                    
                    //                    std::cout << "\t \t \t Nominal Scale: " << m_abs_escale * m_rel_escale[AD] << "; Nominal Offset: " << m_abs_eoffset + m_rel_eoffset[AD] << std::endl;
                }
            }
            break;
        default:
            std::cout << " NO NON-LINEARITY FUNCTION" << std::endl;
            exit(EXIT_FAILURE);
            break;
    }
    
#ifdef PrintEps
    TCanvas* NLC = new TCanvas("NLC","NLC");
    
    NLF->Draw();
    
    NLC->Print(("./Images/"+ AnalysisString+ "/Detector/Nonlinearity.eps").c_str(),".eps");
    delete NLC;
#endif
    
    Interpolation(NLF);//Interpolate the non-linearity function
    

}

void OscillationReactor :: LoadIavCorrection()// From Bryce Littlejohn's results. //I don't use different IAV matrix for each AD since the analysis it's been done for just 1 of them, assume identical ADs.
{
    
    TFile* f;
    TH2D* Correction;
    if(DataSet==2)//P12
    {
        f = new TFile("./IavDistortion/IAVDistortion.root");//Bryce
        Correction = (TH2D*)f->Get("Correction");
    }
    else//P14
    {
        f = new TFile("./IavDistortion/iavMatrix_P14A.root");//Update by Soeren
        Correction = (TH2D*)f->Get("Correction_LS");
    }
    std::cout << "\t \t \t Reading IAV correction file" << std::endl;
    
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
    delete Correction;
    f->Close();
}

// First order true to visible function.
//(x-(Mn-Mp))*(1 - x/Mn*(1.0 - [0]*sqrt(1 - Me*Me/(x-(Mn-Mp))/(x-(Mn-Mp)))))- ((Mn-Mp)*(Mn-Mp) - Me*Me)/(2*Mn)-Me]

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Formula is different from previous version (Mp instead of Mn in 1st order, check which is correct)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t OscillationReactor :: VisibleEnergy1F(Double_t* energ, Double_t* par)
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
}

Double_t OscillationReactor :: VisibleEnergy0F(Double_t* energ, Double_t* par)
{
    return energ[0]-par[0];
}

Double_t OscillationReactor :: ResolutionF(Double_t* energ, Double_t* par)
{
    Double_t e_orig = energ[0];
    Double_t e_sigma = 1.0;
    
    if (e_orig > 0)//To avoid NaNs.
    {
        e_sigma = TMath::Sqrt(par[0]*par[0] + par[1]*par[1]/e_orig + par[2]*par[2]/e_orig/e_orig);
    }
    
    return e_sigma;
}

Double_t OscillationReactor :: NLBCWF(Double_t* energ, Double_t* par)
{
    //Input error
    Double_t e_positron_true = energ[0];
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
Double_t OscillationReactor :: NLLBNLF(Double_t * energ, Double_t * par)
{
    //par[0]: flat energy scale shift for electron
    //par[1]: size of energy scale shift propotional to exp(-1.5*eVis) for electron
    //par[2]: size of energy scale shift due to Ge68 calibration point
    
    Double_t e_positron_true = energ[0];
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

Double_t OscillationReactor :: NLUnifiedF(Double_t * energ, Double_t * par)
{
    Double_t e_positron_true = energ[0];
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
    Double_t random_nl_fac = scinti_nl_fac;//if err = 0 the energy would be the nominal NL function, otherwise energy = nl_nominal_energy(1+Σerrors)
    
    for (Int_t ierr = 0; ierr < m_num_unified_nl_pars; ierr++)
    {
        random_nl_fac += par[ierr]*err[ierr];
    }
    
    Double_t visibleE = random_nl_fac * e_positron_true;
    
    return visibleE * escale_par + escale_offset;
}

void OscillationReactor :: SetSystematic()
{
    if(RelativeEnergyScaleMatrix)
    {
        this->RandomRelativeEnergyScaleMatrix();
    }
    if(IAVMatrixb)
    {
        this->RandomIAVMatrix();
    }
    if(ResolutionMatrix)
    {
        this->RandomResolutionMatrix();
    }
    if(EfficiencyMatrix)
    {
        this->RandomEfficiency();
    }
    if(AllDetectorSystematicsMatrix)
    {
            this->RandomRelativeEnergyScaleMatrix();
    
            this->RandomIAVMatrix();
     
            this->RandomResolutionMatrix();
  
            this->RandomEfficiency();
        
            NLMatrix = 1; //NL case is handled by SetNLParameters();
    }
    //NL case is handled by SetNLParameters();
}

void OscillationReactor :: LoadNominalTotalOscillatedSpectrum()
{
    TFile* InputFile = new TFile(("./RootOutputs/"+ AnalysisString+ "/NominalOutputs/Oscillation.root").c_str());
    Char_t histnameTotalget[100];//47

    InputFile->cd("Total AD Spectra after oscillation");
    for(Int_t AD = 0; AD<NADs; AD++)
    {
        if (Nweeks == 1)
        {
            sprintf(histnameTotalget,"Total spectrum after oscillation at AD%i", AD+1);
        }
        else
        {
            sprintf(histnameTotalget,"Total spectrum after oscillation at AD%i, week%i", AD+1, week+1);
        }

        NominalTotalOscillatedSpectrumAD[AD]= (TH1D*) gDirectory->Get(histnameTotalget);
        
        IntegralNominalTotalOscillatedSpectrumAD[AD] = NominalTotalOscillatedSpectrumAD[AD]->Integral();
        std::cout << "EVENTS IN TOTAL NOMINAL OSCILLATED SPECTRUM : " << NominalTotalOscillatedSpectrumAD[AD]->Integral() << std::endl;
        
    }
    
    delete InputFile;
    
}