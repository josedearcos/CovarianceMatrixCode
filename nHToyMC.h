/*
 The Toy will generate the signals with raw visible energy event by event.
 The output file "roofile_toy.root" includes all information from inital inputs to final outputs step by step.
 
 2014.11.04, the Epositron is flat spectrum, need a random sampling from "true positron spectrum" --- Jixp
 2014.11.04, how to run ---> (1) root (2) .L Toy.cc+ (3) Toy()
 */
#pragma once
#include<iostream>
using namespace std;

#include<cmath>
#include<fstream>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TString.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TEventList.h"
/// select whether nGd or nH IBD analysis - nGd events are currently rejected in the code.  Time cuts are not currently applied.
// 0==special cuts, 1==hydrogen, 64==gadolinium
#define SpecialCuts
//#define Hydrogen
//#define Gadolinium

const double Rshield = 2259.15;  // position of radial shield [mm]
const double ZreflectorTop = 2104.5;  // position of top reflector ESR [mm]
const double ZreflectorBottom = -2027.45;  // position of bottom reflector ESR [mm]

//Uncomment to turn on OAV uncertainty:
const double UncertaintyDeltaRoav = 8.0;  // +-8.0 mm variation in diameter within ADs 1 & 2
const double UncertaintyDeltaZoav = 3.0;  // +-3.0 mm discrepancy between specification and measurement of ADs 1 & 2

//const double deltaRoav = 0;  // +-8.0 mm variation in diameter within ADs 1 & 2
//const double deltaZoav = 0;  // +-3.0 mm discrepancy between specification and measurement of ADs 1 & 2

//#ifdef SpecialCuts  /// Time cuts are not currently applied.
//    const double EpromptCutLow=0.0, EpromptCutHigh=12.0, EdelayCutLow=0.0, EdelayCutHigh=3.3;// 1.5, 12.0, 1.5, 3.3
//    const double TcapCutLow=500.0, TcapCutHigh=400000.0, distCut=500.0;
//    const double radialCut=0.0;//3.9;  // adjust OAV radius [m^2]: nominal=4.0
//#elif def Hydrogen
//    const double EpromptCutLow=1.5, EpromptCutHigh=12.0, EdelayCutLow=1.7955, EdelayCutHigh=2.6535;
//    const double TcapCutLow=500.0, TcapCutHigh=400000.0, distCut=500.0;
//    const double radialCut=0.0;  // adjust OAV radius [m^2]: nominal=4.0
//#elif def Gadolinium
//    const double EpromptCutLow=1.5, EpromptCutHigh=12.0, EdelayCutLow=6.0, EdelayCutHigh=12.0;
//    const double TcapCutLow=500.0, TcapCutHigh=200000.0, distCut=0.0;
//    const double radialCut=0.0;  // adjust OAV radius [m^2]: nominal=4.0
//#endif

/// Time cuts are not currently applied.
#ifdef SpecialCuts  /// Time cuts are not currently applied.
    const double EpromptCutLow=0.0, EpromptCutHigh=12.0, EdelayCutLow=0.0, EdelayCutHigh=12.0; // 1.5, 12.0, 1.5, 3.3
    const double TcapCutLow=0.0, TcapCutHigh=0.0, distCut=500.0;
    const double radialCut=0.0;  // adjust OAV radius^2 [m^2]: nominal=3.964 (MC), 3.935 (data) - unc. 5.1 mm (data).
#elif def Hydrogen
    const double EpromptCutLow=1.5, EpromptCutHigh=12.0, EdelayCutLow=1.7955, EdelayCutHigh=2.6535; // update delayed cut values
    const double TcapCutLow=1000.0, TcapCutHigh=400000.0, distCut=500.0;
    const double radialCut=0.0;  // adjust OAV radius^2 [m^2]: nominal=4.0
#elif def Gadolinium
    const double EpromptCutLow=0.7, EpromptCutHigh=12.0, EdelayCutLow=6.0, EdelayCutHigh=12.0;
    const double TcapCutLow=1000.0, TcapCutHigh=200000.0, distCut=0.0;
    const double radialCut=0.0;  // adjust OAV radius^2 [m^2]: nominal=4.0
#endif

// --- ---- --- ---- --- ---- --- ---- --- ---- --- Global
#define FullEnergyResolution
//#define LoadToyTree

#define LoadTree // To load Eprompt tree, then get the effective entries in a smaller tree (roofile_etree.root)
#define SaveTree // To save toy tree in roofile_toy.root

//#define ReactorShapeinToy //To produce the toy with the reactor shape included
///
TString roostr;

//double hEp_low = 0;
//double hEp_hgh = 12;
//int    hEp_bin = 120;
//
//double hEd_low = 0;
//double hEd_hgh = 3.3; //12;
//int    hEd_bin = 120; //1200;
//
//const int MaxColumnNum = 11; //11 = (R2_binnum+1)

//const Int_t unified_nl_pars = 4;//Number of curves used in the NL error calculation
const Int_t unified_nl_pars = 1;//After the last nl-update, not using marginal curves, instead 1 sigma band
const Int_t InitialSystematic = 0;// Normally 0 to run all the systematics, otherwise select initial and final by changing NSystematic
const Int_t NSystematic = 7;//Nominal,IAV, OAV, NL, Reso, Relative Energy Scale, AllSystematics
//Efficiency(If only global eff this is not applied here actually), if a difference between Data-ToyMC efficiency map is substantial we could apply a gaussian efficiency variation per cell.

class nHToyMC
{
private:
    
    //Logan Plots:
    
//    TH1D *hEp_cc[MaxCellNum];
//    TH1D *hEd_cc[MaxCellNum];
//    TH1D *hEp_cc_clone[MaxCellNum];
//    TH1D *hEd_cc_clone[MaxCellNum];
//    
//    TH1D *hEp_cl[MaxColumnNum];
//    TH1D *hEd_cl[MaxColumnNum];
    
    Double_t GausRelative,GausIAV,GausReso,GausResoCorr,GausOAV;
//    Double_t GausEff;
    Double_t GausNL[unified_nl_pars];
    Double_t ResolutionError,ResolutionErrorUncorrelated;
    
    Double_t deltaRoav;  // +-8.0 mm variation in diameter within ADs 1 & 2
    Double_t deltaZoav;  // +-3.0 mm discrepancy between specification and measurement of ADs 1 & 2

    //Binning parameters:
    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    Double_t evis_bins[MaxNbins+1]; // Single bins between 1.5 and 1.6 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.
    
    Int_t n_evis_bins;
    Int_t n_etrue_bins;
    
    Double_t IAVNominalError; // relative uncertainty of the IAV thickness for each AD.
    Double_t RelativeNominalError; // relative uncertainty of attenuation length for each AD.
    //Systematic parameters;
    bool RelativeEnergyScaleMatrix;
    bool IAVMatrix;
    bool OAVMatrix;
    bool NLMatrix;
    bool ResolutionMatrix;
//    bool EfficiencyMatrix;
    bool AllMatrix;
    string SystematicS;

    //AD configuration parameters:
    Int_t NADs;

    /// find cell
    TH2D *hist_findbin;
    Int_t global_bin_num;
    Int_t local_xbin    ;
    Int_t local_ybin    ;
    Int_t local_zbin    ;
    Double_t usr_r2_P   ;
    Double_t usr_z_P    ;
    Double_t usr_r2_Ng  ;
    Double_t usr_z_Ng   ;
    
    /// non-linearity: NL
    TGraph *graph_electron_LY;
    TGraph *graph_gamma_LY;
    TGraph *graph_electronic;
    
    TH1D* ErrorBand;
    
    TGraph *g_unified_positron_nl;
    TGraph *g_unified_positron_nl_pulls[unified_nl_pars];
    
    /// non-uniformity: NU
    TH2D *hist_map_attenuation;
    TH2D *hist_map_pmt_coverage;
    
    ///
#ifdef ReactorShapeinToy
    TH1D *h_Ev_normal;
#endif
    TRandom3 *gRandom3;
    TRandom3 *RandomSysUncorr;
    TRandom3 *RandomSysCorr;
    
    ///
    TF1 *roofunc_EnergyResolution;
    
    Double_t func_EnergyResolution(double*,double*);
    void func_initialization();
    Int_t  RootCellToVisCell(Int_t RootCell);
    TH1D* TruePredictionH[NSystematic][MaxDetectors][VolumeX][VolumeY];
    TH1D* VisiblePredictionH[NSystematic][MaxDetectors][VolumeX][VolumeY];
    TH1D* DelayedVisiblePredictionH[NSystematic][MaxDetectors][VolumeX][VolumeY];
    TH2D* TransMatrixH[NSystematic][MaxDetectors][VolumeX][VolumeY];
    TH2D* HighResoTransMatrixH[NSystematic][MaxDetectors][VolumeX][VolumeY];
    TH2D* MatrixH[NSystematic][MaxDetectors][VolumeX][VolumeY];
public:
    nHToyMC(NominalData*);//constructor
    ~nHToyMC();//destructor
    void Toy(bool);//main
//    TH2D* LoadnHMatrix(Int_t,Int_t,Int_t);
//    Double_t GetEventsByCell(Int_t,Int_t,Int_t);
    TH2D* HighResoMatrixH[NSystematic][MaxDetectors][VolumeX][VolumeY];
    Double_t PercentualEvents[NSystematic][MaxDetectors][VolumeX][VolumeY];
};

// constructor:
nHToyMC :: nHToyMC(NominalData* Data)
{
    /// find cell
    global_bin_num = 0;
    local_xbin     = 0;
    local_ybin     = 0;
    local_zbin     = 0;
    usr_r2_P    = 0;
    usr_z_P     = 0;
    usr_r2_Ng   = 0;
    usr_z_Ng    = 0;
    
    InitialEnergy = Data->GetEmin();
    FinalEnergy = Data->GetEmax();
    InitialVisibleEnergy = Data->GetEVisMin();
    FinalVisibleEnergy = Data->GetEVisMax();
    
    n_evis_bins = Data->GetVisibleBins();
    
    IAVNominalError = Data->GetIAVError();
    RelativeNominalError = Data->GetRelativeEnergyError();
    
    for (Int_t i = 0; i <= n_evis_bins; i++)
    {
        evis_bins[i] = Data->GetVisibleBinningArray(i);
    }
    
    n_etrue_bins = Data->GetTrueBins();
    
    for (Int_t i = 0; i <= n_etrue_bins; i++)
    {
        enu_bins[i] = Data->GetTrueBinningArray(i);
    }
    
    //Systematic matrix selection:
//    RelativeEnergyScaleMatrix = Data->GetRelativeEnergyScaleMatrix();
//    //    RelativeEnergyOffsetMatrix = Data->GetRelativeEnergyOffsetMatrix();
//    //    AbsoluteEnergyScaleMatrix = Data->GetAbsoluteEnergyScaleMatrix();
//    //    AbsoluteEnergyOffsetMatrix = Data->GetAbsoluteEnergyOffsetMatrix();
//    IAVMatrix = Data->GetIAVMatrix();
//    OAVMatrix = Data->GetOAVMatrix();
//    NLMatrix = Data->GetNLMatrix();
//    ResolutionMatrix = Data->GetResolutionMatrix();
//    EfficiencyMatrix = Data->GetEfficiencyMatrix();
    
    //Errors in parameters:
    ResolutionError = Data->GetResolutionError();
    ResolutionErrorUncorrelated = Data->GetResoUncorrelatedError();
    
    NADs = Data->GetADs();
    
    for(Int_t AD = 0; AD < NADs; AD++)
    {
        for(Int_t Systematic = InitialSystematic; Systematic<NSystematic;Systematic++)
        {
            for(int idx=0; idx<VolumeX; idx++)
            {
                for(int idy=0; idy<VolumeY; idy++)
                {
                    HighResoMatrixH[Systematic][AD][idx][idy] = new TH2D(Form("Fine_nHResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),Form("Fine_nHResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),MatrixBins,0,FinalVisibleEnergy,MatrixBins,0,FinalVisibleEnergy);//from true to visible
                    
                    MatrixH[Systematic][AD][idx][idy] = new TH2D(Form("nHResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),Form("nHResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),n_etrue_bins,enu_bins,n_evis_bins,evis_bins);//from true to visible
                    
                    TransMatrixH[Systematic][AD][idx][idy] = new TH2D(Form("nHTransResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),Form("nHTransResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),n_evis_bins,evis_bins,n_etrue_bins,enu_bins);//from visible to true
                    
                    HighResoTransMatrixH[Systematic][AD][idx][idy] = new TH2D(Form("Fine_nHTransResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),Form("Fine_nHTransResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),MatrixBins,0,FinalVisibleEnergy,MatrixBins,0,FinalVisibleEnergy);//from visible to true
                    
                    TruePredictionH[Systematic][AD][idx][idy] = new TH1D(Form("Pred_AD%i, Cell%i,%i",AD+1,idx,idy),Form("Pred_AD%i, Cell%i,%i",AD+1,idx,idy), Nbins, InitialEnergy, FinalEnergy);
                    
                    VisiblePredictionH[Systematic][AD][idx][idy] = new TH1D(Form("VisPred_AD%i, Cell%i,%i",AD+1,idx,idy),Form("VisPred_AD%i, Cell%i,%i",AD+1,idx,idy), MatrixBins, InitialVisibleEnergy, FinalVisibleEnergy);
                    
                    DelayedVisiblePredictionH[Systematic][AD][idx][idy] = new TH1D(Form("DelayedVisPred_AD%i, Cell%i,%i",AD+1,idx,idy),Form("DelayedVisPred_AD%i, Cell%i,%i",AD+1,idx,idy), MatrixBins, InitialVisibleEnergy, FinalVisibleEnergy);
                }
            }
        }
    }
}

nHToyMC :: ~nHToyMC()
{
    for(Int_t Systematic = InitialSystematic; Systematic<NSystematic;Systematic++)
    {
        for(Int_t AD = 0; AD < NADs; AD++)
        {
            for(int idx=0; idx<VolumeX; idx++)
            {
                for(int idy=0; idy<VolumeY; idy++)
                {
                    delete TruePredictionH[Systematic][AD][idx][idy];
                    delete HighResoMatrixH[Systematic][AD][idx][idy];
                    delete MatrixH[Systematic][AD][idx][idy];
                    delete TransMatrixH[Systematic][AD][idx][idy];
                    delete VisiblePredictionH[Systematic][AD][idx][idy];
                    delete DelayedVisiblePredictionH[Systematic][AD][idx][idy];
                }
            }
        }
    }
}
///
Double_t nHToyMC :: func_EnergyResolution(Double_t *x, Double_t *par)// from Logan
{
    
    Double_t ResolutionBias;
    
    Double_t corr_bias = ResolutionError * GausResoCorr;
    
    ResolutionBias = corr_bias + ResolutionErrorUncorrelated * GausReso;
    
#ifdef FullEnergyResolution
    Double_t E = x[0];
    Double_t R = par[0];
   // Double_t res = sqrt(pow(0.004*E,2) + E*(pow(0.082,2)+R*pow(0.031,2)) + pow((0.028),2))+ResolutionBias; // 7.5/sqrt(E) + 0.9 + 0.865*(R-0.98);
    
    Double_t res = sqrt( 0.004*0.004*E*E + E*(0.082*0.082+R*0.031*0.031) + 0.028*0.028 )+ResolutionBias;  // 7.5/sqrt(E) + 0.9 + 0.865*(R-0.98);

    //res = E*res;
#else
    Double_t P0 = 0.1127;
    Double_t P1 = 0.0192;
    Double_t  R = x[0];
    Double_t res = P0 + P1*R;
#endif
    
    return res;
}

// --- ---- --- ---- --- ---- --- ---- --- ---- --- Function declaration

Int_t nHToyMC :: RootCellToVisCell(Int_t RootCell)
{
    // 01 02 03 04 05 06 07 08 09 10
    // 11 12 13 14 15 16 17 18 19 20
    // 21 22 23 24 25 26 27 28 29 30
    // ...
    // 91 92 93 94 95 96 97 98 99 100
    
    Int_t cell = 0;
    
    TH2D *hist_findbin_temp = new TH2D("hist_findbin_temp","", R2_binnum,R2_lower,R2_upper, Z_binnum,Z_lower,Z_upper);
    hist_findbin_temp->GetBinXYZ(RootCell, local_xbin, local_ybin, local_zbin);
    
    if( local_zbin==0 && local_xbin>=1 && local_xbin<=R2_binnum && local_ybin>=1 && local_ybin<=Z_binnum )
    {
        cell = R2_binnum*Z_binnum - local_ybin*R2_binnum + local_xbin;
    }
    else
    {
        cell = (R2_binnum+2)*(Z_binnum+2);
    }
    
    delete hist_findbin_temp;
    return cell;
}

// --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- ---
// --- ---- --- ---- --- ---- --- ---- ---      --- ---- --- ---- --- ---- --- ---- --- ---- ---
// --- ---- --- ---- --- ---- --- ---- --- MAIN --- ---- --- ---- --- ---- --- ---- --- ---- ---
// --- ---- --- ---- --- ---- --- ---- ---      --- ---- --- ---- --- ---- --- ---- --- ---- ---
// --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- ---

void nHToyMC :: Toy(bool NormalMode)//Normal mode selects between producing systematic variations with a determined number of samples at a time or the whole random sample tree. If we use NormalMode = 0, use FLATSPECTRUM, since it is assumed that there are about same number of energy events per bin
{
        func_initialization();
        
//        Double_t cap_Ev;             // neutrino energy
//        Double_t cap_Target;         // neutron captrure target. e.g. 1 for H
//        Double_t X_Ng, Y_Ng, Z_Ng;   // vertex of gamma(s) induced by neutron capture
//        Double_t Edep_Ng;            // deposited energy of gamma(s) induced by neutron capture
//        Double_t EtotScint_Ng;       // the initial energy of gamma(s) going into scintillator(LS or GdLS)
//        Double_t Etot_Ng;            // the total energy of gamma(s) when generated
//        Double_t EdepE_Ng;           // the deposited energy of secondary electrons induced by the gamma(s)
//        Double_t EdepEA_Ng;
//        Double_t EtotE_Ng;           // the total energy of electrons induced by the gamma(s)
//        Double_t Etot_N;             // the kinetic energy of neutron generated from the IBD
//        
//        Double_t X_P, Y_P, Z_P;  // vetex of positron
//        Double_t Etot_P;         // the initial energy of positron when generated
//        Double_t Etot_Pg1;       // the initial energy of one gamma generated from IBD
//        Double_t Etot_Pg2;       // the initial energy of another gamma generated from IBD
//        Double_t EtotScint_P;    // the initial energy of positron going into scintillator(LS or GdLS)
//        Double_t EtotScint_Pg;   // the initial energy of gamma(s) induced by the positron before annihilation
//        Double_t EtotScint_Pg1;  // the initial energy of one gamma going into scintillator(LS or GdLS)
//        Double_t EtotScint_Pg2;  // ... the meaning of variables bellow can be refered to the neutron capture information
//        Double_t Edep_P;
//        Double_t EdepAcrylic_P;
//        Double_t Edep_Pg1;
//        Double_t Edep_Pg2;
//        Double_t EdepE_P;
//        Double_t EdepEA_P;
//        Double_t EtotE_P;
//        Double_t EdepE_Pg1;
//        Double_t EdepEA_Pg1;
//        Double_t EtotE_Pg1;
//        Double_t EdepE_Pg2;
//        Double_t EdepEA_Pg2;
//        Double_t EtotE_Pg2;
        
        // define variables to get from input file
        double Ev;
        //double cap_Target;
        //int IBDvolume;
        double Tdiff;
        double X_Ng, Y_Ng, Z_Ng;
        //double Etot_Ng;
        double EtotScint_Ng[10];
        //double Edep_Ng;
        double EdepScint_Ng[10];
        //double EdepE_Ng;
        //double EdepEA_Ng;
        //double EtotE_Ng;
        double Etot_N;
        double EtotScint_N;
        
        double X_P, Y_P, Z_P;
        double Etot_P;
        //double Etot_Pg1;
        //double Etot_Pg2;
        double EtotScint_P;
        //double EtotScint_Pg;
        double EtotScint_Pg1;
        double EtotScint_Pg2;
        //double Edep_P;
        double EdepScint_P;
        double EdepAcrylic_P;  // generate a map for EdepAcrylic_P (two energy ranges), rejecting events > R2=3,500,000 and associating sum of GdLS region to Bryce's uncertainty.
        //double Edep_Pg1;
        double EdepScint_Pg1;
        //double Edep_Pg2;
        double EdepScint_Pg2;
        //double EdepE_P;
        //double EdepEA_P;
        //double EtotE_P;
        //double EdepE_Pg1;
        //double EdepEA_Pg1;
        //double EtotE_Pg1;
        //double EdepE_Pg2;
        //double EdepEA_Pg2;
        //double EtotE_Pg2;
        
        Int_t seed_generator = 1;
        Int_t seed_generator_corr = 1;
        Int_t seed_generator_uncorr = 1;
#ifdef LoadTree
        
        ////////////
        
        TFile *roofile_input = new TFile("./Inputs/HInputs/Data/Eprompt.root", "read");
        TTree *wtree = (TTree*)roofile_input->Get("wtree");
        long entries_wtree = wtree->GetEntries();
        cout<<" ---> input entries: "<<entries_wtree<<endl;
        
//        wtree->SetBranchAddress("cap_Target",   &cap_Target);
//        wtree->SetBranchAddress("X_Ng",         &X_Ng);
//        wtree->SetBranchAddress("Y_Ng",         &Y_Ng);
//        wtree->SetBranchAddress("Z_Ng",         &Z_Ng);
//        wtree->SetBranchAddress("Edep_Ng",      &Edep_Ng);
//        wtree->SetBranchAddress("EtotScint_Ng", &EtotScint_Ng);
//        wtree->SetBranchAddress("Etot_Ng",      &Etot_Ng);
//        wtree->SetBranchAddress("EdepE_Ng",     &EdepE_Ng);
//        wtree->SetBranchAddress("EdepEA_Ng",    &EdepEA_Ng);
//        wtree->SetBranchAddress("EtotE_Ng",     &EtotE_Ng);
//        wtree->SetBranchAddress("Etot_N",       &Etot_N);
//        wtree->SetBranchAddress("X_P",          &X_P);
//        wtree->SetBranchAddress("Y_P",          &Y_P);
//        wtree->SetBranchAddress("Z_P",          &Z_P);
//        wtree->SetBranchAddress("Etot_P",       &Etot_P);
//        wtree->SetBranchAddress("Etot_Pg1",     &Etot_Pg1);
//        wtree->SetBranchAddress("Etot_Pg2",     &Etot_Pg2);
//        wtree->SetBranchAddress("EtotScint_P",  &EtotScint_P);
//        wtree->SetBranchAddress("EtotScint_Pg", &EtotScint_Pg);
//        wtree->SetBranchAddress("EtotScint_Pg1", &EtotScint_Pg1);
//        wtree->SetBranchAddress("EtotScint_Pg2", &EtotScint_Pg2);
//        wtree->SetBranchAddress("Edep_P",        &Edep_P);
//        wtree->SetBranchAddress("EdepAcrylic_P", &EdepAcrylic_P);
//        wtree->SetBranchAddress("Edep_Pg1",      &Edep_Pg1);
//        wtree->SetBranchAddress("Edep_Pg2",      &Edep_Pg2);
//        wtree->SetBranchAddress("EdepE_P",       &EdepE_P);
//        wtree->SetBranchAddress("EdepEA_P",      &EdepEA_P);
//        wtree->SetBranchAddress("EtotE_P",       &EtotE_P);
//        wtree->SetBranchAddress("EdepE_Pg1",     &EdepE_Pg1);
//        wtree->SetBranchAddress("EdepEA_Pg1",    &EdepEA_Pg1);
//        wtree->SetBranchAddress("EtotE_Pg1",     &EtotE_Pg1);
//        wtree->SetBranchAddress("EdepE_Pg2",     &EdepE_Pg2);
//        wtree->SetBranchAddress("EdepEA_Pg2",    &EdepEA_Pg2);
//        wtree->SetBranchAddress("EtotE_Pg2",     &EtotE_Pg2);
        
        //wtree->SetBranchAddress("cap_Target",   &cap_Target);
        //wtree->SetBranchAddress("IBDvolume",    &IBDvolume);
        wtree->SetBranchAddress("Tdiff",        &Tdiff);
        wtree->SetBranchAddress("X_Ng",         &X_Ng);
        wtree->SetBranchAddress("Y_Ng",         &Y_Ng);
        wtree->SetBranchAddress("Z_Ng",         &Z_Ng);
        //wtree->SetBranchAddress("Edep_Ng",      &Edep_Ng);
        wtree->SetBranchAddress("EdepScint_Ng", &EdepScint_Ng);
        wtree->SetBranchAddress("EtotScint_Ng", &EtotScint_Ng);
        //wtree->SetBranchAddress("Etot_Ng",      &Etot_Ng);
        //wtree->SetBranchAddress("EdepE_Ng",     &EdepE_Ng);
        //wtree->SetBranchAddress("EdepEA_Ng",    &EdepEA_Ng);
        //wtree->SetBranchAddress("EtotE_Ng",     &EtotE_Ng);
        wtree->SetBranchAddress("Etot_N",       &Etot_N);
        wtree->SetBranchAddress("EtotScint_N",  &EtotScint_N);
        wtree->SetBranchAddress("X_P",          &X_P);
        wtree->SetBranchAddress("Y_P",          &Y_P);
        wtree->SetBranchAddress("Z_P",          &Z_P);
        wtree->SetBranchAddress("Etot_P",       &Etot_P);
        //wtree->SetBranchAddress("Etot_Pg1",     &Etot_Pg1);
        //wtree->SetBranchAddress("Etot_Pg2",     &Etot_Pg2);
        wtree->SetBranchAddress("EtotScint_P",  &EtotScint_P);
        //wtree->SetBranchAddress("EtotScint_Pg", &EtotScint_Pg);
        wtree->SetBranchAddress("EtotScint_Pg1", &EtotScint_Pg1);
        wtree->SetBranchAddress("EtotScint_Pg2", &EtotScint_Pg2);
        //wtree->SetBranchAddress("Edep_P",        &Edep_P);
        wtree->SetBranchAddress("EdepScint_P",   &EdepScint_P);
        wtree->SetBranchAddress("EdepAcrylic_P", &EdepAcrylic_P);
        //wtree->SetBranchAddress("Edep_Pg1",      &Edep_Pg1);
        wtree->SetBranchAddress("EdepScint_Pg1", &EdepScint_Pg1);
        //wtree->SetBranchAddress("Edep_Pg2",      &Edep_Pg2);
        wtree->SetBranchAddress("EdepScint_Pg2", &EdepScint_Pg2);
        //wtree->SetBranchAddress("EdepE_P",       &EdepE_P);
        //wtree->SetBranchAddress("EdepEA_P",      &EdepEA_P);
        //wtree->SetBranchAddress("EtotE_P",       &EtotE_P);
        //wtree->SetBranchAddress("EdepE_Pg1",     &EdepE_Pg1);
        //wtree->SetBranchAddress("EdepEA_Pg1",    &EdepEA_Pg1);
        //wtree->SetBranchAddress("EtotE_Pg1",     &EtotE_Pg1);
        //wtree->SetBranchAddress("EdepE_Pg2",     &EdepE_Pg2);
        //wtree->SetBranchAddress("EdepEA_Pg2",    &EdepEA_Pg2);
        //wtree->SetBranchAddress("EtotE_Pg2",     &EtotE_Pg2);

        
        ////////////
        
//        TH1D *h_Ev_toyMC_input = new TH1D("h_Ev_toyMC_input","h_Ev_toyMC_input",480,0,12);
        
        float EupperLim = 20.0;  // MeV
        TH1D *h_Ev_toyMC_input = new TH1D("h_Ev_toyMC_input","h_Ev_toyMC_input", (int)EupperLim*40, 0.0, EupperLim);
        
        TFile *roofile_etree = new TFile("./Inputs/HInputs/Data/roofile_etree.root", "recreate");
        TTree *etree = new TTree("etree", "effective entries of wtree");
        
//        etree->Branch("cap_Ev",       &cap_Ev,       "cap_Ev/D");
//        etree->Branch("cap_Target",   &cap_Target,   "cap_Target/D");
//        etree->Branch("X_Ng",         &X_Ng,         "X_Ng/D");
//        etree->Branch("Y_Ng",         &Y_Ng,         "Y_Ng/D");
//        etree->Branch("Z_Ng",         &Z_Ng,         "Z_Ng/D");
//        etree->Branch("Edep_Ng",      &Edep_Ng,      "Edep_Ng/D");
//        etree->Branch("EtotScint_Ng", &EtotScint_Ng, "EtotScint_Ng/D");
//        etree->Branch("Etot_Ng",      &Etot_Ng,      "Etot_Ng/D");
//        etree->Branch("EdepE_Ng",     &EdepE_Ng,     "EdepE_Ng/D");
//        etree->Branch("EdepEA_Ng",    &EdepEA_Ng,    "EdepEA_Ng/D");
//        etree->Branch("EtotE_Ng",     &EtotE_Ng,     "EtotE_Ng/D");
//        etree->Branch("Etot_N",       &Etot_N,       "Etot_N/D");
//        etree->Branch("X_P",          &X_P,          "X_P/D");
//        etree->Branch("Y_P",          &Y_P,          "Y_P/D");
//        etree->Branch("Z_P",          &Z_P,          "Z_P/D");
//        etree->Branch("Etot_P",       &Etot_P,       "Etot_P/D");
//        etree->Branch("Etot_Pg1",     &Etot_Pg1,     "Etot_Pg1/D");
//        etree->Branch("Etot_Pg2",     &Etot_Pg2,     "Etot_Pg2/D");
//        etree->Branch("EtotScint_P",  &EtotScint_P,  "EtotScint_P/D");
//        etree->Branch("EtotScint_Pg", &EtotScint_Pg, "EtotScint_Pg/D");
//        etree->Branch("EtotScint_Pg1", &EtotScint_Pg1, "EtotScint_Pg1/D");
//        etree->Branch("EtotScint_Pg2", &EtotScint_Pg2, "EtotScint_Pg2/D");
//        etree->Branch("Edep_P",        &Edep_P,        "Edep_P/D");
//        etree->Branch("EdepAcrylic_P", &EdepAcrylic_P, "EdepAcrylic_P/D");
//        etree->Branch("Edep_Pg1",      &Edep_Pg1,      "Edep_Pg1/D");
//        etree->Branch("Edep_Pg2",      &Edep_Pg2,      "Edep_Pg2/D");
//        etree->Branch("EdepE_P",       &EdepE_P,       "EdepE_P/D");
//        etree->Branch("EdepEA_P",      &EdepEA_P,      "EdepEA_P/D");
//        etree->Branch("EtotE_P",       &EtotE_P,       "EtotE_P/D");
//        etree->Branch("EdepE_Pg1",     &EdepE_Pg1,     "EdepE_Pg1/D");
//        etree->Branch("EdepEA_Pg1",    &EdepEA_Pg1,    "EdepEA_Pg1/D");
//        etree->Branch("EtotE_Pg1",     &EtotE_Pg1,     "EtotE_Pg1/D");
//        etree->Branch("EdepE_Pg2",     &EdepE_Pg2,     "EdepE_Pg2/D");
//        etree->Branch("EdepEA_Pg2",    &EdepEA_Pg2,    "EdepEA_Pg2/D");
//        etree->Branch("EtotE_Pg2",     &EtotE_Pg2,     "EtotE_Pg2/D");
        
        etree->Branch("Ev",           &Ev,           "Ev/D");
        //etree->Branch("cap_Target",   &cap_Target,   "cap_Target/D");
        //etree->Branch("IBDvolume",    &IBDvolume,    "IBDvolume/I");
        etree->Branch("Tdiff",        &Tdiff,        "Tdiff/D");
        etree->Branch("X_Ng",         &X_Ng,         "X_Ng/D");
        etree->Branch("Y_Ng",         &Y_Ng,         "Y_Ng/D");
        etree->Branch("Z_Ng",         &Z_Ng,         "Z_Ng/D");
        //etree->Branch("Edep_Ng",      &Edep_Ng,      "Edep_Ng/D");
        etree->Branch("EdepScint_Ng", EdepScint_Ng, "EdepScint_Ng[10]/D");
        etree->Branch("EtotScint_Ng", EtotScint_Ng, "EtotScint_Ng[10]/D");
        //etree->Branch("Etot_Ng",      &Etot_Ng,      "Etot_Ng/D");
        //etree->Branch("EdepE_Ng",     &EdepE_Ng,     "EdepE_Ng/D");
        //etree->Branch("EdepEA_Ng",    &EdepEA_Ng,    "EdepEA_Ng/D");
        //etree->Branch("EtotE_Ng",     &EtotE_Ng,     "EtotE_Ng/D");
        etree->Branch("Etot_N",       &Etot_N,       "Etot_N/D");
        etree->Branch("EtotScint_N",  &EtotScint_N,  "EtotScint_N/D");
        etree->Branch("X_P",          &X_P,          "X_P/D");
        etree->Branch("Y_P",          &Y_P,          "Y_P/D");
        etree->Branch("Z_P",          &Z_P,          "Z_P/D");
        etree->Branch("Etot_P",       &Etot_P,       "Etot_P/D");
        //etree->Branch("Etot_Pg1",     &Etot_Pg1,     "Etot_Pg1/D");
        //etree->Branch("Etot_Pg2",     &Etot_Pg2,     "Etot_Pg2/D");
        etree->Branch("EtotScint_P",  &EtotScint_P,  "EtotScint_P/D");
        //etree->Branch("EtotScint_Pg", &EtotScint_Pg, "EtotScint_Pg/D");
        etree->Branch("EtotScint_Pg1", &EtotScint_Pg1, "EtotScint_Pg1/D");
        etree->Branch("EtotScint_Pg2", &EtotScint_Pg2, "EtotScint_Pg2/D");
        //etree->Branch("Edep_P",        &Edep_P,        "Edep_P/D");
        etree->Branch("EdepScint_P",   &EdepScint_P,   "EdepScint_P/D");
        etree->Branch("EdepAcrylic_P", &EdepAcrylic_P, "EdepAcrylic_P/D");
        //etree->Branch("Edep_Pg1",      &Edep_Pg1,      "Edep_Pg1/D");
        etree->Branch("EdepScint_Pg1", &EdepScint_Pg1, "EdepScint_Pg1/D");
        //etree->Branch("Edep_Pg2",      &Edep_Pg2,      "Edep_Pg2/D");
        etree->Branch("EdepScint_Pg2", &EdepScint_Pg2, "EdepScint_Pg2/D");
        //etree->Branch("EdepE_P",       &EdepE_P,       "EdepE_P/D");
        //etree->Branch("EdepEA_P",      &EdepEA_P,      "EdepEA_P/D");
        //etree->Branch("EtotE_P",       &EtotE_P,       "EtotE_P/D");
        //etree->Branch("EdepE_Pg1",     &EdepE_Pg1,     "EdepE_Pg1/D");
        //etree->Branch("EdepEA_Pg1",    &EdepEA_Pg1,    "EdepEA_Pg1/D");
        //etree->Branch("EtotE_Pg1",     &EtotE_Pg1,     "EtotE_Pg1/D");
        //etree->Branch("EdepE_Pg2",     &EdepE_Pg2,     "EdepE_Pg2/D");
        //etree->Branch("EdepEA_Pg2",    &EdepEA_Pg2,    "EdepEA_Pg2/D");
        //etree->Branch("EtotE_Pg2",     &EtotE_Pg2,     "EtotE_Pg2/D");
    
        for(long ientry=0; ientry<entries_wtree; ientry++)
        {
            wtree->GetEntry(ientry);
            
            cout.precision(4);
            if(ientry%1000000==0)
                cout<<" ---> processing MC spectrum "<<ientry*100./entries_wtree<<"%"<<endl;
            
            //////
//            
//            if(cap_Target<0) continue;
//            if(Edep_Ng<0) continue;
//            if(EtotScint_Ng<0) continue;
//            if(Etot_Ng<0) continue;
//            if(EdepE_Ng<0) continue;
//            if(EdepEA_Ng<0) continue;
//            if(EtotE_Ng<0) continue;
//            if(Etot_N<0) continue;
//            if(Etot_P<0) continue;
//            if(Etot_Pg1<0) continue;
//            if(Etot_Pg2<0) continue;
//            if(EtotScint_P<0) continue;
//            if(EtotScint_Pg<0) continue;
//            if(EtotScint_Pg1<0) continue;
//            if(EtotScint_Pg2<0) continue;
//            if(Edep_P<0) continue;
//            if(EdepAcrylic_P<0) continue;
//            if(Edep_Pg1<0) continue;
//            if(Edep_Pg2<0) continue;
//            if(EdepE_P<0) continue;
//            if(EdepEA_P<0) continue;
//            if(EtotE_P<0) continue;
//            if(EdepE_Pg1<0) continue;
//            if(EdepEA_Pg1<0) continue;
//            if(EtotE_Pg1<0) continue;
//            if(EdepE_Pg2<0) continue;
//            if(EdepEA_Pg2<0) continue;
//            if(EtotE_Pg2<0) continue;
            
            //////
            
            //if(cap_Target<0) continue;
            //if(IBDvolume<0) continue;
            if(EdepScint_Ng[0]<0) continue;
            if(EtotScint_Ng[0]<0) continue;
            //if(Etot_Ng<0) continue;
            //if(EdepE_Ng<0) continue;
            //if(EdepEA_Ng<0) continue;
            //if(EtotE_Ng<0) continue;
            if(Etot_N<0) continue;
            if(EtotScint_N<0) continue;
            if(Etot_P<0) continue;
            //if(Etot_Pg1<0) continue;
            //if(Etot_Pg2<0) continue;
            if(EtotScint_P<0) continue;
            //if(EtotScint_Pg<0) continue;
            if(EtotScint_Pg1<0) continue;
            if(EtotScint_Pg2<0) continue;
            if(EdepScint_P<0) continue;
            if(EdepAcrylic_P<0) continue;
            if(EdepScint_Pg1<0) continue;
            if(EdepScint_Pg2<0) continue;
            //if(EdepE_P<0) continue;
            //if(EdepEA_P<0) continue;
            //if(EtotE_P<0) continue;
            //if(EdepE_Pg1<0) continue;
            //if(EdepEA_Pg1<0) continue;
            //if(EtotE_Pg1<0) continue;
            //if(EdepE_Pg2<0) continue;
            //if(EdepEA_Pg2<0) continue;
            //if(EtotE_Pg2<0) continue;

            //            cap_Ev = Etot_P +Etot_N + DeltaM;
            //
            //            h_Ev_toyMC_input->Fill(cap_Ev);

            
            Ev = Etot_P + Etot_N + IBDthreshold;
            
            h_Ev_toyMC_input->Fill(Ev);
            
#ifndef ReactorShapeinToy
            etree->Fill();//To fill a flat spectrum to produce a response matrix with nearly equal statistics in every column and about 3 times more events in total
#endif//ReactorShapeInToy
        }
        
#ifdef ReactorShapeinToy        ///////// Uncomment #ReactorShapeinToy to include reactor shape in the toy:
        
//        Double_t max_h_Ev_toyMC_input = h_Ev_toyMC_input->GetBinContent( h_Ev_toyMC_input->GetMaximumBin() );
//        Double_t max_h_Ev_normal = h_Ev_normal->GetBinContent( h_Ev_normal->GetMaximumBin() );
//        h_Ev_normal->Scale(max_h_Ev_toyMC_input/max_h_Ev_normal);
        
        double max_h_Ev_normal = h_Ev_normal->GetBinContent( h_Ev_normal->GetMaximumBin() );
        h_Ev_normal->Scale(1.0/max_h_Ev_normal);
        
        int ibin=0, seed_rand=0;
        double prob=0, rand=0;
        
        for(long ientry=0; ientry<entries_wtree; ientry++)
        {
            wtree->GetEntry(ientry);
            
            cout.precision(4);
            if(ientry%200000==0)
                cout<<" ---> processing reactor shape spectrum "<<ientry*100./entries_wtree<<"%"<<endl;
            
            //////
            
//            if(cap_Target<0) continue;
//            if(Edep_Ng<0) continue;
//            if(EtotScint_Ng<0) continue;
//            if(Etot_Ng<0) continue;
//            if(EdepE_Ng<0) continue;
//            if(EdepEA_Ng<0) continue;
//            if(EtotE_Ng<0) continue;
//            if(Etot_N<0) continue;
//            if(Etot_P<0) continue;
//            if(Etot_Pg1<0) continue;
//            if(Etot_Pg2<0) continue;
//            if(EtotScint_P<0) continue;
//            if(EtotScint_Pg<0) continue;
//            if(EtotScint_Pg1<0) continue;
//            if(EtotScint_Pg2<0) continue;
//            if(Edep_P<0) continue;
//            if(EdepAcrylic_P<0) continue;
//            if(Edep_Pg1<0) continue;
//            if(Edep_Pg2<0) continue;
//            if(EdepE_P<0) continue;
//            if(EdepEA_P<0) continue;
//            if(EtotE_P<0) continue;
//            if(EdepE_Pg1<0) continue;
//            if(EdepEA_Pg1<0) continue;
//            if(EtotE_Pg1<0) continue;
//            if(EdepE_Pg2<0) continue;
//            if(EdepEA_Pg2<0) continue;
//            if(EtotE_Pg2<0) continue;
            
            //if(cap_Target<0) continue;
            //if(IBDvolume<0) continue;
            if(EdepScint_Ng[0]<0) continue;
            if(EtotScint_Ng[0]<0) continue;
            //if(Etot_Ng<0) continue;
            //if(EdepE_Ng<0) continue;
            //if(EdepEA_Ng<0) continue;
            //if(EtotE_Ng<0) continue;
            if(Etot_N<0) continue;
            if(EtotScint_N<0) continue;
            if(Etot_P<0) continue;
            //if(Etot_Pg1<0) continue;
            //if(Etot_Pg2<0) continue;
            if(EtotScint_P<0) continue;
            //if(EtotScint_Pg<0) continue;
            if(EtotScint_Pg1<0) continue;
            if(EtotScint_Pg2<0) continue;
            if(EdepScint_P<0) continue;
            if(EdepAcrylic_P<0) continue;
            if(EdepScint_Pg1<0) continue;
            if(EdepScint_Pg2<0) continue;
            //if(EdepE_P<0) continue;
            //if(EdepEA_P<0) continue;
            //if(EtotE_P<0) continue;
            //if(EdepE_Pg1<0) continue;
            //if(EdepEA_Pg1<0) continue;
            //if(EtotE_Pg1<0) continue;
            //if(EdepE_Pg2<0) continue;
            //if(EdepEA_Pg2<0) continue;
            //if(EtotE_Pg2<0) continue;
            
            //////
            
//            cap_Ev = Etot_P +Etot_N +DeltaM;
//            if( cap_Ev>12 ) continue;
//            
//            Int_t ibin = h_Ev_toyMC_input->FindBin( cap_Ev );
//            Double_t prob = h_Ev_normal->GetBinContent(ibin);
            
            Ev = Etot_P + Etot_N + IBDthreshold;
            
            if( Ev>EupperLim ) {
                cout<<"Determined neutrino energy > "<< EupperLim <<" MeV: "<< Ev <<" at "<< ientry*100./entries_wtree <<"%"<<endl;
                continue;
            }
            
            ibin = h_Ev_toyMC_input->FindBin( Ev );
            prob = h_Ev_normal->GetBinContent(ibin);
            
            Int_t seed_rand = 0;
            
            if( ientry%100==0 )
            {
                seed_generator += 1;
                seed_rand = 2 *seed_generator +1;
                
                gRandom3->SetSeed(seed_rand);
            }
            
            Double_t rand = gRandom3->Uniform(0,1.0);
            if( rand>prob ) continue;
            
            etree->Fill();
        }
#endif//ReactorShapeInToy
        etree->Write();
        roofile_etree->Close();
        
#endif//LoadTree
        
        //////////////////////////////////////////////////////////////////////////
        // The etree file has to have been generated before running this code:
        //////////////////////////////////////////////////////////////////////////
        //    const long MaxEntries = 50;
        
        //    Double_t X_NG[MaxEntries];
        //    Double_t Y_NG[MaxEntries];
        //    Double_t Z_NG[MaxEntries];
        //
        //    Double_t Edep_NG[MaxEntries];
        //    Double_t Etot_NG[MaxEntries];
        //    Double_t EdepEA_NG[MaxEntries];
        //    Double_t EtotE_NG[MaxEntries];
        //    Double_t Etot_Neutron[MaxEntries];
        //    Double_t NeutrinoEnergy[MaxEntries];//Is this too big? Having an array allow us to run the tree reading spectrum only once, instead of once per AD.//cap_Ev
        //    Double_t EdepE_NG[MaxEntries];
        //    Double_t EtotScint_NG[MaxEntries];
        //
        //
        //    Double_t X_Pos[MaxEntries];
        //    Double_t Y_Pos[MaxEntries];
        //    Double_t Z_Pos[MaxEntries];
        //
        //    Double_t Edep_Pos[MaxEntries];
        //    Double_t Edep_PG1[MaxEntries];
        //    Double_t Edep_PG2[MaxEntries];
        //
        //    Double_t EtotScint_Pos[MaxEntries];
        //    Double_t EtotScint_PG[MaxEntries];
        //    Double_t EtotScint_PG1[MaxEntries];
        //    Double_t EtotScint_PG2[MaxEntries];
        //
        //    Double_t EdepE_Pos[MaxEntries];
        //    Double_t EdepE_PG1[MaxEntries];
        //    Double_t EdepE_PG2[MaxEntries];
        //
        //    Double_t EdepAcrylic_Pos[MaxEntries];
        //    Double_t EdepEA_Pos[MaxEntries];
        //    Double_t EdepEA_PG1[MaxEntries];
        //    Double_t EdepEA_PG2[MaxEntries];
        //
        //    Double_t EtotE_Pos[MaxEntries];
        //    Double_t EtotE_PG1[MaxEntries];
        //    Double_t EtotE_PG2[MaxEntries];
        //
        //    Double_t Etot_Pos[MaxEntries];
        //    Double_t Etot_PG1[MaxEntries];
        //    Double_t Etot_PG2[MaxEntries];
#ifndef LoadToyTree

    TFile *roofile_etree_read = new TFile("./Inputs/HInputs/Data/roofile_etree.root", "read");

    TTree *etree_read = (TTree*)roofile_etree_read->Get("etree");
    
    long entries_etree_read = etree_read->GetEntries();
    
    cout<<" ---> filtered entries after neutrino spectrum selection: "<<entries_etree_read<<endl;
    
//        etree_read->SetBranchAddress("cap_Ev",       &cap_Ev);
//        etree_read->SetBranchAddress("cap_Target",   &cap_Target);
//        etree_read->SetBranchAddress("X_Ng",         &X_Ng);
//        etree_read->SetBranchAddress("Y_Ng",         &Y_Ng);
//        etree_read->SetBranchAddress("Z_Ng",         &Z_Ng);
//        etree_read->SetBranchAddress("Edep_Ng",      &Edep_Ng);
//        etree_read->SetBranchAddress("EtotScint_Ng", &EtotScint_Ng);
//        etree_read->SetBranchAddress("Etot_Ng",      &Etot_Ng);
//        etree_read->SetBranchAddress("EdepE_Ng",     &EdepE_Ng);
//        etree_read->SetBranchAddress("EdepEA_Ng",    &EdepEA_Ng);
//        etree_read->SetBranchAddress("EtotE_Ng",     &EtotE_Ng);
//        etree_read->SetBranchAddress("Etot_N",       &Etot_N);
//        etree_read->SetBranchAddress("X_P",          &X_P);
//        etree_read->SetBranchAddress("Y_P",          &Y_P);
//        etree_read->SetBranchAddress("Z_P",          &Z_P);
//        etree_read->SetBranchAddress("Etot_P",       &Etot_P);
//        etree_read->SetBranchAddress("Etot_Pg1",     &Etot_Pg1);
//        etree_read->SetBranchAddress("Etot_Pg2",     &Etot_Pg2);
//        etree_read->SetBranchAddress("EtotScint_P",  &EtotScint_P);
//        etree_read->SetBranchAddress("EtotScint_Pg", &EtotScint_Pg);
//        etree_read->SetBranchAddress("EtotScint_Pg1", &EtotScint_Pg1);
//        etree_read->SetBranchAddress("EtotScint_Pg2", &EtotScint_Pg2);
//        etree_read->SetBranchAddress("Edep_P",        &Edep_P);
//        etree_read->SetBranchAddress("EdepAcrylic_P", &EdepAcrylic_P);
//        etree_read->SetBranchAddress("Edep_Pg1",      &Edep_Pg1);
//        etree_read->SetBranchAddress("Edep_Pg2",      &Edep_Pg2);
//        etree_read->SetBranchAddress("EdepE_P",       &EdepE_P);
//        etree_read->SetBranchAddress("EdepEA_P",      &EdepEA_P);
//        etree_read->SetBranchAddress("EtotE_P",       &EtotE_P);
//        etree_read->SetBranchAddress("EdepE_Pg1",     &EdepE_Pg1);
//        etree_read->SetBranchAddress("EdepEA_Pg1",    &EdepEA_Pg1);
//        etree_read->SetBranchAddress("EtotE_Pg1",     &EtotE_Pg1);
//        etree_read->SetBranchAddress("EdepE_Pg2",     &EdepE_Pg2);
//        etree_read->SetBranchAddress("EdepEA_Pg2",    &EdepEA_Pg2);
//        etree_read->SetBranchAddress("EtotE_Pg2",     &EtotE_Pg2);
        
        etree_read->SetBranchAddress("Ev",           &Ev);
        //etree_read->SetBranchAddress("cap_Target",   &cap_Target);
        //etree_read->SetBranchAddress("IBDvolume",    &IBDvolume);
        etree_read->SetBranchAddress("Tdiff",        &Tdiff);
        etree_read->SetBranchAddress("X_Ng",         &X_Ng);
        etree_read->SetBranchAddress("Y_Ng",         &Y_Ng);
        etree_read->SetBranchAddress("Z_Ng",         &Z_Ng);
        etree_read->SetBranchAddress("EdepScint_Ng", &EdepScint_Ng);
        etree_read->SetBranchAddress("EtotScint_Ng", &EtotScint_Ng);
        //etree_read->SetBranchAddress("Etot_Ng",      &Etot_Ng);
        //etree_read->SetBranchAddress("EdepE_Ng",     &EdepE_Ng);
        //etree_read->SetBranchAddress("EdepEA_Ng",    &EdepEA_Ng);
        //etree_read->SetBranchAddress("EtotE_Ng",     &EtotE_Ng);
        etree_read->SetBranchAddress("Etot_N",       &Etot_N);
        etree_read->SetBranchAddress("EtotScint_N",  &EtotScint_N);
        etree_read->SetBranchAddress("X_P",          &X_P);
        etree_read->SetBranchAddress("Y_P",          &Y_P);
        etree_read->SetBranchAddress("Z_P",          &Z_P);
        etree_read->SetBranchAddress("Etot_P",       &Etot_P);
        //etree_read->SetBranchAddress("Etot_Pg1",     &Etot_Pg1);
        //etree_read->SetBranchAddress("Etot_Pg2",     &Etot_Pg2);
        etree_read->SetBranchAddress("EtotScint_P",  &EtotScint_P);
        //etree_read->SetBranchAddress("EtotScint_Pg", &EtotScint_Pg);
        etree_read->SetBranchAddress("EtotScint_Pg1", &EtotScint_Pg1);
        etree_read->SetBranchAddress("EtotScint_Pg2", &EtotScint_Pg2);
        etree_read->SetBranchAddress("EdepScint_P",   &EdepScint_P);
        etree_read->SetBranchAddress("EdepAcrylic_P", &EdepAcrylic_P);
        etree_read->SetBranchAddress("EdepScint_Pg1", &EdepScint_Pg1);
        etree_read->SetBranchAddress("EdepScint_Pg2", &EdepScint_Pg2);
        //etree_read->SetBranchAddress("EdepE_P",       &EdepE_P);
        //etree_read->SetBranchAddress("EdepEA_P",      &EdepEA_P);
        //etree_read->SetBranchAddress("EtotE_P",       &EtotE_P);
        //etree_read->SetBranchAddress("EdepE_Pg1",     &EdepE_Pg1);
        //etree_read->SetBranchAddress("EdepEA_Pg1",    &EdepEA_Pg1);
        //etree_read->SetBranchAddress("EtotE_Pg1",     &EtotE_Pg1);
        //etree_read->SetBranchAddress("EdepE_Pg2",     &EdepE_Pg2);
        //etree_read->SetBranchAddress("EdepEA_Pg2",    &EdepEA_Pg2);
        //etree_read->SetBranchAddress("EtotE_Pg2",     &EtotE_Pg2);
        
        ///////////////////////////////
        ///////////////////////////////
        ///////////////////////////////
        ///////////////////////////////
#endif
    
//        Double_t LY_E_P;
//        Double_t LY_E_Pg1;
//        Double_t LY_E_Pg2;
//        Double_t LY_E_N;
//        
//        Double_t LY_E_P_Sum;
//        Double_t Opt_E_P_Sum;
//        Double_t FEE_E_P_Sum;
//        Double_t Scale_E_P_Sum;
//        Double_t Res_E_P_Sum;
//        Double_t Eff_E_P_Sum;
//        Double_t Erec_P;

        //        Double_t LY_E_Ng;
        //        Double_t Opt_E_Ng;
        //        Double_t FEE_E_Ng;
        //        Double_t Scale_E_Ng;
        //        Double_t Res_E_Ng;
        //        Double_t Eff_E_Ng;
        //        Double_t Erec_Ng;
        
        double LY_E_P;
        double LY_E_Pg1;
        double LY_E_Pg2;
        double LY_E_N;
        
        double LY_E_P_Sum;
        double Opt_E_P_Sum;
        double FEE_E_P_Sum;
        double Scale_E_P_Sum;
        double Res_E_P_Sum;
        double Erec_P;
        
        double LY_E_Ng = 0;
        double Opt_E_Ng;
        double FEE_E_Ng;
        double Scale_E_Ng;
        double Res_E_Ng;
        double Erec_Ng;

        //Load tree data
        //    for(long ientry=0; ientry<entries_etree_read; ientry++)
        //    {
        //
        //        X_NG[ientry] = X_Ng;
        //        Y_NG[ientry] = Y_Ng;
        //        Z_NG[ientry] = Z_Ng;
        //
        //        Edep_NG[ientry] = Edep_Ng;
        //        Etot_NG[ientry] = Etot_Ng;
        //        EdepEA_NG[ientry] = EdepEA_Ng;
        //        EtotE_NG[ientry] = EtotE_Ng;
        //        Etot_Neutron[ientry]= Etot_N;
        //        NeutrinoEnergy[ientry]=cap_Ev;
        //        EtotScint_NG[ientry] = EtotScint_Ng;
        //        EdepE_NG[ientry] = EdepE_Ng;
        //
        //        X_Pos[ientry] = X_P;
        //        Y_Pos[ientry] = Y_P;
        //        Z_Pos[ientry] = Z_P;
        //
        //        Edep_Pos[ientry] = Edep_P;
        //        Edep_PG1[ientry] = Edep_Pg1;
        //        Edep_PG2[ientry] = Edep_Pg2;
        //
        //        EtotScint_Pos[ientry] =EtotScint_P;
        //        EtotScint_PG[ientry]=EtotScint_Pg;
        //        EtotScint_PG1[ientry]=EtotScint_Pg1;
        //        EtotScint_PG2[ientry]=EtotScint_Pg2;
        //
        //        EdepE_Pos[ientry] = EdepE_P;
        //        EdepE_PG1[ientry]=EdepE_Pg1;
        //        EdepE_PG2[ientry]=EdepE_Pg2;
        //
        //        EdepAcrylic_Pos[ientry]=EdepAcrylic_P;
        //        EdepEA_Pos[ientry]=EdepEA_P;
        //        EdepEA_PG1[ientry]=EdepEA_Pg1;
        //        EdepEA_PG2[ientry]=EdepEA_Pg2;
        //
        //        EtotE_Pos[ientry]=EtotE_P;
        //        EtotE_PG1[ientry]=EtotE_Pg1;
        //        EtotE_PG2[ientry]=EtotE_Pg2;
        //
        //        Etot_Pos[ientry]=Etot_P;
        //        Etot_PG1[ientry]=Etot_Pg1;
        //        Etot_PG2[ientry]=Etot_Pg2;
        //    }
        
        Double_t EendScint_P=0, EendScint_Pg1=0, EendScint_Pg2=0, EendScint_Ng=0;
//        Double_t AdSimple_P=0, AdSimple_Ng=0;
        Double_t usr_opt_attenuation_P=0, usr_pmt_coverage_P=0, usr_opt_attenuation_Ng=0, usr_pmt_coverage_Ng=0;
        Double_t energy_sigma=0,R_average=0;
        //Double_t R2_AA=0, R2_BB=0,
        
        //Efficiency maps:
        TFile *roofile_h2d_ep_ratio2center = new TFile("./Inputs/HInputs/Data/cell_eff/h2d_ep_ratio2center.root", "read");
        TH2D *h2d_Ep_ratio2center = (TH2D*)roofile_h2d_ep_ratio2center->Get("h2d_ep_ratio2center");
        
        TFile *roofile_h2d_ed_ratio2center = new TFile("./Inputs/HInputs/Data/cell_eff/h2d_ed_ratio2center.root", "read");
        TH2D *h2d_Ed_ratio2center = (TH2D*)roofile_h2d_ed_ratio2center->Get("h2d_ed_ratio2center");
    
//        for(int idx=0; idx<MaxColumnNum; idx++)
//        {
//            roostr = TString::Format("hEp_cl_%03d", idx);
//            hEp_cl[idx] = new TH1D(roostr, roostr, hEp_bin, hEp_low, hEp_hgh);
//            
//            roostr = TString::Format("hEd_cl_%03d", idx);
//            hEd_cl[idx] = new TH1D(roostr, roostr, hEd_bin, hEd_low, hEd_hgh);
//        }
//        
//        for(int idx=0; idx<MaxCellNum; idx++)
//        {
//            roostr = TString::Format("hEp_cc_%03d",idx);
//            hEp_cc[idx] = new TH1D(roostr, roostr, hEp_bin, hEp_low, hEp_hgh);
//            
//            roostr = TString::Format("hEd_cc_%03d",idx);
//            hEd_cc[idx] = new TH1D(roostr, roostr, hEd_bin, hEd_low, hEd_hgh);
//            
//            
//            roostr = TString::Format("hEp_cc_clone_%03d",idx);
//            hEp_cc_clone[idx] = new TH1D(roostr, roostr, hEp_bin, hEp_low, hEp_hgh);
//            
//            roostr = TString::Format("hEd_cc_clone_%03d",idx);
//            hEd_cc_clone[idx] = new TH1D(roostr, roostr, hEd_bin, hEd_low, hEd_hgh);
//        }
//        
//        for(int idx=0; idx<MaxColumnNum; idx++)
//        {
//            roostr = TString::Format("hEp_cl_%03d", idx);
    //            hEp_cl[idx] = new TH1D(roostr, roostr, hEp_bin, hEp_low, hEp_hgh);
    //
    //            roostr = TString::Format("hEd_cl_%03d", idx);
    //            hEd_cl[idx] = new TH1D(roostr, roostr, hEd_bin, hEd_low, hEd_hgh);
    //        }
    
    Double_t IntegralEvents[NSystematic][NADs][VolumeX][VolumeY];
    Double_t ADIntegralEvents[NSystematic][NADs];
    
    Double_t distPD=0;
    
    cout.precision(4);
   
    //From Logan NL discussion:
    
//    Since most gammas are 0.51 MeV, I assume that the gammas have all the same uncertainties and are fully correlated.  I assume that the positron is also fully correlated because the gamma curve is derived from the beta.  These assumptions give the following standard deviation due to nonlinearity divided among the three particles
//    
//    sigma_NL = sigma_e + 2*sigma_g
//    
//    Assuming equal contributions from the Boron12 and gamma data to the overall uncertainty and that they correspond directly to sigma_e and sigma_g, gives
//    
//    sigma_e = sigma_NL / 2
//    sigma_g = sigma_NL / 4
//    
//    The above was ignoring electronics.  Assuming equal contributions from both scintillator and electronics gives
//    
//    sigma_e = sigma_NL / (2*sqrt(2))
//    sigma_g = sigma_NL / (4*sqrt(2))
//    
//    for each component (scintillator and FEE), giving two identical pairs of expressions.
//        
//        This result should not be conservative (the only uncorrelated items are the scintillator and FEE components), but I think it is the best estimate of reality.  In addition, we can compare the result with a conservative estimate of
//        
//        sigma_e = sigma_NL / (sqrt(2)*sqrt(2))
//        sigma_g = sigma_NL / (sqrt(8)*sqrt(2))
    
    Double_t Sigma_Electron = 2*TMath::Sqrt(2);
    Double_t Sigma_Gamma = 4*TMath::Sqrt(2);
    
    Int_t VolumeXbin, VolumeYbin;
    Double_t RandomEdepScint_P,RandomEdepScint_Pg1,RandomEdepScint_Pg2;
    
    //Logic to detect to translate cell bins to volume bins (check if the cell is inside the GdLs or the Ls volume)
    VolumeYbin = 1;
    
    Double_t BinWidth = (FinalVisibleEnergy-InitialVisibleEnergy)/(1.*MatrixBins);
    
#ifdef LoadToyTree
    TFile *roofile_etree_read = new TFile("./Inputs/HInputs/Data/Nominalroofile_toy.root", "read");
    TTree *etree_read = (TTree*)roofile_etree_read->Get("toy");
    
    //        etree_read->SetBranchAddress("cap_Ev",       &cap_Ev);
    //        etree_read->SetBranchAddress("cap_Target",   &cap_Target);
    //        etree_read->SetBranchAddress("X_Ng",         &X_Ng);
    //        etree_read->SetBranchAddress("Y_Ng",         &Y_Ng);
    //        etree_read->SetBranchAddress("Z_Ng",         &Z_Ng);
    //        etree_read->SetBranchAddress("Edep_Ng",      &Edep_Ng);
    //        etree_read->SetBranchAddress("EtotScint_Ng", &EtotScint_Ng);
    //        etree_read->SetBranchAddress("Etot_Ng",      &Etot_Ng);
    //        etree_read->SetBranchAddress("EdepE_Ng",     &EdepE_Ng);
    //        etree_read->SetBranchAddress("EdepEA_Ng",    &EdepEA_Ng);
    //        etree_read->SetBranchAddress("EtotE_Ng",     &EtotE_Ng);
    //        etree_read->SetBranchAddress("Etot_N",       &Etot_N);
    //        etree_read->SetBranchAddress("X_P",          &X_P);
    //        etree_read->SetBranchAddress("Y_P",          &Y_P);
    //        etree_read->SetBranchAddress("Z_P",          &Z_P);
    //        etree_read->SetBranchAddress("Etot_P",       &Etot_P);
    //        etree_read->SetBranchAddress("Etot_Pg1",     &Etot_Pg1);
    //        etree_read->SetBranchAddress("Etot_Pg2",     &Etot_Pg2);
    //        etree_read->SetBranchAddress("EtotScint_P",  &EtotScint_P);
    //        etree_read->SetBranchAddress("EtotScint_Pg", &EtotScint_Pg);
    //        etree_read->SetBranchAddress("EtotScint_Pg1", &EtotScint_Pg1);
    //        etree_read->SetBranchAddress("EtotScint_Pg2", &EtotScint_Pg2);
    //        etree_read->SetBranchAddress("Edep_P",        &Edep_P);
    //        etree_read->SetBranchAddress("EdepAcrylic_P", &EdepAcrylic_P);
    //        etree_read->SetBranchAddress("Edep_Pg1",      &Edep_Pg1);
    //        etree_read->SetBranchAddress("Edep_Pg2",      &Edep_Pg2);
    //        etree_read->SetBranchAddress("EdepE_P",       &EdepE_P);
    //        etree_read->SetBranchAddress("EdepEA_P",      &EdepEA_P);
    //        etree_read->SetBranchAddress("EtotE_P",       &EtotE_P);
    //        etree_read->SetBranchAddress("EdepE_Pg1",     &EdepE_Pg1);
    //        etree_read->SetBranchAddress("EdepEA_Pg1",    &EdepEA_Pg1);
    //        etree_read->SetBranchAddress("EtotE_Pg1",     &EtotE_Pg1);
    //        etree_read->SetBranchAddress("EdepE_Pg2",     &EdepE_Pg2);
    //        etree_read->SetBranchAddress("EdepEA_Pg2",    &EdepEA_Pg2);
    //        etree_read->SetBranchAddress("EtotE_Pg2",     &EtotE_Pg2);
    
    etree_read->SetBranchAddress("Ev",           &Ev);
    //etree_read->SetBranchAddress("cap_Target",   &cap_Target);
    //        etree_read->SetBranchAddress("IBDvolume",    &IBDvolume);
    etree_read->SetBranchAddress("Tdiff",        &Tdiff);
    etree_read->SetBranchAddress("X_Ng",         &X_Ng);
    etree_read->SetBranchAddress("Y_Ng",         &Y_Ng);
    etree_read->SetBranchAddress("Z_Ng",         &Z_Ng);
    etree_read->SetBranchAddress("EdepScint_Ng", &EdepScint_Ng);
    etree_read->SetBranchAddress("EtotScint_Ng", &EtotScint_Ng);
    //etree_read->SetBranchAddress("Etot_Ng",      &Etot_Ng);
    //etree_read->SetBranchAddress("EdepE_Ng",     &EdepE_Ng);
    //etree_read->SetBranchAddress("EdepEA_Ng",    &EdepEA_Ng);
    //etree_read->SetBranchAddress("EtotE_Ng",     &EtotE_Ng);
    etree_read->SetBranchAddress("Etot_N",       &Etot_N);
    etree_read->SetBranchAddress("EtotScint_N",  &EtotScint_N);
    etree_read->SetBranchAddress("X_P",          &X_P);
    etree_read->SetBranchAddress("Y_P",          &Y_P);
    etree_read->SetBranchAddress("Z_P",          &Z_P);
    etree_read->SetBranchAddress("Etot_P",       &Etot_P);
    //etree_read->SetBranchAddress("Etot_Pg1",     &Etot_Pg1);
    //etree_read->SetBranchAddress("Etot_Pg2",     &Etot_Pg2);
    etree_read->SetBranchAddress("EtotScint_P",  &EtotScint_P);
    //etree_read->SetBranchAddress("EtotScint_Pg", &EtotScint_Pg);
    etree_read->SetBranchAddress("EtotScint_Pg1", &EtotScint_Pg1);
    etree_read->SetBranchAddress("EtotScint_Pg2", &EtotScint_Pg2);
    etree_read->SetBranchAddress("EdepScint_P",   &EdepScint_P);
    etree_read->SetBranchAddress("EdepAcrylic_P", &EdepAcrylic_P);
    etree_read->SetBranchAddress("EdepScint_Pg1", &EdepScint_Pg1);
    etree_read->SetBranchAddress("EdepScint_Pg2", &EdepScint_Pg2);
    //etree_read->SetBranchAddress("EdepE_P",       &EdepE_P);
    //etree_read->SetBranchAddress("EdepEA_P",      &EdepEA_P);
    //etree_read->SetBranchAddress("EtotE_P",       &EtotE_P);
    //etree_read->SetBranchAddress("EdepE_Pg1",     &EdepE_Pg1);
    //etree_read->SetBranchAddress("EdepEA_Pg1",    &EdepEA_Pg1);
    //etree_read->SetBranchAddress("EtotE_Pg1",     &EtotE_Pg1);
    //etree_read->SetBranchAddress("EdepE_Pg2",     &EdepE_Pg2);
    //etree_read->SetBranchAddress("EdepEA_Pg2",    &EdepEA_Pg2);
    //etree_read->SetBranchAddress("EtotE_Pg2",     &EtotE_Pg2);
    
    ///////////////////////////////
    ///////////////////////////////
    ///////////////////////////////
    ///////////////////////////////

#endif
    
    //Process tree
    for(Int_t Systematic = InitialSystematic; Systematic<NSystematic;Systematic++)
    {
        long entries_etree_read = etree_read->GetEntries();
        
        cout<<" ---> filtered entries after neutrino spectrum selection: "<<entries_etree_read<<endl;

        RelativeEnergyScaleMatrix = 0;
        IAVMatrix = 0;
        OAVMatrix = 0;
        NLMatrix = 0;
        ResolutionMatrix = 0;
        //                EfficiencyMatrix = 0;
        AllMatrix = 0;
        
        switch (Systematic)
        {
            case 0:
                SystematicS = "Nominal";
                //std::cout << " Nominal Matrix " << std::endl;
                break;
            case 1:
                IAVMatrix = 1;
                SystematicS = "IAV";
                //std::cout << " IAV Matrix " << std::endl;
                break;
            case 2:
                OAVMatrix = 1;
                SystematicS = "OAV";
                //std::cout << " OAV Matrix " << std::endl;
                break;
            case 3:
                NLMatrix = 1;
                SystematicS = "NL";
                //std::cout << " NL Matrix " << std::endl;
                break;
            case 4:
                ResolutionMatrix = 1;
                SystematicS = "Resolution";
                //std::cout << " Reso Matrix " << std::endl;
                break;
            case 5:
                RelativeEnergyScaleMatrix = 1;
                SystematicS = "RelativeEnergyScale";
                //std::cout << " RelativeEnergyScale Matrix " << std::endl;
                break;
            case 6:
                AllMatrix = 1;
                SystematicS = "AllDetectorSystematics";
                //std::cout << " All Detector Systematics Matrix " << std::endl;
                break;
                //                    case 6:
                //                        EfficiencyMatrix = 1;
                //                        //std::cout << " Efficiency Matrix " << std::endl;
                //                        break;
            default:
                std::cout << " DEFAULT SHOULDN'T HAPPEN" << std::endl;
                exit(EXIT_FAILURE);
                break;
        }
        
        TTree *toy;
        TFile *roofile_toy;
#ifdef SaveTree
        
        NADs = 1;//Save only 1 AD in the tree
        
//        Int_t pcell = 0;
//        Int_t dcell = 0;
        if(Systematic==0)//Save only Nominal Tree, if we want to save other systematics, then comment this line after every #ifdef SaveTree
        {
            
        roofile_toy = new TFile(("./Inputs/HInputs/Data/"+SystematicS+"roofile_toy.root").c_str(), "recreate");

        toy = new TTree("toy", "toyMC result");
        
        toy->Branch("Ev",           &Ev,           "Ev/D");
        //toy->Branch("cap_Target",   &cap_Target,   "cap_Target/D");
       // toy->Branch("IBDvolume",    &IBDvolume,    "IBDvolume/I");
       toy->Branch("Tdiff",        &Tdiff,        "Tdiff/D");
//        toy->Branch("pcell",        &pcell,        "pcell/I");
//        toy->Branch("dcell",        &dcell,        "dcell/I");
        toy->Branch("X_Ng",         &X_Ng,         "X_Ng/D");
        toy->Branch("Y_Ng",         &Y_Ng,         "Y_Ng/D");
        toy->Branch("Z_Ng",         &Z_Ng,         "Z_Ng/D");
        toy->Branch("EdepScint_Ng", &EdepScint_Ng, "EdepScint_Ng/D");
        toy->Branch("EtotScint_Ng", &EtotScint_Ng, "EtotScint_Ng/D");
        //toy->Branch("Etot_Ng",      &Etot_Ng,      "Etot_Ng/D");
        //toy->Branch("EdepE_Ng",     &EdepE_Ng,     "EdepE_Ng/D");
        //toy->Branch("EdepEA_Ng",    &EdepEA_Ng,    "EdepEA_Ng/D");
        //toy->Branch("EtotE_Ng",     &EtotE_Ng,     "EtotE_Ng/D");
        toy->Branch("Etot_N",       &Etot_N,       "Etot_N/D");
        toy->Branch("EtotScint_N",  &EtotScint_N,  "EtotScint_N/D");
        toy->Branch("X_P",          &X_P,          "X_P/D");
        toy->Branch("Y_P",          &Y_P,          "Y_P/D");
        toy->Branch("Z_P",          &Z_P,          "Z_P/D");
        toy->Branch("Etot_P",       &Etot_P,       "Etot_P/D");
        //toy->Branch("Etot_Pg1",     &Etot_Pg1,     "Etot_Pg1/D");
        //toy->Branch("Etot_Pg2",     &Etot_Pg2,     "Etot_Pg2/D");
        toy->Branch("EtotScint_P",  &EtotScint_P,  "EtotScint_P/D");
        //toy->Branch("EtotScint_Pg", &EtotScint_Pg, "EtotScint_Pg/D");
        toy->Branch("EtotScint_Pg1", &EtotScint_Pg1, "EtotScint_Pg1/D");
        toy->Branch("EtotScint_Pg2", &EtotScint_Pg2, "EtotScint_Pg2/D");
        toy->Branch("EdepScint_P",   &EdepScint_P,   "EdepScint_P/D");
        toy->Branch("EdepAcrylic_P", &EdepAcrylic_P, "EdepAcrylic_P/D");
        toy->Branch("EdepScint_Pg1", &EdepScint_Pg1, "EdepScint_Pg1/D");
        toy->Branch("EdepScint_Pg2", &EdepScint_Pg2, "EdepScint_Pg2/D");
        //toy->Branch("EdepE_P",       &EdepE_P,       "EdepE_P/D");
        //toy->Branch("EdepEA_P",      &EdepEA_P,      "EdepEA_P/D");
        //toy->Branch("EtotE_P",       &EtotE_P,       "EtotE_P/D");
        //toy->Branch("EdepE_Pg1",     &EdepE_Pg1,     "EdepE_Pg1/D");
        //toy->Branch("EdepEA_Pg1",    &EdepEA_Pg1,    "EdepEA_Pg1/D");
        //toy->Branch("EtotE_Pg1",     &EtotE_Pg1,     "EtotE_Pg1/D");
        //toy->Branch("EdepE_Pg2",     &EdepE_Pg2,     "EdepE_Pg2/D");
        //toy->Branch("EdepEA_Pg2",    &EdepEA_Pg2,    "EdepEA_Pg2/D");
        //toy->Branch("EtotE_Pg2",     &EtotE_Pg2,     "EtotE_Pg2/D");
        
        //toy->Branch("LY_E_P",        &LY_E_P,        "LY_E_P/D");
        //toy->Branch("LY_E_Pg1",      &LY_E_Pg1,      "LY_E_Pg1/D");
        //toy->Branch("LY_E_Pg2",      &LY_E_Pg2,      "LY_E_Pg2/D");
        //toy->Branch("LY_E_N",        &LY_E_N,        "LY_E_N/D");
        //toy->Branch("LY_E_P_Sum",    &LY_E_P_Sum,    "LY_E_P_Sum/D");
        //toy->Branch("Opt_E_P_Sum",   &Opt_E_P_Sum,   "Opt_E_P_Sum/D");
        //toy->Branch("FEE_E_P_Sum",   &FEE_E_P_Sum,   "FEE_E_P_Sum/D");
        //toy->Branch("Scale_E_P_Sum", &Scale_E_P_Sum, "Scale_E_P_Sum/D");
//        toy->Branch("Res_E_P_Sum",   &Res_E_P_Sum,   "Res_E_P_Sum/D");
//        toy->Branch("Erec_P",   &Erec_P,   "Erec_P/D");
        
        //toy->Branch("LY_E_Ng",    &LY_E_Ng,    "LY_E_Ng/D");
        //toy->Branch("Opt_E_Ng",   &Opt_E_Ng,   "Opt_E_Ng/D");
        //toy->Branch("FEE_E_Ng",   &FEE_E_Ng,   "FEE_E_Ng/D");
        //toy->Branch("Scale_E_Ng", &Scale_E_Ng, "Scale_E_Ng/D");
//        toy->Branch("Res_E_Ng",   &Res_E_Ng,   "Res_E_Ng/D");
//        toy->Branch("Erec_Ng",   &Erec_Ng,   "Erec_Ng/D");
        }
        
#endif
        
        seed_generator_uncorr = 2863311530; //=(2*4294967295/3) for uncorrelated systematics, chose 2/3*(maxseed) to make it different to 'seed_generator' as a precaution so I don't use the same seeds.
       
        Int_t MaxTrueEnegyLoop = 1;
        
        if(!NormalMode)
        {
           MaxTrueEnegyLoop = MatrixBins;
        }
        for(Int_t AD = 0; AD<NADs;AD++)
        {
            seed_generator = 1;//this way all matrices will be random but have a common nominal spectrum
            seed_generator_corr = 1431655765; //=(4294967295/3) for correlated systematics, chose maxseed/3 to make it different to 'seed_generator' as a precaution so I don't use the same seeds.
           
            long NumberOfEntries = 0;
            bool firstB = 1;
            
            Int_t seed_rand = 0;
            Int_t seed_corr =0;
            Int_t seed_uncorr = 0;
            
            for(Int_t TrueEnergy = 0; TrueEnergy<MaxTrueEnegyLoop; TrueEnergy++)
            {
                if(TrueEnergy<0)
                {
                    TrueEnergy=0;
                }
                Double_t MinTrueEnergy = 1.*TrueEnergy*BinWidth;
                Double_t MaxTrueEnergy = (TrueEnergy+1)*1.*BinWidth;
                
//                std::cout << " MIN ENERGY : " << MinTrueEnergy << " AND MAX ENERGY : " << MaxTrueEnergy << std::endl;
                TEventList* ListEvents;

                if(!NormalMode)
                {
                    //Need a method to load a list of entries from the tree with the energy selected
                    
                    etree_read->SetEventList(0);//Reset the lists
                    
                    etree_read->Draw(">>List",Form("Ev >= %f && Ev < %f",MinTrueEnergy,MaxTrueEnergy));
                   
                    ListEvents = (TEventList*)gDirectory->Get("List");
                    
                    etree_read->SetEventList(ListEvents);
                    if(firstB)
                    {
                        NumberOfEntries = ListEvents->GetN();//It happens that the first non-empty bin has 33993 events and that is the absolute minimum, let's make all the energy bins to have the same number of events in order to minimize the statistical fluctuations due to random sampling. Everytime I skip an event I will add 1 to the loop max limit so the events will keep being the same even under random systematic behaviour.
                        if(NumberOfEntries>0)
                        {
                            firstB=0;
                        }
                    }
                    std::cout << " LIST ENTRIES: " << NumberOfEntries << " FOR ENERGIES " << MinTrueEnergy << " TO " <<  MaxTrueEnergy << std::endl;
                }
                else
                {
                    NumberOfEntries = entries_etree_read;
                    std::cout << " TREE ENTRIES: " << NumberOfEntries << std::endl;
                }
                
                for(long ientry=0; ientry<NumberOfEntries; ientry++)
                {
                    long index = 0;
                    
                    if(!NormalMode)
                    {
                        index = ListEvents->GetEntry(ientry);
                    }
                    else
                    {
                        index = ientry;
                    }
                    
                    etree_read->GetEntry(index);
                    
                    if(ientry%20000==0)
                    {
                        cout<<" ---> processing response matrix in AD" << AD << " " <<ientry*100./entries_etree_read<<"%"<<endl;
                    }
                    
                    if( ientry%100==0 )
                    {
                        seed_generator += 1;
                        seed_generator_corr +=1;
                        seed_generator_uncorr+=1;
                        
                        seed_rand = 2 *seed_generator +1;
                        seed_corr = 2 *seed_generator_corr +1;
                        seed_uncorr = 2 *seed_generator_uncorr +1;
                        //
                        //                    std::cout << "uncorrelated seed: " << seed_generator_uncorr << "- AD: " << AD << "- entry: " << ientry << std::endl;
                        //                    std::cout << "correlated seed: " << seed_generator_corr << "- AD: " << AD << "- entry: " << ientry << std::endl;
                        
                        gRandom3->SetSeed(seed_rand);
                        RandomSysUncorr->SetSeed(seed_uncorr);
                        RandomSysCorr->SetSeed(seed_corr);
                    }
                    
                    //each event has a different systematic error:
                    
                    //Nominal values
                    
                    GausRelative = 0;
                    GausIAV= 0;
                    GausOAV = 0;
                    
                    for (Int_t ierr = 0; ierr < unified_nl_pars; ierr++)
                    {
                        GausNL[ierr]= 0;
                    }
                    GausReso=0;
                    GausResoCorr=0;
                    //                GausEff= 0;
                    //                if(EfficiencyMatrix)
                    //                {
                    //                    GausEff = RandomSysUncorr->Gaus(0,1);
                    //                }
                    if(RelativeEnergyScaleMatrix)
                    {
                        GausRelative = RandomSysUncorr->Gaus(0,1);
                    }
                    else if(IAVMatrix)
                    {
                        GausIAV = RandomSysUncorr->Gaus(0,1);
                    }
                    else if(OAVMatrix)
                    {
                        GausOAV = RandomSysUncorr->Gaus(0,1);;
                    }
                    else if(NLMatrix)
                    {
                        for (Int_t ierr = 0; ierr < unified_nl_pars; ierr++)
                        {
                            GausNL[ierr] = RandomSysUncorr->Gaus(0,1);
                        }
                    }
                    else if(ResolutionMatrix)
                    {
                        GausReso = RandomSysUncorr->Gaus(0,1);
                        GausResoCorr = RandomSysCorr->Gaus(0,1);
                    }
                    else if(AllMatrix)
                    {
                        GausRelative = RandomSysUncorr->Gaus(0,1);
                        
                        GausIAV = RandomSysUncorr->Gaus(0,1);
                        
                        GausOAV =  RandomSysUncorr->Gaus(0,1);
                        
                        for (Int_t ierr = 0; ierr < unified_nl_pars; ierr++)
                        {
                            GausNL[ierr] = RandomSysUncorr->Gaus(0,1);
                        }
                        
                        GausReso = RandomSysUncorr->Gaus(0,1);
                        
                        GausResoCorr = RandomSysCorr->Gaus(0,1);
                        
                        //                    GausEff = RandomSysUncorr->Gaus(0,1);
                        
                    }
                    
                    //factor that is multiplied by the random error so it's added to the nominal by using:  Value = NominalValue + gRandom3->Gauss(0,1)*1sigmaError;
#ifdef DisplayInfo
                    if(ientry%1000000==0)//show only a few
                    {
                        std::cout << "IAV uncorrelated should be the different for different ads and same event: " << GausIAV << "ad: " << AD << " - event: " << ientry <<  std::endl;
                        
                        std::cout << "OAV uncorrelated should be the different for different ads and same event: " << GausOAV << "ad: " << AD << " - event: " << ientry <<  std::endl;
                        
                        std::cout << "gaus reso correlated should be the same for different ads and same event: " << GausResoCorr << "ad: " << AD << " - event: " << ientry <<  std::endl;
                        std::cout << "gaus reso uncorrelated should be different for different ads and same event" << GausReso << "ad: " << AD << " - event: " << ientry <<  std::endl;
                        
                        std::cout << "NL uncorrelated should be different for different ads and same event" << GausNL[0] << "ad: " << AD << " - event: " << ientry <<  std::endl;
                        
                        std::cout << "Relative Energy uncorrelated should be different for different ads and same event" << GausRelative << "ad: " << AD << " - event: " << ientry <<  std::endl;
                        
                        //                    std::cout << "Efficiency Energy uncorrelated should be different for different ads and same event" << GausEff << "ad: " << AD << " - event: " << ientry <<  std::endl;
                        //
                    }
#endif
                    //to check the random numbers are properly generated:
                    //            std::cout << " AD: " << AD << " Gaus RELATIVE: " << GausRelative[AD] << " , Gaus IAV: " << GausIAV[AD] << " , Gaus NL: " << GausNL[AD] << " , Gaus RESO: " << GausReso[AD] << " , Gaus EFF: " << GausEff[AD] << std::endl;
                    
                    ///////////////////////////// LY
                    ///////////////////////////// LY
                    /*
                     /// case01: default ?
                     /// prompt
                     LY_E_P = EtotScint_P * graph_electron_LY->Eval( EtotScint_P, 0, "s" )
                     - (EtotScint_P-Edep_P) * graph_electron_LY->Eval( EtotScint_P-Edep_P, 0, "s" )
                     - (EtotE_P-EdepE_P) * graph_gamma_LY->Eval( EtotE_P-EdepE_P, 0, "s" );
                     
                     LY_E_Pg1 = EtotScint_Pg1 * graph_gamma_LY->Eval( EtotScint_Pg1, 0, "s" )
                     - (EtotScint_Pg1-EdepE_Pg1) * graph_gamma_LY->Eval( EtotScint_Pg1-EdepE_Pg1, 0, "s" );
                     
                     LY_E_Pg2 = EtotScint_Pg2 * graph_gamma_LY->Eval( EtotScint_Pg2, 0, "s" )
                     - (EtotScint_Pg2-EdepE_Pg2) * graph_gamma_LY->Eval( EtotScint_Pg2-EdepE_Pg2, 0, "s" );
                     
                     LY_E_P_Sum = LY_E_P + LY_E_Pg1 + LY_E_Pg2;
                     
                     /// delayed
                     LY_E_Ng = EtotScint_Ng * graph_gamma_LY->Eval( EtotScint_Ng, 0, "s" )
                     - (EtotScint_Ng-EdepE_Ng) * graph_gamma_LY->Eval( EtotScint_Ng-EdepE_Ng, 0, "s" );
                     */
                    /*
                     /// case02: no leakage
                     /// prompt
                     LY_E_P = EtotScint_P * graph_electron_LY->Eval( EtotScint_P, 0, "s" );
                     
                     LY_E_Pg1 = EtotScint_Pg1 * graph_gamma_LY->Eval( EtotScint_Pg1, 0, "s" );
                     
                     LY_E_Pg2 = EtotScint_Pg2 * graph_gamma_LY->Eval( EtotScint_Pg2, 0, "s" );
                     
                     LY_E_P_Sum = LY_E_P + LY_E_Pg1 + LY_E_Pg2;
                     
                     /// delayed
                     LY_E_Ng = EtotScint_Ng * graph_gamma_LY->Eval( EtotScint_Ng, 0, "s" );
                     */
                    
                    /*
                     /// case03: edit
                     /// prompt
                     LY_E_P = EtotScint_P * graph_electron_LY->Eval( EtotScint_P, 0, "s" )
                     - (EtotE_P-EdepE_P) * graph_gamma_LY->Eval( EtotE_P-EdepE_P, 0, "s" );
                     
                     LY_E_Pg1 = EtotScint_Pg1 * graph_gamma_LY->Eval( EtotScint_Pg1, 0, "s" )
                     - (EtotE_Pg1-EdepE_Pg1) * graph_electron_LY->Eval( EtotE_Pg1-EdepE_Pg1, 0, "s" );
                     
                     LY_E_Pg2 = EtotScint_Pg2 * graph_gamma_LY->Eval( EtotScint_Pg2, 0, "s" )
                     - (EtotE_Pg2-EdepE_Pg2) * graph_electron_LY->Eval( EtotE_Pg1-EdepE_Pg1, 0, "s" );
                     
                     LY_E_P_Sum = LY_E_P + LY_E_Pg1 + LY_E_Pg2;
                     
                     /// delayed
                     LY_E_Ng = EtotScint_Ng * graph_gamma_LY->Eval( EtotScint_Ng, 0, "s" )
                     - (EtotE_Ng-EdepE_Ng) * graph_electron_LY->Eval( EtotE_Ng-EdepE_Ng, 0, "s" );
                     */
                    
                    //Iav uncertainty considered for prompt signal. 0.1% bin-to-bin uncorrelated.
                    RandomEdepScint_P = EdepScint_P*(1+GausIAV*IAVNominalError);
                    RandomEdepScint_Pg1 = EdepScint_Pg1*(1+GausIAV*IAVNominalError);
                    RandomEdepScint_Pg2 = EdepScint_Pg2*(1+GausIAV*IAVNominalError);
                    //here I am applying the same variation to the three particles at once, one variation for each event.
                    
                    if( GausIAV*IAVNominalError < -1)
                    {
                        continue;//Very unlikely, but just in case.
                    }
                    
                    /// case02: default
                    EendScint_P = EtotScint_P-RandomEdepScint_P;//Edep_P;
                    EendScint_Pg1 = EtotScint_Pg1-RandomEdepScint_Pg1;//Edep_Pg1;
                    EendScint_Pg2 = EtotScint_Pg2-RandomEdepScint_Pg2;//Edep_Pg2;
                    //EendScint_Ng = EtotScint_Ng-EdepScint_Ng;//Edep_Ng;
                    
                    //                P_Scint = EtotScint_P-RandomEdep_P;
                    //                Pg1_Scint = EtotScint_Pg1-RandomEdep_Pg1;
                    //                Pg2_Scint = EtotScint_Pg2-RandomEdep_Pg2;
                    //                Ng_Scint = EtotScint_Ng-Edep_Ng;
                    
                    /// prompt
                    //                LY_E_P = EtotScint_P * graph_electron_LY->Eval( EtotScint_P, 0, "s" ) - P_Scint * graph_electron_LY->Eval( P_Scint, 0, "s" );
                    //                LY_E_Pg1 = EtotScint_Pg1 * graph_gamma_LY->Eval( EtotScint_Pg1, 0, "s" ) - Pg1_Scint * graph_gamma_LY->Eval( Pg1_Scint, 0, "s" );
                    //                LY_E_Pg2 = EtotScint_Pg2 * graph_gamma_LY->Eval( EtotScint_Pg2, 0, "s" ) - Pg2_Scint * graph_gamma_LY->Eval( Pg2_Scint, 0, "s" );
                    //                LY_E_N = Etot_N * ( 0.186 + exp(-1.142-126.4*Etot_N) + exp(-1.217-17.90*Etot_N) ) * graph_electron_LY->Eval( 0.2, 0, "s" );
                    LY_E_P = 0;
                    LY_E_Pg1 = 0;
                    LY_E_Pg2 = 0;
                    LY_E_N = 0;
                    
                    LY_E_P += EtotScint_P * graph_electron_LY->Eval( EtotScint_P, 0, "s" )*(1+ (ErrorBand->Interpolate(EtotScint_P)*GausNL[0]/Sigma_Electron));
                    
                    if( EendScint_P>0.000001 ) // calculate leakage effect only when significant
                    {
                        LY_E_P -= EendScint_P * graph_electron_LY->Eval( EendScint_P, 0, "s" )*(1+ (ErrorBand->Interpolate(EendScint_P)*GausNL[0]/Sigma_Electron));
                    }
                    
                    LY_E_Pg1 = EtotScint_Pg1 * graph_gamma_LY->Eval( EtotScint_Pg1, 0, "s" )*(1+ (ErrorBand->Interpolate(EtotScint_Pg1)*GausNL[0]/Sigma_Gamma));
                    
                    if( EendScint_Pg1>0.000001 ) // calculate leakage effect only when significant
                    {
                        LY_E_Pg1 -= EendScint_Pg1 * graph_gamma_LY->Eval( EendScint_Pg1, 0, "s" )*(1+ (ErrorBand->Interpolate(EendScint_Pg1)*GausNL[0]/Sigma_Gamma));
                    }
                    
                    LY_E_Pg2 = EtotScint_Pg2 * graph_gamma_LY->Eval( EtotScint_Pg2, 0, "s" )*(1+ (ErrorBand->Interpolate(EtotScint_Pg2)*GausNL[0]/Sigma_Gamma));
                    
                    if( EendScint_Pg2>0.000001 ) // calculate leakage effect only when significant
                    {
                        LY_E_Pg2 -= EendScint_Pg2 * graph_gamma_LY->Eval( EendScint_Pg2, 0, "s" )*(1+ (ErrorBand->Interpolate(EendScint_Pg2)*GausNL[0]/Sigma_Gamma));
                    }
                    
                    /// delayed - loop over gammas
                    
                    LY_E_Ng = 0;
                    for( int k=0; k<10; k++ )
                    {
                        if( EtotScint_Ng[k]<=0 )
                        {
                            continue;  // most array elements will be zero
                        }
                        
                        LY_E_Ng += EtotScint_Ng[k] * graph_gamma_LY->Eval( EtotScint_Ng[k], 0, "s" )*(1+ (ErrorBand->Interpolate(EtotScint_Ng[k])*GausNL[0]/Sigma_Gamma));
                        
                        EendScint_Ng = EtotScint_Ng[k]-EdepScint_Ng[k];//Edep_Ng[k];
                        
                        if( EendScint_Ng>0.000001 ) // calculate leakage effect only when significant
                        {
                            LY_E_Ng -= EendScint_Ng * graph_gamma_LY->Eval( EendScint_Ng, 0, "s" )*(1+ (ErrorBand->Interpolate(EendScint_Ng)*GausNL[0]/Sigma_Gamma));
                        }
                        
                        //if( EendScint_Ng<-0.000001 ) {
                        //  cout<<" ERROR: EtotScint_Ng["<< k <<"] - EdepScint_Ng["<< k <<"] = "<< EendScint_Ng <<".  Skipped gamma "<< k <<endl;
                        //  continue;
                        //}
                    }
                    
                    LY_E_N = EtotScint_N * ( 0.186 + exp(-1.142-126.4*EtotScint_N) + exp(-1.217-17.90*EtotScint_N) ) * graph_electron_LY->Eval( 0.2, 0, "s" )*(1+ (ErrorBand->Interpolate(0.2)*GausNL[0]/Sigma_Electron));
                    
                    
                    //               LY_E_Ng = EtotScint_Ng * graph_gamma_LY->Eval( EtotScint_Ng, 0, "s" ) - Ng_Scint * graph_gamma_LY->Eval( Ng_Scint, 0, "s" );
                    
                    //Random NL: // In each event the NL curve is varied //if GausNL = 0 the energy would be the nominal NL function, otherwise energy = nl_nominal_energy(1+Σerrors)
                    //We don't have pull curves for each term (electron, gamma and electronic), so far use the nominal error in each of the curves and divide by square root of 3 asigning to each particle a proportial side of the band (e+ = 2*E(gamma) + e-)/TMath::Sqrt(3)
                    
                    
                    // This code was design for marginal curves (unified_nl_pars dependency)
                    //
                    //                for (Int_t ierr = 0; ierr < unified_nl_pars; ierr++)
                    //                {
                    //                    LY_E_P += EtotScint_P * graph_electron_LY->Eval( EtotScint_P, 0, "s" )*(1+ (ErrorBand->Interpolate(EtotScint_P)*GausNL[ierr]/TMath::Sqrt(3)));
                    ////                    - P_Scint * graph_electron_LY->Eval( P_Scint, 0, "s" )*(1+ (ErrorBand->Interpolate(P_Scint)*GausNL[ierr]/TMath::Sqrt(3)));
                    //
                    //                    if( EendScint_P>0.000001 ) // calculate leakage effect only when significant
                    //                    {
                    //                        LY_E_P -= EendScint_P * graph_electron_LY->Eval( EendScint_P, 0, "s" )*(1+ (ErrorBand->Interpolate(EendScint_P)*GausNL[ierr]/TMath::Sqrt(3)));
                    //                    }
                    //
                    //                    LY_E_Pg1 += EtotScint_Pg1 * graph_gamma_LY->Eval( EtotScint_Pg1, 0, "s" )*(1+ (ErrorBand->Interpolate(EtotScint_Pg1)*GausNL[ierr]/TMath::Sqrt(3)))
                    //                    - Pg1_Scint * graph_gamma_LY->Eval( Pg1_Scint, 0, "s" )*(1+ (ErrorBand->Interpolate(Pg1_Scint)*GausNL[ierr]/TMath::Sqrt(3)));
                    //
                    //                    LY_E_Pg2 += EtotScint_Pg2 * graph_gamma_LY->Eval( EtotScint_Pg2, 0, "s" )*(1+ (ErrorBand->Interpolate(EtotScint_Pg2)*GausNL[ierr]/TMath::Sqrt(3)))
                    //                    - Pg2_Scint * graph_gamma_LY->Eval( Pg2_Scint, 0, "s" )*(1+ (ErrorBand->Interpolate(Pg2_Scint)*GausNL[ierr]/TMath::Sqrt(3)));
                    //
                    //                    // ------------------------ //
                    ////                    /// neutron - determined from NuWa for 0-0.18 MeV (does not consider leakage)
                    //                   LY_E_N += EtotScint_N * ( 0.186 + exp(-1.142-126.4*EtotScint_N) + exp(-1.217-17.90*EtotScint_N) ) * graph_electron_LY->Eval( 0.2, 0, "s" )*(1+ (ErrorBand->Interpolate(0.2)*GausNL[ierr]/TMath::Sqrt(3)));
                    //
                    ////         Uncomment to use marginal curves, also remove Sqrt(3);
                    //
                    ////                    LY_E_P += EtotScint_P * (graph_electron_LY->Eval( EtotScint_P, 0, "s" ) - g_unified_positron_nl_pulls[ierr]->Eval( EtotScint_P, 0, "s" ))*(GausNL[ierr]/TMath::Sqrt(3)
                    ////                    - P_Scint * (graph_electron_LY->Eval( P_Scint, 0, "s" ) - g_unified_positron_nl_pulls[ierr]->Eval( P_Scint, 0, "s" ))*GausNL[ierr]/TMath::Sqrt(3);
                    ////
                    ////
                    ////                    LY_E_Pg1 += EtotScint_Pg1 * (graph_gamma_LY->Eval( EtotScint_Pg1, 0, "s" ) - g_unified_positron_nl_pulls[ierr]->Eval( EtotScint_Pg1, 0, "s" ))*GausNL[ierr]/TMath::Sqrt(3)
                    ////                    - Pg1_Scint * (graph_gamma_LY->Eval( Pg1_Scint, 0, "s" ) - g_unified_positron_nl_pulls[ierr]->Eval( Pg1_Scint, 0, "s" ))*GausNL[ierr]/TMath::Sqrt(3);
                    ////
                    ////                    LY_E_Pg2 += EtotScint_Pg2 * (graph_gamma_LY->Eval( EtotScint_Pg2, 0, "s" ) - g_unified_positron_nl_pulls[ierr]->Eval( EtotScint_Pg2, 0, "s" ))*GausNL[ierr]/TMath::Sqrt(3)
                    ////                    - Pg2_Scint * (graph_gamma_LY->Eval( Pg2_Scint, 0, "s" ) - g_unified_positron_nl_pulls[ierr]->Eval( Pg2_Scint, 0, "s" ))*GausNL[ierr]/TMath::Sqrt(3);
                    ////
                    ////                    // ------------------------ //
                    ////                    /// neutron - determined from NuWa for 0-0.18 MeV (does not consider leakage)
                    ////                    LY_E_N += Etot_N * ( 0.186 + exp(-1.142-126.4*Etot_N) + exp(-1.217-17.90*Etot_N) ) * (graph_electron_LY->Eval( 0.2, 0, "s" ) - g_unified_positron_nl_pulls[ierr]->Eval( 0.2, 0, "s" ))*GausNL[ierr]/TMath::Sqrt(3);
                    ////
                    //                    /// delayed
                    ////                    LY_E_Ng += EtotScint_Ng * (graph_gamma_LY->Eval( EtotScint_Ng, 0, "s" ) - g_unified_positron_nl_pulls[ierr]->Eval( EtotScint_Ng, 0, "s" ))*GausNL[ierr]/TMath::Sqrt(3) - Ng_Scint * (graph_gamma_LY->Eval( Ng_Scint, 0, "s" ) - g_unified_positron_nl_pulls[ierr]->Eval( Ng_Scint, 0, "s" ))*GausNL[ierr]/TMath::Sqrt(3);
                    //
                    //                    LY_E_Ng += EtotScint_Ng * graph_gamma_LY->Eval( EtotScint_Ng, 0, "s" )*(1+ (ErrorBand->Interpolate(EtotScint_Ng)*GausNL[ierr]/TMath::Sqrt(3))) - Ng_Scint * graph_gamma_LY->Eval( Ng_Scint, 0, "s" )*(1+ (ErrorBand->Interpolate(Ng_Scint)*GausNL[ierr]/TMath::Sqrt(3)));
                    //                }
                    //
                    //                // ------------------------ //
                    //                /// neutron - determined from NuWa for 0-0.18 MeV (does not consider leakage)
                    //                LY_E_N = EtotScint_N * ( 0.186 + exp(-1.142-126.4*EtotScint_N) + exp(-1.217-17.90*EtotScint_N) ) * graph_electron_LY->Eval( 0.2, 0, "s" );
                    //
                    //
                    
                    /// prompt SUM
                    LY_E_P_Sum = LY_E_P + LY_E_Pg1 + LY_E_Pg2 + LY_E_N;
                    //------------------------ //
                    
                    
                    ///////////////////////////// Optical NU
                    /*
                     usr_opt_attenuation_P = 0;
                     usr_pmt_coverage_P  = 0;
                     usr_opt_attenuation_Ng = 0;
                     usr_pmt_coverage_Ng  = 0;
                     */
                    /// prompt
                    usr_r2_P = (X_P*X_P+Y_P*Y_P) * 1e-6;  // mm2 ---> m2
                    usr_z_P  = Z_P * 1e-3;  // mm ---> m
                    
                    global_bin_num = hist_findbin->FindBin(usr_r2_P, usr_z_P);
                    hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
                    
                    if( local_zbin != 0)
                    {
                        continue;
                    }
                    
                    usr_opt_attenuation_P = hist_map_attenuation->GetBinContent(local_xbin, local_ybin);//Relative energy scale
                    usr_pmt_coverage_P  = hist_map_pmt_coverage->GetBinContent(local_xbin, local_ybin);
                    
                    if(RelativeEnergyScaleMatrix)
                    {
                        usr_opt_attenuation_P = usr_opt_attenuation_P*(1+RelativeNominalError*GausRelative);
                    }
                    
                    Opt_E_P_Sum = LY_E_P_Sum * usr_opt_attenuation_P * usr_pmt_coverage_P;  // nH analysis
                    
                    // AdSimple_P = ( (7.84628 * (1 + 3.41294e-02*usr_r2_P) * (1 - 1.21750e-02*usr_z_P - 1.64275e-02*usr_z_P*usr_z_P + 7.33006e-04*pow(usr_z_P,3)))/8.05 );  // AdSimple
                    
                    // Doc7334(old function), updated from http://dayabay.ihep.ac.cn/tracs/dybsvn/browser/dybgaudi/trunk/Reconstruction/QsumEnergy/src/components/QsumEnergyTool.cc
                    
                    //Opt_E_P_Sum = LY_E_P_Sum * AdSimple_P;
                    
                    usr_r2_Ng = (X_Ng*X_Ng+Y_Ng*Y_Ng) * 1e-6;  // mm2 ---> m2
                    usr_z_Ng  = Z_Ng * 1e-3;  // mm ---> m
                    
                    global_bin_num = hist_findbin->FindBin(usr_r2_Ng, usr_z_Ng);
                    hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
                    
                    if( local_zbin != 0)
                    {
                        continue;
                    }
                    
                    usr_opt_attenuation_Ng = hist_map_attenuation->GetBinContent(local_xbin, local_ybin);
                    usr_pmt_coverage_Ng  = hist_map_pmt_coverage->GetBinContent(local_xbin, local_ybin);
                    
                    Opt_E_Ng = LY_E_Ng * usr_opt_attenuation_Ng * usr_pmt_coverage_Ng;  // nH analysis
                    
                    // AdSimple_Ng = ( (7.84628 * (1 + 3.41294e-02*usr_r2_Ng) * (1 - 1.21750e-02*usr_z_Ng - 1.64275e-02*usr_z_Ng*usr_z_Ng + 7.33006e-04*pow(usr_z_Ng,3)))/8.05 );  // AdSimple
                    // Doc7334(old function), updated from http://dayabay.ihep.ac.cn/tracs/dybsvn/browser/dybgaudi/trunk/Reconstruction/QsumEnergy/src/components/QsumEnergyTool.cc
                    
                    //   Opt_E_Ng = LY_E_Ng * AdSimple_Ng;
                    
                    ///////////////////////////// FEE
                    /// prompt
                    //if( Opt_E_P_Sum!=Opt_E_P_Sum ) continue;//Maybe eliminate it, it produces problems. See April 16th Logan's email
                    
                    /// prompt
                    if( Opt_E_P_Sum!=Opt_E_P_Sum || Opt_E_P_Sum<0 )
                    {
                        cout<<" Opt_E_P_Sum = "<< Opt_E_P_Sum <<endl;
                        continue;
                    }
                    
                    FEE_E_P_Sum = 0;
                    
                    FEE_E_P_Sum = Opt_E_P_Sum * graph_electronic->Eval( Opt_E_P_Sum, 0, "s" )*(1+ (ErrorBand->Interpolate(Opt_E_P_Sum)*GausNL[0]/Sigma_Electron));
                    
                    if( FEE_E_P_Sum!=FEE_E_P_Sum || FEE_E_P_Sum<0 ) {
                        cout<<" FEE_E_P_Sum = "<< FEE_E_P_Sum <<endl;
                        continue;
                    }
                    
                    //Random NL: // In each event the NL curve is varied //if GausNL = 0 the energy would be the nominal NL function, otherwise energy = nl_nominal_energy(1+Σerrors)
                    //We don't have pull curves for each term (electron, gamma and electronic), so far use the nominal error in each of the curves.
                    
                    
                    //                for (Int_t ierr = 0; ierr < unified_nl_pars; ierr++)
                    //                {
                    ////                    FEE_E_P_Sum += Opt_E_P_Sum * (graph_electronic->Eval( Opt_E_P_Sum, 0, "s" ) - g_unified_positron_nl_pulls[ierr]->Eval( Opt_E_P_Sum, 0, "s" ))*GausNL[ierr]/TMath::Sqrt(3);
                    //
                    //                    FEE_E_P_Sum += Opt_E_P_Sum * graph_electronic->Eval( Opt_E_P_Sum, 0, "s" )*(1+ (ErrorBand->Interpolate(Opt_E_P_Sum)*GausNL[ierr]/TMath::Sqrt(3)));
                    //
                    //
                    //                }
                    
                    ///////////////////////////// EnergyScale
                    /// prompt
                    Scale_E_P_Sum = FEE_E_P_Sum *EnergyScale;
                    
                    if( Scale_E_P_Sum!=Scale_E_P_Sum || Scale_E_P_Sum<0 )
                    {
                        cout<<" Scale_E_P_Sum = "<< Scale_E_P_Sum <<endl;
                        continue;
                    }
                    ///////////////////////////// FEE
                    /// delayed
                    //if( Opt_E_Ng!=Opt_E_Ng ) continue;//Maybe eliminate it, it produces problems. See April 16th Logan's email
                    
                    
                    /// delayed
                    if( Opt_E_Ng!=Opt_E_Ng || Opt_E_Ng<0 ) {
                        cout<<" Opt_E_Ng = "<< Opt_E_Ng <<endl;
                        continue;
                    }
                    
                    FEE_E_Ng = 0;
                    
                    FEE_E_Ng = Opt_E_Ng * graph_electronic->Eval( Opt_E_Ng, 0, "s" )*(1+ (ErrorBand->Interpolate(Opt_E_Ng)*GausNL[0]/Sigma_Gamma));
                    
                    if( FEE_E_Ng!=FEE_E_Ng || FEE_E_Ng<0 ) {
                        cout<<" FEE_E_Ng = "<< FEE_E_Ng <<", Opt_E_Ng = "<< Opt_E_Ng <<endl;
                        continue;
                    }
                    
                    //                for (Int_t ierr = 0; ierr < unified_nl_pars; ierr++)
                    //                {
                    ////                    FEE_E_Ng += Opt_E_Ng * (graph_electronic->Eval( Opt_E_Ng, 0, "s" ) - g_unified_positron_nl_pulls[ierr]->Eval( Opt_E_Ng, 0, "s" ))*GausNL[ierr]/TMath::Sqrt(3);
                    //
                    //                    FEE_E_Ng += Opt_E_Ng * graph_electronic->Eval( Opt_E_Ng, 0, "s" )*(1+ (ErrorBand->Interpolate(Opt_E_Ng)*GausNL[ierr]/TMath::Sqrt(3)));
                    //
                    //                }
                    ///////////////////////////// EnergyScale
                    /// delayed
                    Scale_E_Ng = FEE_E_Ng *EnergyScale;
                    
                    if( Scale_E_Ng!=Scale_E_Ng || Scale_E_Ng<0 ) {
                        cout<<" Scale_E_Ng = "<< Scale_E_Ng <<endl;
                        continue;
                    }
                    ///////////////////////////// Resolution
                    ///////////////////////////// Resolution
                    
                    //                energy_sigma = 0;
                    //                R2_AA        = 0;
                    //                R2_BB        = 0;
                    //                R_average    = 0;
                    
                    /// prompt
                    global_bin_num = hist_findbin->FindBin(usr_r2_P, usr_z_P);
                    hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
                    
                    if( local_zbin != 0)
                    {
                        continue;
                    }
                    /*
                     R2_AA = (local_xbin-1) * R2_binwidth;
                     R2_BB = local_xbin * R2_binwidth;
                     R_average = sqrt( (R2_AA+R2_BB)/2. );  // volume-weighted average radius of cell
                     */
                    R_average = sqrt( usr_r2_P );
                    
                    roofunc_EnergyResolution->SetParameter( 0, R_average );
                    energy_sigma = roofunc_EnergyResolution->Eval( Scale_E_P_Sum );
                    Res_E_P_Sum = gRandom3->Gaus( Scale_E_P_Sum, energy_sigma );
                    if( Res_E_P_Sum < 0 )  // some low-energy events are smeared to below zero
                    {
                        Res_E_P_Sum = -Res_E_P_Sum;
                    }
                    if( Res_E_P_Sum!=Res_E_P_Sum )
                    {
                        cout<<" Res_E_P_Sum = "<< Res_E_P_Sum <<", R_average="<< R_average <<", Scale_E_P_Sum="<< Scale_E_P_Sum <<", energy_sigma="<< energy_sigma <<endl;
                        //continue;
                    }
                    /// Estimate an "Erec" on which to apply cuts
                    ///////////////////////////// Reconstructed Energy
                    /// prompt
                    Erec_P = Res_E_P_Sum / (usr_opt_attenuation_P *usr_pmt_coverage_P);  // nH analysis
                    //Erec_P = Res_E_P_Sum / AdSimple_P;  // AdSimple
                    
                    /// delayed
                    global_bin_num = hist_findbin->FindBin(usr_r2_Ng, usr_z_Ng);
                    hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
                    
                    if( local_zbin != 0)
                    {
                        continue;
                    }
                    /*
                     R2_AA = (local_xbin-1) * R2_binwidth;
                     R2_BB = local_xbin * R2_binwidth;
                     R_average = sqrt( (R2_AA+R2_BB)/2. );  // volume-weighted average radius of cell
                     */
                    R_average = sqrt( usr_r2_Ng );
                    
                    roofunc_EnergyResolution->SetParameter( 0, R_average );
                    energy_sigma = roofunc_EnergyResolution->Eval( Scale_E_Ng );
                    Res_E_Ng = gRandom3->Gaus( Scale_E_Ng, energy_sigma );
                    if( Res_E_Ng < 0 )  // some low-energy events are smeared to below zero
                        Res_E_Ng = -Res_E_Ng;
                    if( Res_E_Ng!=Res_E_Ng ) {
                        cout<<" Res_E_Ng = "<< Res_E_Ng <<", R_average="<< R_average <<", Scale_E_Ng="<< Scale_E_Ng <<", energy_sigma="<< energy_sigma <<endl;
                        continue;
                    }
                    
                    /// Estimate an "Erec" on which to apply cuts
                    ///////////////////////////// Reconstructed Energy
                    /// delayed
                    Erec_Ng = Res_E_Ng / (usr_opt_attenuation_Ng * usr_pmt_coverage_Ng);  // nH analysis
                    //Erec_Ng = Res_E_Ng / AdSimple_Ng;  // AdSimple
                    
                    ///////////////////////////// Cell
                    ///////////////////////////// Cell
                    
                    ////// Cuts:
                    //////
                    //                if( cap_Target==64 ) continue;  // reject nGd events because not properly considered in toyMC input sample
                    
                    if( Res_E_P_Sum!=Res_E_P_Sum || Res_E_Ng!=Res_E_Ng ) {
                        cout<<" Res_E_P_Sum = "<< Res_E_P_Sum <<", Res_E_Ng = "<< Res_E_Ng << endl;
                        continue;
                    }
                    
                    /// Analysis cuts
                    if( Erec_Ng<EdelayCutLow || Erec_Ng>EdelayCutHigh )
                    {
                        continue;
                    }
                    if( TcapCutHigh>0 && (Tdiff>TcapCutHigh || Tdiff<TcapCutLow) )
                    {
                        continue;
                    }
                    if( Erec_P<EpromptCutLow || Erec_P>EpromptCutHigh )
                    {
                        continue;
                    }
                    if( distCut>0 ) {
                        distPD = sqrt( (X_Ng-X_P)*(X_Ng-X_P) + (Y_Ng-Y_P)*(Y_Ng-Y_P) + (Z_Ng-Z_P)*(Z_Ng-Z_P) );
                        if( distPD>distCut )
                        {
                            continue;
                        }
                    }
                    
                    deltaRoav = (0 + GausOAV*UncertaintyDeltaRoav);
                    deltaZoav = (0 + GausOAV*UncertaintyDeltaZoav);
                    
                    // cut events due to variation in OAV dimensions
                    if( deltaRoav<0 && sqrt(X_P*X_P+Y_P*Y_P)<-deltaRoav )
                    {
                        continue;
                    }
                    if( deltaRoav>0 && sqrt(X_P*X_P+Y_P*Y_P)>Rshield-deltaRoav )
                    {
                        continue;
                    }
                    if( deltaZoav<0 && Z_P<ZreflectorBottom-deltaZoav )
                    {
                        continue;
                    }
                    if( deltaZoav>0 && Z_P>ZreflectorTop-deltaZoav )
                    {
                        continue;
                    }
                    
                    //Efficiency correction:
                    Double_t Eff_p_content = 0;
                    
                    ////// prompt
                    usr_r2_P = pow(deltaRoav + sqrt(X_P*X_P+Y_P*Y_P),2) * 1e-6;  // mm2 ---> m2
                    //if( radialCut>0 && usr_r2>radialCut ) continue;
                    usr_z_P  = (deltaZoav + Z_P) * 1e-3;  // mm ---> m
                    global_bin_num = hist_findbin->FindBin(usr_r2_P, usr_z_P);
                    hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
                    
                    if( local_zbin==0 && local_xbin>=1 && local_xbin<=R2_binnum && local_ybin>=1 && local_ybin<=Z_binnum )
                    {
                        Eff_p_content = h2d_Ep_ratio2center->GetBinContent(local_xbin, local_ybin);
                    }
                    //Efficiency correction:
                    Double_t Eff_d_content = 0;
                    
                    ////// delayed
                    usr_r2_Ng = pow(deltaRoav + sqrt(X_Ng*X_Ng+Y_Ng*Y_Ng),2) * 1e-6;  // mm2 ---> m2
                    //if( radialCut>0 && usr_r2>radialCut ) continue;
                    usr_z_Ng  = (deltaZoav + Z_Ng) * 1e-3;  // mm ---> m
                    global_bin_num = hist_findbin->FindBin(usr_r2_Ng, usr_z_Ng);
                    hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
                    
                    if( local_zbin==0 && local_xbin>=1 && local_xbin<=R2_binnum && local_ybin>=1 && local_ybin<=Z_binnum )
                    {
                        Eff_d_content = h2d_Ed_ratio2center->GetBinContent(local_xbin, local_ybin);
                    }
                    
#ifdef UseVolumes
                    if((local_ybin==1||local_ybin==(Z_binnum)||local_xbin>(R2_binnum-4)))//nH LS
                    {
                        VolumeXbin = 2;
                    }
                    else//GdLS
                    {
                        VolumeXbin = 1;
                    }
#else
                    VolumeXbin = local_xbin;
                    VolumeYbin = local_ybin;
#endif
                    
#ifndef EREC_COMPARISON
                    
                    //Fill histograms: Apply efficiency of the cell to each event, and then add that event in the matrix of the corresponding volume after calculating if the cell is inside the Ls or Gd-Ls
                    
                    Double_t FlatRandomSample = 1;//Normal mode
#ifdef LoadToyTree
                    FlatRandomSample = Double_t((1.*etree_read->GetEntries()/NumberOfEntries));
                    //std::cout << " FLAT RANDOM SAMPLE : " << FlatRandomSample << std::endl;
#endif
                    
                    TruePredictionH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Ev,Eff_p_content*FlatRandomSample);
                    VisiblePredictionH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Res_E_P_Sum,Eff_p_content*FlatRandomSample);
                    
                    DelayedVisiblePredictionH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Res_E_Ng,Eff_d_content*FlatRandomSample);
                    HighResoMatrixH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Ev,Res_E_P_Sum,Eff_p_content*FlatRandomSample);//Fine grid
                    MatrixH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Ev,Res_E_P_Sum,Eff_p_content*FlatRandomSample);//neutrino energy vs visible energy
                    TransMatrixH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Res_E_P_Sum,Ev,Eff_p_content*FlatRandomSample);
                    HighResoTransMatrixH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Res_E_P_Sum,Ev,Eff_p_content*FlatRandomSample);//Fine grid
#else
                    //Fill with Erec for comparison:
                    TruePredictionH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Ev,Eff_p_content*FlatRandomSample);
                    VisiblePredictionH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Erec_P,Eff_p_content*FlatRandomSample);
                    
                    DelayedVisiblePredictionH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Erec_Ng,Eff_d_content*FlatRandomSample);
                    HighResoMatrixH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Ev,Erec_P,Eff_p_content*FlatRandomSample);//Fine grid
                    MatrixH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Ev,Erec_P,Eff_p_content*FlatRandomSample);//neutrino energy vs visible energy
                    TransMatrixH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Erec_P,Ev,Eff_p_content*FlatRandomSample);
                    HighResoTransMatrixH[Systematic][AD][VolumeXbin-1][VolumeYbin-1]->Fill(Erec_P,Ev,Eff_p_content*FlatRandomSample);//Fine grid
#endif
                    
                    //Save number of events in each cell to properly scale the cells afterwards:
                    
                    ADIntegralEvents[Systematic][AD] = 0;//initialize
                    
                    for(Int_t idx=0; idx<VolumeX; idx++)
                    {
                        for(Int_t idy=0; idy<VolumeY; idy++)
                        {
                            IntegralEvents[Systematic][AD][idx][idy] = VisiblePredictionH[Systematic][AD][idx][idy]->Integral();
                            ADIntegralEvents[Systematic][AD] = ADIntegralEvents[Systematic][AD] + IntegralEvents[Systematic][AD][idx][idy];
                        }
                    }
#ifdef SaveTree
                    if(Systematic==0)
                    {
                        toy->Fill();//Save only nominal tree for 1 AD
                    }
#endif
                }//Events
            }//True energy Loop
        }//AD
        
#ifdef SaveTree
        if(Systematic==0)//Uncomment to save other systematic trees
        {
            roofile_toy->cd();
            
            toy->Write();//Save only Nominal Tree
            
            roofile_toy->Close();
        }
#endif
        
        //Save event ratios:
        
        IAVMatrix = 0;
        OAVMatrix = 0;
        NLMatrix = 0;
        ResolutionMatrix = 0;
        RelativeEnergyScaleMatrix = 0;
        //EfficiencyMatrix = 0;
        AllMatrix = 0;
        
        
#ifndef EREC_COMPARISON
        
        //Save txt file in case we want to use it externally:
        FILE *f = fopen(("./Inputs/HInputs/"+SystematicS+"ToyMCEventRatio.txt").c_str(), "w");
#else
        FILE *f = fopen(("./Inputs/HInputs/E_REC_"+SystematicS+"ToyMCEventRatio.txt").c_str(), "w");
        
#endif
        if (f == NULL)
        {
            printf("Error opening file!\n");
            exit(1);
        }
        
        
        for(Int_t AD = 0; AD<NADs; AD++)
        {
            fprintf(f, "#ToyEvents in each AD\n");
            
            fprintf(f, "%d %f\n", AD, ADIntegralEvents[Systematic][AD]);
            
            fprintf(f, "#Events in each cell for AD %d\n", AD);
            
            for(int idx=0; idx<VolumeX; idx++)
            {
                for(int idy=0; idy<VolumeY; idy++)
                {
                    fprintf(f, "%d %d %d %f\n", AD, idx, idy, IntegralEvents[Systematic][AD][idx][idy] );
                }
            }
        }
        
        fclose(f);
        
        //Draw it:
#ifndef UseVolume
        TH2D* MapEvents_ratio2center = new TH2D("MapEventsRatio2center","MapEventsRatio2center",VolumeX,VolumeX_lower,VolumeX_upper,VolumeY,VolumeY_lower,VolumeY_upper);
        
        for(Int_t AD = 0; AD<NADs; AD++)
        {
            for(Int_t idx = 0; idx<VolumeX; idx++)
            {
                for(Int_t idy = 0; idy<VolumeY; idy++)
                {
                    MapEvents_ratio2center->SetBinContent(idx+1,idy+1,IntegralEvents[Systematic][AD][idx][idy]/IntegralEvents[Systematic][AD][0][Int_t(VolumeY/2)]);
                }
            }
        }
        
        TCanvas* MapEventC = new TCanvas("MapEventRatio2Center","MapEventRatio2Center");
        
        MapEventC->cd(1);
        
        MapEvents_ratio2center->Draw("colz");
        
        MapEventC->Print(("./Images/Hydrogen/Detector/"+SystematicS+"MapEventRatio2Center.eps").c_str());
        
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
                    PercentualEvents[Systematic][AD][idx][idy] = IntegralEvents[Systematic][AD][idx][idy]/ADIntegralEvents[Systematic][AD];
                    
                    MapEvents_ratio2total->SetBinContent(idx+1,idy+1,PercentualEvents[Systematic][AD][idx][idy]);
                }
            }
        }
        
        TCanvas* MapEventTotalC = new TCanvas("MapEventRatio2Total","MapEventRatio2Total");
        
        MapEventTotalC->cd(1);
        
        MapEvents_ratio2total->Draw("colz");
        
        MapEventTotalC->Print(("./Images/Hydrogen/Detector/"+SystematicS+"MapEventRatio2Total.eps").c_str());
        
        delete MapEventTotalC;
        
        TFile* SaveFile = new TFile("./Inputs/HInputs/MapEventRatio2Total.root","recreate");
        {
            MapEvents_ratio2total->Write();
        }
        
        delete SaveFile;
        
        delete MapEvents_ratio2total;
        
        //        for(int idx=0; idx<MaxCellNum; idx++)
        //        {
        //            delete hEp_cc[idx];
        //
        //            delete hEd_cc[idx];
        //
        //            delete hEp_cc_clone[idx];
        //
        //            delete hEd_cc_clone[idx];
        //        }
        //        for(int idx=0; idx<MaxColumnNum; idx++)
        //        {
        //            delete hEp_cl[idx];
        //
        //            delete hEd_cl[idx];
        //        }
        
        
#ifdef PrintEps
        TCanvas* NNFineMatrixC = new TCanvas("NNF","NNF");
        NNFineMatrixC->Divide(NADs/2,2);
        
        TCanvas* NNMatrixC = new TCanvas("NNM","NNM");
        NNMatrixC->Divide(NADs/2,2);
        
        TCanvas* NNTransC = new TCanvas("NNT","NNT");
        NNTransC->Divide(NADs/2,2);
        
        TCanvas* NNFineTransC = new TCanvas("NNFT","NNFT");
        NNFineTransC->Divide(NADs/2,2);
        
        for(Int_t i = 0; i<NADs;i++)
        {
            NNMatrixC->cd(i+1);
            MatrixH[Systematic][i][0][0]->Draw("colz");
            
            NNFineMatrixC->cd(i+1);
            HighResoMatrixH[Systematic][i][0][0]->Draw("colz");
            
            NNTransC->cd(i+1);
            TransMatrixH[Systematic][i][0][0]->Draw("colz");
            
            NNFineTransC->cd(i+1);
            HighResoTransMatrixH[Systematic][i][0][0]->Draw("colz");
        }
        
        NNFineMatrixC->Print(("./Images/Hydrogen/ResponseMatrices/"+SystematicS+"NoNormFineHydrogenResponseMatrix.eps").c_str());
        NNMatrixC->Print(("./Images/Hydrogen/ResponseMatrices/"+SystematicS+"NoNormHydrogenResponseMatrix.eps").c_str());
        NNTransC->Print(("./Images/Hydrogen/ResponseMatrices/"+SystematicS+"NoNormTransposeHydrogenMatrix.eps").c_str());
        NNFineTransC->Print(("./Images/Hydrogen/ResponseMatrices/"+SystematicS+"NoNormFineTransposeHydrogenMatrix.eps").c_str());
        
        delete NNMatrixC;
        delete NNFineMatrixC;
        delete NNTransC;
        delete NNFineTransC;
#endif
        
        Double_t Norma[n_etrue_bins];//true->vis
        Double_t HighResoNorma[MatrixBins];//true->vis
        Double_t NormaTrans[n_evis_bins];//vis->true
        Double_t HighResoNormaTrans[MatrixBins];//vis->true
        //Normalize the Matrix
        for(int idx=0; idx<VolumeX; idx++)
        {
            for(int idy=0; idy<VolumeY; idy++)
            {
                for(Int_t AD = 0; AD<NADs;AD++)
                {
                    for(Int_t i=0;i<n_etrue_bins;i++)
                    {
                        Norma[i]=0;
                        for(Int_t j=0;j<n_evis_bins;j++)
                        {
                            Norma[i] = Norma[i]+MatrixH[Systematic][AD][idx][idy]->GetBinContent(i+1,j+1);// true->vis
                        }
                    }
                    
                    for(Int_t i=0;i<n_evis_bins;i++)
                    {
                        NormaTrans[i]=0;
                        
                        for(Int_t j=0;j<n_etrue_bins;j++)
                        {
                            NormaTrans[i] = NormaTrans[i]+TransMatrixH[Systematic][AD][idx][idy]->GetBinContent(i+1,j+1);// vis->true
                        }
                    }
                    
                    for (Int_t i = 0; i < n_etrue_bins; i++)
                    {
                        for (Int_t j = 0; j < n_evis_bins; j++)
                        {
                            if(Norma[i]!=0)
                            {
                                MatrixH[Systematic][AD][idx][idy]->SetBinContent(i+1,j+1,MatrixH[Systematic][AD][idx][idy]->GetBinContent(i+1,j+1)/Norma[i]);//true->vis
                            }
                        }
                    }
                    
                    for (Int_t i = 0; i < n_evis_bins; i++)
                    {
                        for (Int_t j = 0; j < n_etrue_bins; j++)
                        {
                            if(NormaTrans[i]!=0)
                            {
                                TransMatrixH[Systematic][AD][idx][idy]->SetBinContent(i+1,j+1,TransMatrixH[Systematic][AD][idx][idy]->GetBinContent(i+1,j+1)/NormaTrans[i]);
                            }
                            
                        }
                    }
                    
                    for(Int_t i=0;i<MatrixBins;i++)
                    {
                        HighResoNorma[i]=0;
                        HighResoNormaTrans[i]=0;
                        
                        for(Int_t j=0;j<MatrixBins;j++)
                        {
                            HighResoNorma[i] = HighResoNorma[i]+HighResoMatrixH[Systematic][AD][idx][idy]->GetBinContent(i+1,j+1);// true->vis
                            HighResoNormaTrans[i] = HighResoNormaTrans[i]+HighResoTransMatrixH[Systematic][AD][idx][idy]->GetBinContent(i+1,j+1);// vis->true
                            
                        }
                    }
                    
                    for(Int_t i=0;i<MatrixBins;i++)
                    {
                        for(Int_t j=0;j<MatrixBins;j++)
                        {
                            if(HighResoNorma[i]!=0)
                            {
                                HighResoMatrixH[Systematic][AD][idx][idy]->SetBinContent(i+1,j+1,HighResoMatrixH[Systematic][AD][idx][idy]->GetBinContent(i+1,j+1)/HighResoNorma[i]);//true->vis
                            }
                            
                            if(HighResoNormaTrans[i]!=0)
                            {
                                HighResoTransMatrixH[Systematic][AD][idx][idy]->SetBinContent(i+1,j+1,HighResoTransMatrixH[Systematic][AD][idx][idy]->GetBinContent(i+1,j+1)/HighResoNormaTrans[i]);
                            }
                        }
                    }
                }
#ifdef PrintEps
                if(idx==0&&idy==0)
                {
                    TCanvas* FineMatrixC = new TCanvas("F","F");
                    FineMatrixC->Divide(NADs/2,2);
                    
                    TCanvas* MatrixC = new TCanvas("","");
                    MatrixC->Divide(NADs/2,2);
                    
                    TCanvas* PredictionC = new TCanvas("P","P");
                    PredictionC->Divide(NADs/2,2);
                    
                    TCanvas* VisPredictionC = new TCanvas("VP","VP");
                    VisPredictionC->Divide(NADs/2,2);
                    
                    TCanvas* DelayedVisPredictionC = new TCanvas("DVP","DVP");
                    DelayedVisPredictionC->Divide(NADs/2,2);
                    
                    TCanvas* TransC = new TCanvas("T","T");
                    TransC->Divide(NADs/2,2);
                    
                    TCanvas* FineTransC = new TCanvas("FT","FT");
                    FineTransC->Divide(NADs/2,2);
                    
                    for(Int_t i = 0; i<NADs;i++)
                    {
                        PredictionC->cd(i+1);
                        TruePredictionH[Systematic][i][idx][idy]->Draw();
                        
                        VisPredictionC->cd(i+1);
                        VisiblePredictionH[Systematic][i][idx][idy]->Draw();
                        
                        DelayedVisPredictionC->cd(i+1);
                        DelayedVisiblePredictionH[Systematic][i][idx][idy]->Draw();
                        
                        MatrixC->cd(i+1);
                        MatrixH[Systematic][i][idx][idy]->Draw("colz");
                        
                        FineMatrixC->cd(i+1);
                        HighResoMatrixH[Systematic][i][idx][idy]->Draw("colz");
                        
                        TransC->cd(i+1);
                        TransMatrixH[Systematic][i][idx][idy]->Draw("colz");
                        
                        FineTransC->cd(i+1);
                        HighResoTransMatrixH[Systematic][i][idx][idy]->Draw("colz");
                    }
                    
                    PredictionC->Print(("./Images/Hydrogen/ToyMCOutputs/"+SystematicS+"ToyMCTrueHydrogenPrediction.eps").c_str());
                    VisPredictionC->Print(("./Images/Hydrogen/ToyMCOutputs/"+SystematicS+"ToyMCVisibleHydrogenPrediction.eps").c_str());
                    DelayedVisPredictionC->Print(("./Images/Hydrogen/ToyMCOutputs/"+SystematicS+"ToyMCDelayedVisibleHydrogenPrediction.eps").c_str());
                    FineMatrixC->Print(("./Images/Hydrogen/ResponseMatrices/"+SystematicS+"FineHydrogenResponseMatrix.eps").c_str());
                    MatrixC->Print(("./Images/Hydrogen/ResponseMatrices/"+SystematicS+"HydrogenResponseMatrix.eps").c_str());
                    TransC->Print(("./Images/Hydrogen/ResponseMatrices/"+SystematicS+"TransposeHydrogenMatrix.eps").c_str());
                    FineTransC->Print(("./Images/Hydrogen/ResponseMatrices/"+SystematicS+"FineTransposeHydrogenMatrix.eps").c_str());
                    
                    delete PredictionC;
                    delete VisPredictionC;
                    delete DelayedVisPredictionC;
                    delete MatrixC;
                    delete FineMatrixC;
                    delete TransC;
                    delete FineTransC;
                }
#endif
            }
        }
        
        
        //
        //        //Load event ratios:
        //
        //        std::string line;
        //
        //        ifstream mainfile(("./Inputs/HInputs/"+SystematicS+"ToyMCEventRatio.txt").c_str());
        //
        //        Int_t linenum=0;//<---caution: only increments for lines that do not begin with #
        //
        
        //
        //        while(!mainfile.eof())
        //        {
        //            std::getline(mainfile,line);
        //            std::string firstchar = line.substr(0,1);
        //
        //            if(firstchar=="#") continue;//<-- ignore lines with comments
        //
        //            std::istringstream iss(line);
        //
        //            if(linenum == 0)
        //            {
        //                Int_t AD = 0;
        //
        //                iss >> AD >> ADIntegral[Systematic][AD];
        //
        //                std::cout << "line " << linenum << " the integral in AD " << AD << " is " << ADIntegral[Systematic][AD] << std::endl;
        //
        //                linenum++;
        //            }
        //            else if(linenum <=VolumeX*VolumeY)
        //            {
        //                Int_t AD = 0;
        //                Int_t idx = 0;
        //                Int_t idy = 0;
        //
        //                iss >> AD >> idx >> idy >> CellIntegral[Systematic][AD][idx][idy];
        //                //                if(column==0) AD=atoi(firstchar.c_str());
        //                //
        //                //                if(column==1) idx=atoi(sub.c_str());
        //                //
        //                //                if(column==2) idy=atoi(sub.c_str());
        //                //
        //                //                if(column==3)
        //
        //                std::cout << "line " << linenum << " the integral in AD " << AD << " cell " << idx << " , " << idy << " is " << CellIntegral[Systematic][AD][idx][idy] << std::endl;
        //
        //                linenum++;
        //
        //                if(linenum>VolumeX*VolumeY)
        //                {
        //                    linenum = 0;//reset counter
        //                }
        //            }
        //        }
        
#ifndef EREC_COMPARISON
        
        TFile* nHMatrixF = new TFile(("./ResponseMatrices/Hydrogen/"+SystematicS+Form("ResponseMatrix%i_%i.root",n_evis_bins,n_etrue_bins)).c_str(),"recreate");
#else
        TFile* nHMatrixF = new TFile(("./ResponseMatrices/Hydrogen/E_REC_"+SystematicS+Form("ResponseMatrix%i_%i.root",n_evis_bins,n_etrue_bins)).c_str(),"recreate");
#endif
        //TFile* SaveMatrix = new TFile(Form("./ResponseMatrices/Hydrogen/RandomResponseMatrix%i_%i.root",n_evis_bins,n_etrue_bins),"recreate");
        
        for(int idx=0; idx<VolumeX; idx++)
        {
            for(int idy=0; idy<VolumeY; idy++)
            {
                for(Int_t AD = 0; AD<NADs; AD++)//Save different AD matrices to produce the covariance matrices.
                {
                    HighResoMatrixH[Systematic][AD][idx][idy]->Write(Form("FineEvisEnu%i,Cell%i,%i",AD+1,idx,idy));//Save the matrices.
                    MatrixH[Systematic][AD][idx][idy]->Write(Form("EvisEnu%i,Cell%i,%i",AD+1,idx,idy));//Save the matrices.
                    TransMatrixH[Systematic][AD][idx][idy]->Write(Form("EnuEvis%i,Cell%i,%i",AD+1,idx,idy));
                }
                
            }
        }
        delete nHMatrixF;
        
    }//Systematics
    
    
    delete h2d_Ep_ratio2center;
    delete h2d_Ed_ratio2center;
    
    delete gRandom3;
    delete RandomSysUncorr;
    delete RandomSysCorr;
    
}

// --- ---- --- ---- --- ---- --- ---- --- ---- --- Function definition
// --- ---- --- ---- --- ---- --- ---- --- ---- ---
// --- ---- --- ---- --- ---- --- ---- --- ---- ---

void nHToyMC :: func_initialization()
{
    //////
    hist_findbin = new TH2D("hist_findbin","", R2_binnum,R2_lower,R2_upper, Z_binnum,Z_lower,Z_upper);
    
    if( MaxCellNum<(R2_binnum+2)*(Z_binnum+2) ) {
        printf("\n\n MaxCellNum = %d, but must be >= %d\n\n", MaxCellNum, (R2_binnum+2)*(Z_binnum+2) );
        return;
    }
    
    ////// ~/WORK/jixp/hapy/usr_job/work_mc_Ep/Production/process/plot_obj.cc
#ifdef ReactorShapeinToy
    TFile *file_h_Ev = new TFile("./Inputs/HInputs/Data/file_h_Ev.root", "read");

    h_Ev_normal = (TH1D*)file_h_Ev->Get("h_Ev");
#endif
    ////// NL
    //TFile* file_IHEP_NL_models = new TFile("./Inputs/HInputs/Data/IHEP_NL_models/Model1.root", "read");
    TFile* file_IHEP_NL_models = new TFile("./Inputs/unified_nl_data/energyModel_march2015.root", "read");

    graph_electron_LY = (TGraph*)file_IHEP_NL_models->Get("electronScintNL");// positron = electron + gamma1 + gamma2
    graph_gamma_LY    = (TGraph*)file_IHEP_NL_models->Get("gammaScintNL");
    graph_electronic  = (TGraph*)file_IHEP_NL_models->Get("electronicsNL");
    
    g_unified_positron_nl = (TGraph*)file_IHEP_NL_models->Get("nominal")->Clone();
    
    //Extract sigma band:
    const Int_t NPoints = g_unified_positron_nl->GetN();
    Double_t X_axis[NPoints];
    Double_t Y_axis[NPoints];

    for(Int_t point = 0; point < NPoints; point++)
    {
        g_unified_positron_nl->GetPoint(point,X_axis[point],Y_axis[point]);
    }
    
    ErrorBand = new TH1D("NLErrorBand", "NLErrorBand", NPoints, X_axis[0],X_axis[NPoints-1]);//Save in TH1D to interpolate later

    for(Int_t point = 0; point < NPoints; point++)
    {
        ErrorBand->SetBinContent(point+1, g_unified_positron_nl->GetErrorY(point));//No difference between GetErrorYhigh, GetErrorY and GetErrorYLow (symmetric)
    }
    
    //Print:
    
#ifdef PrintEps
    TCanvas* ErrorC = new TCanvas("ErrorC","ErrorC");

    ErrorBand->Draw();
    
    ErrorC->Print("./Images/Hydrogen/Detector/NL_Error_band.eps");
    
    delete ErrorC;
#endif
    
    for (Int_t i = 0; i < unified_nl_pars; i++)
    {
        g_unified_positron_nl_pulls[i] =  (TGraph*)file_IHEP_NL_models->Get(Form("pull%d",i))->Clone();
    }
    
    //NL-Gadolinium models:
    
   // TFile* unifiedF_NL_models = new TFile("./Inputs/unified_nl_data/nl_models_final.root", "read");
//    g_unified_positron_nl = (TGraph*)unifiedF_NL_models->Get("positron_0")->Clone();
//    g_unified_positron_nl_pulls[0] =  (TGraph*)unifiedF_NL_models->Get(Form("positron_%d",1))->Clone();
//    g_unified_positron_nl_pulls[1] =  (TGraph*)unifiedF_NL_models->Get(Form("positron_%d",2))->Clone();
//    g_unified_positron_nl_pulls[2] =  (TGraph*)unifiedF_NL_models->Get(Form("positron_%d",3))->Clone();
//    g_unified_positron_nl_pulls[3] =  (TGraph*)unifiedF_NL_models->Get(Form("positron_%d",4))->Clone();
    
//    delete unifiedF_NL_models;
    
    ////// NU
    TFile *file_hist_attenuation_rel = new TFile("./Inputs/HInputs/Data/uniformity/hist_attenuation_rel.root", "read");
    hist_map_attenuation = (TH2D*)file_hist_attenuation_rel->Get("hist_attenuation_rel");
    hist_map_attenuation->SetName("hist_map_attenuation");
    delete file_hist_attenuation_rel;
    
    TFile *file_hist_coverage_rel = new TFile("./Inputs/HInputs/Data/uniformity/hist_coverage_rel.root", "read");
    hist_map_pmt_coverage = (TH2D*)file_hist_coverage_rel->Get("hist_coverage_rel");
    hist_map_pmt_coverage->SetName("hist_map_pmt_coverage");
    delete file_hist_coverage_rel;
    ///
    gRandom3 = new TRandom3();
    RandomSysUncorr = new TRandom3();
    RandomSysCorr = new TRandom3();
    
    /// energy resolution
    roofunc_EnergyResolution = new TF1("roofunc_EnergyResolution", this, &nHToyMC::func_EnergyResolution, 0,20, 1,"nHToyMC","func_EnergyResolution");
}

//TH2D* nHToyMC :: LoadnHMatrix(Int_t AD,Int_t idx, Int_t idy)
//{
//    return HighResoMatrixH[AD][idx][idy];
//}
//
//Double_t nHToyMC :: GetEventsByCell(Int_t AD,Int_t idx, Int_t idy)
//{
//    return PercentualEvents[AD][idx][idy];
//}

//Double_t nHToyMC :: GetMarginanNL(TGraph* Full_NL, TGraph* NLscint, TGraph* NLcerenkov, TGraph* NLelectronics, Double_t Energy, Int_t mode)
//{
//    Double_t result;
//    
//    //Full NL = a(f_q + k_c*f_c)*f_e
//    
//    if(mode==0)//Extract scintillator
//    {
//        result = Full_NL->Eval( Energy, 0, "s" )/(a*NLelectronics->Eval( Energy, 0, "s" )) - k_c*NLcerenkov->Eval( Energy, 0, "s");
//    }
//    else if(mode==1)//Extract cerenkov
//    {
//        result =  (Full_NL->Eval( Energy, 0, "s" )/(a*NLelectronics->Eval( Energy, 0, "s" )) - NLscint->Eval( Energy, 0, "s" ))/k_c;
//    }
//    else//Extract electronics
//    {
//        result =  Full_NL->Eval( Energy, 0, "s" )/a*(NLscint->Eval( Energy, 0, "s" )+k_c*NLcerenkov->Eval( Energy, 0, "s" ));
//    }
//    
//    return result;
//}