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

/// select whether nGd or nH IBD analysis - nGd events are currently rejected in the code.  Time cuts are not currently applied.
// 0==special cuts, 1==hydrogen, 64==gadolinium
#define SpecialCuts
//#define Hydrogen
//#define Gadolinium

#ifdef SpecialCuts  /// Time cuts are not currently applied.
    const double EpromptCutLow=0.0, EpromptCutHigh=12.0, EdelayCutLow=0.0, EdelayCutHigh=3.3;// 1.5, 12.0, 1.5, 3.3
    const double TcapCutLow=500.0, TcapCutHigh=400000.0, distCut=500.0;
    const double radialCut=0.0;//3.9;  // adjust OAV radius [m^2]: nominal=4.0
#elif def Hydrogen
    const double EpromptCutLow=1.5, EpromptCutHigh=12.0, EdelayCutLow=1.7955, EdelayCutHigh=2.6535;
    const double TcapCutLow=500.0, TcapCutHigh=400000.0, distCut=500.0;
    const double radialCut=0.0;  // adjust OAV radius [m^2]: nominal=4.0
#elif def Gadolinium
    const double EpromptCutLow=1.5, EpromptCutHigh=12.0, EdelayCutLow=6.0, EdelayCutHigh=12.0;
    const double TcapCutLow=500.0, TcapCutHigh=200000.0, distCut=0.0;
    const double radialCut=0.0;  // adjust OAV radius [m^2]: nominal=4.0
#endif

// --- ---- --- ---- --- ---- --- ---- --- ---- --- Global
#define ReDoNominal // if you change binning use this to remake the nominal matrix
#define FullEnergyResolution
#define UseDelayInformation
//#define SaveTree // To save toy tree
//#define LoadTree // To load eprompt tree
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
    
    Double_t GausRelative,GausIAV,GausNL,GausReso,GausResoCorr,GausEff;
    Double_t ResolutionError,ResolutionErrorUncorrelated;
    
    //Binning parameters:
    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.
    
    Int_t n_evis_bins;
    Int_t n_etrue_bins;
    
    Double_t IAVNominalError; // relative uncertainty of the IAV thickness for each AD.
    Double_t RelativeNominalError; // relative uncertainty of attenuation length for each AD.
    //Systematic parameters;
    bool RelativeEnergyScaleMatrix;
    bool IAVMatrix;
    bool NLMatrix;
    bool ResolutionMatrix;
    bool EfficiencyMatrix;
    
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
    
    TGraph *g_unified_positron_nl;
    TGraph *g_unified_positron_nl_pulls[4];
    
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
    TH1D* TruePredictionH[MaxDetectors][VolumeX][VolumeY];
    TH1D* VisiblePredictionH[MaxDetectors][VolumeX][VolumeY];
    TH1D* DelayedVisiblePredictionH[MaxDetectors][VolumeX][VolumeY];
    TH2D* TransMatrixH[MaxDetectors][VolumeX][VolumeY];
    TH2D* HighResoTransMatrixH[MaxDetectors][VolumeX][VolumeY];
    TH2D* MatrixH[MaxDetectors][VolumeX][VolumeY];
public:
    nHToyMC(NominalData*);//constructor
    ~nHToyMC();//destructor
    void Toy(bool);//main
    TH2D* LoadnHMatrix(Int_t,Int_t,Int_t);
    void GeneratenHResponseMatrix();
    TH2D* HighResoMatrixH[MaxDetectors][VolumeX][VolumeY];
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
    RelativeEnergyScaleMatrix = Data->GetRelativeEnergyScaleMatrix();
    //    RelativeEnergyOffsetMatrix = Data->GetRelativeEnergyOffsetMatrix();
    //    AbsoluteEnergyScaleMatrix = Data->GetAbsoluteEnergyScaleMatrix();
    //    AbsoluteEnergyOffsetMatrix = Data->GetAbsoluteEnergyOffsetMatrix();
    IAVMatrix = Data->GetIAVMatrix();
    NLMatrix = Data->GetNLMatrix();
    ResolutionMatrix = Data->GetResolutionMatrix();
    EfficiencyMatrix = Data->GetEfficiencyMatrix();
    
    //Errors in parameters:
    ResolutionError = Data->GetResolutionError();
    ResolutionErrorUncorrelated = Data->GetResoUncorrelatedError();
    
    NADs = Data->GetADs();
    
    for(Int_t AD = 0; AD < NADs; AD++)
    {
        for(int idx=0; idx<VolumeX; idx++)
        {
            for(int idy=0; idy<VolumeY; idy++)
            {
                HighResoMatrixH[AD][idx][idy] = new TH2D(Form("Fine_nHResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),Form("Fine_nHResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),MatrixBins,0,FinalVisibleEnergy,MatrixBins,0,FinalVisibleEnergy);//from true to visible
                
                MatrixH[AD][idx][idy] = new TH2D(Form("nHResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),Form("nHResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),n_etrue_bins,enu_bins,n_evis_bins,evis_bins);//from true to visible
                
                TransMatrixH[AD][idx][idy] = new TH2D(Form("nHTransResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),Form("nHTransResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),n_evis_bins,evis_bins,n_etrue_bins,enu_bins);//from visible to true
                
                HighResoTransMatrixH[AD][idx][idy] = new TH2D(Form("Fine_nHTransResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),Form("Fine_nHTransResponseMatrix_AD%i, Cell%i,%i",AD+1,idx,idy),MatrixBins,0,FinalVisibleEnergy,MatrixBins,0,FinalVisibleEnergy);//from visible to true
                
                TruePredictionH[AD][idx][idy] = new TH1D(Form("Pred_AD%i, Cell%i,%i",AD+1,idx,idy),Form("Pred_AD%i, Cell%i,%i",AD+1,idx,idy), Nbins, InitialEnergy, FinalEnergy);
                
                VisiblePredictionH[AD][idx][idy] = new TH1D(Form("VisPred_AD%i, Cell%i,%i",AD+1,idx,idy),Form("VisPred_AD%i, Cell%i,%i",AD+1,idx,idy), MatrixBins, InitialVisibleEnergy, FinalVisibleEnergy);
                
                DelayedVisiblePredictionH[AD][idx][idy] = new TH1D(Form("DelayedVisPred_AD%i, Cell%i,%i",AD+1,idx,idy),Form("DelayedVisPred_AD%i, Cell%i,%i",AD+1,idx,idy), MatrixBins, InitialVisibleEnergy, FinalVisibleEnergy);
            }
        }
        
        
    }
    
    
}

nHToyMC :: ~nHToyMC()
{
    for(Int_t AD = 0; AD < NADs; AD++)
    {
        for(int idx=0; idx<VolumeX; idx++)
        {
            for(int idy=0; idy<VolumeY; idy++)
            {
                delete TruePredictionH[AD][idx][idy];
                delete HighResoMatrixH[AD][idx][idy];
                delete MatrixH[AD][idx][idy];
                delete TransMatrixH[AD][idx][idy];
                delete VisiblePredictionH[AD][idx][idy];
                delete DelayedVisiblePredictionH[AD][idx][idy];
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
    Double_t res = (7.5/sqrt(E) + 0.9 + 0.865*(R-0.98))*0.01+ResolutionBias;
    res = E*res;
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

void nHToyMC :: Toy(bool mode)
{

    
#ifndef ReDoNominal
    if(mode!=0)
    {
#endif
        
        func_initialization();
        
        Double_t cap_Ev;             // neutrino energy
        Double_t cap_Target;         // neutron captrure target. e.g. 1 for H
        Double_t X_Ng, Y_Ng, Z_Ng;   // vertex of gamma(s) induced by neutron capture
        Double_t Edep_Ng;            // deposited energy of gamma(s) induced by neutron capture
        Double_t EtotScint_Ng;       // the initial energy of gamma(s) going into scintillator(LS or GdLS)
        Double_t Etot_Ng;            // the total energy of gamma(s) when generated
        Double_t EdepE_Ng;           // the deposited energy of secondary electrons induced by the gamma(s)
        Double_t EdepEA_Ng;
        Double_t EtotE_Ng;           // the total energy of electrons induced by the gamma(s)
        Double_t Etot_N;             // the kinetic energy of neutron generated from the IBD
        
        Double_t X_P, Y_P, Z_P;  // vetex of positron
        Double_t Etot_P;         // the initial energy of positron when generated
        Double_t Etot_Pg1;       // the initial energy of one gamma generated from IBD
        Double_t Etot_Pg2;       // the initial energy of another gamma generated from IBD
        Double_t EtotScint_P;    // the initial energy of positron going into scintillator(LS or GdLS)
        Double_t EtotScint_Pg;   // the initial energy of gamma(s) induced by the positron before annihilation
        Double_t EtotScint_Pg1;  // the initial energy of one gamma going into scintillator(LS or GdLS)
        Double_t EtotScint_Pg2;  // ... the meaning of variables bellow can be refered to the neutron capture information
        Double_t Edep_P;
        Double_t EdepAcrylic_P;
        Double_t Edep_Pg1;
        Double_t Edep_Pg2;
        Double_t EdepE_P;
        Double_t EdepEA_P;
        Double_t EtotE_P;
        Double_t EdepE_Pg1;
        Double_t EdepEA_Pg1;
        Double_t EtotE_Pg1;
        Double_t EdepE_Pg2;
        Double_t EdepEA_Pg2;
        Double_t EtotE_Pg2;
        
        Int_t seed_generator = 1;
        Int_t seed_generator_corr = 1;
        Int_t seed_generator_uncorr = 1;
#ifdef LoadTree
        
        ////////////
        
        TFile *roofile_input = new TFile("./Inputs/HInputs/Data/Eprompt.root", "read");
        TTree *wtree = (TTree*)roofile_input->Get("wtree");
        long entries_wtree = wtree->GetEntries();
        cout<<" ---> input entries: "<<entries_wtree<<endl;
        
        wtree->SetBranchAddress("cap_Target",   &cap_Target);
        wtree->SetBranchAddress("X_Ng",         &X_Ng);
        wtree->SetBranchAddress("Y_Ng",         &Y_Ng);
        wtree->SetBranchAddress("Z_Ng",         &Z_Ng);
        wtree->SetBranchAddress("Edep_Ng",      &Edep_Ng);
        wtree->SetBranchAddress("EtotScint_Ng", &EtotScint_Ng);
        wtree->SetBranchAddress("Etot_Ng",      &Etot_Ng);
        wtree->SetBranchAddress("EdepE_Ng",     &EdepE_Ng);
        wtree->SetBranchAddress("EdepEA_Ng",    &EdepEA_Ng);
        wtree->SetBranchAddress("EtotE_Ng",     &EtotE_Ng);
        wtree->SetBranchAddress("Etot_N",       &Etot_N);
        wtree->SetBranchAddress("X_P",          &X_P);
        wtree->SetBranchAddress("Y_P",          &Y_P);
        wtree->SetBranchAddress("Z_P",          &Z_P);
        wtree->SetBranchAddress("Etot_P",       &Etot_P);
        wtree->SetBranchAddress("Etot_Pg1",     &Etot_Pg1);
        wtree->SetBranchAddress("Etot_Pg2",     &Etot_Pg2);
        wtree->SetBranchAddress("EtotScint_P",  &EtotScint_P);
        wtree->SetBranchAddress("EtotScint_Pg", &EtotScint_Pg);
        wtree->SetBranchAddress("EtotScint_Pg1", &EtotScint_Pg1);
        wtree->SetBranchAddress("EtotScint_Pg2", &EtotScint_Pg2);
        wtree->SetBranchAddress("Edep_P",        &Edep_P);
        wtree->SetBranchAddress("EdepAcrylic_P", &EdepAcrylic_P);
        wtree->SetBranchAddress("Edep_Pg1",      &Edep_Pg1);
        wtree->SetBranchAddress("Edep_Pg2",      &Edep_Pg2);
        wtree->SetBranchAddress("EdepE_P",       &EdepE_P);
        wtree->SetBranchAddress("EdepEA_P",      &EdepEA_P);
        wtree->SetBranchAddress("EtotE_P",       &EtotE_P);
        wtree->SetBranchAddress("EdepE_Pg1",     &EdepE_Pg1);
        wtree->SetBranchAddress("EdepEA_Pg1",    &EdepEA_Pg1);
        wtree->SetBranchAddress("EtotE_Pg1",     &EtotE_Pg1);
        wtree->SetBranchAddress("EdepE_Pg2",     &EdepE_Pg2);
        wtree->SetBranchAddress("EdepEA_Pg2",    &EdepEA_Pg2);
        wtree->SetBranchAddress("EtotE_Pg2",     &EtotE_Pg2);
        
        ////////////
        
        TH1D *h_Ev_toyMC_input = new TH1D("h_Ev_toyMC_input","h_Ev_toyMC_input",480,0,12);
        
        TFile *roofile_etree = new TFile("./Inputs/HInputs/Data/roofile_etree.root", "recreate");
        TTree *etree = new TTree("etree", "effective entries of wtree");
        
        etree->Branch("cap_Ev",       &cap_Ev,       "cap_Ev/D");
        etree->Branch("cap_Target",   &cap_Target,   "cap_Target/D");
        etree->Branch("X_Ng",         &X_Ng,         "X_Ng/D");
        etree->Branch("Y_Ng",         &Y_Ng,         "Y_Ng/D");
        etree->Branch("Z_Ng",         &Z_Ng,         "Z_Ng/D");
        etree->Branch("Edep_Ng",      &Edep_Ng,      "Edep_Ng/D");
        etree->Branch("EtotScint_Ng", &EtotScint_Ng, "EtotScint_Ng/D");
        etree->Branch("Etot_Ng",      &Etot_Ng,      "Etot_Ng/D");
        etree->Branch("EdepE_Ng",     &EdepE_Ng,     "EdepE_Ng/D");
        etree->Branch("EdepEA_Ng",    &EdepEA_Ng,    "EdepEA_Ng/D");
        etree->Branch("EtotE_Ng",     &EtotE_Ng,     "EtotE_Ng/D");
        etree->Branch("Etot_N",       &Etot_N,       "Etot_N/D");
        etree->Branch("X_P",          &X_P,          "X_P/D");
        etree->Branch("Y_P",          &Y_P,          "Y_P/D");
        etree->Branch("Z_P",          &Z_P,          "Z_P/D");
        etree->Branch("Etot_P",       &Etot_P,       "Etot_P/D");
        etree->Branch("Etot_Pg1",     &Etot_Pg1,     "Etot_Pg1/D");
        etree->Branch("Etot_Pg2",     &Etot_Pg2,     "Etot_Pg2/D");
        etree->Branch("EtotScint_P",  &EtotScint_P,  "EtotScint_P/D");
        etree->Branch("EtotScint_Pg", &EtotScint_Pg, "EtotScint_Pg/D");
        etree->Branch("EtotScint_Pg1", &EtotScint_Pg1, "EtotScint_Pg1/D");
        etree->Branch("EtotScint_Pg2", &EtotScint_Pg2, "EtotScint_Pg2/D");
        etree->Branch("Edep_P",        &Edep_P,        "Edep_P/D");
        etree->Branch("EdepAcrylic_P", &EdepAcrylic_P, "EdepAcrylic_P/D");
        etree->Branch("Edep_Pg1",      &Edep_Pg1,      "Edep_Pg1/D");
        etree->Branch("Edep_Pg2",      &Edep_Pg2,      "Edep_Pg2/D");
        etree->Branch("EdepE_P",       &EdepE_P,       "EdepE_P/D");
        etree->Branch("EdepEA_P",      &EdepEA_P,      "EdepEA_P/D");
        etree->Branch("EtotE_P",       &EtotE_P,       "EtotE_P/D");
        etree->Branch("EdepE_Pg1",     &EdepE_Pg1,     "EdepE_Pg1/D");
        etree->Branch("EdepEA_Pg1",    &EdepEA_Pg1,    "EdepEA_Pg1/D");
        etree->Branch("EtotE_Pg1",     &EtotE_Pg1,     "EtotE_Pg1/D");
        etree->Branch("EdepE_Pg2",     &EdepE_Pg2,     "EdepE_Pg2/D");
        etree->Branch("EdepEA_Pg2",    &EdepEA_Pg2,    "EdepEA_Pg2/D");
        etree->Branch("EtotE_Pg2",     &EtotE_Pg2,     "EtotE_Pg2/D");
        
        for(long ientry=0; ientry<entries_wtree; ientry++)
        {
            wtree->GetEntry(ientry);
            
            cout.precision(3);
            if(ientry%20000==0)
                cout<<" ---> processing MC spectrum "<<ientry*100./entries_wtree<<"%"<<endl;
            
            //////
            
            if(cap_Target<0) continue;
            if(Edep_Ng<0) continue;
            if(EtotScint_Ng<0) continue;
            if(Etot_Ng<0) continue;
            if(EdepE_Ng<0) continue;
            if(EdepEA_Ng<0) continue;
            if(EtotE_Ng<0) continue;
            if(Etot_N<0) continue;
            if(Etot_P<0) continue;
            if(Etot_Pg1<0) continue;
            if(Etot_Pg2<0) continue;
            if(EtotScint_P<0) continue;
            if(EtotScint_Pg<0) continue;
            if(EtotScint_Pg1<0) continue;
            if(EtotScint_Pg2<0) continue;
            if(Edep_P<0) continue;
            if(EdepAcrylic_P<0) continue;
            if(Edep_Pg1<0) continue;
            if(Edep_Pg2<0) continue;
            if(EdepE_P<0) continue;
            if(EdepEA_P<0) continue;
            if(EtotE_P<0) continue;
            if(EdepE_Pg1<0) continue;
            if(EdepEA_Pg1<0) continue;
            if(EtotE_Pg1<0) continue;
            if(EdepE_Pg2<0) continue;
            if(EdepEA_Pg2<0) continue;
            if(EtotE_Pg2<0) continue;
            
            //////
            
            cap_Ev = Etot_P +Etot_N + DeltaM;
            
            h_Ev_toyMC_input->Fill(cap_Ev);
#ifndef ReactorShapeinToy
            etree->Fill();//To fill a flat spectrum to produce a response matrix with nearly equal statistics in every column and about 3 times more events in total
#endif
        }
        
#ifdef ReactorShapeinToy        ///////// Uncomment #ReactorShapeinToy to include reactor shape in the toy:
        
        Double_t max_h_Ev_toyMC_input = h_Ev_toyMC_input->GetBinContent( h_Ev_toyMC_input->GetMaximumBin() );
        Double_t max_h_Ev_normal = h_Ev_normal->GetBinContent( h_Ev_normal->GetMaximumBin() );
        h_Ev_normal->Scale(max_h_Ev_toyMC_input/max_h_Ev_normal);
        
        for(long ientry=0; ientry<entries_wtree; ientry++)
        {
            wtree->GetEntry(ientry);
            
            cout.precision(3);
            if(ientry%20000==0)
                cout<<" ---> processing reactor shape spectrum "<<ientry*100./entries_wtree<<"%"<<endl;
            
            //////
            
            if(cap_Target<0) continue;
            if(Edep_Ng<0) continue;
            if(EtotScint_Ng<0) continue;
            if(Etot_Ng<0) continue;
            if(EdepE_Ng<0) continue;
            if(EdepEA_Ng<0) continue;
            if(EtotE_Ng<0) continue;
            if(Etot_N<0) continue;
            if(Etot_P<0) continue;
            if(Etot_Pg1<0) continue;
            if(Etot_Pg2<0) continue;
            if(EtotScint_P<0) continue;
            if(EtotScint_Pg<0) continue;
            if(EtotScint_Pg1<0) continue;
            if(EtotScint_Pg2<0) continue;
            if(Edep_P<0) continue;
            if(EdepAcrylic_P<0) continue;
            if(Edep_Pg1<0) continue;
            if(Edep_Pg2<0) continue;
            if(EdepE_P<0) continue;
            if(EdepEA_P<0) continue;
            if(EtotE_P<0) continue;
            if(EdepE_Pg1<0) continue;
            if(EdepEA_Pg1<0) continue;
            if(EtotE_Pg1<0) continue;
            if(EdepE_Pg2<0) continue;
            if(EdepEA_Pg2<0) continue;
            if(EtotE_Pg2<0) continue;
            
            //////
            
            cap_Ev = Etot_P +Etot_N +DeltaM;
            if( cap_Ev>12 ) continue;
            
            Int_t ibin = h_Ev_toyMC_input->FindBin( cap_Ev );
            Double_t prob = h_Ev_normal->GetBinContent(ibin);
            
            Int_t seed_rand = 0;
            
            if( ientry%100==0 )
            {
                seed_generator += 1;
                seed_rand = 2 *seed_generator +1;
                
                gRandom3->SetSeed(seed_rand);
            }
            
            Double_t rand = gRandom3->Uniform(0,max_h_Ev_toyMC_input);
            if( rand>prob ) continue;
            
            etree->Fill();
        }
#endif
        etree->Write();
        roofile_etree->Close();
        
#endif
        
        //////////////////////////////////////////////////////////////////////////
        // The etree file has to have been generated before running this code:
        //////////////////////////////////////////////////////////////////////////
        //    const long MaxEntries = 50;
        
        //#ifdef UseDelayInformation
        //    Double_t X_NG[MaxEntries];
        //    Double_t Y_NG[MaxEntries];
        //    Double_t Z_NG[MaxEntries];
        //
        //    Double_t Edep_NG[MaxEntries];
        //    Double_t Etot_NG[MaxEntries];
        //    Double_t EdepEA_NG[MaxEntries];
        //    Double_t EtotE_NG[MaxEntries];
        //    Double_t Etot_Neutron[MaxEntries];
        //#endif
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
        
        TFile *roofile_etree_read = new TFile("./Inputs/HInputs/Data/roofile_etree.root", "read");
        TTree *etree_read = (TTree*)roofile_etree_read->Get("etree");
        long entries_etree_read = etree_read->GetEntries();
        
        etree_read->SetBranchAddress("cap_Ev",       &cap_Ev);
        etree_read->SetBranchAddress("cap_Target",   &cap_Target);
        etree_read->SetBranchAddress("X_Ng",         &X_Ng);
        etree_read->SetBranchAddress("Y_Ng",         &Y_Ng);
        etree_read->SetBranchAddress("Z_Ng",         &Z_Ng);
        etree_read->SetBranchAddress("Edep_Ng",      &Edep_Ng);
        etree_read->SetBranchAddress("EtotScint_Ng", &EtotScint_Ng);
        etree_read->SetBranchAddress("Etot_Ng",      &Etot_Ng);
        etree_read->SetBranchAddress("EdepE_Ng",     &EdepE_Ng);
        etree_read->SetBranchAddress("EdepEA_Ng",    &EdepEA_Ng);
        etree_read->SetBranchAddress("EtotE_Ng",     &EtotE_Ng);
        etree_read->SetBranchAddress("Etot_N",       &Etot_N);
        etree_read->SetBranchAddress("X_P",          &X_P);
        etree_read->SetBranchAddress("Y_P",          &Y_P);
        etree_read->SetBranchAddress("Z_P",          &Z_P);
        etree_read->SetBranchAddress("Etot_P",       &Etot_P);
        etree_read->SetBranchAddress("Etot_Pg1",     &Etot_Pg1);
        etree_read->SetBranchAddress("Etot_Pg2",     &Etot_Pg2);
        etree_read->SetBranchAddress("EtotScint_P",  &EtotScint_P);
        etree_read->SetBranchAddress("EtotScint_Pg", &EtotScint_Pg);
        etree_read->SetBranchAddress("EtotScint_Pg1", &EtotScint_Pg1);
        etree_read->SetBranchAddress("EtotScint_Pg2", &EtotScint_Pg2);
        etree_read->SetBranchAddress("Edep_P",        &Edep_P);
        etree_read->SetBranchAddress("EdepAcrylic_P", &EdepAcrylic_P);
        etree_read->SetBranchAddress("Edep_Pg1",      &Edep_Pg1);
        etree_read->SetBranchAddress("Edep_Pg2",      &Edep_Pg2);
        etree_read->SetBranchAddress("EdepE_P",       &EdepE_P);
        etree_read->SetBranchAddress("EdepEA_P",      &EdepEA_P);
        etree_read->SetBranchAddress("EtotE_P",       &EtotE_P);
        etree_read->SetBranchAddress("EdepE_Pg1",     &EdepE_Pg1);
        etree_read->SetBranchAddress("EdepEA_Pg1",    &EdepEA_Pg1);
        etree_read->SetBranchAddress("EtotE_Pg1",     &EtotE_Pg1);
        etree_read->SetBranchAddress("EdepE_Pg2",     &EdepE_Pg2);
        etree_read->SetBranchAddress("EdepEA_Pg2",    &EdepEA_Pg2);
        etree_read->SetBranchAddress("EtotE_Pg2",     &EtotE_Pg2);
        
        ///////////////////////////////
        ///////////////////////////////
        ///////////////////////////////
        ///////////////////////////////
        
        Double_t LY_E_P;
        Double_t LY_E_Pg1;
        Double_t LY_E_Pg2;
        Double_t LY_E_N;
        
        Double_t LY_E_P_Sum;
        Double_t Opt_E_P_Sum;
        Double_t FEE_E_P_Sum;
        Double_t Scale_E_P_Sum;
        Double_t Res_E_P_Sum;
        Double_t Eff_E_P_Sum;
        Double_t Erec_P;

#ifdef UseDelayInformation
        Double_t LY_E_Ng;
        Double_t Opt_E_Ng;
        Double_t FEE_E_Ng;
        Double_t Scale_E_Ng;
        Double_t Res_E_Ng;
        Double_t Eff_E_Ng;
        Double_t Erec_Ng;
#endif
        
#ifdef SaveTree
        
        Int_t pcell = 0;
        Int_t dcell = 0;
        
        TFile *roofile_toy = new TFile("roofile_toy.root", "recreate");
        TTree *toy = new TTree("toy", "toyMC result");
        
        toy->Branch("cap_Ev",       &cap_Ev,       "cap_Ev/D");
        toy->Branch("cap_Target",   &cap_Target,   "cap_Target/D");
        toy->Branch("pcell",        &pcell,        "pcell/I");
        toy->Branch("dcell",        &dcell,        "dcell/I");
        toy->Branch("X_Ng",         &X_Ng,         "X_Ng/D");
        toy->Branch("Y_Ng",         &Y_Ng,         "Y_Ng/D");
        toy->Branch("Z_Ng",         &Z_Ng,         "Z_Ng/D");
        toy->Branch("Edep_Ng",      &Edep_Ng,      "Edep_Ng/D");
        toy->Branch("EtotScint_Ng", &EtotScint_Ng, "EtotScint_Ng/D");
        toy->Branch("Etot_Ng",      &Etot_Ng,      "Etot_Ng/D");
        toy->Branch("EdepE_Ng",     &EdepE_Ng,     "EdepE_Ng/D");
        toy->Branch("EdepEA_Ng",    &EdepEA_Ng,    "EdepEA_Ng/D");
        toy->Branch("EtotE_Ng",     &EtotE_Ng,     "EtotE_Ng/D");
        toy->Branch("Etot_N",       &Etot_N,       "Etot_N/D");
        toy->Branch("X_P",          &X_P,          "X_P/D");
        toy->Branch("Y_P",          &Y_P,          "Y_P/D");
        toy->Branch("Z_P",          &Z_P,          "Z_P/D");
        toy->Branch("Etot_P",       &Etot_P,       "Etot_P/D");
        toy->Branch("Etot_Pg1",     &Etot_Pg1,     "Etot_Pg1/D");
        toy->Branch("Etot_Pg2",     &Etot_Pg2,     "Etot_Pg2/D");
        toy->Branch("EtotScint_P",  &EtotScint_P,  "EtotScint_P/D");
        toy->Branch("EtotScint_Pg", &EtotScint_Pg, "EtotScint_Pg/D");
        toy->Branch("EtotScint_Pg1", &EtotScint_Pg1, "EtotScint_Pg1/D");
        toy->Branch("EtotScint_Pg2", &EtotScint_Pg2, "EtotScint_Pg2/D");
        toy->Branch("Edep_P",        &Edep_P,        "Edep_P/D");
        toy->Branch("EdepAcrylic_P", &EdepAcrylic_P, "EdepAcrylic_P/D");
        toy->Branch("Edep_Pg1",      &Edep_Pg1,      "Edep_Pg1/D");
        toy->Branch("Edep_Pg2",      &Edep_Pg2,      "Edep_Pg2/D");
        toy->Branch("EdepE_P",       &EdepE_P,       "EdepE_P/D");
        toy->Branch("EdepEA_P",      &EdepEA_P,      "EdepEA_P/D");
        toy->Branch("EtotE_P",       &EtotE_P,       "EtotE_P/D");
        toy->Branch("EdepE_Pg1",     &EdepE_Pg1,     "EdepE_Pg1/D");
        toy->Branch("EdepEA_Pg1",    &EdepEA_Pg1,    "EdepEA_Pg1/D");
        toy->Branch("EtotE_Pg1",     &EtotE_Pg1,     "EtotE_Pg1/D");
        toy->Branch("EdepE_Pg2",     &EdepE_Pg2,     "EdepE_Pg2/D");
        toy->Branch("EdepEA_Pg2",    &EdepEA_Pg2,    "EdepEA_Pg2/D");
        toy->Branch("EtotE_Pg2",     &EtotE_Pg2,     "EtotE_Pg2/D");
        
        toy->Branch("LY_E_P",        &LY_E_P,        "LY_E_P/D");
        toy->Branch("LY_E_Pg1",      &LY_E_Pg1,      "LY_E_Pg1/D");
        toy->Branch("LY_E_Pg2",      &LY_E_Pg2,      "LY_E_Pg2/D");
        toy->Branch("LY_E_N",    &LY_E_N,    "LY_E_N/D");
        toy->Branch("LY_E_P_Sum",    &LY_E_P_Sum,    "LY_E_P_Sum/D");
        toy->Branch("Opt_E_P_Sum",   &Opt_E_P_Sum,   "Opt_E_P_Sum/D");
        toy->Branch("FEE_E_P_Sum",   &FEE_E_P_Sum,   "FEE_E_P_Sum/D");
        toy->Branch("Scale_E_P_Sum", &Scale_E_P_Sum, "Scale_E_P_Sum/D");
        toy->Branch("Res_E_P_Sum",   &Res_E_P_Sum,   "Res_E_P_Sum/D");
        toy->Branch("Eff_E_P_Sum",   &Eff_E_P_Sum,   "Eff_E_P_Sum/D");
        toy->Branch("Erec_P",   &Erec_P,   "Erec_P/D");

#ifdef UseDelayInformation
        toy->Branch("LY_E_Ng",    &LY_E_Ng,    "LY_E_Ng/D");`
        toy->Branch("Opt_E_Ng",   &Opt_E_Ng,   "Opt_E_Ng/D");
        toy->Branch("FEE_E_Ng",   &FEE_E_Ng,   "FEE_E_Ng/D");
        toy->Branch("Scale_E_Ng", &Scale_E_Ng, "Scale_E_Ng/D");
        toy->Branch("Res_E_Ng",   &Res_E_Ng,   "Res_E_Ng/D");
        toy->Branch("Eff_E_Ng",   &Eff_E_Ng,   "Eff_E_Ng/D");
        toy->Branch("Erec_Ng",   &Erec_Ng,   "Erec_Ng/D");

#endif
        
#endif
        
        //Load tree data
        //    for(long ientry=0; ientry<entries_etree_read; ientry++)
        //    {
        //
        //#ifdef UseDelayInformation
        //        X_NG[ientry] = X_Ng;
        //        Y_NG[ientry] = Y_Ng;
        //        Z_NG[ientry] = Z_Ng;
        //
        //        Edep_NG[ientry] = Edep_Ng;
        //        Etot_NG[ientry] = Etot_Ng;
        //        EdepEA_NG[ientry] = EdepEA_Ng;
        //        EtotE_NG[ientry] = EtotE_Ng;
        //        Etot_Neutron[ientry]= Etot_N;
        //#endif
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
        
        Double_t P_Scint=0, Pg1_Scint=0, Pg2_Scint=0, Ng_Scint=0;
        Double_t AdSimple_P=0, AdSimple_Ng=0;
        Double_t usr_opt_attenuation_P=0, usr_pmt_coverage_P=0, usr_opt_attenuation_Ng=0, usr_pmt_coverage_Ng=0;
        Double_t energy_sigma=0, R2_AA=0, R2_BB=0, R_average=0;
        
        //Efficiency maps:
        TFile *roofile_h2d_ep_ratio2center = new TFile("./Inputs/HInputs/Data/cell_eff/h2d_ep_ratio2center.root", "read");
        TH2D *h2d_Ep_ratio2center = (TH2D*)roofile_h2d_ep_ratio2center->Get("h2d_ep_ratio2center");
        
        TFile *roofile_h2d_ed_ratio2center = new TFile("./Inputs/HInputs/Data/cell_eff/h2d_ed_ratio2center.root", "read");
        TH2D *h2d_Ed_ratio2center = (TH2D*)roofile_h2d_ed_ratio2center->Get("h2d_ed_ratio2center");
        
        seed_generator_uncorr = 2863311530; //=(2*4294967295/3) for uncorrelated systematics, chose 2/3*(maxseed) to make it different to 'seed_generator' as a precaution so I don't use the same seeds.
        
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
        
        Double_t IntegralEvents[NADs][VolumeX][VolumeY];
        Double_t ADIntegralEvents[NADs];
        
        //Use the loaded tree data in each AD in an independent way:
        for(Int_t AD = 0; AD<NADs;AD++)
        {
            seed_generator = 1;//this way all matrices will be random but have a common nominal spectrum
            seed_generator_corr = 1431655765; //=(4294967295/3) for correlated systematics, chose maxseed/3 to make it different to 'seed_generator' as a precaution so I don't use the same seeds.
           
            Double_t distPD=0;

            //Process tree
            for(long ientry=0; ientry<entries_etree_read; ientry++)
            {
                etree_read->GetEntry(ientry);
                
                cout.precision(3);
                if(ientry%20000==0)
                    cout<<" ---> processing response matrix in AD" << AD << " " <<ientry*100./entries_etree_read<<"%"<<endl;
                
                Int_t seed_rand = 0;
                Int_t seed_corr =0;
                Int_t seed_uncorr = 0;
          
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
                GausNL= 0;
                GausReso=0;
                GausResoCorr=0;
                GausEff= 0;
                
                 if(mode!=0)//If not nominal values, vary systematic for each entry and for each AD (uncorrelated)
                {
                    if(RelativeEnergyScaleMatrix)
                    {
                        GausRelative = RandomSysUncorr->Gaus(0,1);
                    }
                    if(IAVMatrix)
                    {
                        GausIAV = RandomSysUncorr->Gaus(0,1);
                    }
                    if(NLMatrix)
                    {
                        GausNL = RandomSysUncorr->Gaus(0,1);
                    }
                    if(ResolutionMatrix)
                    {
                        GausReso = RandomSysUncorr->Gaus(0,1);
                        GausResoCorr = RandomSysCorr->Gaus(0,1);
                    }
                    if(EfficiencyMatrix)
                    {
                        GausEff = RandomSysUncorr->Gaus(0,1);
                    }

                    //factor that is multiplied by the random error so it's added to the nominal by using:  Value = NominalValue + gRandom3->Gauss(0,1)*1sigmaError;
                }
                if(ientry%1000000==0)//show only a few
                {
                    std::cout << "IAV uncorrelated should be the different for different ads and same event: " << GausIAV << "ad: " << AD << " - event: " << ientry <<  std::endl;

                    std::cout << "gaus reso correlated should be the same for different ads and same event: " << GausResoCorr << "ad: " << AD << " - event: " << ientry <<  std::endl;
                    std::cout << "gaus reso uncorrelated should be different for different ads and same event" << GausReso << "ad: " << AD << " - event: " << ientry <<  std::endl;
                }
                
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
                
                /// case04: default ? edit
                /// case03: default
                Double_t RandomEdep_P,RandomEdep_Pg1,RandomEdep_Pg2;
                //Iav uncertainty considered for prompt signal. 0.1% bin-to-bin uncorrelated.
                RandomEdep_P = Edep_P*(1+GausIAV*IAVNominalError);
                RandomEdep_Pg1 = Edep_Pg1*(1+GausIAV*IAVNominalError);
                RandomEdep_Pg2 = Edep_Pg2*(1+GausIAV*IAVNominalError);
                //here I am applying the same variation to the three particles at once, one variation for each event.

                if( GausIAV*IAVNominalError < -1) continue;//Very unlikely, but just in case.

                P_Scint = EtotScint_P-RandomEdep_P;
                Pg1_Scint = EtotScint_Pg1-RandomEdep_Pg1;
                Pg2_Scint = EtotScint_Pg2-RandomEdep_Pg2;
                Ng_Scint = EtotScint_Ng-Edep_Ng;
                
                /// prompt
                LY_E_P = EtotScint_P * graph_electron_LY->Eval( EtotScint_P, 0, "s" )
                - P_Scint * graph_electron_LY->Eval( P_Scint, 0, "s" );
                
                LY_E_Pg1 = EtotScint_Pg1 * graph_gamma_LY->Eval( EtotScint_Pg1, 0, "s" )
                - Pg1_Scint * graph_gamma_LY->Eval( Pg1_Scint, 0, "s" );
                
                LY_E_Pg2 = EtotScint_Pg2 * graph_gamma_LY->Eval( EtotScint_Pg2, 0, "s" )
                - Pg2_Scint * graph_gamma_LY->Eval( Pg2_Scint, 0, "s" );
                
                // ------------------------ //
                /// neutron - determined from NuWa for 0-0.18 MeV (does not consider leakage)
                LY_E_N = Etot_N * ( 0.186 + exp(-1.142-126.4*Etot_N) + exp(-1.217-17.90*Etot_N) ) * graph_electron_LY->Eval( 0.2, 0, "s" );
                
                /// prompt SUM
                LY_E_P_Sum = LY_E_P + LY_E_Pg1 + LY_E_Pg2 + LY_E_N;
                //------------------------ //
                
#ifdef UseDelayInformation
                
                /// delayed
                LY_E_Ng = EtotScint_Ng * graph_gamma_LY->Eval( EtotScint_Ng, 0, "s" )
                - Ng_Scint * graph_gamma_LY->Eval( Ng_Scint, 0, "s" );
#endif
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
                 
                 if( local_zbin != 0) continue;
                
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
                
#ifdef UseDelayInformation
                 usr_r2_Ng = (X_Ng*X_Ng+Y_Ng*Y_Ng) * 1e-6;  // mm2 ---> m2
                 usr_z_Ng  = Z_Ng * 1e-3;  // mm ---> m
                
                 global_bin_num = hist_findbin->FindBin(usr_r2_Ng, usr_z_Ng);
                 hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
                 
                 if( local_zbin != 0) continue;
                 
                 usr_opt_attenuation_Ng = hist_map_attenuation->GetBinContent(local_xbin, local_ybin);
                 usr_pmt_coverage_Ng  = hist_map_pmt_coverage->GetBinContent(local_xbin, local_ybin);
                 
                 Opt_E_Ng = LY_E_Ng * usr_opt_attenuation_Ng * usr_pmt_coverage_Ng;  // nH analysis
                 
                // AdSimple_Ng = ( (7.84628 * (1 + 3.41294e-02*usr_r2_Ng) * (1 - 1.21750e-02*usr_z_Ng - 1.64275e-02*usr_z_Ng*usr_z_Ng + 7.33006e-04*pow(usr_z_Ng,3)))/8.05 );  // AdSimple
                // Doc7334(old function), updated from http://dayabay.ihep.ac.cn/tracs/dybsvn/browser/dybgaudi/trunk/Reconstruction/QsumEnergy/src/components/QsumEnergyTool.cc
                
                //   Opt_E_Ng = LY_E_Ng * AdSimple_Ng;
                
#endif
                ///////////////////////////// FEE
                /// prompt
                //if( Opt_E_P_Sum!=Opt_E_P_Sum ) continue;//Maybe eliminate it, it produces problems. See April 16th Logan's email
                FEE_E_P_Sum = Opt_E_P_Sum * graph_electronic->Eval( Opt_E_P_Sum, 0, "s" );
                
                ///////////////////////////// EnergyScale
                /// prompt
                Scale_E_P_Sum = FEE_E_P_Sum *EnergyScale;
#ifdef UseDelayInformation
                ///////////////////////////// FEE
                /// delayed
                //if( Opt_E_Ng!=Opt_E_Ng ) continue;//Maybe eliminate it, it produces problems. See April 16th Logan's email
                FEE_E_Ng = Opt_E_Ng * graph_electronic->Eval( Opt_E_Ng, 0, "s" );
                
                ///////////////////////////// EnergyScale
                /// delayed
                Scale_E_Ng = FEE_E_Ng *EnergyScale;
#endif

                ///////////////////////////// Resolution
                ///////////////////////////// Resolution
                
                Double_t energy_sigma = 0;
                Double_t R2_AA        = 0;
                Double_t R2_BB        = 0;
                Double_t R_average    = 0;
                
                /// prompt
                global_bin_num = hist_findbin->FindBin(usr_r2_P, usr_z_P);
                hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
                
                if( local_zbin != 0) continue;
                /*
                 R2_AA = (local_xbin-1) * R2_binwidth;
                 R2_BB = local_xbin * R2_binwidth;
                 R_average = sqrt( (R2_AA+R2_BB)/2. );  // volume-weighted average radius of cell
                 */
                R_average = sqrt( usr_r2_P );
                
                roofunc_EnergyResolution->SetParameter( 0, R_average );
                energy_sigma = roofunc_EnergyResolution->Eval( Scale_E_P_Sum );
                Res_E_P_Sum = gRandom3->Gaus( Scale_E_P_Sum, energy_sigma );
                
                /// Estimate an "Erec" on which to apply cuts
                ///////////////////////////// Reconstructed Energy
                /// prompt
                Erec_P = Res_E_P_Sum / (usr_opt_attenuation_P *usr_pmt_coverage_P);  // nH analysis
                //Erec_P = Res_E_P_Sum / AdSimple_P;  // AdSimple
                
#ifdef UseDelayInformation
                /// delayed
                global_bin_num = hist_findbin->FindBin(usr_r2_Ng, usr_z_Ng);
                hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
                
                if( local_zbin != 0) continue;
                /*
                 R2_AA = (local_xbin-1) * R2_binwidth;
                 R2_BB = local_xbin * R2_binwidth;
                 R_average = sqrt( (R2_AA+R2_BB)/2. );  // volume-weighted average radius of cell
                 */
                R_average = sqrt( usr_r2_Ng );
                
                roofunc_EnergyResolution->SetParameter( 0, R_average );
                energy_sigma = roofunc_EnergyResolution->Eval( Scale_E_Ng );
                Res_E_Ng = gRandom3->Gaus( Scale_E_Ng, energy_sigma );
                
                /// Estimate an "Erec" on which to apply cuts
                ///////////////////////////// Reconstructed Energy
                /// delayed
                Erec_Ng = Res_E_Ng / (usr_opt_attenuation_Ng * usr_pmt_coverage_Ng);  // nH analysis
                //Erec_Ng = Res_E_Ng / AdSimple_Ng;  // AdSimple
#endif
                ///////////////////////////// Cell
                ///////////////////////////// Cell
                
                ////// Cuts:
                //////
                if( cap_Target==64 ) continue;  // reject nGd events because not properly considered in toyMC input sample
                /// Analysis cuts
                
                if( Res_E_P_Sum!=Res_E_P_Sum) continue;
                
#ifdef UseDelayInformation
                if(Res_E_Ng!=Res_E_Ng ) continue;
                
                if( Erec_Ng<EdelayCutLow || Erec_Ng>EdelayCutHigh ) continue;
#endif
                //if( TcapCutHigh>0 && (?>TcapCutHigh || ?<TcapCutLow) ) continue;  // time information is currently not saved in the toyMC input sample
                if( Erec_P<EpromptCutLow || Erec_P>EpromptCutHigh ) continue;
                if( distCut>0 ) {
                    distPD = sqrt( (X_Ng-X_P)*(X_Ng-X_P) + (Y_Ng-Y_P)*(Y_Ng-Y_P) + (Z_Ng-Z_P)*(Z_Ng-Z_P) );
                    if( distPD>distCut ) continue;
                }
                //Efficiency correction:
                Double_t Eff_p_content;
                ////// prompt
                usr_r2_P = (X_P*X_P+Y_P*Y_P) * 1e-6;  // mm2 ---> m2
                if( radialCut>0 && usr_r2_P>radialCut ) continue;
                usr_z_P  = Z_P * 1e-3;  // mm ---> m
                global_bin_num = hist_findbin->FindBin(usr_r2_P, usr_z_P);
                hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
               // hEp_cc[global_bin_num]->Fill( Res_E_P_Sum );
                if( local_zbin==0 && local_xbin>=1 && local_xbin<=R2_binnum && local_ybin>=1 && local_ybin<=Z_binnum )
                {
                    Eff_p_content = h2d_Ep_ratio2center->GetBinContent(local_xbin, local_ybin);
                }
#ifdef UseDelayInformation
                //Efficiency correction:
                Double_t Eff_d_content;
                ////// delayed
                usr_r2_Ng = (X_Ng*X_Ng+Y_Ng*Y_Ng) * 1e-6;  // mm2 ---> m2
                if( radialCut>0 && usr_r2_Ng>radialCut ) continue;
                usr_z_Ng  = Z_Ng * 1e-3;  // mm ---> m
                global_bin_num = hist_findbin->FindBin(usr_r2_Ng, usr_z_Ng);
                hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
               // hEd_cc[global_bin_num]->Fill( Res_E_Ng );
                if( local_zbin==0 && local_xbin>=1 && local_xbin<=R2_binnum && local_ybin>=1 && local_ybin<=Z_binnum )
                {
                    Eff_d_content = h2d_Ed_ratio2center->GetBinContent(local_xbin, local_ybin);
                }
#endif
                Int_t VolumeXbin, VolumeYbin;

#ifdef UseVolumes
                //Logic to detect to translate cell bins to volume bins (check if the cell is inside the GdLs or the Ls volume)
                VolumeYbin = 1;

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
                
                //Fill histograms:
                TruePredictionH[AD][VolumeXbin-1][VolumeYbin-1]->Fill(cap_Ev,Eff_p_content);
                VisiblePredictionH[AD][VolumeXbin-1][VolumeYbin-1]->Fill(Res_E_P_Sum,Eff_p_content);
                
#ifdef UseDelayInformation
                DelayedVisiblePredictionH[AD][VolumeXbin-1][VolumeYbin-1]->Fill(Res_E_Ng,Eff_d_content);
#endif
                HighResoMatrixH[AD][VolumeXbin-1][VolumeYbin-1]->Fill(cap_Ev,Res_E_P_Sum,Eff_p_content);//Fine grid
                MatrixH[AD][VolumeXbin-1][VolumeYbin-1]->Fill(cap_Ev,Res_E_P_Sum,Eff_p_content);//neutrino energy vs visible energy
                TransMatrixH[AD][VolumeXbin-1][VolumeYbin-1]->Fill(Res_E_P_Sum,cap_Ev,Eff_p_content);
                HighResoTransMatrixH[AD][VolumeXbin-1][VolumeYbin-1]->Fill(Res_E_P_Sum,cap_Ev,Eff_p_content);//Fine grid
                
#ifdef SaveTree
                toy->Fill();
#endif
            }//events
            
//            //Efficiency map:
//            
//            double entries=0, Effcontent=0;
//            
//            for(int idx=0; idx<MaxCellNum; idx++)
//            {
//                entries = 0;
//                entries = hEp_cc[idx]->GetEntries();
//                if( entries>10 )
//                {
//                    // hEp_cc[idx]->Write();
//                    hEp_cc_clone[idx]->Add(hEp_cc[idx]);
//                }
//                
//                entries = hEd_cc[idx]->GetEntries();
//                if( entries>10 )
//                {
//                    // hEd_cc[idx]->Write();
//                    hEd_cc_clone[idx]->Add(hEd_cc[idx]);
//                }
//                
//                //////
//                hist_findbin->GetBinXYZ(idx, local_xbin, local_ybin, local_zbin);
//                if( local_zbin==0 && local_xbin>=1 && local_xbin<=R2_binnum && local_ybin>=1 && local_ybin<=Z_binnum )
//                {
//                    Effcontent = 0;
//                    
//                    //if(local_ybin<=1 || local_ybin>=Z_binnum)  continue;  // exclude bottom and top rows from columns
//                    
//                    entries = hEp_cc_clone[idx]->Integral(1,hEp_bin);
//                    Effcontent = h2d_Ep_ratio2center->GetBinContent(local_xbin, local_ybin);
//                    hEp_cc_clone[idx]->Scale(Effcontent/entries);
//                    hEp_cl[local_xbin]->Add(hEp_cc_clone[idx]);//Efficiency corrected histogram
//                    
//                    entries = hEd_cc_clone[idx]->Integral(1,hEd_bin);
//                    Effcontent = h2d_Ed_ratio2center->GetBinContent(local_xbin, local_ybin);
//                    hEd_cc_clone[idx]->Scale(Effcontent/entries);
//                    hEd_cl[local_xbin]->Add(hEd_cc_clone[idx]);//Efficiency corrected histogram
//                }
//            }
//            
//            //Reset histograms
//            for(int idx=0; idx<MaxCellNum; idx++)
//            {
//                hEp_cc[idx]->Reset();
//                hEd_cc[idx]->Reset();
//                hEd_cc_clone[idx]->Reset();
//                hEp_cc_clone[idx]->Reset();
//            }
//            for(int idx=0; idx<MaxColumnNum; idx++)
//            {
//                hEp_cl[idx]->Reset();
//                hEd_cl[idx]->Reset();
//            }
            
            
            //Save number of events in each cell to properly scale the cells afterwards:
            
            ADIntegralEvents[AD] = 0;//initialize
            
            for(int idx=0; idx<VolumeX; idx++)
            {
                for(int idy=0; idy<VolumeY; idy++)
                {
                    IntegralEvents[AD][idx][idy] = VisiblePredictionH[AD][idx][idy]->Integral();
                    ADIntegralEvents[AD] = ADIntegralEvents[AD] + IntegralEvents[AD][idx][idy];
                }
            }
        }//ADs
        
        delete h2d_Ep_ratio2center;
        delete h2d_Ed_ratio2center;
        
        //Save event ratios:
        
        string outputS;
    
        if(mode==0)
        {
            outputS= "./Inputs/HInputs/NominalToyMCEventRatio.txt";
        }
        else
        {
             if(RelativeEnergyScaleMatrix)
             {
                 outputS= "./Inputs/HInputs/RelativeEnergyScaleToyMCEventRatio.txt";
             }
             if(IAVMatrix)
             {
                 outputS= "./Inputs/HInputs/IAVToyMCEventRatio.txt";
                 
             }
             else if(NLMatrix)
             {
                 outputS= "./Inputs/HInputs/NLToyMCEventRatio.txt";
                 
             }
             else if(ResolutionMatrix)
             {
                 outputS= "./Inputs/HInputs/ResolutionToyMCEventRatio.txt";
                 
             }
             else if(EfficiencyMatrix)
             {
                 outputS= "./Inputs/HInputs/EfficiencyToyMCEventRatio.txt";
                 
             }
        }
        
        FILE *f = fopen(outputS.c_str(), "w");

        if (f == NULL)
        {
            printf("Error opening file!\n");
            exit(1);
        }
        
        
        for(Int_t AD = 0; AD<NADs; AD++)
        {
            fprintf(f, "#ToyEvents in each AD\n");

            fprintf(f, "%d %f\n", AD, ADIntegralEvents[AD]);

            fprintf(f, "#Events in each cell for AD %d\n", AD);

            for(int idx=0; idx<VolumeX; idx++)
            {
                for(int idy=0; idy<VolumeY; idy++)
                {
                    fprintf(f, "%d %d %d %f\n", AD, idx, idy, IntegralEvents[AD][idx][idy] );
                }
            }
        }

        fclose(f);

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
        
        if(mode!=0)
        {
            delete gRandom3;
            delete RandomSysUncorr;
            delete RandomSysCorr;
        }
        
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
                MatrixH[i][0][0]->Draw("colz");
                
                NNFineMatrixC->cd(i+1);
                HighResoMatrixH[i][0][0]->Draw("colz");
                
                NNTransC->cd(i+1);
                TransMatrixH[i][0][0]->Draw("colz");
                
                NNFineTransC->cd(i+1);
                HighResoTransMatrixH[i][0][0]->Draw("colz");
            }
            
            NNFineMatrixC->Print("./Images/Hydrogen/ResponseMatrices/NoNormFineHydrogenResponseMatrix.eps");
            NNMatrixC->Print("./Images/Hydrogen/ResponseMatrices/NoNormHydrogenResponseMatrix.eps");
            NNTransC->Print("./Images/Hydrogen/ResponseMatrices/NoNormTransposeHydrogenMatrix.eps");
            NNFineTransC->Print("./Images/Hydrogen/ResponseMatrices/NoNormFineTransposeHydrogenMatrix.eps");
            
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
                            Norma[i] = Norma[i]+MatrixH[AD][idx][idy]->GetBinContent(i+1,j+1);// true->vis
                        }
                    }
                    
                    for(Int_t i=0;i<n_evis_bins;i++)
                    {
                        NormaTrans[i]=0;
                        
                        for(Int_t j=0;j<n_etrue_bins;j++)
                        {
                            NormaTrans[i] = NormaTrans[i]+TransMatrixH[AD][idx][idy]->GetBinContent(i+1,j+1);// vis->true
                        }
                    }
                    
                    for (Int_t i = 0; i < n_etrue_bins; i++)
                    {
                        for (Int_t j = 0; j < n_evis_bins; j++)
                        {
                            if(Norma[i]!=0)
                            {
                                MatrixH[AD][idx][idy]->SetBinContent(i+1,j+1,MatrixH[AD][idx][idy]->GetBinContent(i+1,j+1)/Norma[i]);//true->vis
                            }
                        }
                    }
                    
                    for (Int_t i = 0; i < n_evis_bins; i++)
                    {
                        for (Int_t j = 0; j < n_etrue_bins; j++)
                        {
                            if(NormaTrans[i]!=0)
                            {
                                TransMatrixH[AD][idx][idy]->SetBinContent(i+1,j+1,TransMatrixH[AD][idx][idy]->GetBinContent(i+1,j+1)/NormaTrans[i]);
                            }
                            
                        }
                    }
                    
                    for(Int_t i=0;i<MatrixBins;i++)
                    {
                        HighResoNorma[i]=0;
                        HighResoNormaTrans[i]=0;
                        
                        for(Int_t j=0;j<MatrixBins;j++)
                        {
                            HighResoNorma[i] = HighResoNorma[i]+HighResoMatrixH[AD][idx][idy]->GetBinContent(i+1,j+1);// true->vis
                            HighResoNormaTrans[i] = HighResoNormaTrans[i]+HighResoTransMatrixH[AD][idx][idy]->GetBinContent(i+1,j+1);// vis->true
                            
                        }
                    }
                    
                    for(Int_t i=0;i<MatrixBins;i++)
                    {
                        for(Int_t j=0;j<MatrixBins;j++)
                        {
                            if(HighResoNorma[i]!=0)
                            {
                                HighResoMatrixH[AD][idx][idy]->SetBinContent(i+1,j+1,HighResoMatrixH[AD][idx][idy]->GetBinContent(i+1,j+1)/HighResoNorma[i]);//true->vis
                            }
                            
                            if(HighResoNormaTrans[i]!=0)
                            {
                                HighResoTransMatrixH[AD][idx][idy]->SetBinContent(i+1,j+1,HighResoTransMatrixH[AD][idx][idy]->GetBinContent(i+1,j+1)/HighResoNormaTrans[i]);
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
                        TruePredictionH[i][idx][idy]->Draw();
                        
                        VisPredictionC->cd(i+1);
                        VisiblePredictionH[i][idx][idy]->Draw();
                        
                        DelayedVisPredictionC->cd(i+1);
                        DelayedVisiblePredictionH[i][idx][idy]->Draw();
                        
                        MatrixC->cd(i+1);
                        MatrixH[i][idx][idy]->Draw("colz");
                        
                        FineMatrixC->cd(i+1);
                        HighResoMatrixH[i][idx][idy]->Draw("colz");
                        
                        TransC->cd(i+1);
                        TransMatrixH[i][idx][idy]->Draw("colz");
                        
                        FineTransC->cd(i+1);
                        HighResoTransMatrixH[i][idx][idy]->Draw("colz");
                    }
                    
                    PredictionC->Print("./Images/Hydrogen/ToyMCOutputs/ToyMCTrueHydrogenPrediction.eps");
                    VisPredictionC->Print("./Images/Hydrogen/ToyMCOutputs/ToyMCVisibleHydrogenPrediction.eps");
                    DelayedVisPredictionC->Print("./Images/Hydrogen/ToyMCOutputs/ToyMCDelayedVisibleHydrogenPrediction.eps");
                    FineMatrixC->Print("./Images/Hydrogen/ResponseMatrices/FineHydrogenResponseMatrix.eps");
                    MatrixC->Print("./Images/Hydrogen/ResponseMatrices/HydrogenResponseMatrix.eps");
                    TransC->Print("./Images/Hydrogen/ResponseMatrices/TransposeHydrogenMatrix.eps");
                    FineTransC->Print("./Images/Hydrogen/ResponseMatrices/FineTransposeHydrogenMatrix.eps");
                    
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
        if(mode==0)
        {
            TFile* SaveMatrix = new TFile("./ResponseMatrices/Hydrogen/NominalResponseMatrix.root","recreate");
            
            for(int idx=0; idx<VolumeX; idx++)
            {
                for(int idy=0; idy<VolumeY; idy++)
                {
                    for(Int_t AD = 0; AD<NADs; AD++)//Save different AD matrices to produce the covariance matrices.
                    {
                        HighResoMatrixH[AD][idx][idy]->Write(Form("FineEvisEnu%i,Cell%i,%i",AD+1,idx,idy));//Save the matrices.
                        MatrixH[AD][idx][idy]->Write(Form("EvisEnu%i,Cell%i,%i",AD+1,idx,idy));//Save the matrices.
                        TransMatrixH[AD][idx][idy]->Write(Form("EnuEvis%i,Cell%i,%i",AD+1,idx,idy));
                    }
                    
                }
            }
            delete SaveMatrix;
        }
        else
        {
            TFile* SaveMatrix = new TFile("./ResponseMatrices/Hydrogen/RandomResponseMatrix.root","recreate");
           
            for(int idx=0; idx<VolumeX; idx++)
            {
                for(int idy=0; idy<VolumeY; idy++)
                {
                    for(Int_t AD = 0; AD<NADs; AD++)//Save different AD matrices to produce the covariance matrices.
                    {
                        HighResoMatrixH[AD][idx][idy]->Write(Form("FineEvisEnu%i,Cell%i,%i",AD+1,idx,idy));//Save the matrices.
                        MatrixH[AD][idx][idy]->Write(Form("EvisEnu%i,Cell%i,%i",AD+1,idx,idy));//Save the matrices.
                        TransMatrixH[AD][idx][idy]->Write(Form("EnuEvis%i,Cell%i,%i",AD+1,idx,idy));
                    }
                    
                }
            }
            delete SaveMatrix;
        }
#ifdef SaveTree
        toy->Write();
        roofile_toy->Close();
#endif

#ifndef ReDoNominal
    }
    else//if mode == 0
    {
        TFile* NominalMatrixF = new TFile("./ResponseMatrices/Hydrogen/NominalResponseMatrix.root");
        
        for(int idx=0; idx<VolumeX; idx++)
        {
            for(int idy=0; idy<VolumeY; idy++)
            {
                for(Int_t AD = 0; AD<NADs; AD++)//Save different AD matrices to produce the covariance matrices.
                {
                    HighResoMatrixH[AD][idx][idy] = (TH2D*)NominalMatrixF->Get(Form("FineEvisEnu%i,Cell%i,%i",AD+1,idx,idy));//Save the matrices.
                }
            }
        }
        delete NominalMatrixF;
    }
    
#endif
}

void nHToyMC :: GeneratenHResponseMatrix()
{
    
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
    TFile* file_IHEP_NL_models = new TFile("./Inputs/HInputs/Data/IHEP_NL_models/Model1.root", "read");
    graph_electron_LY = (TGraph*)file_IHEP_NL_models->Get("electronScint");// positron = electron + gamma1 + gamma2
    graph_gamma_LY    = (TGraph*)file_IHEP_NL_models->Get("gammaScint");
    graph_electronic  = (TGraph*)file_IHEP_NL_models->Get("electronics");
    
    //NL-Gadolinium models:
    
    TFile* unifiedF_NL_models = new TFile("./Inputs/unified_nl_data/nl_models_final.root", "read");
    g_unified_positron_nl = (TGraph*)unifiedF_NL_models->Get("positron_0")->Clone();
    g_unified_positron_nl_pulls[0] =  (TGraph*)unifiedF_NL_models->Get(Form("positron_%d",1))->Clone();
    g_unified_positron_nl_pulls[1] =  (TGraph*)unifiedF_NL_models->Get(Form("positron_%d",2))->Clone();
    g_unified_positron_nl_pulls[2] =  (TGraph*)unifiedF_NL_models->Get(Form("positron_%d",3))->Clone();
    g_unified_positron_nl_pulls[3] =  (TGraph*)unifiedF_NL_models->Get(Form("positron_%d",4))->Clone();
    
    delete unifiedF_NL_models;
    
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

TH2D* nHToyMC :: LoadnHMatrix(Int_t AD,Int_t idx, Int_t idy)
{
    return HighResoMatrixH[AD][idx][idy];
}