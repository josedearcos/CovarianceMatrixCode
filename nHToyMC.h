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

// --- ---- --- ---- --- ---- --- ---- --- ---- --- Global

#define FullEnergyResolution
//#define UseDelayInformation
//#define SaveTree // To save toy tree
//#define LoadTree // To load eprompt tree
///
TString roostr;

///
const Double_t DeltaM = 1.80433; // Mneutron+Mpositron-Mproton
const Double_t EnergyScale = 0.982; // data_centercell/toy_centercell = 0.982 for nH gamma

class nHToyMC
{
private:
    
    //Binning parameters:
    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    bool LinearBinning;
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.
    
    Int_t n_evis_bins;
    Int_t n_etrue_bins;
    
    //Systematic parameters;
    bool RelativeEnergyScaleMatrix;
    bool IAVMatrix;
    bool NLMatrix;
    bool ResolutionMatrix;
    bool EfficiencyMatrix;
    
    //AD configuration parameters:
    Int_t NADs;
    
    /// vetex region
    Int_t    R2_binnum  ;                      // ---> set option: divide the volume to sub-regions
    Double_t R2_lower   ;// m2                  // ---> set option
    Double_t R2_upper   ;// m2                  // ---> set option
    Double_t R2_binwidth;                   // ---> set option
    
    Int_t    Z_binnum  ;                      // ---> set option
    Double_t Z_lower   ;// m                  // ---> set option
    Double_t Z_upper   ;// m                   // ---> set option
    
    /// find cell
    TH2D *hist_findbin;
    Int_t global_bin_num;
    Int_t local_xbin    ;
    Int_t local_ybin    ;
    Int_t local_zbin    ;
    Double_t usr_r2     ;
    Double_t usr_z      ;
    
    /// non-linearity: NL
    TGraph *graph_electron_LY;
    TGraph *graph_gamma_LY;
    TGraph *graph_electronic;
    
    /// non-uniformity: NU
    TH2D *hist_map_attenuation;
    TH2D *hist_map_pmt_coverage;
    
    ///
    TH1D *h_Ev_normal;
    TRandom3 *gRandom3;
    TRandom3 *RandomSys;

    ///
    TF1 *roofunc_EnergyResolution;
    
    Double_t func_EnergyResolution(double*,double*);
    void func_initialization();
    Int_t  RootCellToVisCell(Int_t RootCell);
    TH1D* PredictionH[MaxDetectors];
    TH1D* VisiblePredictionH[MaxDetectors];
    TH2D* TransMatrixH[MaxDetectors];
    TH2D* HighResoTransMatrixH[MaxDetectors];
    TH2D* MatrixH[MaxDetectors];
public:
    nHToyMC(NominalData*);//constructor
    ~nHToyMC();//destructor
    void Toy(bool);//main
    TH2D* LoadnHMatrix(Int_t);
    void GeneratenHResponseMatrix();
    TH2D* HighResoMatrixH[MaxDetectors];
};

// constructor:
nHToyMC :: nHToyMC(NominalData* Data)
{
    R2_binnum = 10;                      // ---> set option: divide the volume to sub-regions
    R2_lower  = 0;// m2                  // ---> set option
    R2_upper  = 4;// m2                  // ---> set option
    R2_binwidth = 0.4;                   // ---> set option
    
    Z_binnum  = 10;                      // ---> set option
    Z_lower   = -2;// m                  // ---> set option
    Z_upper   = 2;// m                   // ---> set option
    
    /// find cell
    global_bin_num = 0;
    local_xbin     = 0;
    local_ybin     = 0;
    local_zbin     = 0;
    usr_r2      = 0;
    usr_z       = 0;
    
    
    InitialEnergy = Data->GetEmin();
    FinalEnergy = Data->GetEmax();
    InitialVisibleEnergy = Data->GetEVisMin();
    FinalVisibleEnergy = Data->GetEVisMax();
    
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
    
    RelativeEnergyScaleMatrix = Data->GetRelativeEnergyScaleMatrix();
    //    RelativeEnergyOffsetMatrix = Data->GetRelativeEnergyOffsetMatrix();
    //    AbsoluteEnergyScaleMatrix = Data->GetAbsoluteEnergyScaleMatrix();
    //    AbsoluteEnergyOffsetMatrix = Data->GetAbsoluteEnergyOffsetMatrix();
    IAVMatrix = Data->GetIAVMatrix();
    NLMatrix = Data->GetNLMatrix();
    ResolutionMatrix = Data->GetResolutionMatrix();
    EfficiencyMatrix = Data->GetEfficiencyMatrix();
    
    NADs = Data->GetADs();
    for(Int_t AD = 0; AD < NADs; AD++)
    {
        HighResoMatrixH[AD] = new TH2D(Form("Fine_nHResponseMatrix_AD%d",AD+1),Form("Fine_nHResponseMatrix_AD%d",AD+1),MatrixBins,0,FinalVisibleEnergy,MatrixBins,0,FinalVisibleEnergy);//from true to visible
      
        MatrixH[AD] = new TH2D(Form("nHResponseMatrix_AD%d",AD+1),Form("nHResponseMatrix_AD%d",AD+1),n_etrue_bins,enu_bins,n_evis_bins,evis_bins);//from true to visible
        
        TransMatrixH[AD] = new TH2D(Form("nHTransResponseMatrix_AD%d",AD+1),Form("nHTransResponseMatrix_AD%d",AD+1),n_evis_bins,evis_bins,n_etrue_bins,enu_bins);//from visible to true
        
        HighResoTransMatrixH[AD] = new TH2D(Form("Fine_nHTransResponseMatrix_AD%d",AD+1),Form("Fine_nHTransResponseMatrix_AD%d",AD+1),MatrixBins,0,FinalVisibleEnergy,MatrixBins,0,FinalVisibleEnergy);//from visible to true
        
        PredictionH[AD] = new TH1D(Form("PredAD%d",AD),Form("PredAD%d",AD), Nbins, InitialEnergy, FinalEnergy);
        
        VisiblePredictionH[AD] = new TH1D(Form("VisPredAD%d",AD),Form("VisPredAD%d",AD), MatrixBins, InitialVisibleEnergy, FinalVisibleEnergy);
    }
}

nHToyMC :: ~nHToyMC()
{
    for(Int_t AD = 0; AD < NADs; AD++)
    {
        delete PredictionH[AD];
        delete HighResoMatrixH[AD];
        delete MatrixH[AD];
        delete TransMatrixH[AD];
        delete VisiblePredictionH[AD];
    }
}
///
Double_t nHToyMC :: func_EnergyResolution(Double_t *x, Double_t *par)// from Logan
{
    
#ifdef FullEnergyResolution
    Double_t E = x[0];
    Double_t R = par[0];
    Double_t res = 7.5/sqrt(E) + 0.9 + 0.865*(R-0.98);
    res = E*res*0.01;
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
    hist_findbin->GetBinXYZ(RootCell, local_xbin, local_ybin, local_zbin);
    
    if( local_zbin==0 && local_xbin>=1 && local_xbin<=10 && local_ybin>=1 && local_ybin<=10 )
    {
        cell = R2_binnum*Z_binnum - local_ybin*R2_binnum + local_xbin;
    }
    else
    {
        cell = 200;
    }
    
    return cell;
}

// --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- ---
// --- ---- --- ---- --- ---- --- ---- ---      --- ---- --- ---- --- ---- --- ---- --- ---- ---
// --- ---- --- ---- --- ---- --- ---- --- MAIN --- ---- --- ---- --- ---- --- ---- --- ---- ---
// --- ---- --- ---- --- ---- --- ---- ---      --- ---- --- ---- --- ---- --- ---- --- ---- ---
// --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- --- ---- ---

void nHToyMC :: Toy(bool mode)
{
    
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
        if(ientry%10000==0)
            cout<<" ---> processing00 "<<ientry*100./entries_wtree<<"%"<<endl;
        
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
    }
    
    /////////
    
    Double_t max_h_Ev_toyMC_input = h_Ev_toyMC_input->GetBinContent( h_Ev_toyMC_input->GetMaximumBin() );
    Double_t max_h_Ev_normal = h_Ev_normal->GetBinContent( h_Ev_normal->GetMaximumBin() );
    h_Ev_normal->Scale(max_h_Ev_toyMC_input/max_h_Ev_normal);
    
    for(long ientry=0; ientry<entries_wtree; ientry++)
    {
        wtree->GetEntry(ientry);
        
        cout.precision(3);
        if(ientry%10000==0)
            cout<<" ---> processing01 "<<ientry*100./entries_wtree<<"%"<<endl;
        
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
    
    Double_t LY_E_P_Sum;
    Double_t Opt_E_P_Sum;
    Double_t FEE_E_P_Sum;
    Double_t Scale_E_P_Sum;
    Double_t Res_E_P_Sum;
    
#ifdef UseDelayInformation
    Double_t LY_Etot_Ng;
    Double_t Opt_Etot_Ng;
    Double_t FEE_Etot_Ng;
    Double_t Scale_Etot_Ng;
    Double_t Res_Etot_Ng;
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
    toy->Branch("LY_E_P_Sum",    &LY_E_P_Sum,    "LY_E_P_Sum/D");
    toy->Branch("Opt_E_P_Sum",   &Opt_E_P_Sum,   "Opt_E_P_Sum/D");
    toy->Branch("FEE_E_P_Sum",   &FEE_E_P_Sum,   "FEE_E_P_Sum/D");
    toy->Branch("Scale_E_P_Sum", &Scale_E_P_Sum, "Scale_E_P_Sum/D");
    toy->Branch("Res_E_P_Sum",   &Res_E_P_Sum,   "Res_E_P_Sum/D");
    
    toy->Branch("LY_Etot_Ng",    &LY_Etot_Ng,    "LY_Etot_Ng/D");
    toy->Branch("Opt_Etot_Ng",   &Opt_Etot_Ng,   "Opt_Etot_Ng/D");
    toy->Branch("FEE_Etot_Ng",   &FEE_Etot_Ng,   "FEE_Etot_Ng/D");
    toy->Branch("Scale_Etot_Ng", &Scale_Etot_Ng, "Scale_Etot_Ng/D");
    toy->Branch("Res_Etot_Ng",   &Res_Etot_Ng,   "Res_Etot_Ng/D");
#endif
    
    Double_t GausRelative[MaxDetectors],GausIAV[MaxDetectors],GausNL[MaxDetectors],GausReso[MaxDetectors],GausEff[MaxDetectors];
    
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
    
    //Use the loaded tree data in each AD in an independent way:
    for(Int_t AD = 0; AD<NADs;AD++)
    {
        seed_generator = 1;//this way all matrices will be random but have a common nominal spectrum
        
        //Process tree
        for(long ientry=0; ientry<entries_etree_read; ientry++)
        {
            etree_read->GetEntry(ientry);

            cout.precision(3);
            if(ientry%10000==0)
                cout<<" ---> processing02 "<<ientry*100./entries_etree_read<<"%"<<endl;
            
            Int_t seed_rand = 0;
            if( ientry%100==0 )
            {
                seed_generator += 1;
                seed_rand = 2 *seed_generator +1;
                
                gRandom3->SetSeed(seed_rand);
            }
            
            //each event has a different systematic error:
            
            if(mode==0)//Nominal values
            {
                GausRelative[AD] = 0;
                GausIAV[AD]= 0;
                GausNL[AD]= 0;
                GausReso[AD]=0;
                GausEff[AD]= 0;
            }
            else//vary systematic for each entry.
            {
                GausRelative[AD] = RandomSys->Gaus(0,1);
                GausIAV[AD]= RandomSys->Gaus(0,1);
                GausNL[AD]= RandomSys->Gaus(0,1);
                GausReso[AD]= RandomSys->Gaus(0,1);
                GausEff[AD]= RandomSys->Gaus(0,1);
                
                //factor that is multiplied by the random error so it's added to the nominal by using:  Value = NominalValue + gRandom3->Gauss(0,1)*1sigmaError;
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
             LY_Etot_Ng = EtotScint_Ng * graph_gamma_LY->Eval( EtotScint_Ng, 0, "s" )
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
             LY_Etot_Ng = EtotScint_Ng * graph_gamma_LY->Eval( EtotScint_Ng, 0, "s" );
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
             LY_Etot_Ng = EtotScint_Ng * graph_gamma_LY->Eval( EtotScint_Ng, 0, "s" )
             - (EtotE_Ng-EdepE_Ng) * graph_electron_LY->Eval( EtotE_Ng-EdepE_Ng, 0, "s" );
             */
            
            /// case04: default ? edit
            Double_t P_Scint = EtotScint_P-Edep_P;
            Double_t P_totElectron = EtotE_P-EdepE_P;
            
            Double_t Pg1_totElectron = EtotScint_Pg1-EdepE_Pg1;
            Double_t Pg2_totElectron = EtotScint_Pg2-EdepE_Pg2;
            
            Double_t Ng_totElectron = EtotScint_Ng-EdepE_Ng;
            
            if( P_Scint<0 )         P_Scint = 0;
            if( Pg1_totElectron<0 ) continue;
            if( Pg2_totElectron<0 ) continue;
            if( Ng_totElectron<0 )  continue;
            
            /// prompt
            LY_E_P = EtotScint_P * graph_electron_LY->Eval( EtotScint_P, 0, "s" )
            - (P_Scint) * graph_electron_LY->Eval( P_Scint, 0, "s" )
            - (P_totElectron) * graph_gamma_LY->Eval( P_totElectron, 0, "s" );
            
            LY_E_Pg1 = EtotScint_Pg1 * graph_gamma_LY->Eval( EtotScint_Pg1, 0, "s" )
            - (Pg1_totElectron) * graph_gamma_LY->Eval( Pg1_totElectron, 0, "s" );
            
            LY_E_Pg2 = EtotScint_Pg2 * graph_gamma_LY->Eval( EtotScint_Pg2, 0, "s" )
            - (Pg2_totElectron) * graph_gamma_LY->Eval( Pg2_totElectron, 0, "s" );
            
            LY_E_P_Sum = LY_E_P + LY_E_Pg1 + LY_E_Pg2;
            
            /// delayed
#ifdef UseDelayInformation
            LY_Etot_Ng = EtotScint_Ng * graph_gamma_LY->Eval( EtotScint_Ng, 0, "s" )
            - (Ng_totElectron) * graph_gamma_LY->Eval( Ng_totElectron, 0, "s" );
#endif
            ///////////////////////////// Optical NU
            ///////////////////////////// Optical NU
            
            Double_t usr_opt_attention = 0;
            Double_t usr_pmt_coverage  = 0;
            
            /// prompt
            Double_t VX_P = X_P;
            Double_t VY_P = Y_P;
            Double_t VZ_P = Z_P;
            
            usr_r2 = (VX_P*VX_P+VY_P*VY_P) * 1e-6;// mm2 ---> m2
            usr_z  = VZ_P * 1e-3;// mm ---> m
            
            global_bin_num = hist_findbin->FindBin(usr_r2, usr_z);
            hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
            
            if( local_zbin != 0) continue;
            
            usr_opt_attention = hist_map_attenuation->GetBinContent(local_xbin, local_ybin);
            usr_pmt_coverage  = hist_map_pmt_coverage->GetBinContent(local_xbin, local_ybin);
            
            Opt_E_P_Sum = LY_E_P_Sum *usr_opt_attention *usr_pmt_coverage;
            
            /// delayed
#ifdef UseDelayInformation
            
            Double_t VX_Ng = X_Ng;
            Double_t VY_Ng = Y_Ng;
            Double_t VZ_Ng = Z_Ng;
            
            usr_r2 = (VX_Ng*VX_Ng+VY_Ng*VY_Ng) * 1e-6;// mm2 ---> m2
            usr_z  = VZ_Ng * 1e-3;// mm ---> m
            
            global_bin_num = hist_findbin->FindBin(usr_r2, usr_z);
            hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
            
            if( local_zbin != 0) continue;
            
            usr_opt_attention = hist_map_attenuation->GetBinContent(local_xbin, local_ybin);
            usr_pmt_coverage  = hist_map_pmt_coverage->GetBinContent(local_xbin, local_ybin);
            
            Opt_Etot_Ng = LY_Etot_Ng *usr_opt_attention *usr_pmt_coverage;
            
#endif
            ///////////////////////////// FEE
            ///////////////////////////// FEE
            
            /// prompt
            FEE_E_P_Sum = Opt_E_P_Sum *graph_electronic->Eval( Opt_E_P_Sum, 0, "s" );
            
            /// delayed
#ifdef UseDelayInformation
            FEE_Etot_Ng = Opt_Etot_Ng *graph_electronic->Eval( Opt_Etot_Ng, 0, "s" );
#endif
            ///////////////////////////// EnergyScale
            ///////////////////////////// EnergyScale
            
            /// prompt
            Scale_E_P_Sum = FEE_E_P_Sum *EnergyScale;
            
            /// delayed
#ifdef UseDelayInformation
            Scale_Etot_Ng = FEE_Etot_Ng *EnergyScale;
#endif
            ///////////////////////////// Resolution
            ///////////////////////////// Resolution
            
            Double_t energy_sigma = 0;
            Double_t R2_AA        = 0;
            Double_t R2_BB        = 0;
            Double_t R_average    = 0;
            
            /// prompt
            usr_r2 = (VX_P*VX_P+VY_P*VY_P) * 1e-6;// mm2 ---> m2
            usr_z  = VZ_P * 1e-3;// mm ---> m
            
            global_bin_num = hist_findbin->FindBin(usr_r2, usr_z);
            hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
            
            if( local_zbin != 0) continue;
            
            R2_AA = (local_xbin-1) * R2_binwidth;
            R2_BB = local_xbin * R2_binwidth;
            R_average = sqrt( (R2_AA+R2_BB)/2. );// average radius from volume weighted
            
            roofunc_EnergyResolution->SetParameter( 0, R_average );
            energy_sigma = roofunc_EnergyResolution->Eval( Scale_E_P_Sum );
            Res_E_P_Sum = gRandom3->Gaus( Scale_E_P_Sum, energy_sigma );
            
            /// delayed
#ifdef UseDelayInformation
            usr_r2 = (VX_Ng*VX_Ng+VY_Ng*VY_Ng) * 1e-6;// mm2 ---> m2
            usr_z  = VZ_Ng * 1e-3;// mm ---> m
            
            global_bin_num = hist_findbin->FindBin(usr_r2, usr_z);
            hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
            
            if( local_zbin != 0) continue;
            
            R2_AA = (local_xbin-1) * R2_binwidth;
            R2_BB = local_xbin * R2_binwidth;
            R_average = sqrt( (R2_AA+R2_BB)/2. );// average radius from volume weighted
            
            roofunc_EnergyResolution->SetParameter( 0, R_average );
            energy_sigma = roofunc_EnergyResolution->Eval( Scale_Etot_Ng );
            Res_Etot_Ng = gRandom3->Gaus( Scale_Etot_Ng, energy_sigma );
#endif
            ///////////////////////////// Cell
            ///////////////////////////// Cell
            
            //Fill histograms:
            PredictionH[AD]->Fill(cap_Ev);
            VisiblePredictionH[AD]->Fill(Res_E_P_Sum);
            HighResoMatrixH[AD]->Fill(cap_Ev,Res_E_P_Sum);//Fine grid
            MatrixH[AD]->Fill(cap_Ev,Res_E_P_Sum);//neutrino energy vs visible energy
            TransMatrixH[AD]->Fill(Res_E_P_Sum,cap_Ev);
            HighResoTransMatrixH[AD]->Fill(Res_E_P_Sum,cap_Ev);//Fine grid

#ifdef SaveTree
            toy->Fill();
#endif
        }//events
    }//ADs
    
    
    Double_t Norma[n_etrue_bins];//true->vis
    Double_t HighResoNorma[MatrixBins];//true->vis
    Double_t NormaTrans[n_evis_bins];//vis->true
    Double_t HighResoNormaTrans[MatrixBins];//vis->true
    //Normalize the Matrix
    for(Int_t AD = 0; AD<NADs;AD++)
    {
        for(Int_t i=0;i<n_etrue_bins;i++)
        {
            Norma[i]=0;
            for(Int_t j=0;j<n_evis_bins;j++)
            {
                Norma[i] = Norma[i]+MatrixH[AD]->GetBinContent(i+1,j+1);// true->vis
            }
        }
        
        for(Int_t i=0;i<n_evis_bins;i++)
        {
            NormaTrans[i]=0;

            for(Int_t j=0;j<n_etrue_bins;j++)
            {
                NormaTrans[i] = NormaTrans[i]+TransMatrixH[AD]->GetBinContent(i+1,j+1);// vis->true
            }
        }

        for (Int_t i = 0; i < n_etrue_bins; i++)
        {
            for (Int_t j = 0; j < n_evis_bins; j++)
            {
                if(Norma[i]!=0)
                {
                    MatrixH[AD]->SetBinContent(i+1,j+1,MatrixH[AD]->GetBinContent(i+1,j+1)/Norma[i]);//true->vis
                }
            }
        }
        
        for (Int_t i = 0; i < n_evis_bins; i++)
        {
            for (Int_t j = 0; j < n_etrue_bins; j++)
            {
                if(NormaTrans[i]!=0)
                {
                    TransMatrixH[AD]->SetBinContent(i+1,j+1,TransMatrixH[AD]->GetBinContent(i+1,j+1)/NormaTrans[i]);
                }

            }
        }
        
        for(Int_t i=0;i<MatrixBins;i++)
        {
            HighResoNorma[i]=0;
            HighResoNormaTrans[i]=0;
            
            for(Int_t j=0;j<MatrixBins;j++)
            {
                HighResoNorma[i] = HighResoNorma[i]+HighResoMatrixH[AD]->GetBinContent(i+1,j+1);// true->vis
                HighResoNormaTrans[i] = HighResoNormaTrans[i]+HighResoTransMatrixH[AD]->GetBinContent(i+1,j+1);// vis->true
                
            }
        }
        
        for(Int_t i=0;i<MatrixBins;i++)
        {
            for(Int_t j=0;j<MatrixBins;j++)
            {
                if(HighResoNorma[i]!=0)
                {
                    HighResoMatrixH[AD]->SetBinContent(i+1,j+1,HighResoMatrixH[AD]->GetBinContent(i+1,j+1)/HighResoNorma[i]);//true->vis
                }
                
                if(HighResoNormaTrans[i]!=0)
                {
                    HighResoTransMatrixH[AD]->SetBinContent(i+1,j+1,HighResoTransMatrixH[AD]->GetBinContent(i+1,j+1)/HighResoNormaTrans[i]);
                }
            }
        }
    }
    
    if(Print)
    {
        TCanvas* FineMatrixC = new TCanvas("F","F");
        FineMatrixC->Divide(NADs/2,2);
        
        TCanvas* MatrixC = new TCanvas("","");
        MatrixC->Divide(NADs/2,2);
        
        TCanvas* PredictionC = new TCanvas("P","P");
        PredictionC->Divide(NADs/2,2);
        
        TCanvas* VisPredictionC = new TCanvas("VP","VP");
        VisPredictionC->Divide(NADs/2,2);
        
        TCanvas* TransC = new TCanvas("T","T");
        TransC->Divide(NADs/2,2);
        
        TCanvas* FineTransC = new TCanvas("FT","FT");
        FineTransC->Divide(NADs/2,2);
        
        for(Int_t i = 0; i<NADs;i++)
        {
            PredictionC->cd(i+1);
            PredictionH[i]->Draw();
            
            VisPredictionC->cd(i+1);
            VisiblePredictionH[i]->Draw();
            
            MatrixC->cd(i+1);
            MatrixH[i]->Draw("colz");
            
            FineMatrixC->cd(i+1);
            HighResoMatrixH[i]->Draw("colz");
  
            TransC->cd(i+1);
            TransMatrixH[i]->Draw("colz");
            
            FineTransC->cd(i+1);
            HighResoTransMatrixH[i]->Draw("colz");

        }
        
        PredictionC->Print("./Images/TrueHydrogenPrediction.eps");
        VisPredictionC->Print("./Images/VisibleHydrogenPrediction.eps");
        FineMatrixC->Print("./Images/FineHydrogenResponseMatrix.eps");
        MatrixC->Print("./Images/HydrogenResponseMatrix.eps");
        TransC->Print("./Images/TransposeHydrogenMatrix.eps");
        FineTransC->Print("./Images/FineTransposeHydrogenMatrix.eps");

        delete PredictionC;
        delete VisPredictionC;

        delete MatrixC;
        delete FineMatrixC;
        delete TransC;
        delete FineTransC;
    }
    
    if(mode==0)
    {
        TFile* SaveMatrix = new TFile("./ResponseMatrices/Hydrogen/NominalResponseMatrix.root","recreate");
        
        for(Int_t AD = 0; AD<NADs; AD++)//Save different AD matrices to produce the covariance matrices.
        {
            HighResoMatrixH[AD]->Write(Form("FineEvisEnu%d",AD+1));//Save the matrices.
            MatrixH[AD]->Write(Form("EvisEnu%d",AD+1));//Save the matrices.
        }
        
        TransMatrixH[0]->Write("EnuEvis");//All ADs are "identical", only need 1

        delete SaveMatrix;
        
    }
    
#ifdef SaveTree
    toy->Write();
    roofile_toy->Close();
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
    
    ////// ~/WORK/jixp/hapy/usr_job/work_mc_Ep/Production/process/plot_obj.cc
    TFile *file_h_Ev = new TFile("./Inputs/HInputs/Data/file_h_Ev.root", "read");
    h_Ev_normal = (TH1D*)file_h_Ev->Get("h_Ev");
    
    ////// NL
    TFile* file_IHEP_NL_models = new TFile("./Inputs/HInputs/Data/IHEP_NL_models/Model1.root", "read");
    graph_electron_LY = (TGraph*)file_IHEP_NL_models->Get("electronScint");// positron = electron + gamma1 + gamma2
    graph_gamma_LY    = (TGraph*)file_IHEP_NL_models->Get("gammaScint");
    graph_electronic  = (TGraph*)file_IHEP_NL_models->Get("electronics");
    
    ////// NU
    TFile *file_hist_attenuation_rel = new TFile("./Inputs/HInputs/Data/uniformity/hist_attenuation_rel.root", "read");
    hist_map_attenuation = (TH2D*)file_hist_attenuation_rel->Get("hist_attenuation_rel");
    hist_map_attenuation->SetName("hist_map_attenuation");
    
    TFile *file_hist_coverage_rel = new TFile("./Inputs/HInputs/Data/uniformity/hist_coverage_rel.root", "read");
    hist_map_pmt_coverage = (TH2D*)file_hist_coverage_rel->Get("hist_coverage_rel");
    hist_map_pmt_coverage->SetName("hist_map_pmt_coverage");
    
    ///
    gRandom3 = new TRandom3();
    RandomSys = new TRandom3(0);
    /// energy resolution
    roofunc_EnergyResolution = new TF1("roofunc_EnergyResolution", this,&nHToyMC::func_EnergyResolution, 0,20, 1,"nHToyMC","func_EnergyResolution");
}

TH2D* nHToyMC :: LoadnHMatrix(Int_t AD)
{
    return HighResoMatrixH[AD];
}