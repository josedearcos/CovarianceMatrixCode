/*
Usage: root
       .L Toy.cc+
       Toy()
Input is a flat MC IBD spectrum, then shaped to an input neutrino spectrum via random sampling, then applies scintillator non-linearity, detector non-uniformity, electronics non-linearity, and energy resolution.  
Code written by Xiangpan Ji up to 2014/11/04.  
Edited by Logan Lebanowski up to 2014/11/19.
2015/04/04: Logan added the new NL models and a new resolution function.
2015/05/14: Logan added Tdiff, IBDvolume, EtotScint_N, and EtotScint_Ng[10] and EdepScint_Ng[10] as arrays to handle nGd (multi-gamma) events
*/

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

const double IBDthreshold = 1.80607;  //( (Mn+Me)^2-Mp^2 ) / (2Mp) = 
 // = ( (939.565378+0.510998928)^2 - 938.272046^2 ) / (2*938.272046)
 //1.80433;  // Mneutron+Mpositron-Mproton
const double EnergyScale = 0.982;  // data_centercell/toy_centercell = 0.982 for nH gamma toy non-uniformity and AdSimple data
//const double EnergyScale = 1.023;  // data_centercell/toy_centercell = 1.02 when using same non-uniformity in both data and toy

/// vetex region
int    R2_binnum = 10;                      // ---> set option: divide the volume to sub-regions
double R2_lower  = 0;  // m2                  // ---> set option
double R2_upper  = 4;  // m2                  // ---> set option
double R2_binwidth = 0.4;                   // ---> set option

int    Z_binnum  = 10;                      // ---> set option
double Z_lower   = -2;  // m                  // ---> set option
double Z_upper   = 2;  // m                   // ---> set option

/// find cell
TH2D *hist_findbin;
int global_bin_num = 0;
int local_xbin     = 0;
int local_ybin     = 0;
int local_zbin     = 0;
double usr_r2_P    = 0;
double usr_z_P     = 0;
double usr_r2_Ng   = 0;
double usr_z_Ng    = 0;

/// non-linearity: NL
TGraph *graph_electron_LY;
TGraph *graph_gamma_LY;
TGraph *graph_electronic;

/// non-uniformity: NU
TH2D *hist_map_attenuation;
TH2D *hist_map_pmt_coverage;

TH1D *h_Ev_normal;
TRandom3 *gRandom3;

double func_EnergyResolution(double *x, double *par)
{
  double E = x[0];
  double R = par[0];
  double res = sqrt( 0.004*0.004*E*E + E*(0.082*0.082+R*0.031*0.031) + 0.028*0.028 );  // 7.5/sqrt(E) + 0.9 + 0.865*(R-0.98);
  return res;
}

TF1 *roofunc_EnergyResolution;

// --- ---- --- ---- --- ---- --- ---- --- ---- --- Function declaration

void func_initialization();

int RootCellToVisCell(int RootCell)
{
  // 01 02 03 04 05 06 07 08 09 10
  // 11 12 13 14 15 16 17 18 19 20
  // 21 22 23 24 25 26 27 28 29 30
  // ...
  // 91 92 93 94 95 96 97 98 99 100

  int cell = 0;
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

void Toy()
{
  func_initialization();  

  ////////////

  TFile *roofile_input = new TFile("~/WORK2/logan/TWin/Vertex/ExamTruth/aileron/output/Eprompt.root", "read");  //"./DATA/Eprompt.root"
  TTree *wtree = (TTree*)roofile_input->Get("wtree");
  long entries_wtree = wtree->GetEntries();
  cout<<" ---> input entries: "<<entries_wtree<<endl;

  // define variables to get from input file
  double Ev;
  double cap_Target;
  int IBDvolume;
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

  wtree->SetBranchAddress("cap_Target",   &cap_Target);
  wtree->SetBranchAddress("IBDvolume",    &IBDvolume);
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

  float EupperLim = 20.0;  // MeV
  TH1D *h_Ev_toyMC_input = new TH1D("h_Ev_toyMC_input","h_Ev_toyMC_input", (int)EupperLim*40, 0.0, EupperLim);
  int seed_generator = 1;
  // /*
  TFile *roofile_etree = new TFile("roofile_etree.root", "recreate");
  TTree *etree = new TTree("etree", "effective entries of wtree");
  
  etree->Branch("Ev",           &Ev,           "Ev/D");
  etree->Branch("cap_Target",   &cap_Target,   "cap_Target/D");
  etree->Branch("IBDvolume",    &IBDvolume,    "IBDvolume/I");
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

    cout.precision(3);
    if((ientry+1)%100000==0)
      cout<<" ---> Filtering input: "<<(ientry+1)*100./entries_wtree<<"%"<<endl;

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
    
    //////

    Ev = Etot_P + Etot_N + IBDthreshold;

    h_Ev_toyMC_input->Fill(Ev);
  }
  // */
  /////////

  double max_h_Ev_normal = h_Ev_normal->GetBinContent( h_Ev_normal->GetMaximumBin() );
  h_Ev_normal->Scale(1.0/max_h_Ev_normal);

  int ibin=0, seed_rand=0;
  double prob=0, rand=0;
  // /*
  for(long ientry=0; ientry<entries_wtree; ientry++)
  {
    wtree->GetEntry(ientry);

    cout.precision(3);
    if((ientry+1)%100000==0)
      cout<<" ---> processing01 "<<(ientry+1)*100./entries_wtree<<"%"<<endl;

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

    //////

    Ev = Etot_P + Etot_N + IBDthreshold;
    if( Ev>EupperLim ) {
      cout<<"Determined neutrino energy > "<< EupperLim <<" MeV: "<< Ev <<" at "<< ientry*100./entries_wtree <<"%"<<endl;
      continue;
    }

    ibin = h_Ev_toyMC_input->FindBin( Ev );
    prob = h_Ev_normal->GetBinContent(ibin);
    
    seed_rand = 0;
    if( ientry%100==0 ) 
    {
      seed_generator += 1;
      seed_rand = 2 *seed_generator +1;

      gRandom3->SetSeed(seed_rand);
    }
    
    rand = gRandom3->Uniform(0.0,1.0);
    if( rand>prob ) continue;

    etree->Fill();
  }
  etree->Write();
  roofile_etree->Close();
  // */
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  TFile *roofile_etree_read = new TFile("roofile_etree.root", "read");
  TTree *etree_read = (TTree*)roofile_etree_read->Get("etree");
  long entries_etree_read = etree_read->GetEntries();
  cout<<" ---> filtered entries after neutrino spectrum selection: "<<entries_etree_read<<endl;

  etree_read->SetBranchAddress("Ev",           &Ev);
  etree_read->SetBranchAddress("cap_Target",   &cap_Target);
  etree_read->SetBranchAddress("IBDvolume",    &IBDvolume);
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

  int pcell = 0;
  int dcell = 0;

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

  TFile *roofile_toy = new TFile("roofile_toy.root", "recreate");
  TTree *toy = new TTree("toy", "toyMC result");

  toy->Branch("Ev",           &Ev,           "Ev/D");
  toy->Branch("cap_Target",   &cap_Target,   "cap_Target/D");
  toy->Branch("IBDvolume",    &IBDvolume,    "IBDvolume/I");
  toy->Branch("Tdiff",        &Tdiff,        "Tdiff/D");
  toy->Branch("pcell",        &pcell,        "pcell/I");
  toy->Branch("dcell",        &dcell,        "dcell/I");
  toy->Branch("X_Ng",         &X_Ng,         "X_Ng/D");
  toy->Branch("Y_Ng",         &Y_Ng,         "Y_Ng/D");
  toy->Branch("Z_Ng",         &Z_Ng,         "Z_Ng/D");
  //toy->Branch("EdepScint_Ng", &EdepScint_Ng, "EdepScint_Ng/D");
  //toy->Branch("EtotScint_Ng", &EtotScint_Ng, "EtotScint_Ng/D");
  //toy->Branch("Etot_Ng",      &Etot_Ng,      "Etot_Ng/D");
  //toy->Branch("EdepE_Ng",     &EdepE_Ng,     "EdepE_Ng/D");
  //toy->Branch("EdepEA_Ng",    &EdepEA_Ng,    "EdepEA_Ng/D");
  //toy->Branch("EtotE_Ng",     &EtotE_Ng,     "EtotE_Ng/D");
  //toy->Branch("Etot_N",       &Etot_N,       "Etot_N/D");
  //toy->Branch("EtotScint_N",  &EtotScint_N,  "EtotScint_N/D");
  toy->Branch("X_P",          &X_P,          "X_P/D");
  toy->Branch("Y_P",          &Y_P,          "Y_P/D");
  toy->Branch("Z_P",          &Z_P,          "Z_P/D");
  //toy->Branch("Etot_P",       &Etot_P,       "Etot_P/D");
  //toy->Branch("Etot_Pg1",     &Etot_Pg1,     "Etot_Pg1/D");
  //toy->Branch("Etot_Pg2",     &Etot_Pg2,     "Etot_Pg2/D");
  //toy->Branch("EtotScint_P",  &EtotScint_P,  "EtotScint_P/D");
  //toy->Branch("EtotScint_Pg", &EtotScint_Pg, "EtotScint_Pg/D");
  //toy->Branch("EtotScint_Pg1", &EtotScint_Pg1, "EtotScint_Pg1/D");
  //toy->Branch("EtotScint_Pg2", &EtotScint_Pg2, "EtotScint_Pg2/D");
  //toy->Branch("EdepScint_P",   &EdepScint_P,   "EdepScint_P/D");
  //toy->Branch("EdepAcrylic_P", &EdepAcrylic_P, "EdepAcrylic_P/D");
  //toy->Branch("EdepScint_Pg1", &EdepScint_Pg1, "EdepScint_Pg1/D");
  //toy->Branch("EdepScint_Pg2", &EdepScint_Pg2, "EdepScint_Pg2/D");
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
  toy->Branch("Res_E_P_Sum",   &Res_E_P_Sum,   "Res_E_P_Sum/D");
  toy->Branch("Erec_P",   &Erec_P,   "Erec_P/D");
 
  //toy->Branch("LY_E_Ng",    &LY_E_Ng,    "LY_E_Ng/D"); 
  //toy->Branch("Opt_E_Ng",   &Opt_E_Ng,   "Opt_E_Ng/D"); 
  //toy->Branch("FEE_E_Ng",   &FEE_E_Ng,   "FEE_E_Ng/D"); 
  //toy->Branch("Scale_E_Ng", &Scale_E_Ng, "Scale_E_Ng/D"); 
  toy->Branch("Res_E_Ng",   &Res_E_Ng,   "Res_E_Ng/D"); 
  toy->Branch("Erec_Ng",   &Erec_Ng,   "Erec_Ng/D"); 

  double EendScint_P=0, EendScint_Pg1=0, EendScint_Pg2=0, EendScint_Ng=0;
  //double AdSimple_P=0, AdSimple_Ng=0;
  double usr_opt_attenuation_P=0, usr_pmt_coverage_P=0, usr_opt_attenuation_Ng=0, usr_pmt_coverage_Ng=0;
  double energy_sigma=0, R2_AA=0, R2_BB=0, R_average=0;

  for(long ientry=0; ientry<entries_etree_read; ientry++)
  {
    etree_read->GetEntry(ientry);

    cout.precision(3);
    if((ientry+1)%20000==0)
      cout<<" ---> processing02 "<<(ientry+1)*100./entries_etree_read<<"%"<<endl;

    seed_rand = 0;
    if( ientry%100==0 ) 
    {
      seed_generator += 1;
      seed_rand = 2 *seed_generator +1;

      gRandom3->SetSeed(seed_rand);
    }


    ///////////////////////////// LY
    ///////////////////////////// LY
    /*
    /// case01: most realistic (not recommended now due to issues with beta energies)
    /// prompt

    if( EtotScint_Pg-(EtotE_P-EdepE_P) < -0.00001 ) {
      cout<<"EtotScint_Pg-(EtotE_P-EdepE_P) = "<< EtotScint_Pg-(EtotE_P-EdepE_P) <<" < 0."<<endl;
      continue;
    }
    if( EtotE_Pg1-EdepE_Pg1 < -0.00001 ) {
      cout<<"EtotE_Pg1-EdepE_Pg1 = "<< EtotE_Pg1-EdepE_Pg1 <<" < 0."<<endl;
      continue;
    }
    if( EtotE_Pg2-EdepE_Pg2 < -0.00001 ) {
      cout<<"EtotE_Pg2-EdepE_Pg2 = "<< EtotE_Pg2-EdepE_Pg2 <<" < 0."<<endl;
      continue;
    }
    if( EtotE_Ng-EdepE_Ng < -0.00001 ) {
      cout<<"EtotE_Ng-EdepE_Ng = "<< EtotE_Ng-EdepE_Ng <<" < 0."<<endl;
      continue;
    }

    LY_E_P = EtotScint_P * graph_electron_LY->Eval( EtotScint_P, 0, "s" )
	- (EtotScint_P-EdepScint_P - (EtotScint_Pg-(EtotE_P-EdepE_P))) * graph_electron_LY->Eval( EtotScint_P-EdepScint_P - (EtotScint_Pg-(EtotE_P-EdepE_P)), 0, "s" );

    LY_E_Pg1 = EtotE_Pg1 * graph_gamma_LY->Eval( EtotE_Pg1, 0, "s" )
	- (EtotE_Pg1-EdepE_Pg1) * graph_gamma_LY->Eval( EtotE_Pg1-EdepE_Pg1, 0, "s" );

    LY_E_Pg2 = EtotE_Pg2 * graph_gamma_LY->Eval( EtotE_Pg2, 0, "s" )
	- (EtotE_Pg2-EdepE_Pg2) * graph_gamma_LY->Eval( EtotE_Pg2-EdepE_Pg2, 0, "s" );

    /// delayed
    LY_E_Ng = EtotE_Ng * graph_gamma_LY->Eval( EtotE_Ng, 0, "s" )
	- (EtotE_Ng-EdepE_Ng) * graph_gamma_LY->Eval( EtotE_Ng-EdepE_Ng, 0, "s" );
    */ 
    /// case02: default
    EendScint_P = EtotScint_P-EdepScint_P;//Edep_P;
    EendScint_Pg1 = EtotScint_Pg1-EdepScint_Pg1;//Edep_Pg1;
    EendScint_Pg2 = EtotScint_Pg2-EdepScint_Pg2;//Edep_Pg2;
    //EendScint_Ng = EtotScint_Ng-EdepScint_Ng;//Edep_Ng;

    /// prompt
    LY_E_P = EtotScint_P * graph_electron_LY->Eval( EtotScint_P, 0, "s" );
    if( EendScint_P>0.000001 ) // calculate leakage effect only when significant
       LY_E_P -= EendScint_P * graph_electron_LY->Eval( EendScint_P, 0, "s" );

    LY_E_Pg1 = EtotScint_Pg1 * graph_gamma_LY->Eval( EtotScint_Pg1, 0, "s" );
    if( EendScint_Pg1>0.000001 ) // calculate leakage effect only when significant
      LY_E_Pg1 -= EendScint_Pg1 * graph_gamma_LY->Eval( EendScint_Pg1, 0, "s" );

    LY_E_Pg2 = EtotScint_Pg2 * graph_gamma_LY->Eval( EtotScint_Pg2, 0, "s" );
    if( EendScint_Pg2>0.000001 ) // calculate leakage effect only when significant
      LY_E_Pg2 -= EendScint_Pg2 * graph_gamma_LY->Eval( EendScint_Pg2, 0, "s" );

    /// delayed - loop over gammas
    LY_E_Ng = 0;
    for( int k=0; k<10; k++ ) {
      if( EtotScint_Ng[k]<=0 ) continue;  // most array elements will be zero
      LY_E_Ng += EtotScint_Ng[k] * graph_gamma_LY->Eval( EtotScint_Ng[k], 0, "s" );
      EendScint_Ng = EtotScint_Ng[k]-EdepScint_Ng[k];//Edep_Ng[k];
      if( EendScint_Ng>0.000001 ) // calculate leakage effect only when significant
	 LY_E_Ng -= EendScint_Ng * graph_gamma_LY->Eval( EendScint_Ng, 0, "s" );
      //if( EendScint_Ng<-0.000001 ) {
      //  cout<<" ERROR: EtotScint_Ng["<< k <<"] - EdepScint_Ng["<< k <<"] = "<< EendScint_Ng <<".  Skipped gamma "<< k <<endl;
      //  continue;
      //}
    }


    // ------------------------ //
    /// neutron - determined from NuWa for 0-0.18 MeV (does not consider leakage)
    LY_E_N = EtotScint_N * ( 0.186 + exp(-1.142-126.4*EtotScint_N) + exp(-1.217-17.90*EtotScint_N) ) * graph_electron_LY->Eval( 0.2, 0, "s" );

    /// prompt SUM
    LY_E_P_Sum = LY_E_P + LY_E_Pg1 + LY_E_Pg2 + LY_E_N;
    //------------------------ //


    ///////////////////////////// Optical NU
    ///////////////////////////// Optical NU
    /*
    usr_opt_attenuation_P = 0;
    usr_pmt_coverage_P  = 0;
    usr_opt_attenuation_Ng = 0;
    usr_pmt_coverage_Ng  = 0;
    */
    /// prompt
    usr_z_P  = Z_P * 1e-3;  // mm ---> m
    usr_r2_P = (X_P*X_P+Y_P*Y_P) * 1e-6;  // mm2 ---> m2
    // /*
    global_bin_num = hist_findbin->FindBin(usr_r2_P, usr_z_P);
    hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);

    if( local_zbin != 0) continue;

    usr_opt_attenuation_P = hist_map_attenuation->GetBinContent(local_xbin, local_ybin);
    usr_pmt_coverage_P  = hist_map_pmt_coverage->GetBinContent(local_xbin, local_ybin);

    Opt_E_P_Sum = LY_E_P_Sum * usr_opt_attenuation_P * usr_pmt_coverage_P;  // Tsinghua model
    // */
    //AdSimple_P = ( (7.84628 * (1 + 3.41294e-02*usr_r2_P) * (1 - 1.21750e-02*usr_z_P - 1.64275e-02*usr_z_P*usr_z_P + 7.33006e-04*pow(usr_z_P,3)))/8.05 );  // AdSimple
    // Doc7334(old function), updated from http://dayabay.ihep.ac.cn/tracs/dybsvn/browser/dybgaudi/trunk/Reconstruction/QsumEnergy/src/components/QsumEnergyTool.cc

    //Opt_E_P_Sum = LY_E_P_Sum * AdSimple_P;

    /// delayed
    usr_r2_Ng = (X_Ng*X_Ng+Y_Ng*Y_Ng) * 1e-6;  // mm2 ---> m2
    usr_z_Ng  = Z_Ng * 1e-3;  // mm ---> m
    // /*
    global_bin_num = hist_findbin->FindBin(usr_r2_Ng, usr_z_Ng);
    hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);

    if( local_zbin != 0) continue;
      
    usr_opt_attenuation_Ng = hist_map_attenuation->GetBinContent(local_xbin, local_ybin);
    usr_pmt_coverage_Ng  = hist_map_pmt_coverage->GetBinContent(local_xbin, local_ybin);

    Opt_E_Ng = LY_E_Ng * usr_opt_attenuation_Ng * usr_pmt_coverage_Ng;  // Tsinghua model
    // */
    //AdSimple_Ng = ( (7.84628 * (1 + 3.41294e-02*usr_r2_Ng) * (1 - 1.21750e-02*usr_z_Ng - 1.64275e-02*usr_z_Ng*usr_z_Ng + 7.33006e-04*pow(usr_z_Ng,3)))/8.05 );  // AdSimple
    // Doc7334(old function), updated from http://dayabay.ihep.ac.cn/tracs/dybsvn/browser/dybgaudi/trunk/Reconstruction/QsumEnergy/src/components/QsumEnergyTool.cc

    //Opt_E_Ng = LY_E_Ng * AdSimple_Ng;
    

    ///////////////////////////// FEE
    ///////////////////////////// FEE

    /// prompt
    if( Opt_E_P_Sum!=Opt_E_P_Sum || Opt_E_P_Sum<0 ) {
      cout<<" Opt_E_P_Sum = "<< Opt_E_P_Sum <<endl;
      continue;
    }
    FEE_E_P_Sum = Opt_E_P_Sum * graph_electronic->Eval( Opt_E_P_Sum, 0, "s" );
    if( FEE_E_P_Sum!=FEE_E_P_Sum || FEE_E_P_Sum<0 ) {
      cout<<" FEE_E_P_Sum = "<< FEE_E_P_Sum <<endl;
      continue;
    }
    
    /// delayed
    if( Opt_E_Ng!=Opt_E_Ng || Opt_E_Ng<0 ) {
      cout<<" Opt_E_Ng = "<< Opt_E_Ng <<endl;
      continue;
    }
    FEE_E_Ng = Opt_E_Ng * graph_electronic->Eval( Opt_E_Ng, 0, "s" );
    if( FEE_E_Ng!=FEE_E_Ng || FEE_E_Ng<0 ) {
      cout<<" FEE_E_Ng = "<< FEE_E_Ng <<", Opt_E_Ng = "<< Opt_E_Ng <<endl;
      continue;
    }


    ///////////////////////////// EnergyScale
    ///////////////////////////// EnergyScale

    /// prompt
    Scale_E_P_Sum = FEE_E_P_Sum *EnergyScale;
    if( Scale_E_P_Sum!=Scale_E_P_Sum || Scale_E_P_Sum<0 ) {
      cout<<" Scale_E_P_Sum = "<< Scale_E_P_Sum <<endl;
      continue;
    }

    /// delayed
    Scale_E_Ng = FEE_E_Ng *EnergyScale;
    if( Scale_E_Ng!=Scale_E_Ng || Scale_E_Ng<0 ) {
      cout<<" Scale_E_Ng = "<< Scale_E_Ng <<endl;
      continue;
    }


    ///////////////////////////// Resolution
    ///////////////////////////// Resolution

    energy_sigma = 0;
    R2_AA        = 0;
    R2_BB        = 0;
    R_average    = 0;

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
    if( Res_E_P_Sum < 0 )  // some low-energy events are smeared to below zero 
      Res_E_P_Sum = -Res_E_P_Sum;
    if( Res_E_P_Sum!=Res_E_P_Sum ) {
      cout<<" Res_E_P_Sum = "<< Res_E_P_Sum <<", R_average="<< R_average <<", Scale_E_P_Sum="<< Scale_E_P_Sum <<", energy_sigma="<< energy_sigma <<endl;
      //continue;
    }

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
    if( Res_E_Ng < 0 )  // some low-energy events are smeared to below zero
      Res_E_Ng = -Res_E_Ng;
    if( Res_E_Ng!=Res_E_Ng ) {
      cout<<" Res_E_Ng = "<< Res_E_Ng <<", R_average="<< R_average <<", Scale_E_Ng="<< Scale_E_Ng <<", energy_sigma="<< energy_sigma <<endl;
      continue;
    }


    ///////////////////////////// Reconstructed Energy
    ///////////////////////////// Reconstructed Energy
   /// Estimate an "Erec" on which to apply cuts

    /// prompt
    Erec_P = Res_E_P_Sum / (usr_opt_attenuation_P *usr_pmt_coverage_P);  // Tsinghua model
    //Erec_P = Res_E_P_Sum / AdSimple_P;  // AdSimple

    /// delayed
    Erec_Ng = Res_E_Ng / (usr_opt_attenuation_Ng * usr_pmt_coverage_Ng);  // Tsinghua model
    //Erec_Ng = Res_E_Ng / AdSimple_Ng;  // AdSimple


    ///////////////////////////// Cell
    ///////////////////////////// Cell

    
    toy->Fill();

  }

  toy->Write();
  roofile_toy->Close();

}

// --- ---- --- ---- --- ---- --- ---- --- ---- --- Function definition
// --- ---- --- ---- --- ---- --- ---- --- ---- --- 
// --- ---- --- ---- --- ---- --- ---- --- ---- --- 

void func_initialization()
{
  //////
  hist_findbin = new TH2D("hist_findbin","", R2_binnum,R2_lower,R2_upper, Z_binnum,Z_lower,Z_upper);

  TFile *file_h_Ev = new TFile("./DATA/file_h_Ev.root", "read");
  h_Ev_normal = (TH1D*)file_h_Ev->Get("h_Ev");

  ////// NL
  TFile* file_IHEP_NL_models = new TFile("./DATA/IHEP_NL_models/energyModel_march2015.root", "read");
  graph_electron_LY = (TGraph*)file_IHEP_NL_models->Get("electronScintNL");  // positron = electron + gamma1 + gamma2
  graph_gamma_LY    = (TGraph*)file_IHEP_NL_models->Get("gammaScintNL");
  graph_electronic  = (TGraph*)file_IHEP_NL_models->Get("electronicsNL");  
  
  ////// NU
  TFile *file_hist_attenuation_rel = new TFile("./DATA/uniformity/hist_attenuation_rel.root", "read");
  hist_map_attenuation = (TH2D*)file_hist_attenuation_rel->Get("hist_attenuation_rel");
  hist_map_attenuation->SetName("hist_map_attenuation");

  TFile *file_hist_coverage_rel = new TFile("./DATA/uniformity/hist_coverage_rel.root", "read");
  hist_map_pmt_coverage = (TH2D*)file_hist_coverage_rel->Get("hist_coverage_rel");
  hist_map_pmt_coverage->SetName("hist_map_pmt_coverage");

  ///
  gRandom3 = new TRandom3();
  
  /// energy resolution
  roofunc_EnergyResolution = new TF1("roofunc_EnergyResolution", func_EnergyResolution, 0, 20, 1);
}
