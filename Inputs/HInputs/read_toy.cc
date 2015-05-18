/*
Usage: root -b read_toy.cc
Uses output from Toy.cc and applies cuts and efficiency map to obtain final spectra (also saves plots by cell)
Code written by Xiangpan Ji up to 2014/11/04.  
2014/11/15: Added analysis cuts and Erec.  Logan Lebanowski
2015/03/22: Added variation of OAV dimensions (does not preserve event number --> invalid at center & edges) .  Logan Lebanowski
2015/05/14: Enabled Tdiff cuts and removed rejection of nGd events.  Logan Lebanowski
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

#include<map>

// --- ---- --- ---- --- ---- --- ---- --- ---- --- Global

///
TString roostr;

/// select whether nGd or nH IBD analysis - nGd events are currently rejected in the code.  Time cuts are not currently applied.
int analysis = 0;  // 0==special cuts, 1==hydrogen, 64==gadolinium

const double Rshield = 2259.15;  // position of radial shield [mm]
const double ZreflectorTop = 2104.5;  // position of top reflector ESR [mm]
const double ZreflectorBottom = -2027.45;  // position of bottom reflector ESR [mm]
const double deltaRoav = 8.0;  // +-8.0 mm variation in diameter within ADs 1 & 2
const double deltaZoav = 3.0;  // +-3.0 mm discrepancy between specification and measurement of ADs 1 & 2

/// vetex region
const int R2_binnum = 10;              // ---> set option: divide the volume to sub-regions
double R2_lower  = 0;  // m2                // ---> set option
double R2_upper  = 4;  // m2                // ---> set option
double R2_binwidth = 4.0/R2_binnum;         // ---> set option

const int Z_binnum = 10;                       // ---> set option
double Z_lower   = -2;  // m                // ---> set option
double Z_upper   = 2;  // m                 // ---> set option

/// find cell
int global_bin_num = 0;
int local_xbin     = 0;
int local_ybin     = 0;
int local_zbin     = 0;  
double usr_r2      = 0;
double usr_z       = 0;

// --- ---- --- ---- --- ---- --- ---- --- ---- --- Function declaration

int RootCellToVisCell(int RootCell)
{
  // 01 02 03 04 05 06 07 08 09 10
  // 11 12 13 14 15 16 17 18 19 20
  // 21 22 23 24 25 26 27 28 29 30
  // ...
  // 91 92 93 94 95 96 97 98 99 100

  int cell = 0;
  
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

void read_toy()
{
  //gROOT->ProcessLine(".x ./lhcbStyle.C");  

  /// Time cuts are not currently applied.
  if(analysis==0) {  // special - for comparison with NuWa MC
    const double EpromptCutLow=0.0, EpromptCutHigh=12.0, EdelayCutLow=0.0, EdelayCutHigh=12.0; // 1.5, 12.0, 1.5, 3.3
    const double TcapCutLow=0.0, TcapCutHigh=0.0, distCut=500.0;
    const double radialCut=0.0;  // adjust OAV radius^2 [m^2]: nominal=3.964 (MC), 3.935 (data) - unc. 5.1 mm (data).
  } else if(analysis==1) {  // nH
    const double EpromptCutLow=1.5, EpromptCutHigh=12.0, EdelayCutLow=1.7955, EdelayCutHigh=2.6535; // update delayed cut values
    const double TcapCutLow=1000.0, TcapCutHigh=400000.0, distCut=500.0;
    const double radialCut=0.0;  // adjust OAV radius^2 [m^2]: nominal=4.0
  } else if(analysis==64) {  // nGd
    const double EpromptCutLow=0.7, EpromptCutHigh=12.0, EdelayCutLow=6.0, EdelayCutHigh=12.0;
    const double TcapCutLow=1000.0, TcapCutHigh=200000.0, distCut=0.0;
    const double radialCut=0.0;  // adjust OAV radius^2 [m^2]: nominal=4.0
  } else {
    printf("\n\nANALYSIS TYPE IS INVALID.  END.\n\n");
    return;
  }

  //double Ev = 0;
  double cap_Target;
  int IBDvolume;  // -1 undefined, 0 other, 1 GdLS, 2 LS, 3 MO, 4 IAV, 5 OAV
  double Tdiff;
  double X_Ng, Y_Ng, Z_Ng;
  //double Edep_Ng;
  //double EtotScint_Ng;
  //double Etot_Ng;
  ///double EdepE_Ng;
  ///double EdepEA_Ng;
  ///double EtotE_Ng;
  //double Etot_N;

  double X_P, Y_P, Z_P;
  //double Etot_P;
  //double Etot_Pg1;
  //double Etot_Pg2;
  //double EtotScint_P;
  ///double EtotScint_Pg;
  //double EtotScint_Pg1;
  //double EtotScint_Pg2;
  //double Edep_P;
  //double EdepAcrylic_P;
  //double Edep_Pg1;
  //double Edep_Pg2;
  ///double EdepE_P;
  ///double EdepEA_P;
  ///double EtotE_P;
  ///double EdepE_Pg1;
  ///double EdepEA_Pg1;
  ///double EtotE_Pg1;
  ///double EdepE_Pg2;
  ///double EdepEA_Pg2;
  ///double EtotE_Pg2;

  int pcell = 0;
  int dcell = 0;

  //double LY_E_Ng;
  //double Opt_E_Ng;
  //double FEE_E_Ng;
  //double Scale_E_Ng;
  double Res_E_Ng;  
  double Erec_Ng;

  //double LY_E_P;
  //double LY_E_Pg1;
  //double LY_E_Pg2;
  //double LY_E_N;
  //double LY_E_P_Sum;
  //double Opt_E_P_Sum;
  //double FEE_E_P_Sum;
  //double Scale_E_P_Sum;
  double Res_E_P_Sum;
  double Erec_P;


  TFile *toyFile = new TFile("roofile_toy.root", "read");
  TTree *toyTree = (TTree*)toyFile->Get("toy");

  //toyTree->SetBranchAddress("Ev",           &Ev);
  toyTree->SetBranchAddress("cap_Target",   &cap_Target);
  toyTree->SetBranchAddress("IBDvolume",    &IBDvolume);
  toyTree->SetBranchAddress("Tdiff",        &Tdiff);
  toyTree->SetBranchAddress("X_Ng",         &X_Ng);
  toyTree->SetBranchAddress("Y_Ng",         &Y_Ng);
  toyTree->SetBranchAddress("Z_Ng",         &Z_Ng);
  //toyTree->SetBranchAddress("Edep_Ng",      &Edep_Ng);
  //toyTree->SetBranchAddress("EtotScint_Ng", &EtotScint_Ng);
  //toyTree->SetBranchAddress("Etot_Ng",      &Etot_Ng);
  ///toyTree->SetBranchAddress("EdepE_Ng",     &EdepE_Ng);
  ///toyTree->SetBranchAddress("EdepEA_Ng",    &EdepEA_Ng);
  ///toyTree->SetBranchAddress("EtotE_Ng",     &EtotE_Ng);
  //toyTree->SetBranchAddress("Etot_N",       &Etot_N);
  toyTree->SetBranchAddress("X_P",          &X_P);
  toyTree->SetBranchAddress("Y_P",          &Y_P);
  toyTree->SetBranchAddress("Z_P",          &Z_P);
  //toyTree->SetBranchAddress("Etot_P",       &Etot_P);
  //toyTree->SetBranchAddress("Etot_Pg1",     &Etot_Pg1);
  //toyTree->SetBranchAddress("Etot_Pg2",     &Etot_Pg2);
  //toyTree->SetBranchAddress("EtotScint_P",  &EtotScint_P);
  ///toyTree->SetBranchAddress("EtotScint_Pg", &EtotScint_Pg);
  //toyTree->SetBranchAddress("EtotScint_Pg1", &EtotScint_Pg1);
  //toyTree->SetBranchAddress("EtotScint_Pg2", &EtotScint_Pg2);
  //toyTree->SetBranchAddress("Edep_P",        &Edep_P);
  //toyTree->SetBranchAddress("EdepAcrylic_P", &EdepAcrylic_P);
  //toyTree->SetBranchAddress("Edep_Pg1",      &Edep_Pg1);
  //toyTree->SetBranchAddress("Edep_Pg2",      &Edep_Pg2);
  ///toyTree->SetBranchAddress("EdepE_P",       &EdepE_P);
  ///toyTree->SetBranchAddress("EdepEA_P",      &EdepEA_P);
  ///toyTree->SetBranchAddress("EtotE_P",       &EtotE_P);
  ///toyTree->SetBranchAddress("EdepE_Pg1",     &EdepE_Pg1);
  ///toyTree->SetBranchAddress("EdepEA_Pg1",    &EdepEA_Pg1);
  ///toyTree->SetBranchAddress("EtotE_Pg1",     &EtotE_Pg1);
  ///toyTree->SetBranchAddress("EdepE_Pg2",     &EdepE_Pg2);
  ///toyTree->SetBranchAddress("EdepEA_Pg2",    &EdepEA_Pg2);
  ///toyTree->SetBranchAddress("EtotE_Pg2",     &EtotE_Pg2);

  //toyTree->SetBranchAddress("LY_E_P",        &LY_E_P);
  //toyTree->SetBranchAddress("LY_E_Pg1",      &LY_E_Pg1);
  //toyTree->SetBranchAddress("LY_E_Pg2",      &LY_E_Pg2);
  //toyTree->SetBranchAddress("LY_E_N",        &LY_E_N);
  //toyTree->SetBranchAddress("LY_E_P_Sum",    &LY_E_P_Sum);
  //toyTree->SetBranchAddress("Opt_E_P_Sum",   &Opt_E_P_Sum);
  //toyTree->SetBranchAddress("FEE_E_P_Sum",   &FEE_E_P_Sum);
  //toyTree->SetBranchAddress("Scale_E_P_Sum", &Scale_E_P_Sum);
  toyTree->SetBranchAddress("Res_E_P_Sum",   &Res_E_P_Sum);
  toyTree->SetBranchAddress("Erec_P",   &Erec_P);
 
  //toyTree->SetBranchAddress("LY_E_Ng",    &LY_E_Ng);
  //toyTree->SetBranchAddress("Opt_E_Ng",   &Opt_E_Ng);
  //toyTree->SetBranchAddress("FEE_E_Ng",   &FEE_E_Ng);
  //toyTree->SetBranchAddress("Scale_E_Ng", &Scale_E_Ng);
  toyTree->SetBranchAddress("Res_E_Ng",   &Res_E_Ng);
  toyTree->SetBranchAddress("Erec_Ng",   &Erec_Ng);

  int entries_toyTree = toyTree->GetEntries();

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////

  TFile *roofile_h2d_ep_ratio2center = new TFile("./DATA/cell_eff/h2d_ep_ratio2center.root", "read");
  TH2D *h2d_Ep_ratio2center = (TH2D*)roofile_h2d_ep_ratio2center->Get("h2d_ep_ratio2center");

  TFile *roofile_h2d_ed_ratio2center = new TFile("./DATA/cell_eff/h2d_ed_ratio2center.root", "read");
  TH2D *h2d_Ed_ratio2center = (TH2D*)roofile_h2d_ed_ratio2center->Get("h2d_ed_ratio2center");

  TH2D *hist_findbin = new TH2D("hist_findbin","", R2_binnum,R2_lower,R2_upper, Z_binnum,Z_lower,Z_upper);

  const int MaxCellNum = 401;
  if( MaxCellNum<(R2_binnum+2)*(Z_binnum+2) ) {
    printf("\n\n MaxCellNum = %d, but must be >= %d\n\n", MaxCellNum, (R2_binnum+2)*(Z_binnum+2) );
    return;
  }
  TH1D *hEp_cc[MaxCellNum];
  TH1D *hEd_cc[MaxCellNum];
  TH1D *hEp_cc_clone[MaxCellNum];
  TH1D *hEd_cc_clone[MaxCellNum];

  const int MaxColumnNum = R2_binnum+1; //11
  TH1D *hEp_cl[MaxColumnNum];
  TH1D *hEd_cl[MaxColumnNum];

  double hEp_low = 0;
  double hEp_hgh = 12;
  int    hEp_bin = 120;  // any value works except for comparison with MC, which uses 120 bins
  const double hEp_low_mc = 1.5;  // formula used in ->Integral() assumes this is >= hEp_low and <hEp_hgh
  const double hEp_hgh_mc = 13.5;
  const int    hEp_bin_mc = 120;

  double hEd_low = 0;
  double hEd_hgh = 3.3; //12;
  int    hEd_bin = 120; //1200;  // any value works except for comparison with MC, which uses 120 bins
  const double hEd_low_mc = 0.0;
  const double hEd_hgh_mc = 3.3;
  const int    hEd_bin_mc = 120;

  TH1D *hEp = new TH1D("hEp", "hEp", hEp_bin, hEp_low, hEp_hgh);
  TH1D *hEd = new TH1D("hEd", "hEd", hEd_bin, hEd_low, hEd_hgh);

  for(int idx=0; idx<MaxCellNum; idx++)
  {
    roostr = TString::Format("hEp_cc_%03d",idx);
    hEp_cc[idx] = new TH1D(roostr, roostr, hEp_bin, hEp_low, hEp_hgh);

    roostr = TString::Format("hEd_cc_%03d",idx);
    hEd_cc[idx] = new TH1D(roostr, roostr, hEd_bin, hEd_low, hEd_hgh);


    roostr = TString::Format("hEp_cc_clone_%03d",idx);
    hEp_cc_clone[idx] = new TH1D(roostr, roostr, hEp_bin, hEp_low, hEp_hgh);

    roostr = TString::Format("hEd_cc_clone_%03d",idx);
    hEd_cc_clone[idx] = new TH1D(roostr, roostr, hEd_bin, hEd_low, hEd_hgh);
  }

  for(int idx=0; idx<MaxColumnNum; idx++)
  {
    roostr = TString::Format("hEp_cl_%03d", idx);
    hEp_cl[idx] = new TH1D(roostr, roostr, hEp_bin, hEp_low, hEp_hgh);

    roostr = TString::Format("hEd_cl_%03d", idx);
    hEd_cl[idx] = new TH1D(roostr, roostr, hEd_bin, hEd_low, hEd_hgh);
  }

  ////////////

  double distPD=0;
  for(int ientry=0; ientry<entries_toyTree; ientry++)
  {
    toyTree->GetEntry(ientry);

    cout.precision(3);
    if((ientry+1)%20000==0) cout<<" ---> processing00 "<<(ientry+1)*100./entries_toyTree<<"%"<<endl;

    //////
    //if( cap_Target==64 ) continue;  // reject nGd events
    if( Res_E_P_Sum!=Res_E_P_Sum || Res_E_Ng!=Res_E_Ng ) {
      cout<<" Res_E_P_Sum = "<< Res_E_P_Sum <<", Res_E_Ng = "<< Res_E_Ng << endl;
      continue;
    }

    /// Analysis cuts
    if( Erec_Ng<EdelayCutLow || Erec_Ng>EdelayCutHigh ) continue;
    if( TcapCutHigh>0 && (Tdiff>TcapCutHigh || Tdiff<TcapCutLow) ) continue;
    if( Erec_P<EpromptCutLow || Erec_P>EpromptCutHigh ) continue;
    if( distCut>0 ) {
      distPD = sqrt( (X_Ng-X_P)*(X_Ng-X_P) + (Y_Ng-Y_P)*(Y_Ng-Y_P) + (Z_Ng-Z_P)*(Z_Ng-Z_P) );
      if( distPD>distCut ) continue;
    }

    // cut events due to variation in OAV dimensions
    if( deltaRoav<0 && sqrt(X_P*X_P+Y_P*Y_P)<-deltaRoav ) continue;
    if( deltaRoav>0 && sqrt(X_P*X_P+Y_P*Y_P)>Rshield-deltaRoav ) continue;
    if( deltaZoav<0 && Z_P<ZreflectorBottom-deltaZoav ) continue;
    if( deltaZoav>0 && Z_P>ZreflectorTop-deltaZoav ) continue;

    ////// prompt
    usr_r2 = pow(deltaRoav + sqrt(X_P*X_P+Y_P*Y_P),2) * 1e-6;  // mm2 ---> m2
    //if( radialCut>0 && usr_r2>radialCut ) continue;
    usr_z  = (deltaZoav + Z_P) * 1e-3;  // mm ---> m
    global_bin_num = hist_findbin->FindBin(usr_r2, usr_z);
    hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);
    
    hEp_cc[global_bin_num]->Fill( Res_E_P_Sum );

    ////// delayed
    usr_r2 = pow(deltaRoav + sqrt(X_Ng*X_Ng+Y_Ng*Y_Ng),2) * 1e-6;  // mm2 ---> m2
    //if( radialCut>0 && usr_r2>radialCut ) continue;
    usr_z  = (deltaZoav + Z_Ng) * 1e-3;  // mm ---> m
    global_bin_num = hist_findbin->FindBin(usr_r2, usr_z);
    hist_findbin->GetBinXYZ(global_bin_num, local_xbin, local_ybin, local_zbin);

    hEd_cc[global_bin_num]->Fill( Res_E_Ng );  
    
  }

  cout.precision(7);
  cout<<"\nNumber of nH IBD events = "<< numTotalnH <<endl;
  cout<<"Number of events cut = "<< numCut <<endl;
  cout<<"Number of events cut (modified up cut) = "<< numCutModUp <<", (modified down cut) = "<< numCutModDown <<endl;
  cout<<"Efficiency = "<< 1.0-numCut/numTotalnH <<",  EfficiencyModUp = "<< 1.0-numCutModUp/numTotalnH <<",  EfficiencyModDown = "<< 1.0-numCutModDown/numTotalnH <<endl;
  cout<<"deltaEffUp = "<< (numCut-numCutModUp)/numTotalnH <<",  deltaEffDown = "<< (numCut-numCutModDown)/numTotalnH <<",  deltaEffAve = "<< (numCutModUp-numCutModDown)/(2.0*numTotalnH) <<"\n"<<endl;

  ////////////////////

  double entries=0, Effcontent=0;
  TFile *roofile_rawE = new TFile("roofile_rawE.root", "recreate");
  for(int idx=0; idx<MaxCellNum; idx++)
  {
    entries = 0;
    entries = hEp_cc[idx]->GetEntries();
    if( entries>10 )  
    {
      hEp_cc[idx]->Write();
      hEp_cc_clone[idx]->Add(hEp_cc[idx]);
    }

    entries = hEd_cc[idx]->GetEntries();
    if( entries>10 ) 
    {
      hEd_cc[idx]->Write();
      hEd_cc_clone[idx]->Add(hEd_cc[idx]);  
    }

    //////
    hist_findbin->GetBinXYZ(idx, local_xbin, local_ybin, local_zbin);
    if( local_zbin==0 && local_xbin>=1 && local_xbin<=R2_binnum && local_ybin>=1 && local_ybin<=Z_binnum )
    {
      Effcontent = 0;

      //if(local_ybin<=1 || local_ybin>=Z_binnum)  continue;  // exclude bottom and top rows from columns

      // normalize prompt events over energy range that matches that of MC (for comparison)
      entries = hEp_cc_clone[idx]->Integral( 1+ceil(hEp_bin*(hEp_low_mc-hEp_low)/(hEp_hgh-hEp_low)), hEp_bin);
      Effcontent = h2d_Ep_ratio2center->GetBinContent(local_xbin, local_ybin);
      hEp_cc_clone[idx]->Scale(Effcontent/entries);
      hEp_cl[local_xbin]->Add(hEp_cc_clone[idx]);
      hEp->Add(hEp_cc_clone[idx]);

      entries = hEd_cc_clone[idx]->Integral(1, hEd_bin);
      Effcontent = h2d_Ed_ratio2center->GetBinContent(local_xbin, local_ybin);
      hEd_cc_clone[idx]->Scale(Effcontent/entries);
      hEd_cl[local_xbin]->Add(hEd_cc_clone[idx]);
      hEd->Add(hEd_cc_clone[idx]);
    }
  }

  for(int idx=1; idx<MaxColumnNum; idx++)
  {
    hEp_cl[idx]->Write();
    hEd_cl[idx]->Write();
  }
  hEp->Write();
  hEd->Write();

  roofile_rawE->Close();

  ////////////////////

  TCanvas *canv_Ep_cc = new TCanvas("canv_Ep_cc", "canv_Ep_cc", 8000,8000);
  canv_Ep_cc->Divide(R2_binnum,Z_binnum,0,0);

  for(int idx=0; idx<MaxCellNum; idx++)
  {
    entries = 0;
    entries = hEp_cc[idx]->GetEntries();
    if( entries>10 )  
    {
      canv_Ep_cc->cd( RootCellToVisCell(idx) );

      hEp_cc[idx]->SetLineColor(2);
      hEp_cc[idx]->Draw();
    }
  }
  
  canv_Ep_cc->SaveAs("canv_Ep_cc.png");

  ////////////////////

  TCanvas *canv_Ed_cc = new TCanvas("canv_Ed_cc", "canv_Ed_cc", 8000,8000);
  canv_Ed_cc->Divide(R2_binnum,Z_binnum,0,0);

  for(int idx=0; idx<MaxCellNum; idx++)
  {
    entries = 0;

    entries = hEd_cc[idx]->GetEntries();
    if( entries>10 )  
    {
      canv_Ed_cc->cd( RootCellToVisCell(idx) );

      //hEd_cc[idx]->Rebin(2);
      hEd_cc[idx]->GetXaxis()->SetRangeUser(0.75*EdelayCutLow, min(1.3*EdelayCutHigh, 3.5) );  //(0, 3.5);
      hEd_cc[idx]->SetLineColor(2);
      hEd_cc[idx]->Draw();
    }
  }

  canv_Ed_cc->SaveAs("canv_Ed_cc.png");

  ////////////////////

  TCanvas *canv_Ep_cl = new TCanvas("canv_Ep_cl", "canv_Ep_cl", 4000,4000);
  canv_Ep_cl->Divide( ceil(sqrt((float)R2_binnum)), floor(sqrt((float)Z_binnum)), 0, 0 ); //(4,3,0,0);

  for(int idx=1; idx<MaxColumnNum; idx++)
  {
    canv_Ep_cl->cd( idx );
      
    //hEp_cl[idx]->Rebin(2);
    hEp_cl[idx]->SetLineColor(2);
    hEp_cl[idx]->Draw();
  }
  
  canv_Ep_cl->SaveAs("canv_Ep_cl.png");

  ////////////////////

  TCanvas *canv_Ed_cl = new TCanvas("canv_Ed_cl", "canv_Ed_cl", 4000,4000);
  canv_Ed_cl->Divide( ceil(sqrt((float)R2_binnum)), floor(sqrt((float)Z_binnum)), 0, 0 ); //(4,3,0,0);

  for(int idx=1; idx<MaxColumnNum; idx++)
  {
    canv_Ed_cl->cd( idx );

    //hEd_cl[idx]->Rebin(2);
    hEd_cl[idx]->GetXaxis()->SetRangeUser(0.75*EdelayCutLow, min(1.3*EdelayCutHigh, 3.5) ); //(0, 3.5);
    hEd_cl[idx]->SetLineColor(2);
    hEd_cl[idx]->Draw();
  }
  
  canv_Ed_cl->SaveAs("canv_Ed_cl.png");

  //////////////////////////////////////////////////////// mc_Eraw
  ////////////////////////////////////////////////////////
  
  TFile *roofile_mc_Eraw = new TFile("./DATA/mc_Eraw.root", "read");
  TH1D *hEp_mc_cc[MaxCellNum] = {NULL};
  TH1D *hEd_mc_cc[MaxCellNum] = {NULL};
  TH1D *hEp_mc_cl[MaxColumnNum] = {NULL};
  TH1D *hEd_mc_cl[MaxColumnNum] = {NULL};

  TH1D *hEp_mc = new TH1D("hEp_mc", "hEp_mc", hEp_bin_mc, hEp_low_mc, hEp_hgh_mc);
  TH1D *hEd_mc = new TH1D("hEd_mc", "hEd_mc", hEd_bin_mc, hEd_low_mc, hEd_hgh_mc);

  for(int idx=0; idx<MaxColumnNum; idx++)
  {
    roostr = TString::Format("hEp_mc_cl_%03d", idx);
    hEp_mc_cl[idx] = new TH1D(roostr, roostr, hEp_bin_mc, hEp_low_mc, hEp_hgh_mc);

    roostr = TString::Format("hEd_mc_cl_%03d", idx);
    hEd_mc_cl[idx] = new TH1D(roostr, roostr, hEd_bin_mc, hEd_low_mc, hEd_hgh_mc);
  }

  for(int idx=0; idx<MaxCellNum; idx++)
  {
    hist_findbin->GetBinXYZ(idx, local_xbin, local_ybin, local_zbin);

    entries    = 0;
    Effcontent = 0;

    roostr = TString::Format("hist_ep_sig_cell_ad_%03d_1", idx);
    if( roofile_mc_Eraw->Get(roostr)!=NULL )
    {
      hEp_mc_cc[idx] = (TH1D*)roofile_mc_Eraw->Get(roostr);
      entries = hEp_mc_cc[idx]->Integral(1, hEp_bin_mc);
      Effcontent = h2d_Ep_ratio2center->GetBinContent(local_xbin, local_ybin);

      hEp_mc_cc[idx]->Scale(Effcontent/entries);
      hEp_mc_cl[local_xbin]->Add(hEp_mc_cc[idx]);
      hEp_mc->Add(hEp_mc_cc[idx]);
    }

    roostr = TString::Format("hist_ed_sig_cell_ad_%03d_1", idx);
    if( roofile_mc_Eraw->Get(roostr)!=NULL )
    {
      hEd_mc_cc[idx] = (TH1D*)roofile_mc_Eraw->Get(roostr);
      entries = hEd_mc_cc[idx]->Integral(1, hEd_bin_mc);
      Effcontent = h2d_Ed_ratio2center->GetBinContent(local_xbin, local_ybin);

      hEd_mc_cc[idx]->Scale(Effcontent/entries);
      hEd_mc_cl[local_xbin]->Add(hEd_mc_cc[idx]);
      hEd_mc->Add(hEd_mc_cc[idx]);
    }
  }
 
  ////////////////////
  
  TCanvas *canv_Ep_mc_cc = new TCanvas("canv_Ep_mc_cc", "canv_Ep_mc_cc", 8000,8000);
  canv_Ep_mc_cc->Divide(R2_binnum,Z_binnum,0,0);

  for(int idx=0; idx<MaxCellNum; idx++)
  {
    if( hEp_mc_cc[idx]!=NULL )  
    {
      canv_Ep_mc_cc->cd( RootCellToVisCell(idx) );

      hEp_mc_cc[idx]->GetXaxis()->SetRangeUser(0.9*EpromptCutLow, 1.1*EpromptCutHigh);  //(1.35, 13.2);
      hEp_mc_cc[idx]->SetLineColor(2);
      hEp_mc_cc[idx]->Draw();
    }
  }
  
  canv_Ep_mc_cc->SaveAs("canv_Ep_mc_cc.png");

  ////////////////////
  
  TCanvas *canv_Ed_mc_cc = new TCanvas("canv_Ed_mc_cc", "canv_Ed_mc_cc", 8000,8000);
  canv_Ed_mc_cc->Divide(R2_binnum,Z_binnum,0,0);

  for(int idx=0; idx<MaxCellNum; idx++)
  {
    if( hEd_mc_cc[idx]!=NULL )  
    {
      canv_Ed_mc_cc->cd( RootCellToVisCell(idx) );

      //hEd_mc_cc[idx]->Rebin(2);
      hEd_mc_cc[idx]->GetXaxis()->SetRangeUser(0.75*EdelayCutLow, 1.3*EdelayCutHigh);  //(0, 3.5);
      hEd_mc_cc[idx]->SetLineColor(2);
      hEd_mc_cc[idx]->Draw();
    }
  }
  
  canv_Ed_mc_cc->SaveAs("canv_Ed_mc_cc.png");

  ////////////////////

  TCanvas *canv_Ep_mc_cl = new TCanvas("canv_Ep_mc_cl", "canv_Ep_mc_cl", 4000,4000);
  canv_Ep_mc_cl->Divide( ceil(sqrt((float)R2_binnum)), floor(sqrt((float)Z_binnum)), 0, 0 ); //(4,3,0,0);

  for(int idx=1; idx<MaxColumnNum; idx++)
  {
    canv_Ep_mc_cl->cd( idx );
      
    //hEp_mc_cl[idx]->Rebin(2);
    hEp_mc_cl[idx]->GetXaxis()->SetRangeUser(0.9*EpromptCutLow, 1.1*EpromptCutHigh);  //(1.35, 13.2);
    hEp_mc_cl[idx]->SetLineColor(2);
    hEp_mc_cl[idx]->Draw();
  }
  
  canv_Ep_mc_cl->SaveAs("canv_Ep_mc_cl.png");

  ////////////////////

  TCanvas *canv_Ed_mc_cl = new TCanvas("canv_Ed_mc_cl", "canv_Ed_mc_cl", 4000,4000);
  canv_Ed_mc_cl->Divide( ceil(sqrt((float)R2_binnum)), floor(sqrt((float)Z_binnum)), 0, 0 ); //(4,3,0,0);

  for(int idx=1; idx<MaxColumnNum; idx++)
  {
    canv_Ed_mc_cl->cd( idx );

    //hEd_mc_cl[idx]->Rebin(2);
    hEd_mc_cl[idx]->GetXaxis()->SetRangeUser(0.75*EdelayCutLow, 1.3*EdelayCutHigh);  //(0, 3.5);
    hEd_mc_cl[idx]->SetLineColor(2);
    hEd_mc_cl[idx]->Draw();
  }
  
  canv_Ed_mc_cl->SaveAs("canv_Ed_mc_cl.png");
  

  //////////////////////////////////////////////////////// mc_toy_cmp
  ////////////////////////////////////////////////////////

  ////////////////////

  TCanvas *canv_Ep_mc_toy_cl = new TCanvas("canv_Ep_mc_toy_cl", "canv_Ep_mc_toy_cl", 4000,4000);
  canv_Ep_mc_toy_cl->Divide( ceil(sqrt((float)R2_binnum)), floor(sqrt((float)Z_binnum)), 0, 0 ); //(4,3,0,0);

  for(int idx=1; idx<MaxColumnNum; idx++)
  {
    canv_Ep_mc_toy_cl->cd( idx );
    
    hEp_mc_cl[idx]->SetLineColor(1);
    hEp_mc_cl[idx]->Draw();

    hEp_cl[idx]->SetLineColor(2);
    hEp_cl[idx]->Draw("same");
  }
  
  canv_Ep_mc_toy_cl->SaveAs("canv_Ep_mc_toy_cl.png");

  // ---------------- Full spectrum
  TCanvas *canv_Ep_mc_toy = new TCanvas("canv_Ep_mc_toy", "canv_Ep_mc_toy", 4000,4000);
  
  hEp_mc->SetLineColor(1);
  hEp_mc->SetLineWidth(3);
  hEp_mc->Draw();

  hEp->SetLineColor(2);
  hEp->SetLineWidth(3);
  hEp->Draw("same");
  
  canv_Ep_mc_toy->SaveAs("canv_Ep_mc_toy.png");

  ////////////////////

  TCanvas *canv_Ed_mc_toy_cl = new TCanvas("canv_Ed_mc_toy_cl", "canv_Ed_mc_toy_cl", 4000,4000);
  canv_Ed_mc_toy_cl->Divide( ceil(sqrt((float)R2_binnum)), floor(sqrt((float)Z_binnum)), 0, 0 ); //(4,3,0,0);

  for(int idx=1; idx<MaxColumnNum; idx++)
  {
    canv_Ed_mc_toy_cl->cd( idx );
    
    hEd_mc_cl[idx]->SetLineColor(1);
    hEd_mc_cl[idx]->Draw();

    hEd_cl[idx]->SetLineColor(2);
    hEd_cl[idx]->Draw("same");
  }
  
  canv_Ed_mc_toy_cl->SaveAs("canv_Ed_mc_toy_cl.png");

  // ---------------- Full spectrum
  TCanvas *canv_Ed_mc_toy = new TCanvas("canv_Ed_mc_toy", "canv_Ed_mc_toy", 4000,4000);
  
  hEd_mc->SetLineColor(1);
  hEd_mc->SetLineWidth(3);
  hEd_mc->Draw();

  hEd->SetLineColor(2);
  hEd->SetLineWidth(3);
  hEd->Draw("same");
  
  canv_Ed_mc_toy->SaveAs("canv_Ed_mc_toy.png");


}
