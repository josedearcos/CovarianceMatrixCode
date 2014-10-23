#include "TCint.h"
#include "TGaxis.h"
#include <iostream>
#include <string>
#include <ctime>
#include <stdlib.h>
#include <TMath.h>
#include "TF1.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLeaf.h"
#include "Riostream.h"
#include "TTimeStamp.h"
#include "TH1D.h"
#include "TGraphErrors.h"

void AddIAVCorrection(){

  gStyle->SetOptStat(0);
  Double_t xval;
  Double_t yval;
  double junk1, junk2, junk3, junk4, junk5, junk6;

  //Read in reactor flux
  TString fullname;
  fullname = Form("./InputIBD.dat");
  cout<<fullname<<endl;
  ifstream datafile;

  ////Read in IAV Correction
  TFile *file = new TFile("./IAVDistortion.root");
  TH1F *Correciton = file->Get("Correction");
  TH1F *Uncorr = file->Get("UncUncorr");
  
  ////Make histogram of reactor data
  TH1F *Reactor = new TH1F("Reactor","Reactor",240,0,12);
  datafile.open(fullname, ios::in);
  for (Int_t j=0; j<240 ; j++) {
    datafile >> junk1 >> junk2 >> junk3 >> yval >> xval;
    //cout<<xval+0.01<<" "<<yval<<endl;
    Reactor->Fill(xval+0.01,yval);
  }
  datafile.close();

  //////////////////////////////////////////////////////////
  ///////////Add IAV Correction to Histogram: IMPORTANT PART
  //////////////////////////////////////////////////////////
  double CorrectedPercent=0;
  double CorrectedPercentUnc=0;
  double TotalVal=0;
  double StatUnc=0;
  TH1F *ReactorCorr = new TH1F("ReactorCorr","ReactorCorr",240,0,12);
  TH1F *Corr = new TH1F("Corr","Corr",240,0,12);
  for(int i=21;i<241;i++){
    //Total events in input spectrum bin i
    TotalVal = Reactor->GetBinContent(i);
    for(int j=0;j<i+1;j++){
      //Percent of events in bin i moved to bin j
      CorrectedPercent = Correction->GetBinContent(i,j)/Correction->Integral(i,i,0,-1);
      CorrectedPercentUnc = TMath::Sqrt(Correction->GetBinContent(i,j))/Correction->Integral(i,i,0,-1);
      //Add percent of total of input spectrum bin i to corrected spectrum bin j
      ReactorCorr->AddBinContent(j,TotalVal*CorrectedPercent);
      
      //Add statistical uncertainty to bin j
      //Get Existing Uncertainty in bin j
      StatUnc = TMath::Sqrt(TMath::Power(ReactorCorr->GetBinError(j),2)+TMath::Power(TotalVal*CorrectedPercentUnc,2));      
      ReactorCorr->SetBinError(j,StatUnc);
    }    
  }
  
  //Make correlated uncertainty histogram
  for(int i=0;i<241;i++)
    Corr->SetBinContent(i,ReactorCorr->GetBinError(i)/ReactorCorr->GetBinContent(i));
  

  //Draw Histograms
  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,500);
  TCanvas *c2 = new TCanvas("c2","c2",0,0,600,500);
  TCanvas *c3 = new TCanvas("c3","c3",0,0,600,500);
  TLegend *leg = new TLegend(0.5,0.5,0.8,0.7);
  c1->cd();
  c1->SetLogy();
  //ReactorCorr->Divide(Reactor);
  //for(int i=0;i<241;i++)
  //  ReactorCorr->SetBinError(i,0);
  Reactor->Draw("");  
  ReactorCorr->SetLineColor(kRed);
  ReactorCorr->Draw("Same");
  leg->AddEntry(Reactor,"Input IBD Spectrum","EL");
  leg->AddEntry(ReactorCorr,"After IAV Correction","EL");
  leg->Draw();
  c2->cd();
  Corr->Draw("");
  c3->cd();
  Uncorr->Draw("");
  
  //cout<<Reactor->Integral(0,-1)<<endl;
  //cout<<ReactorCorr->Integral(0,-1)<<endl;
  //cout<<(Reactor->Integral(0,-1)-ReactorCorr->Integral(0,-1))/Reactor->Integral(0,-1)<<endl;
}
