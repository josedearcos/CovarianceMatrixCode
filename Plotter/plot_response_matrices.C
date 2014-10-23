#include <iostream>
#include "TCanvas.h"
#include "TH1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TColor.h"
#include "TH2D.h"
#include "TFile.h"
#include <string>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2.h"
#include "TMath.h"
#include <math.h>
#include <TArrayD.h>
#include <TMatrixD.h>
#include <TDecompChol.h>

int plot_response_matrices(void)
{
    TH1::AddDirectory(kFALSE);
    
    TH2D* NominalResponseH;
    TH2D* NominalTransResponseH;
    
    gStyle->SetErrorX(0.0001);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.05,"x");
    gStyle->SetTitleSize(0.05,"y");
    gStyle->SetLabelSize(0.04,"x");
    gStyle->SetLabelSize(0.04,"y");
    
    //To draw using a better palette:
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    
    //////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// ///////////
    ////////////                                            Nominal Response Matrix and Transverse
    //////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// ///////////
    TFile* PlotF = TFile::Open("../ResponseMatrices/NominalResponseMatrix.root");
    
    NominalResponseH=(TH2D*)gDirectory->Get("NoNormalizedEvisEnu0");
    NominalTransResponseH=(TH2D*)gDirectory->Get("NoNormalizedEnuEvis0");
    
    TCanvas* c3 = new TCanvas("NominalResponseMatrix","NominalResponseMatrix",500,500);
    gPad->SetLogz();
    NominalResponseH->GetXaxis()->SetTitle("E_{true} (MeV)");
    NominalResponseH->GetYaxis()->SetTitle("E_{vis} (MeV)");
    
    NominalResponseH->SetTitle("Evis-Enu Matrix");
    NominalResponseH->Draw("colz");
    c3->Print("../Images/NominalResponseMatrix.gif", "gif");
    
    TCanvas* c4 = new TCanvas("TransResponseMatrix","TransResponseMatrix",500,500);
    gPad->SetLogz();
    
    NominalTransResponseH->GetXaxis()->SetTitle("E_{vis} (MeV)");
    NominalTransResponseH->GetXaxis()->SetTitle("E_{true} (MeV)");
    
    NominalTransResponseH->SetTitle("Enu-Evis Matrix");
    NominalTransResponseH->Draw("colz");
    
    c4->Print("../Images/TransResponseMatrix.gif", "gif");
    
    PlotF->Close();
}