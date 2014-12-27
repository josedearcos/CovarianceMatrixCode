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

int plot_flux_fractions()
{
    TH1::AddDirectory(kFALSE);
    
    TH1D* FluxH[4][6][20];
    TH1D* ExtraH[4][4][6];
      
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
    //////////// Factor Studies
    //////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// ///////////
    
    std::string labels[6]={"D1","D2","L1","L2","L3","L4"};
    Double_t colors[6] = {1,2,4,8,6,11};
    Char_t name[20];
    TFile* PlotF1 = TFile::Open("../RootOutputs/Gadolinium/FactorStudies.root");
    
    TCanvas *c1 = new TCanvas("FluxFractions","FluxFractions",1000,300);
    c1->Divide(3,1);
    TLegend *leg = new TLegend(0.1,0.1,0.4,0.5);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    for(Int_t week=0;week<1;++week)
    {
        for(Int_t near=0;near<3;++near)
        {
            c1->cd(near+1);
            for(Int_t reactor=0;reactor<6;++reactor)
            {
                FluxH[near][reactor][week]=(TH1D*)gDirectory->Get(Form("FluxFraction factor, Near AD%i, Reactor%i, Week%i", near+1, reactor+1, week+1));
                
                FluxH[near][reactor][week]->SetTitle("Flux Fractions");
                
                FluxH[near][reactor][week]->SetLineColor(colors[reactor]);
                FluxH[near][reactor][week]->GetXaxis()->SetTitle("E_{true} (MeV)");
                FluxH[near][reactor][week]->GetXaxis()->SetTitleSize(0.04);
                FluxH[near][reactor][week]->GetYaxis()->SetTitleSize(0.045);
                //            FluxH[near][reactor][week]->GetYaxis()->SetTitleOffset(1.35);
                FluxH[near][reactor][week]->GetXaxis()->SetLabelSize(0.045);
                FluxH[near][reactor][week]->GetYaxis()->SetLabelSize(0.045);
                FluxH[near][reactor][week]->SetTitle("");
                
                sprintf(name,"f_{%i,j}",near+1);
                FluxH[near][reactor][week]->GetYaxis()->SetTitle(name);
                if(reactor==0)
                {
                    FluxH[near][reactor][week]->GetYaxis()->SetRangeUser(0,0.5);
                    FluxH[near][reactor][week]->Draw();
                }
                FluxH[near][reactor][week]->Draw("same");
                if(near==0)
                {
                    leg->AddEntry(FluxH[near][reactor][week],labels[reactor].c_str(),"l");
                    leg->Draw("same");
                }
            }
        }
    }
    c1->Print("../Images/Gadolinium/FluxFractions.eps");
    
    PlotF1->Close();
    
    TFile* PlotF2 = TFile::Open("../RootOutputs/Gadolinium/FactorStudies.root");
    
    TCanvas *c2 = new TCanvas("Extrapolation","Extrapolation",900,900);
    c2->Divide(3,3);
    
    int cont=0;
    
    for(Int_t far=0;far<3;++far)
    {
        for(Int_t near=0;near<3;++near)
        {
            ++cont;
            c2->cd(cont);
            
            for(Int_t reactor=0;reactor<6;++reactor)
            {
                ExtraH[near][far][reactor]=(TH1D*)gDirectory->Get(Form("Extrapolation factor, Near AD%i, Far AD%i, Reactor%i", near+1, far+1, reactor+1));
                ExtraH[near][far][reactor]->SetTitle("Extrapolation factors");
                
                ExtraH[near][far][reactor]->SetLineColor(colors[reactor]);
                ExtraH[near][far][reactor]->GetXaxis()->SetTitle("E_{true} (MeV)");
                ExtraH[near][far][reactor]->GetXaxis()->SetTitleSize(0.040);
                ExtraH[near][far][reactor]->GetYaxis()->SetTitleSize(0.045);
                //                ExtraH[near][far][reactor]->GetYaxis()->SetTitleOffset(1.15);
                ExtraH[near][far][reactor]->GetXaxis()->SetLabelSize(0.045);
                ExtraH[near][far][reactor]->GetYaxis()->SetLabelSize(0.045);
                ExtraH[near][far][reactor]->SetTitle("");
                
                sprintf(name,"e_{%ij,%i}",near+1,far+1);
                
                ExtraH[near][far][reactor]->GetYaxis()->SetTitle(name);
                
                if(reactor==0)
                {
                    ExtraH[near][far][reactor]->GetYaxis()->SetRangeUser(0,0.885);
                    ExtraH[near][far][reactor]->Draw();
                }
                ExtraH[near][far][reactor]->Draw("same");
            }
        }
    }
    c2->Print("../Images/Gadolinium/ExtrapolationFactors.eps");
    
    PlotF2->Close();
}