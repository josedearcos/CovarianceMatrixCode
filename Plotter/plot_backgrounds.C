#include <iostream>
#include "TCanvas.h"
#include "TH1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TColor.h"
#include "TH2F.h"
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

const Int_t NSamplesShown = 20;

int plot_backgrounds(void)
{
    
    
    
    string AnalysisString = "Gadolinium";
    
    
    
    TH1::AddDirectory(kFALSE);
    TH1D* VAccidentalVariations[NSamplesShown];
    TH1D* VLiHeVariations[NSamplesShown];
    TH1D* VAmCVariations[NSamplesShown];
    TH1D* VFNVariations[NSamplesShown];
    TH1D* DLiHeVariations[NSamplesShown];
    TH1D* DAmCVariations[NSamplesShown];
    TH1D* DFNVariations[NSamplesShown];
    TH1D* VAccNominal;
    TH1D* VFNNominal;
    TH1D* VLiHeNominal;
    TH1D* VAmCNominal;
    TH1D* DAmCNominal;
    TH1D* DLiHeNominal;
    TH1D* DFNNominal;
    TH1D* VAccRatio;
    TH1D* VAmCRatio;
    TH1D* VFNRatio;
    TH1D* VLiHeRatio;
    TH1D* DAmCRatio;
    TH1D* DFNRatio;
    TH1D* DLiHeRatio;
    
    gStyle->SetErrorX(0.0001);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.05,"x");
    gStyle->SetTitleSize(0.05,"y");
    gStyle->SetLabelSize(0.05,"x");
    gStyle->SetLabelSize(0.05,"y");
    //To draw using a better palette:
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    
    TFile* PlotF8 = TFile::Open(("../RootOutputs/"+AnalysisString+"/Backgrounds/Combine2/AccidentalVariations.root").c_str());
    
    VAccNominal=(TH1D*)gDirectory->Get("Nominal Accidental0_0 period0");
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        VAccidentalVariations[Sample]=(TH1D*)gDirectory->Get(Form("Accidental0_%d period0",Sample));
    }
    PlotF8->Close();
    
    TCanvas* c18 = new TCanvas("VAccidental","VAccidental",500,500);
    VAccNominal->SetTitle("Accidental Variations");
    VAccNominal->SetLineColor(kRed);
    VAccNominal->SetStats(0);
    VAccNominal->GetXaxis()->SetRange(1,36);
    VAccNominal->Draw();

    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        VAccidentalVariations[Sample]->SetStats(0);
        VAccidentalVariations[Sample]->Draw("same");
    }
    VAccNominal->Draw("same");
    
    c18->Print(("../Images/"+AnalysisString+"/BackgroundVariations/VAccidentalVariations.eps").c_str(), ".eps");
    delete c18;
    
    TFile* PlotF9 = TFile::Open(("../RootOutputs/"+AnalysisString+"/Backgrounds/Combine2/AmCVariations.root").c_str());
    VAmCNominal=(TH1D*)gDirectory->Get("Nominal AmCH0_0 period0");
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        VAmCVariations[Sample]=(TH1D*)gDirectory->Get(Form("AmCAD0_%d period0",Sample));
    }
    PlotF9->Close();
    
    TCanvas* c19 = new TCanvas("VAmC","VAmC",500,500);
    
    VAmCNominal->SetTitle("AmC Variations");
    VAmCNominal->SetLineColor(kRed);
    VAmCNominal->SetStats(0);
    VAmCNominal->GetXaxis()->SetRange(1,36);
    VAmCNominal->Draw();
    
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        VAmCVariations[Sample]->SetStats(0);
        VAmCVariations[Sample]->Draw("same");
    }
    VAmCNominal->Draw("same");
    
    c19->Print(("../Images/"+AnalysisString+"/BackgroundVariations/VAmCVariations.eps").c_str(), "eps");
    delete c19;
    
    
    TFile* PlotF10 = TFile::Open(("../RootOutputs/"+AnalysisString+"/Backgrounds/Combine2/LiHeVariations.root").c_str());
    
    VLiHeNominal=(TH1D*)gDirectory->Get("Nominal LiHe0_0 period0");
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        VLiHeVariations[Sample]=(TH1D*)gDirectory->Get(Form("LiHe0_%d period0",Sample));
    }
    PlotF10->Close();
    
    TCanvas* c20 = new TCanvas("VLiHe","VLiHe",500,500);

    VLiHeNominal->SetTitle("LiHe Variations");
    VLiHeNominal->SetLineColor(kRed);
    VLiHeNominal->SetStats(0);
    VLiHeNominal->GetYaxis()->SetRangeUser(0,20);
    VLiHeNominal->GetXaxis()->SetRange(1,36);
    VLiHeNominal->Draw();
    
    
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        VLiHeVariations[Sample]->SetStats(0);
        VLiHeVariations[Sample]->Draw("same");
    }
    VLiHeNominal->Draw("same");
    
    c20->Print(("../Images/"+AnalysisString+"/BackgroundVariations/VLiHeVariations.eps").c_str(), "eps");
    delete c20;
    TFile* PlotF11 = TFile::Open(("../RootOutputs/"+AnalysisString+"/Backgrounds/Combine2/FNVariations.root").c_str());
    TCanvas* c21 = new TCanvas("VFN","VFN",500,500);
    
    VFNNominal=(TH1D*)gDirectory->Get("Nominal FN0_0 period0");
    VFNNominal->SetTitle("FN Variations");
    VFNNominal->SetLineColor(kRed);
    VFNNominal->SetStats(0);
    VFNNominal->GetYaxis()->SetRangeUser(0.0,5);
    VFNNominal->GetXaxis()->SetRange(1,36);
    VFNNominal->Draw();
    
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        VFNVariations[Sample]=(TH1D*)gDirectory->Get(Form("FNAD0_%d period0",Sample));
        VFNVariations[Sample]->SetStats(0);
        VFNVariations[Sample]->Draw("same");
    }
    VFNNominal->Draw("same");
    c21->Print(("../Images/"+AnalysisString+"/BackgroundVariations/VFNVariations.eps").c_str(), "eps");
    delete c21;
    
    PlotF11->Close();
    
    
    TFile* PlotF12 = TFile::Open(("../RootOutputs/"+AnalysisString+"/Backgrounds/Combine2/AmCDistortions.root").c_str());
    TCanvas* c22 = new TCanvas("DAmC","DAmC",500,500);
    
    DAmCNominal=(TH1D*)gDirectory->Get("Nominal AmC0_0 period0");
    DAmCNominal->SetTitle("AmC Distortions");
    DAmCNominal->SetLineColor(kRed);
    DAmCNominal->SetStats(0);
    DAmCNominal->GetXaxis()->SetRange(1,36);
    DAmCNominal->Draw();
    
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        DAmCVariations[Sample]=(TH1D*)gDirectory->Get(Form("AmC0_%d period0",Sample));
        DAmCVariations[Sample]->SetStats(0);
        DAmCVariations[Sample]->Draw("same");
    }
    DAmCNominal->Draw("same");
    
    c22->Print(("../Images/"+AnalysisString+"/BackgroundVariations/DAmCVariations.eps").c_str(), "eps");
    delete c22;
    PlotF12->Close();
    
    TFile* PlotF13 = TFile::Open(("../RootOutputs/"+AnalysisString+"/Backgrounds/Combine2/LiHeDistortions.root").c_str());
    TCanvas* c23 = new TCanvas("DLiHe","DLiHe",500,500);
    
    DLiHeNominal=(TH1D*)gDirectory->Get("Nominal LiHe0_0 period0");
    DLiHeNominal->SetTitle("LiHe Distortions");
    DLiHeNominal->SetLineColor(kRed);
    DLiHeNominal->GetXaxis()->SetRange(1,36);
    DLiHeNominal->GetYaxis()->SetRangeUser(0,20);
    DLiHeNominal->SetStats(0);
    DLiHeNominal->Draw();
    
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        DLiHeVariations[Sample]=(TH1D*)gDirectory->Get(Form("LiHeAD0_%d period0",Sample));
        DLiHeVariations[Sample]->SetStats(0);
        DLiHeVariations[Sample]->Draw("same");
    }
    DLiHeNominal->Draw("same");
    
    c23->Print(("../Images/"+AnalysisString+"/BackgroundVariations/DLiHeVariations.eps").c_str(), "eps");
    delete c23;
    PlotF13->Close();
    
    TFile* PlotF14 = TFile::Open(("../RootOutputs/"+AnalysisString+"/Backgrounds/Combine2/FNDistortions.root").c_str());
    TCanvas* c24 = new TCanvas("DFN","DFN",500,500);
    
    DFNNominal=(TH1D*)gDirectory->Get("Nominal FN0_0 period0");
    DFNNominal->SetTitle("FN Distortions");
    DFNNominal->SetLineColor(kRed);
    DFNNominal->SetStats(0);
    DFNNominal->GetYaxis()->SetRangeUser(0,5);
    DFNNominal->GetXaxis()->SetRange(1,36);
    DFNNominal->Draw();
    
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        DFNVariations[Sample]=(TH1D*)gDirectory->Get(Form("FN0_%d period0",Sample));
        DFNVariations[Sample]->SetStats(0);
        DFNVariations[Sample]->Draw("same");
    }
    DFNNominal->Draw("same");
    
    c24->Print(("../Images/"+AnalysisString+"/BackgroundVariations/DFN.eps").c_str(), "eps");
    delete c24;
    PlotF14->Close();
    
    
    TCanvas* c25 = new TCanvas("VAccRatio","VAccRatio",500,500);
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        VAccRatio=(TH1D*)VAccidentalVariations[Sample]->Clone();
        VAccRatio->SetTitle("VAcc Ratio");
        VAccRatio->Add(VAccNominal,-1);
        VAccRatio->Divide(VAccNominal);
        if(Sample==0)
        {
            VAccRatio->GetYaxis()->SetRangeUser(-2,2);
            VAccRatio->SetStats(0);
            VAccRatio->Draw();
        }
        VAccRatio->Draw("same");
    }
    VAccRatio->Reset();
    VAccRatio->SetLineColor(kRed);
    VAccRatio->Draw("same");
    c25->Print(("../Images/"+AnalysisString+"/BackgroundVariations/VaccRatio.eps").c_str(), "eps");
    delete c25;
    TCanvas* c26 = new TCanvas("VAmCRatio","VAmCRatio",500,500);
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        VAmCRatio=(TH1D*)VAmCVariations[Sample]->Clone();
        VAmCRatio->SetTitle("VAmC Ratio");
        VAmCRatio->Add(VAmCNominal,-1);
        VAmCRatio->Divide(VAmCNominal);
        if(Sample==0)
        {
            VAmCRatio->GetYaxis()->SetRangeUser(-2,2);
            VAmCRatio->SetStats(0);
            VAmCRatio->Draw();
        }
        VAmCRatio->Draw("same");
    }
    VAmCRatio->Reset();
    VAmCRatio->SetLineColor(kRed);
    VAmCRatio->Draw("same");
    c26->Print(("../Images/"+AnalysisString+"/BackgroundVariations/VAmCRatio.eps").c_str(), "eps");
    delete c26;
    TCanvas* c27 = new TCanvas("VFNRatio","VFNRatio",500,500);
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        VFNRatio=(TH1D*)VFNVariations[Sample]->Clone();
        VFNRatio->SetTitle("VFN Ratio");
        VFNRatio->Add(VFNNominal,-1);
        VFNRatio->Divide(VFNNominal);
        if(Sample==0)
        {
            VFNRatio->GetYaxis()->SetRangeUser(-2,2);
            VFNRatio->SetStats(0);
            VFNRatio->Draw();
        }
        VFNRatio->Draw("same");
    }
    VFNRatio->Reset();
    VFNRatio->SetLineColor(kRed);
    VFNRatio->Draw("same");
    c27->Print(("../Images/"+AnalysisString+"/BackgroundVariations/VFNRatio.eps").c_str(), "eps");
    delete c27;
    TCanvas* c28 = new TCanvas("VLiHeRatio","VLiHeRatio",500,500);
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        VLiHeRatio=(TH1D*)VLiHeVariations[Sample]->Clone();
        VLiHeRatio->SetTitle("VLiHe Ratio");
        VLiHeRatio->Add(VLiHeNominal,-1);
        VLiHeRatio->Divide(VLiHeNominal);
        if(Sample==0)
        {
            VLiHeRatio->GetYaxis()->SetRangeUser(-2,2);
            VLiHeRatio->SetStats(0);
            VLiHeRatio->Draw();
        }
        VLiHeRatio->Draw("same");
    }
    VLiHeRatio->Reset();
    VLiHeRatio->SetLineColor(kRed);
    VLiHeRatio->Draw("same");
    c28->Print(("../Images/"+AnalysisString+"/BackgroundVariations/VLiHeRatio.eps").c_str(), "eps");
    delete c28;
    TCanvas* c29 = new TCanvas("DAmCRatio","DAmCRatio",500,500);
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        DAmCRatio=(TH1D*)DAmCVariations[Sample]->Clone();
        DAmCRatio->SetTitle("DAmC Ratio");
        DAmCRatio->Add(DAmCNominal,-1);
        DAmCRatio->Divide(DAmCNominal);
        if(Sample==0)
        {
            DAmCRatio->GetYaxis()->SetRangeUser(-10,10);
            DAmCRatio->SetStats(0);
            DAmCRatio->Draw();
        }
        DAmCRatio->Draw("same");
    }
    DAmCRatio->Reset();
    DAmCRatio->SetLineColor(kRed);
    DAmCRatio->Draw("same");
    c29->Print(("../Images/"+AnalysisString+"/BackgroundVariations/DAmCRatio.eps").c_str(), "eps");
    delete c29;
    TCanvas* c30 = new TCanvas("DFNRatio","DFNRatio",500,500);
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        DFNRatio=(TH1D*)DFNVariations[Sample]->Clone();
        DFNRatio->SetTitle("DFN Ratio");
        DFNRatio->Add(DFNNominal,-1);
        DFNRatio->Divide(DFNNominal);
        if(Sample==0)
        {
            DFNRatio->GetYaxis()->SetRangeUser(-2,2);
            DFNRatio->SetStats(0);
            DFNRatio->Draw();
        }
        DFNRatio->Draw("same");
    }
    DFNRatio->Reset();
    DFNRatio->SetLineColor(kRed);
    DFNRatio->Draw("same");
    c30->Print(("../Images/"+AnalysisString+"/BackgroundVariations/DFNRatio.eps").c_str(), "eps");
    delete c30;
    TCanvas* c31 = new TCanvas("DLiHeRatio","DLiHeRatio",500,500);
    for(Int_t Sample=0;Sample<NSamplesShown;Sample++)
    {
        DLiHeRatio=(TH1D*)DLiHeVariations[Sample]->Clone();
        DLiHeRatio->SetTitle("DLiHe Ratio");
        DLiHeRatio->Add(DLiHeNominal,-1);
        DLiHeRatio->Divide(DLiHeNominal);
        if(Sample==0)
        {
            DLiHeRatio->GetYaxis()->SetRangeUser(-2,2);
            DLiHeRatio->SetStats(0);
            DLiHeRatio->Draw();
        }
        DLiHeRatio->Draw("same");
    }
    DLiHeRatio->Reset();
    DLiHeRatio->SetLineColor(kRed);
    DLiHeRatio->Draw("same");
    c31->Print(("../Images/"+AnalysisString+"/BackgroundVariations/DLiHeRatio.eps").c_str(), "eps");
    delete c31;
}