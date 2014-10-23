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
#include "TH1F.h"
#include "TH2.h"
#include "TMath.h"
#include <math.h>
#include <TArrayD.h>
#include <TMatrixD.h>
#include <TDecompChol.h>

const Int_t MaxSamples = 20;

int plot_systematics(void)
{
    
    TH1::AddDirectory(kFALSE);

    Int_t n_evis_bins =37;
    Int_t n_etrue_bins = 39;
    
    TH2F* PosResponseH;
    TH2F* NLResponseH;
    TH2F* IAVResponseH;
    TH2F* ResoResponseH;
    
    TH1F* PosSpectrumH;
    TH1F* IAVSpectrumH;
    TH1F* NLSpectrumH;
    TH1F* ResoSpectrumH;
    
    TH2F* FinePosResponseH;
    TH2F* FineNLResponseH;
    TH2F* FineIAVResponseH;
    TH2F* FineResoResponseH;
    
    TH1F* FinePosSpectrumH;
    TH1F* FineIAVSpectrumH;
    TH1F* FineNLSpectrumH;
    TH1F* FineResoSpectrumH;
    
    TH1F* NominalSpectrumH;
    TH1F* NominalTrueSpectrumH;
  
    TH1F* FineNominalSpectrumH;
    TH1F* FineNominalTrueSpectrumH;
    
    TH1F* NLVariations[MaxSamples];
    TH1F* ResoVariations[MaxSamples];
    TH1F* IAVVariations[MaxSamples];
    TH1F* IsotopeVariations[MaxSamples];
    TH1F* PowerVariations[MaxSamples];
    TH1F* RelativeEnergyVariations[MaxSamples];

    TH1F* ResoNominal;
    TH1F* IAVNominal;
    TH1F* NLNominal;
    TH1F* IsotopeNominal;
    TH1F* PowerNominal;
    TH1F* ResoRatio;
    TH1F* IAVRatio;
    TH1F* NLRatio;
    TH1F* IsotopeRatio;
    TH1F* PowerRatio;
    TH1F* RelativeEnergyRatio;
    
    Double_t enu_bins[39+1];
    Double_t evis_bins[37+1];
    
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
    
    
    n_evis_bins=37;
    n_etrue_bins=39;
    
    for (Int_t i = 0; i <= n_etrue_bins; i++)
    {
        enu_bins[i] = 0.2 * i + 1.8;
    }
    
    evis_bins[0] = 0.7;
    for (Int_t i = 0; i < n_evis_bins-1; i++)
    {
        evis_bins[i+1] = 0.2 * i + 1.0;
    }
    evis_bins[n_evis_bins] = 12;
    
    ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// //
    //////////////////                                                      Systematic Variations
    ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// //
    
    TFile* PlotF3 = TFile::Open("../CovarianceMatrices/Combine1/Spectrum/NL.root");
    
    NLNominal=(TH1F*)gDirectory->Get("Nominal Far AD0 prediction from Near AD0 VisH Sample0");
    NLNominal->GetXaxis()->SetTitle("E_{vis} (MeV)");
    NLNominal->SetTitle("NL Variations");
    NLNominal->SetLineColor(kRed);
    TCanvas* c5 = new TCanvas("NL","NL",400,400);
    
    for(Int_t Sample=0;Sample<MaxSamples;Sample++)
    {
        NLNominal->Draw("same");
        NLVariations[Sample]=(TH1F*)gDirectory->Get(Form("Varied Far AD0 prediction from Near AD0 VisH Sample%i",Sample));
        NLVariations[Sample]->Draw("same");
        NLNominal->Draw("same");
    }
    c5->Print("../Images/NLVariations.eps", "png");
    
    PlotF3->Close();
    
    TFile* PlotF4 = TFile::Open("../CovarianceMatrices/Combine1/Spectrum/Resolution.root");
    
    ResoNominal=(TH1F*)gDirectory->Get("Nominal Far AD0 prediction from Near AD0 VisH Sample0");
    ResoNominal->GetXaxis()->SetTitle("E_{vis} (MeV)");
    ResoNominal->SetTitle("Resolution Variations");
    ResoNominal->SetLineColor(kRed);
    
    TCanvas* c6 = new TCanvas("Reso","Reso",400,400);
    
    for(Int_t Sample=0;Sample<MaxSamples;Sample++)
    {
        ResoNominal->Draw("same");
        ResoVariations[Sample]=(TH1F*)gDirectory->Get(Form("Varied Far AD0 prediction from Near AD0 VisH Sample%i",Sample));
        ResoVariations[Sample]->Draw("same");
        ResoNominal->Draw("same");
    }
    c6->Print("../Images/ResoVariations.eps", "png");
    
    PlotF4->Close();
    
    TFile* PlotF5 = TFile::Open("../CovarianceMatrices/Combine1/Spectrum/IAV.root");
    
    IAVNominal=(TH1F*)gDirectory->Get("Nominal Far AD0 prediction from Near AD0 VisH Sample0");
    IAVNominal->GetXaxis()->SetTitle("E_{vis} (MeV)");
    IAVNominal->SetTitle("IAV Variations");
    IAVNominal->SetLineColor(kRed);
    
    TCanvas* c7 = new TCanvas("IAV","IAV",400,400);
    for(Int_t Sample=0;Sample<MaxSamples;Sample++)
    {
        IAVNominal->Draw("same");
        IAVVariations[Sample]=(TH1F*)gDirectory->Get(Form("Varied Far AD0 prediction from Near AD0 VisH Sample%i",Sample));
        IAVVariations[Sample]->Draw("same");
        IAVNominal->Draw("same");
    }
    c7->Print("../Images/IAVVariations.eps", "png");
    
    PlotF5->Close();
    
    TFile* PlotF6 = TFile::Open("../CovarianceMatrices/Combine1/Spectrum/Isotope.root");
    
    IsotopeNominal=(TH1F*)gDirectory->Get("Nominal Far AD0 prediction from Near AD0 VisH Sample0");
    IsotopeNominal->GetXaxis()->SetTitle("E_{vis} (MeV)");
    IsotopeNominal->SetTitle("Isotope Variations");
    IsotopeNominal->SetLineColor(kRed);
    
    TCanvas* c11 = new TCanvas("Isotope","Isotope",400,400);
    for(Int_t Sample=0;Sample<MaxSamples;Sample++)
    {
        IsotopeNominal->Draw("same");
        IsotopeVariations[Sample]=(TH1F*)gDirectory->Get(Form("Varied Far AD0 prediction from Near AD0 VisH Sample%i",Sample));
        IsotopeVariations[Sample]->Draw("same");
        IsotopeNominal->Draw("same");
    }
    c11->Print("../Images/IsotopeVariations.eps", "png");
    
    PlotF6->Close();
    
    TFile* PlotF7 = TFile::Open("../CovarianceMatrices/Combine1/Spectrum/Power.root");
    
    PowerNominal=(TH1F*)gDirectory->Get("Nominal Far AD0 prediction from Near AD0 VisH Sample0");
    PowerNominal->GetXaxis()->SetTitle("E_{vis} (MeV)");
    PowerNominal->SetTitle("Reactor Power Variations");
    PowerNominal->SetLineColor(kRed);
    
    TCanvas* c12 = new TCanvas("Power","Power",400,400);
    for(Int_t Sample=0;Sample<MaxSamples;Sample++)
    {
        PowerNominal->Draw("same");
        PowerVariations[Sample]=(TH1F*)gDirectory->Get(Form("Varied Far AD0 prediction from Near AD0 VisH Sample%i",Sample));
        PowerVariations[Sample]->Draw("same");
        PowerNominal->Draw("same");
    }
    c12->Print("../Images/PowerVariations.eps", "png");
    
    PlotF7->Close();
    
    TFile* PlotF8 = TFile::Open("../CovarianceMatrices/Combine1/Spectrum/RelativeEnergyScale.root");
    
    RelativeEnergyNominal=(TH1F*)gDirectory->Get("Nominal Far AD0 prediction from Near AD0 VisH Sample0");
    RelativeEnergyNominal->GetXaxis()->SetTitle("E_{vis} (MeV)");
    RelativeEnergyNominal->SetTitle("RelativeEnergy Variations");
    RelativeEnergyNominal->SetLineColor(kRed);
    TCanvas* c13 = new TCanvas("RelativeEnergy","RelativeEnergy",400,400);
    
    for(Int_t Sample=0;Sample<MaxSamples;Sample++)
    {
        RelativeEnergyNominal->Draw("same");
        RelativeEnergyVariations[Sample]=(TH1F*)gDirectory->Get(Form("Varied Far AD0 prediction from Near AD0 VisH Sample%i",Sample));
        RelativeEnergyVariations[Sample]->Draw("same");
        RelativeEnergyNominal->Draw("same");
    }
    c13->Print("../Images/RelativeEnergyVariations.eps", "png");
    
    PlotF8->Close();
    
    ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// //
    //////////////////                                                     Systematic Ratios
    ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// /// ////////////////// //
    
    TCanvas* c113 = new TCanvas("ResoRatio","ResoRatio",400,400);
    for(Int_t Sample=0;Sample<MaxSamples;Sample++)
    {
        ResoRatio=(TH1F*)ResoVariations[Sample]->Clone();
        ResoRatio->SetTitle("Resolution Ratio");
        ResoRatio->Add(ResoNominal,-1);
        ResoRatio->Divide(ResoNominal);
        if(Sample==0)
        {
            ResoRatio->GetYaxis()->SetRangeUser(-0.1,0.1);
            ResoRatio->Draw();
        }
        ResoRatio->Draw("same");
    }
    ResoRatio->Reset();
    ResoRatio->SetLineColor(kRed);
    ResoRatio->Draw("same");
    c113->Print("../Images/ResoRatio.eps", "png");
    
    TCanvas* c14 = new TCanvas("NLRatio","NLRatio",400,400);
    for(Int_t Sample=0;Sample<MaxSamples;Sample++)
    {
        NLRatio=(TH1F*)NLVariations[Sample]->Clone();
        NLRatio->SetTitle("NL Ratio");
        NLRatio->Add(NLNominal,-1);
        NLRatio->Divide(NLNominal);
        if(Sample==0)
        {
            NLRatio->GetYaxis()->SetRangeUser(-0.1,0.1);
            NLRatio->Draw();
        }
        NLRatio->Draw("same");
    }
    NLRatio->Reset();
    NLRatio->SetLineColor(kRed);
    NLRatio->Draw("same");
    
    c14->Print("../Images/NLRatio.eps", "png");
    
    TCanvas* c15 = new TCanvas("IAVRatio","IAVRatio",400,400);
    for(Int_t Sample=0;Sample<MaxSamples;Sample++)
    {
        IAVRatio=(TH1F*)IAVVariations[Sample]->Clone();
        IAVRatio->SetTitle("IAV Ratio");
        IAVRatio->Add(IAVNominal,-1);
        IAVRatio->Divide(IAVNominal);
        if(Sample==0)
        {
            IAVRatio->GetYaxis()->SetRangeUser(-0.1,0.1);
            IAVRatio->Draw();
        }
        IAVRatio->Draw("same");
    }
    IAVRatio->Reset();
    IAVRatio->SetLineColor(kRed);
    IAVRatio->Draw("same");
    c15->Print("../Images/IAVRatio.eps", "png");
    
    TCanvas* c16= new TCanvas("IsotopeRatio","IsotopeRatio",400,400);
    for(Int_t Sample=0;Sample<MaxSamples;Sample++)
    {
        IsotopeRatio=(TH1F*)IsotopeVariations[Sample]->Clone();
        IsotopeRatio->SetTitle("Isotope Ratio");
        IsotopeRatio->Add(IAVNominal,-1);
        IsotopeRatio->Divide(IAVNominal);
        if(Sample==0)
        {
            IsotopeRatio->GetYaxis()->SetRangeUser(-0.1,0.1);
            IsotopeRatio->Draw();
        }
        IsotopeRatio->Draw("same");
    }
    IsotopeRatio->Reset();
    IsotopeRatio->SetLineColor(kRed);
    IsotopeRatio->Draw("same");
    c16->Print("../Images/IsotopeRatio.eps", "png");
    
    TCanvas* c17 = new TCanvas("PowerRatio","PowerRatio",400,400);
    for(Int_t Sample=0;Sample<MaxSamples;Sample++)
    {
        PowerRatio=(TH1F*)PowerVariations[Sample]->Clone();
        PowerRatio->SetTitle("Power Ratio");
        PowerRatio->Add(PowerNominal,-1);
        PowerRatio->Divide(PowerNominal);
        if(Sample==0)
        {
            PowerRatio->GetYaxis()->SetRangeUser(-0.1,0.1);
            PowerRatio->Draw();
        }
        PowerRatio->Draw("same");
    }
    PowerRatio->Reset();
    PowerRatio->SetLineColor(kRed);
    PowerRatio->Draw("same");
    c17->Print("../Images/PowerRatio.eps", "png");
    
    TCanvas* c18 = new TCanvas("RelativeEnergyRatio","RelativeEnergyRatio",400,400);
    for(Int_t Sample=0;Sample<MaxSamples;Sample++)
    {
        RelativeEnergyRatio=(TH1F*)RelativeEnergyVariations[Sample]->Clone();
        RelativeEnergyRatio->SetTitle("RelativeEnergy Ratio");
        RelativeEnergyRatio->Add(RelativeEnergyNominal,-1);
        RelativeEnergyRatio->Divide(RelativeEnergyNominal);
        if(Sample==0)
        {
            RelativeEnergyRatio->GetYaxis()->SetRangeUser(-0.1,0.1);
            RelativeEnergyRatio->Draw();
        }
        RelativeEnergyRatio->Draw("same");
    }
    RelativeEnergyRatio->Reset();
    RelativeEnergyRatio->SetLineColor(kRed);
    RelativeEnergyRatio->Draw("same");
    c18->Print("../Images/RelativeEnergyRatio.eps", "png");
    
    //////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// ///////////
    ////////////                                                Varied Response Matrices
    //////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// ///////////

    TFile* PlotF101 = TFile::Open("../ResponseMatrices/ResponseMatrix.root");
    
    FinePosResponseH=(TH2F*)gDirectory->Get("FineEvisEnuPos0");
    FineIAVResponseH=(TH2F*)gDirectory->Get("FineEvisEnuIAV0");
    FineNLResponseH=(TH2F*)gDirectory->Get("FineEvisEnuNL0");
    FineResoResponseH=(TH2F*)gDirectory->Get("FineEvisEnuReso0");

    PosResponseH=(TH2F*)gDirectory->Get("EvisEnuPos0");
    IAVResponseH=(TH2F*)gDirectory->Get("EvisEnuIAV0");
    NLResponseH=(TH2F*)gDirectory->Get("EvisEnuNL0");
    ResoResponseH=(TH2F*)gDirectory->Get("EvisEnuReso0");
    
    TH2F* NoFinePosResponseH=(TH2F*)gDirectory->Get("FinePosNoNormalizedEvisEnu0");
    TH2F* NoFineIAVResponseH=(TH2F*)gDirectory->Get("FineIAVNoNormalizedEvisEnu0");
    TH2F* NoFineNLResponseH=(TH2F*)gDirectory->Get("FineNLNoNormalizedEvisEnu0");
    TH2F* NoFineResoResponseH=(TH2F*)gDirectory->Get("FineResoNoNormalizedEvisEnu0");
    
    TH2F* NoPosResponseH=(TH2F*)gDirectory->Get("PosNoNormalizedEvisEnu0");
    TH2F* NoIAVResponseH=(TH2F*)gDirectory->Get("IAVNoNormalizedEvisEnu0");
    TH2F* NoNLResponseH=(TH2F*)gDirectory->Get("NLNoNormalizedEvisEnu0");
    TH2F* NoResoResponseH=(TH2F*)gDirectory->Get("ResoNoNormalizedEvisEnu0");

    PlotF101->Close();
    
    TFile* PlotF102 = TFile::Open("../RootOutputs/NominalOutputs/Oscillation.root");
        PlotF102->cd("Total AD Spectra after oscillation");
        FineNominalTrueSpectrumH=(TH1F*)gDirectory->Get("Total spectrum after oscillation at AD1");
        FineNominalTrueSpectrumH->SetLineColor(kRed);
    PlotF102->Close();
    
    PosSpectrumH = new TH1F("Pos Spectrum","Pos Spectrum", n_evis_bins,evis_bins);
    IAVSpectrumH = new TH1F("IAV Spectrum","IAV Spectrum", n_evis_bins,evis_bins);
    NLSpectrumH = new TH1F("NL Spectrum","NL Spectrum", n_evis_bins,evis_bins);
    ResoSpectrumH = new TH1F("Reso Spectrum","Reso Spectrum", n_evis_bins,evis_bins);
    
    FinePosSpectrumH = new TH1F("Fine Pos Spectrum","Fine Pos Spectrum", 240, 0,12);
    FineIAVSpectrumH = new TH1F("Fine IAV Spectrum","Fine IAV Spectrum", 240, 0,12);
    FineNLSpectrumH = new TH1F("Fine NL Spectrum","Fine NL Spectrum", 240, 0,12);
    FineResoSpectrumH = new TH1F("Fine Reso Spectrum","Fine Reso Spectrum", 240, 0,12);
    
    TH1F* FlatSpectrumH = new TH1F("Flat Spectrum","Flat Spectrum",240,0,12);
    for(Int_t i=1; i<=240; i++)
    {
        FlatSpectrumH->SetBinContent(i,1);
    }
//    for(Int_t i=1; i<=240; i++)
//    {
//        if(i>36)
//        {
//            FineNominalTrueSpectrumH->SetBinContent(i,FineNominalSpectrumH->GetBinContent(i-36));
//        }
//        else
//        {
//            FineNominalTrueSpectrumH->SetBinContent(i,0);
//        }
//    }
    NominalTrueSpectrumH=(TH1F*)FineNominalTrueSpectrumH->Rebin(n_etrue_bins,"Total spectrum after oscillation at AD1",enu_bins);
//
    for(Int_t i=1; i<=240; i++)
    {
        for(Int_t j=1; j<=240; j++)
        {

            FinePosSpectrumH->SetBinContent(i, FinePosSpectrumH->GetBinContent(i)+FinePosResponseH->GetBinContent(i,j)* FineNominalTrueSpectrumH->GetBinContent(j));
            FineIAVSpectrumH->SetBinContent(i, FineIAVSpectrumH->GetBinContent(i)+FineIAVResponseH->GetBinContent(i,j)* FineNominalTrueSpectrumH->GetBinContent(j));
            FineNLSpectrumH->SetBinContent(i, FineNLSpectrumH->GetBinContent(i)+FineNLResponseH->GetBinContent(i,j)*FineNominalTrueSpectrumH->GetBinContent(j));
            FineResoSpectrumH->SetBinContent(i, FineResoSpectrumH->GetBinContent(i)+FineResoResponseH->GetBinContent(i,j)* FineNominalTrueSpectrumH->GetBinContent(j));
        }
    }
    
    for(Int_t i=1; i<=n_evis_bins; i++)
    {
        for(Int_t j=1; j<=n_etrue_bins; j++)
        {
            PosSpectrumH->SetBinContent(i, PosSpectrumH->GetBinContent(i)+PosResponseH->GetBinContent(i,j)* NominalTrueSpectrumH->GetBinContent(j));
            IAVSpectrumH->SetBinContent(i, IAVSpectrumH->GetBinContent(i)+IAVResponseH->GetBinContent(i,j)* NominalTrueSpectrumH->GetBinContent(j));
            NLSpectrumH->SetBinContent(i, NLSpectrumH->GetBinContent(i)+NLResponseH->GetBinContent(i,j)*NominalTrueSpectrumH->GetBinContent(j));
            ResoSpectrumH->SetBinContent(i, ResoSpectrumH->GetBinContent(i)+ResoResponseH->GetBinContent(i,j)* NominalTrueSpectrumH->GetBinContent(j));
        }
    }

    
//    for(Int_t i=1; i<=240; i++)
//    {
//        for(Int_t j=1; j<=240; j++)
//        {
//            FinePosSpectrumH->SetBinContent(i, FinePosSpectrumH->GetBinContent(i)+NoFinePosResponseH->GetBinContent(i,j)* FlatSpectrumH->GetBinContent(j));
//            FineIAVSpectrumH->SetBinContent(i, FineIAVSpectrumH->GetBinContent(i)+NoFineIAVResponseH->GetBinContent(i,j)* FlatSpectrumH->GetBinContent(j));
//            FineNLSpectrumH->SetBinContent(i, FineNLSpectrumH->GetBinContent(i)+NoFineNLResponseH->GetBinContent(i,j)* FlatSpectrumH->GetBinContent(j));
//            FineResoSpectrumH->SetBinContent(i, FineResoSpectrumH->GetBinContent(i)+NoFineResoResponseH->GetBinContent(i,j)* FlatSpectrumH->GetBinContent(j));
//        }
//    }
//    
//    for(Int_t i=1; i<=n_evis_bins; i++)
//    {
//        for(Int_t j=1; j<=n_etrue_bins; j++)
//        {
//            PosSpectrumH->SetBinContent(i, PosSpectrumH->GetBinContent(i)+NoPosResponseH->GetBinContent(i,j)* FlatSpectrumH->GetBinContent(j));
//            IAVSpectrumH->SetBinContent(i, IAVSpectrumH->GetBinContent(i)+NoIAVResponseH->GetBinContent(i,j)* FlatSpectrumH->GetBinContent(j));
//            NLSpectrumH->SetBinContent(i, NLSpectrumH->GetBinContent(i)+NoNLResponseH->GetBinContent(i,j)* FlatSpectrumH->GetBinContent(j));
//            ResoSpectrumH->SetBinContent(i, ResoSpectrumH->GetBinContent(i)+NoResoResponseH->GetBinContent(i,j)* FlatSpectrumH->GetBinContent(j));
//        }
//    }

    
    FinePosSpectrumH->SetTitle("Pos Spectrum");
    FineIAVSpectrumH->SetTitle("IAV Spectrum");
    FineNLSpectrumH->SetTitle("NL Spectrum");
    FineResoSpectrumH->SetTitle("Reso Spectrum");
    
    FinePosSpectrumH->SetLineColor(kBlue);
    FineIAVSpectrumH->SetLineColor(kMagenta);
    FineNLSpectrumH->SetLineColor(kBlack);
    FineResoSpectrumH->SetLineColor(kGreen);
//
//    TCanvas* a9 = new TCanvas("Fine Positron effect","Fine Positron effect",400,400);
//    FinePosSpectrumH->Draw();
//    a9->Print("../Images/FinePosResponseSpectrum.eps","png");
//    
//    TCanvas* a10 = new TCanvas("Fine IAV effect","Fine IAV effect",400,400);
//    FineIAVSpectrumH->Draw("same");
//    a10->Print("../Images/FineIAVResponseSpectrum.eps", "png");
//    
//    TCanvas* a11 = new TCanvas("Fine NL effect","Fine NL effect",400,400);
//    FineNLSpectrumH->Draw("same");
//    a11->Print("../Images/FineNLResponseSpectrum.eps", "png");
//    
//    TCanvas* a12 = new TCanvas("Fine Reso effect","Fine Reso effect",400,400);
//    FineResoSpectrumH->Draw("same");
//    a12->Print("../Images/FineResoResponseSpectrum.eps", "png");
    
    TCanvas* a8 = new TCanvas("Fine Comparison neutrino - pos","Fine Comparison neutrino - pos",400,400);
    FinePosSpectrumH->Draw("same");
    FineNominalTrueSpectrumH->Draw("same");
    
    TCanvas* a5 = new TCanvas("Fine Comparison Pos - IAV","Fine Comparison Pos - IAV",400,400);
    FinePosSpectrumH->Draw();
    FineIAVSpectrumH->Draw("same");
    
    TCanvas* a6 = new TCanvas("Fine Comparison IAV - NL","Fine Comparison IAV - NL",400,400);
    FineIAVSpectrumH->Draw();
    FineNLSpectrumH->Draw("same");

    TCanvas* a7 = new TCanvas("Fine Comparison NL - Reso","Fine Comparison NL - Reso",400,400);
    FineNLSpectrumH->Draw();
    FineResoSpectrumH->Draw("same");
    
    PosSpectrumH->SetTitle("Pos Spectrum");
    IAVSpectrumH->SetTitle("IAV Spectrum");
    NLSpectrumH->SetTitle("NL Spectrum");
    ResoSpectrumH->SetTitle("Reso Spectrum");
    
    PosSpectrumH->SetLineColor(kBlue);
    IAVSpectrumH->SetLineColor(kMagenta);
    NLSpectrumH->SetLineColor(kBlack);
    ResoSpectrumH->SetLineColor(kGreen);
//    
//    TCanvas* a1 = new TCanvas("Positron effect","Positron effect",400,400);
//    PosSpectrumH->Draw("same");
//    a1->Print("../Images/PosResponseSpectrum.eps","png");
//    
//    TCanvas* a2 = new TCanvas("IAV effect","IAV effect",400,400);
//    IAVSpectrumH->Draw("same");
//    a2->Print("../Images/IAVResponseSpectrum.eps", "png");
//    
//    TCanvas* a3 = new TCanvas("NL effect","NL effect",400,400);
//    NLSpectrumH->Draw("same");
//    a3->Print("../Images/NLResponseSpectrum.eps", "png");
//    
//    TCanvas* a4 = new TCanvas("Reso effect","Reso effect",400,400);
//    ResoSpectrumH->Draw("same");
//    a4->Print("../Images/ResoResponseSpectrum.eps", "png");

    TCanvas* a13 = new TCanvas("Comparison neutrino - pos","Comparison neutrino - pos",400,400);
    PosSpectrumH->Draw("same");
    NominalTrueSpectrumH->Draw("same");
    
    TCanvas* a14 = new TCanvas("Comparison Pos - IAV","Comparison Pos - IAV",400,400);
    PosSpectrumH->Draw();
    IAVSpectrumH->Draw("same");
    
    TCanvas* a15 = new TCanvas("Comparison IAV - NL","Comparison IAV - NL",400,400);
    IAVSpectrumH->Draw();
    NLSpectrumH->Draw("same");
    
    TCanvas* a16 = new TCanvas("Comparison NL - Reso","Comparison NL - Reso",400,400);
    NLSpectrumH->Draw();
    ResoSpectrumH->Draw("same");
    
    return 0;
    
}