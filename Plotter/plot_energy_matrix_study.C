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
#include "TGraph.h"
#include "TH1D.h"
#include "TH2.h"

int plot_energy_matrix_study(void)
{
    enum Type{FineColumnFlat, FineColumnSpectrum, FineRowFlat,FineRowSpectrum,RebinColumnFlat,RebinColumnSpectrum,RebinRowFlat,RebinRowSpectrum,NoMatrixFlat,NoMatrixSpectrum};//Add here new NL models
    Type TypeE;
    
    Int_t far = 0;
    Int_t near = 0;
    TH1::AddDirectory(kFALSE);
    TH1D* PredictionH[10];
    TH1D* PredictionRatioH[50];
    gStyle->SetErrorX(0.0001);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.05,"x");
    gStyle->SetTitleSize(0.05,"y");
    
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
    //////////// Energy Matrix Studies
    //////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// ///////////
    
    TFile* PlotF2 = TFile::Open("../RootOutputs/Spectra/EnergyMatrixSpectrumStudy/Row0_Rebin0_Flat1.root");
    TCanvas *c2 = new TCanvas("FineColumnFlat Spectrum","FineColumnFlat Spectrum",300,300);
    
    PredictionH[FineColumnFlat]=(TH1D*)gDirectory->Get(Form("Vis Prediction AD%i from AD%i", far+1, near+1));
    PredictionH[FineColumnFlat]->SetTitle("Fine/Column/Flat");
    PredictionH[FineColumnFlat]->GetXaxis()->SetTitle("E_{vis} (MeV)");
    PredictionH[FineColumnFlat]->GetXaxis()->SetTitleSize(0.05);
    PredictionH[FineColumnFlat]->GetYaxis()->SetTitleSize(0.05);
    PredictionH[FineColumnFlat]->Draw();
    c2->Print("../Images/EnergyMatrixStudy/FineColumnFlat.eps", "png");
    PlotF2->Close();
    
    TFile* PlotF3 = TFile::Open("../RootOutputs/Spectra/EnergyMatrixSpectrumStudy/Row0_Rebin0_Flat0.root");
    TCanvas *c3 = new TCanvas("FineColumnSpectrum","FineColumnSpectrum",300,300);
    
    PredictionH[FineColumnSpectrum]=(TH1D*)gDirectory->Get(Form("Vis Prediction AD%i from AD%i", far+1, near+1));
    PredictionH[FineColumnSpectrum]->SetTitle("FineColumnSpectrum");
    PredictionH[FineColumnSpectrum]->GetXaxis()->SetTitle("E_{vis} (MeV)");
    PredictionH[FineColumnSpectrum]->GetXaxis()->SetTitleSize(0.05);
    PredictionH[FineColumnSpectrum]->GetYaxis()->SetTitleSize(0.05);
    PredictionH[FineColumnSpectrum]->Draw();
    c3->Print("../Images/EnergyMatrixStudy/FineColumnSpectrum.eps", "png");
    PlotF3->Close();

    TFile* PlotF4 = TFile::Open("../RootOutputs/Spectra/EnergyMatrixSpectrumStudy/Row1_Rebin0_Flat1.root");
    TCanvas *c4 = new TCanvas("FineRowFlat Spectrum","FineRowFlat Spectrum",300,300);
    
    PredictionH[FineRowFlat]=(TH1D*)gDirectory->Get(Form("Vis Prediction AD%i from AD%i", far+1, near+1));
    PredictionH[FineRowFlat]->SetTitle("Fine/Row/Flat");
    PredictionH[FineRowFlat]->GetXaxis()->SetTitle("E_{vis} (MeV)");
    PredictionH[FineRowFlat]->GetXaxis()->SetTitleSize(0.05);
    PredictionH[FineRowFlat]->GetYaxis()->SetTitleSize(0.05);
    PredictionH[FineRowFlat]->Draw();
    c4->Print("../Images/EnergyMatrixStudy/FineRowFlat.eps", "png");
    PlotF4->Close();

    TFile* PlotF5 = TFile::Open("../RootOutputs/Spectra/EnergyMatrixSpectrumStudy/Row1_Rebin0_Flat0.root");
    TCanvas *c5 = new TCanvas("FineRowSpectrum","FineRowSpectrum",300,300);
    
    PredictionH[FineRowSpectrum]=(TH1D*)gDirectory->Get(Form("Vis Prediction AD%i from AD%i", far+1, near+1));
    PredictionH[FineRowSpectrum]->SetTitle("Fine/Row/Spectrum");
    PredictionH[FineRowSpectrum]->GetXaxis()->SetTitle("E_{vis} (MeV)");
    PredictionH[FineRowSpectrum]->GetXaxis()->SetTitleSize(0.05);
    PredictionH[FineRowSpectrum]->GetYaxis()->SetTitleSize(0.05);
    PredictionH[FineRowSpectrum]->Draw();
    c5->Print("../Images/EnergyMatrixStudy/FineRowSpectrum.eps", "png");
    PlotF5->Close();

    TFile* PlotF6 = TFile::Open("../RootOutputs/Spectra/EnergyMatrixSpectrumStudy/Row0_Rebin1_Flat1.root");
    TCanvas *c6 = new TCanvas("RebinColumnFlat Spectrum","RebinColumnFlat Spectrum",300,300);
    
    PredictionH[RebinColumnFlat]=(TH1D*)gDirectory->Get(Form("Vis Prediction AD%i from AD%i", far+1, near+1));
    PredictionH[RebinColumnFlat]->SetTitle("Rebin/Column/Flat");
    PredictionH[RebinColumnFlat]->GetXaxis()->SetTitle("E_{vis} (MeV)");
    PredictionH[RebinColumnFlat]->GetXaxis()->SetTitleSize(0.05);
    PredictionH[RebinColumnFlat]->GetYaxis()->SetTitleSize(0.05);
    PredictionH[RebinColumnFlat]->Draw();
    c6->Print("../Images/EnergyMatrixStudy/RebinColumnFlat.eps", "png");
    PlotF6->Close();

    TFile* PlotF8 = TFile::Open("../RootOutputs/Spectra/EnergyMatrixSpectrumStudy/Row0_Rebin1_Flat0.root");
    TCanvas *c8 = new TCanvas("RebinColumnSpectrum","RebinColumnSpectrum",300,300);
    
    PredictionH[RebinColumnSpectrum]=(TH1D*)gDirectory->Get(Form("Vis Prediction AD%i from AD%i", far+1, near+1));
    PredictionH[RebinColumnSpectrum]->SetTitle("Rebin/Column/Spectrum");
    PredictionH[RebinColumnSpectrum]->GetXaxis()->SetTitle("E_{vis} (MeV)");
    PredictionH[RebinColumnSpectrum]->GetXaxis()->SetTitleSize(0.05);
    PredictionH[RebinColumnSpectrum]->GetYaxis()->SetTitleSize(0.05);
    PredictionH[RebinColumnSpectrum]->Draw();
    c8->Print("../Images/EnergyMatrixStudy/RebinColumnSpectrum.eps", "png");
    PlotF8->Close();

    TFile* PlotF9 = TFile::Open("../RootOutputs/Spectra/EnergyMatrixSpectrumStudy/Row1_Rebin1_Flat1.root");
    TCanvas *c9 = new TCanvas("RebinRowFlat Spectrum","RebinRowFlat Spectrum",300,300);
    
    PredictionH[RebinRowFlat]=(TH1D*)gDirectory->Get(Form("Vis Prediction AD%i from AD%i", far+1, near+1));
    PredictionH[RebinRowFlat]->SetTitle("Rebin/Row/Flat");
    PredictionH[RebinRowFlat]->GetXaxis()->SetTitle("E_{vis} (MeV)");
    PredictionH[RebinRowFlat]->GetXaxis()->SetTitleSize(0.05);
    PredictionH[RebinRowFlat]->GetYaxis()->SetTitleSize(0.05);
    PredictionH[RebinRowFlat]->Draw();
    c9->Print("../Images/EnergyMatrixStudy/RebinRowFlat.eps", "png");
    PlotF9->Close();

    TFile* PlotF12 = TFile::Open("../RootOutputs/Spectra/EnergyMatrixSpectrumStudy/Row1_Rebin1_Flat0.root");
    TCanvas *c12 = new TCanvas("RebinRowSpectrum","RebinRowSpectrum",300,300);
    
    PredictionH[RebinRowSpectrum]=(TH1D*)gDirectory->Get(Form("Vis Prediction AD%i from AD%i", far+1, near+1));
    PredictionH[RebinRowSpectrum]->SetTitle("Rebin/Row/Spectrum");
    PredictionH[RebinRowSpectrum]->GetXaxis()->SetTitle("E_{vis} (MeV)");
    PredictionH[RebinRowSpectrum]->GetXaxis()->SetTitleSize(0.05);
    PredictionH[RebinRowSpectrum]->GetYaxis()->SetTitleSize(0.05);
    PredictionH[RebinRowSpectrum]->Draw();
    c12->Print("../Images/EnergyMatrixStudy/RebinRowSpectrum.eps", "png");
    PlotF12->Close();

    TFile* PlotF10 = TFile::Open("../RootOutputs/Spectra/EnergyMatrixSpectrumStudy/NoMatrixSpectrum.root");
    TCanvas *c10 = new TCanvas("NoMatrixSpectrum","NoMatrixSpectrum",300,300);
    
    PredictionH[NoMatrixSpectrum]=(TH1D*)gDirectory->Get(Form("Vis Prediction AD%i from AD%i", far+1, near+1));
    PredictionH[NoMatrixSpectrum]->SetTitle("NoMatrixSpectrum");
    PredictionH[NoMatrixSpectrum]->GetXaxis()->SetTitle("E_{vis} (MeV)");
    PredictionH[NoMatrixSpectrum]->GetXaxis()->SetTitleSize(0.05);
    PredictionH[NoMatrixSpectrum]->GetYaxis()->SetTitleSize(0.05);
    PredictionH[NoMatrixSpectrum]->Draw();
    c10->Print("../Images/EnergyMatrixStudy/NoMatrixSpectrum.eps", "png");
    PlotF10->Close();

    TFile* PlotF11 = TFile::Open("../RootOutputs/Spectra/EnergyMatrixSpectrumStudy/NoMatrixFlat.root");
    TCanvas *c11 = new TCanvas("NoMatrixFlat","NoMatrixFlat",300,300);
    
    PredictionH[NoMatrixFlat]=(TH1D*)gDirectory->Get(Form("Vis Prediction AD%i from AD%i", far+1, near+1));
    PredictionH[NoMatrixFlat]->SetTitle("NoMatrixFlat");
    PredictionH[NoMatrixFlat]->GetXaxis()->SetTitle("E_{vis} (MeV)");
    PredictionH[NoMatrixFlat]->GetXaxis()->SetTitleSize(0.05);
    PredictionH[NoMatrixFlat]->GetYaxis()->SetTitleSize(0.05);
    PredictionH[NoMatrixFlat]->Draw();
    c11->Print("../Images/EnergyMatrixStudy/NoMatrixFlat.eps", "png");
    PlotF11->Close();
    
    //  Ratios
    TCanvas *c13 = new TCanvas("Flat Spectrum Ratio","Flat Spectrum Ratio",600,600);
    c13->Divide(2,2);
    c13->cd(1);
         PredictionRatioH[0]=(TH1D*)PredictionH[FineColumnFlat]->Clone("FineColumn(FlatVsSpectrum)");
         PredictionRatioH[0]->SetTitle("FineColumn(FlatVsSpectrum)");
         PredictionRatioH[0]->Add(PredictionH[FineColumnSpectrum],-1);
         PredictionRatioH[0]->Divide(PredictionH[FineColumnFlat]);
         PredictionRatioH[0]->GetYaxis()->SetTitle("(Flat-Spectrum)/Flat");
         PredictionRatioH[0]->Draw();

    c13->cd(2);
        PredictionRatioH[1]=(TH1D*)PredictionH[FineRowFlat]->Clone("FineRow(FlatVsSpectrum)");
        PredictionRatioH[1]->SetTitle("FineRow(FlatVsSpectrum)");
        PredictionRatioH[1]->Add(PredictionH[FineRowSpectrum],-1);
        PredictionRatioH[1]->Divide(PredictionH[FineRowFlat]);
        PredictionRatioH[1]->GetYaxis()->SetTitle("(Flat-Spectrum)/Flat");
        PredictionRatioH[1]->Draw();
    
    c13->cd(3);
        PredictionRatioH[2]=(TH1D*)PredictionH[RebinColumnFlat]->Clone("RebinColumn(FlatVsSpectrum)");
        PredictionRatioH[2]->SetTitle("RebinColumn(FlatVsSpectrum)");
        PredictionRatioH[2]->Add(PredictionH[RebinRowSpectrum],-1);
        PredictionRatioH[2]->Divide(PredictionH[RebinColumnFlat]);
        PredictionRatioH[2]->GetYaxis()->SetTitle("(Flat-Spectrum)/Flat");
        PredictionRatioH[2]->Draw();
    
    c13->cd(4);
        PredictionRatioH[3]=(TH1D*)PredictionH[RebinRowFlat]->Clone("RebinRow(FlatVsSpectrum)");
        PredictionRatioH[3]->SetTitle("RebinRow(FlatVsSpectrum)");
        PredictionRatioH[3]->Add(PredictionH[RebinRowSpectrum],-1);
        PredictionRatioH[3]->Divide(PredictionH[RebinRowFlat]);
        PredictionRatioH[3]->GetYaxis()->SetTitle("(Flat-Spectrum)/Flat");
        PredictionRatioH[3]->Draw();
    
    c13->Print("../Images/FlatVsSpectrum.eps", "png");

    TCanvas *c14 = new TCanvas("Rebin Fine Ratio","Rebin Fine Ratio",600,600);
    c14->Divide(2,2);
    c14->cd(1);
    PredictionRatioH[4]=(TH1D*)PredictionH[FineColumnFlat]->Clone("FlatColumn(FineVsRebin)");
    PredictionRatioH[4]->SetTitle("FlatColumn(FineVsRebin)");
    PredictionRatioH[4]->Add(PredictionH[RebinColumnFlat],-1);
    PredictionRatioH[4]->Divide(PredictionH[FineColumnFlat]);
    PredictionRatioH[4]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[4]->Draw();
    
    c14->cd(2);
    PredictionRatioH[5]=(TH1D*)PredictionH[FineRowFlat]->Clone("FlatRow(FineVsRebin)");
    PredictionRatioH[5]->SetTitle("FlatRow(FineVsRebin)");
    PredictionRatioH[5]->Add(PredictionH[RebinRowFlat],-1);
    PredictionRatioH[5]->Divide(PredictionH[FineRowFlat]);
    PredictionRatioH[5]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[5]->Draw();
    
    c14->cd(3);
    PredictionRatioH[6]=(TH1D*)PredictionH[FineColumnSpectrum]->Clone("SpectrumColumn(FineVsRebin)");
    PredictionRatioH[6]->SetTitle("SpectrumColumn(FineVsRebin)");
    PredictionRatioH[6]->Add(PredictionH[RebinColumnSpectrum],-1);
    PredictionRatioH[6]->Divide(PredictionH[FineColumnSpectrum]);
    PredictionRatioH[6]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[6]->Draw();
    
    c14->cd(4);
    PredictionRatioH[7]=(TH1D*)PredictionH[FineRowSpectrum]->Clone("SpectrumRow(FineVsRebin)");
    PredictionRatioH[7]->SetTitle("SpectrumRow(FineVsRebin)");
    PredictionRatioH[7]->Add(PredictionH[RebinRowSpectrum],-1);
    PredictionRatioH[7]->Divide(PredictionH[FineRowSpectrum]);
    PredictionRatioH[7]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[7]->Draw();
    
    c14->Print("../Images/FineVsRebin.eps", "png");
    
    TCanvas *c15 = new TCanvas("Column Row Ratio","Columnd Row Ratio",600,600);
    c15->Divide(2,2);
    c15->cd(1);
    PredictionRatioH[4]=(TH1D*)PredictionH[FineColumnFlat]->Clone("FlatFine(ColumnVsRow)");
    PredictionRatioH[4]->SetTitle("FlatFine(ColumnVsRow)");
    PredictionRatioH[4]->Add(PredictionH[FineRowFlat],-1);
    PredictionRatioH[4]->Divide(PredictionH[FineColumnFlat]);
    PredictionRatioH[4]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[4]->Draw();
    
    c15->cd(2);
    PredictionRatioH[5]=(TH1D*)PredictionH[RebinColumnFlat]->Clone("FlatRebin(ColumnVsRow)");
    PredictionRatioH[5]->SetTitle("FlatRebin(ColumnVsRow)");
    PredictionRatioH[5]->Add(PredictionH[RebinRowFlat],-1);
    PredictionRatioH[5]->Divide(PredictionH[RebinColumnFlat]);
    PredictionRatioH[5]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[5]->Draw();
    
    c15->cd(3);
    PredictionRatioH[6]=(TH1D*)PredictionH[FineColumnSpectrum]->Clone("SpectrumFine(ColumnVsRow)");
    PredictionRatioH[6]->SetTitle("SpectrumFine(ColumnVsRow)");
    PredictionRatioH[6]->Add(PredictionH[FineRowSpectrum],-1);
    PredictionRatioH[6]->Divide(PredictionH[FineColumnSpectrum]);
    PredictionRatioH[6]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[6]->Draw();
    
    c15->cd(4);
    PredictionRatioH[7]=(TH1D*)PredictionH[RebinColumnSpectrum]->Clone("SpectrumRebin(ColumnVsRow)");
    PredictionRatioH[7]->SetTitle("SpectrumRebin(ColumnVsRow)");
    PredictionRatioH[7]->Add(PredictionH[RebinRowSpectrum],-1);
    PredictionRatioH[7]->Divide(PredictionH[RebinColumnSpectrum]);
    PredictionRatioH[7]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[7]->Draw();
    
    c15->Print("../Images/ColumnVsRow.eps", "png");
//    THIS IS USELESS:
//    TCanvas *c16 = new TCanvas("Matrix Vs No Matrix Flat Ratio","Matrix Vs No Matrix Flat Ratio",600,1200);
//    c16->Divide(2,4);
//    c16->cd(1);
//    PredictionRatioH[8]=(TH1D*)PredictionH[NoMatrixFlat]->Clone("FlatFineRow(MatrixVsFlatNoMatrix)");
//    PredictionRatioH[8]->SetTitle("FlatFineRow(MatrixVsNoMatrix)");
//    PredictionRatioH[8]->Add(PredictionH[FineRowFlat],-1);
//    PredictionRatioH[8]->Divide(PredictionH[NoMatrixFlat]);
//    PredictionRatioH[8]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
//    PredictionRatioH[8]->Draw();
//    
//    c16->cd(2);
//    PredictionRatioH[9]=(TH1D*)PredictionH[NoMatrixFlat]->Clone("FlatRebinRow(MatrixVsFlatNoMatrix)");
//    PredictionRatioH[9]->SetTitle("FlatRebinRow(MatrixVsNoMatrix)");
//    PredictionRatioH[9]->Add(PredictionH[RebinRowFlat],-1);
//    PredictionRatioH[9]->Divide(PredictionH[NoMatrixFlat]);
//    PredictionRatioH[9]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
//    PredictionRatioH[9]->Draw();
//    
//    c16->cd(3);
//    PredictionRatioH[10]=(TH1D*)PredictionH[NoMatrixFlat]->Clone("SpectrumFineRow(MatrixVsFlatNoMatrix)");
//    PredictionRatioH[10]->SetTitle("SpectrumFineRow(MatrixVsNoMatrix)");
//    PredictionRatioH[10]->Add(PredictionH[FineRowSpectrum],-1);
//    PredictionRatioH[10]->Divide(PredictionH[NoMatrixFlat]);
//    PredictionRatioH[10]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
//    PredictionRatioH[10]->Draw();
//    
//    c16->cd(4);
//    PredictionRatioH[11]=(TH1D*)PredictionH[NoMatrixFlat]->Clone("SpectrumRebinRow(MatrixVsFlatNoMatrix)");
//    PredictionRatioH[11]->SetTitle("SpectrumRebinRow(MatrixVsNoMatrix)");
//    PredictionRatioH[11]->Add(PredictionH[RebinRowSpectrum],-1);
//    PredictionRatioH[11]->Divide(PredictionH[NoMatrixFlat]);
//    PredictionRatioH[11]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
//    PredictionRatioH[11]->Draw();
//    
//    c16->cd(5);
//    PredictionRatioH[12]=(TH1D*)PredictionH[NoMatrixFlat]->Clone("FlatFineColumn(MatrixVsFlatNoMatrix)");
//    PredictionRatioH[12]->SetTitle("FlatFineColumn(MatrixVsNoMatrix)");
//    PredictionRatioH[12]->Add(PredictionH[FineColumnFlat],-1);
//    PredictionRatioH[12]->Divide(PredictionH[NoMatrixFlat]);
//    PredictionRatioH[12]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
//    PredictionRatioH[12]->Draw();
//    
//    c16->cd(6);
//    PredictionRatioH[13]=(TH1D*)PredictionH[NoMatrixFlat]->Clone("FlatRebinColumn(MatrixVsFlatNoMatrix)");
//    PredictionRatioH[13]->SetTitle("FlatRebinColumn(MatrixVsNoMatrix)");
//    PredictionRatioH[13]->Add(PredictionH[RebinColumnFlat],-1);
//    PredictionRatioH[13]->Divide(PredictionH[NoMatrixFlat]);
//    PredictionRatioH[13]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
//    PredictionRatioH[13]->Draw();
//    
//    c16->cd(7);
//    PredictionRatioH[14]=(TH1D*)PredictionH[NoMatrixFlat]->Clone("SpectrumFineColumn(MatrixVsFlatNoMatrix)");
//    PredictionRatioH[14]->SetTitle("SpectrumFineColumn(MatrixVsNoMatrix)");
//    PredictionRatioH[14]->Add(PredictionH[FineColumnSpectrum],-1);
//    PredictionRatioH[14]->Divide(PredictionH[NoMatrixFlat]);
//    PredictionRatioH[14]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
//    PredictionRatioH[14]->Draw();
//    
//    c16->cd(8);
//    PredictionRatioH[15]=(TH1D*)PredictionH[NoMatrixFlat]->Clone("SpectrumRebinColumn(MatrixVsFlatNoMatrix)");
//    PredictionRatioH[15]->SetTitle("SpectrumRebinColumn(MatrixVsNoMatrix)");
//    PredictionRatioH[15]->Add(PredictionH[RebinColumnSpectrum],-1);
//    PredictionRatioH[15]->Divide(PredictionH[NoMatrixFlat]);
//    PredictionRatioH[15]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
//    PredictionRatioH[15]->Draw();
//    
//    c16->Print("../Images/FlatNoMatrixVsAll.eps", "png");
//    
    
    TCanvas *c17 = new TCanvas("Matrix Vs No Matrix Spectrum Ratio","Matrix Vs No Matrix Spectrum Ratio",1200,600);
    c17->Divide(4,2);
    c17->cd(1);
    PredictionRatioH[16]=(TH1D*)PredictionH[NoMatrixSpectrum]->Clone("FlatFineRow(MatrixVsSpectrumNoMatrix)");
    PredictionRatioH[16]->SetTitle("FlatFineRow(MatrixVsNoMatrix)");
    PredictionRatioH[16]->Add(PredictionH[FineRowFlat],-1);
    PredictionRatioH[16]->Divide(PredictionH[NoMatrixSpectrum]);
    PredictionRatioH[16]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[16]->Draw();
    
    c17->cd(2);
    PredictionRatioH[17]=(TH1D*)PredictionH[NoMatrixSpectrum]->Clone("FlatRebinRow(MatrixVsSpectrumNoMatrix)");
    PredictionRatioH[17]->SetTitle("FlatRebinRow(MatrixVsNoMatrix)");
    PredictionRatioH[17]->Add(PredictionH[RebinRowFlat],-1);
    PredictionRatioH[17]->Divide(PredictionH[NoMatrixSpectrum]);
    PredictionRatioH[17]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[17]->Draw();
    
    c17->cd(3);
    PredictionRatioH[18]=(TH1D*)PredictionH[NoMatrixSpectrum]->Clone("SpectrumFineRow(MatrixVsSpectrumNoMatrix)");
    PredictionRatioH[18]->SetTitle("SpectrumFineRow(MatrixVsNoMatrix)");
    PredictionRatioH[18]->Add(PredictionH[FineRowSpectrum],-1);
    PredictionRatioH[18]->Divide(PredictionH[NoMatrixSpectrum]);
    PredictionRatioH[18]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[18]->Draw();
    
    c17->cd(4);
    PredictionRatioH[19]=(TH1D*)PredictionH[NoMatrixSpectrum]->Clone("SpectrumRebinRow(MatrixVsSpectrumNoMatrix)");
    PredictionRatioH[19]->SetTitle("SpectrumRebinRow(MatrixVsNoMatrix)");
    PredictionRatioH[19]->Add(PredictionH[RebinRowSpectrum],-1);
    PredictionRatioH[19]->Divide(PredictionH[NoMatrixSpectrum]);
    PredictionRatioH[19]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[19]->Draw();
    
    c17->cd(5);
    PredictionRatioH[20]=(TH1D*)PredictionH[NoMatrixSpectrum]->Clone("FlatFineColumn(MatrixVsSpectrumNoMatrix)");
    PredictionRatioH[20]->SetTitle("FlatFineColumn(MatrixVsNoMatrix)");
    PredictionRatioH[20]->Add(PredictionH[FineColumnFlat],-1);
    PredictionRatioH[20]->Divide(PredictionH[NoMatrixSpectrum]);
    PredictionRatioH[20]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[20]->Draw();
    
    c17->cd(6);
    PredictionRatioH[21]=(TH1D*)PredictionH[NoMatrixSpectrum]->Clone("FlatRebinColumn(MatrixVsSpectrumNoMatrix)");
    PredictionRatioH[21]->SetTitle("FlatRebinColumn(MatrixVsNoMatrix)");
    PredictionRatioH[21]->Add(PredictionH[RebinColumnFlat],-1);
    PredictionRatioH[21]->Divide(PredictionH[NoMatrixSpectrum]);
    PredictionRatioH[21]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[21]->Draw();
    
    c17->cd(7);
    PredictionRatioH[22]=(TH1D*)PredictionH[NoMatrixSpectrum]->Clone("SpectrumFineColumn(MatrixVsSpectrumNoMatrix)");
    PredictionRatioH[22]->SetTitle("SpectrumFineColumn(MatrixVsNoMatrix)");
    PredictionRatioH[22]->Add(PredictionH[FineColumnSpectrum],-1);
    PredictionRatioH[22]->Divide(PredictionH[NoMatrixSpectrum]);
    PredictionRatioH[22]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[22]->Draw();
    
    c17->cd(8);
    PredictionRatioH[23]=(TH1D*)PredictionH[NoMatrixSpectrum]->Clone("SpectrumRebinColumn(MatrixVsSpectrumNoMatrix)");
    PredictionRatioH[23]->SetTitle("SpectrumRebinColumn(MatrixVsNoMatrix)");
    PredictionRatioH[23]->Add(PredictionH[RebinColumnSpectrum],-1);
    PredictionRatioH[23]->Divide(PredictionH[NoMatrixSpectrum]);
    PredictionRatioH[23]->GetYaxis()->SetTitle("(Fine-Rebin)/Fine");
    PredictionRatioH[23]->Draw();
    
    c17->Print("../Images/SpectrumNoMatrixVsAll.eps", "png");
//    //  Ratios
//    TCanvas *c14 = new TCanvas("Prediction Spectrum Ratio","Prediction Spectrum Ratio",1200,1200);
//    c14->Divide(4,4);
//    c14->cd(1);
//    
//    PredictionRatioH[5]=(TH1D*)PredictionRatioH[0]->Clone("Column(FineVsRebin)(FlatVsSpectrum/FlatVsSpectrum)");
//    PredictionRatioH[5]->SetTitle("Column(FineVsRebin)");
//    PredictionRatioH[5]->Add(PredictionRatioH[2],-1);
//    PredictionRatioH[5]->Divide(PredictionRatioH[0]);
//    PredictionRatioH[5]->GetYaxis()->SetTitle("(FineColumn(Flat-Spectrum)/Flat - RebinColumn(Flat-Spectrum)/Flat) / FineColumn(Flat-Spectrum)/Flat");
//    PredictionRatioH[5]->Draw();
//
//    
//    PredictionRatioH[4]=(TH1D*)PredictionRatioH[1]->Clone("Row(FineVsRebin)(FlatVsSpectrum/FlatVsSpectrum)");
//    PredictionRatioH[4]->SetTitle("Row(FineVsRebin)");
//    PredictionRatioH[4]->Add(PredictionRatioH[3],-1);
//    PredictionRatioH[4]->Divide(PredictionRatioH[1]);
//    PredictionRatioH[4]->GetYaxis()->SetTitle("(FineRow(Flat-Spectrum)/Flat - RebinRow(Flat-Spectrum)/Flat) / FineRow(Flat-Spectrum)/Flat");
//    PredictionRatioH[4]->Draw();
//
//
//
//
//    PredictionRatioH[6]=(TH1D*)PredictionRatioH[5]->Clone("FlatVsSpectrumColumn");
//    PredictionRatioH[6]->SetTitle("FlatVsSpectrumColumn");
//    PredictionRatioH[6]->Add(PredictionRatioH[4],-1);
//    PredictionRatioH[6]->Divide(PredictionRatioH[5]);
//    PredictionRatioH[6]->GetYaxis()->SetTitle("Flat(NoMatrix-ColumnMatrix)/NoMatrix - Spectrum(NoMatrix-ColumnMatrix)/NoMatrix / Flat(NoMatrix-ColumnMatrix)/NoMatrix");
//    PredictionRatioH[6]->Draw();
//    
//    c14->cd(2);
//    
//    PredictionRatioH[7]=(TH1D*)PredictionH[FineColumnFlat]->Clone("FineFlat(ColumndvsRow)");
//    PredictionRatioH[7]->SetTitle("FineFlat(ColumndvsRow)");
//    PredictionRatioH[7]->Add(PredictionH[FineRowFlat],-1);
//    PredictionRatioH[7]->Divide(PredictionH[FineColumnFlat]);
//    PredictionRatioH[7]->GetYaxis()->SetTitle("(Column-Row)/Column");
//    PredictionRatioH[7]->Draw();
//    
//    c14->cd(3);
//    
//    PredictionRatioH[8]=(TH1D*)PredictionH[FineRowFlat]->Clone("FlatRow(FinevsRebin)");
//    PredictionRatioH[8]->SetTitle("FlatRow(FinevsRebin)");
//    PredictionRatioH[8]->Add(PredictionH[RebinRowFlat],-1);
//    PredictionRatioH[8]->Divide(PredictionH[FineRowFlat]);
//    PredictionRatioH[8]->GetYaxis()->SetTitle("Flat(Fine-Rebin)/Fine");
//    PredictionRatioH[8]->Draw();
//    
//    c14->cd(4);
//    
//    PredictionRatioH[9]=(TH1D*)PredictionH[FineRowSpectrum]->Clone("SpectrumRow(FinevsRebin)");
//    PredictionRatioH[9]->SetTitle("SpectrumRow(FinevsRebin)");
//    PredictionRatioH[9]->Add(PredictionH[RebinRowSpectrum],-1);
//    PredictionRatioH[9]->Divide(PredictionH[FineRowSpectrum]);
//    PredictionRatioH[9]->GetYaxis()->SetTitle("Spectrum(Fine-Rebin)/Fine");
//    PredictionRatioH[9]->Draw();
//    
//    
//    c14->cd(5);
//    
//    PredictionRatioH[10]=(TH1D*)PredictionH[NoMatrixFlat]->Clone("FlatNoMatrixVsColumnMatrix");
//    PredictionRatioH[10]->SetTitle("FlatNoMatrixVsColumnMatrix");
//    PredictionRatioH[10]->Add(PredictionH[FineColumnFlat],-1);
//    PredictionRatioH[10]->Divide(PredictionH[NoMatrixFlat]);
//    PredictionRatioH[10]->GetYaxis()->SetTitle("Spectrum(NoMatrix-Matrix)/NoColumnMatrix");
//    PredictionRatioH[10]->Draw();
//    
//    
//    c14->cd(6);
//    
//    PredictionRatioH[11]=(TH1D*)PredictionH[NoMatrixSpectrum]->Clone("SpectrumNoMatrixVsColumnMatrix");
//    PredictionRatioH[11]->SetTitle("SpectrumNoMatrixVsColumnMatrix");
//    PredictionRatioH[11]->Add(PredictionH[FineColumnSpectrum],-1);
//    PredictionRatioH[11]->Divide(PredictionH[NoMatrixSpectrum]);
//    PredictionRatioH[11]->GetYaxis()->SetTitle("Spectrum(NoMatrix-ColumnMatrix)/NoMatrix");
//    PredictionRatioH[11]->Draw();
//    
//    c14->cd(7);
//    
//    PredictionRatioH[12]=(TH1D*)PredictionH[NoMatrixFlat]->Clone("FlatNoMatrixVsRowMatrix");
//    PredictionRatioH[12]->SetTitle("FlatNoMatrixVsRowMatrix");
//    PredictionRatioH[12]->Add(PredictionH[FineRowFlat],-1);
//    PredictionRatioH[12]->Divide(PredictionH[NoMatrixFlat]);
//    PredictionRatioH[12]->GetYaxis()->SetTitle("Flat(NoMatrix-RowMatrix)/NoMatrix");
//    PredictionRatioH[12]->Draw();
//    
//    
//    c14->cd(8);
//    
//    PredictionRatioH[13]=(TH1D*)PredictionH[FineRowSpectrum]->Clone("SpectrumNoMatrixVsRowMatrix");
//    PredictionRatioH[13]->SetTitle("SpectrumNoMatrixVsRowMatrix");
//    PredictionRatioH[13]->Add(PredictionH[NoMatrixSpectrum],-1);
//    PredictionRatioH[13]->Divide(PredictionH[FineRowSpectrum]);
//    PredictionRatioH[13]->GetYaxis()->SetTitle("Spectrum(NoMatrix-RowMatrix)/NoMatrix");
//    PredictionRatioH[13]->Draw();
//    
//    c14->Print("../Images/FlatVsSpectrum.eps", "png");
}
