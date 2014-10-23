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
#include "TGraph.h"
#include "TH1D.h"
#include "TH2.h"

int plot_spectrum(void)
{
    TH1::AddDirectory(kFALSE);
    
    TH1D* AntineutrinoH[6][1];
    TH1D* PredictionH[3][3][1];
    TH1D* ReactorPredictionH[6][1];
    TH1D* OscPredictionH[6][1];
    TH1D* SpectrumRatio[6][1];
    
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
    //////////// Factor Studies
    //////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// /////////// ///////////
    
    Char_t name[20];
    TFile* PlotF1 = TFile::Open("../RootOutputs/Spectra/AntineutrinoSpectrum.root");
    
    TCanvas *c1 = new TCanvas("Antineutrino Spectrum","Antineutrino Spectrum",1200,600);
    c1->Divide(3,2);

    for(Int_t week=0;week<1;++week)
    {
            for(Int_t reactor=0;reactor<6;reactor++)
            {
                c1->cd(reactor+1);

                AntineutrinoH[reactor][week]=(TH1D*)gDirectory->Get(Form("AntineutrinoSpectrumFromReactor%i", reactor+1));
                AntineutrinoH[reactor][week]->SetTitle("Antineutrino Spectrum");
                
                AntineutrinoH[reactor][week]->GetXaxis()->SetTitle("E_{true} (MeV)");
                AntineutrinoH[reactor][week]->GetXaxis()->SetTitleSize(0.05);
                AntineutrinoH[reactor][week]->GetYaxis()->SetTitleSize(0.05);
                AntineutrinoH[reactor][week]->SetTitle("");
                
                sprintf(name,"AntineutrinoReactor%i",reactor+1);
                AntineutrinoH[reactor][week]->GetYaxis()->SetTitle(name);
                AntineutrinoH[reactor][week]->Draw();
            }
        
    }
    c1->Print("../Images/AntineutrinoSpectrum.eps", "png");
    
    PlotF1->Close();
    
    TFile* PlotF2 = TFile::Open("../RootOutputs/Spectra/PredictedSpectrum.root");
    
    TCanvas *c2 = new TCanvas("Prediction Spectrum","Prediction Spectrum",1000,1000);
    c2->Divide(3,3);
    Int_t cont=0;

    for(Int_t week=0;week<1;++week)
    {
        for(Int_t near=0;near<3;near++)
        {
            for(Int_t far=0;far<3;far++)
            {
                ++cont;
                c2->cd(cont);
                PredictionH[far][near][week]=(TH1D*)gDirectory->Get(Form("Vis Prediction AD%i from AD%i", far+1, near+1));
                PredictionH[far][near][week]->SetTitle("Prediction Spectrum");
                
                PredictionH[far][near][week]->GetXaxis()->SetTitle("E_{vis} (MeV)");
                PredictionH[far][near][week]->GetXaxis()->SetTitleSize(0.05);
                PredictionH[far][near][week]->GetYaxis()->SetTitleSize(0.05);
                PredictionH[far][near][week]->SetTitle("");
                
                sprintf(name,"PredictionNear%iFar%i",near+1,far+1);
                PredictionH[far][near][week]->GetYaxis()->SetTitle(name);
                PredictionH[far][near][week]->Draw();
            }
        }
        
    }
    c2->Print("../Images/PredictionNearFar.eps", "png");
    
    PlotF2->Close();
    
    TFile* PlotF3 = TFile::Open("../RootOutputs/NominalOutputs/Oscillation.root");
    
    TCanvas *c3 = new TCanvas("ReactorPrediction Spectrum","ReactorPrediction Spectrum",1200,600);
    c3->Divide(3,2);
    
    for(Int_t week=0;week<1;++week)
    {
        for(Int_t ad=0;ad<6;ad++)
        {
            c3->cd(ad+1);
            PlotF3->cd("Total AD Spectra after oscillation");
            ReactorPredictionH[ad][week]=(TH1D*)gDirectory->Get(Form("Total spectrum after oscillation at AD%i", ad+1));
            sprintf(name,"ReactorPredictionAD%i",ad+1);

            ReactorPredictionH[ad][week]->SetTitle("ReactorPrediction Spectrum");
            
            ReactorPredictionH[ad][week]->GetXaxis()->SetTitle("E_{true} (MeV)");
            ReactorPredictionH[ad][week]->GetXaxis()->SetTitleSize(0.05);
            ReactorPredictionH[ad][week]->GetYaxis()->SetTitleSize(0.05);
            ReactorPredictionH[ad][week]->SetTitle("");
            
            ReactorPredictionH[ad][week]->GetYaxis()->SetTitle(name);
            ReactorPredictionH[ad][week]->Draw();
            
        }
    }
    c3->Print("../Images/ReactorPredictionTrue.eps", "png");
    
    PlotF3->Close();
    
    TFile* PlotF4 = TFile::Open("../RootOutputs/Spectra/NearSpectrumFraction.root");
    
    TCanvas *c4 = new TCanvas("Oscillation Prediction Spectrum","Oscillation Prediction Spectrum",1200,600);
    c4->Divide(3,2);

    for(Int_t week=0;week<1;++week)
    {
        for(Int_t ad=0;ad<6;ad++)
        {
            c4->cd(ad+1);

            if(ad<3)
            {

                OscPredictionH[ad][week]=(TH1D*)gDirectory->Get(Form("AD%i Near Prediction", ad+1));

            }
            else
            {
                if(ad==3)
                {
                    PlotF4->Close();
                    TFile* PlotF5 = TFile::Open("../RootOutputs/Spectra/FarSpectrumFraction.root");
                }
                
                OscPredictionH[ad][week]=(TH1D*)gDirectory->Get(Form("AD%i Far Spectrum prediction from near AD%i", ad-3+1,1));
            }
            
            OscPredictionH[ad][week]->SetTitle("Prediction Spectrum");
            OscPredictionH[ad][week]->GetXaxis()->SetTitle("E_{true} (MeV)");
            OscPredictionH[ad][week]->GetXaxis()->SetTitleSize(0.05);
            OscPredictionH[ad][week]->GetYaxis()->SetTitleSize(0.05);
            OscPredictionH[ad][week]->SetTitle("");

            sprintf(name,"OscPredictionAD%i",ad+1,1);
            OscPredictionH[ad][week]->GetYaxis()->SetTitle(name);
            OscPredictionH[ad][week]->Draw();
            
        }
    }
    c4->Print("../Images/OscPredictionTrue.eps", "png");
    PlotF5->Close();
    
    //  Difference:
    
    TCanvas *c6 = new TCanvas("Prediction Spectrum Ratio","Prediction Spectrum Ratio",1200,600);
    c6->Divide(3,2);
    
    for(Int_t week=0;week<1;++week)
    {
        for(Int_t ad=0;ad<6;ad++)
        {
            c6->cd(ad+1);
            sprintf(name,"Spectrum Ratio AD%i",ad+1);

            ReactorPredictionH[ad][week]->SetTitle(name);
            SpectrumRatio[ad][week]=(TH1D*)ReactorPredictionH[ad][week]->Clone();
            SpectrumRatio[ad][week]->Add(OscPredictionH[ad][week],-1);
            SpectrumRatio[ad][week]->Divide(ReactorPredictionH[ad][week]);
            
            SpectrumRatio[ad][week]->GetYaxis()->SetTitle("(OscModel-RelativeOscModel)/OscModel");
            SpectrumRatio[ad][week]->Draw();
        }
    }
    c6->Print("../Images/OscillationModelRatios.eps", "png");
}
