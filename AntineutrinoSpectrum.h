#pragma once
#include "NominalData.h"
#include <math.h>

//#define UseZhesCrossSection
//
//  AntineutrinoSpectrum.h
//  I multiply the reactour spectrum by the neutrino cross section to get the True Energy Spectrum without oscillation.
//
//  Created by Jose de Arcos	 on 3/13/13.
//
//
//const Int_t NReactors=6;//Defined in NominalData
const Double_t conv_m2_per_cm2 = 1.0e-4;
const Double_t conv_s_per_day = 24*60*60;
const bool TxTFile=1;
const Int_t MaxSamples = 1020;
const bool FlatSpectrumCheck = 0;

class AntineutrinoSpectrum
{
private:
    bool Analysis;
    bool Mode;
    NominalData* Data;
    std::string AnalysisString;
    bool IsotopeMatrix;
    bool ReactorPowerMatrix;
    
    TH1D* TotalReactorSpectrumH[NReactors];
    TH1D* AntineutrinoSpectrumH[NReactors];
    TH1D* CrossSectionH;

    Double_t EnergyM[MaxSamples];
    Double_t CrossSection[MaxSamples];
    
    Double_t InitialEnergy;
    Double_t FinalVisibleEnergy;

//    Double_t DetectorProtons;

    NominalData* Nom;
    
//    void DrawAntineutrinoSpectrum(TH1D*, TH1D*, TH1D*);

    Char_t filenameAntineutrino[200];
    Char_t filenameReactorFile[200];
    void LoadCrossSectionFile();

public:
    
    AntineutrinoSpectrum();
    AntineutrinoSpectrum(NominalData*);
    void AntineutrinoSpectrumMain(bool);
};

AntineutrinoSpectrum :: AntineutrinoSpectrum()
{
    std::cout << " the antineutrino default constructor shouldn't be called" << std::endl;
    
    exit(EXIT_FAILURE);
    
    Nom = new NominalData(0,2);
//    DetectorProtons = Nom->GetDetectorProtons(0);//  AD1
    InitialEnergy = Nom->GetEmin();
    FinalVisibleEnergy = Nom->GetEVisMax();
    Analysis = Nom->GetAnalysis();//  Gd or H data
    if(Analysis)
    {
        AnalysisString = "Hydrogen";
    }
    else
    {
        AnalysisString = "Gadolinium";
    }
    IsotopeMatrix = Nom->GetIsotopeMatrix();
    ReactorPowerMatrix = Nom->GetReactorPowerMatrix();
}

AntineutrinoSpectrum :: AntineutrinoSpectrum(NominalData* Data)
{
//    DetectorProtons = Data->GetDetectorProtons(0);
    InitialEnergy = Data->GetEmin();
    FinalVisibleEnergy = Data->GetEVisMax();
    Analysis = Data->GetAnalysis();//  Gd or H data
    if(Analysis)
    {
        AnalysisString = "Hydrogen";
    }
    else
    {
        AnalysisString = "Gadolinium";
    }
    IsotopeMatrix = Data->GetIsotopeMatrix();
    ReactorPowerMatrix = Data->GetReactorPowerMatrix();
}

void AntineutrinoSpectrum :: LoadCrossSectionFile()
{
        //Read txt file with already calculated nGd cross section
        
            CrossSectionH = new TH1D("Cross Section","Cross Section", Nbins, InitialEnergy,FinalVisibleEnergy);
    
#ifndef UseZhesCrossSection
            std::string line;
            std::ifstream CrossF("./CrossSections/Xsec1_2011.dat");
            
            Int_t Samples=0;
            
            while(Samples<MaxSamples)
            {
                std::getline(CrossF,line);
                std::string firstchar = line.substr(0,1);
                
                if(firstchar=="#")
                {
                    continue;//<-- ignore lines with comments
                }
                
                CrossF >> EnergyM[Samples] >> CrossSection[Samples];
                CrossSection[Samples]*=conv_s_per_day*conv_m2_per_cm2;// We will do the calculations in terms of days and m.
                //            std::cout <<EnergyM[Samples] <<  CrossSection[Samples] << std::endl;
                
                Samples++;
            }
            
            //        std::cout << Samples << "Samples" << std::endl;
            
            //        Double_t BinWidth = ((FinalVisibleEnergy-InitialEnergy)/Nbins);
            //        std::cout << BinWidth << " BinWidth" << std::endl;
            
            for (Int_t i= 0; i< Samples; i++)
            {
                //            std::cout <<  (fmod(EnergyM[i],BinWidth)) << " (fmod(EnergyM[i],BinWidth)" << std::endl;
                
                if((fmod(i,(Samples/Nbins))==0))
                {
                    //                std::cout <<(i/(Samples/Nbins))+1<< "(INDEX)" << std::endl;
                    
                    CrossSectionH->SetBinContent((i/(Samples/Nbins))+1,CrossSection[i]);
                }
            }
#else
        TFile* CrossSectionF = new TFile("./CrossSections/nHCrossSection.root");
        CrossSectionH = (TH1D*)gDirectory->Get("CrossSection");
        CrossSectionH->Scale(10e-42*conv_s_per_day*conv_m2_per_cm2);// We will do the calculations in terms of the # days and m. //DetectorProtons
        delete CrossSectionF;
#endif
    
//    TFile* CrossSecF = new TFile("./CrossSections/CheckCrossSection.root","recreate");
//        CrossSectionH->Write();
//    delete CrossSecF;
    
    #ifdef PrintEps
        TCanvas* crossc = new TCanvas("CrossC","CrossC");
        CrossSectionH->SetStats(0);
        CrossSectionH->Draw();
        crossc->Print(Form(("./Images/"+AnalysisString+"/CrossSection.eps").c_str()),".eps");
        delete crossc;
    #endif
}


void AntineutrinoSpectrum :: AntineutrinoSpectrumMain(bool mode)
{
    Mode = mode;
    
    std::cout <<  "\t ***********************************************************************************************" << std::endl;
    std::cout << "\t Calculating Antineutrino Spectrum" << std::endl;
    
    if((IsotopeMatrix||ReactorPowerMatrix)&&Mode==1)
    {
        sprintf(filenameReactorFile,"./RootOutputs/Reactor/RandomOutputs/ReactorSpectrum_Isotope_%d_Power_%d.root",IsotopeMatrix,ReactorPowerMatrix);
    }
    else
    {
        sprintf(filenameReactorFile,"./RootOutputs/Reactor/NominalOutputs/ReactorSpectrum.root");
    }
    
    TFile* TotalReactorSpectrumF = new TFile(filenameReactorFile);
    
    for(Int_t reactor=0; reactor<NReactors;reactor++)
    {
        TotalReactorSpectrumH[reactor] = (TH1D*)gDirectory->Get(Form("SpectrumFromReactor%i", reactor+1));
        TotalReactorSpectrumH[reactor]->SetStats(0);
        TotalReactorSpectrumH[reactor]->GetXaxis()->SetTitle("E [MeV]");
        TotalReactorSpectrumH[reactor]->GetYaxis()->SetTitle("#bar{#nu_{e}} [MeV^{-1} fissions^{-1}]");
        TotalReactorSpectrumH[reactor]->SetTitle("#nu spectrum");
    }
    
    delete TotalReactorSpectrumF;
    
    LoadCrossSectionFile();
    
    if((IsotopeMatrix||ReactorPowerMatrix)&&Mode==1)
    {
        sprintf(filenameAntineutrino,("./RootOutputs/"+AnalysisString+"/RandomOutputs/AntineutrinoSpectrum_Isotope_%d_Power_%d.root").c_str(),IsotopeMatrix,ReactorPowerMatrix);
    }
    else
    {
        sprintf(filenameAntineutrino,("./RootOutputs/"+AnalysisString+"/NominalOutputs/AntineutrinoSpectrum.root").c_str());
    }
    
    TFile* AntineutrinoSpectrumF = new TFile(filenameAntineutrino,"recreate");
    
    for(Int_t reactor=0; reactor<NReactors;reactor++)
    {
        
        AntineutrinoSpectrumH[reactor] = (TH1D*)TotalReactorSpectrumH[reactor]->Clone(Form("AntineutrinoSpectrumFromReactor%d",reactor+1));
        AntineutrinoSpectrumH[reactor]->Reset();
        AntineutrinoSpectrumH[reactor]->Multiply(CrossSectionH,TotalReactorSpectrumH[reactor]);

        if(FlatSpectrumCheck)
        {
            AntineutrinoSpectrumH[reactor]->Reset();
            for(Int_t i=0;i<Nbins;i++)
            {
                AntineutrinoSpectrumH[reactor]->SetBinContent(i+1,1);
            }
        }
        
        AntineutrinoSpectrumH[reactor]->Write();
        
        // DrawAntineutrinoSpectrum(TotalReactorSpectrumH[reactor], CrossSection, AntineutrinoSpectrumH[reactor]); //Because the pad generates warnings. Uncomment if you want to see the plot of the 3 histograms.
    }
    delete AntineutrinoSpectrumF;

    for(Int_t reactor=0; reactor<NReactors;reactor++)
    {
        delete AntineutrinoSpectrumH[reactor];
        delete TotalReactorSpectrumH[reactor];
    }
    
    delete CrossSectionH;

    
    std::cout << "\t Antineutrino Spectrum calculated" << std::endl;
    
}

//void AntineutrinoSpectrum :: DrawAntineutrinoSpectrum(TH1D* TotalReactorSpectrum, TH1D* CrossSection, TH1D* AntineutrinoSpectrum)
//{
//    TCanvas *c2 = new TCanvas("c2","Antineutrino spectrum",600,400);
//
//    TotalReactorSpectrum->Draw();
//
//    c2->Update();
//
//    Float_t rightmax = 1.1*CrossSection->GetMaximum();
//    Float_t scale = gPad->GetUymax()/rightmax;//?????????????????????????????????????????????????????????????????????????????????????????
//    CrossSection->SetLineColor(kRed);
//    CrossSection->Scale(scale);
//    CrossSection->Draw("same");
//
//    AntineutrinoSpectrum->Scale(5);//?????????????????????????????????????????????????????????????????????????????????????????
//    AntineutrinoSpectrum->SetLineColor(8);
//    AntineutrinoSpectrum->GetXaxis()->SetTitle("E_{#nu} [MeV]");
//    AntineutrinoSpectrum->GetYaxis()->SetTitle("#nu's");
//    
//    AntineutrinoSpectrum->Draw("same");
//    
//    //draw an axis on the right side
//    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
//    axis->SetTitle("IBD cross section [10^{-42} cm^{-2}]");
//    axis->SetLineColor(kRed);
//    axis->SetTextColor(kRed);
//    axis->Draw();
//    
//    c2->Write();
//    
//    c2->Close(); //To avoid showing the canvas.
//    
//}
