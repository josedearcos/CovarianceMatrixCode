#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "NominalData.h"

using namespace std;

class ReactorSpectrum
{

private:
    const char* title;
    Int_t Nbins;
    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t* IsotopeFrac;
    Double_t IsotopeFracError;
    Double_t IsotopeError;
    Double_t* Energy;
    Double_t* Spectrum;
    TRandom3 *rand;
    Double_t Norm;
    NominalData *Nom;
    char* FileName;
    Int_t ADs;
    Int_t Reactors;
    Int_t OriginalNbins;
    
    TH1F* File2Hist(const char* IsotopeName);
    void Plot(TH1F* SpectrumHist);
    void TotalReactorSpectrum(TH1F* SpectrumPu239,TH1F* SpectrumPu241,TH1F* SpectrumU235,TH1F* SpectrumU238);
    void RandomIsotopes();
    
public:
    
    ReactorSpectrum();
    ReactorSpectrum(NominalData Data);
    void SimpleReactorSpectrumMain(bool Random);
  
};


//Constructor

//Default values
ReactorSpectrum :: ReactorSpectrum()
{
    
    Nom = new NominalData();
    
    //Binning variables
    OriginalNbins = 821;
    Nbins = Nom->GetNbins();
    InitialEnergy = Nom->GetEmin();
    FinalEnergy = Nom->GetEmax();
    
    IsotopeFrac = Nom->GetIsotopeFraction(); // Load nominal Isotope fractions
    
    //Isotope errors (5% for all) http://dayabay.ihep.ac.cn/DocDB/0086/008609/002/reactor_technote.pdf
    IsotopeError = 0.05;
    Norm = 0;
    
    ADs = Nom->GetADs();
    Reactors= Nom->GetReactors();
    
    FileName = "./RootOutputs/IsotopeSpectra.root";
    
}

ReactorSpectrum :: ReactorSpectrum(NominalData Data)
{
    
    //Binning variables
    OriginalNbins = 821;
    Nbins = Data.GetNbins();
    InitialEnergy = Data.GetEmin();
    FinalEnergy = Data.GetEmax();
    
    IsotopeFrac = Data.GetIsotopeFraction(); // Load nominal Isotope fractions
    
    //Isotope errors (5% for all) http://dayabay.ihep.ac.cn/DocDB/0086/008609/002/reactor_technote.pdf
    IsotopeError = 0.05;
    Norm = 0;
    
    ADs = Data.GetADs();
    Reactors= Data.GetReactors();
    
    FileName = "./RootOutputs/IsotopeSpectra.root";
    
}

//Use as inputs antiNeuFlux_Pu2392011.dat.txt
//              antiNeuFlux_Pu2412011.dat.txt
//              antiNeuFlux_U2352011.dat.txt
//              antiNeuFlux_U2382011.dat.txt

void ReactorSpectrum :: SimpleReactorSpectrumMain(bool Random)
{
        
    if (Random)
    {
        FileName="./RootOutputs/RandomIsotopeSpectra.root";
        rand = new TRandom3(0);
        RandomIsotopes();
    }

    TH1F* Pu239 =  File2Hist("./ReactorInputs/antiNeuFlux_Pu2392011.dat.txt");
    TH1F* Pu241 =  File2Hist("./ReactorInputs/antiNeuFlux_Pu2412011.dat.txt");
    TH1F* U235 =  File2Hist("./ReactorInputs/antiNeuFlux_U2352011.dat.txt");
    TH1F* U238 =  File2Hist("./ReactorInputs/antiNeuFlux_U2382011.dat.txt");
    
   // TFile outputFile(FileName,"recreate");
    TFile* outputFile = TFile::Open(FileName, "recreate");
    Plot(Pu239);
    Plot(Pu241);
    Plot(U235);
    Plot(U238);
    
    TotalReactorSpectrum(Pu239,Pu241,U235,U238);
    
    outputFile->Close();
    delete[] Energy;
    delete[] Spectrum;

}

TH1F* ReactorSpectrum :: File2Hist (const char* IsotopeName)
{
    
    std::ifstream input(IsotopeName);
    std::string line;
    
    getline( input, line ); //To throw away the first two lines
    getline( input, line );
    
    Energy = (Double_t*)malloc(OriginalNbins*sizeof(Double_t));
    Spectrum = (Double_t*)malloc(OriginalNbins*sizeof(Double_t));
    
    TH1F *SpectrumHist = new TH1F(IsotopeName, IsotopeName, Nbins, InitialEnergy, FinalEnergy);
    
    for(Int_t i=0;i<=820;i++)
    {
        input >> Energy[i] >> Spectrum[i];
        
//        printf("Index is %d \n",i);
//        printf("Energy is %.9f \n",Energy[i]);
//        printf("Spectrum is %.10f \n",Spectrum[i]);     
    }


    for(Int_t i=0;i<=Nbins-1;i++)
    {
        SpectrumHist->SetBinContent(i+1, Spectrum[i*(OriginalNbins-1)/(Nbins-1)]); //821 - 1 is the original size of the data
    }
    
    SpectrumHist->SetDirectory(0);//Solves memory leakying warning
    return SpectrumHist;
}

void ReactorSpectrum :: Plot(TH1F* SpectrumHist)
{
    title = SpectrumHist->GetName();
    SpectrumHist->GetXaxis()->SetTitle("E_{#nu} [MeV]");
    SpectrumHist->GetYaxis()->SetTitle(title);
    SpectrumHist->Write();
    
}

void ReactorSpectrum :: TotalReactorSpectrum(TH1F* SpectrumPu239,TH1F* SpectrumPu241,TH1F* SpectrumU235,TH1F* SpectrumU238)
{
    cout << " Random Isotope 1" << " Fraction is " << IsotopeFrac[0] << "\n";

    SpectrumU235->Scale( IsotopeFrac[0]);
    SpectrumU238->Scale( IsotopeFrac[1]);
    SpectrumPu239->Scale( IsotopeFrac[2]);
    SpectrumPu241->Scale( IsotopeFrac[3]);
    
    TH1F *TotalReactorSpectrum = new TH1F("TotalReactorSpectrum", "Total reactor spectrum", Nbins, InitialEnergy , FinalEnergy);
    TotalReactorSpectrum->GetXaxis()->SetTitle("E_{#nu} [MeV]");
    TotalReactorSpectrum->GetYaxis()->SetTitle("#nu [MeV^{-1} fission^{-1}]");
    TotalReactorSpectrum->SetDirectory(0);
    TotalReactorSpectrum->Add(SpectrumU235);
    TotalReactorSpectrum->Add(SpectrumU238);
    TotalReactorSpectrum->Add(SpectrumPu239);
    TotalReactorSpectrum->Add(SpectrumPu241);
    
    //Now I'm going to multiply by the Power in the core
    TotalReactorSpectrum->Write();
    
}

void ReactorSpectrum :: RandomIsotopes()
{
    for (Int_t i = 0; i<4; i++)
    {
        IsotopeFracError[i] = (IsotopeError * rand->Gaus(0,1));
        IsotopeFrac[i]  = (1 + IsotopeFracError[i]) * IsotopeFrac[i];
        
         Norm = Norm + IsotopeFrac[i];
     // cout << " Random Isotope" << i << " absolute error is " << IsotopeFracError[i] << "\n";
    }
    //Normalize total isotope fraction sum to 1 after being randomized.
    for (Int_t i = 0; i<4; i++)
    {
        IsotopeFrac[i]  = IsotopeFrac[i]/Norm;
        cout << " Random Isotope" << i << " fraction is " << IsotopeFrac[i] << "\n";
    }
    
}