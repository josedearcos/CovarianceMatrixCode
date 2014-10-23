//
//  ReactorSpectrum.C
//  
//
//  Created by Jose de Arcos	 on 2/8/13.
//
//

#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>

void ReactorSpectrum(Int_t Nbins, Double_t InitialEnergy, Double_t FinalEnergy){
    
    //Nominal fractions
    Double_t NmlU235Frac  = 0.64;
    Double_t NmlU238Frac  = 0.08;
    Double_t NmlPu239Frac = 0.25;
    Double_t NmlPu241Frac = 0.03;
   
    
    //For Pu239
  
    std::ifstream input("./SpectrumData/antiNeuFlux_Pu2392011.dat.txt");
    std::string line;

    getline( input, line ); //To throw away the first two lines
    getline( input, line );
    
   // Double_t Energy[Nbins];
   // Double_t Pu239Spectrum[Nbins];
    Double_t *Energy;
    Energy = (Double_t*)malloc(Nbins*sizeof(Double_t));
    Double_t *Pu239Spectrum;
    Pu239Spectrum = (Double_t*)malloc(Nbins*sizeof(Double_t));
    
    TH1F *Pu239 = new TH1F("Pu239", "Pu239 spectrum", Nbins, InitialEnergy, FinalEnergy);
 
    int num = 0;
    for(int i=0;i<=Nbins-1;i++){
        input >> Energy[num] >> Pu239Spectrum[num]; ++num;
        Pu239->SetBinContent(i+1, Pu239Spectrum[i]);

       //  printf("Index is %d \n",num);
       // printf("Energy is %.9f \n",Energy[num-1]);
       // printf("Pu239 Spectrum is %.10f \n",Pu239Spectrum[num-1]);
    }
  
    //For Pu241
    TH1F *Pu241= new TH1F("Pu241", "Pu241 spectrum", Nbins, InitialEnergy, FinalEnergy);

    std::ifstream input("./SpectrumData/antiNeuFlux_Pu2412011.dat.txt");
    std::string line;
    
    getline( input, line ); //To throw away the first two lines
    getline( input, line );
    
    //Double_t Pu241Spectrum[Nbins];
    Double_t *Pu241Spectrum;
    Pu241Spectrum = (Double_t*)malloc(Nbins*sizeof(Double_t));
    
    int num = 0;
    for(int i=0;i<=Nbins-1;i++){
        input >> Energy[num] >> Pu241Spectrum[num]; ++num;
        Pu241->SetBinContent(i+1, Pu241Spectrum[i]);

    }

    //For U235
    TH1F *U235 = new TH1F("U235", "U235 spectrum", Nbins, InitialEnergy, FinalEnergy);

    std::ifstream input("./SpectrumData/antiNeuFlux_U2352011.dat.txt");
    std::string line;
    
    getline( input, line ); //To throw away the first two lines
    getline( input, line );
    
    //Double_t U235Spectrum[Nbins];
    Double_t *U235Spectrum;
    U235Spectrum = (Double_t*)malloc(Nbins*sizeof(Double_t));
    
    int num = 0;
    for(int i=0;i<=Nbins-1;i++){
        input >> Energy[num] >> U235Spectrum[num]; ++num;
        U235->SetBinContent(i+1, U235Spectrum[i]);

    }
    
    //For U238
    TH1F *U238 = new TH1F("U238", "U238 spectrum", Nbins, InitialEnergy, FinalEnergy);

    std::ifstream input("./SpectrumData/antiNeuFlux_U2382011.dat.txt");
    std::string line;
    
    getline( input, line ); //To throw away the first two lines
    getline( input, line );
    
    //Double_t U238Spectrum[Nbins];
    Double_t *U238Spectrum;
    U238Spectrum = (Double_t*)malloc(Nbins*sizeof(Double_t));
    
    int num = 0;
    for(int i=0;i<=Nbins-1;i++){
        input >> Energy[num] >> U238Spectrum[num]; ++num;
        U238->SetBinContent(i+1, U238Spectrum[i]);

    }
    
    //Produce Histograms:


    Pu239->GetXaxis()->SetTitle("E_{#nu} [MeV]");
    Pu239->GetYaxis()->SetTitle("Pu239");
    //Pu239->Draw();
    Pu241->GetXaxis()->SetTitle("E_{#nu} [MeV]");
    Pu241->GetYaxis()->SetTitle("Pu241");
    //Pu241->Draw();

    U235->GetXaxis()->SetTitle("E_{#nu} [MeV]");
    U235->GetYaxis()->SetTitle("U235");
    //U235->Draw();

    U238->GetXaxis()->SetTitle("E_{#nu} [MeV]");
    U238->GetYaxis()->SetTitle("U238");
    //U238->Draw();

    //Save histograms:
  
    TFile f("./RootOutputs/IsotopeSpectra.root","recreate");
    Pu239->Write();
    Pu241->Write();
    U235->Write();
    U238->Write();
       
    U235->Scale(NmlU235Frac);
    U238->Scale(NmlU238Frac);
    Pu239->Scale(NmlPu239Frac);
    Pu241->Scale(NmlPu241Frac);
    
    TH1F *TotalReactorSpectrum = new TH1F("TotalReactorSpectrum", "Total reactor spectrum", Nbins, InitialEnergy , FinalEnergy);
    
    TotalReactorSpectrum->GetXaxis()->SetTitle("E_{#nu} [MeV]");
    TotalReactorSpectrum->GetYaxis()->SetTitle("#nu [MeV^{-1} fission^{-1}]");

    
    TotalReactorSpectrum->Add(U235);
    TotalReactorSpectrum->Add(U238);
    TotalReactorSpectrum->Add(Pu239);
    TotalReactorSpectrum->Add(Pu241);
    
    TotalReactorSpectrum->Write();
    
}