# include <stdio.h>
# include "NominalData.h"

class ChiSquare
{
    
public:
    Double_t chiSq ;
    Double_t binWidth;
    Double_t Detected;
    Double_t Expected;
    Double_t penalty;
    Double_t s22t13;
    Double_t Nmls22t13;
    Double_t Nmls22t13Error;
    Int_t Nweeks;
    Int_t Samples;
    NominalData Nom;
    
    ChiSquare();
    //ChiSquare(NominalData Data, OscillationReactor Osc);
    ChiSquare(NominalData Data, Oscillation Osc);
    Double_t ChiSquareCalc();

};

ChiSquare :: ChiSquare ()
{
    Nom = new NominalData();
    chiSq = 0;
    Nmls22t13 = Nom.GetSin22t13();
    Nmls22t13Error = Nom.GetSin22t13Error();
}

//ChiSquare :: ChiSquare (NominalData Data, OscillationReactor Osc)
ChiSquare :: ChiSquare (NominalData Data, Oscillation Osc)
{
    chiSq = 0;
    Nmls22t13 = Data.GetSin22t13();
    Nmls22t13Error = Data.GetSin22t13Error();
    Nweeks = Data.GetWeeks();
    s22t13 = Osc.GetSin22t13();//Random s22t13

}

Double_t ChiSquare :: ChiSquareCalc()
{

    //Input Spectrums
    TFile *OscillationF = TFile::Open("./RootOutputs/Oscillation.root");

    OscillationF->cd("EH spectra after oscillation");

    TH1F *TotalOscillatedSpectrumEH3 = (TH1F*)gDirectory->Get("Spectrum after oscillation at EH3");
    TotalOscillatedSpectrumEH3->SetDirectory(0);
    OscillationF->Close();
    
    TFile *RandomOscillationF = TFile::Open("./RootOutputs/RandomOscillation.root");
    //  TFile *RandomOscillationF= TFile::Open("./RootOutputs/Oscillation.root"); //File to check if ChiSq is 0 when both inputs are the same. It works!
    RandomOscillationF->cd("EH spectra after oscillation");

    TH1F *TotalRandomSpectrumEH3 = (TH1F*)gDirectory->Get("Spectrum after oscillation at EH3");
    TotalRandomSpectrumEH3->SetDirectory(0);
    RandomOscillationF->Close();
    
    binWidth = TotalRandomSpectrumEH3->GetBinWidth(1);
    for(int idx = 0; idx<TotalRandomSpectrumEH3->GetNbinsX(); idx++){
        Detected = TotalOscillatedSpectrumEH3->GetBinContent(idx+1);
        Expected = TotalRandomSpectrumEH3->GetBinContent(idx+1);
        chiSq += 2*(Expected-Detected);
        
        if(Detected!=0){
            chiSq += 2*Detected*TMath::Log(Detected/Expected);
        }

    }
    
    //cout << "Sin22t13 " << Nmls22t13 << "\n";
    //cout << "Sin22t13Error is " << Nmls22t13Error << "\n";
    //cout << " chiSq without penalty " <<  chiSq << "\n";
    
    // Add penalty based on current best estimate of sin22t13 value
    penalty = ( s22t13 - Nmls22t13) / Nmls22t13Error ;
    
    //cout << "  penalty " <<  penalty  << "\n";

    chiSq += penalty*penalty;

    //cout << "ChiSquare value is " << chiSq  << "\n";
    
    return chiSq;
}