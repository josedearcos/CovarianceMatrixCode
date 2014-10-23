#include "TFile.h"
#include "TH1.h"

//Fast check finding the peaks of the histograms and calculating Erec/Etrue, we should get a trend similar to the NL function.
void NLCrossCheck()
{
    
    TFile* F = TFile::Open("../ResponseMatrices/NominalResponseMatrix.root");
    TH1F* NLH[240];
    TH1F* IAVH[240];
    Double_t MaxIAV[240];
    Double_t MaxNL[240];
    
    for (Int_t i=0; i<240; i++)
    {
        NLH[i] = (TH1F*)gDirectory->Get(Form("NLSpectrum0,%d",i));
        NLH[i]->SetDirectory(0);
        IAVH[i] = (TH1F*)gDirectory->Get(Form("IAVSpectrum0,%d",i));
        IAVH[i]->SetDirectory(0);
        MaxIAV[i]=IAVH[i]->GetMaximumBin();
        MaxNL[i]=NLH[i]->GetMaximumBin();
    }
    
    F->Close();

    TFile* ResultF = TFile::Open("./NLCheck.root","recreate");
    TH1F* Check = new TH1F("NLCHeck","NLCheck", 240, 0, 12);
    for (Int_t i=0; i<240; i++)
    {
        NLCHeck->SetBinContent(i+1,(MaxNL[i])/MaxIAV[i]);//Erec/Etrue
        cout <<MaxNL[i]/MaxIAV[i] << endl;
    }
    NLCHeck->Write();
    NLCHeck->Draw("");

    
    ResultF->Close();
    
}