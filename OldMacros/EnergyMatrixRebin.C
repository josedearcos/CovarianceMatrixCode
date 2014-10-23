#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TH2.h"

#include "TMath.h"

#include "NominalData.h"
#include <math.h>

const Int_t MaxDetectors = 8;

void EnergyMatrixRebin()
{
    Int_t MatrixBins = 240;
    Int_t NADs=6;
    Int_t Nbins = 51;
    Double_t InitialEnergy = 1.8;
    Double_t FinalEnergy = 12;
    
TH2F* EnergyM[MaxDetectors];
TH2F* EnergyMatrixH[MaxDetectors];
TH1F* ProjectionH[MaxDetectors];
TH2* EnergyRebinH[MaxDetectors];

Int_t RebinFactor = Nbins + InitialEnergy/((FinalEnergy-InitialEnergy)/Nbins);

TFile* EnergyMatrixF = TFile::Open("./RootOutputs/EnergyMatrix1.root");

for (Int_t ad =0; ad<1; ad++)//CHANGE LIMTIS
{
    cout << "LOADING" << endl;
    EnergyMatrixH[ad] = (TH2F*)gDirectory->Get(Form("EvisEnu%i",ad));
//    EnergyMatrixH[ad]->SetDirectory(0);
}

EnergyMatrixF->Close();

    for(int idet=0; idet<1; idet++)
    {

        EnergyRebinH[idet] = EnergyMatrixH[idet]->Rebin2D(MatrixBins/RebinFactor,MatrixBins/RebinFactor);
        
        EnergyRebinH[idet]->Draw("cold");
        TFile* EnergyRebinF = TFile::Open("./RootOutputs/EnergyMatrixRebin.root","recreate");

        EnergyRebinH[idet]->Write();
        cout << "WRITE" << endl;
        //        delete EnergyM[idet];
        EnergyRebinF->Close();

    }

    
}