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

int plot_oscillation_curves(void)
{
    TH1::AddDirectory(kFALSE);
    
TH1F* OscillationH[6][3];

TFile* PlotF7 = TFile::Open("../RootOutputs/OscillationCurves.root");

TCanvas *c7 = new TCanvas("Oscillation Curves","Oscillation Curves",1200,600);
c7->Divide(6,3);
Int_t cont = 0;
for(Int_t near=0;near<3;near++)
{
    for(Int_t reactor=0;reactor<6;reactor++)
    {
        ++cont;
        c7->cd(cont);
        
        OscillationH[reactor][near]=(TH1D*)gDirectory->Get(Form("OscProbFromReactor%itoNear%i",reactor,near));
        OscillationH[reactor][near]->SetTitle(Form("OscillationReactor%iNearAD%i",reactor+1,near+1));
        OscillationH[reactor][near]->GetXaxis()->SetTitle("E_{true} (MeV)");
        OscillationH[reactor][near]->GetXaxis()->SetTitleSize(0.05);
        OscillationH[reactor][near]->GetYaxis()->SetTitleSize(0.05);
        OscillationH[reactor][near]->GetYaxis()->SetTitle("Survival Probability");
        
        OscillationH[reactor][near]->Draw();
        
    }
}
c7->Print("../Images/OscillationCurves.eps", "png");

PlotF7->Close();
}