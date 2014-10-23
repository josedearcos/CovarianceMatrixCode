#include "TFile.h"
#include "TH1.h"
#include <float.h>

const Int_t NReactors = 6;
const Int_t Nweeks = 1;
//I need to make all the bins that are 0 a little bit bigger than 0, otherwise I get NaNs in the fraction and extrapolation factors.

void ReadInclusive()
{

    Char_t filenameReactor[1024];
    Char_t filenameFlux[1024];

    TFile* ReactorDataF = TFile::Open("./WeeklyFlux_20week_unblinded_inclusive.root");
   
    TH1F* ReactorSpectrumH[NReactors][Nweeks];

    for(Int_t reactor = 0; reactor < NReactors; reactor++)
    {
        for(Int_t week = 0; week < Nweeks; week++)
        {
            
            sprintf(filenameReactor,"Week%i/%i", week, reactor);
            sprintf(filenameFlux,"Week%i/%i", week, reactor);
            
            ReactorSpectrumH[reactor][week] = (TH1F*)gDirectory->Get(filenameReactor);
            ReactorSpectrumH[reactor][week]->SetDirectory(0);//Somehow this solves all crashes

        }
    }
    ReactorDataF->Close();
    
    TFile* FluxF = TFile::Open("./FormatedInclusiveFlux.root","recreate");
    
    TH1F* Flux[NReactors];
    for(Int_t reactor = 0; reactor < NReactors; reactor++)
    {
        sprintf(filenameFlux,"ZeroCorrectedInclusiveFlux%i",reactor+1);
        Flux[reactor]= new TH1F(filenameFlux,filenameFlux,41,1.8,10);
        for(Int_t pts =1;pts<=41;pts++){
            if(ReactorSpectrumH[reactor][0]->GetBinContent(pts)==0)
                {
                    
                    Flux[reactor]->SetBinContent(pts,ReactorSpectrumH[reactor][0]->GetBinContent(35));
//                    cout<<Flux[reactor]->GetBinContent(pts)<<"\n";
                }
            else
            {
                Flux[reactor]->SetBinContent(pts,ReactorSpectrumH[reactor][0]->GetBinContent(pts));
            }
        }
        Flux[reactor]->Write();
    }
    FluxF->Close();

}