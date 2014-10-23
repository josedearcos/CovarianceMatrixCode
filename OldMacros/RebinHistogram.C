#include "TH1F.h"

void RebinHistogram{
    const MaxDetectors =6;
    const MaxWeeks = 20;
    TH1F* NearHallSpectrumOriginalH[MaxDetectors][MaxWeeks];
    
    NOriginalBins = 37;
    
    Double_t eoriginal_bins[38]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
    
    evis_bins[0] = 0.7;
    for (Int_t i = 0; i < 36; i++){
        eoriginal_bins[i+1] = 0.2 *i + 1.0;
    }
    eoriginal_bins[37] = 12.0;
    
    Double_t e_bins[41+1]; // 41 bins between 1.8 and 10 MeV
    for (Int_t i = 0; i < 41+1; i++){
        e_bins[i] = 0.2 * i + 1.8;
    }
    
    ADdata = "./Inputs/ibd_eprompt_shapes.root";//IBD Gd data
    for(Int_t week=0;week<Nweeks;++week)
    {
        for(Int_t ad=0;ad<NADs;++ad)
        {
            
    if (Nweeks == 1)
    {
        sprintf(name,"h_ibd_eprompt_inclusive_ad%d",ad+1);
    }
    else
    {
        sprintf(name,"h_ibd_eprompt_week%i_ad%i",week,ad+1);
    }
            NearHallSpectrumOriginalH[ad][week] = (TH1F*)gDirectory->Get(name);

            TH1F *NearHallSpectrumH = new TH1F(name,name, Nbins, InitialEnergy, FinalEnergy);
      
            for(Int_t i = 1; i<=NOriginalBins; i++)
            {
                value = NearHallSpectrumOriginalH->GetBinContent(i);
                
                NearHallSpectrumH
                
            
            }
        }
    }



}