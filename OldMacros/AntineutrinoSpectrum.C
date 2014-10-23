//
//  AntineutrinoSpectrum.C
//  I multiply the reactour spectrum by the neutrino cross section to get the True Energy Spectrum without oscillation.
//
//  Created by Jose de Arcos	 on 2/18/13.
//
//

TH1F* AntineutrinoSpectrum(Int_t Nbins, Double_t InitialEnergy, Double_t FinalEnergy){
    
    TFile *CrossSectionF = TFile::Open("./RootOutputs/CrossSection.root");
    
    TH1F *CrossSection = (TH1F*)gDirectory->Get("CrossSection");
    TFile *TotalReactorSpectrumF = TFile::Open("./RootOutputs/IsotopeSpectra.root");
    TH1F *TotalReactorSpectrum = (TH1F*)gDirectory->Get("TotalReactorSpectrum");
    
    TFile *AntineutrinoSpectrumF = new TFile("./RootOutputs/AntineutrinoSpectrum.root","recreate");
    TH1F *AntineutrinoSpectrum = new TH1F("AntineutrinoSpectrumH", "Antineutrino Spectrum", Nbins, InitialEnergy, FinalEnergy);
    
    AntineutrinoSpectrum->Multiply(CrossSection,TotalReactorSpectrum);
    
    TCanvas *c2 = new TCanvas("c2","Antineutrino spectrum",600,400);
    TotalReactorSpectrum->SetStats(0);
    TotalReactorSpectrum->GetXaxis()->SetTitle("E [MeV]");
    TotalReactorSpectrum->GetYaxis()->SetTitle("#bar{#nu_{e}} [MeV^{-1} fissions^{-1}]");
    
    TotalReactorSpectrum->SetTitle("#nu spectrum");
    TotalReactorSpectrum->Draw();
    c2->Update();
    
    Float_t rightmax = 1.1*CrossSection->GetMaximum();
    Float_t scale = gPad->GetUymax()/rightmax;
    CrossSection->SetLineColor(kRed);
    CrossSection->Scale(scale);
    CrossSection->Draw("same");
    
    AntineutrinoSpectrum->Scale(5);
    AntineutrinoSpectrum->SetLineColor(8);
    AntineutrinoSpectrum->GetXaxis()->SetTitle("E_{#nu} [MeV]");
    AntineutrinoSpectrum->GetYaxis()->SetTitle("#nu's");
    
    AntineutrinoSpectrum->Draw("same");
    
    //draw an axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetTitle("IBD cross section [10^{-42} cm^{-2}]");
    axis->SetLineColor(kRed);
    axis->SetTextColor(kRed);
    axis->Draw();
    
    AntineutrinoSpectrumH->Write();
    c2->Write();
    
    return AntineutrinoSpectrumH;
    
}