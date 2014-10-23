int StatisticalFluctuationCheck()
{
    TFile* DataF = TFile::Open("./RootOutputs/Spectra/PredictedSpectrum.root");
    TH1F* Histo = (TH1F*)gDirectory->Get(Get(Form("Combined Near Prediction AD%d",near));
   DataF->Close();

    
    for(Int_t VisibleEnergyIndex=1;VisibleEnergyIndex<=Histo->GetXaxis()->GetNbins();VisibleEnergyIndex++)
    {
         Histo->SetBinContent(VisibleEnergyIndex,(Double_t)(rand->PoissonD(Histo->GetBinContent(VisibleEnergyIndex)*Histo->GetXaxis()->GetBinWidth(VisibleEnergyIndex))/Histo->GetXaxis()->GetBinWidth(VisibleEnergyIndex)));
    }
}
);