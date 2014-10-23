void EnergyComparison()
{
    
    TFile* PREF = TFile::Open("./Energy51.root");
    
    TH2F* PreH = (TH2F*)gDirectory->Get("EvisEnu2");
    PreH->SetDirectory(0);
    PREF->Close();

    TFile* POSTF = TFile::Open("./Energy240.root");
    
    TH2F* PostH = (TH2F*)gDirectory->Get("EvisEnu2");
    PostH->SetDirectory(0);
    POSTF->Close();
    TH2F* DifferenceH = (TH2F*)PostH->Clone("RelativeError");
    DifferenceH->Add(PreH,-1);
    
    TFile* ResultF = TFile::Open("./EnergyMatrixDifference.root","recreate");
    DifferenceH->SetTitle("Relative Error (Truncated)");
    
    for (Int_t i=0; i<60; i++)
    {
        for (Int_t j=0; j<60; j++)
        {
            if(PreH->GetBinContent(i+1,j+1)>0.1)//Truncates the error, otherwise for really small bin values we got 10^6 relative errors.
            {
                DifferenceH->SetBinContent(i+1,j+1,(DifferenceH->GetBinContent(i+1,j+1)/PreH->GetBinContent(i+1,j+1)));
            }
            else
            {
                DifferenceH->SetBinContent(i+1,j+1,0);
            }
        }
    }
    DifferenceH->Write();
    PreH->Write();
    PostH->Write();
    DifferenceH->Draw("samecolz");
    PostH->Draw("samebox1");

    PreH->Draw("samebox");

    ResultF->Close();

}