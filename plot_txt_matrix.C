void plot_txt_matrix()
{
    //To draw using a better palette:
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;
    
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    
    const Int_t n_evis_bins = 37;
    const Int_t n_etrue_bins = 39;
    Double_t evis_bins[37+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
    Double_t enu_bins[39+1]; // 39 bins between 1.8 and 9.6 MeV
    Double_t Energy[n_etrue_bins][n_evis_bins];
    evis_bins[0] = 0.7;
    for (Int_t i = 0; i < n_evis_bins-1; i++)
    {
        evis_bins[i+1] = 0.2*i + 1.0;
    }
    evis_bins[n_evis_bins] = 12;
    
    for (Int_t i = 0; i < n_etrue_bins+1; i++)
    {
        enu_bins[i] = 0.2 * i + 1.8;
    }
    TFile* File = new TFile("./matrix_evis_to_enu_unified_p12e_unblinded.root","recreate");
    TH2F* Hist = new TH2F("EvisEnu","EvisEnu",n_etrue_bins,enu_bins,n_evis_bins,evis_bins);
    TH2F* Hist2 = new TH2F("EnuEvis","EnuEvis",n_etrue_bins,enu_bins,n_evis_bins,evis_bins);

    std::ifstream input("matrix_evis_to_enu_unified_p12e_unblinded.txt");
    std::string line;

    for(Int_t i=0;i<n_evis_bins;i++)
    {
        for(Int_t j=0;j<n_etrue_bins;j++)
        {
            input >> Energy[j][i];
        }
    }
    for(Int_t i=0;i<n_etrue_bins;i++)
    {
        for(Int_t j=0;j<n_evis_bins;j++)
        {
            Hist->SetBinContent(i+1,j+1,Energy[i][j]);
        }
    }
    
    Double_t Normx[n_etrue_bins]={};
    Double_t Normy[n_evis_bins]={};

    for(Int_t i=0;i<n_etrue_bins;i++)
    {
        for(Int_t j=0;j<n_evis_bins;j++)
        {
            //Check normalization
            Normx[i]=Normx[i]+Hist->GetBinContent(i+1,j+1);
            Normy[j]=Normy[j]+Hist->GetBinContent(i+1,j+1);
            
        }
    }
    
    for(Int_t i=0;i<n_etrue_bins;i++)
    {
        for(Int_t j=0;j<n_evis_bins;j++)
        {
            Hist2->SetBinContent(i+1,j+1,Energy[i][j]/Normx[i]);
        }
    }
    
    
    for(Int_t i=0;i<n_etrue_bins;i++)
    {
        for(Int_t j=0;j<n_evis_bins;j++)
        {
            cout << "Norma in i" << Normx[i]<<endl;
            cout << "Norma in j" << Normy[j]<<endl;//This shows it's normalized using the sum for a fixed y bin.
            
        }
    }
//
   Hist->Write();
    Hist2->Write();

   File->Close();

}