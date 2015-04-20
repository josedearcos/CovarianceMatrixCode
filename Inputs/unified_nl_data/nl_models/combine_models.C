{
  const Int_t nmodels = 5;
  TString input_filename[nmodels]
    = {"Model1.root",
       "flt_glb_marg3.root",
       "eff_model1.root",
       "eff_model6.root",
       "lbnl_model.root"};

  TString gamma_name[nmodels] 
    = {"gammaFull",
       "gammaFull",
       "gamma",
       "gamma",
       "gamma_lbnl"};
  
  TString electron_name[nmodels]
    = {"electronFull",
       "electronFull",
       "electron",
       "electron",
       "electron_lbnl"};
  TString positron_name[nmodels]
    = {"positronFull",
       "positronFull",
       "positron",
       "positron",
       "positron_lbnl"};

  TFile * f_in[nmodels];
  for (Int_t i = 0; i < nmodels; i++){
    f_in[i] = new TFile(input_filename[i].Data());
  }
  
  TFile * fout = new TFile("nl_models_final.root","recreate");
  TGraph * g_gamma[nmodels];
  TGraph * g_electron[nmodels];
  TGraph * g_positron[nmodels];

  for (Int_t i = 0; i < nmodels; i++){
    g_gamma[i] = (TGraph*)f_in[i]->Get(gamma_name[i].Data());
    g_gamma[i]->SetName(Form("gamma_%d",i));
    g_gamma[i]->Write();

    g_electron[i] = (TGraph*)f_in[i]->Get(electron_name[i].Data());
    g_electron[i]->SetName(Form("electron_%d",i));
    g_electron[i]->Write();

    g_positron[i] = (TGraph*)f_in[i]->Get(positron_name[i].Data());
    g_positron[i]->SetName(Form("positron_%d",i));
    g_positron[i]->Write();
  }
  fout->Close();
}
