{
  const Int_t npoints_electron = 300;
  const Int_t npoints_positron = 300;
  const Int_t npoints_gamma = 19;
  Double_t e_electron[npoints_electron];
  Double_t nl_electron[npoints_electron];
  Double_t e_positron[npoints_positron];
  Double_t nl_positron[npoints_positron];
  Double_t e_gamma[npoints_gamma];
  Double_t nl_gamma[npoints_gamma];

  ifstream electron_file("lbnl_electron_nl.txt");
  for (Int_t i = 0; i < npoints_electron; i++){
    electron_file >> e_electron[i] >> nl_electron[i];
  }
  ifstream gamma_file("lbnl_gamma_nl.txt");
  for (Int_t i = 0; i < npoints_gamma; i++){
    gamma_file >> e_gamma[i] >> nl_gamma[i];
  }
  ifstream positron_file("lbnl_positron_nl.txt");
  Double_t dummy;
  for (Int_t i = 0; i < npoints_positron; i++){
    positron_file >> e_positron[i] >> nl_positron[i] >> dummy >> dummy >> dummy >> dummy;
  }

  TGraph * g_electron = new TGraph(npoints_electron,e_electron, nl_electron);
  TGraph * g_gamma = new TGraph(npoints_gamma,e_gamma, nl_gamma);
  TGraph * g_positron = new TGraph(npoints_positron,e_positron, nl_positron);

  g_electron->SetName("electron_lbnl");
  g_gamma->SetName("gamma_lbnl");
  g_positron->SetName("positron_lbnl");

  g_electron->Draw("AL");
  g_gamma->SetLineColor(4);
  g_gamma->Draw("L");
  g_positron->SetLineColor(2);
  g_positron->Draw("L");

  TFile * fout = new TFile("lbnl_model.root","recreate");
  
  g_electron->Write();
  g_gamma->Write();
  g_positron->Write();

  fout->Close();
  


}
