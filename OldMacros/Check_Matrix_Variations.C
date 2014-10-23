void Check_Matrix_Variations()
{
    TH1::AddDirectory(kFALSE);

    Char_t filenameResponse[100];

    TH2F* NominalResponseMatrixH[6];
    TH2F* NLResponseMatrixH[6];
    TH2F* IAVResponseMatrixH[6];
    TH2F* ResoResponseMatrixH[6];
    TH2F* DifferenceNL[6];
    TH2F* DifferenceIAV[6];
    TH2F* DifferenceReso[6];

    sprintf(filenameResponse,"./NominalResponseMatrices/LBNLBinning/NominalResponseBCWModel.root");
    
    TFile* ResponseF = TFile::Open(filenameResponse);
    
    for (Int_t AD =0; AD<6; AD++)
    {
        NominalResponseMatrixH[AD] = (TH2F*)gDirectory->Get(Form("EvisEnu%i;2",AD));
    }
    ResponseF->Close();
    
    TFile* NLMatrix1F = TFile::Open("./CrossChecks/EnergyMatrixComparison/NLalteredResponseMatrix.root");
    for (Int_t AD =0; AD<6; AD++)
    {
        NLResponseMatrixH[AD] = (TH2F*)gDirectory->Get(Form("EvisEnu%i;2",AD));
    }
    NLMatrix1F->Close();
        
    TFile* IAVMatrix1F = TFile::Open("./CrossChecks/EnergyMatrixComparison/IAValteredResponseMatrix.root");
    for (Int_t AD =0; AD<6; AD++)
    {
        IAVResponseMatrixH[AD] = (TH2F*)gDirectory->Get(Form("EvisEnu%i;2",AD));
    }
    IAVMatrix1F->Close();
    
    TFile* ResoMatrix1F = TFile::Open("./CrossChecks/EnergyMatrixComparison/ResoalteredResponseMatrix.root");
    for (Int_t AD =0; AD<6; AD++)
    {
        ResoResponseMatrixH[AD] = (TH2F*)gDirectory->Get(Form("EvisEnu%i;2",AD));
    }
    ResoMatrix1F->Close();
    
    TFile* ResultsF = TFile::Open("./CrossChecks/EnergyMatrixComparison/DifferenceNominalAltered.root","recreate");
    for (Int_t i =0; i<6; i++)
    {
        DifferenceNL[i]=(TH2F*)NominalResponseMatrixH[i]->Clone(Form("NLDifferenceAD%d",i));
        DifferenceNL[i]->Add(NLResponseMatrixH[i],-1);
        DifferenceNL[i]->Write();
        DifferenceIAV[i]=(TH2F*)NominalResponseMatrixH[i]->Clone(Form("IAVDifferenceAD%d",i));
        DifferenceIAV[i]->Add(IAVResponseMatrixH[i],-1);
        DifferenceIAV[i]->Write();
        DifferenceReso[i]=(TH2F*)NominalResponseMatrixH[i]->Clone(Form("ResoDifferenceAD%d",i));
        DifferenceReso[i]->Add(ResoResponseMatrixH[i],-1);
        DifferenceReso[i]->Write();
    }
    
    ResultsF->Close();
}