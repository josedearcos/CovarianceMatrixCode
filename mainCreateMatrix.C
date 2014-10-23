#include "CreateEnergyMatrix.h"
#include "NominalData.h"

int main()
{
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);
    
    const bool Analysis =0;//0 for Gd, 1 for H
    Int_t DataSet=2;
    
    NominalData* Data = new NominalData(Analysis,DataSet);
    Double_t sin22t13 = Data->GetSin22t13();
    Double_t deltaM = Data->GetDm231();
    Double_t week = Data->GetWeeks();
    
    CreateEnergyMatrix* Matrix = new CreateEnergyMatrix(Data);

    Matrix->GenerateEnergyMatrix(sin22t13,deltaM,week);
    
    return 0;
}