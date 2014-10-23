#include "TH2F.h"
#include "TFile.h"
#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2.h"
#include "TMath.h"
#include <math.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TCanvas.h"

int InversionCheck(void)
{
    TH1::AddDirectory(kFALSE);

    TH2F* TotalCovarianceMatrixH;
    TArrayD TotalCovarianceMatrixM(37*37*1);
    Double_t InvTotalCovarianceMatrixM[37][37][1];

    TFile* SaveCovarianceMatricesF = new TFile("./CovarianceMatrices/FitterCovarianceMatrixResultsPeriod0.root");
    TotalCovarianceMatrixH = (TH2F*)gDirectory->Get("Total Cov Matrix");
    SaveCovarianceMatricesF->Close();
    
//    TotalCovarianceMatrixH->Draw("colz");
//    Double_t* CombineMatrices();

    TH2F* InvTotalCovarianceMatrixH[1];
    TH2F* TotalCovarianceMatrixH[1];

    InvTotalCovarianceMatrixH[0] = new TH2F("Inv","Inv",37,0,37,37,0,37);
    TMatrixD* TotalCovarianceMatrix = new TMatrixD(37,37);

        for(Int_t i = 0; i<37; i++)
        {
            for(Int_t j = 0; j<37; j++)
            {
                TotalCovarianceMatrixM[(j+37*i)]=TotalCovarianceMatrixH->GetBinContent(i+1,j+1);
                std::cout << "TOTAL COV BY INVERSIONCHECK.h: " << TotalCovarianceMatrixM[(j+37*i)] << std::endl;

            }
        }
    

    TotalCovarianceMatrix->SetMatrixArray(TotalCovarianceMatrixM.GetArray());

//    Double_t* TotalCovarianceMatrixArray = CombineMatrices();
    //        TotalCovarianceMatrix->SetMatrixArray(TotalCovarianceMatrixArray);
    TotalCovarianceMatrix->Invert();
    
    Double_t* InvTotalCovarianceMatrixArray = TotalCovarianceMatrix->GetMatrixArray();
    
    
    for (Int_t i = 0; i < 37; i++)
    {
        for (Int_t j = 0; j < 37; j++)
        {
            InvTotalCovarianceMatrixM[i][j][0]=InvTotalCovarianceMatrixArray[i*37+j];
            InvTotalCovarianceMatrixH[0]->SetBinContent(i+1,j+1,InvTotalCovarianceMatrixM[i][j][0]);
            
//            std::cout << "INV BY INVCHECK.H: " << InvTotalCovarianceMatrixM[i][j][0] << std::endl;

        }
    }
   InvTotalCovarianceMatrixH[0]->Draw("colz");
//
    TH2F* UnityH[1];
    UnityH[0]=(TH2F*)TotalCovarianceMatrixH->Clone("Unity");
    UnityH[0]->Reset();
    Double_t UnityM[37*37];
    for(Int_t i=0; i<37; i++)
    {
        for(Int_t j=0; j<37; j++)
        {
            UnityM[(i+37*(j))] = 0;
            for(Int_t k=0; k<37; k++)
            {
                UnityM[i*37+j] += TotalCovarianceMatrixM[(i*37+k)] * InvTotalCovarianceMatrixM[k][j][0];
            }
        }
    }
    
    for (Int_t i = 0; i < 37; i++)
    {
        for (Int_t j = 0; j < 37; j++)
        {
            UnityH[0]->SetBinContent(i+1,j+1,UnityM[(i*37+j)]);
        }
    }

//    UnityH[0]->Draw("colz");

    delete TotalCovarianceMatrix;
}

//Double_t* CombineMatrices()
//{
//    Double_t TotalCovarianceMatrixM[37][37][1];
//
//    for(Int_t k = 0; k<1; k++)
//    {
//        for(Int_t i = 0; i<37; i++)
//        {
//            for(Int_t j = 0; j<37; j++)
//            {
//                if(i==j)
//                {
//                    TotalCovarianceMatrixM[i][j][k]=1;
//                }
//                else
//                {
//                    TotalCovarianceMatrixM[i][j][k]=0;
//                }
//            }
//        }
//    }
//    return &TotalCovarianceMatrixM[0][0][0];
//}


