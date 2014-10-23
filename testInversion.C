#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2.h"
#include "TMath.h"
#include <math.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TCanvas.h"
#include "TAxis.h"
#include "TArrayD.h"

int testInversion(void)
{
    Int_t n_evis_bins = 37;
    const Int_t MaxBins = 9*37;
     Int_t MaxNearCombine = 3;
     Int_t MaxFarCombine = 3;
    Double_t TotalCovarianceMatrixM[MaxBins*MaxBins];
    Double_t InvTotalCovarianceMatrixM[MaxBins][MaxBins];
    TH2D* InvTotalCovarianceMatrixH = new TH2D("Inv Total Covariance Matrix","Inv Total Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    TH2D* TotalCovarianceMatrixH = new TH2D(" Total Covariance Matrix"," Total Covariance Matrix",MaxBins, 0,MaxBins,MaxBins, 0,MaxBins);
    
    TMatrixD* TotalCovarianceMatrix = new TMatrixD(MaxBins,MaxBins);
    
    ifstream myfile("/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/CovarianceMatrices/Combine0/TotalCovarianceMatrix.txt");
    
    if(!myfile)
    {
        exit(EXIT_FAILURE);
    }
    
    Double_t current;
    Double_t Values[MaxBins*MaxBins];
    Int_t count = 0;
    setprecision(20);
    
    while(myfile >> current)
    {
         Values[count] = current;
         count++;
    }
    
    std::cout << "count is" << count << std::endl;
    
    Int_t x =0;
    Int_t y =0;
    
    Int_t Ni1=0,Ni2=0,Ni3=0,Ni4=0;
    Int_t Nj1=0,Nj2=0,Nj3=0,Nj4=0;
    Int_t Fi1=0,Fi2=0,Fi3=0,Fi4=0;
    Int_t Fj1=0,Fj2=0,Fj3=0,Fj4=0;
    
    for (Int_t neari=0; neari<MaxNearCombine; neari++)
    {
        //Logic for the 2D matrix index done up to 8 ADs
        if(neari==0){Ni1=1;Ni2=0;Ni3=0;Ni4=0;}
        if(neari==1){Ni2++;}
        if(neari==2){Ni3++;}
        if(neari==3){Ni4++;}
        
        for (Int_t fari=0; fari<MaxFarCombine; fari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(Ni1!=Ni2){Fi1=fari+1;Fi2=0;Fi3=0;Fi4=0;}
            if(Ni1==Ni2&&Ni2!=Ni3){Fi1=MaxFarCombine;Fi2=fari+1;}
            if(Ni2==Ni3&&Ni3!=Ni4){Fi2=MaxFarCombine;Fi3=fari+1;}
            if(Ni3==Ni4&&Ni4==1){Fi3=MaxFarCombine;Fi4=fari+1;}
            
            for (Int_t nearj=0; nearj<MaxNearCombine; nearj++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                if(nearj==1){Nj2++;}
                if(nearj==2){Nj3++;}
                if(nearj==3){Nj4++;}
                
                for (Int_t farj=0; farj<MaxFarCombine; farj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(Nj1!=Nj2){Fj1=farj+1;Fj2=0;Fj3=0;Fj4=0;}
                    if(Nj1==Nj2&&Nj2!=Nj3){Fj1=MaxFarCombine; Fj2=farj+1;}
                    if(Nj2==Nj3&&Nj3!=Nj4){Fj2=MaxFarCombine; Fj3=farj+1;}
                    if(Nj3==Nj4&&Nj4==1){Fj3=MaxFarCombine; Fj4=farj+1;}
                    
                    for (Int_t i = 0; i<n_evis_bins; i++)
                    {//columns
                        x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                        
                        for (Int_t j = 0; j<n_evis_bins; j++)
                        {//rows
                            y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                            
                            TotalCovarianceMatrixM[(j+MaxBins*i)] = Values[(j+MaxBins*i)];
                            
                            TotalCovarianceMatrixH->SetBinContent(x+1,y+1,TotalCovarianceMatrixM[(j+MaxBins*i)]);
                        }
                    }
                }
            }
        }
    }
    
    TotalCovarianceMatrix->SetMatrixArray(&TotalCovarianceMatrixM[0]);
    
    TotalCovarianceMatrix->Invert();
    
    Double_t* InvTotalCovarianceMatrixArray = TotalCovarianceMatrix->GetMatrixArray();
    
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j < MaxBins; j++)
        {
            InvTotalCovarianceMatrixM[i][j]=InvTotalCovarianceMatrixArray[(i*MaxBins+j)];
            InvTotalCovarianceMatrixH->SetBinContent(i+1,j+1,InvTotalCovarianceMatrixM[i][j]);
        }
    }
    
    UnityH=(TH2D*)TotalCovarianceMatrixH->Clone("Unity");
    UnityH->Reset();
    Double_t UnityM[MaxBins*MaxBins];
    
    for(Int_t i=0; i<MaxBins; i++)
    {
        for(Int_t j=0; j<MaxBins; j++)
        {
            UnityM[(j+MaxBins*i)] = 0;
            
            for(Int_t k=0; k<MaxBins; k++)
            {
                UnityM[(j+MaxBins*i)] += TotalCovarianceMatrixM[(k+MaxBins*i)] * InvTotalCovarianceMatrixM[k][j];
            }
        }
    }
    
    for (Int_t i = 0; i < MaxBins; i++)
    {
        for (Int_t j = 0; j < MaxBins; j++)
        {
            UnityH->SetBinContent(i+1,j+1,UnityM[(j+MaxBins*i)]);
        }
    }
    TCanvas* c1 = new TCanvas("c1","c1");
    c1->Divide(3,1);
    c1->cd(1);
    TotalCovarianceMatrixH->Draw("colz");
    c1->cd(2);
    InvTotalCovarianceMatrixH->Draw("colz");
    c1->cd(3);
    UnityH->Draw("colz");
    c1->Update();

}