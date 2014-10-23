#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TMath.h"

#include "NominalData.h"
#include <math.h>
#include "EnergyMatrix1.h"

const Int_t MaxPeriods = 20;
const Int_t MaxDetectors = 8;
const Int_t MaxNearDetectors =4;
const Int_t MaxFarDetectors =4;
const Int_t MaxNbins=51;
const Int_t Halls=3;
const bool AddBackgrounds=1;//To add backgrounds
const bool VariateRate = 1; // To vary backgrounds rate
const bool DistortBackgrounds=1;//To vary backgrounds shape
const bool WriteOutput=0;//To save the covariance matrices in a .txt file.
const Int_t NSamples = 1;//Number of samples used in the generation of covariance matrices

class CovarianceMatrix
{
private:
    
    TRandom3* rand;
    EnergyMatrix1* EnergyM;
    //AD configuration parameters:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    Int_t hall;

    //Binning parameters:
    Int_t Nbins;
    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Int_t Nweeks;
    Int_t InitialBins;
    //Relative errors:
    Double_t HAccidentalError[MaxDetectors];
    Double_t HLiHeError[Halls];
    Double_t HFastNeutronsError[Halls];
    Double_t HAmCError[Halls];
    Double_t DistortAcc;
    Double_t DistortLiHe;
    Double_t DistortFN;
    Double_t DistortAmC;

    Double_t ScaleFactorAccidental[MaxDetectors][MaxPeriods];
    Double_t ScaleFactorLiHe[MaxDetectors][MaxPeriods];
    Double_t ScaleFactorFastNeutrons[MaxDetectors][MaxPeriods];
    Double_t ScaleFactorAmC[MaxDetectors][MaxPeriods];

    //Choose matrix
    bool AccidentalMatrix;
    bool LiHeMatrix;
    bool FastNeutronsMatrix;
    bool AmCMatrix;
    
    char* OutputFileName;
    
    enum MatrixCalculation {AccidentalE, LiHeE, FastNeutronsE, AmCE};
    MatrixCalculation BackgroundE;

    TH1F* PredictionH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    TH1F* RandomPredictionH[MaxNearDetectors][MaxNearDetectors][MaxPeriods];
//    TH1F* RandomPredictionH1[MaxNearDetectors][MaxNearDetectors][MaxPeriods];

    TH1F* NearHallSpectrumH[MaxNearDetectors][MaxPeriods];
    TH1F* RandomNearHallSpectrumH[MaxNearDetectors][MaxPeriods];
    TH1F* BackgroundSpectrumH[MaxDetectors][MaxPeriods];
    TH1F* AccidentalsH[MaxDetectors][MaxPeriods];
    TH1F* LiHeH[MaxDetectors][MaxPeriods];
    TH1F* FastNeutronsH[MaxDetectors][MaxPeriods];
    TH1F* AmCH[MaxDetectors][MaxPeriods];
    TH1F* RandomAccidentalsH[MaxDetectors][MaxPeriods];
    TH1F* RandomLiHeH[MaxDetectors][MaxPeriods];
    TH1F* RandomFastNeutronsH[MaxDetectors][MaxPeriods];
    TH1F* RandomAmCH[MaxDetectors][MaxPeriods];
    TH2F* CovMatrix2H;
    TH2F* StatCov2H;
    TH2F* BkgCov2H;
    TH2F* EnergyMatrixH[MaxDetectors];
//    TH1F* TestGaussianH;

    Double_t Sigma_Near[MaxFarDetectors][MaxNearDetectors][MaxPeriods][MaxNbins];
    Double_t Sigma_Far[MaxFarDetectors][MaxNearDetectors][MaxPeriods][MaxNbins];
    Double_t CovStat[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t NormCovStat[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t CovBkg[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t NormCovBkg[9*MaxNbins][9*MaxNbins][MaxPeriods];

    //Main Data loaded from Theta13-inputs_20week.txt if 20 weeks is used. Otherwise change input file.
    Int_t ObsEvts[MaxDetectors][MaxPeriods];
    Double_t Livetime[MaxDetectors][MaxPeriods];
    Double_t MuonVetoEff[MaxDetectors][MaxPeriods];
    Double_t DMCEff[MaxDetectors][MaxPeriods];
    Double_t TargetMass[MaxDetectors][MaxPeriods];
    Double_t BgEvts[MaxDetectors][MaxPeriods];
    Double_t AccEvts[MaxDetectors][MaxPeriods];
    Double_t AccErr[MaxDetectors][MaxPeriods];
    Double_t Li9Evts[MaxDetectors][MaxPeriods];
    Double_t Li9Err[MaxDetectors][MaxPeriods];
    Double_t FnEvts[MaxDetectors][MaxPeriods];
    Double_t FnErr[MaxDetectors][MaxPeriods];
    Double_t AmcEvts[MaxDetectors][MaxPeriods];
    Double_t AmcErr[MaxDetectors][MaxPeriods];
    Double_t AlnEvts[MaxDetectors][MaxPeriods];
    Double_t AlnErr[MaxDetectors][MaxPeriods];
    Double_t ErrEvts[MaxDetectors][MaxPeriods];
    
    void LoadPredictions();
    void LoadBackgrounds();
    void FluctuateBackgrounds();
    void LoadNearHall();
    void LoadEnergyMatrix();
    TH1F* RebinHistogram(TH1F* Histogram);
    //Functions to vary background shapes. So far the same ones than LBNL.
    TF1* GetDistortionFunction(Double_t amount);
    TF1* GetFastNeutronsDistortionFunction(Double_t amount);
    
public:
    
    CovarianceMatrix();
    CovarianceMatrix(NominalData Data);

    void SetAccidentalMatrix(bool matrix);
    void SetLiHeMatrix(bool matrix);
    void SetFastNeutronsMatrix(bool matrix);
    void SetAmCMatrix(bool matrix);
    
    void GenerateStatisticalCovarianceMatrix();
    void GenerateCovarianceMatrix();
    void CovarianceMatrixMain();
    
};
CovarianceMatrix :: CovarianceMatrix()
{
    Nom = new NominalData();
    rand = new TRandom3();
    EnergyM = new EnergyMatrix1(Nom);

    Nbins = Nom.GetNbins();
    InitialEnergy = Nom.GetEmin();
    FinalEnergy = Nom.GetEmax();
    isH = Nom.GetAnalysis();
    InitialBins = Nom.GetInitialBins();
    Nweeks = Nom.GetWeeks();
    
    NADs = Nom.GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    if(NADs == 8) //ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    BackgroundE=4;//A number greater than the amount of different backgrounds included in the model (Right now Acc, LiHe, FN and AmC are indexes 0,1,2,3 respectively, if C(α,n) is added use 5);
    DistortAcc = 0;
    DistortLiHe= 0;
    DistortFN  = 0;
    DistortAmC = 0;
    
    //To draw using a better palette:
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;
    
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    gStyle->SetOptStat(0);
}

CovarianceMatrix :: CovarianceMatrix(NominalData Data, EnergyMatrix1 Energym)
{
    rand = new TRandom3();
    
    EnergyM = Energym;
    
    Nbins = Data.GetNbins();
    InitialEnergy = Data.GetEmin();
    FinalEnergy = Data.GetEmax();
    InitialBins = Data.GetInitialBins();
    Nweeks = Data.GetWeeks();
    
    NADs = Data.GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    if(NADs == 8) //ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    for(Int_t i=0;i<NADs;i++)
    {
        HAccidentalError[i]=Data.GetHAccidentalError(i);
    }
    
    for(Int_t i=0;i<3;i++)
    {
        HLiHeError[i]=Data.GetHLiHeError(i);
    }
    
    for(Int_t i=0;i<3;i++)
    {
        HFastNeutronsError[i]=Data.GetHFNError(i);
    }

    for(Int_t i=0;i<3;i++)
    {
        HAmCError[i]=Data.GetHAmCError(i);
    }
    
    BackgroundE=4;//A number greater than the amount of different backgrounds included in the model (Right now Acc, LiHe, FN and AmC are indexes 0,1,2,3 respectively, if C(α,n) is added use 5);
    DistortAcc = 0;
    DistortLiHe= 0;
    DistortFN  = 0;
    DistortAmC = 0;
    
    //To draw a better palette:
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;
    
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    gStyle->SetOptStat(0);
}

void CovarianceMatrix :: CovarianceMatrixMain()
{
    cout << "INITIAL BINS IS " << InitialBins << endl;
    LoadBackgrounds();
    LoadPredictions();
    LoadNearHall();
    GenerateStatisticalCovarianceMatrix();

    CovMatrix2H=(TH2F*)StatCov2H->Clone("Covariance Matrix");
    CovMatrix2H->Reset();
    CovMatrix2H->SetDirectory(0);
    
//    TestGaussianH = new TH1F("Test Gaussian","Test Gaussian",100,-50,50);
//    TestGaussianH->Reset();
//    TestGaussianH->SetDirectory(0);
    
    for (Int_t samples = 0; samples<NSamples; samples++)
    {
        FluctuateBackgrounds();
        GenerateCovarianceMatrix();
        CovMatrix2H->Add(BkgCov2H);
    }
    
    CovMatrix2H->Scale(1./(NSamples));
   // CovMatrix2H->Draw("colz");
    SaveCovarianceMatrix();
}

void CovarianceMatrix :: GenerateCovarianceMatrix()
{
    cout<<"ENERGY CALL"<<endl;
    EnergyM->EnergyMatrixMain();
    cout <<"ENERGY OUT"<<endl;
    LoadEnergyMatrix();
    cout << "Energy load"<<endl;
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        //Add nominal backgrounds
        for(Int_t AD=0; AD<NADs; AD++)
        {
            BackgroundSpectrumH[AD][week]=(TH1F*)AccidentalsH[0][0]->Clone();
            BackgroundSpectrumH[AD][week]->Reset();
            BackgroundSpectrumH[AD][week]->SetDirectory(0);

            if(AddBackgrounds)//Nominal backgrounds for Fpred
            {
                BackgroundSpectrumH[AD][week]->Add(AccidentalsH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(LiHeH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(FastNeutronsH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(AmCH[AD][week]);
            }
        }
        for (Int_t near=0; near<(ADsEH1+ADsEH2); near++)
        {
            cout << " NEAR " << near << endl;
            RandomNearHallSpectrumH[near][week]=(TH1F*)NearHallSpectrumH[near][week]->Clone();
            RandomNearHallSpectrumH[near][week]->SetDirectory(0);
            
//            EnergyMatrixH[near] = (TH2F*)EnergyM->GetEnergyMatrix(near);
//            cout << " LOAD SUCCESSFUL " << endl;
//            RandomNearHallSpectrumH[near][week]->Reset();
            
//            //Add detector effects //Multiply Matrix and RandomNearHallSpectrum
//
//            for(Int_t i=0;i<Nbins;i++)
//            {
//                Int_t Sum = 0;
//
//                for(Int_t j=0; j<Nbins; j++)
//                {
//                    Sum = Sum+(EnergyMatrixH[near]->GetBinContent(i+1,j+1)*NearHallSpectrumH[near][week]->GetBinContent(i+1));//I don't use RandomNearHallSpectrumH so I don't reset it
//                    
//                    cout << "SUM " << Sum << endl;
//                    cout << "i, j" << i << j << endl;
//                }
//                
//                RandomNearHallSpectrumH[near][week]->SetBinContent(i+1,Sum);
//                cout << "SetBinContent" << endl;
//
//            }
            
            if(AddBackgrounds)//Add fluctuated backgrounds for Fobs
            {
                cout << "HOLA" << endl;

                RandomNearHallSpectrumH[near][week]->Add(RandomAccidentalsH[near][week]);
                cout << "HOLA2" << endl;

                RandomNearHallSpectrumH[near][week]->Add(RandomLiHeH[near][week]);
                RandomNearHallSpectrumH[near][week]->Add(RandomFastNeutronsH[near][week]);
                RandomNearHallSpectrumH[near][week]->Add(RandomAmCH[near][week]);
            }
//            //Add detector effects //Multiply Matrix and RandomNearHallSpectrum
//
            for (Int_t far=0; far<1; far++)//ADsEH3
            {
//                EnergyMatrixH[far] = EnergyM->GetEnergyMatrix(far);//(TH2F*)
                
                TFile* EnergyF = TFile::Open("EnergyTest.root","recreate");
                EnergyMatrixH[far]->Write();
                EnergyF->Close();
                
                cout << " LOAD SUCCESSFUL " << endl;

                RandomPredictionH[far][near][week]=(TH1F*)PredictionH[far][near][week]->Clone();
                RandomPredictionH[far][near][week]->SetDirectory(0);
                RandomPredictionH[far][near][week]->Reset();

//                RandomPredictionH1[far][near][week]=(TH1F*)PredictionH[far][near][week]->Clone();
//                RandomPredictionH1[far][near][week]->SetDirectory(0);
                
                cout << "VEDIAMO" << endl;
                
                for(Int_t i=0;i<Nbins;i++)
                {
                    Int_t Sum = 0;
                    
                    for(Int_t j=0; j<Nbins; j++)
                    {
                        Sum = Sum+(PredictionH[far][near][week]->GetBinContent(i+1));//I don't use RandomNearHallSpectrumH so I don't reset it
                       // Sum = Sum+(EnergyMatrixH[far]->GetBinContent(i+1,j+1)*PredictionH[far][near][week]->GetBinContent(i+1));//I don't use RandomNearHallSpectrumH so I don't reset it

                        cout << "SUM " << Sum << endl;
                        cout << "i, j" << i << j << endl;
                    }
                    
                    RandomPredictionH[far][near][week]->SetBinContent(i+1,Sum);
                    cout << "SetBinContent" << endl;
                    
                }

                if(AddBackgrounds)
                {
                    RandomPredictionH[far][near][week]->Add(RandomAccidentalsH[ADsEH1+ADsEH2+far][week]);
                    RandomPredictionH[far][near][week]->Add(RandomLiHeH[ADsEH1+ADsEH2+far][week]);
                    RandomPredictionH[far][near][week]->Add(RandomFastNeutronsH[ADsEH1+ADsEH2+far][week]);
                    RandomPredictionH[far][near][week]->Add(RandomAmCH[ADsEH1+ADsEH2+far][week]);
                    cout << "HOLA3" << endl;

                    FluctuateBackgrounds();
                    cout << "HOLA4" << endl;

//                    RandomPredictionH1[far][near][week]->Add(RandomAccidentalsH[ADsEH1+ADsEH2+far][week]);
//                    RandomPredictionH1[far][near][week]->Add(RandomLiHeH[ADsEH1+ADsEH2+far][week]);
//                    RandomPredictionH1[far][near][week]->Add(RandomFastNeutronsH[ADsEH1+ADsEH2+far][week]);
//                    RandomPredictionH1[far][near][week]->Add(RandomAmCH[ADsEH1+ADsEH2+far][week]);
                }
                cout << "HOLA5" << endl;

            }
            cout << "HOLA6" << endl;
//            EnergyMatrixH[near]->Reset();

        }
        cout << "HOLA7" << endl;

    }
    cout << "HOLA 8" << endl;
    TH2F* ABkgCov2H = new TH2F("Background Cov Matrix","Background Covariance Matrix",Nbins*9,0,Nbins*9,Nbins*9,0,Nbins*9);
    BkgCov2H=(TH2F*)ABkgCov2H->Clone("Bkg covariance Matrix");
    BkgCov2H->SetDirectory(0);

    delete ABkgCov2H;

    Int_t x =0;
    Int_t y =0;
//    Int_t icounter = 1;
//    Int_t jcounter = 1;

    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<ADsEH3; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
//                cout <<"\n"<< "Ni" << Ni1 << Ni2 << Ni3 << Ni4 << "\n";
//                cout << "Fi"<< Fi1 << Fi2 << Fi3 << Fi4 << "\n";
//                cout << "====================================" << icounter++ <<"\n";
//                cout << Ni1*Fi1 << Ni2*Fi2 << Ni3*Fi3 << Ni4*Fi4;
                
                for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}

                    for (Int_t farj=0; farj<ADsEH3; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
//                        cout <<"\n"<< "Nj"<< Nj1 << Nj2 << Nj3 << Nj4 << "\n";
//                        cout << "Fj"<< Fj1 << Fj2 << Fj3 << Fj4 << "\n";
//                        cout << "===================================="<< jcounter++ <<"\n";
//                        cout << Nj1*Fj1 << Nj2*Fj2 << Nj3*Fj3 << Nj4*Fj4 <<"\n";
//                        cout <<"\n";
                        
                        for (Int_t i = 0; i<Nbins; i++)
                        {//columns
                            //cout << "x = " << i << "+ ("<< Ni1*Fi1<< "+"<< Ni2*Fi2 << "+"<< Ni3*Fi3 << "+"<< Ni4*Fi4 << "-1)*41"<< " that is "<< i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins << "\n" ;

                            for (Int_t j = 0; j<Nbins; j++)
                            {//rows
                                x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*Nbins;
                               
                                CovBkg[x][y][week]=(RandomPredictionH[fari][neari][week]->GetBinContent(i+1)-PredictionH[fari][neari][week]->GetBinContent(i+1)-BackgroundSpectrumH[fari+ADsEH1+ADsEH2][week]->GetBinContent(i+1))*(RandomPredictionH[farj][nearj][week]->GetBinContent(j+1)-PredictionH[farj][nearj][week]->GetBinContent(j+1)-BackgroundSpectrumH[farj+ADsEH1+ADsEH2][week]->GetBinContent(j+1));
                                
                              //  cout << CovBkg[x][y][week] << "\n";
                            }
                          //  TestGaussianH->Fill(RandomPredictionH[fari][neari][week]->GetBinContent(i+1)-RandomPredictionH1[fari][neari][week]->GetBinContent(i+1));
                          //  TestGaussianH->Fill(RandomPredictionH[fari][neari][week]->GetBinContent(i+1)-PredictionH[fari][neari][week]->GetBinContent(i+1)-BackgroundSpectrumH[fari+ADsEH1+ADsEH2][week]->GetBinContent(i+1));

                        }
                    }
                }
            }
        }
    }
   // TestGaussianH->Draw();
    Int_t x =0;
    Int_t y =0;
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<ADsEH3; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<ADsEH3; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
                        for (Int_t i = 0; i<Nbins; i++)
                        {//columns
                            for (Int_t j = 0; j<Nbins; j++)
                            {//rows
                                x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*Nbins;
                                
                                NormCovBkg[x][y][week]=(CovBkg[x][y][week])/(sqrt(CovBkg[x][x][week]*CovBkg[y][y][week]));
                                if(CovBkg[x][y][week]==0&&((CovBkg[x][x][week]||CovBkg[y][y][week])==0))
                                {
                                    NormCovBkg[x][y][week]=0;//To avoid (0/0) nans when bin contents are practically the same; (this happens when backgrounds are not varied)
                                }
                                BkgCov2H->SetBinContent(x+1,y+1,NormCovBkg[x][y][week]);
                              //  cout <<NormCovBkg[x][y][week];
                               // cout << " " << CovBkg[x][y][week] << " " <<  CovBkg[x][x][week] << " " << CovBkg[y][y][week] << " " << sqrt(CovBkg[x][x][week]*CovBkg[y][y][week]) << "\n";

                            }
                        }
                    }
                }
            }
        }
    }
    
//    BkgCov2H->Draw("colz");
    EnergyM->~EnergyMatrix1();
}

////////////////////////////////////////////////////////////////////////////////////
//// I need to correct backgrounds for efficiencies, and maybe events too (Check)
////////////////////////////////////////////////////////////////////////////////////
void CovarianceMatrix :: GenerateStatisticalCovarianceMatrix()
{
    for (Int_t week = 0; week<Nweeks; week++)
    {
        //Add nominal backgrounds
        for(Int_t AD=0; AD<NADs; AD++)
        {
            BackgroundSpectrumH[AD][week]=(TH1F*)AccidentalsH[0][0]->Clone();
            BackgroundSpectrumH[AD][week]->Reset();
            BackgroundSpectrumH[AD][week]->SetDirectory(0);

            if(AddBackgrounds)
            {
                BackgroundSpectrumH[AD][week]->Add(AccidentalsH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(LiHeH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(FastNeutronsH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(AmCH[AD][week]);
            }
        }
        SaveTotalSpectrum();

            for (Int_t far=0; far<ADsEH3; far++)
            {
                for (Int_t near=0; near<(ADsEH1+ADsEH2); near++)
                {
                    for (Int_t pts = 0; pts < Nbins; pts++)
                    {
                        Sigma_Far[far][near][week][pts]=sqrt(PredictionH[far][near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[ADsEH1+ADsEH2+far][week]->GetBinContent(pts+1));
                        
                        Sigma_Near[far][near][week][pts]=(PredictionH[far][near][week]->GetBinContent(pts+1)/NearHallSpectrumH[near][week]->GetBinContent(pts+1))*sqrt(NearHallSpectrumH[near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[near][week]->GetBinContent(pts+1));
                    }
                }
            }
    }
    
    TH2F* AStatCov2H = new TH2F("Statistical Cov Matrix","Statistical Covariance Matrix",Nbins*9,0,Nbins*9,Nbins*9,0,Nbins*9);
    StatCov2H=(TH2F*)AStatCov2H->Clone("Statistical Covariance Matrix");
    StatCov2H->SetDirectory(0);
    delete AStatCov2H;
    Int_t x =0;
    Int_t y =0;
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<ADsEH3; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<ADsEH3; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
                        for (Int_t i = 0; i<Nbins; i++)
                        {//columns
                            for (Int_t j = 0; j<Nbins; j++)
                            {//rows
                                x= i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*Nbins;
                                
                                //Near component correlated
                                if(neari==nearj && fari!=farj)
                                {
                                    CovStat[x][y][week]=Sigma_Near[fari][neari][week][i]*Sigma_Near[farj][nearj][week][j];
                                }
                                //Far component correlated
                                if(fari==farj && neari!=nearj)
                                {
                                    CovStat[x][y][week]=Sigma_Far[fari][neari][week][i]*Sigma_Far[farj][nearj][week][j];
                                }
                                if(neari==nearj && fari==farj)
                                {
                                //General covariance
                                    CovStat[x][y][week]=(Sigma_Near[fari][neari][week][i]*Sigma_Near[fari][neari][week][j])+(Sigma_Far[farj][nearj][week][i]*Sigma_Far[farj][nearj][week][j]);
                                }
                                //Uncorrelated terms
                                if(neari!=nearj && fari!=farj)
                                {
                                    CovStat[x][y][week]=0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Int_t x =0;
    Int_t y =0;
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<ADsEH3; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<ADsEH3; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
                        for (Int_t i = 0; i<Nbins; i++)
                        {//columns
                            for (Int_t j = 0; j<Nbins; j++)
                            {//rows
                                x= i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*Nbins;
                                
                                NormCovStat[x][y][week]=(CovStat[x][y][week])/(sqrt(CovStat[x][x][week]*CovStat[y][y][week]));
                                StatCov2H->SetBinContent(x+1,y+1,NormCovStat[x][y][week]);
                            }
                        }
                    }
                }
            }
        }
    }
    
  // StatCov2H->Draw("colz");

    if(WriteOutput)
    {
        ofstream outf("CovarianceMatrices/StatisticalCovarianceMatrix.txt");
        Int_t x =0;
        Int_t y =0;
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
                if(neari==1){Ni2++;}
                if(neari==2){Ni3++;}
                if(neari==3){Ni4++;}
                
                for (Int_t fari=0; fari<ADsEH3; fari++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
                    if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
                    if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
                    if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
                    
                    for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
                        if(nearj==1){Nj2++;}
                        if(nearj==2){Nj3++;}
                        if(nearj==3){Nj4++;}
                        
                        for (Int_t farj=0; farj<ADsEH3; farj++)
                        {
                            //Logic for the 2D matrix index done up to 8 ADs
                            if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
                            if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
                            if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
                            if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
                            for (Int_t i = 0; i<Nbins; i++)
                            {//columns
                                for (Int_t j = 0; j<Nbins; j++)
                                {//rows
                                    x= i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins;
                                    y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*Nbins;
                                    
                                    outf << NormCovStat[x][y][week] << " ";
                                }
                            }
                        }
                    }
                }
            }
            outf << endl;
        }
        outf.close();
    }
}

void CovarianceMatrix::LoadBackgrounds()
{
    TFile* BackgroundsF = TFile::Open("./BackgroundSpectrum/Backgrounds.root");
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t AD =0; AD<NADs; AD++)
        {
            AccidentalsH[AD][week]= (TH1F*)gDirectory->Get(Form("Accidentals_AD%i",AD+1));
            AccidentalsH[AD][week]->SetDirectory(0);
            LiHeH[AD][week]=(TH1F*)gDirectory->Get("hist_Bkg_StrongAmC");//Missing LiHe inputs so far
            LiHeH[AD][week]->SetDirectory(0);
            FastNeutronsH[AD][week]=(TH1F*)gDirectory->Get("FN");
            FastNeutronsH[AD][week]->SetDirectory(0);
            AmCH[AD][week]=(TH1F*)gDirectory->Get("hist_Bkg_StrongAmC");
            AmCH[AD][week]->SetDirectory(0);
        }
    }
      BackgroundsF->Close();
}
//Randomize the events (rate) inside errors taking into account the corresponding correlations between backgrounds and ADs. Also random shape variations are included using distortion functions.
void CovarianceMatrix::FluctuateBackgrounds()
{

    for (Int_t week = 0; week<Nweeks; week++)
    {
            for (Int_t AD =0; AD<NADs; AD++)
            {
                ScaleFactorAccidental[AD][week]=1;
                ScaleFactorLiHe[AD][week]=1;
                ScaleFactorFastNeutrons[AD][week]=1;
                ScaleFactorAmC[AD][week]=1;
                
                hall=2;
                if(AD<ADsEH1)
                {
                    hall=1;
                }
                if(AD>=ADsEH1+ADsEH2)
                {
                    hall=3;
                }
                
                if(AccidentalMatrix)
                {
                    rand->SetSeed(0);
                    ScaleFactorAccidental[AD][week]=(1+HAccidentalError[AD]*rand->Gaus(0,1));
                }
                if(LiHeMatrix)
                {
                    rand->SetSeed(0);
                    ScaleFactorLiHe[AD][week]=(1+HLiHeError[hall]*rand->Gaus(0,1));
                }
                if(FastNeutronsMatrix)
                {
                    rand->SetSeed(0);
                    ScaleFactorFastNeutrons[AD][week]=(1+HFastNeutronsError[hall]*rand->Gaus(0,1));
                }
                if(AmCMatrix)
                {
                    rand->SetSeed(0);
                    ScaleFactorAmC[AD][week]=(1+HAmCError[hall]*rand->Gaus(0,1));
                }
            }
            
            //Accidentals uncorrelated
            //LiHe and Fast neutrons correlated by hall
            
            ScaleFactorLiHe[1][week]=ScaleFactorLiHe[0][week];
            ScaleFactorLiHe[4][week]=ScaleFactorLiHe[3][week];
            ScaleFactorLiHe[5][week]=ScaleFactorLiHe[3][week];
            ScaleFactorFastNeutrons[1][week]=ScaleFactorFastNeutrons[0][week];
            ScaleFactorFastNeutrons[4][week]=ScaleFactorFastNeutrons[3][week];
            ScaleFactorFastNeutrons[5][week]=ScaleFactorFastNeutrons[3][week];
            
            //AmC all correlated
            
            ScaleFactorAmC[1][week]=ScaleFactorAmC[0][week];
            ScaleFactorAmC[2][week]=ScaleFactorAmC[0][week];
            ScaleFactorAmC[3][week]=ScaleFactorAmC[0][week];
            ScaleFactorAmC[4][week]=ScaleFactorAmC[0][week];
            ScaleFactorAmC[5][week]=ScaleFactorAmC[0][week];
            
            for (Int_t AD =0; AD<NADs; AD++)
            {
                RandomAccidentalsH[AD][week]=(TH1F*)AccidentalsH[AD][week]->Clone();
                RandomLiHeH[AD][week]=(TH1F*)LiHeH[AD][week]->Clone();
                RandomFastNeutronsH[AD][week]=(TH1F*)FastNeutronsH[AD][week]->Clone();
                RandomAmCH[AD][week]=(TH1F*)AmCH[AD][week]->Clone();
                
                RandomLiHeH[AD][week]->SetDirectory(0);
                RandomAmCH[AD][week]->SetDirectory(0);
                RandomLiHeH[AD][week]->SetDirectory(0);
                RandomAmCH[AD][week]->SetDirectory(0);
                
                if(VariateRate)
                {
                    RandomAccidentalsH[AD][week]->Scale(ScaleFactorAccidental[AD][week]);
                    RandomLiHeH[AD][week]->Scale(ScaleFactorLiHe[AD][week]);
                    RandomFastNeutronsH[AD][week]->Scale(ScaleFactorFastNeutrons[AD][week]);
                    RandomAmCH[AD][week]->Scale(ScaleFactorAmC[AD][week]);
                   // RandomAmCH[AD][week]->Draw("same");
                    //            RandomLiHeH[AD][week]->Draw("same");
                     //           RandomFastNeutronsH[AD][week]->Draw("same");
                    // RandomAccidentalsH[AD][week]->Draw("same");

                }

            }

            if(DistortBackgrounds)
            {
                if(AccidentalMatrix)
                {
                    TF1* func_acc=GetDistortionFunction(DistortAcc);
                    
                    for(Int_t iAD=0;iAD<NADs;iAD++)
                    {
                        RandomAccidentalsH[iAD][week]->Multiply(func_acc);
                        RandomAccidentalsH[iAD][week]->Scale(RandomAccidentalsH[iAD][week]->Integral()/AccidentalsH[iAD][week]->Integral());
                        //      cout<<RandomAccidentalsH[AD][week]->Integral()/AccidentalsH[AD][week]->Integral()<<"\n";
                        //      func_acc->Draw();
                    }
//                    RandomAccidentalsH[0][week]->Draw("same");

                    delete func_acc;
                }
                if(LiHeMatrix)
                {
                    TF1* func_LiHe=GetDistortionFunction(DistortLiHe);
                    for(Int_t iAD=0;iAD<NADs;iAD++)
                    {
                        RandomLiHeH[iAD][week]->Multiply(func_LiHe);
                        RandomLiHeH[iAD][week]->Scale(RandomLiHeH[iAD][week]->Integral()/LiHeH[iAD][week]->Integral());
                    //                    cout<<RandomLiHeH[AD][week]->Integral()/LiHeH[AD][week]->Integral()<<"\n";
                    //      func_LiHe->Draw("same");
                    }
                    delete func_LiHe;
                }
                if(FastNeutronsMatrix)
                {

//FN shape distortions are applied taking into account correlations in each EH
                    TF1* func_FN=GetFastNeutronsDistortionFunction(DistortFN);

                        for(Int_t iAD=0;iAD<ADsEH1;iAD++)
                        {
                            RandomFastNeutronsH[iAD][week]->Multiply(func_FN);
                            RandomFastNeutronsH[iAD][week]->Scale(RandomFastNeutronsH[iAD][week]->Integral()/FastNeutronsH[iAD][week]->Integral());
                        }
                    func_FN=GetFastNeutronsDistortionFunction(DistortFN);

                        for(Int_t iAD=ADsEH1;iAD<ADsEH2+ADsEH1;iAD++)
                        {
                            RandomFastNeutronsH[iAD][week]->Multiply(func_FN);
                            RandomFastNeutronsH[iAD][week]->Scale(RandomFastNeutronsH[iAD][week]->Integral()/FastNeutronsH[iAD][week]->Integral());
                        }
                    func_FN=GetFastNeutronsDistortionFunction(DistortFN);

                        for(Int_t iAD=ADsEH1+ADsEH2;iAD<ADsEH1+ADsEH2+ADsEH3;iAD++)
                        {
                            RandomFastNeutronsH[iAD][week]->Multiply(func_FN);
                            RandomFastNeutronsH[iAD][week]->Scale(RandomFastNeutronsH[iAD][week]->Integral()/FastNeutronsH[iAD][week]->Integral());
                        }
                    
//                    for(Int_t iAD=0;iAD<NADs;iAD++)
//                    {
//                       RandomFastNeutronsH[iAD][week]->Draw("same");
//                    }
                    //  func_FN->Draw("same");
                    //   RandomFastNeutronsH[0][week]->Draw("same");

                    delete func_FN;
                }
                if(AmCMatrix)
                {
                    TF1* func_AmC=GetDistortionFunction(DistortAmC);
                    for(Int_t iAD=0;iAD<NADs;iAD++)
                    {
                        RandomAmCH[iAD][week]->Multiply(func_AmC);
                        RandomAmCH[iAD][week]->Scale(RandomAmCH[iAD][week]->Integral()/AmCH[iAD][week]->Integral());
                       //RandomAmCH[iAD][week]->Draw("same");

                        //     func_AmC->Draw("same");
                    }

                    delete func_AmC;
                }
            }
    }
}

void CovarianceMatrix::LoadPredictions()
{
TFile* FarHallPredictionsF = TFile::Open("./RootOutputs/FarSpectrumFraction.root");
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            for (Int_t far =0; far<ADsEH3; far++)
            {
                PredictionH[far][near][week] = (TH1F*)gDirectory->Get(Form("AD%i Far Spectrum prediction from near AD%i",far+1,near+1));
                PredictionH[far][near][week]->SetDirectory(0);
            }
        }
    }
    FarHallPredictionsF->Close();
}

void CovarianceMatrix::LoadEnergyMatrix()
{
    
    TFile* EnergyMatrixF = TFile::Open("./RootOutputs/EnergyMatrix1.root");
    
        for (Int_t ad =0; ad<1; ad++)//ad<ADsEH3
        {
            EnergyMatrixH[ad] = (TH2F*)gDirectory->Get(Form("EvisEnu%i",ad));
            EnergyMatrixH[ad]->SetDirectory(0);
        }
    
    EnergyMatrixF->Close();
}

void CovarianceMatrix::LoadNearHall()
{
    Char_t filenameNear[1024];
    
    TFile* NearHallDataF = TFile::Open("./RootOutputs/NominalOutputs/Oscillation.root"); //This should be real data. It has to be fixed.
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            sprintf(filenameNear,"Total spectrum after oscillation at AD%i",near+1);
            NearHallDataF->cd("Total AD Spectra after oscillation");
            NearHallSpectrumH[near][week] = (TH1F*)gDirectory->Get(filenameNear);
            NearHallSpectrumH[near][week]->SetDirectory(0);
        }
    }
    NearHallDataF->Close();
}

void CovarianceMatrix::SaveTotalSpectrum()
{
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            NearHallSpectrumH[near][week]->Add(BackgroundSpectrumH[near][week]);

            for (Int_t far =0; far<ADsEH3; far++)
            {
                PredictionH[far][near][week]->Add(BackgroundSpectrumH[ADsEH1+ADsEH2+far][week]);
            }
        }
     }
    
TFile* SaveSpectrumDataF = TFile::Open("./RootOutputs/SpectrumWithBackgrounds.root","recreate");
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            NearHallSpectrumH[near][week]->Write();
            
            for (Int_t far =0; far<ADsEH3; far++)
            {
                PredictionH[far][near][week]->Write();
            }
        }
    }
   SaveSpectrumDataF->Close();
}

void CovarianceMatrix::SaveCovarianceMatrix()
{
    TFile* SaveCovarianceMatrixF = TFile::Open("./CovarianceMatrices/CovarianceMatrixMatrices.root","recreate");
    StatCov2H->Write();

    switch (BackgroundE)
    {
        case 1://Vary Accidentals
            CovMatrix2H->SetName("Accidental Covariance Matrix");
            CovMatrix2H->SetTitle("Accidental Covariance Matrix");
        case 2://Vary LiHe
            CovMatrix2H->SetName("LiHe Covariance Matrix");
            CovMatrix2H->SetTitle("LiHe Covariance Matrix");
        case 3://Vary Fast Neutrons
            CovMatrix2H->SetName("FN Covariance Matrix");
            CovMatrix2H->SetTitle("FN Covariance Matrix");
        case 4://Vary AmC
            CovMatrix2H->SetName("AmC Covariance Matrix");
            CovMatrix2H->SetTitle("AmC Covariance Matrix");
        default://Add nominal backgrounds
            CovMatrix2H->SetName("Nominal Covariance Matrix");
            CovMatrix2H->SetTitle("Nominal Covariance Matrix");
    }
   // CovMatrix2H->Draw("colz");
    CovMatrix2H->Write();
    SaveCovarianceMatrixF->Close();
    
    if(WriteOutput)
    {
        ofstream outf("CovarianceMatrices/BackgroundCovarianceMatrix.txt");
        Int_t x =0;
        Int_t y =0;
            for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
                if(neari==1){Ni2++;}
                if(neari==2){Ni3++;}
                if(neari==3){Ni4++;}
                
                for (Int_t fari=0; fari<ADsEH3; fari++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
                    if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
                    if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
                    if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
                    
                    for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
                        if(nearj==1){Nj2++;}
                        if(nearj==2){Nj3++;}
                        if(nearj==3){Nj4++;}
                        
                        for (Int_t farj=0; farj<ADsEH3; farj++)
                        {
                            //Logic for the 2D matrix index done up to 8 ADs
                            if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
                            if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
                            if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
                            if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
                            
                            for (Int_t i = 0; i<Nbins; i++)
                            {//columns
                                for (Int_t j = 0; j<Nbins; j++)
                                {//rows
                                    x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins;
                                    y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*Nbins;
                                    
                                    outf << CovMatrix2H->GetBinContent(x+1,y+1)<< " ";
                                }
                            }
                        }
                    }
                }
            }
            outf << endl;
            outf.close();
    }
}

TF1* CovarianceMatrix::GetDistortionFunction(Double_t amount)
{
    TF1 *func = new TF1("func","TMath::Abs([0]+[1]*x)",InitialEnergy,FinalEnergy);
    rand->SetSeed(0);
    Double_t slope=amount*rand->Gaus(0,1);
    Double_t anchor_point=3.5;
    //want offset to be set by requiring func at anchor point to be 1
    Double_t offset=(1-slope*anchor_point);
    func->SetParameter(0,offset);
    func->SetParameter(1,slope);
    
    return func;
}

TF1* CovarianceMatrix::GetFastNeutronsDistortionFunction(Double_t amount)
{
    TF1 *func = new TF1("func","[0]/(0.2*pow(x,0.1))+[1]",InitialEnergy,FinalEnergy);
    rand->SetSeed(0);
    Double_t scaling =amount*rand->Gaus(0,1);
    func->SetParameter(0,scaling);
// func->SetParameter(1,0); // see if it works
    //set offset so that func(FinalEnergy)=1;
    func->SetParameter(1,1-1*func->Eval(10));//I'm using 10 instead of 12 because when I fit the FN I use the range 10-100 MeV.
    
    return func;
}

void CovarianceMatrix :: SetAccidentalMatrix(bool AccMatrix)
{
    AccidentalMatrix=AccMatrix;
    if(AccMatrix)
    {
        BackgroundE = AccidentalE;
        DistortAcc=0.2;
    }
}

void CovarianceMatrix :: SetLiHeMatrix(bool LiMatrix)
{
    LiHeMatrix=LiMatrix;
    if(LiMatrix)
    {
        BackgroundE = LiHeE;
        DistortLiHe=0.2;
    }
}

void CovarianceMatrix :: SetFastNeutronsMatrix(bool FNmatrix)
{
    FastNeutronsMatrix=FNmatrix;
    if(FNmatrix)
    {
        BackgroundE = FastNeutronsE;
        DistortFN=0.2;
    }
}

void CovarianceMatrix :: SetAmCMatrix(bool AmCmatrix)
{
    AmCMatrix=AmCmatrix;
    if(AmCMatrix)
    {
        BackgroundE = AmCE;
        DistortAmC=0.2;
    }
}

TH1F* CovarianceMatrix :: RebinHistogram(TH1F* Histogram)
{
    TH1F* NewHistogram = new TH1F(Histogram->GetName(),Histogram->GetNbins()+InitialBins,0,FinalEnergy);
    
    for(Int_t i = 1; i<=Histogram->GetNbins()+InitialBins; i++)
    {
        if(i<InitialBins+1)
        {
            NewHistogram->SetBinContent(i,0);
        }
        else
        {
            NewHistogram->SetBinContent(i,Histogram->GetBinContent(i-InitialBins));
        }
    }
    
    return NewHistogram;
}