#pragma once
#include "NominalData.h"
#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include <TMatrixD.h>
#include <vector>
#include "TRandom3.h"
#include "TTree.h"

#define UseLiHeToyMC
// Originally: Given a reactor antineutrino spectrum I output the respective spectra for each AD after oscillation.
// Updated: Class that is able to calculate the spectra in each AD either from Reactor Data or Near Hall Data. 5/3/13.
//
//  Created by Jose de Arcos	 on 3/10/13.
//
//

// Last version:

// LOAD TOY MC FROM OSCILLATION REACTOR OR LOAD DATA

// THEN MULTIPLY BY ROW RESPONSE MATRIX TO GET TRUE ENERGY

// THEN CORRECT FOR AD EFFICIENCY IN TRUE ENERGY

//        TotalOscillatedSpectrumAD[AD]->Scale(1./(DetectorProtons[AD]*DetectorEfficiency[AD][week]*FullTime[AD][week]));//Correct for efficiencies in each AD divide by # protons.

// THEN OSCILLATE BACK AND FORTH

// IN PREDICTION I WILL INTEGRATE BACK IN VISIBLE ENERGY

//#define CheckSumTrueSpectrum
class Oscillation
{
private:
    
    //Using Pedro's MC for LiHe:
#ifdef UseLiHeToyMC
    TH1F* func_LiHe;
#else
    TF1* func_LiHe;
#endif
    std::string RandomString;
    //External classes used
    NominalData* Nom;
    //Survival probability calculation parameters:
    Double_t Energy;
    Int_t hierarchy;
    Double_t deltam2_32;
    Double_t deltam2_21;
    Double_t deltam2_31;
    Double_t dm2ee;
    
    Int_t WEEK;
    TRandom3* rand;

    bool IsotopeMatrix;
    bool ReactorPowerMatrix;
    bool Sin22t12Matrix;
    
    TF1* NominalAmCF;
    bool VaryAccidentalMatrix;
    bool VaryLiHeMatrix;
    bool VaryFastNeutronsMatrix;
    bool VaryAmCMatrix;
    bool DistortLiHeMatrix;
    bool DistortFastNeutronsMatrix;
    bool DistortAmCMatrix;
    
    bool Mode;
    
    std::string ResponseDirectory;

    Double_t s22t12;
    Double_t s22t13;
    
    Int_t week;
    std::vector<Double_t> ADdistances;
    std::vector<Double_t> Norma;
    
    std::vector<Double_t> FluxFraction;
    std::vector<Double_t> Extrapolation;
    
    std::vector<Double_t> DetectorProtons;
    Double_t FullTime[MaxDetectors][MaxPeriods];
    Double_t DetectorEfficiencyDelayed[MaxDetectors];
    Double_t DetectorEfficiency[MaxDetectors][MaxPeriods];
    Double_t NominalDetectorEfficiency[MaxDetectors][MaxPeriods];
    
    //AD configuration parameters:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    
    //For the weekly reactor inputs
    Int_t Nweeks;
    Int_t NReactorPeriods;
    //For the kind of data
    bool Analysis;
    std::string AnalysisString;
    Int_t DataSet;
    //Binning parameters:
    Int_t n_evis_bins;
    Int_t n_etrue_bins;
    Double_t InitialEnergy;
    Double_t FinalVisibleEnergy;
    Double_t InitialVisibleEnergy;
    
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.
    
    //Background errors:
    Double_t AccidentalError[MaxDetectors][MaxPeriods];
    Double_t LiHeError[MaxDetectors][MaxPeriods];
    Double_t FastNeutronsError[MaxDetectors][MaxPeriods];
    Double_t AmCError[MaxDetectors][MaxPeriods];
    
    //Distortion function factor:
    Double_t DistortLiHe;
    Double_t DistortFN;
    Double_t DistortAmC;
    
    Double_t ScaleFactorAccidental[MaxDetectors];
    Double_t ScaleFactorLiHe[MaxDetectors];
    Double_t ScaleFactorFastNeutrons[MaxDetectors];
    Double_t ScaleFactorAmC[MaxDetectors];
    
    //Histograms
    
    TH1D* BackgroundSpectrumH[MaxDetectors];
    TH1D* NearBackgroundSpectrumH[MaxNearDetectors];
    TH1D* FarBackgroundSpectrumH[MaxDetectors];
    
    TH1D* AccidentalsH[MaxDetectors];
    TH1D* LiHeH[MaxDetectors];
    TH1D* FastNeutronsH[MaxDetectors];
    TH1D* AmCH[MaxDetectors];
    
    TH1D* RandomAccidentalsH[MaxDetectors];
    TH1D* RandomLiHeH[MaxDetectors];
    TH1D* RandomFastNeutronsH[MaxDetectors];
    TH1D* RandomAmCH[MaxDetectors];

    TH1D* NearSpectrumFractionH[NReactors*MaxNearDetectors][MatrixBins];//
    TH1D* NearSpectrumH[MaxNearDetectors][MatrixBins];//
    TH1D* FarHallSpectrumFractionH[MaxFarDetectors*MaxNearDetectors*NReactors][MatrixBins];//
    TH1D* FarHallSpectrumH[MaxFarDetectors*MaxNearDetectors][MatrixBins];//
    
    std::vector<TH1D*> FluxFractionH;//
    std::vector<TH1D*> ExtrapolationH;//
    
    TH1D* TrueADSpectrumH[MaxNearDetectors][MatrixBins];//
    //    std::vector<TH1D*> CheckVisADSpectrumH;
    TH1D* ADSpectrumVisH[MaxNearDetectors];//
    TH1D* FluxH[NReactors][MaxDetectors][MaxPeriods];//
    TH1D* InclusiveFluxH[NReactors][MaxDetectors];//
    std::vector<TH1D*> OscillationFunction;//
    TH1D* FDBH;//
    TH1D* FLOH;//
    
    TH2D* TransEnergyMatrix;
    TH2D* EnergyMatrix;
    TH2D* LBNLEnergyMatrix;
    //Files
    Char_t DistanceFileName[100];
    Char_t name[20];
    const Double_t colors[6]= {1,2,4,8,6,11};
    const std::string labels[6]={"D1","D2","L1","L2","L3","L4"};

    bool flagDelete;
    //Load data
    void ReadDistances(Char_t*);
    void LoadToyHistograms(Int_t);
    void LoadNearData(Int_t);
    
    //Functions to vary background shapes. So far the same ones than LBNL.
#ifdef UseLiHeToyMC
    void GetDistortionFunction(Double_t,TH1F*);
#else
    void GetDistortionFunction(Double_t,TF1*);
#endif
    void GetFastNeutronsDistortionFunction(Double_t,TF1*);
    void FluctuateBackgrounds(Int_t);
    void LoadNominalBackgrounds();
    void LoadBackgrounds(Int_t);
    void SaveBackgrounds();

    //Calculate oscillation
    void CalculateFluxFraction();
    void NearSpectrumFraction();
    void GetExtrapolation();
    void FarSpectrumPrediction();
    void SaveOscillationFunction();
    
    Double_t OscProb(Double_t, Double_t, Double_t,Double_t);
    Double_t PlotOscProb(Double_t*,Double_t*);
    void ApplyResponseMatrix();
public:
    Oscillation();
    Oscillation(NominalData*);
    ~Oscillation();
    void OscillationFromNearHallData(Int_t,bool,bool);
    void SetEfficiency(Double_t[]);
    void SetSin22t12(Double_t);
    void LoadFluxHisto();
    void GenerateFluxHisto();
    Double_t GetSin22t13();
    void SetOscillationParameters(Double_t,Double_t);
    TH1D* GetOscillatedADSpectrum(Int_t,Int_t,Int_t,Int_t);
    void SetDelayedEfficiency(Double_t[]);
    void SetFluxHistograms(TH1D*,Int_t,Int_t,Int_t);
    TH1D* GetFluxHisto(Int_t,Int_t,Int_t);
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//DEFAULT CONSTRUCTOR

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Default values
Oscillation :: Oscillation()
{
    std::cout << " oscillation the default constructor shouldn't be called" << std::endl;
    exit(EXIT_FAILURE);
    
    flagDelete=0;

    Nom = new NominalData(0,2);
    rand = new TRandom3(0);

    ResponseDirectory = Nom->GetResponseDirectory();

    hierarchy=Nom->GetHierarchy();
    deltam2_32=Nom->GetDm232();
    deltam2_21=Nom->GetDm221();
    dm2ee = Nom->GetDm2ee();
    s22t13 = Nom->GetSin22t13();
    s22t12 = Nom->GetSin22t12();
    
    Analysis = Nom->GetAnalysis();
    if(Analysis)
    {
        AnalysisString = "Hydrogen";
    }
    else
    {
        AnalysisString = "Gadolinium";
    }
    DataSet = Nom->GetDataSet();
    
    Nweeks = Nom->GetWeeks();
    NReactorPeriods = Nom->GetNReactorPeriods();

    NADs = Nom->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
#ifdef BlindedAnalysis
    sprintf(DistanceFileName,"./Distances/blinded_baseline.txt");
    std::cout << " NEED TO ADD A BLINDED DISTANCE FILE" << std::endl;
    exit(EXIT_FAILURE);
#else
    sprintf(DistanceFileName,"./Distances/unblinded_baseline.txt");
#endif
    
    if(NADs == 8) //ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
#ifdef BlindedAnalysis
        sprintf(DistanceFileName,"./Distances/blinded_baseline.txt");
        std::cout << " NEED TO ADD A BLINDED 8AD DISTANCE FILE" << std::endl;
        exit(EXIT_FAILURE);
#else
        sprintf(DistanceFileName,"./Distances/unblinded_baseline8ADs.txt");//change file to calculate distances to a file that has the information for the 8 ADs
        std::cout << "NEED TO ADD A TXT FILE WITH THE 8AD DISTANCES" << std::endl;
        exit(EXIT_FAILURE);
#endif
    }
    
    DetectorProtons.resize(NADs);
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        //        if(Analysis)
        //        {
        //            //careful with this, we need to distinguish between the different cross nH/nGd ratio events in each volume
        //            DetectorProtons[AD] = Nom->GetDetectorProtonsLs(AD)+Nom->GetDetectorProtonsGdLs(AD);//both volumes
        //        }
        //        else
        //        {
        DetectorProtons[AD] = Nom->GetDetectorProtonsGdLs(AD);
        //        }
        
        for (Int_t week = 0; week <NReactorPeriods; week++)
        {
            FullTime[AD][week] = Nom->GetFullTime(AD,week);
            NominalDetectorEfficiency[AD][week] = Nom->GetDetectorEfficiency(AD,week);
            DetectorEfficiency[AD][week]=NominalDetectorEfficiency[AD][week];
        }
    }
    
    InitialEnergy = Nom->GetEmin();
    InitialVisibleEnergy = Nom->GetEVisMin();
    FinalVisibleEnergy = Nom->GetEVisMax();
    
    n_evis_bins = Nom->GetVisibleBins();
    
    for (Int_t i = 0; i <= n_evis_bins; i++)
    {
        evis_bins[i] = Nom->GetVisibleBinningArray(i);
    }
    
    n_etrue_bins = Nom->GetTrueBins();
    
    for (Int_t i = 0; i <= n_etrue_bins; i++)
    {
        enu_bins[i] = Nom->GetTrueBinningArray(i);
    }
    
    IsotopeMatrix = Nom->GetIsotopeMatrix();
    ReactorPowerMatrix = Nom->GetReactorPowerMatrix();
    Sin22t12Matrix = Nom->GetSin22t12Matrix();
    
    VaryAccidentalMatrix = Nom->GetVaryAccidentalMatrix();
    VaryLiHeMatrix = Nom->GetVaryLiHeMatrix();
    VaryFastNeutronsMatrix = Nom->GetVaryFastNeutronsMatrix();
    VaryAmCMatrix = Nom->GetVaryAmCMatrix();
    DistortLiHeMatrix = Nom->GetDistortLiHeMatrix();
    DistortFastNeutronsMatrix = Nom->GetDistortFastNeutronsMatrix();
    DistortAmCMatrix = Nom->GetDistortAmCMatrix();
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            AccidentalError[AD][week]=Nom->GetAccidentalError(AD,week);
            LiHeError[AD][week]=Nom->GetLiHeError(AD,week);
            FastNeutronsError[AD][week]=Nom->GetFNError(AD,week);
            AmCError[AD][week]=Nom->GetAmCError(AD,week);
        }
    }
    
    delete Nom;
    
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//USEFUL CONSTRUCTOR
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Oscillation :: Oscillation(NominalData* OData)
{
    flagDelete=0;
    ResponseDirectory = OData->GetResponseDirectory();
    rand = new TRandom3(0);

    hierarchy=OData->GetHierarchy();
    deltam2_32=OData->GetDm232();
    deltam2_21=OData->GetDm221();
    dm2ee = OData->GetDm2ee();
    s22t13 = OData->GetSin22t13();
    //s22t12=0;//To check fractions and extrapolation factors
    s22t12 = OData->GetSin22t12();
    
    Analysis = OData->GetAnalysis();
    if(Analysis)
    {
        AnalysisString = "Hydrogen";
    }
    else
    {
        AnalysisString = "Gadolinium";
    }
    DataSet = OData->GetDataSet();
    
    Nweeks = OData->GetWeeks();
    NReactorPeriods = OData->GetNReactorPeriods();
    
    NADs = OData->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
#ifdef BlindedAnalysis
    sprintf(DistanceFileName,"./Distances/blinded_baseline.txt");
    std::cout << " NEED TO ADD A BLINDED DISTANCE FILE" << std::endl;
    exit(EXIT_FAILURE);
#else
    sprintf(DistanceFileName,"./Distances/unblinded_baseline.txt");
#endif
    
    if(NADs == 8) //ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
#ifdef BlindedAnalysis
        sprintf(DistanceFileName,"./Distances/blinded_baseline.txt");
        std::cout << " NEED TO ADD A BLINDED 8AD DISTANCE FILE" << std::endl;
        exit(EXIT_FAILURE);
#else
        sprintf(DistanceFileName,"./Distances/unblinded_baseline8ADs.txt");//change file to calculate distances to a file that has the information for the 8 ADs
        std::cout << "NEED TO ADD A TXT FILE WITH THE 8AD DISTANCES" << std::endl;
        exit(EXIT_FAILURE);
#endif
    }
    
    DetectorProtons.resize(NADs);
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        if(Analysis)
        {
            //careful with this, we need to distinguish between the different cross nH/nGd ratio events in each volume
            DetectorProtons[AD] = OData->GetDetectorProtonsLs(AD)+OData->GetDetectorProtonsGdLs(AD);//both volumes
        }
        else
        {
            DetectorProtons[AD] = OData->GetDetectorProtonsGdLs(AD);
        }
        
        for (Int_t week = 0; week <NReactorPeriods; week++)
        {
            FullTime[AD][week]  = OData->GetFullTime(AD,week);
            NominalDetectorEfficiency[AD][week] = OData->GetDetectorEfficiency(AD,week);
            DetectorEfficiency[AD][week] =NominalDetectorEfficiency[AD][week] ;
        }
    }
    
    InitialEnergy = OData->GetEmin();
    InitialVisibleEnergy = OData->GetEVisMin();
    FinalVisibleEnergy = OData->GetEVisMax();
    
    n_evis_bins = OData->GetVisibleBins();
    
    for (Int_t i = 0; i <= n_evis_bins; i++)
    {
        evis_bins[i] = OData->GetVisibleBinningArray(i);
    }
    
    n_etrue_bins = OData->GetTrueBins();
    
    for (Int_t i = 0; i <= n_etrue_bins; i++)
    {
        enu_bins[i] = OData->GetTrueBinningArray(i);
    }
    
    IsotopeMatrix = OData->GetIsotopeMatrix();
    ReactorPowerMatrix = OData->GetReactorPowerMatrix();
    Sin22t12Matrix = OData->GetSin22t12Matrix();
    
    VaryAccidentalMatrix = OData->GetVaryAccidentalMatrix();
    VaryLiHeMatrix = OData->GetVaryLiHeMatrix();
    VaryFastNeutronsMatrix = OData->GetVaryFastNeutronsMatrix();
    VaryAmCMatrix = OData->GetVaryAmCMatrix();
    DistortLiHeMatrix = OData->GetDistortLiHeMatrix();
    DistortFastNeutronsMatrix = OData->GetDistortFastNeutronsMatrix();
    DistortAmCMatrix = OData->GetDistortAmCMatrix();
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            AccidentalError[AD][week]=OData->GetAccidentalError(AD,week);
            LiHeError[AD][week]=OData->GetLiHeError(AD,week);
            FastNeutronsError[AD][week]=OData->GetFNError(AD,week);
            AmCError[AD][week]=OData->GetAmCError(AD,week);
        }
    }
}


Oscillation :: ~Oscillation()
{
    if(flagDelete)
    {
        for (Int_t ad = 0; ad<NADs; ad++)//far
        {
            for (Int_t reactor = 0; reactor<NReactors; reactor++)
            {
                if(Nweeks==1)
                {
                    delete InclusiveFluxH[reactor][ad];
                }
                else
                {
                    delete FluxH[reactor][ad][week];
                }
            }
        }
    }
    else
    {
        for(Int_t AD=0;AD<NADs;AD++)
        {
            delete BackgroundSpectrumH[AD];
            delete AccidentalsH[AD];
            delete LiHeH[AD];
            delete AmCH[AD];
            delete FastNeutronsH[AD];
        }

        delete NominalAmCF;
        for (Int_t ad = 0; ad<NADs; ad++)//far
        {
            for (Int_t reactor = 0; reactor<NReactors; reactor++)
            {
                delete FluxH[reactor][ad][week];
            }
        }
        
        for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
        {
            for(Int_t j = 0; j < n_evis_bins; j ++)
            {
                delete NearSpectrumH[near][j];
                delete TrueADSpectrumH[near][j];
                
                for (Int_t reactor = 0; reactor < NReactors; reactor++)
                {
                    delete NearSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)][j];
                    
                    for (Int_t far = 0; far<(ADsEH3); far++)//far
                    {
                        delete FarHallSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors][j];
                    }
                }
            }
        }
        
        for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
        {
            for (Int_t reactor = 0; reactor < NReactors; reactor++)
            {
                delete FluxFractionH[near+reactor*(ADsEH1+ADsEH2)];
                delete OscillationFunction[near*NReactors+reactor];
                for (Int_t far = 0; far<(ADsEH3); far++)//far
                {
                    delete ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors];
                }
            }
        }
        
        for (Int_t far = 0; far<(ADsEH3); far++)//far
        {
            delete FarBackgroundSpectrumH[far];
        }
        for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)//near
        {
            delete NearBackgroundSpectrumH[near];
        }
        for (Int_t far = 0; far<(ADsEH3); far++)//far
        {
            for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)//near
            {
                for(Int_t j = 0; j < n_evis_bins; j ++)
                {
                    delete FarHallSpectrumH[far*(ADsEH1+ADsEH2)+near][j];
                }
                //                FarHallSpectrumH[far*(ADsEH1+ADsEH2)+near]=NULL;
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//GETTERS AND SETTERS

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Oscillation:: SetOscillationParameters(Double_t sin22t13,Double_t Dm2_ee)
{
    s22t13=sin22t13;
    dm2ee=Dm2_ee;
}

void Oscillation:: SetEfficiency(Double_t Efficiency[])
{
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week < Nweeks; week++)
        {
            DetectorEfficiency[AD][week] = Efficiency[week+AD*Nweeks];
            
            std::cout << "\t \t \t Detector Efficiency in Oscillation: " << DetectorEfficiency[AD][week] << std::endl;
        }
    }
}

void Oscillation :: SetDelayedEfficiency(Double_t Delayed[])
{
    for(Int_t AD = 0; AD<NADs; AD++)
    {
        DetectorEfficiencyDelayed[AD] = Delayed[AD];
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//MAINSPECTRUM SHOULD CALL THIS IF NEAR HALL DATA IS USED

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Main script to produce the histograms in each AD after oscillation
void Oscillation :: OscillationFromNearHallData(Int_t Week,bool ToyMC,bool mode)
{
    flagDelete = 0;
    
    week = Week;

    Mode = mode;

    if(Mode)
    {
        RandomString="Random";
    }
    else
    {
        RandomString="Nominal";
    }
    
    std::cout <<  "\t ***********************************************************************************************" << std::endl;
    std::cout << "\t Calculating Oscillation" << std::endl;
    
    ReadDistances(DistanceFileName);
    LoadBackgrounds(week);

    if(ToyMC)//Toy MC
    {
        LoadToyHistograms(week);//Used for prediction (Covariance matrix and tests)
    }
    else
    {
        LoadNearData(week);//Used for data
    }
    
    //Apply Matrix:
    ApplyResponseMatrix();
    
    //Oscillate back and forth:
    CalculateFluxFraction();
    NearSpectrumFraction();
    GetExtrapolation();
    FarSpectrumPrediction();
    
    //    //Save oscillation functions for demonstration purposes.
    SaveOscillationFunction();
    std::cout << "\t Oscillation calculated" << std::endl;
    
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//LOAD DISTANCES BETWEEN ADS AND CORES

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Reads baseline distances from file
void Oscillation :: ReadDistances(Char_t* distanceFileName)
{
    ADdistances.resize(NADs*NReactors);
    std::ifstream infile(distanceFileName);
    std::string line;
    
    getline(infile,line); //To throw away the first line
    
    //Get baselines from text file into a Matrix
    for(Int_t i=0;i<NADs;i++)
    {
        infile >> ADdistances[i*NReactors+0] >> ADdistances[i*NReactors+1] >> ADdistances[i*NReactors+2] >> ADdistances[i*NReactors+3] >> ADdistances[i*NReactors+4] >> ADdistances[i*NReactors+5];
        //                printf("Baseline AD%d to D1 is: %f \n", i+1, ADdistances[i*NReactors+0]);
        //                printf("Baseline AD%d to D2 is: %f \n", i+1, ADdistances[i*NReactors+1]);
        //                printf("Baseline AD%d to L1 is: %f \n", i+1, ADdistances[i*NReactors+2]);
        //                printf("Baseline AD%d to L2 is: %f \n", i+1, ADdistances[i*NReactors+3]);
        //                printf("Baseline AD%d to L3 is: %f \n", i+1, ADdistances[i*NReactors+4]);
        //                printf("Baseline AD%d to L4 is: %f \n", i+1, ADdistances[i*NReactors+5]);
        //        for(int j=0;j<NReactors;j++)
        //        {
        //            ADdistances[i*NReactors+j]=ADdistances[i*NReactors+j]*km; //Km instead of meters
        //        }
    }
    
    //Test:
    if(TestAllTheSame)
    {
        for(Int_t i=0;i<NADs;i++)
        {
            for(Int_t j=0;j<NReactors;j++)
            {
                ADdistances[i*NReactors+j]=1000;
            }
        }
    }
    
    infile.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//ANTINEUTRINO SURVIVAL PROBABILITY

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Calculates oscillation probability

Double_t Oscillation :: OscProb(Double_t L, Double_t E, Double_t S22t13,Double_t Dm2ee)
{
    Double_t theta13=TMath::ASin(sqrt(S22t13))*0.5;
    Double_t theta12=TMath::ASin(sqrt(s22t12))*0.5;
    Double_t dm231 = Dm2ee-hierarchy*(5.21e-5+deltam2_21);
//    std::cout << "dm231 is: " << dm231 << "calculated from dmee: " << Dm2ee << std::endl;
  
    if(!DeltaMee)
    {
        //This calculation follows equation 9 in the TDR
        Double_t term1=s22t12*pow(cos(theta13),4)*pow(sin(1.267*deltam2_21*L/E),2);
        Double_t term2=pow(cos(theta12)*sin(2*theta13)*sin(1.267*dm231*L/E),2);
        Double_t term3=pow(sin(theta12)*sin(2*theta13)*sin(1.267*deltam2_32*L/E),2);
        return 1-(term1+term2+term3);
    }
    else
    {
        //Δm_ee aproximation
        Double_t term1=s22t12*pow(cos(theta13),4)*pow(sin(1.267*deltam2_21*L/E),2);
        Double_t term2=pow(sin(2*theta13)*sin(1.267*Dm2ee*L/E),2);
        return 1-(term1+term2);
    }
}

Double_t Oscillation :: PlotOscProb(Double_t* L, Double_t*par)
{
    Double_t E = 4;//Plot for 1MeV antineutrino
    Double_t dm231 = dm2ee-hierarchy*(5.21e-5+deltam2_21);
    Double_t S22t13 = 0.09;
    
    Double_t theta13=TMath::ASin(sqrt(S22t13))*0.5;
    Double_t theta12=TMath::ASin(sqrt(s22t12))*0.5;

    if(!DeltaMee)
    {
        //This calculation follows equation 9 in the TDR
        Double_t term1=s22t12*pow(cos(theta13),4)*pow(sin(1.267*deltam2_21*L[0]/E),2);
        Double_t term2=pow(cos(theta12)*sin(2*theta13)*sin(1.267*dm231*L[0]/E),2);
        Double_t term3=pow(sin(theta12)*sin(2*theta13)*sin(1.267*deltam2_32*L[0]/E),2);
        return 1-(term1+term2+term3);
    }
    else
    {
        //Δm_ee aproximation
        Double_t term1=s22t12*pow(cos(theta13),4)*pow(sin(1.267*deltam2_21*L[0]/E),2);
        Double_t term2=pow(sin(2*theta13)*sin(1.267*dm2ee*L[0]/E),2);
        return 1-(term1+term2);
    }
}
void Oscillation :: SaveOscillationFunction()
{
    TFile* OscillationF = new TFile("./RootOutputs/OscillationCurves.root","recreate");
    
    for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
    {
        for (Int_t k = 0; k<NReactors; k++)
        {
            OscillationFunction[near*NReactors+k]->GetXaxis()->SetTitle("E_{true} (MeV)");
            OscillationFunction[near*NReactors+k]->GetXaxis()->SetTitleSize(0.04);
            OscillationFunction[near*NReactors+k]->GetYaxis()->SetTitleSize(0.04);
            OscillationFunction[near*NReactors+k]->GetYaxis()->SetTitle("Survival Probability");
            OscillationFunction[near*NReactors+k]->Write();
        }
    }
    
    #ifdef PrintEps
        TCanvas* Curves = new TCanvas("OscCurves","OscCurves");
        
        TF1* OscCurve = new TF1("OscillationCurve",this,&Oscillation::PlotOscProb,0,4000,0,"Oscillation","PlotOscProb");
        OscCurve->GetXaxis()->SetTitle("Baseline (m)");
        OscCurve->GetXaxis()->SetTitleSize(0.04);
        OscCurve->GetYaxis()->SetTitleSize(0.04);
        OscCurve->GetYaxis()->SetTitle("Survival Probability");
        OscCurve->Draw();
        Curves->Print("./Images/OscillationCurves.eps",".eps");
        delete OscCurve;
        delete Curves;
    #endif
    delete OscillationF;
}

void Oscillation :: CalculateFluxFraction()
{
    std::cout << "\t \t \t Sin22t12 in Oscillation.h is: " << s22t12 << std::endl;
    std::cout << "\t \t \t Sin22t13 in Oscillation.h is: " << s22t13 << std::endl;
    std::cout << "\t \t \t Dm2_ee in Oscillation.h is: " << dm2ee << std::endl;
    
    //TFile* FactorF = new TFile(("./RootOutputs/"+ AnalysisString+ "/FactorStudies.root").c_str(),"recreate");
    
    Char_t filenameFraction[100];
    
    Norma.resize(n_etrue_bins*(ADsEH1+ADsEH2));
    OscillationFunction.resize((ADsEH1+ADsEH2)*NReactors);
    FluxFraction.resize(n_etrue_bins*(ADsEH1+ADsEH2)*NReactors);
    FluxFractionH.resize(NReactors*(ADsEH1+ADsEH2));
    
    for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
    {
        for (Int_t k = 0; k<NReactors; k++)
        {
            OscillationFunction[near*NReactors+k] = new TH1D(Form("OscProbFromReactor%itoNear%i",k,near),Form("OscProbFromReactor%itoNear%i",k,near),n_etrue_bins,enu_bins[0],enu_bins[n_etrue_bins]);
        }
        
        for(Int_t pts=0;pts<n_etrue_bins;pts++)
        {
            Norma[near+pts*(ADsEH1+ADsEH2)]=0.0;
            
            Energy = TrueADSpectrumH[0][0]->GetXaxis()->GetBinCenter(pts+1);
            
            for (Int_t reactor = 0; reactor<NReactors; reactor++)
            {
                Norma[near+pts*(ADsEH1+ADsEH2)] = Norma[near+pts*(ADsEH1+ADsEH2)] + ((FluxH[reactor][near][week]->GetBinContent(pts+1)*OscProb(ADdistances[near*NReactors+reactor],Energy,s22t13,dm2ee)/(ADdistances[near*NReactors+reactor]*ADdistances[near*NReactors+reactor])));//  Σ Over Cores of all fractions
                //            std::cout <<  Norma[near+pts*(ADsEH1+ADsEH2)] << "NORMA" << std::endl;
                
                OscillationFunction[near*NReactors+reactor]->SetBinContent(pts+1,OscProb(ADdistances[near*NReactors+reactor],Energy,s22t13,dm2ee));
            }
        }
        
        for (Int_t reactor = 0; reactor<NReactors; reactor++)
        {
            sprintf(filenameFraction,"FluxFraction factor, Near AD%i, Reactor%i, Week%i", near+1, reactor+1, week+1);
            
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)] = (TH1D*)FluxH[0][0][0]->Clone(filenameFraction);
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->Reset();
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->SetTitle(filenameFraction);
            
            for(Int_t pts=0;pts<n_etrue_bins;pts++)
            {
                Energy = TrueADSpectrumH[0][0]->GetXaxis()->GetBinCenter(pts+1);
                
                if(Norma[near+pts*(ADsEH1+ADsEH2)]!=0.0)
                {
                    FluxFraction[near+reactor*(ADsEH1+ADsEH2)+pts*(ADsEH1+ADsEH2)*NReactors] = (FluxH[reactor][near][week]->GetBinContent(pts+1)*OscProb(ADdistances[near*NReactors+reactor],Energy,s22t13,dm2ee))/(ADdistances[near*NReactors+reactor]*ADdistances[near*NReactors+reactor]*Norma[near+pts*(ADsEH1+ADsEH2)]);
                    
                   // std::cout << "bin: " << pts << " flux fraction: " << FluxFraction[near+reactor*(ADsEH1+ADsEH2)+pts*(ADsEH1+ADsEH2)*NReactors] << " flux: " << FluxH[reactor][near][week]->GetBinContent(pts+1) << " norm: " << Norma[near+pts*(ADsEH1+ADsEH2)] << " osc prob: " << OscProb(ADdistances[near*NReactors+reactor],Energy,s22t13,dm2ee) << std::endl;
                }
                else
                {
                    std::cout << " NORMA IS 0, THIS MUST BE AN ERROR" << std::endl;
                    
                    exit(EXIT_FAILURE);
                    
                    FluxFraction[near+reactor*(ADsEH1+ADsEH2)+pts*(ADsEH1+ADsEH2)*NReactors]=0.0;
                }
                
                FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->SetBinContent(pts+1,FluxFraction[near+reactor*(ADsEH1+ADsEH2)+pts*(ADsEH1+ADsEH2)*NReactors]);
            }
            
            //FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->Write();
        }
    }
    
    //delete FactorF;
    
    #ifdef PrintEps
    TCanvas* FluxFracC = new TCanvas("FluxFractions","FluxFractions",1000,300);

    FluxFracC->Divide(3,1);
    TLegend *leg = new TLegend(0.1,0.1,0.4,0.5);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    
    for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
    {
        FluxFracC->cd(near+1);
        
        for (Int_t reactor = 0; reactor<NReactors; reactor++)
        {
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->SetTitle("Flux Fractions");
            
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->SetLineColor(colors[reactor]);
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->GetXaxis()->SetTitle("E_{true} (MeV)");
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->GetXaxis()->SetTitleSize(0.04);
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->GetYaxis()->SetTitleSize(0.045);
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->GetXaxis()->SetLabelSize(0.045);
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->GetYaxis()->SetLabelSize(0.045);
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->SetTitle("");
            
            sprintf(name,"f_{%i,j}",near+1);
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->GetYaxis()->SetTitle(name);
            if(reactor==0)
            {
                FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->GetYaxis()->SetRangeUser(0,0.5);
                FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->Draw();
            }
            FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->Draw("same");
            if(near==0)
            {
                leg->AddEntry(FluxFractionH[near+reactor*(ADsEH1+ADsEH2)],labels[reactor].c_str(),"l");
                leg->Draw("same");
            }
            
        }
    }
    
    FluxFracC->Print(("./Images/"+AnalysisString+"/FluxFractions.eps").c_str());
    
    delete FluxFracC;
    #endif
}

void Oscillation :: NearSpectrumFraction()
{
    TFile* NearSpectrumF = new TFile(("./RootOutputs/"+ AnalysisString+ "/Spectra/NearSpectrumFraction.root").c_str(),"recreate");
    
    Char_t filenameNear[100];
    Char_t filenameFrac[100];
#ifdef CheckSumTrueSpectrum
    TH1D* SumNearSpectrumH[(ADsEH1+ADsEH2)*NReactors];

    for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
    {
        for (Int_t reactor = 0; reactor<NReactors; reactor++)
        {
            SumNearSpectrumH[near+reactor*(ADsEH1+ADsEH2)] = new TH1D(Form("Near_AD%d_Spectrum_from_reactor%d",near+1,reactor+1),Form("Near_AD%d_Spectrum_from_reactor%d",near+1,reactor+1),n_etrue_bins,enu_bins);
        }
    }
#endif
    for(Int_t j = 0; j < n_evis_bins; j++)
    {
#ifdef CheckSumTrueSpectrum

        TCanvas* NearSpectrumFC = new TCanvas("NearSpectrumFC","NearSpectrumFC");

        NearSpectrumFC->Divide(NReactors,(ADsEH1+ADsEH2));
#endif

        for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
        {
            if(Nweeks==1)
            {
                sprintf(filenameNear,"AD%i Near Prediction Vis%i", near+1,j+1);
            }
            else
            {
                sprintf(filenameNear,"AD%i Near Prediction, Week%i Vis%i", near+1,week+1,j+1);
            }
            
            NearSpectrumH[near][j] = (TH1D*) TrueADSpectrumH[near][j]->Clone(filenameNear);
            NearSpectrumH[near][j]->Reset();
            
            for (Int_t reactor = 0; reactor<NReactors; reactor++)
            {
                if(Nweeks==1)
                {
                    sprintf(filenameFrac,"AD%i Near Spectrum fraction from Reactor%i Vis%i", near+1, reactor+1,j+1);
                }
                else
                {
                    sprintf(filenameFrac,"AD%i Near Spectrum fraction from Reactor%i, Week%i Vis%i", near+1, reactor+1, week+1,j+1);
                }
                
                NearSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)][j] = (TH1D*)TrueADSpectrumH[near][j]->Clone(filenameFrac);
                NearSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)][j]->SetTitle(filenameFrac);
                
                for(Int_t pts=0;pts<TrueADSpectrumH[near][0]->GetXaxis()->GetNbins();pts++)
                {
                    Double_t SpectrumScaled = (TrueADSpectrumH[near][j]->GetBinContent(pts+1))*FluxFraction[near+reactor*(ADsEH1+ADsEH2)+pts*(ADsEH1+ADsEH2)*NReactors];
                    
 //                   std::cout << pts << " " << j << " " << TrueADSpectrumH[near][j]->GetBinContent(pts+1) << " " << FluxFraction[near+reactor*(ADsEH1+ADsEH2)+pts*(ADsEH1+ADsEH2)*NReactors] << std::endl;
                    
                    NearSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)][j]->SetBinContent(pts+1, SpectrumScaled);
                }
                
                NearSpectrumH[near][j]->Add(NearSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)][j]);//sum over cores
                
                NearSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)][j]->Write();
                
#ifdef CheckSumTrueSpectrum
                SumNearSpectrumH[near+reactor*(ADsEH1+ADsEH2)]->Add(NearSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)][j]);
                
                #ifdef PrintEps
                    NearSpectrumFC->cd(near+reactor*(ADsEH1+ADsEH2)+1);
                    
                    SumNearSpectrumH[near+reactor*(ADsEH1+ADsEH2)]->Draw();
                    
                    TrueADSpectrumH[near][j]->Draw();
                #endif
#endif
            }
            
            NearSpectrumH[near][j]->SetTitle(filenameNear);
            NearSpectrumH[near][j]->Write();
        }
#ifdef CheckSumTrueSpectrum

        NearSpectrumFC->Print(Form("./Images/SpectrumFractions/NearTrueSpectrumSumVisible%d.eps",j+1));

        delete NearSpectrumFC;
#endif
    }
    
    delete NearSpectrumF;
    
#ifdef CheckSumTrueSpectrum
    for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
    {
        for (Int_t reactor = 0; reactor<NReactors; reactor++)
        {
            delete SumNearSpectrumH[near+reactor*(ADsEH1+ADsEH2)];
        }
    }
#endif
}
void Oscillation :: GetExtrapolation()
{
  //  TFile* FactorF = new TFile(("./RootOutputs/"+ AnalysisString+ "/FactorStudies.root").c_str(),"update");
    
    Char_t filenameExtrapolation[100];
    
    Extrapolation.resize(n_etrue_bins*ADsEH3*NReactors*(ADsEH1+ADsEH2));
    ExtrapolationH.resize(ADsEH3*NReactors*(ADsEH1+ADsEH2));
#ifdef PrintEps
    TCanvas *c2 = new TCanvas("Extrapolation","Extrapolation",900,900);
    c2->Divide(3,3);
    Int_t cont=0;
#endif
    
    for (Int_t far = 0; far<(ADsEH3); far++)//far
    {
        for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)//near
        {
            #ifdef PrintEps
            ++cont;
            c2->cd(cont);
            #endif
            for (Int_t reactor = 0; reactor<NReactors; reactor++)
            {
                if(Nweeks==1)
                {
                    sprintf(filenameExtrapolation,"Extrapolation factor, Near AD%i, Far AD%i, Reactor%i", near+1, far+1, reactor+1);
                }
                else
                {
                    sprintf(filenameExtrapolation,"Extrapolation factor, Near AD%i, Far AD%i, Reactor%i, Week%i", near+1, far+1, reactor+1, week+1);
                }
                
                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors] = (TH1D*)TrueADSpectrumH[0][0]->Clone(filenameExtrapolation);
                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->SetTitle(filenameExtrapolation);
                
                for(Int_t pts=0;pts<n_etrue_bins;pts++)
                {
                    Energy = TrueADSpectrumH[near][0]->GetXaxis()->GetBinCenter(pts+1);
                    
                    Extrapolation[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors+pts*(ADsEH1+ADsEH2)*NReactors*ADsEH3]=
                    (FluxH[reactor][far+(ADsEH1+ADsEH2)][week]->GetBinContent(pts+1)*
                     //(DetectorEfficiency[(far+(ADsEH1+ADsEH2))*Nweeks+week]*FullTime[(far+(ADsEH1+ADsEH2))][week]*DetectorProtons[far+(ADsEH1+ADsEH2)])*
                     OscProb(ADdistances[(far+(ADsEH1+ADsEH2))*NReactors+reactor],Energy,s22t13,dm2ee)*ADdistances[near*NReactors+reactor]*ADdistances[near*NReactors+reactor])/
                    (FluxH[reactor][near][week]->GetBinContent(pts+1)
                    // *(DetectorEfficiency[near*Nweeks+week]*FullTime[near][week]*DetectorProtons[near])
                     *OscProb(ADdistances[near*NReactors+reactor],Energy,s22t13,dm2ee)*ADdistances[(far+(ADsEH1+ADsEH2))*NReactors+reactor]*ADdistances[(far+(ADsEH1+ADsEH2))*NReactors+reactor]);
                    
                    ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->SetBinContent(pts+1,Extrapolation[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors+pts*(ADsEH1+ADsEH2)*NReactors*ADsEH3]);
                }
                #ifdef PrintEps
                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->SetLineColor(colors[reactor]);
                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->GetXaxis()->SetTitle("E_{true} (MeV)");
                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->GetXaxis()->SetTitleSize(0.040);
                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->GetYaxis()->SetTitleSize(0.045);

                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->GetXaxis()->SetLabelSize(0.045);
                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->GetYaxis()->SetLabelSize(0.045);
                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->SetTitle("");
                
                sprintf(name,"e_{%ij,%i}",near+1,far+1);
                
                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->GetYaxis()->SetTitle(name);
                
                if(reactor==0)
                {
                    ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->GetYaxis()->SetRangeUser(0,0.885);
                    ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->Draw();
                }
                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->Draw("same");
                #endif
                //ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->Write();
                
            }
        }
    }
    #ifdef PrintEps
    c2->Print(("./Images/"+AnalysisString+"/ExtrapolationFactors.eps").c_str());
    delete c2;
    #endif
    
    //delete FactorF;
}

void Oscillation :: FarSpectrumPrediction()
{
    Char_t filenameFraction[100];
    Char_t filenameFar[100];
    
    TFile* FarSpectrumFractionF = new TFile(("./RootOutputs/"+ AnalysisString+ "/Spectra/FarSpectrumFraction.root").c_str(),"recreate");
    
    TCanvas* FarSpectrumFractionC;
    
    for(Int_t j = 0; j < n_evis_bins; j ++)
    {
        #ifdef PrintEps
            FarSpectrumFractionC = new TCanvas(Form("FarSpectrumFractionVisible%d",j+1),Form("FarSpectrumFractionVisible%d",j+1));
            FarSpectrumFractionC->Divide(ADsEH1+ADsEH2,ADsEH3);
        #endif
        
        for (Int_t far = 0; far<(ADsEH3); far++)//far
        {
            for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)//near
            {
                if(Nweeks==1)
                {
                    sprintf(filenameFar,"AD%i Far Spectrum prediction from near AD%i Vis %d", far+1, near+1,j+1);
                }
                else
                {
                    sprintf(filenameFar,"AD%i Far Spectrum prediction from near AD%i, Week%i Vis%d", far+1, near+1,week+1,j+1);
                }
                
                FarHallSpectrumH[far*(ADsEH1+ADsEH2)+near][j] = (TH1D*) TrueADSpectrumH[near][0]->Clone(filenameFar);
                FarHallSpectrumH[far*(ADsEH1+ADsEH2)+near][j]->Reset();
                
                for (Int_t reactor = 0; reactor<NReactors; reactor++)
                {
                    
                    sprintf(filenameFraction,"AD%i Far Spectrum fraction from Reactor%i and near AD%i, Week%i Vis%i", far+(ADsEH1+ADsEH2)+1, reactor+1, near+1, week+1,j+1);
                    
                    FarHallSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors][j] = (TH1D*)TrueADSpectrumH[near][j]->Clone(filenameFraction);
                    FarHallSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors][j]->Reset();
                    
                    for(Int_t pts=0;pts<n_etrue_bins;pts++)
                    {
                        Double_t FarValue = NearSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)][j]->GetBinContent(pts+1)*Extrapolation[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors+pts*(ADsEH1+ADsEH2)*NReactors*ADsEH3];
                        
                        FarHallSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors][j]->SetBinContent(pts+1, FarValue);
                        
                    }
                    
                    FarHallSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors][j]->SetTitle(filenameFraction);
                    
                    // Save predictions:
                    
                    FarHallSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors][j]->Write();
                    
                    FarHallSpectrumH[far*(ADsEH1+ADsEH2)+near+week][j]->Add(FarHallSpectrumFractionH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors][j],1);//sum over cores
                }
                
                // Multiply by efficiencies and the corresponding # of protons in each AD
                
                FarHallSpectrumH[far*(ADsEH1+ADsEH2)+near][j]->Scale(DetectorEfficiency[(far+(ADsEH1+ADsEH2))][week]*FullTime[(far+(ADsEH1+ADsEH2))][week]*DetectorProtons[far+(ADsEH1+ADsEH2)]);
                
                FarHallSpectrumH[far*(ADsEH1+ADsEH2)+near][j]->SetTitle(filenameFar);
                FarHallSpectrumH[far*(ADsEH1+ADsEH2)+near][j]->Write();
                #ifdef PrintEps
                FarSpectrumFractionC->cd(near+far*(ADsEH1+ADsEH2)+1);
                FarHallSpectrumH[far*(ADsEH1+ADsEH2)+near][j]->Draw();
                #endif
            }
        }
        #ifdef PrintEps
          //  FarSpectrumFractionC->Print(Form("./Images/SpectrumFractions/FarSpectrumFractionVisible%d.eps",j+1));
            delete FarSpectrumFractionC;
        #endif

    }
    delete FarSpectrumFractionF;
}

TH1D* Oscillation :: GetOscillatedADSpectrum(Int_t near,Int_t far,Int_t j,Int_t week)
{
    return FarHallSpectrumH[far*(ADsEH1+ADsEH2)+near][j];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//                                                            LOAD REACTOR TOY DATA

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Loads near hall data
void Oscillation :: LoadToyHistograms(Int_t week)
{
    std::cout << "\t \t Loading toy reactor model" << std::endl;
    Char_t ADdata[200];
    
    if((Sin22t12Matrix||IsotopeMatrix||ReactorPowerMatrix)&&Mode==1)
    {
        sprintf(ADdata,("./RootOutputs/"+ AnalysisString+ "/RandomOutputs/RandomOscillation_Isotope_%d_Power_%d_Sin22t12_%d_.root").c_str(),IsotopeMatrix,ReactorPowerMatrix,Sin22t12Matrix);
    }
    else
    {
        sprintf(ADdata,("./RootOutputs/"+ AnalysisString+ "/NominalOutputs/Oscillation.root").c_str());
    }
    
    TFile* NearHallDataF = new TFile(ADdata);
    
    for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
    {
        NearHallDataF->cd("Total AD Spectra after oscillation");
        
        ADSpectrumVisH[near] = (TH1D*)gDirectory->Get(Form("Oscillation Prediction AD%d, week%d", near+1,week));
        //  In oscillationreactor ADs have been corrected.

        ADSpectrumVisH[near]->Add(BackgroundSpectrumH[near],-1);
    }
    
    delete NearHallDataF;
    
    #ifdef PrintEps
        TCanvas* cr = new TCanvas("cr","cr", 1200,400);
        cr->Divide(ADsEH1+ADsEH2,1);
        
        for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            cr->cd(near+1);
            ADSpectrumVisH[near]->Draw("HIST");
        }
        
        cr->Print(("./Images/"+ AnalysisString+ "/"+RandomString+"NearReactorPredictionWithoutBackground.eps").c_str(),".eps");
        delete cr;
    #endif
    
    TFile* checkADS = new TFile("./RootOutputs/VisibleBinCheck.root","recreate");
    for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
    {
        ADSpectrumVisH[near]->Write();
    }
    delete checkADS;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//LOAD NEAR HALL DATA

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Oscillation :: LoadBackgrounds(Int_t week)
{
    LoadNominalBackgrounds();
    
    if(Mode)
    {
        if(DistortLiHeMatrix)
        {
            DistortLiHe=0.2;
        }
        if(DistortFastNeutronsMatrix)
        {
            DistortFN=0.2;
        }
        if(DistortAmCMatrix)
        {
            DistortAmC=0.15;
        }
 
        FluctuateBackgrounds(week);//Vary backgrounds
    }
    
    for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
    {
        NearBackgroundSpectrumH[near]=(TH1D*)BackgroundSpectrumH[near*Nweeks+week]->Clone();
    }
    for (Int_t far = 0; far < ADsEH3; far++)
    {
        FarBackgroundSpectrumH[far]=(TH1D*)BackgroundSpectrumH[(far+ADsEH1+ADsEH2)]->Clone();
    }
    
    SaveBackgrounds();
}

void Oscillation :: LoadNearData(Int_t week)
{
    std::cout << " LOADING NEAR DATA " << std::endl;
    Char_t NearData[100];
    Char_t NearDataSpec[100];
    
    if(Analysis)
    {
        std::cout << "\t \t \t \t \t \t IBD DATA FILE NEEDED FOR HYDROGEN ANALYSIS" << std::endl;
        sprintf(NearData,"./Inputs/HInputs/ibd_eprompt_shapes.root");
    }
    else
    {
        if(DataSet==2)//LBNL
        {
            //            sprintf(NearData,"./Inputs/GdInputs/ibd_eprompt_shapes.root");
            sprintf(NearData,"./Inputs/GdInputs/IHEP_data_lbnlbin_6AD.root");
            
        }
        else if(DataSet==1)//P12E
        {
            sprintf(NearData,"./Inputs/GdInputs/file_P12E_nGd_IBD_Acc_spectrum.root");
        }
        else if(DataSet==0)
        {
            std::cout << "ERROR, USING NEAR DATA WITH SIMPLE REACTOR MODEL" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    
    for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
    {
        if(Nweeks == 1)
        {
            if(DataSet==2)//LBNL
            {
                //                sprintf(NearDataSpec,"h_ibd_eprompt_inclusive_ad%i",near+1);
                if(near<ADsEH1)
                {
                    sprintf(NearDataSpec,"h_ibd_eprompt_inclusive_eh1_ad%i",near+1);
                }
                else
                {
                    sprintf(NearDataSpec,"h_ibd_eprompt_inclusive_eh2_ad%i",near-ADsEH1+1);
                }
                
            }
            else if(DataSet==1)//P12E
            {
                if(near<ADsEH1)
                {
                    sprintf(NearDataSpec,"hist_AccSub_%d_EH1",near+1);
                }
                if(near>ADsEH1)
                {
                    sprintf(NearDataSpec,"hist_AccSub_%d_EH2",near-ADsEH1+1);
                }
            }
        }
        else
        {
            if(DataSet==2)//LBNL
            {
                sprintf(NearDataSpec,"h_ibd_eprompt_week%d_ad%d", week, near+1);
            }
        }
        
        TFile* NearDataF = new TFile(NearData);
        
        ADSpectrumVisH[near] = (TH1D*)gDirectory->Get(NearDataSpec);
        
        delete NearDataF;
        
        if(DataSet!=2)//need to rebin files to LBNL binning
        {
            ADSpectrumVisH[near]=(TH1D*)ADSpectrumVisH[near]->Rebin(n_evis_bins,Form("Rebinned Vis Near%d Data Spectrum",near),evis_bins);
        }
        
        //Substract backgrounds:
        ADSpectrumVisH[near]->Add(NearBackgroundSpectrumH[near],-1);//Substract (varied if necessary) backgrounds from data
        
        for(Int_t i = 0; i< ADSpectrumVisH[near]->GetXaxis()->GetNbins();i++)
        {
            if(ADSpectrumVisH[near]->GetBinContent(i+1)<0)
            {
                ADSpectrumVisH[near]->SetBinContent(i+1,0);//Check if using 0 would will produce NaNs! Maybe need to extrapolate
            }
        }
    }
    
    #ifdef PrintEps
        TCanvas* c1 = new TCanvas("c1","Near Data", 1200,400);
        c1->Divide(ADsEH1+ADsEH2,1);
        
        for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            c1->cd(near+1);
            ADSpectrumVisH[near]->Draw();
        }
        
        c1->Print(("./Images/"+ AnalysisString+ "/"+RandomString+"NearData.eps").c_str(),".eps");
        
        if(RandomString!="Nominal")
        {
            std::cout << "REAL DATA SHOULD HAVE NOMINAL BACKGROUNDS SUBTRACTED INSTEAD OF RANDOM BACKGROUNDS, CHECK CASES" << std::endl;
            exit(EXIT_FAILURE);
        }
        delete c1;
    #endif
}

void Oscillation :: ApplyResponseMatrix()
{
    if(!strcmp((ResponseDirectory).c_str(),""))
    {
        TFile* TransEnergyMatrixDataF = new TFile(("./ResponseMatrices/"+ AnalysisString+ "/NominalResponseMatrix.root").c_str());
            TransEnergyMatrix = (TH2D*)gDirectory->Get("EnuEvis");
        
        delete TransEnergyMatrixDataF;
    }
    else
    {
        TFile* TransEnergyMatrixDataF = new TFile((ResponseDirectory).c_str());
            TransEnergyMatrix = (TH2D*)gDirectory->Get("EnuEvis");
        delete TransEnergyMatrixDataF;
    }
    
    //If other binnings are used we should include a Rebin method here to match the binning, or recalculate the matrix before running this.
    
    //        TransEnergyMatrix->Rebin();
    
    for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
    {
        ADSpectrumVisH[near]->Scale(1./(DetectorEfficiency[near][week]*FullTime[near][week]*DetectorProtons[near]));// Remove AD performance included in oscillation reactor and data. This is later included in the far spectrum fraction calculation.
        
        for(Int_t j = 0; j < n_evis_bins; j ++)
        {
            TrueADSpectrumH[near][j] = new TH1D(Form("TrueEnergyData%i",j),Form("TrueEnergyData%i",j),n_etrue_bins,enu_bins);
        }
        
        for(Int_t i=0; i<n_etrue_bins; i++)
        {
            for(Int_t j=0; j<n_evis_bins; j++)
            {
                //[True_0 ... True_n] = [0,0]  [...] [0,m]  * [Vis_0]
                //                      [0,n]  [...] [n,m]     [...]
                //                                            [Vis_m]
                
                // [nx1] = [n x m] x [mx1]
                TrueADSpectrumH[near][j]->SetBinContent(i+1,TrueADSpectrumH[near][j]->GetBinContent(i+1)+TransEnergyMatrix->GetBinContent(j+1,i+1)*ADSpectrumVisH[near]->GetBinContent(j+1));
            }
        }
    }
    if(TrueADSpectrumH[0][0]->GetXaxis()->GetNbins()!=TransEnergyMatrix->GetYaxis()->GetNbins())
    {
        std::cout << "Binning disagreement between true spectrum and response matrix" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    
    if(ADSpectrumVisH[0]->GetXaxis()->GetNbins()!=TransEnergyMatrix->GetXaxis()->GetNbins())
    {
        std::cout << "Binning disagreement between vis spectrum and response matrix" << std::endl;
        exit(EXIT_FAILURE);
    }
    delete TransEnergyMatrix;
}

void Oscillation :: GenerateFluxHisto()
{
    std::cout << "Generating efficiency weighted flux Histogram with " << NReactorPeriods << " periods of data" << std::endl;
    
    flagDelete = 1;
    
    Char_t ReactorData[100];
    Char_t ReactorSpectrum[100];
    
    //Reactor model is not Hydrogen/Gadolinium dependent
    
    if (NReactorPeriods == 1)//testing superhistogram vs inclusive histogram
    {
        //need to use superhistograms? efficiency corrected fluxes in each AD
        if(TestAllTheSame)
        {
            sprintf(ReactorData,("./RootOutputs/"+ AnalysisString+"/NominalOutputs/AntineutrinoSpectrum.root").c_str()); //The reactor model spectrum, with no reactor on/off information, just assumes the nominal spectrum is on all the time.
        }
        else
        {
#ifdef BlindedAnalysis
            sprintf(ReactorData,Form("./ReactorInputs/WeeklyFlux_%dweek_blinded_inclusive.root",NReactorPeriods));
            std::cout << " NEEDED A BLINDED WEEKLY INCLUSIVE FLUX FILE" << std::endl;
            exit(EXIT_FAILURE);
#else
            sprintf(ReactorData,Form("./ReactorInputs/WeeklyFlux_%dweek_unblinded_inclusive.root",NReactorPeriods));
#endif
        }
    }
    else
    {
        
#ifdef BlindedAnalysis
        sprintf(ReactorData, Form("./ReactorInputs/WeeklyFlux_%dweek_blinded.root",NReactorPeriods));
        std::cout << " NEEDED A BLINDED WEEKLY FLUX FILE" << std::endl;
        exit(EXIT_FAILURE);
#else
        sprintf(ReactorData, Form("./ReactorInputs/WeeklyFlux_%dweek_unblinded.root",NReactorPeriods));
#endif
    }
    
#ifdef NoOscillation
    //To show a flat extrapolation factor and flux fraction
    sprintf(ReactorData,("./RootOutputs/"+ AnalysisString+"/NominalOutputs/AntineutrinoSpectrum.root").c_str()); //The reactor model spectrum, with no reactor on/off information, just assumes the nominal spectrum is on all the time.
#endif
    
    TFile* ReactorDataF = new TFile(ReactorData);
    
    for(Int_t period = 0; period < NReactorPeriods; period++)
    {
        for(Int_t reactor = 0; reactor < NReactors; reactor++)
        {
            if(Nweeks==1)
            {
                if(TestAllTheSame)
                {
                    sprintf(ReactorSpectrum,"AntineutrinoSpectrumFromReactor%i",reactor+1);
                    // In this case I don't scale it because the spectrum has not been altered usng the # of protons, technically is spectrum in the reactors not in the ad.
                }
                else
                {
                    sprintf(ReactorSpectrum,"Week%i/%i",period,reactor);
                }
            }
            else
            {
                sprintf(ReactorSpectrum,"Week%i/%i", period, reactor);
            }
            #ifdef NoOscillation
            //To show a flat extrapolation factor and flux fraction
            sprintf(ReactorSpectrum,"AntineutrinoSpectrumFromReactor%i",reactor+1);
            #endif
            TH1D* NonBinnedFluxH = (TH1D*)gDirectory->Get(ReactorSpectrum);
            
            //        NonBinnedFluxH = (TH1D*)NonBinnedFluxH->Rebin(n_etrue_bins,Form("Flux True Spectrum from reactor%d",reactor+1),enu_bins);
            
            for(Int_t AD=0;AD<NADs;AD++)
            {
                FluxH[reactor][AD][period] = new TH1D(Form("FluxH%i_%i_%i",reactor+1,AD+1,period),Form("FluxH%i_%i_%i",reactor+1,AD+1,period),n_etrue_bins,enu_bins);
                
                Int_t lastbin = n_etrue_bins;
                
                for(Int_t i = 0; i<n_etrue_bins; i++)
                {
                    if(0<NonBinnedFluxH->GetBinContent(i+1))
                    {
                        FluxH[reactor][AD][period]->SetBinContent(i+1,NonBinnedFluxH->GetBinContent(i+1));
                        lastbin=i;
                    }
                    else//the flux ends at 8.8 MeV, otherwise the bins are 0 and we get NaNs. This extrapolation produces larger error in the large enery band. Need to do it consistently with other groups.
                    {
                        FluxH[reactor][AD][period]->SetBinContent(i+1,NonBinnedFluxH->GetBinContent(lastbin+1));
                    }
                }
            }
            
            delete NonBinnedFluxH;
        }
    }
    
    //Correct flux for efficiencies:
    
    for(Int_t AD=0;AD<NADs;AD++)
    {
        for(Int_t reactor=0;reactor<NReactors;reactor++)
        {
            InclusiveFluxH[reactor][AD] = new TH1D(Form("InclusiveFlux%i_%i",reactor+1,AD+1),Form("InclusiveFlux%i_%i",reactor+1,AD+1),n_etrue_bins,enu_bins);
            
            Double_t scalesum=0;
            
            for(Int_t week=0;week<NReactorPeriods;week++)
            {
                Double_t scale = FullTime[AD][week]*NominalDetectorEfficiency[AD][week];//DetectorEfficiency[AD][week] to use random efficiencies. LBNL pregenerates the flux so the random changes are not contemplated.
                
                std::cout << " AD: " << AD << "- Week: " << week << " - Scale: " << scale << " - Full time: " << FullTime[AD][week] << " - Efficiency:" << NominalDetectorEfficiency[AD][week] << std::endl;
                
                //Something to consider:
                
                // If we want to apply random detector efficiencies to the flux in different periods, we should considerate their correlation. If they were correlated we should produce the efficiencies in Prediction.h in order to hold the variations in memory and share a common efficiency. So far I haven't attempted to do a multi period fit, but if you do you should hold this carefully for every systematic.
                
                scalesum+=scale;
                
                FluxH[reactor][AD][week]->Scale(scale);
                
                InclusiveFluxH[reactor][AD]->Add(FluxH[reactor][AD][week]);
                
                FluxH[reactor][AD][week]->Scale(1./scale);//if fluxh are used weekly we don't need to correct them for efficiencies

                //Rebin Flux Spectrum to match the data binning
                FluxH[reactor][AD][week]=(TH1D*)FluxH[reactor][AD][week]->Rebin(n_etrue_bins,Form("Rebinned Flux Spectrum from reactor%d with week%d ad%d efficiencies",reactor+1,week+1,AD+1),enu_bins);
                
            }
            
            InclusiveFluxH[reactor][AD]->Scale(1./scalesum);
            //Rebin Flux Spectrum to match the data binning
            InclusiveFluxH[reactor][AD]=(TH1D*)InclusiveFluxH[reactor][AD]->Rebin(n_etrue_bins,Form("Rebinned Flux Spectrum from reactor%d with ad%d efficiencies",reactor+1,AD+1),enu_bins);
        }
    }
    
    #ifdef PrintEps
        for(Int_t week=0;week<NReactorPeriods;week++)
        {
            TCanvas* WeeklyFluxC = new TCanvas("WeeklyFluxC","WeeklyFluxC");
            
            WeeklyFluxC->Divide(NReactors,NADs);
            
            for(Int_t AD=0;AD<NADs;AD++)
            {
                for(Int_t reactor=0;reactor<NReactors;reactor++)
                {
                    WeeklyFluxC->cd(reactor+AD*NReactors+1);
                    FluxH[reactor][AD][week]->Draw();
                }
            }
            
            WeeklyFluxC->Print(Form("./Images/Reactor/Weekly%dFluxes.eps",week));
            
            delete WeeklyFluxC;
        }
        
        TCanvas* FluxC = new TCanvas("FluxC","FluxC");
        
        FluxC->Divide(NReactors,NADs);
        
        for(Int_t AD=0;AD<NADs;AD++)
        {
            for(Int_t reactor=0;reactor<NReactors;reactor++)
            {
                FluxC->cd(reactor+AD*NReactors+1);
                InclusiveFluxH[reactor][AD]->Draw();
            }
        }
        
        FluxC->Print("./Images/Reactor/InclusiveAndWeeklyFluxes.eps");
        
        delete FluxC;
    #endif
    
    delete ReactorDataF;
    
    TFile* FluxF = new TFile(Form("./Inputs/SuperFlux%d.root",NReactorPeriods),"recreate");
    
    for(Int_t AD=0;AD<NADs;AD++)
    {
        for(Int_t reactor=0;reactor<NReactors;reactor++)
        {
            InclusiveFluxH[reactor][AD]->Write();
        }
    }
    
    for(Int_t period = 0; period < NReactorPeriods; period++)
    {
        for(Int_t AD=0;AD<NADs;AD++)
        {
            for(Int_t reactor=0;reactor<NReactors;reactor++)
            {
                FluxH[reactor][AD][period]->Write();
            }
        }
    }
    
    delete FluxF;
}

void Oscillation :: LoadFluxHisto()//weekly data periods
{
    flagDelete = 1;
    
    TFile* FluxF = new TFile(Form("./Inputs/SuperFlux%d.root",NReactorPeriods));
    
    if(Nweeks==1)//Inclusive fit
    {
        for(Int_t AD=0;AD<NADs;AD++)
        {
            for(Int_t reactor=0;reactor<NReactors;reactor++)
            {
                InclusiveFluxH[reactor][AD]=(TH1D*)FluxF->Get(Form("Rebinned Flux Spectrum from reactor%d with ad%d efficiencies",reactor+1,AD+1));
            }
        }
    }
    else//weekly fit
    {
        for(Int_t period = 0; period < Nweeks; period++)
        {
            for(Int_t AD=0;AD<NADs;AD++)
            {
                for(Int_t reactor=0;reactor<NReactors;reactor++)
                {
                    FluxH[reactor][AD][period]=(TH1D*)FluxF->Get(Form("Rebinned Flux Spectrum from reactor%d with week%d ad%d efficiencies",reactor+1,period+1,AD+1));
                }
            }
        }
    }
    
    delete FluxF;
}

TH1D* Oscillation :: GetFluxHisto(Int_t reactor, Int_t AD, Int_t week)
{
    if(Nweeks == 1)
    {
        return InclusiveFluxH[reactor][AD];
    }
    else
    {
        return FluxH[reactor][AD][week];
    }
}

void Oscillation :: SetFluxHistograms(TH1D* fluxH,Int_t reactor,Int_t AD,Int_t week)
{
    FluxH[reactor][AD][week] = (TH1D*)fluxH->Clone();
}

void Oscillation :: SetSin22t12(Double_t S22t12)
{
    s22t12=S22t12;
}

void Oscillation :: SaveBackgrounds(){
    
    std::cout << "\t Saving Backgrounds" << std::endl;

    TFile* BackgroundF1 = new TFile(("./RootOutputs/"+ AnalysisString+ "/Backgrounds/"+RandomString+"Backgrounds.root").c_str(),"recreate");
    
    for(Int_t near=0; near<(ADsEH1+ADsEH2); near++)
    {
        NearBackgroundSpectrumH[near]->Write((Form("Near AD%i ",near)+RandomString+ Form(" Background Period%d",week)).c_str());
    }
    for(Int_t far=0; far<ADsEH3; far++)
    {
        FarBackgroundSpectrumH[far]->Write((Form("Far AD%i ",far)+RandomString+ Form(" Background Period%d",week)).c_str());
    }
    delete BackgroundF1;
}



//Randomize the events (rate) inside errors taking into account the corresponding correlations between backgrounds and ADs. Also random shape variations are included using distortion functions.
void Oscillation :: FluctuateBackgrounds(Int_t week)
{
    for (Int_t AD =0; AD<NADs; AD++)
    {
        RandomAccidentalsH[AD]=(TH1D*)AccidentalsH[AD]->Clone();
        RandomLiHeH[AD]=(TH1D*)LiHeH[AD]->Clone();
        RandomFastNeutronsH[AD]=(TH1D*)FastNeutronsH[AD]->Clone();
        RandomAmCH[AD]=(TH1D*)AmCH[AD]->Clone();
        
        //  Initialize scale factors
        if(VaryAccidentalMatrix||VaryLiHeMatrix||VaryFastNeutronsMatrix||VaryAmCMatrix)
        {
            ScaleFactorAccidental[AD]=1.;
            ScaleFactorLiHe[AD]=1.;
            ScaleFactorFastNeutrons[AD]=1.;
            ScaleFactorAmC[AD]=1.;
        }
        
        if(VaryAccidentalMatrix)
        {
            rand->SetSeed(0);
            ScaleFactorAccidental[AD]=(1.+(AccidentalError[AD][week]*rand->Gaus(0,1)));
            //            std::cout << "\t Accidental error " << AccidentalError[AD][week] << std::endl;
        }
        if(VaryLiHeMatrix)
        {
            rand->SetSeed(0);
            ScaleFactorLiHe[AD]=(1.+(LiHeError[AD][week]*rand->Gaus(0,1)));
            //            std::cout <<"\t LiHe error " << LiHeError[AD][week] << std::endl;
        }
        if(VaryFastNeutronsMatrix)
        {
            rand->SetSeed(0);
            ScaleFactorFastNeutrons[AD]=(1.+(FastNeutronsError[AD][week]*rand->Gaus(0,1)));
            //            std::cout <<"\t Fast Neutron error " << FastNeutronsError[AD][week] << std::endl;
        }
        if(VaryAmCMatrix)
        {
            rand->SetSeed(0);
            ScaleFactorAmC[AD]=(1.+(AmCError[AD][week]*rand->Gaus(0,1)));
            //            std::cout <<"\t AmC error " << AmCError[AD][week] << std::endl;
        }
    }
    
    //Accidentals uncorrelated
    if(VaryLiHeMatrix)
    {   //  LiHe correlated by hall
        ScaleFactorLiHe[1]=ScaleFactorLiHe[0];// EH1
        ScaleFactorLiHe[4]=ScaleFactorLiHe[3];// EH3
        ScaleFactorLiHe[5]=ScaleFactorLiHe[3];// EH3
    }
    if(VaryFastNeutronsMatrix)
    {
        //  Fast neutrons by hall
        ScaleFactorFastNeutrons[1]=ScaleFactorFastNeutrons[0];// EH1
        ScaleFactorFastNeutrons[4]=ScaleFactorFastNeutrons[3];// EH3
        ScaleFactorFastNeutrons[5]=ScaleFactorFastNeutrons[3];// EH3
    }
    if(VaryAmCMatrix)
    {   //  AmC all correlated
        ScaleFactorAmC[1]=ScaleFactorAmC[0];
        ScaleFactorAmC[2]=ScaleFactorAmC[0];
        ScaleFactorAmC[3]=ScaleFactorAmC[0];
        ScaleFactorAmC[4]=ScaleFactorAmC[0];
        ScaleFactorAmC[5]=ScaleFactorAmC[0];
    }
    
    //Since the rate distortion is not accounted here, in case we want to use all the distortions at the same time we have to distort the shape first, and then the rate:
    
    if(DistortLiHeMatrix)// Distort LiHe shape for all ADs in the same way
    {
        if(Analysis)//Hydrogen
        {
            //9Li8He
            TFile* LiF = TFile::Open("./BackgroundSpectrum/HBackground/Li9_beta_neutron.root");
            TH1D* LithiumH = (TH1D*)LiF->Get("h");
            delete LiF;//90%
            
            TFile* HF = TFile::Open("./BackgroundSpectrum/HBackground/He8_beta_neutron.root");
            TH1D* HeliumH = (TH1D*)HF->Get("h");
            delete HF;//10%
            
            Double_t FractionError = 0.2;//This number is provisional
            Double_t DistortLiHeError = FractionError*rand->Gaus(0,1);
           
            if((.9+DistortLiHeError)>1)
            {
                DistortLiHeError = 0.1;//maximum error to avoid positive fractions
            }
            else if((.9+DistortLiHeError)<0)
            {
                DistortLiHeError = 0.9;//maximum error to avoid negative fractions
            }
            
            std::cout << LiHeError << std::endl;
            
            RandomLiHeH[NADs] = new TH1D("LiHe","LiHe",n_evis_bins,evis_bins);
            
            for(Int_t AD=0;AD<NADs;AD++)
            {
                for (Int_t visbin = 1; visbin <= n_evis_bins; visbin++)
                {
                    RandomLiHeH[AD]->SetBinContent(visbin,(.9+DistortLiHeError)*LithiumH->Interpolate(RandomLiHeH[AD]->GetXaxis()->GetBinCenter(visbin))+(.1-DistortLiHeError)*HeliumH->Interpolate(RandomLiHeH[AD]->GetXaxis()->GetBinCenter(visbin)));
                }
            }
            for(Int_t AD=0;AD<NADs;AD++)
            {
                RandomLiHeH[AD]->Scale(LiHeH[AD]->Integral()/RandomLiHeH[AD]->Integral());
            }
            
            TCanvas* LiC = new TCanvas("a","");

                LithiumH->Draw();
            
            LiC->Print("./Images/Li.eps");
            delete LiC;
            
            TCanvas* HeC = new TCanvas("b","");

                HeliumH->Draw();
            
            HeC->Print("./Images/He.eps");
            delete HeC;
            
            TCanvas* LiHeC = new TCanvas("c","");
            LiHeC->Divide(NADs/2,2);
            for(Int_t AD=0;AD<NADs;AD++)
            {
                LiHeC->cd(AD+1);
                RandomLiHeH[AD]->Draw();
            }
            LiHeC->Print("./Images/LiHe.eps");
            delete LiHeC;

            delete LithiumH;
            delete HeliumH;
        }
        else
        {
            //function to obtaing distortion function or histogram:
            GetDistortionFunction(DistortLiHe,func_LiHe);
            
#ifdef UseLiHeToyMC//If TOY MC generated histogram:
            TFile* OriginalLiHeF = TFile::Open("./BackgroundSpectrum/GDBackground/li9_spectrum.root");
            
            TH1F* OriginalLiHeH=(TH1F*)gDirectory->Get("h_li9_smeared_toy");
            
            OriginalLiHeF->Close();
            
            TH1F* func_LiHe_interpolated = (TH1F*)OriginalLiHeH->Clone();//To solve different bin limits error
            
            func_LiHe_interpolated->Reset();
            
            for(Int_t i = 1; i<=func_LiHe_interpolated->GetXaxis()->GetNbins();i++)
            {
                func_LiHe_interpolated->SetBinContent(i,func_LiHe->Interpolate(func_LiHe_interpolated->GetXaxis()->GetBinCenter(i)));
            }
            
            OriginalLiHeH->Multiply(func_LiHe_interpolated);

            for(Int_t AD=0;AD<NADs;AD++)
            {
                RandomLiHeH[AD] = (TH1D*)OriginalLiHeH->Rebin(n_evis_bins,Form("Rebinned LiHe"),evis_bins);//Rebin to visible binning
                
                RandomLiHeH[AD]->Scale(LiHeH[AD]->Integral()/RandomLiHeH[AD]->Integral());
            }
            
            delete func_LiHe_interpolated;
            delete OriginalLiHeH;
#else//If distortion function
            
            func_LiHe = new TF1("func_LiHe","TMath::Abs([0]+[1]*x)",InitialVisibleEnergy,FinalVisibleEnergy);

            for(Int_t AD=0;AD<NADs;AD++)
            {
                RandomLiHeH[AD]->Multiply(func_LiHe);
                RandomLiHeH[AD]->Scale(LiHeH[AD]->Integral()/RandomLiHeH[AD]->Integral());
                
            }
#endif

        }
        
        TFile* SaveLiHe = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/LiHeDistortions.root")).c_str(),"recreate");
        for(Int_t AD=0;AD<NADs;AD++)
        {
            LiHeH[AD]->Write(Form("Nominal LiHe%d period%d",AD,week));
            RandomLiHeH[AD]->Write(Form("LiHeAD%d period%d",AD,week));
        }

        func_LiHe->Write();
      
        delete func_LiHe;
        
        delete SaveLiHe;
        
        
    }
    if(DistortFastNeutronsMatrix)//Distort FN shape for all ADs in the same hall
    {
        TF1* func_FN1;
        TF1* func_FN2;
        TF1* func_FN3;
        
        if(Analysis)//Hydrogen
        {
            func_FN1 = new TF1("func_FN1","pow(x,[0])",InitialVisibleEnergy,FinalVisibleEnergy);
            func_FN2 = new TF1("func_FN2","pow(x,[0])",InitialVisibleEnergy,FinalVisibleEnergy);
            func_FN3 = new TF1("func_FN3","pow(x,[0])",InitialVisibleEnergy,FinalVisibleEnergy);

            Double_t FNError;
            
            FNError = 0.036*rand->Gaus(0,1);
            func_FN1->SetParameter(0,-0.639+FNError);//From Xiang Pan Fast Neutron = E^{-P} where P = 0.639±0.036
            
            FNError = 0.036*rand->Gaus(0,1);
            func_FN2->SetParameter(0,-0.639+FNError);//From Xiang Pan Fast Neutron = E^{-P} where P = 0.639±0.036
            
            FNError = 0.036*rand->Gaus(0,1);
            func_FN3->SetParameter(0,-0.639+FNError);//From Xiang Pan Fast Neutron = E^{-P} where P = 0.639±0.036
            
            for(Int_t AD=0;AD<NADs;AD++)
            {
                if(AD<ADsEH1)
                {
                    RandomFastNeutronsH[AD]->Reset();
                    RandomFastNeutronsH[AD]->Add(func_FN1);
                    RandomFastNeutronsH[AD]->Scale(FastNeutronsH[AD]->Integral()/RandomFastNeutronsH[AD]->Integral());
                }
                else if(AD<(ADsEH1+ADsEH2))
                {
                    RandomFastNeutronsH[AD]->Reset();
                    RandomFastNeutronsH[AD]->Add(func_FN2);
                    RandomFastNeutronsH[AD]->Scale(FastNeutronsH[AD]->Integral()/RandomFastNeutronsH[AD]->Integral());
                }
                else
                {
                    RandomFastNeutronsH[AD]->Reset();
                    RandomFastNeutronsH[AD]->Add(func_FN3);
                    RandomFastNeutronsH[AD]->Scale(FastNeutronsH[AD]->Integral()/RandomFastNeutronsH[AD]->Integral());
                }
            }
        }
        else//Gadolinium
        {
            func_FN1 = new TF1("func_FN1","[0]/(0.2*pow(x,0.1))+[1]",InitialVisibleEnergy,FinalVisibleEnergy);
            
            //  FN shape distortions are applied taking into account correlations in each EH
            GetFastNeutronsDistortionFunction(DistortFN,func_FN1);
            
            for(Int_t AD=0;AD<ADsEH1;AD++)
            {
                RandomFastNeutronsH[AD]->Multiply(func_FN1);
                RandomFastNeutronsH[AD]->Scale(FastNeutronsH[AD]->Integral()/RandomFastNeutronsH[AD]->Integral());
            }
            func_FN2 = new TF1("func_FN2","[0]/(0.2*pow(x,0.1))+[1]",InitialVisibleEnergy,FinalVisibleEnergy);
            
            GetFastNeutronsDistortionFunction(DistortFN,func_FN2);
            
            for(Int_t AD=ADsEH1;AD<(ADsEH1+ADsEH2);AD++)
            {
                RandomFastNeutronsH[AD]->Multiply(func_FN2);
                RandomFastNeutronsH[AD]->Scale(FastNeutronsH[AD]->Integral()/RandomFastNeutronsH[AD]->Integral());
            }
            func_FN3 = new TF1("func_FN3","[0]/(0.2*pow(x,0.1))+[1]",InitialVisibleEnergy,FinalVisibleEnergy);
            
            GetFastNeutronsDistortionFunction(DistortFN,func_FN3);
            
            for(Int_t AD=(ADsEH1+ADsEH2);AD<(ADsEH1+ADsEH2+ADsEH3);AD++)
            {
                RandomFastNeutronsH[AD]->Multiply(func_FN3);
                RandomFastNeutronsH[AD]->Scale(FastNeutronsH[AD]->Integral()/RandomFastNeutronsH[AD]->Integral());
            }
        }
        
        TFile* SaveFN = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/FNDistortions.root")).c_str(),"recreate");
        for(Int_t AD = 0; AD<NADs;AD++)
        {
            FastNeutronsH[AD]->Write(Form("Nominal FN%d period%d",AD,week));
            RandomFastNeutronsH[AD]->Write(Form("FN%d period%d",AD,week));
        }
        func_FN1->Write();
        func_FN2->Write();
        func_FN3->Write();
        
        delete SaveFN;
        
        delete func_FN1;
        delete func_FN2;
        delete func_FN3;
    }
    if(DistortAmCMatrix)//  Distort AmC shape for all ADs in the same way
    {
        TF1* func_AmC = (TF1*)NominalAmCF->Clone("amccopy");
        rand->SetSeed(0);
        Double_t invpar1=1./NominalAmCF->GetParameter(1)*rand->Gaus(1,DistortAmC);
        func_AmC->SetParameter(1,1./invpar1);
        
        for(Int_t AD=0;AD<NADs;AD++)
        {
            RandomAmCH[AD]->Reset();
            RandomAmCH[AD]->Add(func_AmC);
            RandomAmCH[AD]->Scale(AmCH[AD]->Integral()/RandomAmCH[AD]->Integral());
        }
        
        TFile* SaveAmC = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/AmCDistortions.root")).c_str(),"recreate");
        for(Int_t AD=0;AD<NADs;AD++)
        {
            AmCH[AD]->Write(Form("Nominal AmC%d period%d",AD,week));
            RandomAmCH[AD]->Write(Form("AmC%d period%d",AD,week));
        }
        func_AmC->Write();
        delete SaveAmC;
        delete func_AmC;
    }
    
    if(VaryAccidentalMatrix)
    {
        TFile* SaveVaryACCCanvasF = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/AccidentalVariations.root")).c_str(),"recreate");
        
        for (Int_t AD =0; AD<NADs; AD++)
        {
            RandomAccidentalsH[AD]->Scale(1.*ScaleFactorAccidental[AD]);
            AccidentalsH[AD]->Write(Form("Nominal Accidental%d period%d",AD,week));
            RandomAccidentalsH[AD]->Write(Form("Accidental%d period%d",AD,week));
        }
        delete SaveVaryACCCanvasF;
    }
    if(VaryLiHeMatrix)
    {
        
        TFile* SaveVaryLiHeCanvasF = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/LiHeVariations.root")).c_str(),"recreate");
        for (Int_t AD =0; AD<NADs; AD++)
        {
            RandomLiHeH[AD]->Scale(1.*ScaleFactorLiHe[AD]);
            LiHeH[AD]->Write(Form("Nominal LiHe%d period%d",AD,week));
            RandomLiHeH[AD]->Write(Form("LiHe%d period%d",AD,week));
        }
        delete SaveVaryLiHeCanvasF;
    }
    if(VaryFastNeutronsMatrix)
    {
        
        TFile* SaveVaryFNCanvasF = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/FNVariations.root")).c_str(),"recreate");
        for (Int_t AD =0; AD<NADs; AD++)
        {
            RandomFastNeutronsH[AD]->Scale(1.*ScaleFactorFastNeutrons[AD]);
            FastNeutronsH[AD]->Write(Form("Nominal FN%d period%d",AD,week));
            RandomFastNeutronsH[AD]->Write(Form("FNAD%d period%d",AD,week));
        }
        delete SaveVaryFNCanvasF;
    }
    if(VaryAmCMatrix)
    {
        TFile* SaveVaryAmCCanvasF = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/AmCVariations.root")).c_str(),"recreate");
        for (Int_t AD =0; AD<NADs; AD++)
        {
            RandomAmCH[AD]->Scale(1.*ScaleFactorAmC[AD]);
            AmCH[AD]->Write(Form("Nominal AmCH%d period%d",AD,week));
            RandomAmCH[AD]->Write(Form("AmCAD%d period%d",AD,week));
        }
        delete SaveVaryAmCCanvasF;
    }
    
    std::cout <<  "********************************************************************************************************" << std::endl;
    
    for(Int_t AD=0; AD<NADs; AD++)
    {
        if(VaryAccidentalMatrix||VaryLiHeMatrix||VaryFastNeutronsMatrix||VaryAmCMatrix)
        {
            std::cout << "\n \t Background Scale Factors:" << std::endl;
            std::cout << "\t Acc LiHe FN AmC" << std::endl;
            std::cout << "AD: "<<AD << "\t " << ScaleFactorAccidental[AD] << "   " << ScaleFactorLiHe[AD] <<"   " <<  ScaleFactorFastNeutrons[AD] <<"   " <<  ScaleFactorAmC[AD] << std::endl;
        }
        
        //  Now BackgroundSpectrumH holds random backgrounds, either varying the rate or the shape
        BackgroundSpectrumH[AD]=(TH1D*)RandomAccidentalsH[AD]->Clone();
        BackgroundSpectrumH[AD]->Add(RandomLiHeH[AD]);
        BackgroundSpectrumH[AD]->Add(RandomFastNeutronsH[AD]);
        BackgroundSpectrumH[AD]->Add(RandomAmCH[AD]);
        BackgroundSpectrumH[AD]->SetTitle(Form("Random Background Spectrum_%d",AD));
        
        delete RandomAccidentalsH[AD];
        delete RandomLiHeH[AD];
        delete RandomFastNeutronsH[AD];
        delete RandomAmCH[AD];
    }
}

#ifdef UseLiHeToyMC
void Oscillation :: GetDistortionFunction(Double_t amount,TH1F* DistortionFunc)
{
    TH1F* LocalCopy;
    TFile* m_file_distortLi9Bg = new TFile("./Inputs/GdInputs/8he9li_distort_neutron100_alpha100_frac0.1_N250.root","READ");
    TTree* m_tree_distortLi9Bg = (TTree*)m_file_distortLi9Bg->Get("tr_distort");
    m_tree_distortLi9Bg->SetBranchAddress("h_distort",&LocalCopy);
    Int_t m_entries_distortLi9Bg = (Int_t)m_tree_distortLi9Bg->GetEntries();
    std::cout << "The distortion tree has " << m_entries_distortLi9Bg << " entries" << std::endl;
    Int_t entry = rand->Uniform(0,m_entries_distortLi9Bg);
    std::cout << "Reading Li He entry : " << entry << std::endl;
    m_tree_distortLi9Bg->GetEntry(entry);
    func_LiHe = (TH1F*)LocalCopy->Clone();
  //  delete m_file_distortLi9Bg;//deletes the file and the tree
}
#else
void Oscillation :: GetDistortionFunction(Double_t amount,TF1* DistortionFunc)
{
    rand->SetSeed(0);
    Double_t slope = amount*rand->Gaus(0,1);
    
    Double_t anchor_point=3.5;
    //want offset to be set by requiring func at anchor point to be 1
    Double_t offset = (1-slope*anchor_point);
    
    DistortionFunc->SetParameter(0,offset);
    DistortionFunc->SetParameter(1,slope);
}
#endif

void Oscillation :: GetFastNeutronsDistortionFunction(Double_t amount,TF1* DistortionFunc)
{
    rand->SetSeed(0);
    Double_t scaling =amount*rand->Gaus(0,1);
    
    DistortionFunc->SetParameter(0,scaling);
    //set offset so that func(FinalEnergy)=1;
    DistortionFunc->SetParameter(1,1-1*DistortionFunc->Eval(10));//I'm using 10 instead of 12 because when I fit the FN I use the range 10-100 MeV.
}


void Oscillation :: LoadNominalBackgrounds()
{
    Char_t BackgroundFile[100];
    
    if(Analysis)
    {
        sprintf(BackgroundFile, "./BackgroundSpectrum/HBackground/Backgrounds.root");
    }
    else
    {
        sprintf(BackgroundFile, "./BackgroundSpectrum/GDBackground/Backgrounds.root");
    }
    
    TFile* BackgroundsF = new TFile(BackgroundFile);
    
    for (Int_t AD =0; AD<NADs; AD++)
    {
        AccidentalsH[AD]=(TH1D*)gDirectory->Get(Form("Accidentals_AD%i",AD));
        LiHeH[AD]=(TH1D*)gDirectory->Get(Form("LiHe_AD%i",AD));
        FastNeutronsH[AD]=(TH1D*)gDirectory->Get(Form("FN_AD%i",AD));
        AmCH[AD]=(TH1D*)gDirectory->Get(Form("AmC_AD%i",AD));
        NominalAmCF = (TF1*)gDirectory->Get("AmCFunc");
        
        //Nominal background spectrum
        BackgroundSpectrumH[AD] = new TH1D(Form("Background Spectrum_%d",AD),Form("Background Spectrum_%d",AD),n_evis_bins,evis_bins);
        BackgroundSpectrumH[AD]->Add(AccidentalsH[AD]);
        BackgroundSpectrumH[AD]->Add(LiHeH[AD]);
        BackgroundSpectrumH[AD]->Add(FastNeutronsH[AD]);
        BackgroundSpectrumH[AD]->Add(AmCH[AD]);
    }
    
    #ifdef PrintEps
        TCanvas* AccidentalsC = new TCanvas("AccidentalsC","AccidentalsC");
        AccidentalsC->Divide(NADs/2,2);
        for(Int_t AD=0; AD<NADs; AD++)
        {
            AccidentalsC->cd(AD+1);
            AccidentalsH[AD]->Draw();
        }
        AccidentalsC->Update();
        
        AccidentalsC->Print(("./Images/"+AnalysisString+"/BackgroundVariations/NominalAccidentals.eps").c_str(),".eps");
        
        delete AccidentalsC;
        
        TCanvas* LiHeC = new TCanvas("LiHeC","LiHeC");
        LiHeC->Divide(NADs/2,2);
        for(Int_t AD=0; AD<NADs; AD++)
        {
            LiHeC->cd(AD+1);
            LiHeH[AD]->Draw();
            
        }
        LiHeC->Update();
        
        LiHeC->Print(("./Images/"+AnalysisString+"/BackgroundVariations/NominalLiHe.eps").c_str(),".eps");
        
        delete LiHeC;
        
        TCanvas* AmCC = new TCanvas("AmCC","AmCC");
        AmCC->Divide(NADs/2,2);
        for(Int_t AD=0; AD<NADs; AD++)
        {
            AmCC->cd(AD+1);
            AmCH[AD]->Draw();
        }
        AmCC->Update();
        
        AmCC->Print(("./Images/"+AnalysisString+"/BackgroundVariations/NominalAmC.eps").c_str(),".eps");
        
        delete AmCC;
        
        TCanvas* FastNeutronsC = new TCanvas("FastNeutronsC","FastNeutronsC");
        FastNeutronsC->Divide(NADs/2,2);
        for(Int_t AD=0; AD<NADs; AD++)
        {
            FastNeutronsC->cd(AD+1);
            FastNeutronsH[AD]->Draw();
        }
        FastNeutronsC->Update();
        
        FastNeutronsC->Print(("./Images/"+AnalysisString+"/BackgroundVariations/NominalFastNeutrons.eps").c_str(),".eps");
        
        delete FastNeutronsC;
        
    #endif
    delete BackgroundsF;
}





