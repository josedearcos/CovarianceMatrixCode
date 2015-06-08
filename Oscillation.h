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

//Comment the plot options to speed it up.

//#define PlotOscillationFunctions
#define PlotExtrapolationFactors

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
    TH1F* Hist_LiHe;
    TH1F* LocalCopy;
    TFile* m_file_distortLi9Bg;
    TTree* m_tree_distortLi9Bg;
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
    
    bool IHEPReactorModel;
    std::string ResponseDirectory;
    Double_t s22t12;
    Double_t s22t13;
    
    Int_t week;
    Double_t ADdistances[MaxDetectors][NReactors];
    
    std::vector<Double_t> Norma;
    
    Double_t FluxFraction[MaxDetectors/2][NReactors][MaxNbins][VolumeX][VolumeY];
    Double_t Extrapolation[MaxDetectors/2][NReactors][MaxDetectors/2][MaxNbins][VolumeX][VolumeY];

    Double_t DetectorProtons[MaxDetectors][VolumeX][VolumeY];
    Double_t FullTime[MaxDetectors][MaxPeriods];
    Double_t MultiMuonEff[MaxDetectors][MaxPeriods];
    Double_t FluxFullTime[MaxDetectors][MaxPeriods];
    Double_t DetectorEfficiencyDelayed[MaxDetectors];
    Double_t DetectorEfficiency[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    Double_t NominalDetectorEfficiency[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    Double_t FluxNominalDetectorEfficiency[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
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
    std::string FluxInputS;
    Int_t DataSet;
    //Cell parameters:
    Int_t XCellLimit;
    Int_t YCellLimit;
    
    //Binning parameters:
    Int_t n_evis_bins;
    Int_t n_etrue_bins;
    Double_t InitialEnergy;
    Double_t FinalVisibleEnergy;
    Double_t InitialVisibleEnergy;
    
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.
    
    //Background errors:
    Double_t AccidentalError[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    Double_t LiHeError[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    Double_t FastNeutronsError[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    Double_t AmCError[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    
    //Distortion function factor:
    Double_t DistortLiHe;
    Double_t DistortFN;
    Double_t DistortAmC;
    
    Double_t ScaleFactorAccidental[MaxDetectors][VolumeX][VolumeY];
    Double_t ScaleFactorLiHe[MaxDetectors][VolumeX][VolumeY];
    Double_t ScaleFactorFastNeutrons[MaxDetectors][VolumeX][VolumeY];
    Double_t ScaleFactorAmC[MaxDetectors][VolumeX][VolumeY];
    
    //Histograms
    
    TH1D* BackgroundSpectrumH[MaxDetectors][VolumeX][VolumeY];
    TH1D* NearBackgroundSpectrumH[MaxNearDetectors][VolumeX][VolumeY];
    TH1D* FarBackgroundSpectrumH[MaxDetectors][VolumeX][VolumeY];
    
    TH1D* AccidentalsH[MaxDetectors][VolumeX][VolumeY];
    TH1D* LiHeH[MaxDetectors][VolumeX][VolumeY];
    TH1D* FastNeutronsH[MaxDetectors][VolumeX][VolumeY];
    TH1D* AmCH[MaxDetectors][VolumeX][VolumeY];
    
    TH1D* RandomAccidentalsH[MaxDetectors][VolumeX][VolumeY];
    TH1D* RandomLiHeH[MaxDetectors][VolumeX][VolumeY];
    TH1D* RandomFastNeutronsH[MaxDetectors][VolumeX][VolumeY];
    TH1D* RandomAmCH[MaxDetectors][VolumeX][VolumeY];

    TH1D* NearSpectrumH[MaxNearDetectors][MatrixBins];//
    TH1D* FarHallSpectrumH[MaxFarDetectors][MaxNearDetectors][MatrixBins];//

    TH1D* NearSpectrumFractionH[MaxNearDetectors][NReactors][MatrixBins][VolumeX][VolumeY];//
    TH1D* FarHallSpectrumFractionH[MaxNearDetectors][NReactors][MaxFarDetectors][MatrixBins][VolumeX][VolumeY];//
    TH1D* FarHallCellSpectrumH[MaxFarDetectors][MaxNearDetectors][MatrixBins][VolumeX][VolumeY];//

    TH1D* FluxFractionH[MaxNearDetectors][NReactors][VolumeX][VolumeY];//
    std::vector<TH1D*> ExtrapolationH;//
    
    TH1D* TrueADSpectrumH[MaxNearDetectors][MatrixBins][VolumeX][VolumeY];//
    //    std::vector<TH1D*> CheckVisADSpectrumH;
    TH1D* ADSpectrumVisH[MaxNearDetectors][VolumeX][VolumeY];//
    TH1D* FluxH[NReactors][MaxDetectors][MaxPeriods][VolumeX][VolumeY];//
    TH1D* InclusiveFluxH[NReactors][MaxDetectors][VolumeX][VolumeY];//
    std::vector<TH1D*> OscillationFunction;//
    TH1D* FDBH;//
    TH1D* FLOH;//
    
    TH2D* TransEnergyMatrix[MaxDetectors][VolumeX][VolumeY];
    TH2D* EnergyMatrix;
    TH2D* LBNLEnergyMatrix;
    //Files
    Char_t name[20];
    const Double_t colors[6]= {1,2,4,8,6,11};
    const std::string labels[6]={"D1","D2","L1","L2","L3","L4"};

    bool flagDelete;
    //Load data
    void LoadToyHistograms(Int_t);
    void LoadNearData(Int_t);
    
    //Functions to vary background shapes. So far the same ones than LBNL.

    void GetDistortionFunction(Double_t,TF1*);

    void GetFastNeutronsDistortionFunction(Double_t,TF1*);
    void FluctuateBackgrounds(Int_t);
    void LoadNominalBackgrounds();
    void LoadBackgrounds(Int_t);
    void SaveBackgrounds();

    //Calculate oscillation
    void PrintFluxFraction();
    void CalculateFluxFraction(Int_t,Int_t);
    void GetExtrapolation(Int_t,Int_t);

    void NearSpectrumFraction(Int_t,Int_t);
    void FarSpectrumPrediction(Int_t,Int_t);
    void SaveOscillationFunction();
    void CalculateOscillationFunction();
    
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
    TH1D* GetOscillatedADSpectrum(Int_t,Int_t,Int_t);
    void SetDelayedEfficiency(Double_t[]);
    void SetFluxHistograms(TH1D*,Int_t,Int_t,Int_t,Int_t,Int_t);
    TH1D* GetFluxHisto(Int_t,Int_t,Int_t,Int_t,Int_t);
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
        FluxInputS = "HInputs";
        AnalysisString = "Hydrogen";
        XCellLimit = VolumeX;
        YCellLimit = VolumeY;
    }
    else
    {
        FluxInputS="GdInputs";
        AnalysisString = "Gadolinium";
        XCellLimit = 1;
        YCellLimit = 1;
    }
    DataSet = Nom->GetDataSet();
    
    Nweeks = Nom->GetWeeks();
    NReactorPeriods = Nom->GetNReactorPeriods();
    IHEPReactorModel = Nom->UsingIHEPReactorModel();
    NADs = Nom->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for(Int_t reactor = 0; reactor<NReactors; reactor++)
        {
            ADdistances[AD][reactor] = Nom->GetDistances(AD,reactor)/1000;//IN KILOMETERS
        }
        
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
#ifdef UseVolumes //To use 2 volumes, otherwise 100 cells.

                DetectorProtons[AD][idx][idy] = Nom->GetDetectorProtons(AD,idx);
                
#else
                if((idy==0||idy==(YCellLimit-1)||idx>=(XCellLimit-4)))//nH LS
                {
                    if(Analysis)//Hydrogen LS
                    {
                        DetectorProtons[AD][idx][idy] = Nom->GetDetectorProtons(AD,1)/(XCellLimit+YCellLimit+(YCellLimit-2)*4);//52 cells (10+10+(10-2)+(10-2)+(10-2)+(10-2)) Share uniformly the mass
                        
                    }
                    else//nGd only 1 volume
                    {
                        DetectorProtons[AD][idx][idy] = Nom->GetDetectorProtons(AD,0);
                    }
                }
                else//GdLS
                {
                    DetectorProtons[AD][idx][idy] = Nom->GetDetectorProtons(AD,0)/((XCellLimit-4)*(YCellLimit-2));//48 cells, Share uniformly the mass
                }
#endif
                //For fitter:
                for (Int_t week = 0; week <Nweeks; week++)
                {
                    NominalDetectorEfficiency[AD][week][idx][idy] = Nom->GetDetectorEfficiency(AD,week,idx,idy,0);
                    DetectorEfficiency[AD][week][idx][idy] = NominalDetectorEfficiency[AD][week][idx][idy];
                    FullTime[AD][week] = Nom->GetFullTime(AD,week,0);
                    MultiMuonEff[AD][week] = Nom->GetMultiMuonEff(AD,week,0);
                }
                //For flux:
                for (Int_t week = 0; week <NReactorPeriods; week++)
                {
                    FluxFullTime[AD][week] = Nom->GetFullTime(AD,week,1);
                    FluxNominalDetectorEfficiency[AD][week][idx][idy] = Nom->GetMultiMuonEff(AD,week,1);//We only need the weekly dependent efficiencies, we could use GetDetectorEfficiency also since the global efficiencies cancel out when we weight the weekly data, but this is more "correct"
                }
            }
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
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    AccidentalError[AD][week][idx][idy]=Nom->GetAccidentalError(AD,week,idx,idy);
                    LiHeError[AD][week][idx][idy]=Nom->GetLiHeError(AD,week,idx,idy);
                    FastNeutronsError[AD][week][idx][idy]=Nom->GetFNError(AD,week,idx,idy);
                    AmCError[AD][week][idx][idy]=Nom->GetAmCError(AD,week,idx,idy);
                }
            }
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
    IHEPReactorModel = OData->UsingIHEPReactorModel();

    Analysis = OData->GetAnalysis();
    
    if(Analysis)
    {
        FluxInputS = "HInputs";
        AnalysisString = "Hydrogen";
        XCellLimit = VolumeX;
        YCellLimit = VolumeY;
    }
    else
    {
        FluxInputS="GdInputs";
        AnalysisString = "Gadolinium";
        XCellLimit = 1;
        YCellLimit = 1;
    }
    DataSet = OData->GetDataSet();
    
    Nweeks = OData->GetWeeks();
    NReactorPeriods = OData->GetNReactorPeriods();

    std::cout << "PERIODS AND WEEKS BEFORE CONSTRUCTOR: :" <<  NReactorPeriods << " " << Nweeks << std::endl;

    NADs = OData->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for(Int_t reactor = 0; reactor<NReactors; reactor++)
        {
            ADdistances[AD][reactor] = OData->GetDistances(AD,reactor)/1000;//IN KILOMETERS
        }
        
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                
#ifdef UseVolumes //To use 2 volumes, otherwise 100 cells.
                
                DetectorProtons[AD][idx][idy] = OData->GetDetectorProtons(AD,idx);
                
#else
                if((idy==0||idy==(YCellLimit-1)||idx>=(XCellLimit-4)))//nH LS
                {
                    if(Analysis)//Hydrogen LS
                    {
                        DetectorProtons[AD][idx][idy] = OData->GetDetectorProtons(AD,1)/(XCellLimit+YCellLimit+(YCellLimit-2)*4);//52 cells (10+10+(10-2)+(10-2)+(10-2)+(10-2)) Share uniformly the mass
                        
                    }
                    else//nGd only 1 volume
                    {
                        DetectorProtons[AD][idx][idy] = OData->GetDetectorProtons(AD,0);
                    }
                }
                else//GdLS
                {
                    DetectorProtons[AD][idx][idy] = OData->GetDetectorProtons(AD,0)/((XCellLimit-4)*(YCellLimit-2));//48 cells, Share uniformly the mass
                }
#endif
                
                //For fitter:
                for (Int_t week = 0; week <Nweeks; week++)
                {
                    NominalDetectorEfficiency[AD][week][idx][idy] = OData->GetDetectorEfficiency(AD,week,idx,idy,0);
                    DetectorEfficiency[AD][week][idx][idy]=NominalDetectorEfficiency[AD][week][idx][idy];
                    FullTime[AD][week] = OData->GetFullTime(AD,week,0);
                    MultiMuonEff[AD][week] = OData->GetMultiMuonEff(AD,week,0);
                }
                //For flux:
                for (Int_t week = 0; week <NReactorPeriods; week++)
                {
                    FluxNominalDetectorEfficiency[AD][week][idx][idy] = OData->GetMultiMuonEff(AD,week,1);
                    FluxFullTime[AD][week]  = OData->GetFullTime(AD,week,1);

                }
            }
        }

        //Commented out because I created a fake file with weekly data
//        if(Analysis)
//        {
//            std::cout << "REMINDER: NEED A WEEKLY EFFICIENCY, FULL TIME, BACKGROUND ETC INFO FILE FOR THE NH ANALYSIS" << std::endl;
//            NReactorPeriods = 1;//This is temporary, otherwise I cannot run the cross-check, the NReactorPeriods is 1 because I have hardcoded the inclusive fulltime efficiencies etc, I need a file with all the weekly information before deleting this line.
//            Nweeks =1;
//        }


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
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    AccidentalError[AD][week][idx][idy]=OData->GetAccidentalError(AD,week,idx,idy);
                    LiHeError[AD][week][idx][idy]=OData->GetLiHeError(AD,week,idx,idy);
                    FastNeutronsError[AD][week][idx][idy]=OData->GetFNError(AD,week,idx,idy);
                    AmCError[AD][week][idx][idy]=OData->GetAmCError(AD,week,idx,idy);
                }
            }
        }
    }
    
    std::cout << "PERIODS AND WEEKS AFTER CONSTRUCTOR: :" <<  NReactorPeriods << " " << Nweeks << std::endl;
  
}


Oscillation :: ~Oscillation()
{
    if(flagDelete)
    {
        for (Int_t ad = 0; ad<NADs; ad++)//far
        {
            for (Int_t reactor = 0; reactor<NReactors; reactor++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        if(Nweeks==1)
                        {
                            delete InclusiveFluxH[reactor][ad][idx][idy];
                        }
                        else
                        {
                            delete FluxH[reactor][ad][week][idx][idy];
                        }
                    }
                }
            }
        }
    }
    else
    {
        for(Int_t AD=0;AD<NADs;AD++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    delete BackgroundSpectrumH[AD][idx][idy];
                    delete AccidentalsH[AD][idx][idy];
                    delete LiHeH[AD][idx][idy];
                    delete AmCH[AD][idx][idy];
                    delete FastNeutronsH[AD][idx][idy];
                }
            }
        }

        delete NominalAmCF;

        
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                for (Int_t reactor = 0; reactor<NReactors; reactor++)
                {
                    for (Int_t ad = 0; ad<NADs; ad++)
                    {
                        delete FluxH[reactor][ad][week][idx][idy];
                    }
                    for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
                    {
                        delete FluxFractionH[near][reactor][idx][idy];
                    }
                }
            }
        }
        


        for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
        {
            for(Int_t j = 0; j < n_evis_bins; j ++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        delete TrueADSpectrumH[near][j][idx][idy];
                    }
                }
            }
        }
        
        for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
        {
            for (Int_t reactor = 0; reactor < NReactors; reactor++)
            {
#ifdef PlotOscillationFunctions
                delete OscillationFunction[near*NReactors+reactor];
#endif
//                for (Int_t far = 0; far<(ADsEH3); far++)//far
//                {
//                    delete ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors];
//                }
            }
        }
        
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                for (Int_t far = 0; far<(ADsEH3); far++)//far
                {
                    delete FarBackgroundSpectrumH[far][idx][idy];
                }
                for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)//near
                {
                    delete NearBackgroundSpectrumH[near][idx][idy];
                }
            }
        }
        for (Int_t far = 0; far<(ADsEH3); far++)//far
        {
            for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)//near
            {
                for(Int_t j = 0; j < n_evis_bins; j ++)
                {
                    delete FarHallSpectrumH[far][near][j];
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
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    DetectorEfficiency[AD][week][idx][idy] = Efficiency[AD+NADs*week+NADs*Nweeks*idx+NADs*Nweeks*XCellLimit*idy];
                    std::cout << "\t \t \t Detector Efficiency in AD" << AD << " ,week: " << week << " Cell: " << idx << "," << idy << " inside Oscillation.h: " << DetectorEfficiency[AD][week][idx][idy] << std::endl;
                }
            }
            
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
    
    std::cout << "\t \t \t Sin22t12 in Oscillation.h is: " << s22t12 << std::endl;
    std::cout << "\t \t \t Sin22t13 in Oscillation.h is: " << s22t13 << std::endl;
    std::cout << "\t \t \t Dm2_ee in Oscillation.h is: " << dm2ee << std::endl;
    
#ifdef PlotOscillationFunctions
    //Save oscillation functions for demonstration purposes.
    CalculateOscillationFunction();
    SaveOscillationFunction();
#endif
    
    //resize has to be done only once, outside the loops
    Norma.resize(n_etrue_bins*(ADsEH1+ADsEH2));

#ifdef PlotExtrapolationFactors
    ExtrapolationH.resize(ADsEH3*NReactors*(ADsEH1+ADsEH2));
#endif
    
    for(Int_t idx=0; idx<XCellLimit; idx++)
    {
        for(Int_t idy=0; idy<YCellLimit; idy++)
        {
            CalculateFluxFraction(idx,idy);// If we only use muon and multi eff this doesn't depend on the volume, nevertheless the global efficiency of each volume is different, but the definition of superflux removes the constants in time and probably this is superfluous and could be eliminated. Since it doesn't hurt (just 1 more loop run) and it's a waste of time to change back the code, I will leave it like this at the moment, but removing the idx,idy dependency could be a small improvement on performance.
            GetExtrapolation(idx,idy);//Same here
            NearSpectrumFraction(idx,idy);
            FarSpectrumPrediction(idx,idy);
        }
    }
#ifdef PrintEps
    PrintFluxFraction();
#endif
    std::cout << "\t Oscillation calculated" << std::endl;
    
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//LOAD DISTANCES BETWEEN ADS AND CORES

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//ANTINEUTRINO SURVIVAL PROBABILITY

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Calculates oscillation probability

Double_t Oscillation :: OscProb(Double_t L, Double_t E, Double_t S22t13,Double_t Dm2ee)
{
    L = L*1000; //Distances are given in km, the 1.267 factor works for km/gev or m/mev
    
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

void Oscillation :: CalculateOscillationFunction()
{
    OscillationFunction.resize((ADsEH1+ADsEH2)*NReactors);
    
    for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
    {
        for (Int_t k = 0; k<NReactors; k++)
        {
            OscillationFunction[near*NReactors+k] = new TH1D(Form("OscProbFromReactor%itoNear%i",k,near),Form("OscProbFromReactor%itoNear%i",k,near),n_etrue_bins,enu_bins[0],enu_bins[n_etrue_bins]);
            
            for(Int_t pts=0;pts<n_etrue_bins;pts++)
            {
                Energy = TrueADSpectrumH[0][0][0][0]->GetXaxis()->GetBinCenter(pts+1);

                OscillationFunction[near*NReactors+k]->SetBinContent(pts+1,OscProb(ADdistances[near][k],Energy,s22t13,dm2ee));
            }
        }
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

void Oscillation :: CalculateFluxFraction(Int_t idx, Int_t idy)
{
    //TFile* FactorF = new TFile(("./RootOutputs/"+ AnalysisString+ "/FactorStudies.root").c_str(),"recreate");
    
    Char_t filenameFluxFraction[100];
    
    for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
    {
        for(Int_t pts=0;pts<n_etrue_bins;pts++)
        {
            Norma[near+pts*(ADsEH1+ADsEH2)]=0.0;
            
            Energy = TrueADSpectrumH[0][0][0][0]->GetXaxis()->GetBinCenter(pts+1);
            
            for (Int_t reactor = 0; reactor<NReactors; reactor++)
            {
                Norma[near+pts*(ADsEH1+ADsEH2)] = Norma[near+pts*(ADsEH1+ADsEH2)] + ((FluxH[reactor][near][week][idx][idy]->GetBinContent(pts+1)*OscProb(ADdistances[near][reactor],Energy,s22t13,dm2ee)/(ADdistances[near][reactor]*ADdistances[near][reactor])));//  Σ Over Cores of all fractions
                  //std::cout <<  "ENERGY: " << Energy << " POINT " << pts << " NORMA: " << Norma[near+pts*(ADsEH1+ADsEH2)] << std::endl;
            }
        }
        
        for (Int_t reactor = 0; reactor<NReactors; reactor++)
        {
            sprintf(filenameFluxFraction,"FluxFraction factor, Near AD%i, Reactor%i, Week%i, Cell%i,%i", near+1, reactor+1, week+1,idx,idy);
            
            FluxFractionH[near][reactor][idx][idy] = (TH1D*)FluxH[0][0][0][0][0]->Clone(filenameFluxFraction);
            FluxFractionH[near][reactor][idx][idy]->Reset();
            FluxFractionH[near][reactor][idx][idy]->SetTitle(filenameFluxFraction);
            
            for(Int_t pts=0;pts<n_etrue_bins;pts++)
            {
                Energy = TrueADSpectrumH[0][0][0][0]->GetXaxis()->GetBinCenter(pts+1);
                
                if(Norma[near+pts*(ADsEH1+ADsEH2)]!=0.0)
                {
                    FluxFraction[near][reactor][pts][idx][idy] = (FluxH[reactor][near][week][idx][idy]->GetBinContent(pts+1)*OscProb(ADdistances[near][reactor],Energy,s22t13,dm2ee))/(ADdistances[near][reactor]*ADdistances[near][reactor]*Norma[near+pts*(ADsEH1+ADsEH2)]);
                    
                //std::cout << "bin: " << pts << " flux fraction: " << FluxFraction[near][reactor][pts][idx][idy] << " flux: " << FluxH[reactor][near][week][idx][idy]->GetBinContent(pts+1) << " norm: " << Norma[near+pts*(ADsEH1+ADsEH2)] << " osc prob: " << OscProb(ADdistances[near][reactor],Energy,s22t13,dm2ee) << std::endl;
                }
                else
                {
                    std::cout << " NORMA IS 0, THIS MUST BE AN ERROR" << std::endl;
                    
                    exit(EXIT_FAILURE);
                    
                    FluxFraction[near][reactor][pts][idx][idy]=0.0;
                }
                
                FluxFractionH[near][reactor][idx][idy]->SetBinContent(pts+1,FluxFraction[near][reactor][pts][idx][idy]);
            }
            
            //FluxFractionH[near+reactor*(ADsEH1+ADsEH2)]->Write();
        }
    }
    
    std::cout << " CalculateFluxFraction ended " << std::endl;

    //delete FactorF;
}

void Oscillation :: PrintFluxFraction()
{
    TCanvas* FluxFracC = new TCanvas("FluxFractions","FluxFractions",1000,300);
        
    FluxFracC->Divide(ADsEH1+ADsEH2,XCellLimit);
    TLegend *leg = new TLegend(0.1,0.1,0.4,0.5);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    
    for (Int_t idx = 0; idx<XCellLimit; idx++)
    {
        for (Int_t idy = 0; idy<YCellLimit; idy++)
        {
            for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
            {
                FluxFracC->cd(near+idx*(ADsEH1+ADsEH2)+1);
                
                for (Int_t reactor = 0; reactor<NReactors; reactor++)
                {
                    FluxFractionH[near][reactor][idx][idy]->SetTitle("Flux Fractions");
                    
                    FluxFractionH[near][reactor][idx][idy]->SetLineColor(colors[reactor]);
                    FluxFractionH[near][reactor][idx][idy]->GetXaxis()->SetTitle("E_{true} (MeV)");
                    FluxFractionH[near][reactor][idx][idy]->GetXaxis()->SetTitleSize(0.04);
                    FluxFractionH[near][reactor][idx][idy]->GetYaxis()->SetTitleSize(0.045);
                    FluxFractionH[near][reactor][idx][idy]->GetXaxis()->SetLabelSize(0.045);
                    FluxFractionH[near][reactor][idx][idy]->GetYaxis()->SetLabelSize(0.045);
                    FluxFractionH[near][reactor][idx][idy]->SetTitle("");
                    
                    sprintf(name,"f_{%i,j}",near+1);
                    FluxFractionH[near][reactor][idx][idy]->GetYaxis()->SetTitle(name);
                    if(reactor==0)
                    {
                        FluxFractionH[near][reactor][idx][idy]->GetYaxis()->SetRangeUser(0,0.5);
                        FluxFractionH[near][reactor][idx][idy]->Draw();
                    }
                    FluxFractionH[near][reactor][idx][idy]->Draw("same");
                    if(near==0&&idx==0)
                    {
                        leg->AddEntry(FluxFractionH[near][reactor][idx][idy],labels[reactor].c_str(),"l");
                        leg->Draw("same");
                    }
                    
                }
            }
        }
    }
    FluxFracC->Print(("./Images/"+AnalysisString+"/FluxFractions.eps").c_str());
    
    delete FluxFracC;
}

void Oscillation :: NearSpectrumFraction(Int_t idx, Int_t idy)
{
    TFile* NearSpectrumF;
    Char_t filenameNear[100];
    Char_t filenameFrac[100];
    
    if(!Analysis)
    {
        Char_t OptionSave[20];
        if(idx==0&&idy==0)
        {
            sprintf(OptionSave,"recreate");
        }
        else
        {
            sprintf(OptionSave,"update");
        }
        
        NearSpectrumF = new TFile(("./RootOutputs/"+ AnalysisString+ "/Spectra/NearSpectrumFraction.root").c_str(),OptionSave);

    }
#ifdef CheckSumTrueSpectrum
    if(idx==Int_t(XCellLimit/2)&&idy==Int_t(YCellLimit/2))
    {
        TH1D* SumNearSpectrumH[(ADsEH1+ADsEH2)][NReactors];
        
        for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
        {
            for (Int_t reactor = 0; reactor<NReactors; reactor++)
            {
                SumNearSpectrumH[near][reactor] = new TH1D(Form("Near_AD%d_Spectrum_from_reactor%d",near+1,reactor+1),Form("Near_AD%d_Spectrum_from_reactor%d",near+1,reactor+1),n_etrue_bins,enu_bins);
            }
        }
    }
#endif
    for(Int_t j = 0; j < n_evis_bins; j++)
    {
#ifdef CheckSumTrueSpectrum
        if(idx==Int_t(XCellLimit/2)&&idy==Int_t(YCellLimit/2))
        {
            TCanvas* NearSpectrumFC = new TCanvas("NearSpectrumFC","NearSpectrumFC");
            
            NearSpectrumFC->Divide(NReactors,(ADsEH1+ADsEH2));
        }
#endif

        for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
        {
            if(!Analysis)
            {
                if(Nweeks==1)
                {
                    sprintf(filenameNear,"AD%i Near Prediction Vis%i, Cell%i,%i", near+1,j+1,idx,idy);
                }
                else
                {
                    sprintf(filenameNear,"AD%i Near Prediction, Week%i Vis%i, Cell%i,%i", near+1,week+1,j+1,idx,idy);
                }
            }
            NearSpectrumH[near][j] = (TH1D*) TrueADSpectrumH[near][j][idx][idy]->Clone();
            NearSpectrumH[near][j]->Reset();
            
            for (Int_t reactor = 0; reactor<NReactors; reactor++)
            {
                NearSpectrumFractionH[near][reactor][j][idx][idy] = (TH1D*)TrueADSpectrumH[near][j][idx][idy]->Clone();
                
                for(Int_t pts=0;pts<TrueADSpectrumH[near][0][0][0]->GetXaxis()->GetNbins();pts++)
                {
                    Double_t SpectrumScaled = (TrueADSpectrumH[near][j][idx][idy]->GetBinContent(pts+1))*FluxFraction[near][reactor][pts][idx][idy];
                    
 //                   std::cout << pts << " " << j << " " << TrueADSpectrumH[near][j]->GetBinContent(pts+1) << " " << FluxFraction[near][reactor][pts][idx][idy] << std::endl;
                    
                    NearSpectrumFractionH[near][reactor][j][idx][idy]->SetBinContent(pts+1, SpectrumScaled);
                }
                
                NearSpectrumH[near][j]->Add(NearSpectrumFractionH[near][reactor][j][idx][idy]);//sum over cores
                
                if(!Analysis)
                {
                    if(Nweeks==1)
                    {
                        sprintf(filenameFrac,"AD%i Near Spectrum fraction from Reactor%i Vis%i Cell%i,%i", near+1, reactor+1,j+1,idx,idy);
                    }
                    else
                    {
                        sprintf(filenameFrac,"AD%i Near Spectrum fraction from Reactor%i, Week%i Vis%i Cell%i,%i", near+1, reactor+1, week+1,j+1,idx,idy);
                    }
                    
                    NearSpectrumFractionH[near][reactor][j][idx][idy]->SetTitle(filenameFrac);
                    
                    NearSpectrumFractionH[near][reactor][j][idx][idy]->Write(filenameFrac);
                }
                
#ifdef CheckSumTrueSpectrum
                if(idx==Int_t(XCellLimit/2)&&idy==Int_t(YCellLimit/2))
                {
                    SumNearSpectrumH[near][reactor]->Add(NearSpectrumFractionH[near][reactor][j][idx][idy]);
                    
#ifdef PrintEps
                    NearSpectrumFC->cd(near+reactor*(ADsEH1+ADsEH2)+1);
                    
                    SumNearSpectrumH[near][reactor]->Draw();
                    
                    TrueADSpectrumH[near][j][idx][idy]->Draw();
#endif
                }
#endif
            }//reactor
            if(!Analysis)
            {
                NearSpectrumH[near][j]->SetTitle(filenameNear);
                NearSpectrumH[near][j]->Write(filenameNear);
            }
        }//visible bins
#ifdef CheckSumTrueSpectrum
        if(idx==Int_t(XCellLimit/2)&&idy==Int_t(YCellLimit/2))
        {
            NearSpectrumFC->Print(Form("./Images/SpectrumFractions/NearTrueSpectrumSumVisible%d.eps",j+1));
            
            delete NearSpectrumFC;
        }
#endif
    }
    
    for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
    {
        for(Int_t j = 0; j < n_evis_bins; j ++)
        {
            delete NearSpectrumH[near][j];
        }
    }
    if(!Analysis)
    {
        delete NearSpectrumF;
    }

#ifdef CheckSumTrueSpectrum
    if(idx==Int_t(XCellLimit/2)&&idy==Int_t(YCellLimit/2))
    {
        for (Int_t near = 0; near<ADsEH1+ADsEH2; near++)
        {
            for (Int_t reactor = 0; reactor<NReactors; reactor++)
            {
                delete SumNearSpectrumH[near][reactor];
            }
        }
    }
#endif
    
    std::cout << "Near Spectrum calculation ended" << idx << idy << std::endl;
}

void Oscillation :: GetExtrapolation(Int_t idx, Int_t idy)
{
  //  TFile* FactorF = new TFile(("./RootOutputs/"+ AnalysisString+ "/FactorStudies.root").c_str(),"update");
    
    Char_t filenameExtrapolation[100];
    
#ifdef PlotExtrapolationFactors
    Int_t cont;
    TCanvas *c2;

    c2 = new TCanvas("Extrapolation","Extrapolation",900,900);
    c2->Divide(3,3);
    cont=0;
    
#endif

    for (Int_t far = 0; far<(ADsEH3); far++)//far
    {
        for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)//near
        {
#ifdef PlotExtrapolationFactors
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
#ifdef PlotExtrapolationFactors
                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors] = (TH1D*)TrueADSpectrumH[0][0][0][0]->Clone(filenameExtrapolation);
                ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->SetTitle(filenameExtrapolation);
#endif
                for(Int_t pts=0;pts<n_etrue_bins;pts++)
                {
                    Energy = TrueADSpectrumH[near][0][0][0]->GetXaxis()->GetBinCenter(pts+1);
                    
                    Extrapolation[near][reactor][far][pts][idx][idy]=
                    (FluxH[reactor][far+(ADsEH1+ADsEH2)][week][idx][idy]->GetBinContent(pts+1)*
                     //(DetectorEfficiency[(far+(ADsEH1+ADsEH2))*Nweeks+week]*FullTime[(far+(ADsEH1+ADsEH2))][week]*DetectorProtons[far+(ADsEH1+ADsEH2)])*
                     OscProb(ADdistances[(far+(ADsEH1+ADsEH2))][reactor],Energy,s22t13,dm2ee)*ADdistances[near][reactor]*ADdistances[near][reactor])/
                    (FluxH[reactor][near][week][idx][idy]->GetBinContent(pts+1)
                    // *(DetectorEfficiency[near*Nweeks+week]*FullTime[near][week]*DetectorProtons[near])
                     *OscProb(ADdistances[near][reactor],Energy,s22t13,dm2ee)*ADdistances[(far+(ADsEH1+ADsEH2))][reactor]*ADdistances[(far+(ADsEH1+ADsEH2))][reactor]);
                    
                   // std::cout << Extrapolation[near][reactor][far][pts][idx][idy] << near << reactor << far << pts << idx << idy << " Flux: " << FluxH[reactor][far+(ADsEH1+ADsEH2)][week][idx][idy]->GetBinContent(pts+1) << std::endl;
#ifdef PlotExtrapolationFactors

                    ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->SetBinContent(pts+1,Extrapolation[near][reactor][far][pts][idx][idy]);
#endif
                }
#ifdef PlotExtrapolationFactors
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
                //ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors]->Write();
                //delete ExtrapolationH[near+reactor*(ADsEH1+ADsEH2)+far*(ADsEH1+ADsEH2)*NReactors];
#endif
            }
        }
    }
    #ifdef PlotExtrapolationFactors
   
    if(idx==1)
    {
        c2->Print(("./Images/"+AnalysisString+"/ExtrapolationFactorsLsVolume.eps").c_str());

    }
    if(idx==0)
    {
        c2->Print(("./Images/"+AnalysisString+"/ExtrapolationFactorsGdLSVolume.eps").c_str());
    }
    
    //just to check that they are identical, (they should)
    delete c2;
    
    #endif
    
    //delete FactorF;
    
    std::cout << "Get Extrapolation ended" << std::endl;
}

void Oscillation :: FarSpectrumPrediction(Int_t idx, Int_t idy)
{
    Char_t filenameFraction[100];
    Char_t filenameCellFar[100];
    Char_t filenameFar[100];
    TFile* FarSpectrumFractionF;
    
    if(!Analysis)
    {
        Char_t OptionFarSave[20];
        
        if(idx==0&&idy==0)
        {
            sprintf(OptionFarSave,"recreate");
        }
        else
        {
            sprintf(OptionFarSave,"update");
        }
        
        FarSpectrumFractionF = new TFile(("./RootOutputs/"+ AnalysisString+ "/Spectra/FarSpectrumFraction.root").c_str(),OptionFarSave);
    }
#ifdef PrintEps
    TCanvas* FarSpectrumFractionC;
#endif
    for(Int_t j = 0; j < n_evis_bins; j ++)
    {
        #ifdef PrintEps
            FarSpectrumFractionC = new TCanvas(Form("FarSpectrumFractionVisible%d",j+1),Form("FarSpectrumFractionVisible%d",j+1));
            FarSpectrumFractionC->Divide(ADsEH1+ADsEH2,ADsEH3);
        #endif
        
        
        for (Int_t near = 0; near<(ADsEH1+ADsEH2); near++)//near
        {
            for (Int_t far = 0; far<(ADsEH3); far++)//far
            {
                if(idx==0&&idy==0)//only create it and reset it the first time, so when I adds all of the spectrum fractions the information remains.
                {
                    FarHallSpectrumH[far][near][j]= (TH1D*) TrueADSpectrumH[near][0][0][0]->Clone();
                    FarHallSpectrumH[far][near][j]->Reset();
                }
                FarHallCellSpectrumH[far][near][j][idx][idy] = (TH1D*) TrueADSpectrumH[near][0][0][0]->Clone();
                FarHallCellSpectrumH[far][near][j][idx][idy]->Reset();
                
                if(!Analysis)//Write without overloading the process
                {
                    if(Nweeks==1)
                    {
                        sprintf(filenameCellFar,"AD%i Far Spectrum prediction from near AD%i Vis %d, Cell%i,%i", far+1, near+1,j+1,idx,idy);
                    }
                    else
                    {
                        sprintf(filenameCellFar,"AD%i Far Spectrum prediction from near AD%i, Week%i Vis%d, Cell%i,%i", far+1, near+1,week+1,j+1,idx,idy);
                    }
                }
                
                for (Int_t reactor = 0; reactor<NReactors; reactor++)
                {
                    
                    sprintf(filenameFraction,"AD%i Far Spectrum fraction from Reactor%i and near AD%i, Week%i Vis%i", far+(ADsEH1+ADsEH2)+1, reactor+1, near+1, week+1,j+1);
                    
                    FarHallSpectrumFractionH[near][reactor][far][j][idx][idy] = (TH1D*)TrueADSpectrumH[near][j][idx][idy]->Clone(filenameFraction);
                    FarHallSpectrumFractionH[near][reactor][far][j][idx][idy]->Reset();
                    
                    for(Int_t pts=0;pts<n_etrue_bins;pts++)
                    {
                        Double_t FarValue = NearSpectrumFractionH[near][reactor][j][idx][idy]->GetBinContent(pts+1)*Extrapolation[near][reactor][far][pts][idx][idy];
                        
                        FarHallSpectrumFractionH[near][reactor][far][j][idx][idy]->SetBinContent(pts+1, FarValue);
                        
                    }
                    
                    // Save predictions:
                    if(!Analysis)//Write without overloading the process
                    {
                        FarHallSpectrumFractionH[near][reactor][far][j][idx][idy]->SetTitle(filenameFraction);

                        FarHallSpectrumFractionH[near][reactor][far][j][idx][idy]->Write(filenameCellFar);
                    }
                    
                    FarHallCellSpectrumH[far][near][j][idx][idy]->Add(FarHallSpectrumFractionH[near][reactor][far][j][idx][idy]);//sum over cores
                    
                    delete FarHallSpectrumFractionH[near][reactor][far][j][idx][idy];
                    
                }//reactor
                
                // Multiply by efficiencies and the corresponding # of protons in each AD, since the detector protons would make the spectrum huge, normalize with respect to AD1 protons.
                
//                FarHallCellSpectrumH[far][near][j][idx][idy]->Scale(DetectorEfficiency[(far+(ADsEH1+ADsEH2))][week][idx][idy]*FullTime[(far+(ADsEH1+ADsEH2))][week]*DetectorProtons[far+(ADsEH1+ADsEH2)][idx][idy]/DetectorProtons[0][idx][idy]);
                
                //Don't do this anymore, right now we correct data and predictions to not to have any detector effects included, also covariance matrices will be produced like that.
                
                //Sum all volumes/cells:
                
                FarHallSpectrumH[far][near][j]->Add(FarHallCellSpectrumH[far][near][j][idx][idy]);//sum over cells

                
                if(!Analysis)//Write without overloading the process
                {
                    FarHallCellSpectrumH[far][near][j][idx][idy]->SetTitle(filenameCellFar);
                    FarHallCellSpectrumH[far][near][j][idx][idy]->Write();

                    if(Nweeks==1)
                    {
                        sprintf(filenameFar,"AD%i Far Spectrum prediction from near AD%i Vis %d", far+1, near+1,j+1);
                    }
                    else
                    {
                        sprintf(filenameFar,"AD%i Far Spectrum prediction from near AD%i, Week%i Vis%d", far+1, near+1,week+1,j+1);
                    }
                    FarHallSpectrumH[far][near][j]->SetTitle(filenameCellFar);
                    FarHallSpectrumH[far][near][j]->Write(filenameFar);
                }
                
                #ifdef PrintEps
                if(idx==0&&idy==0)
                {
                    FarSpectrumFractionC->cd(near+far*(ADsEH1+ADsEH2)+1);
                    FarHallCellSpectrumH[far][near][j][idx][idy]->Draw();
                }
                #endif
                
                delete FarHallCellSpectrumH[far][near][j][idx][idy];
            }//far
            
            for (Int_t reactor = 0; reactor<NReactors; reactor++)
            {
                delete NearSpectrumFractionH[near][reactor][j][idx][idy];
            }
        }//near
        #ifdef PrintEps
          //  FarSpectrumFractionC->Print(Form("./Images/SpectrumFractions/FarSpectrumFractionVisible%d.eps",j+1));
            delete FarSpectrumFractionC;
        #endif

    }//evis bins
    
    if(!Analysis)
    {
        delete FarSpectrumFractionF;
    }
    std::cout << " Far Spectrum calculation ended" << idx << idy << std::endl;
}

TH1D* Oscillation :: GetOscillatedADSpectrum(Int_t near,Int_t far,Int_t j)
{
    return FarHallSpectrumH[far][near][j];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//                                                            LOAD REACTOR TOY DATA

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Loads near hall data
void Oscillation :: LoadToyHistograms(Int_t week)
{
    std::cout << "\t \t Loading toy reactor model" << std::endl;
    Char_t ADdata[200];
    
    if((Sin22t12Matrix||IsotopeMatrix||ReactorPowerMatrix)&&Mode)
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
        
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                ADSpectrumVisH[near][idx][idy] = (TH1D*)gDirectory->Get(Form("Oscillation Prediction AD%d, week%d, cell %d,%d", near+1,week,idx,idy));//It should be calculated with corrected events
                
                //  In oscillationreactor ADs have backgrounds included.
                
                ADSpectrumVisH[near][idx][idy]->Add(BackgroundSpectrumH[near][idx][idy],-1.);
            }
        }
    }
    
    delete NearHallDataF;
    
    #ifdef PrintEps

    TCanvas* cr = new TCanvas("cr","cr", 1200,400);
    cr->Divide(ADsEH1+ADsEH2,XCellLimit*YCellLimit);
    
    if(Analysis)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
                {
                    cr->cd(near+idy*(ADsEH1+ADsEH2)+(ADsEH1+ADsEH2)*YCellLimit*idx+1);

                    ADSpectrumVisH[near][idx][idy]->Draw();
                }
            }
        }
        
    }
    else
    {
        for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            cr->cd(near+1);
            
            ADSpectrumVisH[near][0][0]->Draw("HIST");
        }
    }
    
    cr->Print(("./Images/"+ AnalysisString+ "/"+RandomString+"NearReactorPredictionWithoutBackground.eps").c_str(),".eps");
    delete cr;
    
    #endif

//    TFile* checkADS = new TFile("./RootOutputs/VisibleBinCheck.root","recreate");
//    for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
//    {
//        ADSpectrumVisH[near]->Write();
//    }
//    delete checkADS;
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
    
    for(Int_t idx=0; idx<XCellLimit; idx++)
    {
        for(Int_t idy=0; idy<YCellLimit; idy++)
        {
            for (Int_t near = 0; near < ADsEH1+ADsEH2; near++)
            {
                NearBackgroundSpectrumH[near][idx][idy]=(TH1D*)BackgroundSpectrumH[near*Nweeks+week][idx][idy]->Clone();
            }
            for (Int_t far = 0; far < ADsEH3; far++)
            {
                FarBackgroundSpectrumH[far][idx][idy]=(TH1D*)BackgroundSpectrumH[(far+ADsEH1+ADsEH2)][idx][idy]->Clone();
            }
        }
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
        
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                ADSpectrumVisH[near][idx][idy] = (TH1D*)gDirectory->Get(NearDataSpec);
                
                ADSpectrumVisH[near][idx][idy]->Scale(1./(MultiMuonEff[near][week]*FullTime[near][week]*DetectorProtons[near][idx][idy]/DetectorProtons[0][idx][idy]));//Correct events for efficiencies and calculate it in days. This is done inclusively if NWeeks = 1, or weekly otherwise.

                ////////////////////////////////////////////////////////////////////////////////////////////////////

                //  need to change neardataspec depending on idx,idy cells for nH analysis data when it is available.
                
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            }
        }
        delete NearDataF;
        
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                if(DataSet!=2)//need to rebin files to LBNL binning
                {
                    ADSpectrumVisH[near][idx][idy]=(TH1D*)ADSpectrumVisH[near][idx][idy]->Rebin(n_evis_bins,Form("Rebinned Vis Near%d Data Spectrum, Cell%i,%i",near,idx,idy),evis_bins);
                }
                
                //Substract backgrounds:
                ADSpectrumVisH[near][idx][idy]->Add(NearBackgroundSpectrumH[near][idx][idy],-1);//Substract (varied if necessary) backgrounds from data, which are already corrected events.
                
//Negative bins don't make any sense for data
                for(Int_t i = 0; i< ADSpectrumVisH[near][idx][idy]->GetXaxis()->GetNbins();i++)
                {
                    if(ADSpectrumVisH[near][idx][idy]->GetBinContent(i+1)<0)
                    {
                        ADSpectrumVisH[near][idx][idy]->SetBinContent(i+1,0);
                    }
                }
            }
        }
    }
    
    #ifdef PrintEps

        TCanvas* c1 = new TCanvas("c1","Near Data", 1200,400);
        c1->Divide(ADsEH1+ADsEH2,XCellLimit*YCellLimit);
        
        for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            for(Int_t idx = 0; idx < XCellLimit; idx++)
            {
                for(Int_t idy = 0; idy < YCellLimit; idy++)
                {
                    c1->cd(near+idx*(ADsEH1+ADsEH2)+1);

                        ADSpectrumVisH[near][idx][idy]->Draw();
                }
            }
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
        
        TFile* TransEnergyMatrixDataF;
#ifndef EREC_COMPARISON
        if(Analysis)
        {
             TransEnergyMatrixDataF = new TFile(("./ResponseMatrices/"+ AnalysisString+Form("/NominalResponseMatrix%i_%i.root",n_evis_bins,n_etrue_bins)).c_str());
        }
        else
        {
             TransEnergyMatrixDataF = new TFile(("./ResponseMatrices/"+ AnalysisString+"/NominalResponseMatrix.root").c_str());
        }
#else
         TransEnergyMatrixDataF = new TFile(("./ResponseMatrices/"+ AnalysisString+Form("/E_REC_NominalResponseMatrix%i_%i.root",n_evis_bins,n_etrue_bins)).c_str());
#endif
        
        for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    TransEnergyMatrix[near][idx][idy] = (TH2D*)gDirectory->Get(Form("EnuEvis%i,Cell%i,%i",near+1,idx,idy));
                }
            }
        }
        delete TransEnergyMatrixDataF;
    }
    else
    {
        TFile* TransEnergyMatrixDataF = new TFile((ResponseDirectory).c_str());
        
        for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    TransEnergyMatrix[near][idx][idy] = (TH2D*)gDirectory->Get(Form("EnuEvis%i,Cell%i,%i",near+1,idx,idy));
                }
            }
        }
        delete TransEnergyMatrixDataF;
    }
    
    //If other binnings are used we should include a Rebin method here to match the binning, or recalculate the matrix before running this.
    
    //        TransEnergyMatrix->Rebin();
    
    for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                
                //The scaling in oscillation reactor is already corrected.
                
           //     ADSpectrumVisH[near][idx][idy]->Scale(1./(DetectorEfficiency[near][week][idx][idy]*FullTime[near][week]*DetectorProtons[near][idx][idy]/DetectorProtons[0][idx][idy]));// Remove AD performance included in oscillation reactor and data. This is later included in the far spectrum fraction calculation.
                
                for(Int_t j = 0; j < n_evis_bins; j ++)
                {
                    TrueADSpectrumH[near][j][idx][idy] = new TH1D(Form("TrueAD%i_Energy%i_Data,Cell%i,%i",near,j,idx,idy),Form("TrueAD%i_Energy%i_Data,Cell%i,%i",near,j,idx,idy),n_etrue_bins,enu_bins);
                }
                
                for(Int_t i=0; i<n_etrue_bins; i++)
                {
                    for(Int_t j=0; j<n_evis_bins; j++)
                    {
                        //[True_0 ... True_n] = [0,0]  [...] [0,m]  * [Vis_0]
                        //                      [0,n]  [...] [n,m]     [...]
                        //                                            [Vis_m]
                        
                        // [nx1] = [n x m] x [mx1]
                        TrueADSpectrumH[near][j][idx][idy]->SetBinContent(i+1,TrueADSpectrumH[near][j][idx][idy]->GetBinContent(i+1)+TransEnergyMatrix[near][idx][idy]->GetBinContent(j+1,i+1)*ADSpectrumVisH[near][idx][idy]->GetBinContent(j+1));
                    }
                }
            }
        }
    }
    if(TrueADSpectrumH[0][0][0][0]->GetXaxis()->GetNbins()!=TransEnergyMatrix[0][0][0]->GetYaxis()->GetNbins())
    {
        std::cout << "Binning disagreement between true spectrum and response matrix" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    
    if(ADSpectrumVisH[0][0][0]->GetXaxis()->GetNbins()!=TransEnergyMatrix[0][0][0]->GetXaxis()->GetNbins())
    {
        std::cout << "Binning disagreement between vis spectrum " << ADSpectrumVisH[0][0][0]->GetXaxis()->GetNbins() << " and response matrix" << TransEnergyMatrix[0][0][0]->GetXaxis()->GetNbins() << std::endl;

        exit(EXIT_FAILURE);
    }
    
    for(Int_t near = 0; near < ADsEH1+ADsEH2; near++)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                delete TransEnergyMatrix[near][idx][idy];
            }
        }
    }
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
#ifdef TestAllTheSame

            sprintf(ReactorData,("./RootOutputs/"+ AnalysisString+"/NominalOutputs/AntineutrinoSpectrum.root").c_str()); //The reactor model spectrum, with no reactor on/off information, just assumes the nominal spectrum is on all the time.
#else
#ifdef BlindedAnalysis
            sprintf(ReactorData,Form("./ReactorInputs/WeeklyFlux_%dweek_blinded_inclusive.root",NReactorPeriods));
            std::cout << " NEEDED A BLINDED WEEKLY INCLUSIVE FLUX FILE" << std::endl;
            exit(EXIT_FAILURE);
#else
            sprintf(ReactorData,Form("./ReactorInputs/WeeklyFlux_%dweek_unblinded_inclusive.root",NReactorPeriods));
            std::cout << " NEEDED AN UNBLINDED WEEKLY INCLUSIVE FLUX FILE" << std::endl;
            exit(EXIT_FAILURE);
#endif//BlindedAnalysis
#endif//TestAllTheSame
    }
    else
    {
        
#ifdef BlindedAnalysis
        sprintf(ReactorData, Form("./ReactorInputs/WeeklyFlux_%dweek_blinded.root",NReactorPeriods));
        std::cout << " NEEDED A BLINDED WEEKLY FLUX FILE" << std::endl;
        exit(EXIT_FAILURE);
#else
        sprintf(ReactorData, Form("./ReactorInputs/WeeklyFlux_%dweek_unblinded.root",NReactorPeriods));//Default one
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
#ifdef TestAllTheSame

                    sprintf(ReactorSpectrum,"AntineutrinoSpectrumFromReactor%i",reactor+1);
                    // In this case I don't scale it because the spectrum has not been altered usng the # of protons, technically is spectrum in the reactors not in the ad.
#else
                    sprintf(ReactorSpectrum,"Week%i/%i",period,reactor);
#endif
            }
            else
            {
                sprintf(ReactorSpectrum,"Week%i/%i", period, reactor);//Default one
            }
            
#ifdef NoOscillation
            //To show a flat extrapolation factor and flux fraction
            sprintf(ReactorSpectrum,"AntineutrinoSpectrumFromReactor%i",reactor+1);
#endif
            
            TH1D* NonBinnedFluxH = (TH1D*)gDirectory->Get(ReactorSpectrum);
            
            //        NonBinnedFluxH = (TH1D*)NonBinnedFluxH->Rebin(n_etrue_bins,Form("Flux True Spectrum from reactor%d",reactor+1),enu_bins);
            
            for(Int_t AD=0;AD<NADs;AD++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        
                        FluxH[reactor][AD][period][idx][idy] = new TH1D(Form("FluxH%i_%i_%i,Cell%i_%i",reactor+1,AD+1,period,idx,idy),Form("FluxH%i_%i_%i,Cell%i_%i",reactor+1,AD+1,period,idx,idy),n_etrue_bins,enu_bins);

                        Int_t lastbin = NonBinnedFluxH->GetNbinsX();
                        
                        Double_t last_energy = NonBinnedFluxH->GetBinCenter(lastbin);
                        
                       // std::cout << " Last energy: " << last_energy << " last bin " << lastbin << std::endl;
                        bool First_time = 1;

                        for(Int_t i = 0; i<n_etrue_bins; i++)
                        {
                            Double_t Energy = FluxH[reactor][AD][period][idx][idy]->GetBinCenter(i+1);
                            if(!Analysis)//Interpolation for gd flux files
                            {
                                if(NonBinnedFluxH->Interpolate(Energy)==0)//0
                                {
                                    if(First_time)
                                    {
                                        last_energy = NonBinnedFluxH->GetBinCenter(i); //Last bin with no zero energy
                                        
                                        First_time = 0;
                                    }
                                    FluxH[reactor][AD][period][idx][idy]->SetBinContent(i+1,NonBinnedFluxH->Interpolate(last_energy));
                                    
//                                    std::cout << NonBinnedFluxH->Interpolate(last_energy) << std::endl;
                                }
                            }
                            else//Hydrogen
                            {
                                if(Energy<last_energy)
                                {
                                    FluxH[reactor][AD][period][idx][idy]->SetBinContent(i+1,NonBinnedFluxH->Interpolate(Energy));
                                }
                                else//the flux ends at 8.8 MeV, otherwise the bins are 0 and we get NaNs. This extrapolation produces larger error in the large enery band. Need to do it consistently with other groups.
                                {
                                    FluxH[reactor][AD][period][idx][idy]->SetBinContent(i+1,NonBinnedFluxH->Interpolate(last_energy));
                                }
                            }
                        }
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
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    InclusiveFluxH[reactor][AD][idx][idy] = new TH1D(Form("InclusiveFlux%i_%i,Cell%i_%i",reactor+1,AD+1,idx,idy),Form("InclusiveFlux%i_%i,Cell%i_%i",reactor+1,AD+1,idx,idy),n_etrue_bins,enu_bins);
                    
                    Double_t scalesum=0;
                    
                    for(Int_t week=0;week<NReactorPeriods;week++)
                    {
                        Double_t scale = FluxFullTime[AD][week]*FluxNominalDetectorEfficiency[AD][week][idx][idy];//DetectorEfficiency[AD][week] to use random efficiencies. LBNL pregenerates the flux so the random changes are not contemplated.
                        
                        std::cout << " AD: " << AD << " - reactor: " << reactor << " - Week: " << week << " - Scale: " << scale << " - Full time: " << FluxFullTime[AD][week] << "Cell: " << idx << ", " << idy << " - Efficiency:" << FluxNominalDetectorEfficiency[AD][week][idx][idy] << std::endl;
                        
                        //Something to consider:
                        
                        // If we want to apply random detector efficiencies to the flux in different periods, we should considerate their correlation. If they were correlated we should produce the efficiencies in Prediction.h in order to hold the variations in memory and share a common efficiency. So far I haven't attempted to do a multi period fit, but if you do you should hold this carefully for every systematic.
                        
                        scalesum+=scale;
                        
                        FluxH[reactor][AD][week][idx][idy]->Scale(scale);
                        
                        InclusiveFluxH[reactor][AD][idx][idy]->Add(FluxH[reactor][AD][week][idx][idy]);
                        
                        if(scale==0)
                        {
                            FluxH[reactor][AD][week][idx][idy]->Scale(0);//Some weeks have 0 eff, livetime etc, be careful with NaNs.
                        }
                        else
                        {
                            FluxH[reactor][AD][week][idx][idy]->Scale(1./scale);//if fluxh are used weekly we don't need to correct their efficiencies
                        }
                        //Rebin Flux Spectrum to match the data binning
                        FluxH[reactor][AD][week][idx][idy]=(TH1D*)FluxH[reactor][AD][week][idx][idy]->Rebin(n_etrue_bins,Form("Rebinned Flux Spectrum from reactor%i with week%i ad%i cell%i,%i efficiencies",reactor+1,week+1,AD+1,idx,idy),enu_bins);
                        
                    }
                    if(scalesum==0)
                    {
                        InclusiveFluxH[reactor][AD][idx][idy]->Scale(0);
                        std::cout << " WEEKLY FLUX MUST BE EMPTY, NOT READING THE FILE? " << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    else
                    {
                        InclusiveFluxH[reactor][AD][idx][idy]->Scale(1./scalesum);
                    }
                    //Rebin Flux Spectrum to match the data binning
                    InclusiveFluxH[reactor][AD][idx][idy]=(TH1D*)InclusiveFluxH[reactor][AD][idx][idy]->Rebin(n_etrue_bins,Form("Rebinned Flux Spectrum from reactor%d with ad%d cell%i,%i efficiencies",reactor+1,AD+1,idx,idy),enu_bins);
                }
            }
        }
    }
    
#ifdef PrintEps
    for(Int_t week=0;week<NReactorPeriods;week++)
    {
        TCanvas* WeeklyFluxC = new TCanvas("WeeklyFluxC","WeeklyFluxC");
        
        WeeklyFluxC->Divide(NReactors*XCellLimit,NADs);
        
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                for(Int_t AD=0;AD<NADs;AD++)
                {
                    for(Int_t reactor=0;reactor<NReactors;reactor++)
                    {
                        WeeklyFluxC->cd(reactor+AD*NReactors+1+NReactors*NADs*idx);
                        
                        FluxH[reactor][AD][week][idx][idy]->Draw();
                    }
                }
            }
        }
        
        WeeklyFluxC->Print(("./Images/Reactor/"+AnalysisString+Form("/Weekly%dFluxes.eps",week)).c_str());
        
        delete WeeklyFluxC;
    }
    
    TCanvas* FluxC = new TCanvas("FluxC","FluxC");
    
    FluxC->Divide(NReactors,NADs);
    
    for(Int_t idx=0; idx<XCellLimit; idx++)
    {
        for(Int_t idy=0; idy<YCellLimit; idy++)
        {
            for(Int_t AD=0;AD<NADs;AD++)
            {
                for(Int_t reactor=0;reactor<NReactors;reactor++)
                {
                    FluxC->cd(reactor+AD*NReactors+1+NReactors*NADs*idx);
                    InclusiveFluxH[reactor][AD][idx][idy]->Draw();
                }
            }
        }
    }
    
    FluxC->Print(("./Images/Reactor/"+AnalysisString+"/InclusiveAndWeeklyFluxes.eps").c_str());
    
    delete FluxC;
#endif
    
    delete ReactorDataF;
    
    TFile* FluxF = new TFile(("./Inputs/"+FluxInputS+Form("/SuperFlux%i.root",NReactorPeriods)).c_str(),"recreate");
    
    for(Int_t AD=0;AD<NADs;AD++)
    {
        for(Int_t reactor=0;reactor<NReactors;reactor++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    InclusiveFluxH[reactor][AD][idx][idy]->Write();
                    
                    for(Int_t period = 0; period < NReactorPeriods; period++)
                    {
                        FluxH[reactor][AD][period][idx][idy]->Write();
                    }
                }
            }
            
        }
    }
    
    delete FluxF;
}

void Oscillation :: LoadFluxHisto()//weekly data periods
{
    flagDelete = 1;
    TFile* FluxF;
//    if(IHEPReactorModel)
//    {
//        FluxF = new TFile("./Inputs/SuperFlux20.root");
//        
//        std::cout << "YOU NEED to create a super flux for the IHEP model with 101 periods. Using 20 period file so far" << std::endl;
//    }
//    else
//    {
        FluxF = new TFile(("./Inputs/"+FluxInputS+Form("/SuperFlux%i.root",NReactorPeriods)).c_str());
//    }
    
    if(Nweeks==1)//Inclusive fit
    {
        for(Int_t AD=0;AD<NADs;AD++)
        {
            for(Int_t reactor=0;reactor<NReactors;reactor++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        InclusiveFluxH[reactor][AD][idx][idy]=(TH1D*)FluxF->Get(Form("Rebinned Flux Spectrum from reactor%d with ad%d cell%i,%i efficiencies",reactor+1,AD+1,idx,idy));

                    }
                }
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
                    for(Int_t idx=0; idx<XCellLimit; idx++)
                    {
                        for(Int_t idy=0; idy<YCellLimit; idy++)
                        {
                            FluxH[reactor][AD][period][idx][idy]=(TH1D*)FluxF->Get(Form("Rebinned Flux Spectrum from reactor%i with week%i ad%i cell%i,%i efficiencies",reactor+1,week+1,AD+1,idx,idy));
                        }
                    }
                }
            }
        }
    }
    
    delete FluxF;
}

TH1D* Oscillation :: GetFluxHisto(Int_t reactor, Int_t AD, Int_t week, Int_t idx, Int_t idy)
{
    if(Nweeks == 1)
    {
        return InclusiveFluxH[reactor][AD][idx][idy];
    }
    else
    {
        return FluxH[reactor][AD][week][idx][idy];
    }
}

void Oscillation :: SetFluxHistograms(TH1D* fluxH,Int_t reactor,Int_t AD,Int_t week,Int_t idx, Int_t idy)
{
    FluxH[reactor][AD][week][idx][idy] = (TH1D*)fluxH->Clone();
}

void Oscillation :: SetSin22t12(Double_t S22t12)
{
    s22t12=S22t12;
}

void Oscillation :: SaveBackgrounds(){
    
    std::cout << "\t Saving Backgrounds" << std::endl;

    TFile* BackgroundF1 = new TFile(("./RootOutputs/"+ AnalysisString+ "/Backgrounds/"+RandomString+"Backgrounds.root").c_str(),"recreate");
    for(Int_t idx=0; idx<XCellLimit; idx++)
    {
        for(Int_t idy=0; idy<YCellLimit; idy++)
        {
            for(Int_t near=0; near<(ADsEH1+ADsEH2); near++)
            {
                NearBackgroundSpectrumH[near][idx][idy]->Write((Form("Near AD%i_Volume%i ",near, idx)+RandomString+ Form(" Background Period%i",week)).c_str());
            }
            for(Int_t far=0; far<ADsEH3; far++)
            {
                FarBackgroundSpectrumH[far][idx][idy]->Write((Form("Far AD%i_Volume%i ",far,idx)+RandomString+ Form(" Background Period%i",week)).c_str());
            }
        }
    }
    delete BackgroundF1;
}



//Randomize the events (rate) inside errors taking into account the corresponding correlations between backgrounds and ADs. Also random shape variations are included using distortion functions.
void Oscillation :: FluctuateBackgrounds(Int_t week)
{
    for (Int_t AD =0; AD<NADs; AD++)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                RandomAccidentalsH[AD][idx][idy]=(TH1D*)AccidentalsH[AD][idx][idy]->Clone();
                RandomLiHeH[AD][idx][idy]=(TH1D*)LiHeH[AD][idx][idy]->Clone();
                RandomFastNeutronsH[AD][idx][idy]=(TH1D*)FastNeutronsH[AD][idx][idy]->Clone();
                RandomAmCH[AD][idx][idy]=(TH1D*)AmCH[AD][idx][idy]->Clone();
                
                //  Initialize scale factors
                if(VaryAccidentalMatrix||VaryLiHeMatrix||VaryFastNeutronsMatrix||VaryAmCMatrix)
                {
                    ScaleFactorAccidental[AD][idx][idy]=1.;
                    ScaleFactorLiHe[AD][idx][idy]=1.;
                    ScaleFactorFastNeutrons[AD][idx][idy]=1.;
                    ScaleFactorAmC[AD][idx][idy]=1.;
                }
                
                if(VaryAccidentalMatrix)
                {
                    rand->SetSeed(0);
                    ScaleFactorAccidental[AD][idx][idy]=(1.+(AccidentalError[AD][week][idx][idy]*rand->Gaus(0,1)));
                    //            std::cout << "\t Accidental error " << AccidentalError[AD][week] << std::endl;
                }
                if(VaryLiHeMatrix)
                {
                    rand->SetSeed(0);
                    ScaleFactorLiHe[AD][idx][idy]=(1.+(LiHeError[AD][week][idx][idy]*rand->Gaus(0,1)));
                    //            std::cout <<"\t LiHe error " << LiHeError[AD][week] << std::endl;
                }
                if(VaryFastNeutronsMatrix)
                {
                    rand->SetSeed(0);
                    ScaleFactorFastNeutrons[AD][idx][idy]=(1.+(FastNeutronsError[AD][week][idx][idy]*rand->Gaus(0,1)));
                    //            std::cout <<"\t Fast Neutron error " << FastNeutronsError[AD][week] << std::endl;
                }
                if(VaryAmCMatrix)
                {
                    rand->SetSeed(0);
                    ScaleFactorAmC[AD][idx][idy]=(1.+(AmCError[AD][week][idx][idy]*rand->Gaus(0,1)));
                    //            std::cout <<"\t AmC error " << AmCError[AD][week] << std::endl;
                }
            }
        }
    }
    
    
    for(Int_t idx=0; idx<XCellLimit; idx++)
    {
        for(Int_t idy=0; idy<YCellLimit; idy++)
        {
            //Accidentals uncorrelated
            
            if(VaryLiHeMatrix)
            {   //  LiHe correlated by hall
             
                if(NADs==MaxDetectors)//8AD
                {
                    ScaleFactorLiHe[1][idx][idy]=ScaleFactorLiHe[0][idx][idy];// EH1
                    ScaleFactorLiHe[3][idx][idy]=ScaleFactorLiHe[2][idx][idy];// EH2

                    ScaleFactorLiHe[5][idx][idy]=ScaleFactorLiHe[4][idx][idy];// EH3
                    ScaleFactorLiHe[6][idx][idy]=ScaleFactorLiHe[4][idx][idy];// EH3
                    ScaleFactorLiHe[7][idx][idy]=ScaleFactorLiHe[4][idx][idy];// EH3
                }
                else//6 AD
                {
                    ScaleFactorLiHe[1][idx][idy]=ScaleFactorLiHe[0][idx][idy];// EH1
                    ScaleFactorLiHe[4][idx][idy]=ScaleFactorLiHe[3][idx][idy];// EH3
                    ScaleFactorLiHe[5][idx][idy]=ScaleFactorLiHe[3][idx][idy];// EH3

                }
            }
            if(VaryFastNeutronsMatrix)
            {
                if(NADs==MaxDetectors)//8AD
                {
                    ScaleFactorFastNeutrons[1][idx][idy]=ScaleFactorFastNeutrons[0][idx][idy];// EH1
                    ScaleFactorFastNeutrons[3][idx][idy]=ScaleFactorFastNeutrons[2][idx][idy];// EH2
                    
                    ScaleFactorFastNeutrons[5][idx][idy]=ScaleFactorFastNeutrons[4][idx][idy];// EH3
                    ScaleFactorFastNeutrons[6][idx][idy]=ScaleFactorFastNeutrons[4][idx][idy];// EH3
                    ScaleFactorFastNeutrons[7][idx][idy]=ScaleFactorFastNeutrons[4][idx][idy];// EH3
                }
                else//6 AD
                {
                    ScaleFactorFastNeutrons[1][idx][idy]=ScaleFactorFastNeutrons[0][idx][idy];// EH1
                    ScaleFactorFastNeutrons[4][idx][idy]=ScaleFactorFastNeutrons[3][idx][idy];// EH3
                    ScaleFactorFastNeutrons[5][idx][idy]=ScaleFactorFastNeutrons[3][idx][idy];// EH3
                }
            }
            if(VaryAmCMatrix)
            {   //  AmC all correlated
                ScaleFactorAmC[1][idx][idy]=ScaleFactorAmC[0][idx][idy];
                ScaleFactorAmC[2][idx][idy]=ScaleFactorAmC[0][idx][idy];
                ScaleFactorAmC[3][idx][idy]=ScaleFactorAmC[0][idx][idy];
                ScaleFactorAmC[4][idx][idy]=ScaleFactorAmC[0][idx][idy];
                ScaleFactorAmC[5][idx][idy]=ScaleFactorAmC[0][idx][idy];
                
                if(NADs==MaxDetectors)//8AD
                {
                    ScaleFactorAmC[6][idx][idy]=ScaleFactorAmC[0][idx][idy];
                    ScaleFactorAmC[7][idx][idy]=ScaleFactorAmC[0][idx][idy];
                }
            }
        }
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
            
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    for(Int_t AD=0;AD<NADs;AD++)
                    {
                        //RandomLiHeH[AD][idx][idy] = new TH1D(Form("LiHe_AD%d_volume%d",AD,idx),Form("LiHe_AD%d_volume%d",AD,idx),n_evis_bins,evis_bins);
                        
                        for (Int_t visbin = 1; visbin <= n_evis_bins; visbin++)
                        {
                            RandomLiHeH[AD][idx][idy]->SetBinContent(visbin,(.9+DistortLiHeError)*LithiumH->Interpolate(RandomLiHeH[AD][idx][idy]->GetXaxis()->GetBinCenter(visbin))+(.1-DistortLiHeError)*HeliumH->Interpolate(RandomLiHeH[AD][idx][idy]->GetXaxis()->GetBinCenter(visbin)));
                        }
                    }
                    for(Int_t AD=0;AD<NADs;AD++)
                    {
                        RandomLiHeH[AD][idx][idy]->Scale(LiHeH[AD][idx][idy]->Integral()/RandomLiHeH[AD][idx][idy]->Integral());
                    }
                }
            }
            TCanvas* LiC = new TCanvas("a","");

                LithiumH->Draw();
            
            LiC->Print(("./Images/"+AnalysisString+"/BackgroundVariations/NominalLi.eps").c_str());
            delete LiC;
            
            TCanvas* HeC = new TCanvas("b","");

                HeliumH->Draw();
            
            HeC->Print(("./Images/"+AnalysisString+"/BackgroundVariations/NominalHe.eps").c_str());
            delete HeC;
            
            TCanvas* LiHeC = new TCanvas("c","");
            LiHeC->Divide(NADs/2,2*XCellLimit*YCellLimit);
            for(Int_t AD=0;AD<NADs;AD++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        LiHeC->cd(AD+idx*NADs+idy*NADs*XCellLimit+1);
                        RandomLiHeH[AD][idx][idy]->Draw();
                    }
                }
            }
            LiHeC->Print(("./Images/"+AnalysisString+"/BackgroundVariations/RandomLiHe.eps").c_str());
            delete LiHeC;

            delete LithiumH;
            delete HeliumH;
        }
        else
        {
        
#ifdef UseLiHeToyMC//If TOY MC generated histogram:
            
            m_file_distortLi9Bg = new TFile("./Inputs/GdInputs/8he9li_distort_neutron100_alpha100_frac0.1_N250.root");
            m_tree_distortLi9Bg = (TTree*)m_file_distortLi9Bg->Get("tr_distort");
            m_tree_distortLi9Bg->SetBranchAddress("h_distort",&LocalCopy);
            Int_t m_entries_distortLi9Bg = (Int_t)m_tree_distortLi9Bg->GetEntries();
            std::cout << "The distortion tree has " << m_entries_distortLi9Bg << " entries" << std::endl;
            Int_t entry = rand->Uniform(0,m_entries_distortLi9Bg);
            std::cout << "Reading Li He entry : " << entry << std::endl;
            m_tree_distortLi9Bg->GetEntry(entry);
            Hist_LiHe = (TH1F*)LocalCopy->Clone();
            delete m_file_distortLi9Bg;//Once the histogram is saved in a local copy we can close the file and the tree.

            TFile* OriginalLiHeF = TFile::Open("./BackgroundSpectrum/GDBackground/li9_spectrum.root");
            
            TH1F* OriginalLiHeH=(TH1F*)gDirectory->Get("h_li9_smeared_toy");
            
            OriginalLiHeF->Close();
            
            TH1F* Hist_LiHe_interpolated = (TH1F*)OriginalLiHeH->Clone();//To solve different bin limits error
            
            Hist_LiHe_interpolated->Reset();
            
            for(Int_t i = 1; i<=Hist_LiHe_interpolated->GetXaxis()->GetNbins();i++)
            {
                Hist_LiHe_interpolated->SetBinContent(i,Hist_LiHe->Interpolate(Hist_LiHe_interpolated->GetXaxis()->GetBinCenter(i)));
            }
            
            OriginalLiHeH->Multiply(Hist_LiHe_interpolated);

            for(Int_t AD=0;AD<NADs;AD++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        RandomLiHeH[AD][idx][idy] = (TH1D*)OriginalLiHeH->Rebin(n_evis_bins,Form("Rebinned LiHe"),evis_bins);//Rebin to visible binning
                        
                        RandomLiHeH[AD][idx][idy]->Scale(LiHeH[AD][idx][idy]->Integral()/RandomLiHeH[AD][idx][idy]->Integral());
                    }
                }
            }
            
            delete Hist_LiHe_interpolated;
            delete OriginalLiHeH;
#else//If distortion function
            
            func_LiHe = new TF1("func_LiHe","TMath::Abs([0]+[1]*x)",InitialVisibleEnergy,FinalVisibleEnergy);

            //function to obtaing distortion function or histogram:
            GetDistortionFunction(DistortLiHe,func_LiHe);
            
            for(Int_t AD=0;AD<NADs;AD++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        RandomLiHeH[AD][idx][idy]->Multiply(func_LiHe);
                        RandomLiHeH[AD][idx][idy]->Scale(LiHeH[AD][idx][idy]->Integral()/RandomLiHeH[AD][idx][idy]->Integral());
                    }
                }
            }
#endif
        }
        
        TFile* SaveLiHe = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/LiHeDistortions.root")).c_str(),"recreate");
        for(Int_t AD=0;AD<NADs;AD++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    LiHeH[AD][idx][idy]->Write(Form("Nominal LiHe%d period%d",AD,week));
                    RandomLiHeH[AD][idx][idy]->Write(Form("LiHeAD%d period%d",AD,week));
                }
            }
        }
#ifdef UseLiHeToyMC//If TOY MC generated histogram:
        if(!Analysis)
        {
            Hist_LiHe->Write();
        }
#else
        func_LiHe->Write();
        delete func_LiHe;
#endif

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
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        if(AD<ADsEH1)
                        {
                            RandomFastNeutronsH[AD][idx][idy]->Reset();
                            RandomFastNeutronsH[AD][idx][idy]->Add(func_FN1);
                            RandomFastNeutronsH[AD][idx][idy]->Scale(FastNeutronsH[AD][idx][idy]->Integral()/RandomFastNeutronsH[AD][idx][idy]->Integral());
                        }
                        else if(AD<(ADsEH1+ADsEH2))
                        {
                            RandomFastNeutronsH[AD][idx][idy]->Reset();
                            RandomFastNeutronsH[AD][idx][idy]->Add(func_FN2);
                            RandomFastNeutronsH[AD][idx][idy]->Scale(FastNeutronsH[AD][idx][idy]->Integral()/RandomFastNeutronsH[AD][idx][idy]->Integral());
                        }
                        else
                        {
                            RandomFastNeutronsH[AD][idx][idy]->Reset();
                            RandomFastNeutronsH[AD][idx][idy]->Add(func_FN3);
                            RandomFastNeutronsH[AD][idx][idy]->Scale(FastNeutronsH[AD][idx][idy]->Integral()/RandomFastNeutronsH[AD][idx][idy]->Integral());
                        }
                    }
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
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        RandomFastNeutronsH[AD][idx][idy]->Multiply(func_FN1);
                        RandomFastNeutronsH[AD][idx][idy]->Scale(FastNeutronsH[AD][idx][idy]->Integral()/RandomFastNeutronsH[AD][idx][idy]->Integral());
                    }
                }
            }
            func_FN2 = new TF1("func_FN2","[0]/(0.2*pow(x,0.1))+[1]",InitialVisibleEnergy,FinalVisibleEnergy);
            
            GetFastNeutronsDistortionFunction(DistortFN,func_FN2);
            
            for(Int_t AD=ADsEH1;AD<(ADsEH1+ADsEH2);AD++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        RandomFastNeutronsH[AD][idx][idy]->Multiply(func_FN2);
                        RandomFastNeutronsH[AD][idx][idy]->Scale(FastNeutronsH[AD][idx][idy]->Integral()/RandomFastNeutronsH[AD][idx][idy]->Integral());
                    }
                }
            }
            func_FN3 = new TF1("func_FN3","[0]/(0.2*pow(x,0.1))+[1]",InitialVisibleEnergy,FinalVisibleEnergy);
            
            GetFastNeutronsDistortionFunction(DistortFN,func_FN3);
            
            for(Int_t AD=(ADsEH1+ADsEH2);AD<(ADsEH1+ADsEH2+ADsEH3);AD++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        RandomFastNeutronsH[AD][idx][idy]->Multiply(func_FN3);
                        RandomFastNeutronsH[AD][idx][idy]->Scale(FastNeutronsH[AD][idx][idy]->Integral()/RandomFastNeutronsH[AD][idx][idy]->Integral());
                    }
                }
            }
        }
        
        TFile* SaveFN = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/FNDistortions.root")).c_str(),"recreate");
        for(Int_t AD = 0; AD<NADs;AD++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    FastNeutronsH[AD][idx][idy]->Write(Form("Nominal FN%d period%d",AD,week));
                    RandomFastNeutronsH[AD][idx][idy]->Write(Form("FN%d period%d",AD,week));
                }
            }
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
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    RandomAmCH[AD][idx][idy]->Reset();
                    RandomAmCH[AD][idx][idy]->Add(func_AmC);
                    RandomAmCH[AD][idx][idy]->Scale(AmCH[AD][idx][idy]->Integral()/RandomAmCH[AD][idx][idy]->Integral());
                }
            }
        }
        
        TFile* SaveAmC = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/AmCDistortions.root")).c_str(),"recreate");
        for(Int_t AD=0;AD<NADs;AD++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    AmCH[AD][idx][idy]->Write(Form("Nominal AmC%d period%d",AD,week));
                    RandomAmCH[AD][idx][idy]->Write(Form("AmC%d period%d",AD,week));
                }
            }
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
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    RandomAccidentalsH[AD][idx][idy]->Scale(1.*ScaleFactorAccidental[AD][idx][idy]);
                    AccidentalsH[AD][idx][idy]->Write(Form("Nominal Accidental%d period%d",AD,week));
                    RandomAccidentalsH[AD][idx][idy]->Write(Form("Accidental%d period%d",AD,week));
                }
            }
        }
        delete SaveVaryACCCanvasF;
    }
    if(VaryLiHeMatrix)
    {
        
        TFile* SaveVaryLiHeCanvasF = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/LiHeVariations.root")).c_str(),"recreate");
        for (Int_t AD =0; AD<NADs; AD++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    RandomLiHeH[AD][idx][idy]->Scale(1.*ScaleFactorLiHe[AD][idx][idy]);
                    LiHeH[AD][idx][idy]->Write(Form("Nominal LiHe%d period%d",AD,week));
                    RandomLiHeH[AD][idx][idy]->Write(Form("LiHe%d period%d",AD,week));
                }
            }
        }
        delete SaveVaryLiHeCanvasF;
    }
    if(VaryFastNeutronsMatrix)
    {
        
        TFile* SaveVaryFNCanvasF = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/FNVariations.root")).c_str(),"recreate");
        for (Int_t AD =0; AD<NADs; AD++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    RandomFastNeutronsH[AD][idx][idy]->Scale(1.*ScaleFactorFastNeutrons[AD][idx][idy]);
                    FastNeutronsH[AD][idx][idy]->Write(Form("Nominal FN%d period%d",AD,week));
                    RandomFastNeutronsH[AD][idx][idy]->Write(Form("FNAD%d period%d",AD,week));
                }
            }
        }
        delete SaveVaryFNCanvasF;
    }
    if(VaryAmCMatrix)
    {
        TFile* SaveVaryAmCCanvasF = TFile::Open(("./RootOutputs/"+AnalysisString+Form("/Backgrounds/AmCVariations.root")).c_str(),"recreate");
        for (Int_t AD =0; AD<NADs; AD++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    RandomAmCH[AD][idx][idy]->Scale(1.*ScaleFactorAmC[AD][idx][idy]);
                    AmCH[AD][idx][idy]->Write(Form("Nominal AmCH%d period%d",AD,week));
                    RandomAmCH[AD][idx][idy]->Write(Form("AmCAD%d period%d",AD,week));
                }
            }
        }
        delete SaveVaryAmCCanvasF;
    }
    
    std::cout <<  "********************************************************************************************************" << std::endl;
    
    for(Int_t AD=0; AD<NADs; AD++)
    {
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                if(VaryAccidentalMatrix||VaryLiHeMatrix||VaryFastNeutronsMatrix||VaryAmCMatrix)
                {
                    std::cout << "\n \t Background Scale Factors:" << std::endl;
                    std::cout << "\t Acc LiHe FN AmC" << std::endl;
                    std::cout << "AD: "<< AD << " Volume: " << idx << "\t " << ScaleFactorAccidental[AD][idx][idy] << "   " << ScaleFactorLiHe[AD][idx][idy] <<"   " <<  ScaleFactorFastNeutrons[AD][idx][idy] <<"   " <<  ScaleFactorAmC[AD][idx][idy] << std::endl;
                }
                
                //  Now BackgroundSpectrumH holds random backgrounds, either varying the rate or the shape
                BackgroundSpectrumH[AD][idx][idy]=(TH1D*)RandomAccidentalsH[AD][idx][idy]->Clone();
                BackgroundSpectrumH[AD][idx][idy]->Add(RandomLiHeH[AD][idx][idy]);
                BackgroundSpectrumH[AD][idx][idy]->Add(RandomFastNeutronsH[AD][idx][idy]);
                BackgroundSpectrumH[AD][idx][idy]->Add(RandomAmCH[AD][idx][idy]);
                BackgroundSpectrumH[AD][idx][idy]->SetTitle(Form("Random Background Spectrum_%d",AD));
                
        delete RandomAccidentalsH[AD][idx][idy];
        delete RandomLiHeH[AD][idx][idy];
        delete RandomFastNeutronsH[AD][idx][idy];
        delete RandomAmCH[AD][idx][idy];
            }
        }
    }
}

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
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                AccidentalsH[AD][idx][idy]=(TH1D*)gDirectory->Get(Form("Accidentals_AD%i_Volume%i",AD,idx));
                LiHeH[AD][idx][idy]=(TH1D*)gDirectory->Get(Form("LiHe_AD%i_Volume%i",AD,idx));
                FastNeutronsH[AD][idx][idy]=(TH1D*)gDirectory->Get(Form("FN_AD%i_Volume%i",AD,idx));
                AmCH[AD][idx][idy]=(TH1D*)gDirectory->Get(Form("AmC_AD%i_Volume%i",AD,idx));
            }
        }
        NominalAmCF = (TF1*)gDirectory->Get("AmCFunc");
        
        //Nominal background spectrum
        for(Int_t idx=0; idx<XCellLimit; idx++)
        {
            for(Int_t idy=0; idy<YCellLimit; idy++)
            {
                BackgroundSpectrumH[AD][idx][idy] = new TH1D(Form("Background Spectrum_AD%i_Volume%i",AD,idx),Form("Background Spectrum_AD%i_Volume%i",AD,idx),n_evis_bins,evis_bins);
                BackgroundSpectrumH[AD][idx][idy]->Add(AccidentalsH[AD][idx][idy]);
                BackgroundSpectrumH[AD][idx][idy]->Add(LiHeH[AD][idx][idy]);
                BackgroundSpectrumH[AD][idx][idy]->Add(FastNeutronsH[AD][idx][idy]);
                BackgroundSpectrumH[AD][idx][idy]->Add(AmCH[AD][idx][idy]);
            }
        }
    }
    
    #ifdef PrintEps
        TCanvas* AccidentalsC = new TCanvas("AccidentalsC","AccidentalsC");
        AccidentalsC->Divide(NADs/2,2*XCellLimit*YCellLimit);
        for(Int_t AD=0; AD<NADs; AD++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    AccidentalsC->cd(AD+idx*NADs+idy*NADs*XCellLimit+1);
                    AccidentalsH[AD][idx][idy]->Draw();
                }
            }
        }
        AccidentalsC->Update();
        
        AccidentalsC->Print(("./Images/"+AnalysisString+"/BackgroundVariations/NominalAccidentals.eps").c_str(),".eps");
        
        delete AccidentalsC;
        
        TCanvas* LiHeC = new TCanvas("LiHeC","LiHeC");
        LiHeC->Divide(NADs/2,2*XCellLimit*YCellLimit);
        for(Int_t AD=0; AD<NADs; AD++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    LiHeC->cd(AD+idx*NADs+idy*NADs*XCellLimit+1);
                    LiHeH[AD][idx][idy]->Draw();
                }
            }
        }
        LiHeC->Update();
        
        LiHeC->Print(("./Images/"+AnalysisString+"/BackgroundVariations/NominalLiHe.eps").c_str(),".eps");
        
        delete LiHeC;
        
        TCanvas* AmCC = new TCanvas("AmCC","AmCC");
        AmCC->Divide(NADs/2,2*XCellLimit*YCellLimit);
        for(Int_t AD=0; AD<NADs; AD++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    AmCC->cd(AD+idx*NADs+idy*NADs*XCellLimit+1);
                    AmCH[AD][idx][idy]->Draw();
                }
            }
        }
        AmCC->Update();
        
        AmCC->Print(("./Images/"+AnalysisString+"/BackgroundVariations/NominalAmC.eps").c_str(),".eps");
        
        delete AmCC;
        
        TCanvas* FastNeutronsC = new TCanvas("FastNeutronsC","FastNeutronsC");
        FastNeutronsC->Divide(NADs/2,2*XCellLimit*YCellLimit);
        for(Int_t AD=0; AD<NADs; AD++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    FastNeutronsC->cd(AD+idx*NADs+idy*NADs*XCellLimit+1);
                    FastNeutronsH[AD][idx][idy]->Draw();
                }
            }
        }
        FastNeutronsC->Update();
        
        FastNeutronsC->Print(("./Images/"+AnalysisString+"/BackgroundVariations/NominalFastNeutrons.eps").c_str(),".eps");
        
        delete FastNeutronsC;
        
    #endif
    delete BackgroundsF;
}
