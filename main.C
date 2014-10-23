#include <stdio.h>
#include "FitBackgrounds2.h"
#include "CovarianceMatrix3.h"
#include "NominalData.h"
#include "TColor.h"
#include "TStyle.h"
#include "CrossSection.h"
#include <stdlib.h>  

int main(void)
{    
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);

    gROOT->SetBatch(kTRUE);//This avoids canvas popping up every time.
    
    const Int_t CovarianceMatrices=15;
    
    //To draw using a better palette:
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;    
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    
    const bool Analysis =0;//0 for Gd, 1 for H
    Int_t DataSet=2;//0 is Simulation, 1 is P12E, 2 is LBNL DataSet

    NominalData* Data = new NominalData(Analysis,DataSet);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // CrossCalc takes a lot of time and it needs to be run only once to produce the root file. Uncomment only if binning is different from standards and you want to produce a new cross section for this new binning.
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    CrossSection* Cross = new CrossSection(Data);
    Cross->CrossCalc();
    delete Cross;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    const bool ToyMC=1;//  1 for Toy MC, 0 for data. To produce covariance matrices we use ToyMC.
    //Parameters of the model
    
    Data->SetNSamples(500);// 500 in the final version.

    Data->SetToyMC(ToyMC);

    Data->SetCombineMode(0); //0 is 9x9, 1 is 1x1 and 2 is 2x2
    
    Data->SetBinning(0);//  0 for LBNL binning or 1 for Linear binning
    Data->SetWeeks(1); //   1 or 20 so far.

    std::cout <<  "********************************************************************************************************" << std::endl;
    std::cout <<  "**********************************            MAIN             *****************************************" << std::endl;
    std::cout <<  "********************************************************************************************************" << std::endl;

    //Chose Data Set
    if(Data->GetAnalysis())//   Hydrogen data
    {
        switch (DataSet)
        {
            case 0://  Simple reactor model used as input data
                std::cout << "\t Loading simple reactor model" << std::endl;
                break;
            case 1://   P12E
                if(1==Data->GetWeeks())
                {
                    std::cout << "\t Loading nH P12E Data" << std::endl;
                    Data->LoadMainData("./Inputs/HInputs/P12E_Inclusive.txt");
                }
                else
                {
                    std::cout << "\t \t \t NO MULTIPLE WEEK P12E DATA IN H ANALYSIS YET " << std::endl;
                    exit(EXIT_FAILURE);
                    Data->LoadMainData("./Inputs/HInputs/P12E_20.txt");
                }
                break;
            case 2:// LBNL
                std::cout << "\t \t \t NO LBNL H ANALYSIS, OPTION NOT VALID " << std::endl;
                exit(EXIT_FAILURE);
                break;
            default:
                break;
        }
    }
    else//  Gd data
    {
        switch (DataSet)
        {
            case 0://  Simple reactor model used as input data
                std::cout << "\t Loading simple reactor model" << std::endl;
                break;
            case 1://   P12E
                if(1==Data->GetWeeks())
                {
                    std::cout << "\t Loading Gd P12E Data" << std::endl;
//                    Data->LoadMainData("./Inputs/GdInputs/P12E_Inclusive.txt");
                }
                else
                {
                    std::cout << "\t \t \t NO MULTIPLE WEEK P12E DATA IN Gd ANALYSIS YET " << std::endl;
                    exit(EXIT_FAILURE);
                    Data->LoadMainData("./Inputs/GdInputs/P12E_20.txt");
                }
                break;
            case 2:// LBNL
                if(1==Data->GetWeeks())
                {
                    std::cout << "\t Loading LBNL Gd Data" << std::endl;
//                    Data->LoadMainData("./Inputs/GdInputs/Theta13-inputs_20week_inclusive.txt");
                }
                else
                {
                    std::cout << "\t Loading 20 weeks LBNL Gd Data" << std::endl;
                    Data->LoadMainData("./Inputs/GdInputs/Theta13-inputs_20week.txt");
                }
                break;
            default:
                break;
        }
    }

    //To choose the NL model set to 1. Activate only one at a time.
    //NL
    bool BCW=0;
    bool LBNL=0;
    bool Unified=1;//Random NL with correlations not functional yet for this one.

    //To check non theta13 oscillation behaviour of the fitter. To avoid any oscillation at all then theta12 has to be set to 0 too.
//    Data->SetSin22t12(0);
//    Data->SetSin22t13(0);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    FitBackgrounds2* Bkg = new FitBackgrounds2(Data);

    if(Analysis)
    {
        Bkg->ReadHBackgrounds();
    }
    else
    {
        Bkg->ReadGdBackgrounds();
    }
    
    delete Bkg;
    
    CovarianceMatrix3* Cov;
    
    Int_t Automatic = 1;
    if (Automatic==1)//1 for producing all matrices, others to choose specific combinations manually
    {//Automatic Script:
        std::cout <<  "********************************************************************************************************" << std::endl;
        std::cout << "\t Automatic Script" << std::endl;

        Data->SetBCWModel(BCW);
        Data->SetLBNLModel(LBNL);
        Data->SetUnifiedModel(Unified);
        //Order is important (Inititialize NL model before randomize it)

        for (Int_t i = 0; i<(CovarianceMatrices); i++)
        {
            Data->SetVaryAccidentalMatrix(0==i);
            Data->SetVaryLiHeMatrix(1==i);
            Data->SetVaryFastNeutronsMatrix(2==i);
            Data->SetVaryAmCMatrix(3==i);
            Data->SetDistortLiHeMatrix(4==i);
            Data->SetDistortFastNeutronsMatrix(5==i);
            Data->SetDistortAmCMatrix(6==i);
            Data->SetIsotopeMatrix(7==i);
            Data->SetReactorPowerMatrix(8==i);
            Data->SetRelativeEnergyScaleMatrix(9==i);
//            Data->SetAbsoluteEnergyScaleMatrix(10==i);
//            Data->SetRelativeEnergyOffsetMatrix(11==i);
//            Data->SetAbsoluteEnergyOffsetMatrix(12==i);
            Data->SetIAVMatrix(10==i);
            Data->SetNLMatrix(11==i);
            Data->SetResolutionMatrix(12==i);
            Data->SetSin22t12Matrix(13==i);
            Data->SetEfficiencyMatrix(14==i);
            
            Cov = new CovarianceMatrix3(Data);
                Cov->CovarianceMatrixMain();
            delete Cov;
        }
    }
    else if(Automatic==0)
    {
        std::cout <<  "********************************************************************************************************" << std::endl;
        std::cout << "\t Manual Script" << std::endl;

        Data->SetBCWModel(BCW);
        Data->SetLBNLModel(LBNL);
        Data->SetUnifiedModel(Unified);
        
        //Manual Control:
        //To calculate Covariance Matrices set to 1. Activate only one at a time.
        //Backgrounds
        bool VaryAccidentalMatrix=0;
        bool VaryLiHeMatrix=0;
        bool VaryFastNeutronsMatrix=0;
        bool VaryAmCMatrix=0;
        bool DistortLiHeMatrix=0;
        bool DistortFastNeutronsMatrix=0;
        bool DistortAmCMatrix=0;
        //Systematics
        bool IsotopeMatrix=0;
        bool ReactorPowerMatrix=0;
        bool EnergyScaleMatrix=0;
//        bool EnergyOffsetMatrix=0;
//        bool AbsoluteScaleMatrix=0;
//        bool AbsoluteOffsetMatrix=0;
        bool IAVMatrix=0;
        bool NLMatrix=0;
        bool ResolutionMatrix=0;
        bool Sin22t12Matrix=0;
        bool EfficiencyMatrix=0;
        
        Data->SetVaryAccidentalMatrix(VaryAccidentalMatrix);
        Data->SetVaryLiHeMatrix(VaryLiHeMatrix);
        Data->SetVaryFastNeutronsMatrix(VaryFastNeutronsMatrix);
        Data->SetVaryAmCMatrix(VaryAmCMatrix);
        Data->SetDistortLiHeMatrix(DistortLiHeMatrix);
        Data->SetDistortFastNeutronsMatrix(DistortFastNeutronsMatrix);
        Data->SetDistortAmCMatrix(DistortAmCMatrix);
        Data->SetIsotopeMatrix(IsotopeMatrix);
        Data->SetReactorPowerMatrix(ReactorPowerMatrix);
        Data->SetRelativeEnergyScaleMatrix(EnergyScaleMatrix);
//        Data->SetAbsoluteEnergyScaleMatrix(AbsoluteScaleMatrix);
//        Data->SetRelativeEnergyOffsetMatrix(EnergyOffsetMatrix);
//        Data->SetAbsoluteEnergyOffsetMatrix(AbsoluteOffsetMatrix);
        Data->SetIAVMatrix(IAVMatrix);
        Data->SetNLMatrix(NLMatrix);
        Data->SetResolutionMatrix(ResolutionMatrix);
        Data->SetSin22t12Matrix(Sin22t12Matrix);
        Data->SetEfficiencyMatrix(EfficiencyMatrix);

        Cov = new CovarianceMatrix3(Data);

        Cov->CovarianceMatrixMain();
        
        delete Cov;
    }

    delete Data;
    
    gROOT->SetBatch(kFALSE);

    return 0;
}
