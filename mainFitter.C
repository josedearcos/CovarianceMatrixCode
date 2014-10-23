#include "Fitter.h"
#include <stdio.h>
#include "TColor.h"
#include "TStyle.h"
#include "NominalData.h"

int main(void)
{
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);

    gROOT->SetBatch(kTRUE);

    const bool Fit2D = 0; // 1 for 2D Fit, 0 for 1D Fit
    const bool FitSin22t13 = 1; //  1 for Sin22t13 Fit, 0 for DM231 Fit.
    
    const bool AutomaticBudget = 0; //  To loop the fitter over all possible systematics
    const bool TurnOnBudget=0; // Generates Turn On Error Budget
    const bool TurnOffBudget=0; // Generates Turn Off Error Budget

    const bool RangeSin = 1;    //  Fit over a range of 20 Sin22t13
    const bool RangeDelta = 0;  //  Fit over a range of 20 Delta2M31

    const bool ToyMC=0;//1 ToyMC, 0 Data. To test the fitter use 1, to fit real data use 0.

    const bool BCW=0;
    const bool LBNL=0;
    const bool Unified=1;//Random NL not functional yet for this one.
    
    //To draw using a better palette:
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    std::cout <<  "********************************************************************************************************" << std::endl;
    std::cout <<  "**********************************           FITTER            *****************************************" << std::endl;
    std::cout <<  "********************************************************************************************************" << std::endl;
    
    const bool Analysis =0;//0 for Gd, 1 for H
    Int_t DataSet=2;
    NominalData* Data = new NominalData(Analysis,DataSet);
    Data->SetDataSet(DataSet);//0 is Simulation, 1 is P12E, 2 is LBNL DataSet

    Data->SetCombineMode(1); //0 is 9x9, 1 is 1x1 and 2 is 2x2
//    Data->SetSin22t12(0);
//    Data->SetSin22t13(0);
    //Parameters of the model
    Data->SetStatisticalFluctuation(0);
    Data->SetAnalysis(Analysis);//  Gd or H data
    Data->SetBinning(0);//  0 for LBNL binning or 1 for Linear binning
    Data->SetNSteps(101);// 101 in the final version.
    Data->SetWeeks(1); //   1 or 20 so far.

    Data->SetBCWModel(BCW);
    Data->SetLBNLModel(LBNL);
    Data->SetUnifiedModel(Unified);
    
//    // To fit we use the nominal set, no variations:
//    Data->SetIsotopeMatrix(1);
//    Data->SetReactorPowerMatrix(1);
//    Data->SetRelativeEnergyScaleMatrix(1);
//    Data->SetAbsoluteEnergyScaleMatrix(1);
//    Data->SetRelativeEnergyOffsetMatrix(1);
//    Data->SetAbsoluteEnergyOffsetMatrix(1);
//    Data->SetIAVMatrix(1);
//    Data->SetNLMatrix(1);
//    Data->SetResolutionMatrix(1);
//    Data->SetEfficiencyMatrix(1);
//    Data->SetSin22t12Matrix(1);

    Data->SetToyMC(ToyMC);
    
    //Chose Data Set
    if(Data->GetAnalysis())//   Hydrogen data
    {
        if(DataSet==1) // P12 production data
        {

            if(1==Data->GetWeeks())
            {
                std::cout << "\t Loading nH P12E Data" << std::endl;
                Data->LoadMainData("./Inputs/HInputs/P12E_Inclusive.txt");
            }
            else
            {
                std::cout << "\t NO MULTIPLE WEEK P12E DATA IN H ANALYSIS YET " << std::endl;
                Data->LoadMainData("./Inputs/HInputs/P12E_20.txt");
            }
        }
        else // Simple reactor model used as input data
        {
            std::cout << "\t Loading simple reactor model" << std::endl;
        }
    }
    else//  Gd data
    {
        if(DataSet==1)// P12E production data
        {

            if(1==Data->GetWeeks())
            {
                std::cout << "\t Loading Gd P12E Data" << std::endl;
//                Data->LoadMainData("./Inputs/GdInputs/P12E_Inclusive.txt");
            }
            else
            {
                std::cout << "\t NO MULTIPLE WEEK P12E DATA IN Gd ANALYSIS YET " << std::endl;
                Data->LoadMainData("./Inputs/GdInputs/P12E_20.txt");
            }
        }
        else if(DataSet==2)
        {
            if(1==Data->GetWeeks())
            {
                std::cout << "\t Loading LBNL Gd Data" << std::endl;
//                Data->LoadMainData("./Inputs/GdInputs/Theta13-inputs_20week_inclusive.txt");
            }
            else
            {
                std::cout << "\t Loading 20 weeks LBNL Gd Data" << std::endl;
                Data->LoadMainData("./Inputs/GdInputs/Theta13-inputs_20week.txt");
            }
        }
        else//  Simple reactor model used as input data
        {
            std::cout << "\t Loading simple reactor model" << std::endl;
        }
    }
    
    // Specific combinations to produce the error budget:

    if (TurnOnBudget||TurnOffBudget)//   Turns off the selected systematic
    {
        Data->SetTurnOffBudget(TurnOffBudget);
        Data->SetTurnOnBudget(TurnOnBudget);
        
        if(AutomaticBudget)
        {
            for (Int_t i = 0; i<=(15); i++) //18+1, the 19th is to produce the normal fit with all systematics in the TurnOff case.
            {
                std::cout << "\t Error Budget run #" << i << std::endl;
                
                Data->SetVaryAccidentalBudget(0==i);
                Data->SetVaryLiHeBudget(1==i);
                Data->SetVaryFastNeutronsBudget(2==i);
                Data->SetVaryAmCBudget(3==i);
                Data->SetDistortLiHeBudget(4==i);
                Data->SetDistortFastNeutronsBudget(5==i);
                Data->SetDistortAmCBudget(6==i);
                Data->SetIsotopeBudget(7==i);
                Data->SetReactorPowerBudget(8==i);
                Data->SetRelativeEnergyScaleBudget(9==i);
//                Data->SetRelativeEnergyOffsetBudget(10==i);
//                Data->SetAbsoluteEnergyScaleBudget(11==i);
//                Data->SetAbsoluteEnergyOffsetBudget(12==i);
                Data->SetIAVBudget(10==i);
                Data->SetNLBudget(11==i);
                Data->SetResolutionBudget(12==i);
                Data->SetSin22t12Budget(13==i);
                Data->SetSystematicBudget(14==i);
                Data->SetBackgroundBudget(15==i);
                
                Fitter* Fit = new Fitter(Data);
                Fit->MainFitter(Fit2D,FitSin22t13);
                
                if (Fit2D)
                {
                    Fit->Save2DFit(i);
                }
                else
                {
                    if (FitSin22t13)
                    {
                        Fit->SaveSin1DFit(i,TurnOnBudget,TurnOffBudget);
                    }
                    else
                    {
                        Fit->SaveDM1DFit(i,TurnOnBudget,TurnOffBudget);
                    }
                }
                
                delete Fit;
            }
        }
        else
        {
            Data->SetVaryAccidentalBudget(0);
            Data->SetVaryLiHeBudget(0);
            Data->SetVaryFastNeutronsBudget(1);
            Data->SetVaryAmCBudget(0);
            Data->SetDistortLiHeBudget(0);
            Data->SetDistortFastNeutronsBudget(0);
            Data->SetDistortAmCBudget(0);
            Data->SetIsotopeBudget(0);
            Data->SetReactorPowerBudget(0);
            Data->SetRelativeEnergyScaleBudget(0);
//            Data->SetRelativeEnergyOffsetBudget(0);
//            Data->SetAbsoluteEnergyScaleBudget(0);
//            Data->SetAbsoluteEnergyOffsetBudget(0);
            Data->SetIAVBudget(0);
            Data->SetNLBudget(0);
            Data->SetResolutionBudget(0);
            Data->SetSin22t12Budget(0);
            Data->SetSystematicBudget(0);
            Data->SetBackgroundBudget(0);
            
            Fitter* Fit = new Fitter(Data);
            Fit->MainFitter(Fit2D,FitSin22t13);
            
            if (Fit2D)
            {
                Fit->Save2DFit(0);
            }
            else
            {
                if (FitSin22t13)
                {
                    Fit->SaveSin1DFit(0,TurnOnBudget,TurnOffBudget);
                }
                else
                {
                    Fit->SaveDM1DFit(0,TurnOnBudget,TurnOffBudget);
                }
            }
            
            delete Fit;
        }
    }
    else
    {
        std::cout << "\t Fit run" << std::endl;

        if(RangeSin)
        {
            std::cout << "\t Sin range run" << std::endl;

            Double_t Sin22t13[21];
            
            for (Int_t i = 0; i<21; i++)
            {
                Sin22t13[i] = i*.2/(20);
            }

            for (Int_t i = 0; i<21; i++)
            {
                Data->SetSin22t13(Sin22t13[i]);
                
                Fitter* Fit = new Fitter(Data);
                Fit->MainFitter(Fit2D,FitSin22t13);
                Fit->SaveSinRangeChiSquare(i,FitSin22t13);
                delete Fit;
            }
        }
        else if(RangeDelta)
        {
            std::cout << "\t Delta range run" << std::endl;

            Double_t Dm2ee[21];

            for (Int_t i = 0; i<21; i++)
            {
                Dm2ee[i] = 0.0015+(i*.002/20);
            }
            
            for (Int_t i = 0; i<21; i++)
            {
                Data->SetDm2ee(Dm2ee[i]);
                
                Fitter* Fit = new Fitter(Data);
                Fit->MainFitter(Fit2D,FitSin22t13);
                Fit->SaveDeltaMRangeChiSquare(i,FitSin22t13);
                delete Fit;
            }
        }
        else
        {
            Fitter* Fit = new Fitter(Data);
            Fit->MainFitter(Fit2D,FitSin22t13);
            
            if (Fit2D)
            {
                Fit->Save2DFit(0);
            }
            else
            {
                if (FitSin22t13)
                {
                    Fit->SaveSin1DFit(0,TurnOnBudget,TurnOffBudget);
                }
                else
                {
                    Fit->SaveDM1DFit(0,TurnOnBudget,TurnOffBudget);
                }
            }
            delete Fit;
        }
    }
    
    gROOT->SetBatch(kFALSE);

    return 0;
}