#pragma once
#include "TH2D.h"
#include "TFile.h"
#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2.h"
#include "TMath.h"
#include <math.h>
#include <TArrayD.h>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include "TCanvas.h"
#include "NominalData.h"
#include "CrossSection.h"
#include "CovarianceMatrix3.h"
#include "Prediction.h"
#include "TMinuit.h"
#include "TDecompLU.h"

const bool FullGrid = 0;//To save all histograms in the 2D grid, in the case of 101x101 it can be immense!.
const Int_t MaxSteps = 501;
const Double_t chi2limit = 2.3;// 1 sigma for 2D

Prediction* Pred;// Not sure how to initialize it with Data, need to check this !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//  Function for minuit minimization
void minuit_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *x, Int_t iflag)
{
    std::cout << " FOR MINUIT PROPERLY TO WORK THIS SHOULD BE RUN SECOND" << std::endl;
    
    Double_t sin22t13 = x[0];
    Double_t dm2_ee = x[1];
    Int_t week = x[2];
    bool ToyMC = x[3];
    
    f = Pred->CalculateChi2(sin22t13, dm2_ee + 7.5e-5, week,ToyMC);
    
    std::cout << "\t Calculate χ2 FCN result: " << f << std::endl;
}

class Fitter
{
private:
    NominalData* Nom;
    
    Double_t s22t13start;
    Double_t s22t13end;
    Double_t dm2_eestart;
    Double_t dm2_eeend;
    
    bool ToyMC;
    bool Mode;
    
    Int_t Nsteps;
    Int_t NstepsM;
    Int_t Combine;
    Int_t Nweeks;
    Int_t hall;
    
    Double_t chi2_min[MaxPeriods];
    Double_t s2t_min[MaxPeriods];
    Double_t s2t_error_min[MaxPeriods];
    Double_t dm2_ee_error_min[MaxPeriods];
    Double_t dm2_ee_min[MaxPeriods];
    Int_t degenerancy;
    Double_t chi2;
    
    bool VaryAccidentalBudget;
    bool VaryLiHeBudget;
    bool VaryFastNeutronsBudget;
    bool VaryAmCBudget;
    bool DistortLiHeBudget;
    bool DistortFastNeutronsBudget;
    bool DistortAmCBudget;
    bool IsotopeBudget;
    bool ReactorPowerBudget;
    bool IAVBudget;
    bool NLBudget;
    bool ResolutionBudget;
    bool RelativeEnergyScaleBudget;
    bool RelativeEnergyOffsetBudget;
    bool AbsoluteEnergyScaleBudget;
    bool AbsoluteEnergyOffsetBudget;
    bool Sin22t12Budget;
    bool EfficiencyBudget;
    bool SystematicBudget;
    bool TotalBudget;
    bool BackgroundBudget;
    
    bool StatisticalFluctuation;
    
    //    Double_t chi2;
    Double_t Sin22t13;
    Double_t Dm2_ee;
    Double_t sin22t13[MaxSteps];
    Double_t Realsin22t13[MaxSteps];
    Double_t dm2_ee[MaxSteps];
    Double_t Realdm2_ee[MaxSteps];

    
    Double_t DeltaWidth;
    Double_t SinWidth;
    Double_t RealSinWidth;
    Double_t RealDeltaWidth;
    
    TH2D* ChiSquare2DH[MaxPeriods];
    TH1D* SinChiSquareH[MaxPeriods][MaxSteps];
    TH1D* DeltaChiSquareH[MaxPeriods][MaxSteps];
    TH1D* ASinChiSquareH[MaxPeriods];
    TH1D* ADeltaChiSquareH[MaxPeriods];
    TH1D* TestSinChiSquareH[MaxPeriods][MaxSteps];
    TH1D* TestDeltaChiSquareH[MaxPeriods][MaxSteps];
    
    Char_t FileName[100];
    bool Analysis;
    string AnalysisString;

    void SaveChiSquare(Int_t);
  public:
    Fitter();
    Fitter(NominalData*);
    void MainFitter(bool,bool,bool,FitterGui*);
    void TestRandomExperimentsChi2(bool,bool,bool,FitterGui*,Int_t,Int_t);
    void SaveSinRangeChiSquare(Int_t,bool);
    void SaveDeltaMRangeChiSquare(Int_t,bool);
    void SaveSin1DFit(Int_t,bool,bool);
    void SaveDM1DFit(Int_t,bool,bool);
    void Save2DFit(Int_t);
    void SaveChiTest(bool,bool,Int_t);
};

Fitter :: Fitter()
{
    std::cout << " the fitter default constructor shouldn't be called" << std::endl;
    
    exit(EXIT_FAILURE);
    
    Nom = new NominalData(0,2);
    
    Pred = new Prediction(Nom);
    Combine = Nom->GetCombineMode();
    Nweeks = Nom->GetWeeks();
    Nsteps = Nom->GetNSteps();
    NstepsM = Nsteps;
    
    Analysis = Nom->GetAnalysis();

    s22t13start= Nom->GetSinStart();
    s22t13end=Nom->GetSinEnd();
    dm2_eestart=Nom->GetDmeeStart();
    dm2_eeend=Nom->GetDmeeEnd();
    
    StatisticalFluctuation = Nom->GetStatisticalFluctuation();
    
    if(Analysis)
    {
        AnalysisString = "Hydrogen";
    }
    else
    {
        AnalysisString = "Gadolinium";
    }
    
    Sin22t13 = Nom->GetSin22t13();
    Dm2_ee = Nom->GetDm2ee();
    
    VaryAccidentalBudget = Nom->GetVaryAccidentalBudget();
    VaryLiHeBudget = Nom->GetVaryLiHeBudget();
    VaryFastNeutronsBudget = Nom->GetVaryFastNeutronsBudget();
    VaryAmCBudget = Nom->GetVaryAmCBudget();
    DistortLiHeBudget = Nom->GetDistortLiHeBudget();
    DistortFastNeutronsBudget = Nom->GetDistortFastNeutronsBudget();
    DistortAmCBudget = Nom->GetDistortAmCBudget();
    
    IsotopeBudget = Nom->GetIsotopeBudget();
    ReactorPowerBudget = Nom->GetReactorPowerBudget();
    RelativeEnergyScaleBudget = Nom->GetRelativeEnergyScaleBudget();
    RelativeEnergyOffsetBudget = Nom->GetRelativeEnergyOffsetBudget();
    AbsoluteEnergyScaleBudget = Nom->GetAbsoluteEnergyScaleBudget();
    AbsoluteEnergyOffsetBudget = Nom->GetAbsoluteEnergyOffsetBudget();
    IAVBudget = Nom->GetIAVBudget();
    NLBudget = Nom->GetNLBudget();
    ResolutionBudget = Nom->GetResolutionBudget();
    Sin22t12Budget = Nom->GetSin22t12Budget();
    EfficiencyBudget = Nom->GetEfficiencyBudget();
    SystematicBudget = Nom->GetSystematicBudget();
    BackgroundBudget = Nom->GetBackgroundBudget();
    TotalBudget = Nom->GetTotalBudget();
    ToyMC = Nom->GetToyMC();
    delete Nom;
}

Fitter :: Fitter(NominalData* Data)
{
    std::cout << " FOR MINUIT PROPERLY TO WORK THIS SHOULD BE RUN FIRST" << std::endl;
    Pred = new Prediction(Data);
    
    Combine = Data->GetCombineMode();
    Nweeks = Data->GetWeeks();
    Nsteps = Data->GetNSteps();
    NstepsM = Nsteps;
    
    Sin22t13 = Data->GetSin22t13();
    Dm2_ee = Data->GetDm2ee();
    
    Analysis = Data->GetAnalysis();
    
    s22t13start= Data->GetSinStart();
    s22t13end=Data->GetSinEnd();
    dm2_eestart=Data->GetDmeeStart();
    dm2_eeend=Data->GetDmeeEnd();
    
    StatisticalFluctuation = Data->GetStatisticalFluctuation();

    if(Analysis)
    {
        AnalysisString = "Hydrogen";
    }
    else
    {
        AnalysisString = "Gadolinium";
    }
    
    VaryAccidentalBudget = Data->GetVaryAccidentalBudget();
    VaryLiHeBudget = Data->GetVaryLiHeBudget();
    VaryFastNeutronsBudget = Data->GetVaryFastNeutronsBudget();
    VaryAmCBudget = Data->GetVaryAmCBudget();
    DistortLiHeBudget = Data->GetDistortLiHeBudget();
    DistortFastNeutronsBudget = Data->GetDistortFastNeutronsBudget();
    DistortAmCBudget = Data->GetDistortAmCBudget();
    IsotopeBudget = Data->GetIsotopeBudget();
    ReactorPowerBudget = Data->GetReactorPowerBudget();
    RelativeEnergyScaleBudget = Data->GetRelativeEnergyScaleBudget();
    RelativeEnergyOffsetBudget = Data->GetRelativeEnergyOffsetBudget();
    AbsoluteEnergyScaleBudget = Data->GetAbsoluteEnergyScaleBudget();
    AbsoluteEnergyOffsetBudget = Data->GetAbsoluteEnergyOffsetBudget();
    IAVBudget = Data->GetIAVBudget();
    NLBudget = Data->GetNLBudget();
    ResolutionBudget = Data->GetResolutionBudget();
    Sin22t12Budget = Data->GetSin22t12Budget();
    EfficiencyBudget = Data->GetEfficiencyBudget();
    SystematicBudget = Data->GetSystematicBudget();
    BackgroundBudget = Data->GetBackgroundBudget();
    TotalBudget = Data->GetTotalBudget();
    
    ToyMC = Data->GetToyMC();
}

void Fitter :: MainFitter(bool Minuit, bool Fit2D, bool FitSin22t13,FitterGui* FitterGui)
{
    std::cout << "\t Scanning parameters" << std::endl;
    
    for (Int_t week = 0; week < Nweeks; week++)
    {
        chi2_min[week] = 1e10;
        s2t_min[week] = 1e10;
        s2t_error_min[week] = 1e10;
        dm2_ee_min[week]= 1e10;
    }
    
    degenerancy = 0;
    Mode=0;//No fluctuation in data by default
    
    DeltaWidth =(dm2_eeend-dm2_eestart)*1./(NstepsM-1);
    SinWidth = (s22t13end-s22t13start)*1./(Nsteps-1);
    
    std::cout << " " << std::endl;
    std::cout << "\t \t \t STEP " << (s22t13end-s22t13start)*1./(Nsteps-1) << std::endl;
    std::cout << "\t \t \t STEPM " <<  (dm2_eeend-dm2_eestart)*1./(NstepsM-1) << std::endl;
    std::cout << " " << std::endl;
    
    if(Minuit)
    {
        const Int_t NPars = 2;//    2   when delta m is also included
        TMinuit* minu = new TMinuit(NPars);
        
        Int_t ierflag;
        Double_t arglist[4];
        arglist[0] = 1.0;
        
        for (Int_t week = 0; week < Nweeks; week++)
        {
            // Fit by Minuit
            Pred->LoadData(week,ToyMC,Nsteps,Mode);
            
            if(ReadTxt)
            {
                Pred->LoadTxtCovarianceMatrices(week);
            }
            else
            {
                Pred->LoadRootCovarianceMatrices(week);
            }
            Pred->GenerateInverseMatrix(Sin22t13, Dm2_ee, week, 0, ToyMC,Nsteps);

            minu->SetFCN(minuit_fcn);
            
            minu->mnexcm("SET ERR",arglist,1,ierflag);
            minu->mnparm(0, "SinSq2Theta13", Sin22t13, (s22t13end-s22t13start)*1./(Nsteps-1), s22t13start,s22t13end,ierflag);
            minu->mnparm(1, "DeltaMSqee",  Dm2_ee, (dm2_eeend-dm2_eestart)*1./(NstepsM-1), dm2_eestart,dm2_eeend,ierflag);
            minu->mnparm(2, "Week", week, 0, 0,MaxPeriods, ierflag);
            minu->mnparm(3, "ToyMC", ToyMC, 0, 0,1, ierflag);
            
            minu->FixParameter(2);
            minu->FixParameter(3);
            
            arglist[0] = 10000;// #of iterations
            arglist[1] = 1.0;// 1 ChiSquare, 0.5 log likelihood
            minu->mnexcm("MIGRAD", arglist ,2,ierflag);//   Executes Migrad command
            
            minu->Release(2);
            minu->Release(3);

            //  Probably equivalent:
            //            minu->SetMaxIterations(10000);// #of iterations
            //            minu->SetErrorDef(1);// 1 ChiSquare, 0.5 log likelihood
            //            minu->Migrad();
            
            // minu->FixParameter(0);
            // minu->FixParameter(1);
            
            std::cout << "================ MIGRAD finished with error flag = " << ierflag  << " ============================" << std::endl;
            
            Double_t fpar,ferr;
            
            minu->GetParameter(0,fpar,ferr);
            s2t_min[week] = fpar;
            
            minu->GetParameter(1,fpar,ferr);
            dm2_ee_min[week] = fpar;
            
            Double_t min_pars[2] = {s2t_min[week],dm2_ee_min[week]};
            Double_t * grad;
            
            minu->Eval(2,grad,chi2_min[week],min_pars,0);
            
            std::cout << "======== fit results (for checking) :" <<  s2t_min[week] << " " << dm2_ee_min[week] << " " << chi2_min[week] << std::endl;
            
            // minu->SetErrorDef(n*n); //Call before Countour() To plot at n sigma
            TGraph* ContourMap = (TGraph*)minu->Contour(10,0,1);
            TCanvas* ContourMapCanvas = new TCanvas("ContourMap","CountourMap");
            ContourMap->Draw("ACP");
            ContourMapCanvas->Print(("./Images/"+AnalysisString+"/MinuitContourMap.eps").c_str(),".eps");
            delete ContourMapCanvas;
        }
        Pred->DeleteData();
        Pred->DeleteMatrices();
        delete minu;
    }
    else
    {
        for (Int_t week = 0; week < Nweeks; week++)
        {
            chi2=0;
            
            Pred->LoadData(week,ToyMC,Nsteps,Mode);
            
            if(ReadTxt)
            {
                Pred->LoadTxtCovarianceMatrices(week);
            }
            else
            {
                Pred->LoadRootCovarianceMatrices(week);
            }
            
            Pred->GenerateInverseMatrix(Sin22t13, Dm2_ee, week, 0, ToyMC,Nsteps);

            if(Fit2D)
            {
                //  Fit 101 x 101 grid in Sin22t13 & DM2_ee phase-space (10,201 Fits!)
                //  Computational demanding process, reduce the number of samples, my computer holds up to 20x20 (400 Fits)
                
                ChiSquare2DH[week] = new TH2D("2DChiDistribution","#chi^{2} Distribution",Nsteps,s22t13start,s22t13end,NstepsM,dm2_eestart,dm2_eeend);
                
                for(Int_t step=0;step < Nsteps;++step)
                {
                    sin22t13[step] = SinWidth*step + s22t13start;
                    
                    std::cout << "\t sin22t13 " << sin22t13[step] << std::endl;
                    
                    for(Int_t stepM = 0; stepM < NstepsM ; stepM++)
                    {
                        dm2_ee[stepM] = DeltaWidth*stepM + dm2_eestart;
                        
                        std::cout << "\t dm2_ee " << dm2_ee[stepM] << std::endl;
                        
                        chi2 = Pred->CalculateChi2(sin22t13[step],dm2_ee[stepM],week,ToyMC);
                        
                        ChiSquare2DH[week]->SetBinContent(step+1,stepM+1,TMath::Abs(chi2));
                        
                        if (TMath::Abs(chi2)<chi2_min[week])
                        {
                            s2t_min[week]=sin22t13[step];
                            dm2_ee_min[week]=dm2_ee[stepM];
                            chi2_min[week]=TMath::Abs(chi2);
                        }
                        else if(TMath::Abs(chi2)==chi2_min[week])
                        {
                            degenerancy++;
                        }
                        
                        if((chi2-chi2_min[week])<=chi2limit)//1sigma cut
                        {
                            s2t_error_min[week] = sin22t13[step];
                            dm2_ee_error_min[week] = dm2_ee[step];
                        }
                        
                        FitterGui->UpdateFitter(step*NstepsM+stepM);
                    }
                }
                
                if(FullGrid)//  Produces 202 (2*MaxSamples) cut histograms
                {
                    for(Int_t stepM=1;stepM <=NstepsM;++stepM)
                    {
                        SinChiSquareH[week][stepM]=(TH1D*)ChiSquare2DH[week]->ProjectionX(Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,dm2_ee[stepM],week),stepM+1,stepM+1);
                    }
                    for(Int_t step = 1; step <= Nsteps ; step++)
                    {
                        DeltaChiSquareH[week][step]=(TH1D*)ChiSquare2DH[week]->ProjectionY(Form("Sin%f_Distribution_DeltaM%f_Period%d",sin22t13[step],Dm2_ee,week),step+1,step+1);
                    }
                }
            }
            else if (FitSin22t13)
            {
                //  Fit 101 different Sin22t13 with a fix Dm2_ee
                ASinChiSquareH[week] = new TH1D(Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,Dm2_ee,week),Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,Dm2_ee,week),Nsteps,s22t13start,s22t13end);
                
                dm2_ee_min[week] = Dm2_ee;
                
                for(Int_t step=0;step < Nsteps;++step)
                {
                    sin22t13[step] = SinWidth*step + s22t13start;
                    
                    chi2 = Pred->CalculateChi2(sin22t13[step],Dm2_ee,week,ToyMC);
                    
                    if (TMath::Abs(chi2)<chi2_min[week])
                    {
                        s2t_min[week]=sin22t13[step];
                        chi2_min[week]=TMath::Abs(chi2);
                    }
                    else if(TMath::Abs(chi2)==chi2_min[week])
                    {
                        degenerancy++;
                    }
                    
                    ASinChiSquareH[week]->SetBinContent(step+1,TMath::Abs(chi2));
                    FitterGui->UpdateFitter(step);
                }
            }
            else
            {
                //  Fit 101 different DM2_ee with a fix Sin22t13
                
                ADeltaChiSquareH[week] = new TH1D(Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,Dm2_ee,week),Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,Dm2_ee,week),NstepsM,dm2_eestart,dm2_eeend);
                
                s2t_min[week] = Sin22t13;
                
                for(Int_t stepM = 0; stepM < NstepsM ; stepM++)
                {
                    dm2_ee[stepM] = DeltaWidth*stepM + dm2_eestart;
                    
                    chi2 = Pred->CalculateChi2(Sin22t13,dm2_ee[stepM],week,ToyMC);
                    
                    if (TMath::Abs(chi2)<chi2_min[week])
                    {
                        dm2_ee_min[week]=dm2_ee[stepM];
                        chi2_min[week]=TMath::Abs(chi2);
                    }
                    else if(TMath::Abs(chi2)==chi2_min[week])
                    {
                        degenerancy++;
                    }
                    
                    ADeltaChiSquareH[week]->SetBinContent(stepM+1,TMath::Abs(chi2));
                    FitterGui->UpdateFitter(stepM);
                }
            }
            
            std::cout << "\t \t \t For Period: " << week << std::endl;
            std::cout << "\t \t \t \t χ2 min" << chi2_min[week] << std::endl;
            std::cout << "\t \t \t \t Sin22t13 min" << s2t_min[week] << std::endl;
            //            std::cout << "\t \t \t \t Sin22t13 error" << s2t_error_min[week] << std::endl;
            std::cout << "\t \t \t \t Δm2_ee_min min" << dm2_ee_min[week] << std::endl;
            //            std::cout << "\t \t \t \t Δm2_ee error" << dm2_ee_error_min[week] << std::endl;
            std::cout << "\t \t \t \t χ2 degenerated " << degenerancy << " times" << std::endl;
            
            Pred->DeleteData();
            Pred->DeleteMatrices();
        }
    }
    
    delete Pred;
    
    std::cout << "\t Finished Scanning parameters" << std::endl;
}

void Fitter :: SaveSinRangeChiSquare(Int_t sample,bool FitSin22t13)
{
    std::cout << "\t SAVING Sin22t13 Range Fit χ2 HISTOGRAMS" << std::endl;
    
    TString optionS("recreate");
    
    if(sample>0)
    {
        optionS.ReplaceAll("recreate","update");
    }
    
    TString fileS(("./ChiSquare/"+AnalysisString+Form("/Combine%d/Sin22t13Sin22t13Range.root",Combine)).c_str());
    
    if(!FitSin22t13)
    {
        fileS.ReplaceAll(("./ChiSquare/"+AnalysisString+Form("/Combine%d/Sin22t13Sin22t13Range.root",Combine)).c_str(),("./ChiSquare/"+AnalysisString+Form("/Combine%d/DeltaMSin22t13Range.root",Combine)).c_str());
    }
    
    TFile* SaveSR = TFile::Open(fileS,optionS);
    
    for (Int_t week = 0; week < Nweeks; week++)
    {
        if(FitSin22t13)
        {
            ASinChiSquareH[week]->Write(Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,Dm2_ee,week));
            
            delete  ASinChiSquareH[week];
        }
        else
        {
            ADeltaChiSquareH[week]->Write(Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,Dm2_ee,week));
            
            delete ADeltaChiSquareH[week];
        }
    }
    
    SaveSR->Close();
}

void Fitter :: SaveDeltaMRangeChiSquare(Int_t sample, bool FitSin22t13)
{
    std::cout << "\t SAVING D2M_13 Range Fit χ2 HISTOGRAMS" << std::endl;
    
    TString optionM("recreate");
    
    if(sample>0)
    {
        optionM.ReplaceAll("recreate","update");
    }
    
    TString fileM(("./ChiSquare/"+AnalysisString+Form("/Combine%d/Sin22t13Delta2MeeRange.root",Combine)).c_str());
    
    if(!FitSin22t13)
    {
        fileM.ReplaceAll(("./ChiSquare/"+AnalysisString+Form("/Combine%d/Sin22t13Delta2MeeRange.root",Combine)).c_str(),("./ChiSquare/"+AnalysisString+Form("/Combine%d/DeltaMDelta2MeeRange.root",Combine)).c_str());
    }
    
    TFile* SaveMR = TFile::Open(Form(fileM,Combine),optionM);
    
    for (Int_t week = 0; week < Nweeks; week++)
    {
        if(FitSin22t13)
        {
            ASinChiSquareH[week]->Write(Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,Dm2_ee,week));
            
            delete  ASinChiSquareH[week];
        }
        else
        {
            ADeltaChiSquareH[week]->Write(Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,Dm2_ee,week));
            
            delete ADeltaChiSquareH[week];
        }
    }
    
    SaveMR->Close();
}

void Fitter :: SaveSin1DFit(Int_t sample, bool TurnOnBudget, bool TurnOffBudget)
{
    std::cout << "\t SAVING Sin22t13 1D Fit χ2 HISTOGRAMS" << std::endl;
    
    TString optionS1("recreate");
    
    if(sample>0)
    {
        optionS1.ReplaceAll("recreate","update");
    }
    
    if(TurnOnBudget)
    {
        if(VaryAccidentalBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/VaryAccidentalChiSquare.root",Combine)).c_str());
        }
        else if(VaryLiHeBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/VaryLiHeChiSquare.root",Combine)).c_str());
        }
        else if(VaryFastNeutronsBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/VaryFastNeutronsChiSquare.root",Combine)).c_str());
        }
        else if(VaryAmCBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/VaryAmCChiSquare.root",Combine)).c_str());
        }
        else if(DistortLiHeBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/DistortLiHeChiSquare.root",Combine)).c_str());
        }
        else if(DistortFastNeutronsBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/DistortFastNeutronsChiSquare.root",Combine)).c_str());
        }
        else if(DistortAmCBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/DistortAmCChiSquare.root",Combine)).c_str());
        }
        else if(IsotopeBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/IsotopeChiSquare.root",Combine)).c_str());
        }
        else if(ReactorPowerBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/ReactorPowerChiSquare.root",Combine)).c_str());
        }
        else if(RelativeEnergyScaleBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/RelativeEnergyScaleChiSquare.root",Combine)).c_str());
        }
        else if(RelativeEnergyOffsetBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/RelativeEnergyOffsetChiSquare.root",Combine)).c_str());
        }
        else if(AbsoluteEnergyScaleBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/AbsoluteEnergyScaleChiSquare.root",Combine)).c_str());
        }
        else if(AbsoluteEnergyOffsetBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/AbsoluteEnergyOffsetChiSquare.root",Combine)).c_str());
        }
        else if(IAVBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/IAVChiSquare.root",Combine)).c_str());
        }
        else if(NLBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/NLChiSquare.root",Combine)).c_str());
        }
        else if(ResolutionBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/ResolutionChiSquare.root",Combine)).c_str());
        }
        else if(Sin22t12Budget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/Sin22t12ChiSquare.root",Combine)).c_str());
        }
        else if(EfficiencyBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/EfficiencyChiSquare.root",Combine)).c_str());
        }
        else if(SystematicBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/SystematicChiSquare.root",Combine)).c_str());
        }
        else if(BackgroundBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/BackgroundChiSquare.root",Combine)).c_str());
        }
        else if(TotalBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/TotalChiSquare.root",Combine)).c_str());
        }
        else
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/S2/OnlyStatChiSquare.root",Combine)).c_str());
        }
    }
    else if(TurnOffBudget)
    {
        if(VaryAccidentalBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/VaryAccidentalChiSquare.root",Combine)).c_str());
        }
        else if(VaryLiHeBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/VaryLiHeChiSquare.root",Combine)).c_str());
        }
        else if(VaryFastNeutronsBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/VaryFastNeutronsChiSquare.root",Combine)).c_str());
        }
        else if(VaryAmCBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/VaryAmCChiSquare.root",Combine)).c_str());
        }
        else if(DistortLiHeBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/DistortLiHeChiSquare.root",Combine)).c_str());
        }
        else if(DistortFastNeutronsBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/DistortFastNeutronsChiSquare.root",Combine)).c_str());
        }
        else if(DistortAmCBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/DistortAmCChiSquare.root",Combine)).c_str());
        }
        else if(IsotopeBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/IsotopeChiSquare.root",Combine)).c_str());
        }
        else if(ReactorPowerBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/ReactorPowerChiSquare.root",Combine)).c_str());
        }
        else if(RelativeEnergyScaleBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/RelativeEnergyScaleChiSquare.root",Combine)).c_str());
        }
        else if(RelativeEnergyOffsetBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/RelativeEnergyOffsetChiSquare.root",Combine)).c_str());
        }
        else if(AbsoluteEnergyScaleBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/AbsoluteEnergyScaleChiSquare.root",Combine)).c_str());
        }
        else if(AbsoluteEnergyOffsetBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/AbsoluteEnergyOffsetChiSquare.root",Combine)).c_str());
        }
        else if(IAVBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/IAVChiSquare.root",Combine)).c_str());
        }
        else if(NLBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/NLChiSquare.root",Combine)).c_str());
        }
        else if(ResolutionBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/ResolutionChiSquare.root",Combine)).c_str());
        }
        else if(Sin22t12Budget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/Sin22t12ChiSquare.root",Combine)).c_str());
        }
        else if(EfficiencyBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/EfficiencyChiSquare.root",Combine)).c_str());
        }
        else if(SystematicBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/SystematicChiSquare.root",Combine)).c_str());
        }
        else if(BackgroundBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/BackgroundChiSquare.root",Combine)).c_str());
        }
        else if(TotalBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/TotalChiSquare.root",Combine)).c_str());
        }
        else
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/S2/EverythingChiSquare.root",Combine)).c_str());
        }
    }
    else
    {
        sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/1DChiSquareFitSin22t13.root",Combine)).c_str());
    }
    
    
    TFile* SaveS = TFile::Open(FileName,optionS1);
    
    for (Int_t week = 0; week < Nweeks; week++)
    {
        TF1* fitpol2 = new TF1("fitpol2","pol9",Sin22t13-(Sin22t13/5),Sin22t13+(Sin22t13/5));
        ASinChiSquareH[week]->Fit("fitpol2","MLLRQO");
        ASinChiSquareH[week]->Write(Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,Dm2_ee,week));
        TF1* fitcurve = (TF1*)ASinChiSquareH[week]->GetFunction("fitpol2");
        s2t_error_min[week] = (fitcurve->GetX(1.,0,s2t_min[week])-s2t_min[week]);
        
        if(Print)
        {
            TCanvas* SinC = new TCanvas("SinC","SinC");
            ASinChiSquareH[week]->SetTitle("#chi^{2} distribution");
            ASinChiSquareH[week]->GetXaxis()->SetTitleSize(0.03);
            ASinChiSquareH[week]->GetXaxis()->SetTitleOffset(1.3);
            ASinChiSquareH[week]->GetXaxis()->SetTitle("Sin^{2}(2#theta_{13})");
            ASinChiSquareH[week]->GetYaxis()->SetTitleSize(0.25);
            ASinChiSquareH[week]->GetYaxis()->SetTitle("#chi^{2}");

            ASinChiSquareH[week]->SetStats(0);
            ASinChiSquareH[week]->GetFunction("fitpol2")->SetBit(TF1::kNotDraw);
            ASinChiSquareH[week]->Draw();
            
            if(ToyMC&&!StatisticalFluctuation)
            {
                SinC->Print(("./Images/"+AnalysisString+Form("/ChiSquareResults/Combine%d_NoFSinChi2Period_%d.eps",Combine,week)).c_str(),".eps");
            }
            else if(ToyMC)
            {
                SinC->Print(("./Images/"+AnalysisString+Form("/ChiSquareResults/Combine%d_FSinChi2Period_%d.eps",Combine,week)).c_str(),".eps");
            }
            else
            {
                SinC->Print(("./Images/"+AnalysisString+Form("/ChiSquareResults/Combine%d_DataSinChi2Period_%d.eps",Combine,week)).c_str(),".eps");
            }
            delete SinC;
        }
        
        delete ASinChiSquareH[week];
        delete fitpol2;
        std::cout << FileName << " sigma: " << s2t_error_min[week] << " min: " << s2t_min[week] << std::endl;
    }
    
    SaveS->Close();

    
    for (Int_t week = 0; week < Nweeks; week++)
    {
        if((TurnOnBudget)||(TurnOffBudget))
        {
            string filename;
            
            if(TurnOnBudget)
            {
                filename = Form("SinTurnOnErrorBudget_%d.txt",week);
            }
            else
            {
                filename = Form("SinTurnOffErrorBudget_%d.txt",week);
            }
            
            std::fstream file;
            // open file named "data.txt" for writing (std::fstream::app lets you add text to the end of the file)
            file.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);
            // could not open file
            if(!file.is_open())
            {
                std::cout << "Couldn't open file." << endl;
            }
            else
            {
                //  Filename - sigma - best fit
                
                file << FileName << " " << s2t_error_min[week] << " " << s2t_min[week] <<"\n";
                file.close(); // close file
            }
        }
        else
        {
            string filename;
            
            filename = Form("SinFitResults_%d.txt",week);
            
            std::fstream file;
            // open file named "data.txt" for writing (std::fstream::app lets you add text to the end of the file)
            file.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);
            // could not open file
            if(!file.is_open())
            {
                std::cout << "Couldn't open file." << endl;
            }
            else
            {
                //  Filename - sigma - best fit

                file << FileName << " " << s2t_error_min[week] << " " << s2t_min[week] <<"\n";
                file.close(); // close file
            }
        }
    }
    
}

void Fitter :: SaveDM1DFit(Int_t sample, bool TurnOnBudget, bool TurnOffBudget)
{
    std::cout << "\t SAVING D2M_13 1D Fit χ2 HISTOGRAMS" << std::endl;
    
    TString optionM1("recreate");
    if(sample>0)
    {
        optionM1.ReplaceAll("recreate","update");
    }
    
    if(TurnOnBudget)
    {
        if(VaryAccidentalBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/VaryAccidentalChiSquare.root",Combine)).c_str());
        }
        else if(VaryLiHeBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/VaryLiHeChiSquare.root",Combine)).c_str());
        }
        else if(VaryFastNeutronsBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/VaryFastNeutronsChiSquare.root",Combine)).c_str());
        }
        else if(VaryAmCBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/VaryAmCChiSquare.root",Combine)).c_str());
        }
        else if(DistortLiHeBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/DistortLiHeChiSquare.root",Combine)).c_str());
        }
        else if(DistortFastNeutronsBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/DistortFastNeutronsChiSquare.root",Combine)).c_str());
        }
        else if(DistortAmCBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/DistortAmCChiSquare.root",Combine)).c_str());
        }
        else if(IsotopeBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/IsotopeChiSquare.root",Combine)).c_str());
        }
        else if(ReactorPowerBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/ReactorPowerChiSquare.root",Combine)).c_str());
        }
        else if(RelativeEnergyScaleBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/RelativeEnergyScaleChiSquare.root",Combine)).c_str());
        }
        else if(RelativeEnergyOffsetBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/RelativeEnergyOffsetChiSquare.root",Combine)).c_str());
        }
        else if(AbsoluteEnergyScaleBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/AbsoluteEnergyScaleChiSquare.root",Combine)).c_str());
        }
        else if(AbsoluteEnergyOffsetBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/AbsoluteEnergyOffsetChiSquare.root",Combine)).c_str());
        }
        else if(IAVBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/IAVChiSquare.root",Combine)).c_str());
        }
        else if(NLBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/NLChiSquare.root",Combine)).c_str());
        }
        else if(ResolutionBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/ResolutionChiSquare.root",Combine)).c_str());
        }
        else if(Sin22t12Budget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/Sin22t12ChiSquare.root",Combine)).c_str());
        }
        else if(EfficiencyBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/EfficiencyChiSquare.root",Combine)).c_str());
        }
        else if(SystematicBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/SystematicChiSquare.root",Combine)).c_str());
        }
        else if(BackgroundBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/BackgroundChiSquare.root",Combine)).c_str());
        }
        else if(TotalBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/TotalChiSquare.root",Combine)).c_str());
        }
        else
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOnBudget/DM/OnlyStatChiSquare.root",Combine)).c_str());
        }
    }
    else if(TurnOffBudget)
    {
        if(VaryAccidentalBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/VaryAccidentalChiSquare.root",Combine)).c_str());
        }
        else if(VaryLiHeBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/VaryLiHeChiSquare.root",Combine)).c_str());
        }
        else if(VaryFastNeutronsBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/VaryFastNeutronsChiSquare.root",Combine)).c_str());
        }
        else if(VaryAmCBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/VaryAmCChiSquare.root",Combine)).c_str());
        }
        else if(DistortLiHeBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/DistortLiHeChiSquare.root",Combine)).c_str());
        }
        else if(DistortFastNeutronsBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/DistortFastNeutronsChiSquare.root",Combine)).c_str());
        }
        else if(DistortAmCBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/DistortAmCChiSquare.root",Combine)).c_str());
        }
        else if(IsotopeBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/IsotopeChiSquare.root",Combine)).c_str());
        }
        else if(ReactorPowerBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/ReactorPowerChiSquare.root",Combine)).c_str());
        }
        else if(RelativeEnergyScaleBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/RelativeEnergyScaleChiSquare.root",Combine)).c_str());
        }
        else if(RelativeEnergyOffsetBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/RelativeEnergyOffsetChiSquare.root",Combine)).c_str());
        }
        else if(AbsoluteEnergyScaleBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/AbsoluteEnergyScaleChiSquare.root",Combine)).c_str());
        }
        else if(AbsoluteEnergyOffsetBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/AbsoluteEnergyOffsetChiSquare.root",Combine)).c_str());
        }
        else if(IAVBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/IAVChiSquare.root",Combine)).c_str());
        }
        else if(NLBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/NLChiSquare.root",Combine)).c_str());
        }
        else if(ResolutionBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/ResolutionChiSquare.root",Combine)).c_str());
        }
        else if(Sin22t12Budget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/Sin22t12ChiSquare.root",Combine)).c_str());
        }
        else if(EfficiencyBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/EfficiencyChiSquare.root",Combine)).c_str());
        }
        else if(SystematicBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/SystematicChiSquare.root",Combine)).c_str());
        }
        else if(BackgroundBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/BackgroundChiSquare.root",Combine)).c_str());
        }
        if(TotalBudget)
        {
            sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/TurnOffBudget/DM/TotalChiSquare.root",Combine)).c_str());
        }
        
    }
    else
    {
        sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/1DChiSquareFitDm.root",Combine)).c_str());
    }
    
    TFile* SaveDM = TFile::Open(FileName,optionM1);
    
    for (Int_t week = 0; week < Nweeks; week++)
    {
        TF1* fitpol1 = new TF1("fitpol1","pol9",Dm2_ee-(Dm2_ee/5),Dm2_ee+(Dm2_ee/5));
        ADeltaChiSquareH[week]->Fit("fitpol1","MLLRQO");
        ADeltaChiSquareH[week]->Write(Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,Dm2_ee,week));
        TF1* fitcurve1 = (TF1*)ADeltaChiSquareH[week]->GetFunction("fitpol1");
        dm2_ee_error_min[week] = (fitcurve1->GetX(1)-dm2_ee_min[week]);

        std::cout << FileName << " sigma: " << dm2_ee_error_min[week] << " min: " << s2t_min[week] << std::endl;
        
        if(Print)
        {
            TCanvas* DmC = new TCanvas("DmC","DmC");
            
            ADeltaChiSquareH[week]->SetTitle("#chi^{2} distribution");
            ADeltaChiSquareH[week]->GetXaxis()->SetTitleSize(0.03);
            ADeltaChiSquareH[week]->GetXaxis()->SetTitleOffset(1.3);
            ADeltaChiSquareH[week]->GetXaxis()->SetTitle("#Delta^{2}m_{23}");
            ADeltaChiSquareH[week]->GetYaxis()->SetTitleSize(0.25);
            ADeltaChiSquareH[week]->GetYaxis()->SetTitle("#chi^{2}");
            ADeltaChiSquareH[week]->SetStats(0);
            ADeltaChiSquareH[week]->GetFunction("fitpol1")->SetBit(TF1::kNotDraw);
            ADeltaChiSquareH[week]->Draw();
            
            if(ToyMC&&!StatisticalFluctuation)
            {
                DmC->Print(("./Images/"+AnalysisString+Form("/ChiSquareResults/Combine%d_NoFDMChi2Period_%d.eps",Combine,week)).c_str(),".eps");
            }
            else if(ToyMC)
            {
                DmC->Print(("./Images/"+AnalysisString+Form("/ChiSquareResults/Combine%d_FDMChi2Period_%d.eps",Combine,week)).c_str(),".eps");
            }
            else
            {
                DmC->Print(("./Images/"+AnalysisString+Form("/ChiSquareResults/Combine%d_DataDMChi2Period_%d.eps",Combine,week)).c_str(),".eps");
            }
            
            delete DmC;
        }
        
        delete ADeltaChiSquareH[week];
        delete fitpol1;
    }
    
    SaveDM->Close();
    
    
    for (Int_t week = 0; week < Nweeks; week++)
    {
        if((TurnOnBudget)||(TurnOffBudget))
        {
            string filename;
            
            if(TurnOnBudget)
            {
                filename = Form("DMTurnOnErrorBudget_%d.txt",week);
            }
            else
            {
                filename = Form("DMTurnOffErrorBudget_%d.txt",week);
            }
            
            std::fstream file;
            // open file named "data.txt" for writing (std::fstream::app lets you add text to the end of the file)
            file.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);
            // could not open file
            if(!file.is_open())
            {
                std::cout << "Couldn't open file." << endl;
            }
            else
            {
                //  Filename - sigma - best fit
                file << FileName << " " << dm2_ee_error_min[week] << " " << dm2_ee_min[week] <<"\n";
                file.close(); // close file
            }
        }
        else
        {
            string filename;
            
            filename = Form("DMFitResults_%d.txt",week);
            
            std::fstream file;
            // open file named "data.txt" for writing (std::fstream::app lets you add text to the end of the file)
            file.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);
            // could not open file
            if(!file.is_open())
            {
                std::cout << "Couldn't open file." << endl;
            }
            else
            {
                //  Filename - sigma - best fit

                file << FileName << " " << dm2_ee_error_min[week] << " " << dm2_ee_min[week] <<"\n";
                file.close(); // close file
            }
        }
    }
}

void Fitter :: Save2DFit(Int_t sample)
{
    
    std::cout << "\t SAVING 2DFit χ2 HISTOGRAMS" << std::endl;
    
    TString option2D("recreate");
    
    if(sample>0)
    {
        option2D.ReplaceAll("recreate","update");
    }
    
    sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/2DChiSquare.root",Combine)).c_str());
    
    TFile* Save2D = TFile::Open(FileName,option2D);
    
    for (Int_t week = 0; week < Nweeks; week++)
    {
        ChiSquare2DH[week]->SetTitle("2D #chi^{2} map");
        ChiSquare2DH[week]->GetXaxis()->SetTitle("Sin^{2}(2#theta_{13})");
        ChiSquare2DH[week]->GetXaxis()->SetTitleSize(0.25);
        ChiSquare2DH[week]->GetYaxis()->SetTitleSize(0.25);
        ChiSquare2DH[week]->GetYaxis()->SetTitle("#Delta m^{2}_{23}");
        
        ChiSquare2DH[week]->Write();
        
        if(Print)
        {
            TCanvas* Fit2DC = new TCanvas("2DC","2DC");

            ChiSquare2DH[week]->SetStats(0);
            ChiSquare2DH[week]->Draw("colz");
            
            
            if(ToyMC&&!StatisticalFluctuation)
            {
                Fit2DC->Print(("./Images/"+AnalysisString+Form("/ChiSquareResults/Combine%d_NoF2DChi2Period_%d.eps",Combine,week)).c_str(),".eps");
            }
            else if(ToyMC)
            {
                Fit2DC->Print(("./Images/"+AnalysisString+Form("/ChiSquareResults/Combine%d_F2DChi2Period_%d.eps",Combine,week)).c_str(),".eps");
            }
            else
            {
                Fit2DC->Print(("./Images/"+AnalysisString+Form("/ChiSquareResults/Combine%d_Data2DChi2Period_%d.eps",Combine,week)).c_str(),".eps");
            }
            
            delete Fit2DC;
        }
        
        delete ChiSquare2DH[week];

        if(FullGrid)
        {
            for(Int_t stepM=0;stepM < NstepsM;++stepM)
            {
                SinChiSquareH[week][stepM]->Write(Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,stepM*DeltaWidth+dm2_eestart,week));
                delete SinChiSquareH[week][stepM];
            }
            for(Int_t step = 0; step < Nsteps ; step++)
            {
                DeltaChiSquareH[week][step]->Write(Form("Sin%f_Distribution_DeltaM%f_Period%d",step*SinWidth+s22t13start,Dm2_ee,week));
                delete DeltaChiSquareH[week][step];
            }
        }
    }
    
    Save2D->Close();
}

void Fitter::TestRandomExperimentsChi2(bool Minuit, bool Fit2D, bool FitSin22t13,FitterGui* FitterGui,Int_t FakeExperiments,Int_t GridSteps)
{
    Mode=1;//random experiments.
    // run 1000 experiments in 21 grid points
    NominalData* ndata = Pred->GetNominalData();
    ndata->SetAllRandomSystematics(1);//all systematics random
    Pred->SetNominalData(ndata);
    string SaveTestsS;
    
    for (Int_t week = 0; week < Nweeks; week++)
    {
        chi2_min[week] = 1e10;
        s2t_min[week] = 1e10;
        s2t_error_min[week] = 1e10;
        dm2_ee_min[week]= 1e10;
    }
    
    degenerancy = 0;
    
    std::cout << " " << std::endl;
    std::cout << "\t \t \t STEP " << (s22t13end-s22t13start)*1./(Nsteps-1) << std::endl;
    std::cout << "\t \t \t STEPM " <<  (dm2_eeend-dm2_eestart)*1./(NstepsM-1) << std::endl;
    std::cout << " " << std::endl;
    
    DeltaWidth =(dm2_eeend-dm2_eestart)*1./(NstepsM-1);
    SinWidth = (s22t13end-s22t13start)*1./(Nsteps-1);
    
    RealSinWidth = (s22t13end-s22t13start)*1./(GridSteps-1);
    RealDeltaWidth = (dm2_eeend-dm2_eestart)*1./(GridSteps-1);

    for (Int_t week = 0; week < Nweeks; week++)
    {
        chi2=0;
        
        if(ReadTxt)
        {
            Pred->LoadTxtCovarianceMatrices(week);
        }
        else
        {
            Pred->LoadRootCovarianceMatrices(week);
        }
        
        if(Fit2D)
        {
            //  Fit 101 x 101 grid in Sin22t13 & DM2_ee phase-space (10,201 Fits!)
            //  Generating and saving Toy MC samples in a tree once avoids to calculate this large number of fits every time I test the fitter.
            
            SaveTestsS = ("./ChiSquare/"+AnalysisString+Form("/Combine%d/Tests/Test2D.root",Combine));
            
            ChiSquare2DH[week] = new TH2D("2DChiDistribution","#chi^{2} Distribution",Nsteps,s22t13start,s22t13end,NstepsM,dm2_eestart,dm2_eeend);
            
            for(Int_t step=0;step < Nsteps;++step)
            {
                sin22t13[step] = SinWidth*step + s22t13start;
                
                std::cout << "\t sin22t13 " << sin22t13[step] << std::endl;
                
                for(Int_t stepM = 0; stepM < NstepsM ; stepM++)
                {
                    dm2_ee[stepM] = DeltaWidth*stepM + dm2_eestart;
                    
                    std::cout << "\t dm2_ee " << dm2_ee[stepM] << std::endl;
                    
                    chi2 = Pred->CalculateChi2(sin22t13[step],dm2_ee[stepM],week,ToyMC);
                    
                    ChiSquare2DH[week]->SetBinContent(step+1,stepM+1,TMath::Abs(chi2));
                    
                    if (TMath::Abs(chi2)<chi2_min[week])
                    {
                        s2t_min[week]=sin22t13[step];
                        dm2_ee_min[week]=dm2_ee[stepM];
                        chi2_min[week]=TMath::Abs(chi2);
                    }
                    else if(TMath::Abs(chi2)==chi2_min[week])
                    {
                        degenerancy++;
                    }
                    
                    if((chi2-chi2_min[week])<=chi2limit)//1sigma cut
                    {
                        s2t_error_min[week] = sin22t13[step];
                        dm2_ee_error_min[week] = dm2_ee[step];
                    }
                    
                    FitterGui->UpdateFitter(step*NstepsM+stepM);
                }
            }
            
            if(FullGrid)//  Produces 202 (2*MaxSamples) cut histograms
            {
                for(Int_t stepM=1;stepM <=NstepsM;++stepM)
                {
                    SinChiSquareH[week][stepM]=(TH1D*)ChiSquare2DH[week]->ProjectionX(Form("Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,dm2_ee[stepM],week),stepM+1,stepM+1);
                }
                for(Int_t step = 1; step <= Nsteps ; step++)
                {
                    DeltaChiSquareH[week][step]=(TH1D*)ChiSquare2DH[week]->ProjectionY(Form("Sin%f_Distribution_DeltaM%f_Period%d",sin22t13[step],Dm2_ee,week),step+1,step+1);
                }
            }
        }
        else if (FitSin22t13)
        {
            //  Fit 101 different Sin22t13 with a fix Dm2_ee
            
            SaveTestsS = ("./ChiSquare/"+AnalysisString+Form("/Combine%d/Tests/TestSin22t13.root",Combine));

            dm2_ee_min[week] = Dm2_ee;
            
            for(Int_t gridstep=0;gridstep < GridSteps;++gridstep)
            {
                TestSinChiSquareH[week][gridstep] = new TH1D(Form("#chi^{2}_Distribution_Period%d_Step%d",week,gridstep),Form("#chi^{2}_Distribution_Period%d_Step%d",week,gridstep),gridstep,s22t13start,s22t13end);
                
                Realsin22t13[gridstep] = RealSinWidth*gridstep + s22t13start;

                Pred->SetDM213(Dm2_ee);
                Pred->SetSin22t13(Realsin22t13[gridstep]);
                
                Pred->SetExperiment(0); // load different far data for each experiment
                
                Pred->LoadData(week,1,GridSteps,Mode);// load near data a first time for generating the inverse matrix
                
                Pred->GenerateInverseMatrix(Realsin22t13[gridstep], Dm2_ee, week, 0, ToyMC,GridSteps);
                
                std::cout << " Sin inside varied data tree " << Realsin22t13[gridstep] << " and Δm2ee: " << Dm2_ee << std::endl;

                for(Int_t experiment = 0;experiment<FakeExperiments;experiment++)
                {
                    Pred->SetExperiment(experiment);
                    
                    Pred->LoadData(week,1,GridSteps,Mode);//always ToyMC
                    
                    degenerancy = 0;

                    for(Int_t step=0;step < Nsteps;++step)
                    {
                        sin22t13[step] = SinWidth*step + s22t13start;

                        chi2 = Pred->CalculateChi2(sin22t13[step],Dm2_ee,week,1);
                        
                        if (TMath::Abs(chi2)<chi2_min[week])
                        {
                            s2t_min[week]=sin22t13[step];
                            chi2_min[week]=TMath::Abs(chi2);
                        }
                        else if(TMath::Abs(chi2)==chi2_min[week])
                        {
                            degenerancy++;
                        }
                        FitterGui->UpdateFitter(step);
                    }
                    
                    TestSinChiSquareH[week][gridstep]->Fill(s2t_min[week]);
                }
                
                Pred->DeleteData();
            }
        }
        else
        {
            //  Fit 101 different DM2_ee with a fix Sin22t13
            
            SaveTestsS = ("./ChiSquare/"+AnalysisString+Form("/Combine%d/Tests/TestDM2ee.root",Combine));

            s2t_min[week] = Sin22t13;
            
            for(Int_t gridstep=0;gridstep < GridSteps;++gridstep)
            {
                TestDeltaChiSquareH[week][gridstep] = new TH1D(Form("#chi^{2}_Distribution_Period%d_Step%d",week,gridstep),Form("#chi^{2}_Distribution_Period%d_Step%d",week,gridstep),gridstep,s22t13start,s22t13end);
                
                Realdm2_ee[gridstep] = RealDeltaWidth*gridstep + dm2_eestart;
                
                Pred->SetDM213(Realdm2_ee[gridstep]);
                Pred->SetSin22t13(Sin22t13);
                
                Pred->SetExperiment(0); // load different far data for each experiment

                Pred->LoadData(week,1,GridSteps,Mode);// load near data a first time for generating the inverse matrix

                Pred->GenerateInverseMatrix(Sin22t13, Realdm2_ee[gridstep], week, 0, ToyMC,GridSteps);

                for(Int_t experiment = 0;experiment<FakeExperiments;experiment++)
                {
                    Pred->SetExperiment(experiment); // load different far data for each experiment

                    Pred->LoadData(week,1,GridSteps,Mode);//always ToyMC with variations

                    degenerancy = 0;

                    for(Int_t stepM = 0; stepM < NstepsM ; stepM++)
                    {
                        dm2_ee[stepM] = DeltaWidth*stepM + dm2_eestart;
                        
                        chi2 = Pred->CalculateChi2(Sin22t13,dm2_ee[stepM],week,1);
                        
                        if (TMath::Abs(chi2)<chi2_min[week])
                        {
                            dm2_ee_min[week]=dm2_ee[stepM];
                            chi2_min[week]=TMath::Abs(chi2);
                        }
                        else if(TMath::Abs(chi2)==chi2_min[week])
                        {
                            degenerancy++;
                        }
                        FitterGui->UpdateFitter(stepM);
                    }
                    TestDeltaChiSquareH[week][gridstep]->Fill(dm2_ee_min[week]);
                }
                
                Pred->DeleteData();
            }
        }
        
        std::cout << "\t \t \t For Period: " << week << std::endl;
        std::cout << "\t \t \t \t χ2 min" << chi2_min[week] << std::endl;
        std::cout << "\t \t \t \t Sin22t13 min" << s2t_min[week] << std::endl;
        //            std::cout << "\t \t \t \t Sin22t13 error" << s2t_error_min[week] << std::endl;
        std::cout << "\t \t \t \t Δm2_ee_min min" << dm2_ee_min[week] << std::endl;
        //            std::cout << "\t \t \t \t Δm2_ee error" << dm2_ee_error_min[week] << std::endl;
        std::cout << "\t \t \t \t χ2 degenerated " << degenerancy << " times" << std::endl;
        
        Pred->DeleteMatrices();
    }

    delete Pred;
    
    TFile* SaveTests = new TFile(SaveTestsS.c_str(),"recreate");
    
    for(Int_t gridstep=0;gridstep < GridSteps;++gridstep)
    {
        for (Int_t week = 0; week < Nweeks; week++)
        {
            if(Fit2D)
            {
                ChiSquare2DH[week]->Write();
            }
            else if(FitSin22t13)
            {
                TestSinChiSquareH[week][gridstep]->Write();
                
            }
            else
            {
                TestDeltaChiSquareH[week][gridstep]->Write();
            }
        }
    }
    delete SaveTests;
    std::cout << "\t Finished Scanning parameters" << std::endl;

}

void Fitter :: SaveChiTest(bool Fit2D, bool FitSin22t13,Int_t steps)
{
    for (Int_t week = 0; week < Nweeks; week++)
    {
        for(Int_t gridstep=0;gridstep < steps;++gridstep)
        {
            if(Fit2D)
            {
                //Not coded yet.
            }
            else if (FitSin22t13)
            {
                TestSinChiSquareH[week][gridstep]->Write();
            }
            else
            {
                TestDeltaChiSquareH[week][gridstep]->Write();
            }
        }
    }
}