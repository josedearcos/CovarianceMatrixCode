#pragma once
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include <math.h>
#include <iostream>
#include "NominalData.h"
#include "TROOT.h"

const Int_t BGNDs = 9;

class FitBackgrounds2
{
private:
    NominalData* Nom;
    std::string AnalysisString;
    std::string SaveString;
    TF1* HAmCFunc;
    TF1* GdAmCFunc;
    bool isH;
    Int_t XCellLimit;
    Int_t YCellLimit;
    
    Int_t Nweeks;
    Int_t NADs;
    Int_t hall;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    
    //Histograms
    TH1D* AccidentalsH[MaxDetectors];
    TH1D* AmCH;
    TH1D* BackgroundsH[BGNDs];
    
    TH1D* FinalAccidentalsH[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    TH1D* FinalLiHeH[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    TH1D* FinalFastNeutronsH[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    TH1D* FinalAmCH[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    
    TH1D* NearBackgroundSpectrumH[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    TH1D* FarBackgroundSpectrumH[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    
    //  Background rates (with AD effects included)
    Double_t ScaleAcc[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    Double_t ScaleLiHe[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    Double_t ScaleFN[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    Double_t ScaleAmC[MaxDetectors][MaxPeriods][VolumeX][VolumeY];
    
    //Binning parameters:
    Int_t n_evis_bins;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    
    void FitAmc();
    void FitFastNeutrons();
    void PrintBackgrounds();
    void SaveAndDelete();
    void ReadHBackgrounds();
    void ReadGdBackgrounds();
public:
    
    FitBackgrounds2();
    FitBackgrounds2(NominalData*);
    void ProduceNominalBackgroundHistograms();

};

FitBackgrounds2 :: FitBackgrounds2()
{
    Nom = new NominalData(0,2);
    Nweeks = Nom->GetWeeks();
    NADs = Nom->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    isH = Nom->GetAnalysis();
    
    if(isH)
    {
        AnalysisString = "Hydrogen";
        SaveString = "HBackground";
        XCellLimit = VolumeX;
        YCellLimit = VolumeY;
    }
    else
    {
        AnalysisString = "Gadolinium";
        SaveString = "GDBackground";
        XCellLimit = 1;
        YCellLimit = 1;
    }
    if(NADs == 8)//    ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    //  Binning variables
    InitialVisibleEnergy = Nom->GetEVisMin();
    FinalVisibleEnergy =  Nom->GetEVisMax();
    
    n_evis_bins = Nom->GetVisibleBins();
    
    for (Int_t i = 0; i <= n_evis_bins; i++)
    {
        evis_bins[i] = Nom->GetVisibleBinningArray(i);
    }
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    ScaleAcc[AD][week][idx][idy]=Nom->GetAccidentalEvents(AD,week,idx,idy);//AD efficiencies included
                    ScaleLiHe[AD][week][idx][idy]=Nom->GetLiHeEvents(AD,week,idx,idy);
                    ScaleFN[AD][week][idx][idy]=Nom->GetFNEvents(AD,week,idx,idy);
                    ScaleAmC[AD][week][idx][idy]=Nom->GetAmCEvents(AD,week,idx,idy);
                }
            }
        }
    }
    delete Nom;
}

FitBackgrounds2 :: FitBackgrounds2(NominalData* Data)
{
    Nweeks = Data->GetWeeks();
    
    NADs = Data->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    isH = Data->GetAnalysis();
    if(isH)
    {
        AnalysisString = "Hydrogen";
        SaveString = "HBackground";
        XCellLimit = VolumeX;
        YCellLimit = VolumeY;
    }
    else
    {
        AnalysisString = "Gadolinium";
        SaveString = "GDBackground";
        XCellLimit = 1;
        YCellLimit = 1;
    }
    if(NADs == 8)//    ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    n_evis_bins = Data->GetVisibleBins();
    
    for (Int_t i = 0; i <= n_evis_bins; i++)
    {
        evis_bins[i] = Data->GetVisibleBinningArray(i);
    }
    
    InitialVisibleEnergy = Data->GetEVisMin();
    FinalVisibleEnergy = Data->GetEVisMax();
    
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    ScaleAcc[AD][week][idx][idy]=Data->GetAccidentalEvents(AD,week,idx,idy);//AD efficiencies included
                    ScaleLiHe[AD][week][idx][idy]=Data->GetLiHeEvents(AD,week,idx,idy);
                    ScaleFN[AD][week][idx][idy]=Data->GetFNEvents(AD,week,idx,idy);
                    ScaleAmC[AD][week][idx][idy]=Data->GetAmCEvents(AD,week,idx,idy);
                }
            }
        }
    }
}

//This works for the current given Backgrounds and 6 ADs, will have to be changed in the case of 8 ADs or diverse inputs.
void FitBackgrounds2 :: ReadHBackgrounds()
{
    std::cout << "Reading Hydrogen Backgrounds" << std::endl;
    
    Int_t j=0;
    
    for(Int_t AD = 0; AD < NADs; AD++)
    {
        TFile* AccidentalsF = TFile::Open(("./BackgroundSpectrum/"+SaveString+Form("/Bkg_Accidental_AD%i.root", AD+1)).c_str());
        
        if(AD==2)
        {
            AccidentalsH[AD]=(TH1D*)gDirectory->Get("h1dAccBkg22_1");
        }
        else if(AD<2)
        {
            AccidentalsH[AD]=(TH1D*)gDirectory->Get(Form("h1dAccBkg22_%i",j+1));
            j++;
        }
        else if(AD>2)
        {
            AccidentalsH[AD]=(TH1D*)gDirectory->Get(Form("h1dAccBkg22_%i",j-1));
            j++;
        }
        
        if(LoganBinning)//240 but accidentals files have 120, need to interpolate!
        {
            BackgroundsH[AD] = new TH1D(Form("Accidentals_AD%i",AD+1),Form("Accidentals_AD%i",AD+1),n_evis_bins,evis_bins);
            
            BackgroundsH[AD]->Reset();
            
            //         cout<<"NPoints "<<AccidentalsH[AD]->GetXaxis()->GetNbins();
            for (Int_t visbin = 1; visbin <= n_evis_bins; visbin++)
            {
                BackgroundsH[AD]->SetBinContent(visbin,AccidentalsH[AD]->Interpolate(BackgroundsH[AD]->GetXaxis()->GetBinCenter(visbin)));
                
//                std::cout << "Vis Bin: " << visbin << " - Center Bin: " << BackgroundsH[AD]->GetXaxis()->GetBinCenter(visbin) << std::endl;
//                std::cout << "interpolation value: " << AccidentalsH[AD]->Interpolate(BackgroundsH[AD]->GetXaxis()->GetBinCenter(visbin)) << std::endl;
//                std::cout << "evis_bins: " << evis_bins[visbin] << std::endl;

            }
        }
        else
        {
            BackgroundsH[AD]=(TH1D*)AccidentalsH[AD]->Clone(Form("Accidentals_AD%i",AD+1));
        }
        AccidentalsF->Close();
    }
    
    FitFastNeutrons();
    
    FitAmc();//Exponential fit
    
    //9Li8He
    TFile* LiF = TFile::Open("./BackgroundSpectrum/HBackground/Li9_beta_neutron.root");
    TH1D* LithiumH = (TH1D*)LiF->Get("h");
    delete LiF;//90%
    
    TFile* HF = TFile::Open("./BackgroundSpectrum/HBackground/He8_beta_neutron.root");
    TH1D* HeliumH = (TH1D*)HF->Get("h");
    delete HF;//10%
    
  //  TH1D* LiHeBackgroundH = (TH1D*)LithiumH->Clone();
   // LiHeBackgroundH->Scale(.9);
    

    
   // LiHeBackgroundH->Add(HeliumH,.1);
    
    BackgroundsH[NADs] = new TH1D("LiHe","LiHe",n_evis_bins,evis_bins);
    TH1D* LinearBackgroundsH = new TH1D("LinearLiHe","LiHe",240,0,12);
    for (Int_t visbin = 1; visbin <= 240; visbin++)
    {
        LinearBackgroundsH->SetBinContent(visbin,(.9)*LithiumH->Interpolate(LinearBackgroundsH->GetXaxis()->GetBinCenter(visbin))+(.1)*HeliumH->Interpolate(HeliumH->GetXaxis()->GetBinCenter(visbin)));
    }
    
    BackgroundsH[NADs] = (TH1D*)LinearBackgroundsH->Rebin(n_evis_bins,"LiHe",evis_bins);

    delete LinearBackgroundsH;
//    for (Int_t visbin = 1; visbin <= n_evis_bins; visbin++)
//    {
//        BackgroundsH[NADs]->SetBinContent(visbin,(.9)*LithiumH->Interpolate(BackgroundsH[NADs]->GetXaxis()->GetBinCenter(visbin))+(.1)*HeliumH->Interpolate(HeliumH->GetXaxis()->GetBinCenter(visbin)));
//
//        
//      //  BackgroundsH[NADs]->SetBinContent(visbin,LiHeBackgroundH->Interpolate(BackgroundsH[NADs]->GetXaxis()->GetBinCenter(visbin)));
//        
//        //                std::cout << "Vis Bin: " << visbin << " - Center Bin: " << BackgroundsH[AD]->GetXaxis()->GetBinCenter(visbin) << std::endl;
//        //                std::cout << "interpolation value: " << AccidentalsH[AD]->Interpolate(BackgroundsH[AD]->GetXaxis()->GetBinCenter(visbin)) << std::endl;
//        //                std::cout << "evis_bins: " << evis_bins[visbin] << std::endl;
//        
//    }
    
    delete LithiumH;
    delete HeliumH;
   // delete LiHeBackgroundH;
    
    #ifdef PrintEps
    TCanvas* BGNDC = new TCanvas("","");
    BGNDC->Divide(3,3);
    
    for(Int_t AD = 0; AD < NADs+3; AD++)
    {
        BGNDC->cd(AD+1);
        // AccidentalsH[AD]->Draw();
        BackgroundsH[AD]->Draw();
    }
    BGNDC->Print("./Images/Hydrogen/BackgroundVariations/HydrogenBackgrounds.eps");
    
    delete BGNDC;
    #endif
}

void FitBackgrounds2 :: ReadGdBackgrounds()
{
    std::cout << "Reading Gadollnium Backgrounds" << std::endl;
    
    Char_t filenameBackgrounds[100];
    
    if(ADSimple) //P12E ADSimple
    {
        sprintf(filenameBackgrounds,("./BackgroundSpectrum/"+SaveString+"/accidental_eprompt_shapes.root").c_str());//ADSimple
    }
    else //P14 ADScaled
    {
        sprintf(filenameBackgrounds,("./BackgroundSpectrum/"+SaveString+"/IHEP_accidental_lbnlbin_6AD.root").c_str());//ADScaled
    }
    
    TFile* BackgroundsF = TFile::Open(filenameBackgrounds);
    
    for(Int_t AD = 0; AD < NADs; AD++)
    {
        if(ADSimple)//P12E ADSimple
        {
            BackgroundsH[AD]=(TH1D*)gDirectory->Get(Form("h_accidental_eprompt_fine_inclusive_ad%i",AD+1));//ADSimple
        }
        else//P14 ADScaled
        {
            if(AD<ADsEH1)
            {
                BackgroundsH[AD]=(TH1D*)gDirectory->Get(Form("h_accidental_eprompt_fine_inclusive_eh1_ad%i",AD+1));//ADScaled
            }
            else if(AD<(ADsEH1+ADsEH2))
            {
                BackgroundsH[AD]=(TH1D*)gDirectory->Get(Form("h_accidental_eprompt_fine_inclusive_eh2_ad%i",AD-ADsEH1+1));//ADScaled
                
            }
            else
            {
                BackgroundsH[AD]=(TH1D*)gDirectory->Get(Form("h_accidental_eprompt_fine_inclusive_eh3_ad%i",AD-ADsEH1-ADsEH2+1));//ADScaled
            }
        }
    }
    
    BackgroundsF->Close();
    
    TFile* BackgroundsF1 = TFile::Open(("./BackgroundSpectrum/"+SaveString+"/li9_spectrum.root").c_str());
    
    BackgroundsH[NADs]=(TH1D*)gDirectory->Get("h_li9_smeared_toy");
    
    BackgroundsF1->Close();
    
    TFile* BackgroundsF2 = TFile::Open(("./BackgroundSpectrum/"+SaveString+"/fn_spectrum.root").c_str());
    
    BackgroundsH[NADs+1]=(TH1D*)gDirectory->Get("h_toy");
    
    BackgroundsF2->Close();
    
    TFile* BackgroundsF3 = TFile::Open(("./BackgroundSpectrum/"+SaveString+"/amc_spectrum.root").c_str());
    GdAmCFunc=(TF1*)gDirectory->Get("expo");
    GdAmCFunc->SetRange(InitialVisibleEnergy,FinalVisibleEnergy);//InitialVisibleEnergy = 0.7 before
    GdAmCFunc->SetName("AmCFunc");
    BackgroundsH[NADs+2]=(TH1D*)gDirectory->Get("h_toy");
    BackgroundsH[NADs+2]->Reset();
    BackgroundsH[NADs+2]->Add(GdAmCFunc);
    BackgroundsF3->Close();
    
    #ifdef PrintEps
        TCanvas* GdAmCC = new TCanvas ("AmCC","AmCC");
        GdAmCFunc->SetTitle("AmC Background");
        GdAmCFunc->SetLineColor(kBlue);
        GdAmCFunc->Draw();
        GdAmCC->Print("./Images/Gadolinium/BackgroundVariations/LinearBinningAmC.eps");
        delete GdAmCC;
        
        TCanvas* GdFNC = new TCanvas ("FNC","FNC");
        BackgroundsH[NADs+1]->Draw();
        GdFNC->Print("./Images/Gadolinium/BackgroundVariations/LinearBinningFN.eps");
        delete GdFNC;
        
        TCanvas* GdLiHeC = new TCanvas ("LiHeC","LiHeC");
        BackgroundsH[NADs]->Draw();
        GdLiHeC->Print("./Images/Gadolinium/BackgroundVariations/LinearBinningLiHe.eps");
        delete GdLiHeC;
    #endif
    
    std::cout << "Finished Reading Gadolinium Backgrounds" << std::endl;
    
}

void FitBackgrounds2::SaveAndDelete()
{
    //Save Backgrounds
    for(Int_t BGND = 0; BGND < BGNDs; BGND++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for (Int_t AD = 0; AD<NADs; AD++)
            {
                for(Int_t idx=0; idx<XCellLimit; idx++)
                {
                    for(Int_t idy=0; idy<YCellLimit; idy++)
                    {
                        if(BGND<NADs)//Accidentals
                        {
                            FinalAccidentalsH[AD][week][idx][idy]=(TH1D*)BackgroundsH[BGND]->Clone(Form("Accidentals_AD%i_Volume%i",AD,idx));
                            FinalAccidentalsH[AD][week][idx][idy]=(TH1D*)FinalAccidentalsH[AD][week][idx][idy]->Rebin(n_evis_bins,Form("Accidentals_AD%i_Volume%i",AD,idx),evis_bins);
                            
                            FinalAccidentalsH[AD][week][idx][idy]->Scale(ScaleAcc[AD][week][idx][idy]/FinalAccidentalsH[AD][week][idx][idy]->Integral());
                            FinalAccidentalsH[AD][week][idx][idy]->SetTitle(Form("Accidentals_AD%i_Volume%i",AD,idx));
                          //  std::cout << " SCALE ACCIDENTALS: " << ScaleAcc[AD][week][idx][idy]  << " " <<  FinalAccidentalsH[AD][week][idx][idy]->Integral() << std::endl;
                        }
                        else if(BGND==NADs)//LiHe
                        {
                            
                            FinalLiHeH[AD][week][idx][idy]=(TH1D*)BackgroundsH[BGND]->Clone(Form("LiHe_AD%i_Volume%i",AD,idx));
                            FinalLiHeH[AD][week][idx][idy]=(TH1D*)FinalLiHeH[AD][week][idx][idy]->Rebin(n_evis_bins,Form("LiHe_AD%i_Volume%i",AD,idx),evis_bins);
                            FinalLiHeH[AD][week][idx][idy]->Scale(ScaleLiHe[AD][week][idx][idy]/FinalLiHeH[AD][week][idx][idy]->Integral());
                            FinalLiHeH[AD][week][idx][idy]->SetTitle("LiHe");
                           // std::cout << " SCALE LI-HE: " << ScaleLiHe[AD][week][idx][idy]  << " " <<  FinalLiHeH[AD][week][idx][idy]->Integral() << std::endl;
                            
                        }
                        else if(BGND==NADs+1)//FN
                        {
                            FinalFastNeutronsH[AD][week][idx][idy]=(TH1D*)BackgroundsH[BGND]->Clone(Form("FN_AD%i_Volume%i",AD,idx));
                            FinalFastNeutronsH[AD][week][idx][idy]=(TH1D*)FinalFastNeutronsH[AD][week][idx][idy]->Rebin(n_evis_bins,Form("FN_AD%i_Volume%i",AD,idx),evis_bins);
                            FinalFastNeutronsH[AD][week][idx][idy]->Scale(ScaleFN[AD][week][idx][idy]/FinalFastNeutronsH[AD][week][idx][idy]->Integral());
                            
                            FinalFastNeutronsH[AD][week][idx][idy]->SetTitle("FN");
                           // std::cout << " SCALE FN : " << ScaleFN[AD][week][idx][idy]  << " " <<  FinalFastNeutronsH[AD][week][idx][idy]->Integral() << std::endl;
                            
                        }
                        else//AmC
                        {
                            FinalAmCH[AD][week][idx][idy]=(TH1D*)BackgroundsH[BGND]->Clone(Form("AmC_AD%i_Volume%i",AD,idx));
                            FinalAmCH[AD][week][idx][idy]=(TH1D*)FinalAmCH[AD][week][idx][idy]->Rebin(n_evis_bins,Form("AmC_AD%i_Volume%i",AD,idx),evis_bins);
                            FinalAmCH[AD][week][idx][idy]->Scale(ScaleAmC[AD][week][idx][idy]/FinalAmCH[AD][week][idx][idy]->Integral());
                            FinalAmCH[AD][week][idx][idy]->SetTitle("AmC");
                           // std::cout << " SCALE AMC : " << ScaleAmC[AD][week][idx][idy]  << " " <<  FinalAmCH[AD][week][idx][idy]->Integral() << std::endl;
                            
                        }
                    }
                }
            }
        }
    }
    
    TFile* BackgroundsF4 = TFile::Open(("./BackgroundSpectrum/"+SaveString+"/Backgrounds.root").c_str(),"recreate");
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    if(AD<ADsEH1+ADsEH2)
                    {
                        NearBackgroundSpectrumH[AD][week][idx][idy] = new TH1D(Form("Near AD%i Nominal Background Period%i_Volume%i",AD,week,idx),Form("Near AD%i Nominal Background Period%i_Volume%i",AD,week,idx),n_evis_bins,evis_bins);
                        NearBackgroundSpectrumH[AD][week][idx][idy]->Add(FinalAccidentalsH[AD][week][idx][idy]);
                        NearBackgroundSpectrumH[AD][week][idx][idy]->Add(FinalLiHeH[AD][week][idx][idy]);
                        NearBackgroundSpectrumH[AD][week][idx][idy]->Add(FinalFastNeutronsH[AD][week][idx][idy]);
                        NearBackgroundSpectrumH[AD][week][idx][idy]->Add(FinalAmCH[AD][week][idx][idy]);
                        NearBackgroundSpectrumH[AD][week][idx][idy]->Write();
                    }
                    else
                    {
                        FarBackgroundSpectrumH[AD][week][idx][idy] = new TH1D(Form("Far AD%i Nominal Background Period%i_Volume%i",AD,week,idx),Form("Far AD%i Nominal Background Period%i_Volume%i",AD,week,idx),n_evis_bins,evis_bins);
                        FarBackgroundSpectrumH[AD][week][idx][idy]->Add(FinalAccidentalsH[AD][week][idx][idy]);
                        FarBackgroundSpectrumH[AD][week][idx][idy]->Add(FinalLiHeH[AD][week][idx][idy]);
                        FarBackgroundSpectrumH[AD][week][idx][idy]->Add(FinalFastNeutronsH[AD][week][idx][idy]);
                        FarBackgroundSpectrumH[AD][week][idx][idy]->Add(FinalAmCH[AD][week][idx][idy]);
                        FarBackgroundSpectrumH[AD][week][idx][idy]->Write();
                    }
                    
                    FinalAccidentalsH[AD][week][idx][idy]->Write();
                    FinalLiHeH[AD][week][idx][idy]->Write();
                    FinalFastNeutronsH[AD][week][idx][idy]->Write();
                    FinalAmCH[AD][week][idx][idy]->Write();
                }
            }
        }
    }
    
    if(isH)
    {
        HAmCFunc->Write();
        delete HAmCFunc;
    }
    else
    {
        GdAmCFunc->Write();
        delete GdAmCFunc;
    }

    BackgroundsF4->Close();
    
    //print
    #ifdef PrintEps
        PrintBackgrounds();
    #endif
    //delete
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        if(isH)
        {
            delete AccidentalsH[AD];
        }
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    delete FinalAccidentalsH[AD][week][idx][idy];
                    delete FinalLiHeH[AD][week][idx][idy];
                    delete FinalFastNeutronsH[AD][week][idx][idy];
                    delete FinalAmCH[AD][week][idx][idy];
                    
                    if(AD<ADsEH1+ADsEH2)
                    {
                        delete NearBackgroundSpectrumH[AD][week][idx][idy];
                    }
                    else
                    {
                        delete FarBackgroundSpectrumH[AD][week][idx][idy];
                    }
                }
            }
        }
    }
    
    for(Int_t BGND = 0; BGND<BGNDs; BGND++)
    {
        delete BackgroundsH[BGND];
    }
}

void FitBackgrounds2::FitFastNeutrons()
{
    bool FlatNeutron = 0;
    if(FlatNeutron)
    {
        TFile* FastNeutronF = TFile::Open(("./BackgroundSpectrum/"+SaveString+"/Bkg_FastNeutron_EH1_MC.root").c_str());
        TH1D* FastNeutronH=(TH1D*)gDirectory->Get("hist_FastNeutron");
        FastNeutronF->Close();
        
        FastNeutronH->Fit("pol0","Q0","",10,100);
        TF1* FitFNResult = FastNeutronH->GetFunction("pol0");
        Double_t AveragePol0 = FitFNResult->GetParameter(0);
        
        FastNeutronH->Fit("pol1","+Q0","",10,100);
        TF1* polynomial = FastNeutronH->GetFunction("pol1");
        
        Double_t ConstantPol1 = polynomial->GetParameter(0);
        Double_t SlopePol1 = polynomial->GetParameter(1);
        Double_t AveragePol1=ConstantPol1+SlopePol1*55; //Value of the function in the middle point of the range 55 = ((100-10)/2)
        
        Double_t TotalAverage = (AveragePol0 + AveragePol1)/2;
        
        BackgroundsH[NADs+1] = new TH1D("FN","FN",n_evis_bins,evis_bins);//60 so the bins are 0.2
        
        TH1D* LinearFN = new TH1D("LinearFN","FN",240,0,12);
        
        for(Int_t i=1;i<=240;i++)
        {
            LinearFN->SetBinContent(i,TotalAverage);
        }
        
        BackgroundsH[NADs+1] = (TH1D*)LinearFN->Rebin(n_evis_bins,"FN",evis_bins);
        
        delete LinearFN;

        delete FastNeutronH;
    }
    else
    {
        TF1* FNCurve = new TF1("FastNeutron","pow(x,[0])",0,12);
        FNCurve->SetParameter(0,-0.639);//From Xiang Pan Fast Neutron = E^{-P} where P = 0.639±0.036
        
        TH1D* LinearFN = new TH1D("LinearFN","FN",240,0,12);
        
        LinearFN->Add(FNCurve);
        
        BackgroundsH[NADs+1] = (TH1D*)LinearFN->Rebin(n_evis_bins,"FN",evis_bins);
        
        delete LinearFN;
        
        delete FNCurve;
    }
//    delete polynomial;
//    delete FitFNResult;
}

void FitBackgrounds2 :: FitAmc()
{
    TFile* AmCF1 = TFile::Open(("./BackgroundSpectrum/"+SaveString+"/Bkg_StrongAmC.root").c_str());
    TH1D* OriginalAmCH=(TH1D*)gDirectory->Get("hist_Bkg_StrongAmC");
    AmCF1->Close();
    
    if(LoganBinning)
    {
        AmCH = new TH1D("hist_Bkg_StrongAmC","hist_Bkg_StrongAmC",n_evis_bins,evis_bins);
        
        for(Int_t i=1;i<=n_evis_bins;i++)
        {
            AmCH->SetBinContent(i,(OriginalAmCH->Interpolate(AmCH->GetXaxis()->GetBinCenter(i))));
            // std::cout << AmCH->GetBinContent(i) << std::endl;
            if((AmCH->GetBinContent(i))<0)
            {
                AmCH->SetBinContent(i,-1*AmCH->GetBinContent(i)); //Invert negative bins??
                //                          OR
                // AmCH->SetBinContent(i,0);//Set negative bins to 0?
            }
        }
    }
    else
    {
        AmCH = (TH1D*)OriginalAmCH->Clone();
        
        for(Int_t i=1;i<=AmCH->GetXaxis()->GetNbins();i++)
        {
            if((AmCH->GetBinContent(i))<0)
            {
                AmCH->SetBinContent(i,-1*(AmCH->GetBinContent(i))); //Invert negative bins??
                //                          OR
                // AmCH->SetBinContent(i,0);//Set negative bins to 0?
            }
        }
    }
    
    delete OriginalAmCH;
    
    BackgroundsH[NADs+2]=(TH1D*)AmCH->Clone();
    BackgroundsH[NADs+2]->Reset();
    
    //    I need to study what are the implications. Right now I'm not applying any correction to the original histogram.
    
    Double_t IntegralOriginal = AmCH->Integral();
    //    cout << IntegralOriginal<<"\n";
    
    // FIT from 1.5-12 (original limits)
    HAmCFunc = new TF1("AmCFunc","exp([0]+x*[1])",InitialVisibleEnergy,FinalVisibleEnergy);
    AmCH->Fit("AmCFunc","Q0");//"Q0 is quiet mode and no draw
    
    //Extend background to 0-12
    TF1* ExtendedHAmCFunc = new TF1("AmCFunc","exp([0]+x*[1])",0,FinalVisibleEnergy);
    ExtendedHAmCFunc->SetParameter(0,HAmCFunc->GetParameter(0));
    ExtendedHAmCFunc->SetParameter(1,HAmCFunc->GetParameter(1));
    
    ///////////////////////////////////////////////////
    //Should I fit this with or without errors? (CHECK) (Original histogram (errors+content) vs only GetBinContent)???? What about negative bins?
    ///////////////////////////////////////////////////
    
    BackgroundsH[NADs+2]->Add(ExtendedHAmCFunc);
    
    delete ExtendedHAmCFunc;
    
    Double_t IntegralExponential = BackgroundsH[NADs+2]->Integral();//  Check that this works correctly, otherwise write explicitly limits in terms of the bins (InitialVisibleEnergy/binwidth, Final...)
    //    cout << IntegralExponential <<"\n";
    
    //Normalize histogram with respect to integral of the original graph (not accurate since the limits are different, anyway when backgrounds are used they are normalized using data to match the expected number of AmC events)
    
    BackgroundsH[NADs+2]->Scale(IntegralOriginal/IntegralExponential);
    // BackgroundsH[NADs]->Draw("same");
    //cout << IntegralOriginal/IntegralExponential <<"\n";
    delete AmCH;
}

void FitBackgrounds2::PrintBackgrounds()
{
        TCanvas* AccidentalC = new TCanvas("AccidentalC","Accidental Backgrounds", NADs*400/2,400*2);
        for (Int_t week = 0; week<Nweeks; week++)
        {
            AccidentalC->Clear();
            AccidentalC->Divide(NADs/2,2*XCellLimit*YCellLimit);
            
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    for (Int_t AD = 0; AD<NADs; AD++)
                    {
                        AccidentalC->cd(AD+NADs*idy+NADs*YCellLimit*idx+1);
                        FinalAccidentalsH[AD][week][idx][idy]->SetStats(0);
                        FinalAccidentalsH[AD][week][idx][idy]->SetTitle(Form("Accidentals AD%i Volume%i",AD+1,idx));
                        if(!LoganBinning)
                        {
                            FinalAccidentalsH[AD][week][idx][idy]->GetXaxis()->SetRange(1,n_evis_bins-1);
                        }
                        FinalAccidentalsH[AD][week][idx][idy]->Draw("HIST");
                        AccidentalC->Modified();
                    }
                }
            }
            AccidentalC->Print(("./Images/"+AnalysisString+Form("/BackgroundVariations/AccidentalBackgroundsPeriod%i.eps",week)).c_str(),".eps");
        }
        delete AccidentalC;
        
        TCanvas* LiHeC = new TCanvas("LiHeC","LiHe Backgrounds", NADs*400/2,400*2);
        for (Int_t week = 0; week<Nweeks; week++)
        {
            LiHeC->Clear();
            LiHeC->Divide(NADs/2,2*XCellLimit*YCellLimit);
            
          
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    for (Int_t AD = 0; AD<NADs; AD++)
                    {
                        LiHeC->cd(AD+NADs*idy+NADs*YCellLimit*idx+1);
                        FinalLiHeH[AD][week][idx][idy]->SetStats(0);
                        FinalLiHeH[AD][week][idx][idy]->SetTitle(Form("Li He AD%i Volume%i",AD+1,idx));
                        if(!LoganBinning)
                        {
                            FinalLiHeH[AD][week][idx][idy]->GetXaxis()->SetRange(1,n_evis_bins-1);
                        }
                        FinalLiHeH[AD][week][idx][idy]->Draw("HIST");
                        LiHeC->Modified();
                    }
                }
            }
            LiHeC->Print(("./Images/"+AnalysisString+Form("/BackgroundVariations/LiHeBackgroundsPeriod%i.eps",week)).c_str(),".eps");
        }
        delete LiHeC;
        
        TCanvas* FastC = new TCanvas("FastC","Fast Backgrounds", NADs*400/2,400*2);
        for (Int_t week = 0; week<Nweeks; week++)
        {
            FastC->Clear();
            FastC->Divide(NADs/2,2*XCellLimit*YCellLimit);
            
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    for (Int_t AD = 0; AD<NADs; AD++)
                    {
                        FastC->cd(AD+NADs*idy+NADs*YCellLimit*idx+1);
                        FinalFastNeutronsH[AD][week][idx][idy]->SetStats(0);
                        FinalFastNeutronsH[AD][week][idx][idy]->SetTitle(Form("Fast Neutrons AD%i Volume%i",AD+1,idx));
                        if(!LoganBinning)
                        {
                            FinalFastNeutronsH[AD][week][idx][idy]->GetXaxis()->SetRange(1,n_evis_bins-1);
                        }
                        FinalFastNeutronsH[AD][week][idx][idy]->Draw("HIST");
                        FastC->Modified();
                    }
                }
            }
            FastC->Print(("./Images/"+AnalysisString+Form("/BackgroundVariations/FastNeutronBackgroundsPeriod%i.eps",week)).c_str(),".eps");
        }
        delete FastC;
        
        TCanvas* AmCC = new TCanvas("AmCC","AmC Backgrounds", NADs*400/2,400*2);
        for (Int_t week = 0; week<Nweeks; week++)
        {
            AmCC->Clear();
            AmCC->Divide(NADs/2,2*XCellLimit*YCellLimit);
            
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    for (Int_t AD = 0; AD<NADs; AD++)
                    {
                        AmCC->cd(AD+NADs*idy+NADs*YCellLimit*idx+1);
                        FinalAmCH[AD][week][idx][idy]->SetStats(0);
                        FinalAmCH[AD][week][idx][idy]->SetTitle(Form("AmC AD%i Volume%i",AD+1,idx));
                        if(!LoganBinning)
                        {
                            FinalAmCH[AD][week][idx][idy]->GetXaxis()->SetRange(1,n_evis_bins-1);
                        }
                        FinalAmCH[AD][week][idx][idy]->Draw("HIST");
                        AmCC->Modified();
                    }
                }
            }
            AmCC->Print(("./Images/"+AnalysisString+Form("/BackgroundVariations/AmCBackgroundsPeriod%i.eps",week)).c_str(),".eps");
        }
        delete AmCC;

        TCanvas* NearBackgroundC = new TCanvas("NearBackgroundC","Near Backgrounds", (ADsEH1+ADsEH2)*400,400);
        TCanvas* FarBackgroundC = new TCanvas("FarBackgroundC","Far Backgrounds", (ADsEH3)*400,400);
        
        NearBackgroundC->Divide(ADsEH1+ADsEH2,XCellLimit*YCellLimit);
        FarBackgroundC->Divide(ADsEH3,XCellLimit*YCellLimit);
        
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for(Int_t idx=0; idx<XCellLimit; idx++)
            {
                for(Int_t idy=0; idy<YCellLimit; idy++)
                {
                    for (Int_t AD = 0; AD<NADs; AD++)
                    {
                        if(AD<ADsEH1+ADsEH2)
                        {
                            NearBackgroundC->cd(AD+NADs*idy+NADs*YCellLimit*idx+1);
                            NearBackgroundSpectrumH[AD][week][idx][idy]->Draw("HIST");
                            NearBackgroundC->Modified();
                        }
                        else
                        {
                            FarBackgroundC->cd(AD+NADs*idy+NADs*YCellLimit*idx+1);
                            FarBackgroundSpectrumH[AD][week][idx][idy]->Draw("HIST");
                            FarBackgroundC->Modified();
                        }
                    }
                }
            }
            NearBackgroundC->Print(("./Images/"+AnalysisString+Form("/BackgroundVariations/NearBackgroundsPeriod%i.eps",week)).c_str(),".eps");
            FarBackgroundC->Print(("./Images/"+AnalysisString+Form("/BackgroundVariations/FarBackgroundsPeriod%i.eps",week)).c_str(),".eps");
        }
        
        delete NearBackgroundC;
        delete FarBackgroundC;
    
}

void FitBackgrounds2::ProduceNominalBackgroundHistograms()
{
    if(isH)
    {
         ReadHBackgrounds();
    }
    else
    {
         ReadGdBackgrounds();
    }
    
    SaveAndDelete();
    
}
