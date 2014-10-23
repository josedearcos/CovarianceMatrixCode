#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TH2.h"

#include "TMath.h"

#include "NominalData.h"
#include <math.h>

const bool PreRebin = 1; //0 to rebin the Energy Matrix after its calculation, 1 to rebin the IAV matrix before the Energy Matrix calculation.
const Int_t MaxPeriods = 20;
const Int_t MaxDetectors = 8;
const Int_t MaxNearDetectors =4;
const Int_t MaxFarDetectors =4;
const Int_t MaxNbins=41;
const Int_t Halls=3;
const Int_t MatrixBins=240;//Same than IAV matrix for Gd. When the IAV H matrix is produced this can be selected through NominalData.h
const Int_t n_bcw_positron_nl = 1000;

const Double_t Me = 0.510999; // MeV
const Double_t Mn = 939.565; // MeV
const Double_t Mp = 938.272; // MeV

Double_t VisibleEnergy0F(Double_t * x, Double_t * par);// Zeroth order true to visible function.
Double_t VisibleEnergy1F(Double_t * x, Double_t * par);// First order true to visible function.

Double_t enu_to_epositron(Double_t * x, Double_t * par);
Double_t reso_func_bcw(Double_t * x, Double_t * par);

class EnergyMatrix1
{
    private:
    
    TRandom3* rand;
    
    enum Systematic{ReactorE, EnergyE, IAVE, NLE, ResolutionE};
    Systematic SystematicE;
    
    //AD configuration parameters:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    Int_t hall;
    
    TH1F* PredictionH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    TH2F* TrueEnergyToVisibleEnergy;

    //Binning parameters:
    Int_t Nbins;
    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    Int_t Nweeks;
    Int_t InitialBins;
    Int_t RebinFactor;
    Double_t BinWidth;
    Int_t TotalBins;
    
    char* OutputFileName;
    
    Double_t m_energy[MaxDetectors][MatrixBins];
    Double_t m_positronTrueSpectrum[MaxDetectors][MatrixBins];
    Double_t m_positronDetectedSpectrum[MaxDetectors][MatrixBins];
    
    TH1F* PositronTrueSpectrumH[MaxDetectors][MatrixBins];
    TH1F* PositronIAVSpectrumH[MaxDetectors][MatrixBins];
    TH1F* PositronNLSpectrumH[MaxDetectors][MatrixBins];
    TH1F* PositronVisibleSpectrumH[MaxDetectors][MatrixBins];
    TH1D* EnergyRebinH[MaxDetectors];
    TH1D* ProjectionH[MaxDetectors];
    
    TH2F* EnergyM[MaxDetectors];
    TH2F* CopyEnergyM[MaxDetectors];
    TH2* h_evis_vs_enu_ad[MaxDetectors];

    //Non linearity
    Double_t m_bcw_elec_nl_par[5];
    Double_t m_bcw_elec_nl_par_nominal[5];
    Double_t m_bcw_elec_nl_par_error[5];
    
    Double_t m_rel_escale[MaxDetectors];
    Double_t m_rel_escale_nominal[MaxDetectors];
    Double_t m_rel_escale_error[MaxDetectors];
    Double_t m_rel_eoffset[MaxDetectors];
    
    //IAV
    Double_t IAVError; // relative uncertainty of the IAV thickness
    Double_t IAVError[MaxDetectors]; // relative uncertainty of the IAV thickness
    Double_t IAVMatrix[MatrixBins][MatrixBins];
    
    Double_t m_detectorResolution; // Detector energy resolution parameter
    Double_t m_detectorResolution_nominal; // Detector energy resolution parameter
    Double_t m_detectorResolution_error; // Detector energy resolution parameter
    Double_t m_detectorResolution_error_uncorr; // Detector energy resolution parameter
    Double_t m_detectorResolution_bias[MaxDetectors]; // Random biases of the detector resolution.
    
    Double_t m_positronNLSpectrum[MaxDetectors][MatrixBins];
    
    public:
    
    EnergyMatrix1();
    EnergyMatrix1(NominalData Data);
    void EnergyMatrixMain();
    TF1* GetNeutrinoToVisibleFunction();
    void LoadPredictions();
    void GetPositronEnergy();
    void LoadIavCorrection();
    void CreateEnergyMatrix();
    
    void SetReactorMatrix(bool ReactorMatrix);
    void SetEnergyScaleMatrix(bool EnergyMatrix);
    void SetIAVMatrix(bool IAVMatrix);
    void SetNLMatrix(bool NLMatrix);
    void SetResolutionMatrix(bool ResolutionMatrix);

    void RandomReactorMatrix();
    void RandomEnergyScaleMatrix();
    void RandomIAVMatrix();
    void RandomNLMatrix();
    void RandomResolutionMatrix();
    TH2F* GetEnergyMatrix(Int_t Detector);
};

EnergyMatrix1::EnergyMatrix1()
{
    Nom = new NominalData();
    rand = new TRandom3();
    
    Nbins = Nom->GetNbins();
    InitialEnergy = Nom->GetEmin();
    FinalEnergy = Nom->GetEmax();
    InitialVisibleEnergy = Nom->GetEVisMin();
    FinalVisibleEnergy = Nom->GetEVisMax();
    InitialBins = Nom->GetInitialBins();
    RebinFactor = 240/(Nbins + InitialBins);
    cout << "NBINS IS " << Nbins << endl;

    cout << "ORIGINAL REBIN FACTOR IS " << RebinFactor << endl;
    Nweeks = Nom->GetWeeks();
    
    NADs = Nom->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    TotalBins = MatrixBins;

    if(PreRebin)
    {
        TotalBins = (Nbins + InitialBins);
    }
    
    for(Int_t idet = 0; idet< ADsEH1+ADsEH2; idet++)
    {
        CopyEnergyM[idet] = new TH2F(Form("EvisEnu%d",idet),Form("EvisEnu%d",idet),TotalBins,InitialVisibleEnergy,FinalVisibleEnergy,TotalBins,InitialVisibleEnergy,FinalVisibleEnergy);
    }
    
    BinWidth=(FinalVisibleEnergy-InitialVisibleEnergy)/TotalBins;
    
    if(NADs == 8) //ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    IAVError=Nom->GetIAVError();
    
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

EnergyMatrix1 :: EnergyMatrix1 (NominalData Data)
{
    rand = new TRandom3();

    Nbins = Data.GetNbins();
    InitialEnergy = Data.GetEmin();
    FinalEnergy = Data.GetEmax();
    InitialVisibleEnergy =Data.GetEVisMin();
    FinalVisibleEnergy = Data.GetEVisMax();
    InitialBins = Data.GetInitialBins();
    RebinFactor = MatrixBins/(Nbins + InitialBins);

    Nweeks = Data.GetWeeks();
    
    NADs = Data.GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    TotalBins = MatrixBins;

    if(PreRebin)
    {
        TotalBins = (Nbins + InitialBins);
    }
    
    for(Int_t idet = 0; idet< ADsEH1+ADsEH2; idet++)
    {
        CopyEnergyM[idet] = new TH2F(Form("EvisEnu%d",idet),Form("EvisEnu%d",idet),TotalBins,InitialVisibleEnergy,FinalVisibleEnergy,TotalBins,InitialVisibleEnergy,FinalVisibleEnergy);
    }
    
    BinWidth=(FinalVisibleEnergy-InitialVisibleEnergy)/TotalBins;
    
    if(NADs == 8) //ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    IAVError=Data.GetIAVError();
    
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

void EnergyMatrix1 :: EnergyMatrixMain()
{
      LoadPredictions();
      LoadIavCorrection();
      cout << "Systematic" << SystematicE << endl;
      switch (SystematicE)
      {
        case 1://Vary Reactor

        case 2://Vary Energy Scale

        case 3://Vary IAV
              RandomIAVMatrix();
        case 4://Vary NL
              
        case 5://Vary Resolution

        default://Add nominal systematics
   
      }
    
      GetPositronEnergy();
}

void EnergyMatrix1 :: RandomReactorMatrix()
{

}
void EnergyMatrix1 :: RandomEnergyScaleMatrix()
{
    
}
void EnergyMatrix1 :: RandomIAVMatrix()
{
    for (Int_t AD=0; AD<NADs; AD++)
    {
        IAVError[AD]= (1+IAVError*rand->Gaus(0,1)); //Each AD is varied individually.
    }

    for(Int_t i=0; i<TotalBins; i++)
    {
        for(Int_t j=0; i<TotalBins; i++)
        {
            cout<<"IAV ORIGINAL" << IAVMatrix[i][j] << endl;

            IAVMatrix[i][j]=IAVError[AD]*IAVMatrix[i][j];
            
            cout<< "IAV ALTERED" << IAVMatrix[i][j] << endl;
        }
    }
}
void EnergyMatrix1 :: RandomNLMatrix()
{
    
}
void EnergyMatrix1 :: RandomResolutionMatrix()
{
    
}
void EnergyMatrix1 :: CreateEnergyMatrix()
{
    TFile* SaveSpectrumDataF = TFile::Open("./RootOutputs/EnergyMatrix1.root","update");

    // This method should be used with the relative oscillation method
    for(int idet=0; idet<(3); idet++)//(ADsEH1+ADsEH2)
    {
        EnergyM[idet] = new TH2F(Form("EvisEnu%d",idet),Form("EvisEnu%d",idet),TotalBins,InitialVisibleEnergy,FinalVisibleEnergy,TotalBins,InitialVisibleEnergy,FinalVisibleEnergy);
        EnergyM[idet]->SetDirectory(0);
        for(int TrueEnergyIndex=1; TrueEnergyIndex<=TotalBins; TrueEnergyIndex++)
        {
            for(int VisibleBin=1;VisibleBin<=TotalBins;VisibleBin++)
            {
                EnergyM[idet]->SetBinContent(TrueEnergyIndex,VisibleBin,PositronVisibleSpectrumH[idet][TrueEnergyIndex-1]->GetBinContent(VisibleBin));
            }
        }
        
        cout << " HALLELUYA" << endl;
        if(!PreRebin)
        {
            cout << "REBIN FACTOR IS " << RebinFactor;
            EnergyM[idet] = (TH2F*)EnergyM[idet]->Rebin2D(RebinFactor,RebinFactor);
        }
        
        CopyEnergyM[idet] = (TH2F*)EnergyM[idet]->Clone();
        CopyEnergyM[idet]->SetDirectory(0);
        cout << " HALLELUYA" << endl;

        EnergyM[idet]->Write();
        cout << " HALLELUYA" << endl;

            delete EnergyM[idet];
            cout << " HALLELUYA" << endl;
    }
    SaveSpectrumDataF->Close();
}

void EnergyMatrix1 :: GetPositronEnergy()
{
    Double_t m_abs_escale =1;
    Double_t m_abs_escale_nominal =1;
    Double_t m_abs_escale_error = 0.01;
    
    Double_t m_abs_eoffset = 0.0;
    Double_t m_abs_eoffset_error = 0.08;  //(MeV)
    
    Double_t m_rel_eoffset_error = 0.013; //  (MeV)

    //Input data
    ifstream bcw_elec_data("bcw_nl_data/par.dat");
    for (Int_t i = 0; i < 3; i++)
    {
        bcw_elec_data >> m_bcw_elec_nl_par_nominal[i]  >> m_bcw_elec_nl_par_error[i];
        m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];
    }
    bcw_elec_data.close();

    //Energy shift function set up
    VisibleF=GetNeutrinoToVisibleFunction(0);
    
    //NL function set up
    nl_func = new TF1("nl_func",nl_func_bcw,0,12,7);
    
    //Resolution function set up
    reso_func = new TF1("reso_func",reso_func_bcw,0,20,3);
    reso_func->SetParameters(0.022,0.077,0.018); // based on Bryce's TN
    
    //BCW values http://dayabay.ihep.ac.cn/DocDB/0087/008768/013/6AdAnalysis-BCW.pdf are different! 0.13226034931, 0.32604048828, 0.26196488314
   
    for(int ii=0;ii<NADs;ii++)
    {
        m_detectorResolution_bias[ii] = 0;//this will be used to distort it
    }   

    //Energy shift
    for(int idet=0; idet<1; idet++)
    {
        //I define this parameters here to use the for loop only once
        m_rel_escale[idet] = 1.0;
        m_rel_escale_error[idet] = 0.0035; // 0.35%
        m_rel_escale_nominal[idet] = m_rel_escale[idet];
        m_rel_eoffset[idet] = 0.0;
        
        nl_func->SetParameters(m_bcw_elec_nl_par[0],m_bcw_elec_nl_par[1],m_bcw_elec_nl_par[2], m_bcw_elec_nl_par[3], m_bcw_elec_nl_par[4], m_abs_escale * m_rel_escale[idet], m_abs_eoffset+m_rel_eoffset[idet]);
        
        //Set up the histograms
        for(int TrueEnergyIndex=0; TrueEnergyIndex<TotalBins; TrueEnergyIndex++)
        {
            PositronTrueSpectrumH[idet][TrueEnergyIndex]= new TH1F(Form("PositronTrueSpectrum%d,%d",idet,TrueEnergyIndex),Form("PositronTrueSpectrum%d,%d",idet,TrueEnergyIndex), TotalBins,InitialVisibleEnergy,FinalVisibleEnergy);
            PositronTrueSpectrumH[idet][TrueEnergyIndex]->SetDirectory(0);
            PositronTrueSpectrumH[idet][TrueEnergyIndex]->Reset();

            //IAV
            PositronIAVSpectrumH[idet][TrueEnergyIndex]=(TH1F*)PositronTrueSpectrumH[idet][TrueEnergyIndex]->Clone(Form("IAVSpectrum%d,%d",idet,TrueEnergyIndex));
            PositronIAVSpectrumH[idet][TrueEnergyIndex]->SetDirectory(0);
            PositronIAVSpectrumH[idet][TrueEnergyIndex]->SetTitle("IAV Spectrum");
            PositronIAVSpectrumH[idet][TrueEnergyIndex]->Reset();
            //NL
            PositronNLSpectrumH[idet][TrueEnergyIndex] = (TH1F*)PositronIAVSpectrumH[idet][TrueEnergyIndex]->Clone(Form("NLSpectrum%d,%d",idet,TrueEnergyIndex));
            PositronNLSpectrumH[idet][TrueEnergyIndex]->SetDirectory(0);
            PositronNLSpectrumH[idet][TrueEnergyIndex]->SetTitle("NL Spectrum");
            PositronNLSpectrumH[idet][TrueEnergyIndex]->Reset();
            //Resolution
            PositronVisibleSpectrumH[idet][TrueEnergyIndex] = (TH1F*)PositronNLSpectrumH[idet][TrueEnergyIndex]->Clone(Form("VisibleSpectrum%d,%d",idet,TrueEnergyIndex));
            PositronVisibleSpectrumH[idet][TrueEnergyIndex]->SetDirectory(0);
            PositronVisibleSpectrumH[idet][TrueEnergyIndex]->SetTitle("Visible Spectrum");
            PositronVisibleSpectrumH[idet][TrueEnergyIndex]->Reset();
        }
        //Calculate Energy Shift
        for(Int_t TrueEnergyIndex=0; TrueEnergyIndex<TotalBins; TrueEnergyIndex++)
        {
            Double_t e_neutrino = (TrueEnergyIndex)*BinWidth;
            Double_t e_positron = VisibleF->Eval(e_neutrino);
            Int_t e_positron_index=(Int_t)(e_positron/BinWidth);
  
            if (e_positron < 1.015)//1.8-0.78 is about the minimum visible energy (1.017999 for 0th order, 1.01574 for 1st order calculation)
            {
                PositronTrueSpectrumH[idet][TrueEnergyIndex]->Reset();
            }
            else
            {
                PositronTrueSpectrumH[idet][TrueEnergyIndex]->SetBinContent(e_positron_index+1, PredictionH[0][0][0]->GetBinContent(e_positron_index+1)); //Change this to include all ADs.
            }
        //Calculate IAV Shift
            for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TrueEnergyIndex+1; VisibleEnergyIndex++)
            {
                if (e_positron >= 1.015)//1.8-0.78 is about the minimum visible energy (1.017999 for 0th order, 1.01574 for 1st order calculation)
                {
                    PositronIAVSpectrumH[idet][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, PositronIAVSpectrumH[idet][TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1)+ IAVMatrix[e_positron_index][VisibleEnergyIndex] * PositronTrueSpectrumH[idet][TrueEnergyIndex]->GetBinContent(e_positron_index+1));
                    
//                    cout << e_positron_index << " is the e_positron_index" << endl;
//                    cout << VisibleEnergyIndex << " is the VisibleEnergyIndex" << endl;
                }
            }
        //Calculate NL Shift
            for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
            {
                PositronNLSpectrumH[idet][TrueEnergyIndex]->SetBinContent(nl_func->Eval(VisibleEnergyIndex)+1,PositronIAVSpectrumH[idet][TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1));
             //   cout << "positron NL index" << nl_func->Eval(VisibleEnergyIndex)+1 << "with values" << PositronIAVSpectrumH[idet][TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1) << "and true energy" << TrueEnergyIndex << endl;
                //weird behaviour at TrueEnergy = 66; 117 and 167
            }
        //Calculate Resolution effect
            for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
            {
                //Resolution
                Double_t resolutionRange = 8; // sigma // Why 8?
                
                Double_t sigma = (reso_func->Eval(VisibleEnergyIndex*BinWidth) + m_detectorResolution_bias[idet]) * VisibleEnergyIndex*BinWidth;
//                cout << "SIGMA " << sigma<<endl;
                Double_t minDetE = VisibleEnergyIndex*BinWidth - resolutionRange*sigma;
                Double_t maxDetE = VisibleEnergyIndex*BinWidth + resolutionRange*sigma;
                Int_t minDetEIdx = (Int_t)((minDetE)/BinWidth);
                Int_t maxDetEIdx = (Int_t)((maxDetE)/BinWidth);
                
//                cout << "minDetE "<<minDetE<<endl;
//                cout << "maxDetE"<<maxDetE<<endl;
//                cout << "minDetEIdx "<<minDetEIdx<<endl;
//                cout << "maxDetEIdx"<<maxDetEIdx<<endl;
                
                if(minDetEIdx < 0)
                {
                    minDetEIdx = 0;
                }
                if(maxDetEIdx >= TotalBins)
                {
                    maxDetEIdx = TotalBins-1;
                }
                for(Int_t detIdx=minDetEIdx; detIdx<=maxDetEIdx; detIdx++)
                {
                    if(detIdx==0)
                    {
                        continue;
                    }
                    Double_t gausFactor = TMath::Gaus((VisibleEnergyIndex-detIdx)*BinWidth,0,sigma,true);
                    
//                    cout << "detidx " <<detIdx << endl;
//                    cout << "VisibleEnergyIndex " << VisibleEnergyIndex <<endl;
//                    cout << "X gaussian "<< ((VisibleEnergyIndex-detIdx)/BinWidth) << endl;
//                    cout << "gaus factor "<<gausFactor<<endl;

                    PositronVisibleSpectrumH[idet][TrueEnergyIndex]->SetBinContent(detIdx+1, PositronVisibleSpectrumH[idet][TrueEnergyIndex]->GetBinContent(detIdx+1)+(PositronNLSpectrumH[idet][TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1)*gausFactor));
                }
            }
        }
    }
    delete VisibleF;
    delete reso_func;
    delete nl_func;
                                                    
// UNCOMMENT FOLLOWING LINES TO SAVE ENERGY MATRIX SLICES IN EACH PRODUCTION STEP
//    for(int idet=0; idet<1; idet++)
//    {
//        for(int TrueEnergyIndex=0; TrueEnergyIndex<TotalBins; TrueEnergyIndex++)
//        {
//            PositronTrueSpectrumH[idet][TrueEnergyIndex]->Write();
//            PositronIAVSpectrumH[idet][TrueEnergyIndex]->Write();
//            PositronNLSpectrumH[idet][TrueEnergyIndex]->Write();
//            PositronVisibleSpectrumH[idet][TrueEnergyIndex]->Write();
//        }
//    }
    CreateEnergyMatrix();
}
                                                    
TF1* EnergyMatrix1 :: GetNeutrinoToVisibleFunction(Int_t order)
{
    if(order==0)
    {
        Double_t Correction = Mn-Mp-Me;
        TF1 *VisibleF = new TF1("VisibleF",VisibleEnergy0F,InitialEnergy,FinalEnergy,1);
        VisibleF->SetParameter(0,Correction);
        cout<<"Zeroth order"<<endl;
    }
    if(order==1)
    {
        TF1 *VisibleF = new TF1("VisibleF",VisibleEnergy1F,InitialEnergy,FinalEnergy,1);
        VisibleF->SetParameter(0,0);//Maybe I can improve this using a MC simulation of the angular distribution calculated in http://authors.library.caltech.edu/2796/1/VOGprd99.pdf
        cout<<"First order"<<endl;
    }
    return VisibleF;
}

// First order true to visible function.
//(x-(Mn-Mp))*(1 - x/Mn*(1.0 - [0]*sqrt(1 - Me*Me/(x-(Mn-Mp))/(x-(Mn-Mp)))))- ((Mn-Mp)*(Mn-Mp) - Me*Me)/(2*Mn)-Me]
Double_t VisibleEnergy1F(Double_t * x, Double_t * par)
{
    Double_t Enu = x[0];
    Double_t costheta = par[0];

    Double_t Delta = Mn-Mp;
    
    Double_t Ee0 = Enu - Delta;
    
    Double_t gamma0 = Ee0/Me;
    Double_t beta0 = sqrt(1 - 1/gamma0/gamma0);
    Double_t y2 = (Delta*Delta - Me*Me)/2.;
    
    Double_t Ee1 = Ee0 * (1 - Enu/Mn*(1.0 - costheta*beta0)) - y2/Mn;
    
    return Ee1 + Me;
}

Double_t VisibleEnergy0F(Double_t * x, Double_t * par)
{
    return x[0]-par[0];
}

void EnergyMatrix1 :: LoadPredictions()
{
    TFile* FarHallPredictionsF = TFile::Open("./RootOutputs/FarSpectrumFraction240.root");//To match Gd IAV size.
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            for (Int_t far =0; far<ADsEH3; far++)
            {
                PredictionH[far][near][week] = (TH1F*)gDirectory->Get(Form("AD%i Far Spectrum prediction from near AD%i",far+1,near+1));
                PredictionH[far][near][week]->SetDirectory(0);
                if(PreRebin)
                {
                    PredictionH[far][near][week]->Rebin(RebinFactor);
                }
            }
        }
    }
    

    FarHallPredictionsF->Close();
}

Double_t enu_to_epositron(Double_t * x, Double_t * par)
{
    Double_t Enu = x[0];
    Double_t costheta = par[0];
    
    Double_t Me = 0.510999; // MeV
    Double_t Mn = 939.565; // MeV
    Double_t Mp = 938.272; // MeV
    Double_t Delta = Mn-Mp;
    
    Double_t Ee0 = Enu - Delta;
    // if (Ee0 < Me){    // not allowed!
    //   return -1;
    // }
    
    Double_t gamma0 = Ee0/Me;
    Double_t beta0 = sqrt(1 - 1/gamma0/gamma0);
    Double_t y2 = (Delta*Delta - Me*Me)/2.;
    
    Double_t Ee1 = Ee0 * (1 - Enu/Mn*(1.0 - costheta*beta0)) - y2/Mn;

    return Ee1 + Me; // Add Me to include annihiration gamma energy
}

void EnergyMatrix1 :: LoadIavCorrection()//From Bryce Littlejohn's results. //I don't use different IAV matrix for each AD since the analysis it's been done for just 1 of them, assume identical ADs.
{
    
    TFile * f = new TFile("./IavDistortion/IAVDistortion.root");
    TH2F * Correction = (TH2F*)f->Get("Correction");
    
    cout << "reading Iav correction file" << endl;
    
    // I can rebin here, it decreases the calculation time of covariance matrices. But before doing so I have to check whether rebinning before the full calculation affects the final result or not.

    if(PreRebin)
    {
        Correction->Rebin2D(RebinFactor,RebinFactor);
    }
    for(Int_t i=0;i<TotalBins;i++)
    { // i: true positron energy bin; j: distorted energy bin
        //Total events in input spectrum bin i
        if (Correction->Integral(i+1,i+1,0,-1) > 0)
        {
            for(Int_t j=0;j<i+1;j++)
            {
                IAVMatrix[i][j]= Correction->GetBinContent(i+1,j+1)/Correction->Integral(i+1,i+1,0,-1);
            }
        }
        else
        {
            for(Int_t j=0;j<i+1;j++)
            {
                if (i==j)
                {
                    IAVMatrix[i][j] = 1;
                }
                else
                {
                    IAVMatrix[i][j] = 0;
                }
            }
        }
    }
    

    f->Close();
}


Double_t reso_func_bcw(Double_t * x, Double_t * par)
{
    
    Double_t e_orig = x[0];
    Double_t e_sigma = 1.0;
    
    if (e_orig > 0)
    {
        e_sigma = TMath::Sqrt(par[0]*par[0] + par[1]*par[1]/e_orig + par[2]*par[2]/e_orig/e_orig);
    }
    
    return e_sigma;
}

Double_t nl_func_bcw(Double_t * x, Double_t * par)
{
    //Input error
    TGraph* g_bcw_elec_nl_error[2];

    TFile *bcw_ele_err_file = new TFile("bcw_nl_data/ele_err.root");
    g_bcw_elec_nl_error[0] = (TGraph*)bcw_ele_err_file->Get("g_up")->Clone();
    g_bcw_elec_nl_error[1] = (TGraph*)bcw_ele_err_file->Get("g_down")->Clone();
    bcw_ele_err_file->Close();
    
    //BCW NL model (Apr. 1, 2013)
    Double_t m_bcw_positron_nl_e[n_bcw_positron_nl];
    Double_t m_bcw_positron_nl_fac[n_bcw_positron_nl];
    Double_t m_bcw_elec_nl_par[5];
    Double_t m_bcw_elec_nl_par_nominal[5];
    Double_t m_bcw_elec_nl_par_error[5];
    
    ifstream bcw_positron_data("bcw_nl_data/positron.dat");
    
    for (Int_t i = 0; i < n_bcw_positron_nl; i++)
    {
        bcw_positron_data >> m_bcw_positron_nl_e[i] >> m_bcw_positron_nl_fac[i];
    }
    
    bcw_positron_data.close();
    
    Double_t e_positron_true = x[0];
    Double_t escale_par = par[5]; // Add flat energy scale parameter
    Double_t escale_offset = par[6]; // Add fixed energy offset
    
    Double_t scinti_nl_fac = 1;
    
    if (e_positron_true < m_bcw_positron_nl_e[0])
    {
        scinti_nl_fac = m_bcw_positron_nl_fac[0];
    }
    else if (e_positron_true > m_bcw_positron_nl_e[n_bcw_positron_nl-1])
    {
        scinti_nl_fac = m_bcw_positron_nl_fac[n_bcw_positron_nl-1];
    }
    else
    {
        for (Int_t i = 0; i < n_bcw_positron_nl-1; i++)
        {
            if (e_positron_true >= m_bcw_positron_nl_e[i] && e_positron_true < m_bcw_positron_nl_e[i+1])
            {
                scinti_nl_fac = ((m_bcw_positron_nl_e[i+1] - e_positron_true)*m_bcw_positron_nl_fac[i] + (e_positron_true - m_bcw_positron_nl_e[i])*m_bcw_positron_nl_fac[i+1]) / (m_bcw_positron_nl_e[i+1] - m_bcw_positron_nl_e[i]);
                break;
            }
        }
    }
    
    Double_t visibleE = scinti_nl_fac * e_positron_true;
    //  double electronicsCorrection = exp(par[0] + par[1] * visibleE) + exp(par[2] + par[3] * visibleE);
    
    Double_t err_offset = 0;
    Double_t err_shift = 0;
    
    Double_t elec_err_min_x = 3.45; // MeV
    Double_t par3_up = g_bcw_elec_nl_error[0]->Eval(elec_err_min_x) - 1;
    Double_t par3_down = g_bcw_elec_nl_error[1]->Eval(elec_err_min_x) - 1;
    
    if (par[3] > 0)
    {
        err_offset = par[3]*par3_up;
    }
    else
    {
        err_offset = par[3]*par3_down;
    }
    
    Double_t par4_up = 0;
    Double_t par4_down = 0;
    
    if (visibleE > elec_err_min_x)
    {
        par4_up = g_bcw_elec_nl_error[0]->Eval(visibleE)- g_bcw_elec_nl_error[0]->Eval(elec_err_min_x);
        par4_down = g_bcw_elec_nl_error[1]->Eval(visibleE)- g_bcw_elec_nl_error[1]->Eval(elec_err_min_x);
    }
    else
    {
        par4_up = g_bcw_elec_nl_error[1]->Eval(visibleE)- g_bcw_elec_nl_error[1]->Eval(elec_err_min_x);
        par4_down = g_bcw_elec_nl_error[0]->Eval(visibleE)- g_bcw_elec_nl_error[0]->Eval(elec_err_min_x);
    }
    if (par[4] > 0) err_shift = par[4]*par4_up;
    else err_shift = par[4]*par4_down;
    
    Double_t electronicsCorrection = exp(par[0] + par[1] * visibleE) + par[2] + err_offset + err_shift;
    
    Double_t final_energy =  visibleE * electronicsCorrection * escale_par + escale_offset;
    
    return final_energy;
}

void EnergyMatrix1 :: SetReactorMatrix(bool ReactorMatrix)
{
    if(ReactorMatrix)
    {
        RandomReactorMatrix();
        SystematicE = ReactorE;
    }
}

void EnergyMatrix1 :: SetEnergyScaleMatrix(bool EnergyScaleMatrix)
{
    if(EnergyScaleMatrix)
    {
        RandomEnergyScaleMatrix();
        SystematicE = EnergyE;
    }
}

void EnergyMatrix1 :: SetIAVMatrix(bool IAVMatrix)
{
    if(IAVMatrix)
    {
        RandomIAVMatrix();
        SystematicE = IAVE;
    }
}

void EnergyMatrix1 :: SetNLMatrix(bool NLMatrix)
{
    if(NLMatrix)
    {
        RandomNLMatrix();
        SystematicE = NLE;
    }
}

void EnergyMatrix1 :: SetResolutionMatrix(bool ResolutionMatrix)
{
    if(ResolutionMatrix)
    {
        RandomResolutionMatrix();
        SystematicE = ResolutionE;
    }
}

TH2F* EnergyMatrix1 :: GetEnergyMatrix(Int_t Detector)
{
    return CopyEnergyM[Detector];
}
