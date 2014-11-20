#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"

#include "NominalData.h"
#include <math.h>

const Int_t MaxPeriods = 20;
const Int_t MaxDetectors = 8;
const Int_t MaxNearDetectors =4;
const Int_t MaxFarDetectors =4;
const Int_t MaxNbins=41;
const Int_t Halls=3;
const Int_t MatrixBins=240;//Same than IAV matrix for Gd. When the IAV H matrix is produced this can be controlled through NominalData.h

const Double_t Me = 0.510999; // MeV
const Double_t Mn = 939.565; // MeV
const Double_t Mp = 938.272; // MeV
const double m_binWidth=(12.0-0)/MatrixBins;

Double_t VisibleEnergyF(Double_t * x, Double_t * par);// First order true to visible function.
Double_t enu_to_epositron(Double_t * x, Double_t * par);
Double_t reso_func_bcw(Double_t * x, Double_t * par);



class EnergyMatrix
{
    private:
    
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
    
    char* OutputFileName;
    
    Double_t m_energy[MaxDetectors][MatrixBins];
    Double_t m_positronTrueSpectrum[MaxDetectors][MatrixBins];
    Double_t m_positronDetectedSpectrum[MaxDetectors][MatrixBins];
    
    TH1F* PositronTrueSpectrumH[MaxDetectors];
    TH1F* PositronIAVSpectrumH[MaxDetectors];
    TH1F* PositronNLSpectrumH[MaxDetectors];
    TH1F* PositronDetectedSpectrumH[MaxDetectors];
    TH1D* htmp_rebin;
    TH1D* htmp_rebin1;

    TH2F* h_evis_vs_enu_ad[MaxDetectors];
    
    Double_t m_nominal_iav_corr[MatrixBins][MatrixBins];
    Double_t m_nominal_iav_frac[MatrixBins];
    
    Double_t m_iav_error; // relative uncertainty of the IAV thickness
    Double_t m_iav_corr[MaxDetectors][MatrixBins][MatrixBins];
    Double_t m_iav_frac[MaxDetectors][MatrixBins];
    
    Double_t m_positronIavDistortedSpectrum[MaxDetectors][MatrixBins];
    Double_t m_detectorResolution; // Detector energy resolution parameter
    Double_t m_detectorResolution_nominal; // Detector energy resolution parameter
    Double_t m_detectorResolution_error; // Detector energy resolution parameter
    Double_t m_detectorResolution_error_uncorr; // Detector energy resolution parameter
    Double_t m_detectorResolution_bias[MaxDetectors]; // Random biases of the detector resolution.
    
    Double_t m_positronNLSpectrum[MaxDetectors][MatrixBins];
    
    public:
    
    EnergyMatrix();
    EnergyMatrix(NominalData Data);
    void EnergyMatrixMain();
    void NeutrinoToVisible();
    TF1* GetNeutrinoToVisibleFunction();
    void LoadPredictions();
    void updatePositronTrue(Int_t index);
    void LoadIavCorrection();
    void SaveTotalSpectrum();

};

EnergyMatrix::EnergyMatrix()
{
    Nom = new NominalData();
    
    Nbins = Nom.GetNbins();
    InitialEnergy = Nom.GetEmin();
    FinalEnergy = Nom.GetEmax();
    InitialVisibleEnergy = Nom.GetEVisMin();
    FinalVisibleEnergy = Nom.GetEVisMax();
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
    
    BackgroundE=4;
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

EnergyMatrix :: EnergyMatrix (NominalData Data)
{
    Nbins = Data.GetNbins();
    InitialEnergy = Data.GetEmin();
    FinalEnergy = Data.GetEmax();
    InitialVisibleEnergy =Data.GetEVisMin();
    FinalVisibleEnergy = Data.GetEVisMax();
    
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

void EnergyMatrix :: EnergyMatrixMain()
{
      LoadPredictions();
      LoadIavCorrection();
      NeutrinoToVisible();
      SaveTotalSpectrum();
    
    for(Int_t idet =0;idet<NADs;idet++)
    {
        delete h_evis_vs_enu_ad[idet];
    }
//    EnergyMatrix2H=(TH2F*)Energy2H->Clone("Energy Matrix");
//    EnergyMatrix2H->Reset();
//    EnergyMatrix2H->SetDirectory(0);
//    
//    for (Int_t samples = 0; samples<NSamples; samples++)
//    {
//        FluctuateSystematics();
//        GenerateEnergyMatrix();
//        EnergyMatrix2H->Add(BkgCov2H);
//    }
//    
//    EnergyMatrix2H->Scale(1./(NSamples));
//
//    SaveCovarianceMatrix();
}

void EnergyMatrix::NeutrinoToVisible()
{
    for(int ii=0;ii<NADs;ii++){
        m_detectorResolution_bias[ii] = 0;
    }
//    Double_t VisibleE[MatrixBins];
//    GetNeutrinoToVisibleFunction(0);//Gets zeroth order (0) or 1st order (1) energy shift
//
//    TrueEnergyToVisibleEnergy = new TH2F ("TrueEnergyToVisibleEnergy","True Energy To Visible Energy",MatrixBins,InitialEnergy,FinalEnergy,MatrixBins,InitialVisibleEnergy,FinalVisibleEnergy);
//
//    for (Int_t i =0; i<MatrixBins; i++)
//    {
//        VisibleE[i] = VisibleF->Eval(InitialEnergy+(i+1)*((FinalEnergy-InitialEnergy)/MatrixBins));
//      //  cout<<VisibleF->GetX(1);
//      //  cout<<VisibleF->GetX(4);
//
//        //cout <<"\n"<< VisibleE[i]<<"at i" <<i+1;
//        //To check compatibility with LBNl results I should implement their energy array;
//        
//        for(int ivis=0;ivis<MatrixBins;++ivis)
//        {
//           TrueEnergyToVisibleEnergy->Fill(InitialEnergy+(i)*((FinalEnergy-InitialEnergy)/MatrixBins),InitialVisibleEnergy+(ivis)*((FinalVisibleEnergy-InitialVisibleEnergy)/MatrixBins),VisibleE[i]);
//        }
//        TrueEnergyToVisibleEnergy->Draw("colz");
    const Int_t n_evis_bins = 37;
    Double_t evis_bins[100]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
    evis_bins[0] = 0.7;
    for (Int_t i = 0; i < n_evis_bins-1; i++)
    {
        evis_bins[i+1] = 0.2 *i + 1.0;
    //    cout << evis_bins[i+1];
    }
    evis_bins[n_evis_bins] = 12.0;
    
    const Int_t n_enu_bins = 39;
    
    Double_t enu_bins[39+1]; // 39 bins between 1.8 and 9.6 MeV
    
    for (Int_t i = 0; i < 39+1; i++)
    {
        enu_bins[i] = 0.2 * i + 1.8;
    }
    
    for(int idet=0;idet<NADs;++idet)
    {
        h_evis_vs_enu_ad[idet] = new TH2F (Form("EvisEnu%d",idet),Form("EvisEnu%d",idet),MatrixBins,0,12,n_evis_bins,evis_bins);
    }
    
    for(int ibin_enu=0;ibin_enu<MatrixBins;++ibin_enu)
    {
        updatePositronTrue(ibin_enu);
        
        for(int ibin=0;ibin<MatrixBins;++ibin)
        {
            for(int idet=0;idet<NADs;++idet)
            {
                h_evis_vs_enu_ad[idet]->Reset();
                h_evis_vs_enu_ad[idet]->Fill(m_energy[idet][ibin_enu],m_energy[idet][ibin],m_positronDetectedSpectrum[idet][ibin]*m_binWidth);
                //  cout << m_positronTrueSpectrum[idet][ibin] << endl;
                h_evis_vs_enu_ad[idet]->SetDirectory(0);
            }
        }
    }
//    h_evis_vs_enu_ad[2]->Draw();

    //idet loop
    //    for(int idet=0;idet<NADs;++idet)
    //    {
    //    delete h_evis_vs_enu_ad[idet];
    //    }
    for (Int_t iEvisBin = 0; iEvisBin < n_evis_bins; iEvisBin++)
    {
        TH1D * htmp = h_evis_vs_enu_ad[0]->ProjectionX(Form("h_enu_%d",iEvisBin),iEvisBin+1,iEvisBin+1);
        
        htmp_rebin = (TH1D*)htmp->Rebin(n_enu_bins,Form("h_enu_rebin_%d",iEvisBin),enu_bins);
        Double_t norm = htmp_rebin->Integral();
        if(norm==0) norm =1;
        //cout << norm;
        htmp_rebin->Scale(1./norm);
        
        for (Int_t iEnuBin = 0; iEnuBin < n_enu_bins; iEnuBin++)
        {
          //  cout << htmp_rebin->GetBinContent(iEnuBin+1) << endl;
            
        }
//        cout << endl;
//    htmp_rebin->Draw("colz");
        // // for debugging
//    TCanvas * c1 = new TCanvas();
//    htmp_rebin->Draw();
//    c1->Update();
//    Char_t dummy[20];
//    cin >> dummy;
    }
    for (Int_t iEvisBin = 0; iEvisBin < n_evis_bins; iEvisBin++)
    {
        //    h_enu_tmp->Reset();
        TH1D * htmp1 = h_evis_vs_enu_ad[3]->ProjectionX(Form("h_enu_%d",iEvisBin));// add everything together
        
        htmp_rebin1 = (TH1D*)htmp1->Rebin(n_enu_bins,Form("h_enu_rebin_%d",iEvisBin),enu_bins);
        Double_t norm = htmp_rebin1->Integral();
        htmp_rebin1->Scale(1./norm);
     //   htmp_rebin1->Draw("colz");
    }
}

void EnergyMatrix::updatePositronTrue(Int_t index)
{
    enu_to_epositron_func = new TF1("enu_to_epositron_func",enu_to_epositron,0,20,1);
    enu_to_epositron_func->SetParameter(0,0);
    
    reso_func = new TF1("reso_func",reso_func_bcw,0,20,3);
    reso_func->SetParameters(0.022,0.077,0.018); // based on Bryce's TN

    for(int idet=0; idet<NADs; idet++)
    {
        PositronTrueSpectrumH[idet]= new TH1F(Form("PositronTrueSpectrum%d",idet),Form("PositronTrueSpectrum%d",idet), MatrixBins,0,12);
        PositronIAVSpectrumH[idet]= new TH1F(Form("PositronIAVSpectrum%d",idet),Form("PositronIAVSpectrum%d",idet), MatrixBins,0,12);
        PositronNLSpectrumH[idet]= new TH1F(Form("PositronNLSpectrum%d",idet),Form("PositronNLSpectrum%d",idet), MatrixBins,0,12);
        PositronDetectedSpectrumH[idet]= new TH1F(Form("PositronDetectedSpectrum%d",idet),Form("PositronDetectedSpectrum%d",idet), MatrixBins,0,12);
    }
    
    for(int idet=0; idet<NADs; idet++)
    {
        
        for(int idx=0; idx<MatrixBins; idx++)
        {
            m_energy[idet][idx] = ((0.5+idx)*m_binWidth)+0;//minimum energy shows as 0
            //                                  cout << "m_energy is" << (0.5+idx)*m_binWidth+0 << endl; //Checked
            
            double e_positron = m_energy[idet][idx];
            if (e_positron < 1.022)
            {
                // cout << e_positron << endl; Checked
                //                    cout << "idx" << idx << endl;
                m_positronTrueSpectrum[idet][idx] = 0;
                continue;
            }
            double e_nu = enu_to_epositron_func->GetX(e_positron);
            int e_nu_Idx = (int)((e_nu-0)/m_binWidth);
            //                cout << e_nu << endl;
            //                cout<<m_binWidth<<endl;
            //                  cout << "enuidex is " << e_nu_Idx << endl;
            //                cout << "idx" << idx << endl;
            
            //  cout << e_nu << endl; Checked
            
            if(e_nu_Idx<0 || e_nu_Idx>=MatrixBins)
            {
                m_positronTrueSpectrum[idet][idx] = 0;
                //                    cout << idx << "this is over or under limits";
                
                continue;
            }
            double binScaling = enu_to_epositron_func->Derivative(e_nu,0,0.0001);
            //cout << e_nu << " " << e_positron << " " << binScaling << endl; //checked
            double dNdE = 0;
            int sign = 1;
            double dE = e_nu - (0.5 + e_nu_Idx)*m_binWidth;
            //                cout << "dE" << dE << endl;
            if (dE < 0) sign = -1;
            
            if(e_nu_Idx==MatrixBins-1)
            {
                dNdE =PredictionH[0][0][0]->GetBinContent(e_nu_Idx+1);
                //                    cout<< dNdE << " dNdE, last index " << e_nu_Idx << endl; //checked
                //                    cout<< dNdE; //checked
            }
            else
            {
                dNdE =(TMath::Abs(dE)/m_binWidth)*PredictionH[0][0][0]->GetBinContent(e_nu_Idx+sign+1)+(1 - TMath::Abs(dE)/m_binWidth)*PredictionH[0][0][0]->GetBinContent(e_nu_Idx+1);
                //                 cout<< dNdE << " dNdE, index " << e_nu_Idx << endl; //checked
                
                //    cout << PredictionH[0][0][0]->GetBinContent(idx+1) << "bin" << idx+1 << endl;
                
            }
            //        cout << energy_nl << "\t" << energy_true << "\t" << dNdE/binScaling << endl;
            m_positronTrueSpectrum[idet][idx] = dNdE/ binScaling;
            //                cout << "this is correct";
            
            PositronTrueSpectrumH[idet]->SetBinContent(idx+1, m_positronTrueSpectrum[idet][idx]);
            //     PositronTrueSpectrumH[idet]->Draw();
        }
    }
//IAV
    for(int idet=0; idet<NADs; idet++)
    {
        for(int idx=0; idx<MatrixBins; idx++)
        {
            m_positronIavDistortedSpectrum[idet][idx] = 0;
        }
        
        for(int idx=0; idx<MatrixBins; idx++)
        {
            for(int jdx=0; jdx<idx+1; jdx++)
            {
                m_positronIavDistortedSpectrum[idet][jdx]+= m_iav_corr[idet][idx][jdx] * m_positronTrueSpectrum[idet][idx];
            }
            PositronIAVSpectrumH[idet]->SetBinContent(idx+1, m_positronIavDistortedSpectrum[idet][idx]);
        }
    }

//Non linearity
    double alpha = 0.08;
    double beta = 1 * 1;
    for(int idet=0; idet<NADs; idet++)
    {
        if(beta==1.0 && alpha==0.0)
        {
            for(int idx=0; idx<MatrixBins; idx++)
            {
                m_positronNLSpectrum[idet][idx] = m_positronIavDistortedSpectrum[idet][idx];
            }
        }
        else
        {
            for(int idx=0; idx<MatrixBins; idx++)
            {
                double energy_nl = m_energy[idet][idx];
                double energy_true = 0.5*((energy_nl/beta)+ TMath::Sqrt(4*alpha+ ((energy_nl*energy_nl)/(beta*beta))));
                int eTrueIdx = (int)((energy_true-0)/m_binWidth);
                
                if(eTrueIdx<0 || eTrueIdx>=MatrixBins)
                {
                    m_positronNLSpectrum[idet][idx] = 0;
                    continue;
                }
                
                double binScaling = (1 + alpha/(energy_true*energy_true))*beta;
                double dNdE = 0;
                double dEtrue = energy_true - eTrueIdx*m_binWidth;
                
                if(eTrueIdx==MatrixBins-1)
                {
                    dNdE = m_positronIavDistortedSpectrum[idet][eTrueIdx];
                }else
                {
                    // Linear interpolate between points
                    dNdE = (m_positronIavDistortedSpectrum[idet][eTrueIdx]+((dEtrue/m_binWidth)*(m_positronIavDistortedSpectrum[idet][eTrueIdx+1]-m_positronIavDistortedSpectrum[idet][eTrueIdx])));
                }
                m_positronNLSpectrum[idet][idx] = dNdE / binScaling;
                // cout << "NL bin is " << m_positronNLSpectrum[idet][idx];
                PositronNLSpectrumH[idet]->SetBinContent(idx+1,  m_positronNLSpectrum[idet][idx]);
            }
        }
    }


//Resolution
    double resolutionRange = 8; // sigma
    for(int idet=0; idet<NADs; idet++)
    {
        for(int idx=0; idx<MatrixBins; idx++)
        {
            m_positronDetectedSpectrum[idet][idx] = 0;
        }
        
        for(int idx=0; idx<MatrixBins; idx++)
        {
//            if(m_detectorResolution == 0)
//            {
//                m_positronDetectedSpectrum[idet][idx] = m_positronNLSpectrum[idet][idx];
//                continue;
//            }
            if(m_positronNLSpectrum[idet][idx]==0) continue;
            double energy_nl = m_energy[idet][idx];
            //      double sigma = m_detectorResolution * TMath::Sqrt(energy_nl);
            double sigma = (reso_func->Eval(energy_nl) + m_detectorResolution_bias[idet]) * energy_nl;
           // cout << sigma << endl;
            double minDetE = energy_nl - resolutionRange*sigma;
            double maxDetE = energy_nl + resolutionRange*sigma;
            int minDetEIdx = (int)((minDetE-0)/m_binWidth);
            int maxDetEIdx = (int)((maxDetE-0)/m_binWidth);
            if(minDetEIdx < 0) minDetEIdx = 0;
            if(maxDetEIdx >= MatrixBins) maxDetEIdx = MatrixBins-1;
            
            for(int detIdx=minDetEIdx; detIdx<=maxDetEIdx; detIdx++)
            {
                if(detIdx==0) continue;
                double gausFactor = TMath::Gaus((energy_nl-m_energy[idet][detIdx]),0,sigma,true);
                
                m_positronDetectedSpectrum[idet][detIdx] += (m_positronNLSpectrum[idet][idx]*gausFactor* m_binWidth);
            }
            PositronDetectedSpectrumH[idet]->SetBinContent(idx+1, m_positronDetectedSpectrum[idet][idx]);

        }

    }
    
    delete enu_to_epositron_func;
    delete reso_func;

    for(Int_t idet =0;idet<NADs;idet++)
    {
        delete PositronTrueSpectrumH[idet];
        delete PositronIAVSpectrumH[idet];
        delete PositronNLSpectrumH[idet];
        delete PositronDetectedSpectrumH[idet];
    }
    
    return m_positronDetectedSpectrum[idet];
}


TF1* EnergyMatrix::GetNeutrinoToVisibleFunction(Int_t order)
{
    if(order==0)
    {
        Double_t Correction = Mn-Mp-Me;
        TF1 *VisibleF = new TF1("VisibleF","x-[0]",InitialEnergy,FinalEnergy,1);
        VisibleF->SetParameter(0,Correction);
    }
    if(order==1)
    {
        TF1 *VisibleF = new TF1("VisibleF",VisibleEnergyF,InitialEnergy,FinalEnergy,1);
        VisibleF->SetParameter(0,0);//Maybe I can improve this using a MC simulation of the angular distribution calculated in http://authors.library.caltech.edu/2796/1/VOGprd99.pdf
    }
    return VisibleF;
}
// First order true to visible function.
//(x-(Mn-Mp))*(1 - x/Mn*(1.0 - [0]*sqrt(1 - Me*Me/(x-(Mn-Mp))/(x-(Mn-Mp)))))- ((Mn-Mp)*(Mn-Mp) - Me*Me)/(2*Mn)-Me]
Double_t VisibleEnergyF(Double_t * x, Double_t * par)
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

void EnergyMatrix::LoadPredictions()
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
            }
        }
    }
    FarHallPredictionsF->Close();
}

Double_t enu_to_epositron(Double_t * x, Double_t * par){
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

void EnergyMatrix::LoadIavCorrection(){
    
    TFile * f = new TFile("./IavDistortion/IAVDistortion.root");
    TH2F * Correction = (TH2F*)f->Get("Correction");
    
    cout << "reading Iav correction file" << endl;
    for(int i=0;i<MatrixBins;i++)
    { // i: true positron energy bin; j: distorted energy bin
        //Total events in input spectrum bin i
        m_nominal_iav_frac[i] = 0;
        if (Correction->Integral(i+1,i+1,0,-1) > 0)
        {
            for(int j=0;j<i+1;j++)
            {
                m_nominal_iav_corr[i][j]= Correction->GetBinContent(i+1,j+1)/Correction->Integral(i+1,i+1,0,-1);
                if (i!=j)
                {
                    m_nominal_iav_frac[i] += m_nominal_iav_corr[i][j];
                }
            }
        }
        else
        {
            for(int j=0;j<i+1;j++)
            {
                if (i==j)
                {
                    m_nominal_iav_corr[i][j] = 1;
                }
                else
                {
                    m_nominal_iav_corr[i][j] = 0;
                }
            }
            m_nominal_iav_frac[i] = 0;
        }
    }
    
    for(int idet=0;idet<NADs;idet++)
    {
        for(int i=0;i<MatrixBins;i++)
        { // i: true positron energy bin; j: distorted energy bin
            m_iav_frac[idet][i] = m_nominal_iav_frac[i];
            for(int j=0;j<MatrixBins;j++)
            {
                m_iav_corr[idet][i][j] = m_nominal_iav_corr[i][j];
            }
        }
    }
    f->Close();

}


Double_t reso_func_bcw(Double_t * x, Double_t * par){
    
    Double_t e_orig = x[0];
    Double_t e_sigma = 1.0;
    
    if (e_orig > 0)
    {
        e_sigma = TMath::Sqrt(par[0]*par[0] + par[1]*par[1]/e_orig + par[2]*par[2]/e_orig/e_orig);
    }
    
    return e_sigma;
}

void EnergyMatrix::SaveTotalSpectrum()
{
    TFile* SaveSpectrumDataF = TFile::Open("./RootOutputs/EnergyMatrix.root","recreate");
    for (Int_t idet = 0;idet<NADs; idet++)
    {
//            PositronTrueSpectrumH[idet]->Write();
//            PositronIAVSpectrumH[idet]->Write();
//            PositronNLSpectrumH[idet]->Write();
//            PositronDetectedSpectrumH[idet]->Write();
            h_evis_vs_enu_ad[idet]->Write();
        htmp_rebin->Write();
       // htmp_rebin1->Write();


    }
    SaveSpectrumDataF->Close();
}