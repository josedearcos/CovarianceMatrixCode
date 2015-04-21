#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"

int CreateWeeklyFluxRootFile(){
    
    TH1::AddDirectory(kFALSE);
    const Double_t conv_m2_per_cm2 = 1.0e-4;
    const Double_t conv_s_per_day = 24*60*60;
    
    //Load Cross Section:
    TFile* CrossSectionF = new TFile("../../CrossSections/nHCrossSection.root");//Zhe's file
    TH1D* CrossSectionH = (TH1D*)gDirectory->Get("CrossSection");
    CrossSectionH->Scale(10e-42);
    delete CrossSectionF;
    
    const Int_t NReactors=6;
    const Int_t NIHEPReactorBins = 33; // 1.5, 1.75, ..., 9.25, 9.5 in 0.25 MeV steps
    Double_t eMin=1.5;
    Double_t eMax=9.5;
    Double_t binWidth = 0.25;
    const Int_t MaxPeriods = 101;
    TH1D* ReactorSpectrum[NReactors][MaxPeriods];
    TH1D* InclusiveReactorSpectrum[NReactors];

    Double_t m_dNdE_nom[NReactors*(NIHEPReactorBins+1)*MaxPeriods];

    Int_t curPeriod;
    std::fstream fileData;
    Double_t trash;
    
    for (Int_t reactor = 0; reactor < NReactors; reactor++)
    {
        curPeriod = 0;
        
        switch(reactor)
        {
            case 0://Daya Bay A
                fileData.open("./DayaBayA_2011-12-24_2013-11-27.txt",std::fstream::in);
                break;
            case 1://Daya Bay B
                fileData.open("./DayaBayB_2011-12-24_2013-11-27.txt",std::fstream::in);
                break;
            case 2://LingAo IA
                fileData.open("./LingAoIA_2011-12-24_2013-11-27.txt",std::fstream::in);
                break;
            case 3://LingAo IB
                fileData.open("./LingAoIB_2011-12-24_2013-11-27.txt",std::fstream::in);
                break;
            case 4://LingAo IIA
                fileData.open("./LingAoIIA_2011-12-24_2013-11-27.txt",std::fstream::in);
                break;
            case 5://LingAo IIB
                fileData.open("./LingAoIIB_2011-12-24_2013-11-27.txt",std::fstream::in);
                break;
            default:
                std::cout << "Reactor data file not found" << std::endl;
                
                exit(EXIT_FAILURE);
        }
        
        if(!fileData.is_open())
        {
            std::cout << "P14A Reactor file not open" << std::endl;
            
            exit(EXIT_FAILURE);
        }
        
        while(!fileData.eof())
        {
            for(Int_t idr = 0; idr<5; idr++)//first 5 columns are the reactor ID
            {
                fileData >> trash;
                
                //std::cout << trash << std::endl;
                
            }
            for (Int_t bin = 0; bin <=NIHEPReactorBins; bin++)
            {
                fileData >> m_dNdE_nom[reactor+NReactors*bin+curPeriod*NReactors*NIHEPReactorBins];
            }
            curPeriod++;
        }
        curPeriod--;
        
        fileData.close();
    }
    
    static bool linearInterp = true;// By default, use linear interpolation between nearest sample points

    for (Int_t reactor = 0; reactor < NReactors; reactor++)
    {
        InclusiveReactorSpectrum[reactor] = new TH1D(Form("%i",reactor),Form("%i",reactor),NIHEPReactorBins-1,1.5,9.5);
        
        for(Int_t period = 0; period<curPeriod; period++)
        {
            ReactorSpectrum[reactor][period] = new TH1D(Form("%i%i",reactor,period),Form("%i%i",reactor,period),NIHEPReactorBins-1,1.5,9.5);
            
            for (Int_t bin = 0; bin <=NIHEPReactorBins; bin++)
            {
                Int_t binIdxLow = bin;
                Double_t e_nu = eMin + binIdxLow*binWidth + binWidth/2;//center of the bin
                Double_t value = m_dNdE_nom[reactor+NReactors*bin+period*NReactors*NIHEPReactorBins];
                
                if(linearInterp)
                {
                    Int_t binIdxHigh = binIdxLow+1;
                    Double_t binLowE = eMin + binIdxLow*binWidth;
                    Double_t dE = e_nu - binLowE;
                    Double_t slope = (m_dNdE_nom[reactor+NReactors*binIdxHigh+period*NReactors*NIHEPReactorBins]
                                      - m_dNdE_nom[reactor+NReactors*binIdxLow+period*NReactors*NIHEPReactorBins])/binWidth;
                    value += dE*slope;
                }
                ReactorSpectrum[reactor][period]->SetBinContent(bin+1,value*CrossSectionH->Interpolate(e_nu));
            }
            
            InclusiveReactorSpectrum[reactor]->Add(ReactorSpectrum[reactor][period]);
        }
    }
    
    TFile* SaveWeeklyFluxF = new TFile(Form("../WeeklyFlux_%dweek_unblinded.root",curPeriod),"recreate");
    
    for(Int_t period = 0; period<curPeriod; period++)
    {
        TDirectory* TotalDirectory = SaveWeeklyFluxF->mkdir(Form("Week%i",period));

        SaveWeeklyFluxF->cd(Form("Week%i",period));

        for (Int_t reactor = 0; reactor < NReactors; reactor++)
        {
            ReactorSpectrum[reactor][period]->SetName(Form("%i",reactor));
            ReactorSpectrum[reactor][period]->SetTitle(Form("%i",reactor));

            ReactorSpectrum[reactor][period]->Write();
        }
    }
    
    delete SaveWeeklyFluxF;
    
    TFile* SaveInclusiveWeeklyFluxF = new TFile(Form("../WeeklyFlux_%dweek_unblinded_inclusive.root",curPeriod),"recreate");

        TDirectory* TotalDirectory = SaveInclusiveWeeklyFluxF->mkdir(Form("Week%i",0));
        
        SaveInclusiveWeeklyFluxF->cd(Form("Week%i",0));
        
        for (Int_t reactor = 0; reactor < NReactors; reactor++)
        {
            InclusiveReactorSpectrum[reactor]->Write();
        }
    
    delete SaveInclusiveWeeklyFluxF;
    
    return 0;
}