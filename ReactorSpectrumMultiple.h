#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "NominalData.h"
#include <vector>
#include <math.h>
#include <iomanip>
#include <streambuf>
#include <TLegend.h>

//const Int_t NReactors=6; //Already declared in NominalData
//const Int_t NIsotopes=4;
const Double_t conv_mev_per_s_gw = 6.2415e21;
//const Int_t m_bcw_nbins = 26;
const bool ChristineReactorModel = 1;
//const Int_t Nbins=240-36;

class ReactorSpectrumMultiple
{
private:
    //    std::vector<TH1D*> Prueba;
    std::vector<TH1D*> NominalSpectra;
    //    TH2D* CovMatrixM;
    //    TH2D* CovMatrixLM;
    TH1D* SpectrumHist;
    
    bool Mode;
    
    Double_t InitialEnergy;
    Double_t FinalVisibleEnergy;
    Double_t BinWidth;
    
    Double_t m_eMin;
    Double_t m_eMax;
    Int_t m_nSamples;
    Double_t m_binWidth;
    
    Double_t e_nu;
    Double_t eMin;
    Double_t eMax;
    Double_t eNu;
    Double_t dNdE[NReactors];
    Int_t curSample;
    
    Double_t NominalIsotopeFrac[NIsotopes];
    Double_t IsotopeFrac[NIsotopes];
    Double_t IsotopeError[NIsotopes];
    Double_t RandomIsotopeError[NIsotopes];
    
    Double_t NominalReactorPower[NReactors];
    Double_t ReactorPower[NReactors];
    Double_t ReactorPowerError[NReactors];
    Double_t RandomReactorPowerError[NReactors];
    Double_t EnergyPerFission[NIsotopes];
    Double_t EnergyPerFissionError[NIsotopes];
    Double_t SumEnergyPerFission;
    Double_t* Energy;
    Double_t* Spectrum;
    
    bool IsotopeMatrix;
    bool ReactorPowerMatrix;
    
    Double_t m_dNdE_nom[NReactors][MatrixBins];
    Double_t m_dNdE[NReactors][MatrixBins];
    Double_t L[NReactors*MatrixBins][NReactors*MatrixBins]; // lower triangle of the covariance matrix
    
    TRandom3* rand;
    Double_t Norm;
    NominalData* Nom;
    
    Int_t ADs;
    Int_t Reactors;
    Int_t OriginalNbins;
    
    std::vector<TH1D*> Pu239H;
    std::vector<TH1D*> Pu241H;
    std::vector<TH1D*> U235H;
    std::vector<TH1D*> U238H;
    
    std::vector<TH1D*> TotalReactorSpectrumH;
    void File2Hist(const char* IsotopeName);
    void Plot(TH1D*);
    void RandomIsotopes();
    void RandomPower();
    
public:
    
    ReactorSpectrumMultiple();
    ~ReactorSpectrumMultiple();
    ReactorSpectrumMultiple(NominalData* Data);
    void MultipleReactorSpectrumMain(bool);
};

//Default values
ReactorSpectrumMultiple :: ReactorSpectrumMultiple()
{
    std::cout << " the reactor default constructor shouldn't be called" << std::endl;
    
    exit(EXIT_FAILURE);
    
    Nom = new NominalData(0,2);
    rand = new TRandom3(0);
    
    std::setprecision(40);
    
    //Binning variables
    OriginalNbins = 820;
    InitialEnergy = Nom->GetEmin();
    FinalVisibleEnergy = Nom->GetEVisMax();
    BinWidth = (FinalVisibleEnergy-InitialEnergy)/Nbins;
    
    //Isotope Fractions
    for (Int_t isotope=0; isotope<NIsotopes; isotope++)
    {
        NominalIsotopeFrac[isotope] = Nom->GetIsotopeFraction(isotope); // Load nominal Isotope fractions
        IsotopeFrac[isotope]=NominalIsotopeFrac[isotope];
        IsotopeError[isotope]= Nom->GetIsotopeFractionError(isotope); //Isotope errors (5% for all) http://dayabay.ihep.ac.cn/DocDB/0086/008609/002/reactor_technote.pdf
        EnergyPerFission[isotope]= Nom->GetEnergyPerFission(isotope);
        EnergyPerFissionError[isotope]= Nom->GetEnergyPerFissionError(isotope);
    }
    
    //Reactor Power
    for(Int_t r =0;r<NReactors;r++)
    {
        NominalReactorPower[r]=Nom->GetReactorPower(r);
        ReactorPower[r]=NominalReactorPower[r];
        ReactorPowerError[r]=Nom->GetReactorPowerError(r);
    }
    
    ADs = Nom->GetADs();
    
    //reactor:
    
    m_nSamples = Nom->GetReactorSamples();
    m_binWidth = Nom->GetReactorBinWidth();
    m_eMin = Nom->GetReactorEmin();
    m_eMax = Nom->GetReactorEmax();
    
    for (Int_t i = 0; i < m_nSamples * NReactors; i++)
    {
        for (Int_t j = 0; j <  m_nSamples * NReactors; j++)
        {
            L[i][j]=Nom->GetReactorCovMatrix(i,j);
        }
    }
    
    for(Int_t reactor=0; reactor<NReactors;reactor++)
    {
        for(Int_t sample=0; sample<m_nSamples;sample++)
        {
            m_dNdE_nom[reactor][sample] = Nom->GetNominalReactorSpectrum(reactor,sample);
            m_dNdE[reactor][sample]=m_dNdE_nom[reactor][sample];
        }
    }
    
    IsotopeMatrix = Nom->GetIsotopeMatrix();
    ReactorPowerMatrix = Nom->GetReactorPowerMatrix();
    
    delete Nom;
}

ReactorSpectrumMultiple :: ReactorSpectrumMultiple(NominalData* Data)
{
    rand = new TRandom3(0);
    
    //Binning variables
    OriginalNbins = 820;
    InitialEnergy = Data->GetEmin();
    FinalVisibleEnergy = Data->GetEVisMax();
    BinWidth = (FinalVisibleEnergy-InitialEnergy)/Nbins;
    
    //Isotope Fractions
    for (Int_t isotope=0; isotope<NIsotopes; isotope++)
    {
        NominalIsotopeFrac[isotope] = Data->GetIsotopeFraction(isotope); // Load nominal Isotope fractions
        IsotopeFrac[isotope]=NominalIsotopeFrac[isotope];
        IsotopeError[isotope]= Data->GetIsotopeFractionError(isotope); //Isotope errors (5% for all) http://dayabay.ihep.ac.cn/DocDB/0086/008609/002/reactor_technote.pdf
        EnergyPerFission[isotope]= Data->GetEnergyPerFission(isotope);
        EnergyPerFissionError[isotope]= Data->GetEnergyPerFissionError(isotope);
    }
    
    //Reactor Power
    for(Int_t r =0;r<NReactors;r++)
    {
        NominalReactorPower[r]=Data->GetReactorPower(r);
        ReactorPower[r]=NominalReactorPower[r];
        ReactorPowerError[r]=Data->GetReactorPowerError(r);
    }
    
    ADs = Data->GetADs();
    
    //reactor:
    
    m_nSamples = Data->GetReactorSamples();
    m_binWidth = Data->GetReactorBinWidth();
    m_eMin = Data->GetReactorEmin();
    m_eMax = Data->GetReactorEmax();
    
    for (Int_t i = 0; i < m_nSamples * NReactors; i++)
    {
        for (Int_t j = 0; j <  m_nSamples * NReactors; j++)
        {
            L[i][j]=Data->GetReactorCovMatrix(i,j);
        }
    }
    
    for(Int_t reactor=0; reactor<NReactors;reactor++)
    {
        for(Int_t sample=0; sample<m_nSamples;sample++)
        {
            m_dNdE_nom[reactor][sample] = Data->GetNominalReactorSpectrum(reactor,sample);
            m_dNdE[reactor][sample]=m_dNdE_nom[reactor][sample];
        }
    }
    
    IsotopeMatrix = Data->GetIsotopeMatrix();
    ReactorPowerMatrix = Data->GetReactorPowerMatrix();
}

ReactorSpectrumMultiple ::~ReactorSpectrumMultiple()
{
    delete rand;
    
    for(Int_t reactor=0; reactor<NReactors;reactor++)
    {
        delete TotalReactorSpectrumH[reactor];
    }
}
void ReactorSpectrumMultiple :: MultipleReactorSpectrumMain(bool mode)
{
    Mode = mode;
    
    Char_t FileName[200];
    std::cout <<  "\t ***********************************************************************************************" << std::endl;
    std::cout << " \t Calculating reactor spectrum " << std::endl;
    
    if((IsotopeMatrix||ReactorPowerMatrix)&&Mode==1)
    {
        sprintf(FileName,"./RootOutputs/Reactor/RandomOutputs/ReactorSpectrum_Isotope_%d_Power_%d.root",IsotopeMatrix,ReactorPowerMatrix);
    }
    else
    {
        sprintf(FileName,"./RootOutputs/Reactor/NominalOutputs/ReactorSpectrum.root");
    }
    
    if(ChristineReactorModel)
    {
        NominalSpectra.resize(NReactors);
        
        for (Int_t reactor = 0; reactor < NReactors; reactor++)
        {
            NominalSpectra[reactor] = new TH1D(Form("Nominal%d",reactor),Form("Nominal%d",reactor),218,1.85,12.75);
        }
        for(Int_t curSample = 0; curSample<m_nSamples;curSample++)
        {
            for (Int_t reactor = 0; reactor < NReactors; reactor++)
            {
                NominalSpectra[reactor]->SetBinContent(curSample+1,m_dNdE[reactor][curSample] * 1.0e18);
            }
        }
        //        CovMatrixM = new TH2D("COV","COV",m_nSamples*NReactors,0,72,m_nSamples*NReactors,0,72);
        //        CovMatrixLM = new TH2D("COVL","COVL",m_nSamples*NReactors,0,72,m_nSamples*NReactors,0,72);
        
        //  Apply variations to the nominal spectrum:
        if (IsotopeMatrix&&Mode==1)
        {
            Double_t ranvec[NReactors*MatrixBins];
            
            for (Int_t i = 0; i < NReactors*m_nSamples; i++)
            {
                rand->SetSeed(0);
                ranvec[i] = rand->Gaus(0,1);
            }
            
            for (Int_t iCore = 0; iCore < NReactors; iCore++)
            {
                for (Int_t iSample = 0; iSample < m_nSamples; iSample++)
                {
                    m_dNdE[iCore][iSample] = m_dNdE_nom[iCore][iSample];
                    
                    for (Int_t jCore = 0; jCore < NReactors; jCore++)
                    {
                        for (Int_t jSample = 0; jSample < m_nSamples; jSample++)
                        {
                            //  Multiply: Correlated Vector = Lower Decomposed Correlation Matrix * Uncorrelated Gaussian Vector
                            m_dNdE[iCore][iSample] += L[iCore*m_nSamples+iSample][jCore*m_nSamples+jSample] * ranvec[jCore*m_nSamples+jSample];
                        }
                    }
                    //                          std::cout << "dNdE, correlated vector, random vs nominal: " << m_dNdE[iCore][iSample]  << "\t" << m_dNdE_nom[iCore][iSample]  << std::endl;
                }
            }
        }
        
        Int_t binIdxLow;
        Int_t binIdxHigh;
        
        Double_t value_low;
        Double_t value_high;
        Double_t value;
        
        Double_t binLowE;
        Double_t dE;
        Double_t slope;
        
        TotalReactorSpectrumH.resize(NReactors);
        //        Prueba.resize(NReactors);
        
        if(ReactorPowerMatrix&&Mode==1)
        {
            RandomPower();
        }
        
        for (Int_t reactor = 0; reactor < NReactors; reactor++)
        {
            //            Prueba[reactor] = new TH1D(Form("Prueba%d",reactor+1),Form("Prueba%d",reactor+1),m_nSamples,InitialEnergy,FinalVisibleEnergy);
            
            //            for (Int_t i = 0; i < m_nSamples; i++)
            //            {
            //                Prueba[reactor]->SetBinContent(i+1,m_dNdE[reactor][i] * 1.0e18);
            //            }
            
            //  Apply power variations in all reactors
            
            
            TotalReactorSpectrumH[reactor] = new TH1D(Form("SpectrumFromReactor%d",reactor+1),Form("SpectrumFromReactor%d",reactor+1),Nbins,InitialEnergy,FinalVisibleEnergy);
            
            //  Interpolate for missing points and write reactor spectrum
            for (Int_t i = 0; i < Nbins; i++)
            {
                e_nu = i*BinWidth+InitialEnergy;
                
                //                std::cout << " e_nu " << e_nu << std::endl;
                
                // From LBNL's notes:
                // By default, use linear interpolation between nearest sample points
                // Christine's spectra starts from 1.85 MeV, while we want spectra from 1.80 MeV. So changed it so that it can do linear extrapolation below 1.85 MeV.
                
                binIdxLow = 0;
                
                if (e_nu > m_eMin)
                {
                    binIdxLow = (Int_t)((e_nu-m_eMin)/m_binWidth);
                }
                
                binIdxHigh = binIdxLow+1;
                
                value_low = m_dNdE[reactor][binIdxLow];
                value_high = m_dNdE[reactor][binIdxHigh];
                value = value_low;
                
                binLowE = m_eMin + binIdxLow*m_binWidth;
                dE = e_nu - binLowE;
                slope = (value_high - value_low)/m_binWidth;
                value += dE*slope;
                
                if(e_nu>=m_eMax)
                {
                    value = 0.0;
                }
                if (value < 0)
                {
                    value = 0.0;
                }
                
                TotalReactorSpectrumH[reactor]->SetBinContent(i+1, value * ReactorPower[reactor]/NominalReactorPower[reactor] * 1.0e18);
            }
            std::cout << " Reactor factor " <<  reactor << " " << ReactorPower[reactor]/NominalReactorPower[reactor] << std::endl;
            
        }
        
        TFile* ReactorFile = new TFile(FileName, "recreate");
        
        for(Int_t reactor=0; reactor<NReactors;reactor++)
        {
            TotalReactorSpectrumH[reactor]->Write();
            //            Prueba[reactor]->Write();
            NominalSpectra[reactor]->Write();
            //                delete Prueba[reactor];
        }
        delete ReactorFile;

        if(Print)
        {
            TCanvas* SaveReactor = new TCanvas("Reactor","Reactor",1200,1200);
            
            TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
            legend->SetTextFont(72);
            legend->SetTextSize(0.02);
            legend->SetFillColor(0);
            for(Int_t reactor=0; reactor<NReactors;reactor++)
            {
                legend->AddEntry(NominalSpectra[reactor],Form("Core #%d",reactor),"l");

                SaveReactor->cd((NReactors+reactor)+1);
                NominalSpectra[reactor]->SetStats(0);
                NominalSpectra[reactor]->SetTitle("");
                NominalSpectra[reactor]->SetLineColor(reactor+1);
                NominalSpectra[reactor]->GetXaxis()->SetTitle("E_{true} (MeV)");
                NominalSpectra[reactor]->GetXaxis()->SetTitleSize(0.04);
                NominalSpectra[reactor]->GetYaxis()->SetTitleSize(0.04);
                // Units: [neutrinos MeV^-1 fission^-1]

                NominalSpectra[reactor]->GetYaxis()->SetTitle("#nu_{e} MeV^{-1} fission^{-1}");
                NominalSpectra[reactor]->Draw("same");
                legend->Draw("same");

                SaveReactor->Modified();
            }
            SaveReactor->Update();
            
            SaveReactor->Print("./Images/Reactor/ChristineReactorSpectrum.eps",".eps");
            delete SaveReactor;
        }
        
        for(Int_t reactor=0; reactor<NReactors;reactor++)
        {
            delete NominalSpectra[reactor];
        }
        //        CovMatrixM->Write();
        //        CovMatrixLM->Write();
        //        delete CovMatrixM;
        //        delete CovMatrixLM;
    }
    else
    {
        //  A simple spectrum, no SNF, non-eq corrections or burn-up has been applied, neither time information. In a fast comparison to Christine's spectrum, it looks as it has an excess of events at high energies (probably due to the information provided (non updated txt inputs), or extrapolation method below) and the peak is shifted 0.7MeV (around 4.2MeV). Use this just as a provisional model but not as a final version.
        
        TotalReactorSpectrumH.resize(NReactors);
        
        Char_t filenameMultiple[1024];
        Char_t filenamePu239[1024];
        Char_t filenamePu241[1024];
        Char_t filenameU235[1024];
        Char_t filenameU238[1024];
        
        File2Hist("./ReactorInputs/antiNeuFlux_Pu2392011.dat.txt");
        TH1D* Pu239Nominal = (TH1D*)SpectrumHist->Clone();
        delete SpectrumHist;
        File2Hist("./ReactorInputs/antiNeuFlux_Pu2412011.dat.txt");
        TH1D* Pu241Nominal = (TH1D*)SpectrumHist->Clone();
        delete SpectrumHist;
        File2Hist("./ReactorInputs/antiNeuFlux_U2352011.dat.txt");
        TH1D* U235Nominal = (TH1D*)SpectrumHist->Clone();
        delete SpectrumHist;
        File2Hist("./ReactorInputs/antiNeuFlux_U2382011.dat.txt");
        TH1D* U238Nominal = (TH1D*)SpectrumHist->Clone();
        delete SpectrumHist;
        
        U235H.resize(NReactors);
        U238H.resize(NReactors);
        Pu239H.resize(NReactors);
        Pu241H.resize(NReactors);
        
        TFile* outputFile = new TFile(FileName, "recreate");
        
        //First update fission fractions
        
        if (IsotopeMatrix)
        {
            RandomIsotopes();
        }
        if (ReactorPowerMatrix)
        {
            RandomPower();
        }
        
        for(Int_t reactor=0; reactor<NReactors;reactor++)
        {
            sprintf(filenameU235,"U235 Isotope spectrum for reactor %i", reactor+1);
            sprintf(filenameU238,"U238 Isotope spectrum for reactor %i", reactor+1);
            sprintf(filenamePu239,"Pu239 Isotope spectrum for reactor %i", reactor+1);
            sprintf(filenamePu241,"Pu241 Isotope spectrum for reactor %i", reactor+1);
            
            U235H[reactor] = (TH1D*)U235Nominal->Clone(filenameU235);
            U238H[reactor] = (TH1D*)U238Nominal->Clone(filenameU238);
            Pu239H[reactor] = (TH1D*)Pu239Nominal->Clone(filenamePu239);
            Pu241H[reactor] = (TH1D*)Pu241Nominal->Clone(filenamePu241);
            
            U235H[reactor]->Scale(IsotopeFrac[0]);
            U238H[reactor]->Scale(IsotopeFrac[1]);
            Pu239H[reactor]->Scale(IsotopeFrac[2]);
            Pu241H[reactor]->Scale(IsotopeFrac[3]);
            
            Plot(U235H[reactor]);
            Plot(U238H[reactor]);
            Plot(Pu239H[reactor]);
            Plot(Pu241H[reactor]);
            
            sprintf(filenameMultiple,"SpectrumFromReactor%i", reactor+1);
            
            TotalReactorSpectrumH[reactor] = (TH1D*)U235H[reactor]->Clone(filenameMultiple);
            TotalReactorSpectrumH[reactor]->Add(U238H[reactor]);
            TotalReactorSpectrumH[reactor]->Add(Pu239H[reactor]);
            TotalReactorSpectrumH[reactor]->Add(Pu241H[reactor]);
            TotalReactorSpectrumH[reactor]->SetTitle(filenameMultiple);
            
            delete U235H[reactor];
            delete U238H[reactor];
            delete Pu239H[reactor];
            delete Pu241H[reactor];
        }
        
        //Calculate average energy per fission
        SumEnergyPerFission=0;
        for (Int_t isotope=0; isotope<NIsotopes; isotope++)
        {
            SumEnergyPerFission+=IsotopeFrac[isotope]*EnergyPerFission[isotope];
        }
        
        for(Int_t reactor=0; reactor<NReactors;reactor++)
        {
            Double_t AverageEnergyPerFission = (ReactorPower[reactor]/SumEnergyPerFission)*conv_mev_per_s_gw ;
            TotalReactorSpectrumH[reactor]->Scale(AverageEnergyPerFission);
            TotalReactorSpectrumH[reactor]->Write();
        }
        delete outputFile;
        
        if(Print)
        {
            //Scale with reactor power and Energy per fission
            TCanvas* SaveReactor = new TCanvas("Reactor","Reactor");
            
            SaveReactor->Divide(NReactors,2);
            
            for(Int_t reactor=0; reactor<NReactors;reactor++)
            {
                SaveReactor->cd(reactor+1);
                
                TotalReactorSpectrumH[reactor]->Draw();
            }
            
            SaveReactor->Print("./Images/Reactor/ReactorSpectrum.eps",".eps");
            delete SaveReactor;
        }
        
        delete[] Energy;
        delete[] Spectrum;
    }
    
    std::cout << "\t Reactor Spectrum calculated" << std::endl;
}

void ReactorSpectrumMultiple :: File2Hist (const char* IsotopeName)
{
    Energy = (Double_t*)malloc(OriginalNbins*sizeof(Double_t));
    Spectrum = (Double_t*)malloc((OriginalNbins+200)*sizeof(Double_t));
    
    SpectrumHist = new TH1D(IsotopeName, IsotopeName, Nbins, InitialEnergy, FinalVisibleEnergy);
    
    std::string line;
    std::ifstream input(IsotopeName);
    
    std::getline(input,line); //To throw away the first two lines
    std::getline(input,line);
    
    for(Int_t i=0;i<=OriginalNbins;i++)
    {
        input >> Energy[i] >> Spectrum[i];
        
        //       printf("Index is %d \n",i);
        //       printf("Energy is %.9f \n",Energy[i]);
        //       printf("Spectrum is %.10f \n",Spectrum[i]);
    }
    
    //I will extrapolate the reactor spectrum from 10 to 12MeV so it matches the rest of the inputs
    //There are 2 MeV to fill, steps are 0.01 MeV so we need 200 extra bins.
    for(Int_t i=OriginalNbins;i<OriginalNbins+200;i++)
    {
        Spectrum[i+1]=Spectrum[i-1]+2*(Spectrum[i]-Spectrum[i-1]);
        if(Spectrum[i]<=0)
        {
            Spectrum[i]=Spectrum[i-1];//This to extrapolate
            //            Spectrum[i]=0;//This to just pad zeros
        }
    }
    
    for(Int_t i=0;i<Nbins;i++)
    {
        SpectrumHist->SetBinContent(i+1, Spectrum[i*(OriginalNbins+200)/(Nbins+36)]);
    }
    
    
    //    for(Int_t i=0;i<Nbins;i++)
    //    {
    //        SpectrumHist->SetBinContent(i+1, Spectrum[i*(OriginalNbins)/(Nbins)]);
    //    }
}

void ReactorSpectrumMultiple :: Plot(TH1D* Hist)
{
    Hist->GetXaxis()->SetTitle("E_{#nu} [MeV]");
    Hist->GetYaxis()->SetTitle("#nu [MeV^{-1} fission^{-1}]");
    Hist->Write();
}

void ReactorSpectrumMultiple :: RandomIsotopes()
{
    Norm =0;//reset
    for (Int_t isotope=0; isotope<NIsotopes; isotope++)
    {
        rand->SetSeed(0);
        RandomIsotopeError[isotope] = (IsotopeError[isotope] * rand->Gaus(0,1));
        IsotopeFrac[isotope]  = (1 + RandomIsotopeError[isotope]) * NominalIsotopeFrac[isotope];
        
        Norm = Norm + IsotopeFrac[isotope];
        std::cout << " Random Isotope" << isotope+1 << " Error is " << RandomIsotopeError[isotope] << "\n";
        std::cout << " Random Isotope" << isotope+1 << " Isotope Frac is " << IsotopeFrac[isotope]  << "\n";
        
    }
    
    //Normalize total isotope fraction sum to 1 after being randomized.
    for (Int_t i = 0; i<NIsotopes; i++)
    {
        IsotopeFrac[i]  = IsotopeFrac[i]/Norm;
        std::cout << " Random Isotope" << i << " Fraction is " << IsotopeFrac[i] << "\n";
    }
}

void ReactorSpectrumMultiple :: RandomPower()
{
    for(Int_t reactor=0; reactor<NReactors;reactor++)
    {
        rand->SetSeed(0);
        RandomReactorPowerError[reactor] = (ReactorPowerError[reactor] * rand->Gaus(0,1));
        ReactorPower[reactor] = (1 + RandomReactorPowerError[reactor]) * NominalReactorPower[reactor];
        std::cout << " Random Power" << reactor+1 << " is " << ReactorPower[reactor] << "\n";
    }
}
