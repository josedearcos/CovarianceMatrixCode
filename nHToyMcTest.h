#pragma once
#include "TRandom3.h"
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF3.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2.h"
#include "TMath.h"
#include "TCanvas.h"

//Particle masses
const Double_t Me = 0.510999; // MeV
const Double_t Mn = 939.565; // MeV
const Double_t Mp = 938.272; // MeV

//M14A IBD events:
//const Int_t TotalIBDGdLs = 4978445;
//const Int_t TotalIBDLs = 5423625;

const Double_t GdLsnHRatio = 0.1570;//+/- 0.80% Doc-DB 9624 M14A IBD Samples
const Double_t GdLsnGdRatio = 0.8417;//+/- 0.80% Doc-DB 9624 M14A IBD Samples

const Double_t LsnHRatio = 0.9580;

const Int_t MatrixBins = 240;
const Int_t NADs = 6;
const Double_t MaxLsRadius = 4;    //3 m vessel for GdLs and 4 m vessel for Ls
const Double_t MaxGdLsRadius = 3;    //3 m vessel for GdLs and 4 m vessel for Ls
const Double_t MaxTheta = 360;   //degrees
const Double_t MaxZ = 2;//meters
const Double_t LsMeanAbsorptionTime = 209e-6;//   In the Ls region
const Double_t GdLsMeanAbsorptionTime = 29e-6;//   In the GdLs region
const Double_t nHpeak[NADs] = {2.3333,2.3399,2.3399,2.3382,2.3483,2.3415};//  MeV
const Double_t nHSigma[NADs] = {0.1439,0.1451,0.1407,0.1416,0.1402,0.1449};// MeV

const Int_t TimeCut = 400; // μs
const Double_t EnergyCut = 1.5; //MeV
const Int_t DistanceCut = 500;//mm

const Double_t LsDistanceCutEfficiency = 0.7527;//Average +/- 0.0035 corr +/- 0.003 uncorr
const Double_t GdLsDistanceCutEfficiency = 0.6562;

typedef struct
{
    Double_t PromptRad;
    Double_t PromptTheta;
    Double_t PromptZ;
    
    Double_t DelayRad;
    Double_t DelayTheta;
    Double_t DelayZ;
    
    Double_t PromptEnergy;
    Double_t DelayEnergy;
    
    Double_t Time;//Between Prompt and Delay.
    
    bool valid;
} Event;

class nHToyMcTest
{
public:
    void nHToyMcTestMain();
private:
    
    TH2D* DelayPositionMap[NADs];
    TH2D* PromptPositionMap[NADs];

    TH2D* nHDelayPromptMap[NADs];
    TH2D* nGdDelayPromptMap[NADs];
    
    TH2D* nHDelayPromptMapUniformCorrected[NADs];
    TH2D* nGdDelayPromptMapUniformCorrected[NADs];
    
    TH2D* nHDelayPromptMapAfterCuts[NADs];
    TH2D* nGdDelayPromptMapAfterCuts[NADs];
    
    TH1D* nHAntineutrinoSpectrumH[NADs];
    TH1D* nGdAntineutrinoSpectrumH[NADs];
    
    TF1* GetNeutrinoToVisibleFunction(bool);
    Double_t VisibleEnergy1F(Double_t*,Double_t*);
    Double_t VisibleEnergy0F(Double_t*,Double_t*);
//    Double_t CorrectionF(Double_t*,Double_t*);
    TF1* VisibleF;
    Int_t TotalBins;
    Double_t InitialDelayEnergy;
    Double_t FinalDelayEnergy;
    void Save();

};

//nHToyMcTest :: nHToyMcTest()
//{
//
//}
void nHToyMcTest::nHToyMcTestMain()
{
    TRandom3* rand = new TRandom3(0);

//    TF3* CorrF = new TF3("CorrF",this,&nHToyMcTest::CorrectionF,0,20,0,5,-2.5,2.5,6,"nHToyMcTest","CorrectionF");
//    CorrF->SetParameters(8.05,7.84628,0.0341294,0.0121750,0.0164275,0.000733006);//from Doc-DB 9441
    
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        TFile* File = TFile::Open("./RootOutputs/NominalOutputs/Oscillation.root");
        
        Char_t histnameTotal[100];//47
        
        sprintf(histnameTotal,"Total spectrum after oscillation at AD%i", AD+1);
        
        File->cd("Total AD Spectra after oscillation");
        nGdAntineutrinoSpectrumH[AD] = (TH1D*)gDirectory->Get(histnameTotal);//For Gd

        nHAntineutrinoSpectrumH[AD] = (TH1D*)gDirectory->Get(histnameTotal);//For Gd!! Need nH absorption cross section
//        nHAntineutrinoSpectrumH[AD]->Scale((TotalIBDGdLs*GdLsnHRatio+TotalIBDLs*LsnHRatio)/nHAntineutrinoSpectrumH[AD]->Integral());// Normalize shape and then multiply by the number of events
        
        Int_t TotalnGDEvents = nGdAntineutrinoSpectrumH[AD]->Integral();
        Int_t TotalnHEvents = nHAntineutrinoSpectrumH[AD]->Integral();
        
        std::cout << TotalnGDEvents << " "<< TotalnHEvents << std::endl;
        Int_t MaxPromptBins = nHAntineutrinoSpectrumH[AD]->GetXaxis()->GetNbins();
        Double_t InitialPromptEnergy = nHAntineutrinoSpectrumH[AD]->GetXaxis()->GetXmin();
        Double_t FinalPromptEnergy = nHAntineutrinoSpectrumH[AD]->GetXaxis()->GetXmax();
        
        File->Close();

        TotalBins = 240;
        InitialDelayEnergy = 0;
        FinalDelayEnergy = 12;
        
        Double_t BinWidthDelay=(FinalDelayEnergy-InitialDelayEnergy)/TotalBins;
        Double_t BinWidthPrompt=(FinalPromptEnergy-InitialPromptEnergy)/MaxPromptBins;
        
        nHDelayPromptMap[AD] = new TH2D(Form("PromptVsDelay%d",AD+1),Form("Prompt Vs Delay%d",AD+1), MaxPromptBins,InitialPromptEnergy,FinalPromptEnergy,TotalBins,InitialDelayEnergy,FinalDelayEnergy);
        nHDelayPromptMapUniformCorrected[AD] = new TH2D(Form("PromptVsDelayUnifCorrected%d",AD+1),Form("Prompt Vs Delay Unif Corrected%d",AD+1), MaxPromptBins,InitialPromptEnergy,FinalPromptEnergy,TotalBins,InitialDelayEnergy,FinalDelayEnergy);
        nHDelayPromptMapAfterCuts[AD] = new TH2D(Form("PromptVsDelayAfterCuts%d",AD+1),Form("Prompt Vs Delay After Cuts%d",AD+1), MaxPromptBins,InitialPromptEnergy,FinalPromptEnergy,TotalBins,InitialDelayEnergy,FinalDelayEnergy);
        DelayPositionMap[AD] = new TH2D(Form("DelayPositionMap%d",AD+1),Form("DelayPositionMap%d",AD+1), 200, 0,4,200, -2,2);
        PromptPositionMap[AD] = new TH2D(Form("PromptPositionMap%d",AD+1),Form("PromptPositionMap%d",AD+1), 200, 0,4,200, -2,2);

        for(Int_t i = 0; i<TotalBins;i++)
        {
            Int_t NumbernHNeutrinos = nHAntineutrinoSpectrumH[AD]->GetBinContent(i+1);
            
            std::vector<Event> Evento(NumbernHNeutrinos);
            std::vector<Event> UniformCorrected(NumbernHNeutrinos);
            std::vector<Event> KinematicCorrected(NumbernHNeutrinos);
            
            for(Int_t j = 0; j<NumbernHNeutrinos; j++)
            {
                //Maybe I need to generate the distance as a exp distribution with a mean distance, I need to follow a distance distribution that agrees with the cut efficiency!
                
                Double_t TotalGdLsEvents = TotalnGDEvents/GdLsnGdRatio;
                Double_t TotalGdLsnHEvents = GdLsnHRatio*TotalGdLsEvents;
                Double_t TotalLsnHEvents = TotalnHEvents-TotalGdLsnHEvents;
                
                Double_t ProbnHLs =TotalLsnHEvents/TotalnHEvents;
                Double_t ProbnHGdLs =TotalGdLsnHEvents/TotalnHEvents;

                if(rand->Uniform()<ProbnHLs)//LS region
                {
                    if(rand->Integer(7)>2)// random 3-6, give 4/7% chances, same as area covered by the region (3 square meters vs 4 square meters)
                    {
                        Evento[j].PromptRad = MaxLsRadius*rand->Uniform(0.75,1);//Inside LS region
                        Evento[j].PromptZ = MaxZ*rand->Uniform(-1,1);//Inside LS region
                    }
                    else
                    {
                        Evento[j].PromptRad = MaxLsRadius*rand->Uniform(0,0.75);//Inside LS region
                        Evento[j].PromptZ = TMath::Power(-1,Int_t(rand->Integer(2)))*MaxZ*rand->Uniform(0.75,1);//Inside LS region
                    }
                    
                    Evento[j].PromptTheta = MaxTheta*rand->Uniform();
                    Evento[j].DelayTheta = MaxTheta*rand->Uniform();
                    
                    if(rand->Integer(7)>2)
                    {
                        Evento[j].DelayRad = MaxLsRadius*rand->Uniform(0.75,1);//Inside LS region
                        Evento[j].DelayZ = MaxZ*rand->Uniform(-1,1);//Inside LS region
                    }
                    else
                    {
                        Evento[j].DelayRad = MaxLsRadius*rand->Uniform(0,0.75);//Inside LS region
                        Evento[j].DelayZ = TMath::Power(-1,Int_t(rand->Integer(2)))*MaxZ*rand->Uniform(0.75,1);//Inside LS region
                    }
                }
                else
                {
                    Evento[j].PromptRad = MaxLsRadius*rand->Uniform(0,0.75);//Inside GDLS region
                    Evento[j].PromptZ = MaxZ*rand->Uniform(-0.75,0.75);//Inside GDLS region
                    Evento[j].DelayRad = MaxLsRadius*rand->Uniform(0,0.75);//Inside GDLS region
                    Evento[j].DelayZ = MaxZ*rand->Uniform(-0.75,0.75);//Inside GDLS region
                    Evento[j].PromptTheta = MaxTheta*rand->Uniform();
                    Evento[j].DelayTheta = MaxTheta*rand->Uniform();

                }
//                std::cout << Int_t(Evento[j].DelayRad*200/4) << std::endl;
//                std::cout << Int_t((Evento[j].DelayZ+2)*200/4) << std::endl;
//                std::cout << Int_t(Evento[j].PromptRad*200/4) << std::endl;
//                std::cout << Int_t((Evento[j].PromptZ+2)*200/4) << std::endl;

                DelayPositionMap[AD]->SetBinContent(Int_t(Evento[j].DelayRad*200/4)+1,Int_t((Evento[j].DelayZ+2)*200/4)+1,1+DelayPositionMap[AD]->GetBinContent(Int_t(Evento[j].DelayRad*200/4)+1,Int_t((Evento[j].DelayZ+2)*200/4)+1));
                
                PromptPositionMap[AD]->SetBinContent(Int_t((Evento[j].PromptRad)*200/4)+1,Int_t((Evento[j].PromptZ+2)*200/4)+1,1+PromptPositionMap[AD]->GetBinContent(Int_t((Evento[j].PromptRad)*200/4)+1,Int_t((Evento[j].PromptZ+2)*200/4)+1));
                
                Evento[j].PromptEnergy = nHAntineutrinoSpectrumH[AD]->GetBinCenter(i);
                Evento[j].DelayEnergy = rand->Gaus(nHpeak[AD],nHSigma[AD]);//   This should follow the tail shape function for each energy
                
                Evento[j].Time = rand->Exp(LsMeanAbsorptionTime);
                
                Int_t DelayEnergyBin = (Int_t)(Evento[j].DelayEnergy*TotalBins/(FinalDelayEnergy-InitialDelayEnergy));
                
//                std::cout << i << " " << DelayEnergyBin << " " << j << std::endl;
                
                nHDelayPromptMap[AD]->SetBinContent(i+1,DelayEnergyBin+1,1+nHDelayPromptMap[AD]->GetBinContent(i+1,DelayEnergyBin+1));
                
                UniformCorrected[j].PromptRad = Evento[j].PromptRad;
                UniformCorrected[j].PromptZ = Evento[j].PromptZ;
                UniformCorrected[j].DelayRad = Evento[j].DelayRad;
                UniformCorrected[j].DelayZ = Evento[j].DelayZ;
                UniformCorrected[j].Time = Evento[j].Time;
                
//                UniformCorrected[j].PromptEnergy = CorrF->Eval(Evento[j].PromptEnergy,Evento[j].PromptRad,Evento[j].PromptZ);
//                UniformCorrected[j].DelayEnergy = CorrF->Eval(Evento[j].DelayEnergy,Evento[j].DelayRad,Evento[j].DelayZ);
                UniformCorrected[j].valid=1;
                
                Int_t DelayEnergyBinUniformCorrected = (Int_t)(UniformCorrected[j].DelayEnergy/BinWidthDelay);
                
                nHDelayPromptMapUniformCorrected[AD]->SetBinContent(i+1,DelayEnergyBinUniformCorrected+1,1+nHDelayPromptMapUniformCorrected[AD]->GetBinContent(i+1,DelayEnergyBinUniformCorrected+1));
                //  Cuts:
                //  Low Energy Cut for prompt signal
                if(UniformCorrected[j].DelayEnergy<EnergyCut)
                {
                    UniformCorrected[j].valid = 0;
                }
                //  Time cut:
                if(UniformCorrected[j].Time>TimeCut)
                {
                    UniformCorrected[j].valid = 0;
                }
                //  Distance cut:
                //                if(UniformCorrected[j].Distance>DistanceCut)
                //                {
                //                    UniformCorrected[j].valid = 0;
                //                }
                
                //  3σ cut:
                if(TMath::Abs(UniformCorrected[j].DelayEnergy-nHpeak[AD])>3*nHSigma[AD])
                {
                    UniformCorrected[j].valid = 0;
                }
                
                //  Fill histogram after cuts:
                if(UniformCorrected[j].valid == 1)
                {
                    nHDelayPromptMapAfterCuts[AD]->SetBinContent(i+1,DelayEnergyBinUniformCorrected+1,1+nHDelayPromptMapAfterCuts[AD]->GetBinContent(i+1,DelayEnergyBinUniformCorrected+1));
//                    std::cout << DelayEnergyBinUniformCorrected << std::endl;
                }
                
                //  Now apply distortions (IAV, NL, Resolution)
                
                //  Kinematic Energy Shift:
                
                GetNeutrinoToVisibleFunction(1);//0 for 0th order, 1 for 1st order
                
                if(UniformCorrected[j].valid == 1)
                {
                    KinematicCorrected[j].PromptEnergy = VisibleF->Eval(UniformCorrected[j].PromptEnergy);
                    KinematicCorrected[j].DelayEnergy = KinematicCorrected[j].DelayEnergy+(KinematicCorrected[j].PromptEnergy-KinematicCorrected[j].PromptEnergy);//The neutron absorbs the difference
                }
                
                delete VisibleF;
                
                if (UniformCorrected[j].PromptEnergy <= 1.8)
                {
                    KinematicCorrected[j].PromptEnergy = 0;//There's no enough energy for IBD to happen.
                }
                
                //                KinematicCorrectedVector->SetBinContent(j+1,KinematicCorrected[j].PromptEnergy);//TH1D
            }
        }
    }
    
    //  IAV in nH is roughly 3 times the IAV effect in nGd
    //IAV
    
    //Calculate IAV Shift
    
    //    IAVCorrected[j].PromptEnergy +=
    //    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TrueEnergyIndex+1; VisibleEnergyIndex++)
    //    {
    //        OscDeltaIAVSpectrumH[TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, OscDeltaIAVSpectrumH[TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1) + 3*IAVMatrix[far+ADsEH1+ADsEH2][PositronEnergyIndex][VisibleEnergyIndex] * OscDeltaPositronSpectrumH[TrueEnergyIndex]->GetBinContent(PositronEnergyIndex+1));
    //    }
    //
    
    //  Now Sum all events!
    
    //    if(Evento[j].Rad<MaxGdLsRadius)
    //    {
    //        Evento[j].Time = rand->Exp(GdLsMeanAbsorptionTime);
    //    }
    Save();
    
}

void nHToyMcTest :: Save()
{
    TFile* SaveF = new TFile("nHToyTest.root","recreate");
    TCanvas* c = new TCanvas("PromptVsDelay","PromptVsDelay",800,600);
    c->Divide(3,2);
    for (Int_t AD = 0; AD<NADs; AD++)
    {
        c->cd(AD+1);

        nHDelayPromptMap[AD]->Draw("colz");
        nHDelayPromptMap[AD]->Write();
        nHDelayPromptMapUniformCorrected[AD]->Write();
        nHDelayPromptMapAfterCuts[AD]->Write();
        nHAntineutrinoSpectrumH[AD]->Write();
        DelayPositionMap[AD]->Write();
        PromptPositionMap[AD]->Write();
    }
    c->Print("./Images/nHPromptVsDelay.eps", "png");
    c->Close();

        TCanvas* c2 = new TCanvas("EventMap");
        c2->Divide(2,1);
        c2->cd(1);
        DelayPositionMap[0]->Draw("colz");
        c2->cd(2);
        PromptPositionMap[0]->Draw("colz");
        c2->Print("./Images/nHADMAP.eps", "png");
        c2->Close();
    
    SaveF->Close();
}
TF1* nHToyMcTest :: GetNeutrinoToVisibleFunction(bool order)
{
    if(order==0)
    {
        Double_t Correction = Mn-Mp-Me;
        VisibleF = new TF1("VisibleF",this,&nHToyMcTest::VisibleEnergy0F,InitialDelayEnergy,FinalDelayEnergy,1,"nHToyMcTest","VisibleEnergy0F");
        VisibleF->SetParameter(0,Correction);
    }
    if(order==1)
    {
        VisibleF = new TF1("VisibleF",this,&nHToyMcTest::VisibleEnergy1F,InitialDelayEnergy,FinalDelayEnergy,1,"nHToyMcTest","VisibleEnergy1F");
        VisibleF->SetParameter(0,0);//Maybe I can improve this using a MC simulation of the angular distribution calculated in http://authors.library.caltech.edu/2796/1/VOGprd99.pdf
    }
    return VisibleF;
}

Double_t nHToyMcTest :: VisibleEnergy1F(Double_t* energ, Double_t* par)
{
    Double_t Enu = energ[0];
    
    Double_t Delta = Mn-Mp;
    
    Double_t Ee0 = Enu - Delta;
    
    Double_t gamma0 = Ee0/Me;
    Double_t v0 = sqrt(1 - 1/gamma0/gamma0);
    
    Double_t costheta = -0.034*v0+2.4*Enu/Mp;
    
    Double_t y2 = (Delta*Delta - Me*Me)/2.;
    
    Double_t Ee1 = Ee0 * (1 - Enu/Mp*(1.0 - v0*costheta)) - y2/Mp;
    
    return Ee1 + Me;
}

Double_t nHToyMcTest :: VisibleEnergy0F(Double_t* energ, Double_t* par)
{
    return energ[0]-par[0];
}

//Double_t nHToyMcTest :: CorrectionF(Double_t* input,Double_t* par)
//{
//    //Apply Eraw Erec correction here
//    // 8.05/(7.84628*(1+0.0341294*R2)*(1−0.0121750*Z −0.0164275*Z2 +7.33006e−04*Z3))
//    Double_t e_final;
//    Double_t e_orig = input[0];
//    Double_t R = input[1];
//    Double_t Z = input[2];
//    
//    e_final = e_orig*(par[0]/(par[1]*(1+par[2]*TMath::Power(R,2))*(1-par[3]*Z-par[4]*TMath::Power(Z,2)+par[5]*TMath::Power(Z,3))));
//    
//    return e_final;
//}
//generate random energies following antineutrino spectra
//apply non-uniformityusing standarized detector spots
//gaussian around that energy to simulate tail effect
// apply distortions, fee and scintillator non linearities are applied at different stages
// sum all the energies after being distorted to get the final spectrum
