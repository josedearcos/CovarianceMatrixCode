#include "TFile.h"
#include "TH1F.h"
#include <math.h>

#include "NominalData.h"

const Int_t NADs = 6; //Change to 8 if I get to have 8 ADs backgrounds
const Int_t BGNDs = 8; //LiHe is not included yet

using namespace std;

class FitBackgrounds
{
private:
    

    //Histograms
    TH1F* AccidentalsH[NADs];
    TH1F* AmCH;
    TH1F* BackgroundsH[BGNDs];
    TH1F* RebinBackgroundsH[BGNDs];
    TH1F* RebinHist;

    //Files
    TFile* BackgroundsF;
    TFile* AccidentalsF[NADs];
    TFile* FastNeutronF;
    TFile* AmCF;
    
    //Binning parameters:
    Int_t Nbins;
    Int_t InitialBins;
    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    
    void FitAmc();
    void FitFastNeutrons();
    
public:
    
    FitBackgrounds();
    FitBackgrounds(NominalData Data);
    void ReadBackgrounds();
};

FitBackgrounds::FitBackgrounds()
{
    Nom = new NominalData();
    
    //Binning variables
    Nbins = Nom->GetNbins();
    InitialBins = Nom->GetInitialBins();
    InitialEnergy = Nom->GetEmin();
    FinalEnergy = Nom->GetEmax();
    InitialVisibleEnergy = Nom->.GetEVisMin();
    FinalVisibleEnergy = Nom->GetEVisMax();
}

FitBackgrounds::FitBackgrounds(NominalData Data)
{
    ran = new TRandom3();
    
    Nbins = Data.GetNbins();
    InitialBins = Data.GetInitialBins();
    InitialEnergy = Data.GetEmin();
    FinalEnergy = Data.GetEmax();
    InitialVisibleEnergy = Data.GetEVisMin();
    FinalVisibleEnergy = Data.GetEVisMax();
}

//This works for the current given Backgrounds and 6 ADs, will have to be changed in the case of 8 ADs or diverse inputs.
void FitBackgrounds::ReadBackgrounds()
{
    Char_t filenameAccidentals[1024];
    
    Int_t j=0;

    for(Int_t AD = 0; AD < NADs; AD++)
    {
        sprintf(filenameAccidentals,"./BackgroundSpectrum/Bkg_Accidental_AD%i.root", AD+1);

        AccidentalsF[AD] = TFile::Open(filenameAccidentals);
        
        if(AD==2)
        {
            AccidentalsH[2]=(TH1F*)gDirectory->Get("h1dAccBkg22_1");
        }
        if(AD<2)
        {
            AccidentalsH[AD]=(TH1F*)gDirectory->Get(Form("h1dAccBkg22_%i",j+1));
            j++;
        }
        if(AD>2)
        {
            AccidentalsH[AD]=(TH1F*)gDirectory->Get(Form("h1dAccBkg22_%i",j-1));
            j++;
        }
        AccidentalsH[AD]->SetDirectory(0);

        BackgroundsH[AD]=(TH1F*)AccidentalsH[AD]->Clone(Form("Accidentals_AD%i",AD+1));

        BackgroundsH[AD]->Reset();

        // cout<<"NPoints "<<AccidentalsH[AD]->GetXaxis()->GetNbins();
        for (Int_t pts=1; pts <= AccidentalsH[AD]->GetXaxis()->GetNbins(); pts++)
        {
            BackgroundsH[AD]->SetBinContent(pts,AccidentalsH[AD]->GetBinContent(pts));
        }
        BackgroundsH[AD]->SetDirectory(0);

        AccidentalsF[AD]->Close();
    }
    
    FitAmc();//Exponential fit

    FitFastNeutrons();//Pol0 and Pol1 fit
    
    //9Li8He missing
    
    BackgroundsF = TFile::Open("./BackgroundSpectrum/Backgrounds.root","recreate");

    for(Int_t BGND = 0; BGND < BGNDs; BGND++)
    {
        BackgroundsH[BGND]->Rebin(BackgroundsH[BGND]->GetXaxis()->GetNbins()/(Nbins+InitialBins));
        BackgroundsH[BGND]->Write();
    }
    BackgroundsF->Close();
}

void FitBackgrounds::FitAmc()
{
    AmCF = TFile::Open("./BackgroundSpectrum/Bkg_StrongAmC.root");
    AmCH=(TH1F*)gDirectory->Get("hist_Bkg_StrongAmC");
    AmCH->SetDirectory(0);
    AmCF->Close();

    AmCF = TFile::Open("./BackgroundSpectrum/Backgrounds.root","update");
    BackgroundsH[NADs]=(TH1F*)AmCH->Clone();
    BackgroundsH[NADs]->Reset();
    BackgroundsH[NADs]->SetDirectory(0);
    AmCF->Close();

    for(Int_t i=1;i<=AmCH->GetXaxis()->GetNbins();i++)
    {
        if((AmCH->GetBinContent(i))<0)
        {
            AmCH->SetBinContent(i,-1*(AmCH->GetBinContent(i))); //Invert negative bins??
            //                          OR
           // AmCH->SetBinContent(i,0);//Set negative bins to 0?
        }
    }
//    I need to study what are the implications. Right now I'm not applying any correction to the original histogram.
    
    Double_t IntegralOriginal = AmCH->Integral();
//    cout << IntegralOriginal<<"\n";
    
    TF1 *AmCFunc = new TF1("AmCFunc","exp([0]+x*[1])",InitialVisibleEnergy,FinalVisibleEnergy);
    AmCH->Fit("AmCFunc","Q0");//"Q0 is quiet mode and no draw
    // Try to change FIT options from 1.5-12? Maybe improves the output.

    ///////////////////////////////////////////////////
    //Should I fit this with or without errors? (CHECK) (Original histogram (errors+content) vs only GetBinContent)???? What about negative bins?
    ///////////////////////////////////////////////////

    BackgroundsH[NADs]->Add(AmCFunc);
    Double_t IntegralExponential = BackgroundsH[NADs]->Integral(InitialVisibleEnergy,FinalVisibleEnergy);
//    cout << IntegralExponential <<"\n";

    //Normalize histogram with respect to integral of the original graph, I can use data once I have it and normalize it to the expected number of AmC events.
    BackgroundsH[NADs]->Scale(IntegralOriginal/IntegralExponential);
   // BackgroundsH[NADs]->Draw("same");

    //cout << IntegralOriginal/IntegralExponential <<"\n";

}

void FitBackgrounds::FitFastNeutrons()
{
    FastNeutronF = TFile::Open("./BackgroundSpectrum/Bkg_FastNeutron_EH1_MC.root");
    FastNeutronH=(TH1F*)gDirectory->Get("hist_FastNeutron");
    FastNeutronH->SetDirectory(0);
    
    FastNeutronH->Fit("pol0","Q0","",10,100);
    TF1* FitFNResult = FastNeutronH->GetFunction("pol0");

    Double_t AveragePol0 = FitFNResult->GetParameter(0);
    FastNeutronH->Fit("pol1","+Q0","",10,100);
    Double_t ConstantPol1 = pol1->GetParameter(0);
    Double_t SlopePol1 = pol1->GetParameter(1);
    Double_t AveragePol1=ConstantPol1+SlopePol1*55; //Value of the function in the middle point of the range 55 = ((100-10)/2)
    Double_t TotalAverage = (AveragePol0 + AveragePol1)/2;
    
    BackgroundsH[7] = new TH1F("FN","FN",Nbins+InitialBins,InitialVisibleEnergy,FinalVisibleEnergy);//60 so the bins are 0.2
    
    for(Int_t i=1;i<=Nbins+InitialBins;i++)
    {
        BackgroundsH[7]->SetBinContent(i,TotalAverage);
    }
    
    BackgroundsH[7]->SetDirectory(0);
    FastNeutronF->Close();
}
