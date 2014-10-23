// Class to calculate Antineutrino cross section given a number of bins, initial and final energies.
//.L InverseBetaCS.C; before using this.
//gROOT->ProcessLine(".x InverseBetaCS.C");//somehow this doesn't work so it has to be called manually in console

// So I ended up copying InverseBetaCS.h inside this class.
//Calculate cross section giving NBins, Initial Energy and Final Energy

// José de Arcos
#pragma once
#include "stdio.h"
#include "TFile.h"
#include "NominalData.h"
#include "TF1.h"

const Int_t gNbins=240-36;
const Double_t gkMassProton = 938.27203;  // MeV
const Double_t gkMassNeutron = 939.56536; // MeV
const Double_t gkMassElectron = 0.51099892; // MeV

class CrossSection
{
private:
    NominalData* Nom;
    Double_t InitialEnergy;
    Double_t FinalVisibleEnergy;
    
    Double_t fEpos; // positron energy (MeV)
    Double_t fEnu;  // anti neutrino energy (MeV)
    
    //utility constants set in constructor...
    Double_t fF;    // hadronic weak interaction constants
    Double_t fG;    // hadronic weak interaction constants
    Double_t fF2;   // hadronic weak interaction constants
    Double_t fCosCab; // cosine of Cabibbo angle
    Double_t fDelta;  // M_n - M_p (MeV)
    Double_t fYsq;   // y defined in text below eqn 11
    Double_t fMassEsq; // electron mass squared (MeV^2)
    Double_t fSigma0; // eqn 8,9
    Double_t fF2Plus3G2;
    Double_t fF2MinusG2;
    
    // variables that are here to save calculation
    Double_t fE0;  //set in Ee1()
    Double_t fP0;  //set in Ee1();
    Double_t fV0;  //set in Ee1();
    Double_t fE1;  //set in DSigDCosTh
    
    TH1D* fHTotalCrossSection;
    
    void setupTotalCrossSection();
public:
    CrossSection();
    CrossSection(NominalData*);
    void CrossCalc();

    Double_t Ee0(Double_t aEnu) { return (aEnu - fDelta); }  // eqn 6 in Vogel/Beacom
    Double_t Ee1(Double_t aEnu, Double_t aCosTheta); // eqn 11
    Double_t DSigDCosTh(Double_t aEnu, Double_t aCosTheta); // eqn 12
    Double_t GammaTerm(Double_t aCosTheta); // eqn 13
    Double_t SigmaTot(Double_t aEnu); // integration of eqn 12 by Gaussian quadrature
    Double_t PromptEnergyToNeutrinoEnergy(Double_t aEprompt); // incles 2*.511
    TF1* fDifferential; //! here for easy integration of diff cross section
    
    Double_t gDSigmaByDCosTheta(Double_t* x, Double_t* a);
    Double_t gSigmaTotal(Double_t* x, Double_t* a);
};

CrossSection* gCrossSection;

CrossSection :: CrossSection()
{
    std::cout << " the default constructor shouldn't be called" << std::endl;
    
    exit(EXIT_FAILURE);
    
    Nom = new NominalData(0,2);
    
    InitialEnergy = Nom->GetEmin();
    FinalVisibleEnergy = Nom->GetEVisMax();
    
    fF = 1;
    fG = 1.26;
    fF2 = 3.706;
    fDelta = gkMassNeutron - gkMassProton;
    fMassEsq = gkMassElectron*gkMassElectron;
    fYsq = (fDelta*fDelta - fMassEsq)/2;
    fF2Plus3G2 = fF*fF + 3*fG*fG;
    fF2MinusG2 = fF*fF - fG*fG;
    fSigma0 = 0.0952/(fF2Plus3G2); // *10^{-42} cm^2, eqn 9
    fDifferential = new TF1("fDifferential", gCrossSection, &CrossSection::gDSigmaByDCosTheta, -1, 1, 1, "CrossSection","gDSigmaByDCosTheta");
    fDifferential->SetParameter(0,3); // default enu = 3 MeV
    fHTotalCrossSection = NULL;
    
    delete Nom;
}

CrossSection :: CrossSection(NominalData* Data)
{
    InitialEnergy = Data->GetEmin();
    FinalVisibleEnergy = Data->GetEVisMax();
    
    fF = 1;
    fG = 1.26;
    fF2 = 3.706;
    fDelta = gkMassNeutron - gkMassProton;
    fMassEsq = gkMassElectron*gkMassElectron;
    fYsq = (fDelta*fDelta - fMassEsq)/2;
    fF2Plus3G2 = fF*fF + 3*fG*fG;
    fF2MinusG2 = fF*fF - fG*fG;
    fSigma0 = 0.0952/(fF2Plus3G2); // *10^{-42} cm^2, eqn 9
    fDifferential = new TF1("fDifferential", gCrossSection, &CrossSection::gDSigmaByDCosTheta, -1, 1, 1, "CrossSection","gDSigmaByDCosTheta");
    fDifferential->SetParameter(0,3); // default enu = 3 MeV
    fHTotalCrossSection = NULL;
}

void CrossSection :: CrossCalc()
{
    Double_t Energies[2];
    Energies[0]= InitialEnergy;//Initial energy for the X axis
    Energies[1]= FinalVisibleEnergy;//Final energy for the X axis
    Double_t SigmaTotal;
    Double_t a[1]={0.0};
    
    TH1D* CrossSection = new TH1D("CrossSection", "#nu cross section", gNbins, Energies[0], Energies[1]);
    
    //Fill histogram with Inverse Beta Cross Section for all energies
    for (int i=1;i<=gNbins;i++)
    {
        Energies[0]=Energies[0]+((FinalVisibleEnergy-InitialEnergy)/(gNbins));//there's some small error due to double decimal handling.
        SigmaTotal = gSigmaTotal(Energies,a);
        CrossSection->SetBinContent(i, SigmaTotal);
        delete fHTotalCrossSection;
    }
    
    //Write histogram in the root file
    CrossSection->GetXaxis()->SetTitle("E_{#nu} [MeV]");
    CrossSection->GetYaxis()->SetTitle("#sigma  [10^{-42} cm^{-2}]");
    
    TFile f("CrossSections/nHCrossSection.root", "recreate");
        CrossSection->Write();
    f.Close();
    
    delete CrossSection;
    delete fDifferential;
    delete gCrossSection;
}

/// in the calculation of the cross section, this must be called first.
Double_t CrossSection :: Ee1(Double_t aEnu, Double_t aCosTheta) {
    Double_t answer;
    fE0 = Ee0(aEnu);
    if( fE0 <= gkMassElectron )
    { // below threshold
        fE0 = 0;
        fP0 = 0;
        answer = 0;
    }
    else
    {
        fP0 = TMath::Sqrt(fE0*fE0 - fMassEsq);
        fV0 = fP0/fE0;
        Double_t sqBracket = 1 - aEnu*(1-fV0*aCosTheta)/gkMassProton;
        answer = fE0*sqBracket - fYsq/gkMassProton;
    }
    return answer;
}

Double_t CrossSection :: DSigDCosTh(Double_t aEnu, Double_t aCosTheta) {
    fE1 = Ee1(aEnu, aCosTheta);
    Double_t answer;
    if( fE1<gkMassElectron )
    {
        answer = 0;
    }
    else
    {
        Double_t pe1 = TMath::Sqrt(fE1*fE1 - fMassEsq);
        Double_t ve1 = pe1/fE1;
        Double_t firstLine = (fSigma0/2) * (fF2Plus3G2 + fF2MinusG2*ve1*aCosTheta)* fE1*pe1;
        Double_t secondLine = fSigma0 * GammaTerm(aCosTheta) * fE0 * fP0 / (2* gkMassProton);
        answer = firstLine - secondLine;
    }
    return answer;
}

Double_t CrossSection :: GammaTerm(Double_t aCosTheta)
{
    Double_t firstLine = (2*fE0+fDelta)*(1-fV0*aCosTheta)-fMassEsq/fE0;
    firstLine *= (2*(fF+fF2)*fG);
    Double_t secondLine = fDelta*(1+fV0*aCosTheta) +fMassEsq/fE0;
    secondLine *= (fF*fF+fG*fG);
    Double_t thirdLine = fF2Plus3G2 *( (fE0+fDelta)*(1-aCosTheta/fV0) - fDelta);
    Double_t fourthLine = (fF*fF - fG*fG)*fV0*aCosTheta;
    fourthLine *= ( (fE0+fDelta)*(1-aCosTheta/fV0) - fDelta );
    
    Double_t answer = firstLine + secondLine + thirdLine + fourthLine;
    return answer;
}

Double_t CrossSection :: SigmaTot(Double_t aEnu)
{
    //  fDifferential->SetParameter(0, aEnu);
    Double_t answer;
    if( aEnu<0 )
    {
        Warning("SigmaTot(Double_t aEnu)","Tried to calculate cross section for -ve nu energy");
        return 0;
    }
    if( aEnu>10 )
    {
        // table on precalculated to 10 MeV
        answer = fDifferential->Integral(-1.0, 1.0, &aEnu, 1e-9);
        return answer;
    }
    if( fHTotalCrossSection==NULL ) setupTotalCrossSection();
    int bin = fHTotalCrossSection->FindBin(aEnu);
    Double_t answer1 =  fHTotalCrossSection->GetBinContent(bin);
    Double_t e1 = fHTotalCrossSection->GetBinCenter(bin);
    Double_t answer2, e2;
    if( bin!=1000 )
    {
        answer2 =  fHTotalCrossSection->GetBinContent(bin+1);
        e2 = fHTotalCrossSection->GetBinCenter(bin+1);
    }
    else
    {
        answer2 =  fHTotalCrossSection->GetBinContent(bin-1);
        e2 = fHTotalCrossSection->GetBinCenter(bin-1);
    }
    answer = answer1 + (answer2-answer1)*( aEnu-e1)/(e2-e1); // answer1 + slope correction

    return answer;
}

void CrossSection :: setupTotalCrossSection()
{
    fHTotalCrossSection = new TH1D("htotalCross","Total #nu (p,n) e Cross Section", 1000, 0, 10);
    Double_t enu;
    for( int i=1 ; i<=1000 ; i++ )
    {
        enu = fHTotalCrossSection->GetBinCenter(i);
        fHTotalCrossSection->SetBinContent(i, fDifferential->Integral(-1.0, 1.0, &enu, 1e-9) );
    }
}


/// takes given prompt energy and turns it into a neutrino energy
/// inverts eqn 11
Double_t CrossSection :: PromptEnergyToNeutrinoEnergy(Double_t aE_prompt)
{
    Double_t ePos = aE_prompt - gkMassElectron; // paper uses total E so only
    // subtract one mass
    Double_t a = -1/(gkMassProton);
    Double_t b = 1+fDelta/gkMassProton;
    Double_t c = -fYsq/gkMassProton - ePos - fDelta;
    Double_t answer1 = ( -b + sqrt(b*b-4*a*c) )/(2*a);
    Double_t answer2 = ( -b - sqrt(b*b-4*a*c) )/(2*a);
    //  printf("%e\t%e\t%e\t%e\t%e\n", aE_prompt, ePos, answer1, answer2, ePos+fDelta);
    return answer1;
}

///global function to allow calls from a TF1.
///This is necessary to use the gaussian quadrature method
///built into TF1
/// dsig/dcos(theta)
Double_t CrossSection :: gDSigmaByDCosTheta(Double_t* x, Double_t* a)
{
    Double_t cosTheta = x[0];
    Double_t enu = a[0];
    if( gCrossSection == NULL ) gCrossSection = new CrossSection();
    return gCrossSection->DSigDCosTh(enu, cosTheta);
    
}

//a global function to allow the total cross section call be put into a TF1
// a redirect to InverseBetaCS::SigmaTot()
Double_t CrossSection ::  gSigmaTotal(Double_t* x, Double_t* a)
{
    // a not used
    // x[0] = neutrino energy (MeV)
    if( gCrossSection == NULL ) gCrossSection = new CrossSection();
    return gCrossSection->SigmaTot(x[0]);
}

