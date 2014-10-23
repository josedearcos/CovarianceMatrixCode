#include <TApplication.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TRootEmbeddedCanvas.h>
#include "FitterGui.h"
#include "GuiLinkDef.h"
#include <stdio.h>
#include "FitBackgrounds2.h"
#include "CovarianceMatrix3.h"
#include "NominalData.h"
#include "TColor.h"
#include "TStyle.h"
#include "CrossSection.h"
#include <stdlib.h>
#include "Fitter.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TImage.h"
#include "CreateEnergyMatrix.h"
#include "TSystemDirectory.h"
#include "TPaletteAxis.h"
#include "TLegend.h"
#include "TKey.h"
#include <assert.h> // Used to debug, if you don't want it to run, uncomment the following line:

//#define NDEBUG
const bool TestExternalInputs=0;
const bool PrintOnConsole=0;

//To interpolate the chi2 curves
const Int_t InterpolationFactor = 100;
const Int_t SetStats = 1111001;//neimr or 1111001 to show integral in statbox, 0 to don't show stats in the plots

TH1D* RatioH[MaximumSamples][MaxFarDetectors][MaxNearDetectors];
TH1D* SpectrumH[MaximumSamples][MaxFarDetectors][MaxNearDetectors];
TH1D* NominalSpectrumH[MaxFarDetectors][MaxNearDetectors];

streambuf * coutstream;

int main(int argc, char **argv)
{
    
    //To avoid printing on console and make operations faster,
    
    //To draw using a better palette:
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    gStyle->SetOptStat(SetStats);
    gStyle->SetTitleOffset(1.3);
    gStyle->SetTitleSize(0.03);
    
    TApplication theApp("App",&argc,argv);
    new FitterGui(gClient->GetRoot(),200,200);
    
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);
    
    gROOT->SetBatch(kTRUE);//This avoids canvas popping up every time.
    theApp.Run();
    return 0;
}

FitterGui::FitterGui(const TGWindow *p,UInt_t w,UInt_t h)
{
    gBenchmark = new TBenchmark();
    
    //Write "" if you want to generate the inputs
    if(TestExternalInputs)
    {
        ToyMCSampleDirectory = "/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/ToyMCTrees/nGdInputs/toySpectra_allsys_and_stat.root";
        NominalPredictionsDirectory = "/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/ToyMCTrees/nGdInputs/toySpectra_allsys_and_stat.root";
        ResponseMatrixDirectory = "/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/matrix_evis_to_enu_unified_p12e_unblinded.root";
        SysCovarianceMatrixDirectory = "/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/CovarianceMatrices/Berkeley/matrix_sigsys_full.root";
        BkgCovarianceMatrixDirectory = "/Users/royal/DayaBay/jdearcos/CovarianceMatrixCode/CovarianceMatrices/Berkeley/matrix_bgsys_full.root";
    }
    else
    {
        ToyMCSampleDirectory =  "";
        NominalPredictionsDirectory = "";
        ResponseMatrixDirectory = "";
        SysCovarianceMatrixDirectory = "";
        BkgCovarianceMatrixDirectory = "";
    }
    
    Minuit = 0;//Choose between applying Minuit or manual grid fitting.
    Fit2D = 0; // 1 for 2D Fit, 0 for 1D Fit
    FitSin22t13 = 1; //  1 for Sin22t13 Fit, 0 for DM2ee Fit.
    NFits = 101;//101 in the final version
    Period=1;
    NReactorPeriods=20;
    PlotBin=0;//Initial bin
    Binning = 0;
    NADs=6;
    ToyMC=1;//1 ToyMC, 0 Data. To test the fitter use 1, to fit real data use 0.
    deleteFlag=0;
    deleteFlagSpec=0;
    NearTrueIndex = 0;
    flagCombine = 1;//Show Combine Plot as default
    NominalString = "Nominal";
    CovString = "";
    Analysis = 0; //0 for Gd, 1 for H
    if(Analysis)
    {
        AnalysisString = "Hydrogen";
        
    }
    else
    {
        AnalysisString = "Gadolinium";
    }
    NSamples = 10; //change default to 500 once debug is completed
    CombineMode = 2; // change default to 2 once debug is completed
    NL[0]=0;//BCW
    NL[1]=0;//LBNL
    NL[2]=1;//Unified (Default)
    
    Automatic = 0;//Default
    VaryAccidentalMatrix=0;
    VaryLiHeMatrix=0;
    VaryFastNeutronsMatrix=0;
    VaryAmCMatrix=0;
    DistortLiHeMatrix=0;
    DistortFastNeutronsMatrix=0;
    DistortAmCMatrix=0;
    //Systematics
    IsotopeMatrix=0;
    ReactorPowerMatrix=0;
    EnergyScaleMatrix=0;
    //         EnergyOffsetMatrix=0;
    //         AbsoluteScaleMatrix=0;
    //         AbsoluteOffsetMatrix=0;
    IAVMatrix=0;
    NLMatrix=0;
    ResolutionMatrix=0;
    Sin22t12Matrix=0;
    EfficiencyMatrix=0;
    
    UseToyMCTree = 0;
    PlotCovariance= 1;
    //    // Create a main frame
    fMain = new TGMainFrame(p,w,h);
    // Create canvas widget
    fMaincanvas = new TRootEmbeddedCanvas("Maincanvas",fMain,250,250);
    TCanvas* Maincanvas = fMaincanvas->GetCanvas();
    TImage *img = TImage::Open("Images/dyb_logo_shadow.jpg");
    img->SetConstRatio(kFALSE);
    img->Draw();
    Maincanvas->Update();
    fMain->AddFrame(fMaincanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY, 10,10,10,1));
    // Create a horizontal frame widget with buttons
    
    TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain,250,40);
    
    TGTextButton *inputs = new TGTextButton(hframe,"&Load Inputs");
    inputs->Connect("Clicked()","FitterGui",this,"DoReadInputs()");
    hframe->AddFrame(inputs, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    TGTextButton *toyMC = new TGTextButton(hframe,"&ToyMC");
    toyMC->Connect("Clicked()","FitterGui",this,"DoToyMC()");
    hframe->AddFrame(toyMC, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    TGTextButton *fitter = new TGTextButton(hframe,"&Fitter");
    fitter->Connect("Clicked()","FitterGui",this,"DoFitter()");
    hframe->AddFrame(fitter, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    TGTextButton *exit = new TGTextButton(hframe,"&Exit", "gApplication->Terminate(0)");
    hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    
    // Set a name to the main frame
    fMain->SetWindowName(" Main Menu ");
    
    // Map all subwindows of main frame
    fMain->MapSubwindows();
    
    // Initialize the layout algorithm
    fMain->Resize(fMain->GetDefaultSize());
    
    // Map main frame
    fMain->MapWindow();
}

void FitterGui::DoFitter()
{
    // main frame
    TGMainFrame *fRunFitterFrame = new TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame);
    fRunFitterFrame->SetName("fRunFitterFrame");
    fRunFitterFrame->SetLayoutBroken(kTRUE);
    fRunFitterFrame->SetWindowName(" Fitter Menu ");
    
    // composite frame
    TGCompositeFrame *fFitterFrame = new TGCompositeFrame(fRunFitterFrame,849,587,kVerticalFrame);
    fFitterFrame->SetName("fFitterFrame");
    fFitterFrame->SetLayoutBroken(kTRUE);
    
    ULong_t ucolor;        // will reflect user color changes
    gClient->GetColorByName("#ffffff",ucolor);
    
    // composite frame
    fFitterFrame2 = new TGCompositeFrame(fFitterFrame,490,372,kVerticalFrame,ucolor);
    fFitterFrame2->SetName("fFitterFrame2");
    fFitterFrame2->SetLayoutBroken(kTRUE);
    
    TGFont *ufont;         // will reflect user font changes
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    TGGC   *uGC;           // will reflect user GC changes
    // graphics context changes
    GCValues_t vall717;
    vall717.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall717.fForeground);
    gClient->GetColorByName("#e0e0e0",vall717.fBackground);
    vall717.fFillStyle = kFillSolid;
    vall717.fFont = ufont->GetFontHandle();
    vall717.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall717, kTRUE);
    
    gClient->GetColorByName("#ffffff",ucolor);
    TGLabel *fLabel717 = new TGLabel(fFitterFrame2,"nH/nGd",uGC->GetGC(),ufont->GetFontStruct(),kChildFrame,ucolor);
    fLabel717->SetTextJustify(36);
    fLabel717->SetMargins(0,0,0,0);
    fLabel717->SetWrapLength(-1);
    fFitterFrame2->AddFrame(fLabel717, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel717->MoveResize(8,0,121,18);
    
    // Fitter box
    
    FitterBox = new TGButtonGroup(fFitterFrame2,"Fitter Mode",kVerticalFrame,TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),0xffffff);
    FitterBox->SetName("FitterBox");
    fF[0] = new TGRadioButton(FitterBox,"1D Fit",0);
    fF[1] = new TGRadioButton(FitterBox,"2D Fit",1);
    
    fF[0]->ChangeBackground(ucolor);
    fF[1]->ChangeBackground(ucolor);
    
    fF[0]->SetState(kButtonDown);
    
    for (Int_t i = 0; i < 2; ++i)
    {
        fF[i]->Connect("Pressed()", "FitterGui",  this, "DoFitterMode()");
    }
    FitterBox->Show();
    
    fFitterFrame2->AddFrame(FitterBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    FitterBox->MoveResize(182,15,150,60);
    
    // Fitter box
    Fitter1DBox = new TGButtonGroup(fFitterFrame2,"Oscillation Parameter",kVerticalFrame,TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),0xffffff);
    Fitter1DBox->SetName("FitterBox");
    f1DF[0] = new TGRadioButton(Fitter1DBox,"Sin^2(2theta13) Fit",0);
    f1DF[1] = new TGRadioButton(Fitter1DBox,"DM^2ee Fit",1);
    
    f1DF[0]->ChangeBackground(ucolor);
    f1DF[1]->ChangeBackground(ucolor);
    
    f1DF[0]->SetState(kButtonDown);
    
    for (Int_t i = 0; i < 2; ++i)
    {
        f1DF[i]->Connect("Pressed()", "FitterGui",  this, "Do1DFitterMode()");
    }
    fFitterFrame2->AddFrame(Fitter1DBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    Fitter1DBox->MoveResize(182,80,150,60);
    
    
    // NL box
    
    NLBox = new TGButtonGroup(fFitterFrame2,"Non Linearity Model",kVerticalFrame,TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),0xffffff);
    NLBox->SetName("NLBox");
    fR[0] = new TGRadioButton(NLBox,"BCW Model",0);
    fR[1] = new TGRadioButton(NLBox,"LBNL Model",1);
    fR[2] = new TGRadioButton(NLBox,"Unified Model",2);
    
    fR[0]->ChangeBackground(ucolor);
    fR[1]->ChangeBackground(ucolor);
    fR[2]->ChangeBackground(ucolor);
    
    fR[2]->SetState(kButtonDown);
    
    for (Int_t i = 0; i < 3; ++i)
    {
        fR[i]->Connect("Pressed()", "FitterGui",  this, "DoNLModel()");
    }
    NLBox->Show();
    
    fFitterFrame2->AddFrame(NLBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    NLBox->MoveResize(8,56,150,90);
    
    gClient->GetColorByName("#ffffff",ucolor);
    
    // combo box
    HGdBox = new TGComboBox(fFitterFrame2,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    HGdBox->SetName("HGdBox");
    HGdBox->AddEntry("Gadollinium Analysis",0);
    HGdBox->AddEntry("Hydrogen Analysis",1);
    HGdBox->Resize(152,22);
    HGdBox->Connect("Selected(Int_t)", "FitterGui", this, "DoAnalysis()");
    HGdBox->Select(0);//Gd Default right now.
    fFitterFrame2->AddFrame(HGdBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    HGdBox->MoveResize(8,24,152,22);
    
    gClient->GetColorByName("#ffffff",ucolor);
    TGLabel *binninglabel = new TGLabel(fFitterFrame2,"Binning",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
    binninglabel->SetTextJustify(36);
    binninglabel->SetMargins(0,0,0,0);
    binninglabel->SetWrapLength(-1);
    fFitterFrame2->AddFrame(binninglabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    binninglabel->MoveResize(390,108,128,18);
    
    BinningBox = new TGComboBox(fFitterFrame2,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    BinningBox->SetName("BinningBox");
    BinningBox->AddEntry("LBNL Binning",0);
    BinningBox->AddEntry("Linear Binning",1);
    BinningBox->Connect("Selected(Int_t)", "FitterGui", this, "DoBinning()");
    BinningBox->Select(0);//LBNL Binning default for now
    fFitterFrame2->AddFrame(BinningBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    BinningBox->MoveResize(367,126,115,22);
    
    gClient->GetColorByName("#ffffff",ucolor);
    TGLabel *NADsLabel = new TGLabel(fFitterFrame2,"Detectors",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
    NADsLabel->SetTextJustify(36);
    NADsLabel->SetMargins(0,0,0,0);
    NADsLabel->SetWrapLength(-1);
    fFitterFrame2->AddFrame(NADsLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    NADsLabel->MoveResize(390,155,128,18);
    
    NADsBox = new TGComboBox(fFitterFrame2,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    NADsBox->SetName("NADsBox");
    NADsBox->AddEntry("6 ADs",6);
    NADsBox->AddEntry("8 ADs",8);
    NADsBox->Connect("Selected(Int_t)", "FitterGui", this, "DoNADs()");
    NADsBox->Select(6);//6 ADs default for now
    fFitterFrame2->AddFrame(NADsBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    NADsBox->MoveResize(424,175,60,22);
    
    gClient->GetColorByName("#ffffff",ucolor);
    TGLabel *fLabel746 = new TGLabel(fFitterFrame2,"# Fits",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
    fLabel746->SetTextJustify(36);
    fLabel746->SetMargins(0,0,0,0);
    fLabel746->SetWrapLength(-1);
    fFitterFrame2->AddFrame(fLabel746, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel746->MoveResize(390,0,128,18);
    
    ToyMCBox = new TGComboBox(fFitterFrame2,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    ToyMCBox->SetName("ToyMCBox");
    ToyMCBox->AddEntry("Fit Data",0);
    ToyMCBox->AddEntry("Fit ToyMC",1);
    ToyMCBox->Connect("Selected(Int_t)", "FitterGui", this, "ChooseToyMC()");
    ToyMCBox->Select(1);//6 ADs default for now
    fFitterFrame2->AddFrame(ToyMCBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    ToyMCBox->MoveResize(400,250,90,22);
    
    MinuitBox = new TGComboBox(fFitterFrame2,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    MinuitBox->SetName("MinuitBox");
    MinuitBox->AddEntry("Manual grid",0);
    MinuitBox->AddEntry("Minuit",1);
    MinuitBox->Connect("Selected(Int_t)", "FitterGui", this, "ChooseMinuit()");
    MinuitBox->Select(0);//Minuit still needs testing
    fFitterFrame2->AddFrame(MinuitBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    MinuitBox->MoveResize(390,280,100,22);
    
    FluctuationBox = new TGComboBox(fFitterFrame2,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    FluctuationBox->SetName("Statistical Fluctuation");
    FluctuationBox->AddEntry("Stat fluc Off",0);
    FluctuationBox->AddEntry("Stat fluc On",1);
    FluctuationBox->Connect("Selected(Int_t)", "FitterGui", this, "DoStatisticalFluctuation()");
    FluctuationBox->Select(0);
    fFitterFrame2->AddFrame(FluctuationBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    FluctuationBox->MoveResize(390,310,100,22);
    
    ToyMCTreeBox = new TGComboBox(fFitterFrame2,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    ToyMCTreeBox->SetName("ToyMC Samples");
    ToyMCTreeBox->AddEntry("Toy MC",0);
    ToyMCTreeBox->AddEntry("Toy MC Tree",1);
    ToyMCTreeBox->Connect("Selected(Int_t)", "FitterGui", this, "DoUseToyMCSamples()");
    ToyMCTreeBox->Select(0);
    fFitterFrame2->AddFrame(ToyMCTreeBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    ToyMCTreeBox->MoveResize(390,330,100,22);
    
    NFitsBox = new TGNumberEntry(fFitterFrame2, (Int_t) NFits,6,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 0,(TGNumberFormat::ELimit) 2,0,1001);
    NFitsBox->SetName("NFits");
    fFitterFrame2->AddFrame(NFitsBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    NFitsBox->MoveResize(424,24,59,22);
    
    NFitsBox->Connect("ValueSet(Long_t)", "FitterGui", this, "DoFits()");
    
    (NFitsBox->GetNumberEntry())->Connect("ReturnPressed()", "FitterGui", this,"DoFits()");
    
    gClient->GetColorByName("#ffffff",ucolor);
    TGLabel *PeriodLabel = new TGLabel(fFitterFrame2,"Period",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
    PeriodLabel->SetTextJustify(36);
    PeriodLabel->SetMargins(0,0,0,0);
    PeriodLabel->SetWrapLength(-1);
    fFitterFrame2->AddFrame(PeriodLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    PeriodLabel->MoveResize(390,200,128,18);
    
    PeriodBox = new TGNumberEntry(fFitterFrame2, (Int_t) Period,6,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 0,(TGNumberFormat::ELimit) 2,0,500);//It can be any number of weeks,500 is about 10 years of data taking, if you need more change this atribute.
    PeriodBox->SetName("Period");
    fFitterFrame2->AddFrame(PeriodBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    PeriodBox->MoveResize(424,220,59,22);
    
    PeriodBox->Connect("ValueSet(Long_t)", "FitterGui", this, "DoPeriod()");
    
    (PeriodBox->GetNumberEntry())->Connect("ReturnPressed()", "FitterGui", this,"DoPeriod()");
    
    gClient->GetColorByName("#ffffff",ucolor);
    TGLabel *fLabel751 = new TGLabel(fFitterFrame2,"Combine Mode",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
    fLabel751->SetTextJustify(36);
    fLabel751->SetMargins(0,0,0,0);
    fLabel751->SetWrapLength(-1);
    fFitterFrame2->AddFrame(fLabel751, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel751->MoveResize(380,56,120,18);
    
    CombineBox = new TGNumberEntry(fFitterFrame2, (Int_t) CombineMode,6,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 0,(TGNumberFormat::ELimit) 2,0,2);
    CombineBox->SetName("CombineBox");
    fFitterFrame2->AddFrame(CombineBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    CombineBox->MoveResize(424,80,59,22);
    
    CombineBox->Connect("ValueSet(Long_t)", "FitterGui", this, "DoCombine()");
    
    (CombineBox->GetNumberEntry())->Connect("ReturnPressed()", "FitterGui", this,"DoCombine()");
    
    gClient->GetColorByName("#333399",ucolor);
    
    // horizontal frame
    TGHorizontalFrame *PlotFrame = new TGHorizontalFrame(fFitterFrame2,245,144,kHorizontalFrame,ucolor);
    PlotFrame->SetName("PlotFrame");
    PlotFrame->SetLayoutBroken(kTRUE);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valButton758;
    valButton758.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#ffffff",valButton758.fForeground);
    gClient->GetColorByName("#e0e0e0",valButton758.fBackground);
    valButton758.fFillStyle = kFillSolid;
    valButton758.fFont = ufont->GetFontHandle();
    valButton758.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valButton758, kTRUE);
    
    TGTextButton *fTextButton758 = new TGTextButton(PlotFrame,"Generate ToyMC Samples",-1,uGC->GetGC());
    fTextButton758->Connect("Clicked()","FitterGui",this,"GenerateToyMC()");
    fTextButton758->SetTextJustify(36);
    fTextButton758->SetMargins(0,0,0,0);
    fTextButton758->SetWrapLength(-1);
    fTextButton758->Resize(232,24);
    
    fTextButton758->ChangeBackground(ucolor);
    PlotFrame->AddFrame(fTextButton758, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTextButton758->MoveResize(8,24,232,24);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valButton759;
    valButton759.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#ffffff",valButton759.fForeground);
    gClient->GetColorByName("#e0e0e0",valButton759.fBackground);
    valButton759.fFillStyle = kFillSolid;
    valButton759.fFont = ufont->GetFontHandle();
    valButton759.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valButton759, kTRUE);
    
    TGTextButton *fTextButton759 = new TGTextButton(PlotFrame,"Plot Error Budget",-1,uGC->GetGC());
    fTextButton759->Connect("Clicked()","FitterGui",this,"PlotErrorBudget()");
    fTextButton759->SetTextJustify(36);
    fTextButton759->SetMargins(0,0,0,0);
    fTextButton759->SetWrapLength(-1);
    fTextButton759->Resize(232,24);
    
    fTextButton759->ChangeBackground(ucolor);
    PlotFrame->AddFrame(fTextButton759, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTextButton759->MoveResize(8,72,232,24);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valButton760;
    valButton760.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#ffffff",valButton760.fForeground);
    gClient->GetColorByName("#e0e0e0",valButton760.fBackground);
    valButton760.fFillStyle = kFillSolid;
    valButton760.fFont = ufont->GetFontHandle();
    valButton760.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valButton760, kTRUE);
    
    TGTextButton *fTextButton760 = new TGTextButton(PlotFrame,"Plot Spectrum Fraction",-1,uGC->GetGC());
    fTextButton760->Connect("Clicked()","FitterGui",this,"PlotSpectrumFraction()");
    fTextButton760->SetTextJustify(36);
    fTextButton760->SetMargins(0,0,0,0);
    fTextButton760->SetWrapLength(-1);
    fTextButton760->Resize(232,24);
    
    fTextButton760->ChangeBackground(ucolor);
    PlotFrame->AddFrame(fTextButton760, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTextButton760->MoveResize(8,96,232,24);
    
    TGTextButton *ResponseMatrixB = new TGTextButton(PlotFrame,"Plot Response Matrix",-1,uGC->GetGC());
    ResponseMatrixB->Connect("Clicked()","FitterGui",this,"PlotResponseMatrix()");
    ResponseMatrixB->SetTextJustify(36);
    ResponseMatrixB->SetMargins(0,0,0,0);
    ResponseMatrixB->SetWrapLength(-1);
    ResponseMatrixB->Resize(232,24);
    
    ResponseMatrixB->ChangeBackground(ucolor);
    PlotFrame->AddFrame(ResponseMatrixB, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    ResponseMatrixB->MoveResize(8,120,232,24);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valButton762;
    valButton762.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#ffffff",valButton762.fForeground);
    gClient->GetColorByName("#e0e0e0",valButton762.fBackground);
    valButton762.fFillStyle = kFillSolid;
    valButton762.fFont = ufont->GetFontHandle();
    valButton762.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valButton762, kTRUE);
    
    TGTextButton *fTextButton762 = new TGTextButton(PlotFrame,"Plot AD True Spectrum",-1,uGC->GetGC());
    fTextButton762->Connect("Clicked()","FitterGui",this,"PlotADTrue()");
    fTextButton762->SetTextJustify(36);
    fTextButton762->SetMargins(0,0,0,0);
    fTextButton762->SetWrapLength(-1);
    fTextButton762->Resize(232,24);
    
    fTextButton762->ChangeBackground(ucolor);
    PlotFrame->AddFrame(fTextButton762, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTextButton762->MoveResize(8,48,232,24);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valButton757;
    valButton757.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#ffffff",valButton757.fForeground);
    gClient->GetColorByName("#e0e0e0",valButton757.fBackground);
    valButton757.fFillStyle = kFillSolid;
    valButton757.fFont = ufont->GetFontHandle();
    valButton757.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valButton757, kTRUE);
    
    TGTextButton *fTextButton757 = new TGTextButton(PlotFrame,"Plot chi square distribution",-1,uGC->GetGC());
    fTextButton757->Connect("Clicked()","FitterGui",this,"PlotChi()");
    fTextButton757->SetTextJustify(36);
    fTextButton757->SetMargins(0,0,0,0);
    fTextButton757->SetWrapLength(-1);
    fTextButton757->Resize(232,24);
    
    fTextButton757->ChangeBackground(ucolor);
    PlotFrame->AddFrame(fTextButton757, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTextButton757->MoveResize(8,0,232,24);
    
    fFitterFrame2->AddFrame(PlotFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    PlotFrame->MoveResize(120,160,245,144);
    
    TGTextButton *RunFitterButton = new TGTextButton(fFitterFrame2,"Run Fitter");
    RunFitterButton->Connect("Clicked()","FitterGui",this,"RunFitter()");
    RunFitterButton->SetTextJustify(36);
    RunFitterButton->SetMargins(0,0,0,0);
    RunFitterButton->SetWrapLength(-1);
    RunFitterButton->Resize(92,24);
    
    fFitterFrame2->AddFrame(RunFitterButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    RunFitterButton->MoveResize(200,310,92,24);
    
    //    ToyMCProgressBar = new TGHProgressBar(fFitterFrame2,TGProgressBar::kFancy,500);
    
    FitterProgressBar = new TGHProgressBar(fFitterFrame2,TGProgressBar::kFancy,500);
    FitterProgressBar->SetName("ToyMCProgressBar");
    FitterProgressBar->ChangeOptions(kSunkenFrame | kDoubleBorder | kOwnBackground);
    FitterProgressBar->SetBarColor("green");
    //    ToyMCProgressBar->SetForegroundColor(0x000000);//To change text color, black so far
    if(Fit2D)
    {
        FitterProgressBar->SetRange(0,NFits*NFits);
    }
    else
    {
        FitterProgressBar->SetRange(0,NFits);
    }
    FitterProgressBar->ShowPosition(kTRUE,kFALSE,"Fit #%.0f ");
    
    FitterProgressBar->MoveResize(154,342,184,24);
    
    fFitterFrame2->AddFrame(FitterProgressBar, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    
    fFitterFrame->AddFrame(fFitterFrame2, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    fFitterFrame2->MoveResize(0,0,490,372);
    
    fRunFitterFrame->AddFrame(fFitterFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    fFitterFrame->MoveResize(0,0,490,372);
    
    fRunFitterFrame->SetMWMHints(kMWMDecorAll, kMWMFuncAll, kMWMInputModeless);
    fRunFitterFrame->MapSubwindows();
    
    fRunFitterFrame->Resize(fRunFitterFrame->GetDefaultSize());
    fRunFitterFrame->MapWindow();
    fRunFitterFrame->Resize(490,372);
}
void FitterGui::RunFitter()
{
    
    gBenchmark->Start("Fitter");
    std::cout << "Running Fitter" << std::endl;
    if(!PrintOnConsole)
    {
        coutstream = cout.rdbuf(0);//change stream of cout
    }
    const bool RangeSin = 0;    //  Fit over a range of 20 Sin22t13
    const bool RangeDelta = 0;  //  Fit over a range of 20 Delta2Mee
    
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
    
    Int_t DataSet=2;
    NominalData* Data = new NominalData(Analysis,DataSet);
    
    Data->SetToyMCSamplesDirectory(ToyMCSampleDirectory);
    Data->SetPredictionDirectory(NominalPredictionsDirectory);
    Data->SetResponseDirectory(ResponseMatrixDirectory);
    Data->SetBkgCovDirectory(BkgCovarianceMatrixDirectory);
    Data->SetSysCovDirectory(SysCovarianceMatrixDirectory);
    
    //  No variations in the ToyMC Prediction, if desired set to 1:
    Data->SetAllRandomSystematics(0);
    Data->SetStatisticalFluctuation(StatisticalFluctuation);
    
    Data->SetCombineMode(CombineMode); //0 is 9x9, 1 is 1x1 and 2 is 2x2
    Data->SetUseToyMCTree(UseToyMCTree);
    //    Data->SetSin22t12(0);
    //    Data->SetSin22t13(0);
    //Parameters of the model
    Data->SetAnalysis(Analysis);//  Gd or H data
    Data->SetBinning(Binning);//  0 for LBNL binning or 1 for Linear binning
    Data->SetNSteps(NFits);// 101 in the final version.
    Data->SetWeeks(Period);
    Data->SetNReactorPeriods(NReactorPeriods);
    Data->SetBCWModel(NL[0]);
    Data->SetLBNLModel(NL[1]);
    Data->SetUnifiedModel(NL[2]);
    
    Data->SetToyMC(ToyMC);
    
    std::string FluxInputS;
    if(Analysis)
    {
        AnalysisString = "Hydrogen";
        FluxInputS="HInputs";
    }
    else
    {
        AnalysisString = "Gadolinium";
        FluxInputS="GdInputs";
    }
    
    NReactorPeriods=Data->GetNReactorPeriods();

    // Generate SuperHistogram:
    NominalData* FluxData = new NominalData(0,2);//Same for nH and nGd, difference is on efficiencies.
    FluxData->LoadMainData(("./Inputs/"+FluxInputS+Form("/Theta13-inputs_%dweek.txt",NReactorPeriods)).c_str());
    Oscillation* FluxOsc= new Oscillation(FluxData);
    FluxOsc->GenerateFluxHisto();
    delete FluxOsc;
    delete FluxData;
    
    //Chose Data Set
    if(Data->GetAnalysis())//   Hydrogen data
    {
        if(DataSet==1) // P12E production data given by Xiang Pan
        {
            
            if(1==Data->GetWeeks())
            {
                std::cout << "\t Loading nH P12E Data" << std::endl;
                Data->LoadMainData("./Inputs/HInputs/P12E_Inclusive.txt");
            }
            else
            {
                std::cout << "\t NO MULTIPLE WEEK P12E DATA IN H ANALYSIS YET " << std::endl;
                Data->LoadMainData(Form("./Inputs/HInputs/P12E_%d.txt",NReactorPeriods));
            }
        }
        else // Simple reactor model used as input data
        {
            std::cout << "\t Loading simple reactor model" << std::endl;
        }
    }
    else//  Gd data
    {
        if(DataSet==2)// P12E production data given by Xiang Pan
        {
            if(1==Data->GetWeeks())
            {
                std::cout << "\t Loading LBNL Gd P12E Data" << std::endl;
                //                Data->LoadMainData("./Inputs/GdInputs/P12E_Inclusive.txt");
            }
            else
            {
                std::cout << "\t NO MULTIPLE WEEK P12E DATA IN LBNL Gd ANALYSIS YET " << std::endl;
                exit(EXIT_FAILURE);

                Data->LoadMainData(Form("./Inputs/GdInputs/Theta13-inputs_%dweek_inclusive.txt",NReactorPeriods));
            }
        }
        else if(DataSet==1)
        {
            if(1==Data->GetWeeks())
            {
                std::cout << "\t Loading Gd Data by Xiang Pan" << std::endl;
                //                Data->LoadMainData("./Inputs/GdInputs/Theta13-inputs_20week_inclusive.txt");
            }
            else
            {
                std::cout << "\t Loading 20 weeks Gd Data by Xiang Pan" << std::endl;
                Data->LoadMainData(Form("./Inputs/GdInputs/Theta13-inputs_%dweek_inclusive.txt",NReactorPeriods));
            }
        }
        else//  Simple reactor model used as input data
        {
            std::cout << "\t Loading simple reactor model" << std::endl;
        }
    }
    
    // Specific combinations to produce the error budget:
    
    if (TurnOnBudget||TurnOffBudget)
    {
        Data->SetTurnOffBudget(TurnOffBudget);
        Data->SetTurnOnBudget(TurnOnBudget);
        
        for (Int_t week = 0; week < Period; week++)
        {
            //Reset files
            std::ofstream ofs;
            
            ofs.open(Form("DMTurnOnErrorBudget_%d.txt",week), std::ofstream::out | std::ofstream::trunc);
            ofs.close();
            
            ofs.open(Form("DMTurnOffErrorBudget_%d.txt",week), std::ofstream::out | std::ofstream::trunc);
            ofs.close();
            
            ofs.open(Form("SinTurnOnErrorBudget_%d.txt",week), std::ofstream::out | std::ofstream::trunc);
            ofs.close();
            
            ofs.open(Form("SinTurnOffErrorBudget_%d.txt",week), std::ofstream::out | std::ofstream::trunc);
            ofs.close();
        }
        
        if(AutomaticBudget)
        {
            for (Int_t i = 0; i<=(18); i++) //17+1, the 18th is to produce the normal fit with all systematics in the TurnOff case or only statistics if TurnOn is used.
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
                Data->SetEfficiencyBudget(14==i);
                Data->SetSystematicBudget(15==i);
                Data->SetBackgroundBudget(16==i);
                Data->SetTotalBudget(17==i);
                
                Fitter* Fit = new Fitter(Data);
                
                Fit->MainFitter(Minuit,Fit2D,FitSin22t13,this);
                
                if(!Minuit)
                {
                    if (Fit2D)
                    {
                        Fit->Save2DFit(i);
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
                }
                delete Fit;
            }
        }
        else
        {
            Data->SetVaryAccidentalBudget(0);
            Data->SetVaryLiHeBudget(0);
            Data->SetVaryFastNeutronsBudget(0);
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
            Data->SetEfficiencyBudget(0);
            Data->SetSin22t12Budget(0);
            Data->SetSystematicBudget(0);
            Data->SetBackgroundBudget(0);
            Data->SetTotalBudget(1);
            
            Fitter* Fit = new Fitter(Data);
            Fit->MainFitter(Minuit,Fit2D,FitSin22t13,this);
            
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
            Double_t Sin22t13[21];
            
            for (Int_t i = 0; i<21; i++)
            {
                Sin22t13[i] = i*.2/(20);
            }
            
            for (Int_t i = 0; i<21; i++)
            {
                Data->SetSin22t13(Sin22t13[i]);
                
                Fitter* Fit = new Fitter(Data);
                Fit->MainFitter(Minuit,Fit2D,FitSin22t13,this);
                Fit->SaveSinRangeChiSquare(i,FitSin22t13);
                delete Fit;
            }
        }
        else if(RangeDelta)
        {
            Double_t Dm2ee[21];
            
            for (Int_t i = 0; i<21; i++)
            {
                Dm2ee[i] = 0.0015+(i*.002/20);
            }
            
            for (Int_t i = 0; i<21; i++)
            {
                Data->SetDm2ee(Dm2ee[i]);
                
                Fitter* Fit = new Fitter(Data);
                Fit->MainFitter(Minuit,Fit2D,FitSin22t13,this);
                Fit->SaveDeltaMRangeChiSquare(i,FitSin22t13);
                delete Fit;
            }
        }
        else
        {
            Fitter* Fit = new Fitter(Data);
            Fit->MainFitter(Minuit,Fit2D,FitSin22t13,this);
            
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
    cout.rdbuf(coutstream);
    gBenchmark->Show("Fitter");
}

FitterGui::~FitterGui()
{
    // Clean up used widgets: frames, buttons, layouthints
    fMain->Cleanup();
    delete fMain;
}

void FitterGui::DoReadInputs()
{
    // main frame
    TGMainFrame *fReadInputs = new TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame);
    fReadInputs->SetName("fReadInputs");
    fReadInputs->SetLayoutBroken(kTRUE);
    fReadInputs->SetWindowName(" Load Inputs Menu ");
    
    // composite frame
    TGCompositeFrame *fLoadInputs = new TGCompositeFrame(fReadInputs,849,587,kVerticalFrame);
    fLoadInputs->SetName("fLoadInputs");
    fLoadInputs->SetLayoutBroken(kTRUE);
    
    TGFont *ufont;         // will reflect user font changes
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    ULong_t ucolor;        // will reflect user color changes
    gClient->GetColorByName("#ffffff",ucolor);
    
    TGGC   *uGC;           // will reflect user GC changes
    // graphics context changes
    GCValues_t vall717;
    vall717.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall717.fForeground);
    gClient->GetColorByName("#e0e0e0",vall717.fBackground);
    vall717.fFillStyle = kFillSolid;
    vall717.fFont = ufont->GetFontHandle();
    vall717.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall717, kTRUE);
    
    gClient->GetColorByName("#ffffff",ucolor);
    
    // composite frame
    fLoadInputs2 = new TGCompositeFrame(fLoadInputs,490,372,kVerticalFrame,ucolor);
    fLoadInputs2->SetName("fLoadInputs2");
    fLoadInputs2->SetLayoutBroken(kTRUE);
    
    fLoadInputs->AddFrame(fLoadInputs2, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    fLoadInputs2->MoveResize(0,0,490,372);
    
    fReadInputs->AddFrame(fLoadInputs, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    fLoadInputs->MoveResize(0,0,490,372);
    
    TGTextButton *fResponseButton= new TGTextButton(fLoadInputs2,"Load Response Matrix",-1,uGC->GetGC());
    fResponseButton->Connect("Clicked()","FitterGui",this,"LoadResponseMatrix()");
    fResponseButton->SetTextJustify(36);
    fResponseButton->SetMargins(0,0,0,0);
    fResponseButton->SetWrapLength(-1);
    fResponseButton->Resize(232,24);
    
    fLoadInputs2->AddFrame(fResponseButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fResponseButton->MoveResize(8,24,232,24);
    
    TGTextButton *fNominalPredictions= new TGTextButton(fLoadInputs2,"Load Nominal Predictions",-1,uGC->GetGC());
    fNominalPredictions->Connect("Clicked()","FitterGui",this,"LoadNominalPredictions()");
    fNominalPredictions->SetTextJustify(36);
    fNominalPredictions->SetMargins(0,0,0,0);
    fNominalPredictions->SetWrapLength(-1);
    fNominalPredictions->Resize(232,24);
    
    fLoadInputs2->AddFrame(fNominalPredictions, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fNominalPredictions->MoveResize(8,56,232,24);
    
    
    TGTextButton *fToyMCSamples= new TGTextButton(fLoadInputs2,"Load Toy MC Samples",-1,uGC->GetGC());
    fToyMCSamples->Connect("Clicked()","FitterGui",this,"LoadToyMCSamples()");
    fToyMCSamples->SetTextJustify(36);
    fToyMCSamples->SetMargins(0,0,0,0);
    fToyMCSamples->SetWrapLength(-1);
    fToyMCSamples->Resize(232,24);
    
    fLoadInputs2->AddFrame(fToyMCSamples, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fToyMCSamples->MoveResize(8,88,232,24);
    
    fReadInputs->SetMWMHints(kMWMDecorAll, kMWMFuncAll, kMWMInputModeless);
    fReadInputs->MapSubwindows();
    
    fReadInputs->Resize(fReadInputs->GetDefaultSize());
    fReadInputs->MapWindow();
    fReadInputs->Resize(490,372);
}

void FitterGui::LoadToyMCSamples()
{
    TGFileInfo file_info_;
    const char *filetypes[] = {"Root Files", "*.root", 0, 0};
    file_info_.fFileTypes = filetypes;
    file_info_.fIniDir = StrDup(".");
    
    new TGFileDialog(gClient->GetDefaultRoot(), gClient->GetDefaultRoot(), kFDOpen, &file_info_);
    
    ToyMCSampleDirectory = file_info_.fFilename;
    
    std::cout << "Loading Toy MC samples from: " <<  ToyMCSampleDirectory << std::endl;
    
}

void FitterGui::LoadNominalPredictions()
{
    TGFileInfo file_info_;
    const char *filetypes[] = {"Root Files", "*.root", 0, 0};
    file_info_.fFileTypes = filetypes;
    file_info_.fIniDir = StrDup(".");
    
    new TGFileDialog(gClient->GetDefaultRoot(), gClient->GetDefaultRoot(), kFDOpen, &file_info_);
    
    NominalPredictionsDirectory = file_info_.fFilename;
    
    std::cout << "Loading predictions from: " <<  NominalPredictionsDirectory << std::endl;
    
}

void FitterGui::LoadResponseMatrix()
{
    TGFileInfo file_info_;
    const char *filetypes[] = {"Root Files", "*.root", 0, 0};
    file_info_.fFileTypes = filetypes;
    file_info_.fIniDir = StrDup(".");
    
    new TGFileDialog(gClient->GetDefaultRoot(), gClient->GetDefaultRoot(), kFDOpen, &file_info_);
    ResponseMatrixDirectory = file_info_.fFilename;
    
    std::cout << "Loading response matrix from: " <<  ResponseMatrixDirectory << std::endl;
    
}

void FitterGui::DoToyMC()
{
    // main frame
    TGMainFrame *fRunToyMCFrame = new TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame);
    fRunToyMCFrame->SetName("fRunToyMCFrame");
    fRunToyMCFrame->SetLayoutBroken(kTRUE);
    fRunToyMCFrame->SetWindowName(" Toy MC Menu ");
    
    // composite frame
    TGCompositeFrame *fToyMCFrame = new TGCompositeFrame(fRunToyMCFrame,849,587,kVerticalFrame);
    fToyMCFrame->SetName("fToyMCFrame");
    fToyMCFrame->SetLayoutBroken(kTRUE);
    
    ULong_t ucolor;        // will reflect user color changes
    gClient->GetColorByName("#ffffff",ucolor);
    
    // composite frame
    fToyMCFrame2 = new TGCompositeFrame(fToyMCFrame,490,372,kVerticalFrame,ucolor);
    fToyMCFrame2->SetName("fToyMCFrame2");
    fToyMCFrame2->SetLayoutBroken(kTRUE);
    
    TGFont *ufont;         // will reflect user font changes
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    TGGC   *uGC;           // will reflect user GC changes
    // graphics context changes
    GCValues_t vall717;
    vall717.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall717.fForeground);
    gClient->GetColorByName("#e0e0e0",vall717.fBackground);
    vall717.fFillStyle = kFillSolid;
    vall717.fFont = ufont->GetFontHandle();
    vall717.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall717, kTRUE);
    
    gClient->GetColorByName("#ffffff",ucolor);
    TGLabel *fLabel717 = new TGLabel(fToyMCFrame2,"nH/nGd",uGC->GetGC(),ufont->GetFontStruct(),kChildFrame,ucolor);
    fLabel717->SetTextJustify(36);
    fLabel717->SetMargins(0,0,0,0);
    fLabel717->SetWrapLength(-1);
    fToyMCFrame2->AddFrame(fLabel717, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel717->MoveResize(8,0,121,18);
    
    gClient->GetColorByName("#ffffff",ucolor);
    
    TGLabel *ChoseCovarianceMatrixLabel = new TGLabel(fToyMCFrame2,"Choose Covariance Matrix",uGC->GetGC(),ufont->GetFontStruct(),kChildFrame,ucolor);
    ChoseCovarianceMatrixLabel->SetTextJustify(36);
    ChoseCovarianceMatrixLabel->SetMargins(0,0,0,0);
    ChoseCovarianceMatrixLabel->SetWrapLength(-1);
    fToyMCFrame2->AddFrame(ChoseCovarianceMatrixLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    ChoseCovarianceMatrixLabel->MoveResize(182,0,200,18);
    
    // combo box
    
    CovarianceMatrixBox = new TGComboBox(fToyMCFrame2,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    CovarianceMatrixBox->SetName("CovarianceMatrixBox");
    
    CovarianceMatrixBox->AddEntry("Run all Covariance Matrices",0);
    CovarianceMatrixBox->AddEntry("Vary Accidental Matrix",1);
    CovarianceMatrixBox->AddEntry("Vary LiHe Matrix",2);
    CovarianceMatrixBox->AddEntry("Vary FN Matrix",3);
    CovarianceMatrixBox->AddEntry("Vary AmC Matrix",4);
    CovarianceMatrixBox->AddEntry("Distort LiHe Matrix",5);
    CovarianceMatrixBox->AddEntry("Distort FN Matrix",6);
    CovarianceMatrixBox->AddEntry("Distort AmC Matrix",7);
    CovarianceMatrixBox->AddEntry("Reactor Spectrum Matrix",8);
    CovarianceMatrixBox->AddEntry("Reactor Power Matrix",9);
    CovarianceMatrixBox->AddEntry("IAV Matrix",10);
    CovarianceMatrixBox->AddEntry("NL Matrix",11);
    CovarianceMatrixBox->AddEntry("Resolution Matrix",12);
    CovarianceMatrixBox->AddEntry("Efficiency Matrix",13);
    CovarianceMatrixBox->AddEntry("Sin^{2}(2#theta_{12}) Matrix",14);
    CovarianceMatrixBox->AddEntry("Relative Energy Scale Matrix",15);
    
    CovarianceMatrixBox->Connect("Selected(Int_t)", "FitterGui", this, "ChooseCovarianceMatrix()");
    CovarianceMatrixBox->Select(0);//Run all covariance matrices as default
    fToyMCFrame2->AddFrame(CovarianceMatrixBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    CovarianceMatrixBox->MoveResize(182,23,200,22);
    
    // combo box
    
    NLBox = new TGButtonGroup(fToyMCFrame2,"Non Linearity Model",kVerticalFrame,TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),0xffffff);
    NLBox->SetName("NLBox");
    fR[0] = new TGRadioButton(NLBox,"BCW Model",0);
    fR[1] = new TGRadioButton(NLBox,"LBNL Model",1);
    fR[2] = new TGRadioButton(NLBox,"Unified Model",2);
    
    fR[0]->ChangeBackground(ucolor);
    fR[1]->ChangeBackground(ucolor);
    fR[2]->ChangeBackground(ucolor);
    
    fR[2]->SetState(kButtonDown);
    
    for (Int_t i = 0; i < 3; ++i)
    {
        fR[i]->Connect("Pressed()", "FitterGui",  this, "DoNLModel()");
    }
    NLBox->Show();
    
    //    TGComboBox *NLBox = new TGComboBox(fToyMCFrame2,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    //    NLBox->SetName("NLBox");
    //    NLBox->AddEntry("BCW Model",1);
    //    NLBox->AddEntry("LBNL Model",2);
    //    NLBox->AddEntry("Unified Model",3);
    //    NLBox->Resize(102,22);
    //    NLBox->Select(3);
    fToyMCFrame2->AddFrame(NLBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    NLBox->MoveResize(8,56,150,90);
    
    gClient->GetColorByName("#ffffff",ucolor);
    
    // combo box
    HGdBox = new TGComboBox(fToyMCFrame2,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    HGdBox->SetName("HGdBox");
    HGdBox->AddEntry("Gadollinium Analysis",0);
    HGdBox->AddEntry("Hydrogen Analysis",1);
    HGdBox->Resize(152,22);
    HGdBox->Connect("Selected(Int_t)", "FitterGui", this, "DoAnalysis()");
    HGdBox->Select(0);//Gd Default right now.
    fToyMCFrame2->AddFrame(HGdBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    HGdBox->MoveResize(8,24,152,22);
    
    gClient->GetColorByName("#ffffff",ucolor);
    TGLabel *binninglabel = new TGLabel(fToyMCFrame2,"Binning",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
    binninglabel->SetTextJustify(36);
    binninglabel->SetMargins(0,0,0,0);
    binninglabel->SetWrapLength(-1);
    fToyMCFrame2->AddFrame(binninglabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    binninglabel->MoveResize(390,108,128,18);
    
    BinningBox = new TGComboBox(fToyMCFrame2,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    BinningBox->SetName("BinningBox");
    BinningBox->AddEntry("LBNL Binning",0);
    BinningBox->AddEntry("Linear Binning",1);
    BinningBox->Connect("Selected(Int_t)", "FitterGui", this, "DoBinning()");
    BinningBox->Select(0);//LBNL Binning default for now
    fToyMCFrame2->AddFrame(BinningBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    BinningBox->MoveResize(367,126,118,22);
    
    gClient->GetColorByName("#ffffff",ucolor);
    TGLabel *NADsLabel = new TGLabel(fToyMCFrame2,"Detectors",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
    NADsLabel->SetTextJustify(36);
    NADsLabel->SetMargins(0,0,0,0);
    NADsLabel->SetWrapLength(-1);
    fToyMCFrame2->AddFrame(NADsLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    NADsLabel->MoveResize(390,155,128,18);
    
    NADsBox = new TGComboBox(fToyMCFrame2,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    NADsBox->SetName("NADsBox");
    NADsBox->AddEntry("6 ADs",6);
    NADsBox->AddEntry("8 ADs",8);
    NADsBox->Connect("Selected(Int_t)", "FitterGui", this, "DoNADs()");
    NADsBox->Select(6);//6 ADs default for now
    fToyMCFrame2->AddFrame(NADsBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    NADsBox->MoveResize(424,175,60,22);
    
    gClient->GetColorByName("#ffffff",ucolor);
    TGLabel *fLabel746 = new TGLabel(fToyMCFrame2,"# Samples",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
    fLabel746->SetTextJustify(36);
    fLabel746->SetMargins(0,0,0,0);
    fLabel746->SetWrapLength(-1);
    fToyMCFrame2->AddFrame(fLabel746, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel746->MoveResize(390,0,128,18);
    
    NSamplesBox = new TGNumberEntry(fToyMCFrame2, (Int_t) NSamples,6,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 0,(TGNumberFormat::ELimit) 2,0,MaximumSamples);
    NSamplesBox->SetName("NSamples");
    fToyMCFrame2->AddFrame(NSamplesBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    NSamplesBox->MoveResize(424,24,59,22);
    
    NSamplesBox->Connect("ValueSet(Long_t)", "FitterGui", this, "DoSamples()");
    
    (NSamplesBox->GetNumberEntry())->Connect("ReturnPressed()", "FitterGui", this,"DoSamples()");
    
    gClient->GetColorByName("#ffffff",ucolor);
    TGLabel *PeriodLabel = new TGLabel(fToyMCFrame2,"Period",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
    PeriodLabel->SetTextJustify(36);
    PeriodLabel->SetMargins(0,0,0,0);
    PeriodLabel->SetWrapLength(-1);
    fToyMCFrame2->AddFrame(PeriodLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    PeriodLabel->MoveResize(390,200,128,18);
    
    PeriodBox = new TGNumberEntry(fToyMCFrame2, (Int_t) Period,6,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 0,(TGNumberFormat::ELimit) 2,0,500);//It can be any number of weeks,500 is about 10 years of data taking, if you need more change this atribute.
    PeriodBox->SetName("Period");
    fToyMCFrame2->AddFrame(PeriodBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    PeriodBox->MoveResize(424,220,59,22);
    
    PeriodBox->Connect("ValueSet(Long_t)", "FitterGui", this, "DoPeriod()");
    
    (PeriodBox->GetNumberEntry())->Connect("ReturnPressed()", "FitterGui", this,"DoPeriod()");
    
    gClient->GetColorByName("#ffffff",ucolor);
    TGLabel *fLabel751 = new TGLabel(fToyMCFrame2,"Combine Mode",TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),kChildFrame,ucolor);
    fLabel751->SetTextJustify(36);
    fLabel751->SetMargins(0,0,0,0);
    fLabel751->SetWrapLength(-1);
    fToyMCFrame2->AddFrame(fLabel751, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel751->MoveResize(380,56,120,18);
    
    CombineBox = new TGNumberEntry(fToyMCFrame2, (Int_t) CombineMode,6,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 0,(TGNumberFormat::ELimit) 2,0,2);
    CombineBox->SetName("CombineBox");
    fToyMCFrame2->AddFrame(CombineBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    CombineBox->MoveResize(424,80,59,22);
    
    CombineBox->Connect("ValueSet(Long_t)", "FitterGui", this, "DoCombine()");
    
    (CombineBox->GetNumberEntry())->Connect("ReturnPressed()", "FitterGui", this,"DoCombine()");
    
    gClient->GetColorByName("#333399",ucolor);
    
    // horizontal frame
    TGHorizontalFrame *PlotFrame = new TGHorizontalFrame(fToyMCFrame2,245,144,kHorizontalFrame,ucolor);
    PlotFrame->SetName("PlotFrame");
    PlotFrame->SetLayoutBroken(kTRUE);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valButton758;
    valButton758.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#ffffff",valButton758.fForeground);
    gClient->GetColorByName("#e0e0e0",valButton758.fBackground);
    valButton758.fFillStyle = kFillSolid;
    valButton758.fFont = ufont->GetFontHandle();
    valButton758.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valButton758, kTRUE);
    
    TGTextButton *fTextButton758 = new TGTextButton(PlotFrame,"Plot Backgrounds",-1,uGC->GetGC());
    fTextButton758->Connect("Clicked()","FitterGui",this,"PlotBkgd()");
    fTextButton758->SetTextJustify(36);
    fTextButton758->SetMargins(0,0,0,0);
    fTextButton758->SetWrapLength(-1);
    fTextButton758->Resize(232,24);
    
    fTextButton758->ChangeBackground(ucolor);
    PlotFrame->AddFrame(fTextButton758, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTextButton758->MoveResize(8,24,232,24);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valButton759;
    valButton759.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#ffffff",valButton759.fForeground);
    gClient->GetColorByName("#e0e0e0",valButton759.fBackground);
    valButton759.fFillStyle = kFillSolid;
    valButton759.fFont = ufont->GetFontHandle();
    valButton759.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valButton759, kTRUE);
    
    TGTextButton *fTextButton759 = new TGTextButton(PlotFrame,"Plot AD Visible Spectrum",-1,uGC->GetGC());
    fTextButton759->Connect("Clicked()","FitterGui",this,"PlotADVis()");
    fTextButton759->SetTextJustify(36);
    fTextButton759->SetMargins(0,0,0,0);
    fTextButton759->SetWrapLength(-1);
    fTextButton759->Resize(232,24);
    
    fTextButton759->ChangeBackground(ucolor);
    PlotFrame->AddFrame(fTextButton759, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTextButton759->MoveResize(8,72,232,24);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valButton760;
    valButton760.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#ffffff",valButton760.fForeground);
    gClient->GetColorByName("#e0e0e0",valButton760.fBackground);
    valButton760.fFillStyle = kFillSolid;
    valButton760.fFont = ufont->GetFontHandle();
    valButton760.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valButton760, kTRUE);
    
    TGTextButton *fTextButton760 = new TGTextButton(PlotFrame,"Plot Spectrum Fraction",-1,uGC->GetGC());
    fTextButton760->Connect("Clicked()","FitterGui",this,"PlotSpectrumFraction()");
    fTextButton760->SetTextJustify(36);
    fTextButton760->SetMargins(0,0,0,0);
    fTextButton760->SetWrapLength(-1);
    fTextButton760->Resize(232,24);
    
    fTextButton760->ChangeBackground(ucolor);
    PlotFrame->AddFrame(fTextButton760, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTextButton760->MoveResize(8,96,232,24);
    
    TGTextButton *ResponseMatrixB = new TGTextButton(PlotFrame,"Plot Response Matrix",-1,uGC->GetGC());
    ResponseMatrixB->Connect("Clicked()","FitterGui",this,"PlotResponseMatrix()");
    ResponseMatrixB->SetTextJustify(36);
    ResponseMatrixB->SetMargins(0,0,0,0);
    ResponseMatrixB->SetWrapLength(-1);
    ResponseMatrixB->Resize(232,24);
    
    ResponseMatrixB->ChangeBackground(ucolor);
    PlotFrame->AddFrame(ResponseMatrixB, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    ResponseMatrixB->MoveResize(8,120,232,24);
    
    
    TGTextButton *VariationsB = new TGTextButton(PlotFrame,"Plot Spectrum Variations",-1,uGC->GetGC());
    VariationsB->Connect("Clicked()","FitterGui",this,"PlotVariations()");
    VariationsB->SetTextJustify(36);
    VariationsB->SetMargins(0,0,0,0);
    VariationsB->SetWrapLength(-1);
    VariationsB->Resize(232,24);
    
    VariationsB->ChangeBackground(ucolor);
    PlotFrame->AddFrame(VariationsB, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    VariationsB->MoveResize(8,144,232,24);
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valButton762;
    valButton762.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#ffffff",valButton762.fForeground);
    gClient->GetColorByName("#e0e0e0",valButton762.fBackground);
    valButton762.fFillStyle = kFillSolid;
    valButton762.fFont = ufont->GetFontHandle();
    valButton762.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valButton762, kTRUE);
    
    TGTextButton *fTextButton762 = new TGTextButton(PlotFrame,"Plot AD True Spectrum",-1,uGC->GetGC());
    fTextButton762->Connect("Clicked()","FitterGui",this,"PlotADTrue()");
    fTextButton762->SetTextJustify(36);
    fTextButton762->SetMargins(0,0,0,0);
    fTextButton762->SetWrapLength(-1);
    fTextButton762->Resize(232,24);
    
    fTextButton762->ChangeBackground(ucolor);
    PlotFrame->AddFrame(fTextButton762, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTextButton762->MoveResize(8,48,232,24);
    
    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
    
    // graphics context changes
    GCValues_t valButton757;
    valButton757.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#ffffff",valButton757.fForeground);
    gClient->GetColorByName("#e0e0e0",valButton757.fBackground);
    valButton757.fFillStyle = kFillSolid;
    valButton757.fFont = ufont->GetFontHandle();
    valButton757.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valButton757, kTRUE);
    
    TGTextButton *fTextButton757 = new TGTextButton(PlotFrame,"Plot Covariance Matrices",-1,uGC->GetGC());
    fTextButton757->Connect("Clicked()","FitterGui",this,"PlotCov()");
    fTextButton757->SetTextJustify(36);
    fTextButton757->SetMargins(0,0,0,0);
    fTextButton757->SetWrapLength(-1);
    fTextButton757->Resize(232,24);
    
    fTextButton757->ChangeBackground(ucolor);
    PlotFrame->AddFrame(fTextButton757, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fTextButton757->MoveResize(8,0,232,24);
    
    fToyMCFrame2->AddFrame(PlotFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    PlotFrame->MoveResize(120,160,245,168);
    
    TGTextButton *RunToyMCButton = new TGTextButton(fToyMCFrame2,"Run Toy MC");
    RunToyMCButton->Connect("Clicked()","FitterGui",this,"RunToyMC()");
    RunToyMCButton->SetTextJustify(36);
    RunToyMCButton->SetMargins(0,0,0,0);
    RunToyMCButton->SetWrapLength(-1);
    RunToyMCButton->Resize(92,24);
    
    fToyMCFrame2->AddFrame(RunToyMCButton, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    RunToyMCButton->MoveResize(200,80,112,44);
    
    //    ToyMCProgressBar = new TGHProgressBar(fToyMCFrame2,TGProgressBar::kFancy,500);
    
    ToyMCProgressBar = new TGHProgressBar(fToyMCFrame2,TGProgressBar::kFancy,500);
    ToyMCProgressBar->SetName("ToyMCProgressBar");
    ToyMCProgressBar->ChangeOptions(kSunkenFrame | kDoubleBorder | kOwnBackground);
    ToyMCProgressBar->SetBarColor("green");
    //    ToyMCProgressBar->SetForegroundColor(0x000000);//To change text color, black so far
    ToyMCProgressBar->SetRange(0,NSamples);
    ToyMCProgressBar->ShowPosition(kTRUE,kFALSE,"Sample #%.0f ");
    
    ToyMCProgressBar->MoveResize(154,342,184,24);
    
    fToyMCFrame2->AddFrame(ToyMCProgressBar, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    
    fToyMCFrame->AddFrame(fToyMCFrame2, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    fToyMCFrame2->MoveResize(0,0,490,372);
    
    fRunToyMCFrame->AddFrame(fToyMCFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    fToyMCFrame->MoveResize(0,0,490,372);
    
    //    // terminate ROOT session when window is closed
    //    fRunToyMCFrame->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    //    fRunToyMCFrame->DontCallClose();
    
    fRunToyMCFrame->SetMWMHints(kMWMDecorAll, kMWMFuncAll, kMWMInputModeless);
    fRunToyMCFrame->MapSubwindows();
    
    fRunToyMCFrame->Resize(fRunToyMCFrame->GetDefaultSize());
    fRunToyMCFrame->MapWindow();
    fRunToyMCFrame->Resize(490,372);
}

void FitterGui::RunToyMC()
{
    gBenchmark->Start("ToyMC");
    std::cout << "Running Covariance Matrix ToyMC" << std::endl;
    if(!PrintOnConsole)
    {
        coutstream = cout.rdbuf(0);//change stream of cout
    }
    const Int_t CovarianceMatrices=15;
    
    //To draw using a better palette:
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    
    Int_t DataSet=2;//0 is Simulation, 2 is P12E
    
    NominalData* Data = new NominalData(Analysis,DataSet);
    
    Data->SetToyMCSamplesDirectory(ToyMCSampleDirectory);
    Data->SetPredictionDirectory(NominalPredictionsDirectory);
    Data->SetResponseDirectory(ResponseMatrixDirectory);
    Data->SetBkgCovDirectory(BkgCovarianceMatrixDirectory);
    Data->SetSysCovDirectory(SysCovarianceMatrixDirectory);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // CrossCalc takes a lot of time and it needs to be run only once to produce the root file. Uncomment only if binning is different from standards and you want to produce a new cross section for this new binning.
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(Analysis)//nGd cross section is already tabulated, this is to calculate nH
    {
        CrossSection* Cross = new CrossSection(Data);
        Cross->CrossCalc();
        delete Cross;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Parameters of the model
    
    Data->SetNSamples(NSamples);// 500 in the final version.
    
    Data->SetToyMC(1);//  1 for Toy MC, 0 for data. To produce covariance matrices we use ToyMC.
    
    Data->SetCombineMode(CombineMode); //0 is 9x9, 1 is 1x1 and 2 is 2x2
    Data->SetUseToyMCTree(UseToyMCTree);
    Data->SetBinning(Binning);//  0 for LBNL binning or 1 for Linear binning
    Data->SetWeeks(Period);
    Data->SetNReactorPeriods(NReactorPeriods);
    std::cout <<  "********************************************************************************************************" << std::endl;
    std::cout <<  "**********************************            MAIN             *****************************************" << std::endl;
    std::cout <<  "********************************************************************************************************" << std::endl;
    
    //Chose Data Set
    if(Data->GetAnalysis())//   Hydrogen data
    {
        switch (DataSet)
        {
            case 0://  Simple reactor model used as input data
                std::cout << "\t Loading simple reactor model" << std::endl;
                break;
            case 1://   P12E
                if(1==Data->GetWeeks())
                {
                    std::cout << "\t Loading nH P12E Data" << std::endl;
                    Data->LoadMainData("./Inputs/HInputs/P12E_Inclusive.txt");
                }
                else
                {
                    std::cout << "\t \t \t NO MULTIPLE WEEK P12E DATA IN H ANALYSIS YET " << std::endl;
                    exit(EXIT_FAILURE);
                    Data->LoadMainData(Form("./Inputs/HInputs/P12E_%d.txt",NReactorPeriods));
                }
                break;
            case 2:// LBNL
                std::cout << "\t \t \t NO LBNL H ANALYSIS, OPTION NOT VALID " << std::endl;
                exit(EXIT_FAILURE);
                break;
            default:
                break;
        }
    }
    else//  Gd data
    {
        switch (DataSet)
        {
            case 0://  Simple reactor model used as input data
                std::cout << "\t Loading simple reactor model" << std::endl;
                break;
            case 1://   P12E
                if(1==Data->GetWeeks())
                {
                    std::cout << "\t Loading Gd P12E Data" << std::endl;
                    //                    Data->LoadMainData("./Inputs/GdInputs/P12E_Inclusive.txt");
                }
                else
                {
                    std::cout << "\t \t \t NO MULTIPLE WEEK P12E DATA IN Gd ANALYSIS YET " << std::endl;
                    exit(EXIT_FAILURE);
                    Data->LoadMainData(Form("./Inputs/GdInputs/P12E_%d.txt",NReactorPeriods));
                }
                break;
            case 2:// LBNL
                if(1==Data->GetWeeks())
                {
                    std::cout << "\t Loading LBNL Gd Data" << std::endl;
                    //                    Data->LoadMainData("./Inputs/GdInputs/Theta13-inputs_20week_inclusive.txt");
                }
                else
                {
                    std::cout << "\t Loading weekly LBNL Gd Data" << std::endl;
                    Data->LoadMainData(Form("./Inputs/GdInputs/Theta13-inputs_%dweek.txt",NReactorPeriods));
                }
                break;
            default:
                break;
        }
    }
    
    //To check non theta13 oscillation behaviour of the fitter. To avoid any oscillation at all then theta12 has to be set to 0 too.
    //    Data->SetSin22t12(0);
    //    Data->SetSin22t13(0);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    FitBackgrounds2* Bkg = new FitBackgrounds2(Data);
    
    if(Analysis)
    {
        Bkg->ReadHBackgrounds();
    }
    else
    {
        Bkg->ReadGdBackgrounds();
    }
    
    delete Bkg;
    
    CovarianceMatrix3* Cov;
    
    if (Automatic==1)//1 for producing all matrices, others to choose specific combinations manually
    {//Automatic Script:
        std::cout <<  "********************************************************************************************************" << std::endl;
        std::cout << "\t Automatic Script" << std::endl;
        
        Data->SetBCWModel(NL[0]);
        Data->SetLBNLModel(NL[1]);
        Data->SetUnifiedModel(NL[2]);
        //Order is important (Inititialize NL model before randomize it)
        
        for (Int_t i = 0; i<(CovarianceMatrices); i++)
        {
            Data->SetVaryAccidentalMatrix(0==i);
            Data->SetVaryLiHeMatrix(1==i);
            Data->SetVaryFastNeutronsMatrix(2==i);
            Data->SetVaryAmCMatrix(3==i);
            Data->SetDistortLiHeMatrix(4==i);
            Data->SetDistortFastNeutronsMatrix(5==i);
            Data->SetDistortAmCMatrix(6==i);
            Data->SetIsotopeMatrix(7==i);
            Data->SetReactorPowerMatrix(8==i);
            Data->SetRelativeEnergyScaleMatrix(9==i);
            //            Data->SetAbsoluteEnergyScaleMatrix(10==i);
            //            Data->SetRelativeEnergyOffsetMatrix(11==i);
            //            Data->SetAbsoluteEnergyOffsetMatrix(12==i);
            Data->SetIAVMatrix(10==i);
            Data->SetNLMatrix(11==i);
            Data->SetResolutionMatrix(12==i);
            Data->SetSin22t12Matrix(13==i);
            Data->SetEfficiencyMatrix(14==i);
            
            Cov = new CovarianceMatrix3(Data);
            Cov->CovarianceMatrixMain(this);
            delete Cov;
        }
    }
    else if(Automatic==0)
    {
        std::cout <<  "********************************************************************************************************" << std::endl;
        std::cout << "\t Manual Script" << std::endl;
        
        Data->SetBCWModel(NL[0]);
        Data->SetLBNLModel(NL[1]);
        Data->SetUnifiedModel(NL[2]);
        
        Data->SetVaryAccidentalMatrix(VaryAccidentalMatrix);
        Data->SetVaryLiHeMatrix(VaryLiHeMatrix);
        Data->SetVaryFastNeutronsMatrix(VaryFastNeutronsMatrix);
        Data->SetVaryAmCMatrix(VaryAmCMatrix);
        Data->SetDistortLiHeMatrix(DistortLiHeMatrix);
        Data->SetDistortFastNeutronsMatrix(DistortFastNeutronsMatrix);
        Data->SetDistortAmCMatrix(DistortAmCMatrix);
        Data->SetIsotopeMatrix(IsotopeMatrix);
        Data->SetReactorPowerMatrix(ReactorPowerMatrix);
        Data->SetRelativeEnergyScaleMatrix(EnergyScaleMatrix);
        //        Data->SetAbsoluteEnergyScaleMatrix(AbsoluteScaleMatrix);
        //        Data->SetRelativeEnergyOffsetMatrix(EnergyOffsetMatrix);
        //        Data->SetAbsoluteEnergyOffsetMatrix(AbsoluteOffsetMatrix);
        Data->SetIAVMatrix(IAVMatrix);
        Data->SetNLMatrix(NLMatrix);
        Data->SetResolutionMatrix(ResolutionMatrix);
        Data->SetSin22t12Matrix(Sin22t12Matrix);
        Data->SetEfficiencyMatrix(EfficiencyMatrix);
        
        Cov = new CovarianceMatrix3(Data);
        
        Cov->CovarianceMatrixMain(this);
        
        delete Cov;
    }
    
    delete Data;
    cout.rdbuf(coutstream);
    gBenchmark->Show("ToyMC");
    
}

void FitterGui::DoSamples()
{
    NSamples = (Int_t)NSamplesBox->GetNumberEntry()->GetIntNumber();
    std::cout << "NSAMPLES : " << NSamples << std::endl;
    ToyMCProgressBar->SetRange(0,NSamples);
}

void FitterGui::DoFits()
{
    NFits = (Int_t)NFitsBox->GetNumberEntry()->GetIntNumber();
    
    if(Fit2D)
    {
        std::cout << "FIT GRID : " << NFits << " x " << NFits << std::endl;
        FitterProgressBar->SetRange(0,NFits*NFits);
    }
    else
    {
        FitterProgressBar->SetRange(0,NFits);
        std::cout << "FIT RANGE : " << NFits  << std::endl;
    }
}

void FitterGui::DoCombine()
{
    CombineMode = (Int_t)CombineBox->GetNumberEntry()->GetIntNumber();
    std::cout << "COMBINE MODE : " << CombineMode << std::endl;
}

void FitterGui::DoBinning()
{
    Binning = BinningBox->GetSelected();
    std::cout << "BINNING MODE : " << CombineMode << " ( 0 is LBNL binning, 1 is linear binning )" << std::endl;
    
}

void FitterGui::DoStatisticalFluctuation()
{
    StatisticalFluctuation = FluctuationBox->GetSelected();
    std::cout << "STATISTICAL FLUCTUATION : " << StatisticalFluctuation << std::endl;
}

void FitterGui::DoUseToyMCSamples()
{
    UseToyMCTree = ToyMCTreeBox->GetSelected();
    std::cout << "Use Toy MC Tree : " << UseToyMCTree << std::endl;
}

void FitterGui :: DoFitterMode()
{
    TGButton *btn = (TGButton *) gTQSender;
    Int_t id = btn->WidgetId();
    
    for (int i = 0; i < 2; i++)
    {
        FitterMode[i] = 1;
        
        if (fF[i]->WidgetId() != id)
        {
            FitterMode[i] = 0;
            fF[i]->SetState(kButtonUp);
        }
    }
    
    
    if(FitterMode[0]==1)
    {
        Fit2D = 0;
        Fitter1DBox->Show();
        std::cout << "1D FIT" << std::endl;
        FitterProgressBar->SetRange(0,NFits);
        std::cout << "# Fits set :" << NFits << std::endl;
        
    }
    else
    {
        Fit2D = 1;
        Fitter1DBox->Hide();
        std::cout << "2D FIT" << std::endl;
        FitterProgressBar->SetRange(0,NFits*NFits);
        std::cout << "# Fits set :" << NFits*NFits << std::endl;
    }
}

void FitterGui :: Do1DFitterMode()
{
    TGButton *btn = (TGButton *) gTQSender;
    Int_t id = btn->WidgetId();
    
    for (int i = 0; i < 2; i++)
    {
        Fitter1DMode[i] = 1;
        
        if (f1DF[i]->WidgetId() != id)
        {
            Fitter1DMode[i] = 0;
            f1DF[i]->SetState(kButtonUp);
        }
    }
    
    
    if(Fitter1DMode[0]==1)
    {
        FitSin22t13 = 1;
        std::cout << "Sin22theta13 FIT" << std::endl;
    }
    else
    {
        FitSin22t13 = 0;
        std::cout << "DeltaM FIT" << std::endl;
    }
}

void FitterGui :: DoNLModel()
{
    TGButton *btn = (TGButton *) gTQSender;
    Int_t id = btn->WidgetId();
    
    for (int i = 0; i < 3; i++)
    {
        NL[i] = 1;
        
        if (fR[i]->WidgetId() != id)
        {
            NL[i] = 0;
            fR[i]->SetState(kButtonUp);
        }
    }
    
    std::cout << "NL MODEL : BCW " << NL[0] << " - LBNL MODEL " << NL[1] << " - UNIFIED MODEL " << NL[2] << std::endl;
    
}

void FitterGui :: DoNear()
{
    TGButton *btn = (TGButton *) gTQSender;
    Int_t id = btn->WidgetId();
    
    for (int i = 0; i < 3; i++)
    {
        if (fI[i]->WidgetId() != id)
        {
            fI[i]->SetState(kButtonUp);
        }
        else
        {
            NearTrueIndex = i;
        }
    }
    
    PlotFarADTrue();
    
    std::cout << "Plot Far Prediction from Near AD1:" << NearTrueIndex << std::endl;
    
}

void FitterGui :: DoNominal()
{
    
    TGButton *btn = (TGButton *) gTQSender;
    Int_t id = btn->WidgetId();
    
    if (id == 0)
    {
        NominalString = "Nominal";
        fN[id]->SetState(kButtonUp);
    }
    else
    {
        NominalString = "";
        fN[id]->SetState(kButtonUp);
    }
    
    if(flagCombine == 1)//Update the current active window
    {
        PlotCombine();
    }
    else
    {
        PlotAllADVis();
    }
    
    std::cout << "Plot Nominal/Random Spectrum : " << NominalString  << std::endl;
    
}

void FitterGui::DoChangeBin()
{
    PlotBin = hslider->GetPosition();
    
    if(flagNear == 1)//Update the current active window
    {
        PlotNear();
    }
    else
    {
        PlotFar();
    }
    
    std::cout << "BIN: " << PlotBin << std::endl;
}

void FitterGui::DoAnalysis()
{
    Analysis = HGdBox->GetSelected();
    if(Analysis)
    {
        AnalysisString = "Hydrogen";
    }
    else
    {
        AnalysisString = "Gadolinium";
    }
    
    std::cout << "Analysis: " << AnalysisString << std::endl;
}
void FitterGui::Update(Int_t samples)
{
    ToyMCProgressBar->SetPosition(samples+1);
    
    gSystem->ProcessEvents();
}

void FitterGui::UpdateFitter(Int_t samples)
{
    FitterProgressBar->SetPosition(samples+1);
    
    gSystem->ProcessEvents();
}

void FitterGui::PlotVariations()
{
    TGMainFrame *fPlotVariationsFrame = new TGMainFrame(gClient->GetRoot(),1000,800,kMainFrame | kVerticalFrame);
    fPlotVariationsFrame->SetName("fPlotVariationsFrame");
    fPlotVariationsFrame->SetWindowName(" Plot Spectrum Variations ");
    
    fVariationsCanvas = new TRootEmbeddedCanvas("VariationsCanvas",fPlotVariationsFrame,1000,800);
    VCanvas = fVariationsCanvas->GetCanvas();
    VarString = "NL"; //default
    TGHorizontalFrame *hframe = new TGHorizontalFrame(fPlotVariationsFrame,1000,40);
    
    VariationsPlotBox = new TGComboBox(hframe,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    VariationsPlotBox->SetName("VariationsMatrixBox");
    
    VariationsPlotBox->AddEntry("Vary Accidental Variations",1);
    VariationsPlotBox->AddEntry("Vary LiHe Variations",2);
    VariationsPlotBox->AddEntry("Vary FN Variations",3);
    VariationsPlotBox->AddEntry("Vary AmC Variations",4);
    VariationsPlotBox->AddEntry("Distort LiHe Variations",5);
    VariationsPlotBox->AddEntry("Distort FN Variations",6);
    VariationsPlotBox->AddEntry("Distort AmC Variations",7);
    VariationsPlotBox->AddEntry("Reactor Spectrum Variations",8);
    VariationsPlotBox->AddEntry("Reactor Power Variations",9);
    VariationsPlotBox->AddEntry("IAV Variations",10);
    VariationsPlotBox->AddEntry("NL Variations",11);
    VariationsPlotBox->AddEntry("Resolution Variations",12);
    VariationsPlotBox->AddEntry("Efficiency Variations",13);
    VariationsPlotBox->AddEntry("Sin^{2}(2#theta_{12}) Variations",14);
    VariationsPlotBox->AddEntry("Relative Energy Scale Variations",15);
    
    VariationsPlotBox->Connect("Selected(Int_t)", "FitterGui", this, "ChooseVariations()");
    VariationsPlotBox->Select(11);
    VariationsPlotBox->MoveResize(0,0,200,22);
    hframe->AddFrame(VariationsPlotBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    
    VariationsBox = new TGButtonGroup(hframe,"Plot Ratio/Spectrum",kVerticalFrame,TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),0xffffff);
    VariationsBox->SetName("VariationsBox");
    fV[0] = new TGRadioButton(VariationsBox,"Plot Ratio",0);
    fV[1] = new TGRadioButton(VariationsBox,"Plot Spectrum",1);
    
    fV[0]->ChangeBackground(0xffffff);
    fV[1]->ChangeBackground(0xffffff);
    
    fV[0]->SetState(kButtonDown);
    
    for (Int_t i = 0; i < 2; ++i)
    {
        fV[i]->Connect("Pressed()", "FitterGui",  this, "ChoosePlotVariations()");
    }
    
    VariationsBox->Show();
    
    hframe->AddFrame(VariationsBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    VariationsBox->MoveResize(400,800,400,15);
    
    hframe->MoveResize(0,800,800,40);
    
    fPlotVariationsFrame->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    fPlotVariationsFrame->AddFrame(fVariationsCanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY, 10,10,10,1));
    
    
    // Set a name to the main frame
    fPlotVariationsFrame->SetWindowName(" Plot Spectrum Variations ");
    
    // Map all subwindows of main frame
    fPlotVariationsFrame->MapSubwindows();
    
    // Initialize the layout algorithm
    fPlotVariationsFrame->Resize(fPlotVariationsFrame->GetDefaultSize());
    
    // Map main frame
    fPlotVariationsFrame->MapWindow();
    
}
void FitterGui::PlotBkgd()
{
    TGMainFrame *fPlotBkgdFrame = new TGMainFrame(gClient->GetRoot(),1000,800,kMainFrame | kVerticalFrame);
    fPlotBkgdFrame->SetName("fPlotBkgdFrame");
    fPlotBkgdFrame->SetWindowName(" Plot Background ");
    fPlotBkgdFrame->SetLayoutBroken(kTRUE);
    
    fBackgroundcanvas = new TRootEmbeddedCanvas("Backgroundcanvas",fPlotBkgdFrame,1000,800);
    fCanvas = fBackgroundcanvas->GetCanvas();
    fCanvas->Divide(3,2);
    fPlotBkgdFrame->AddFrame(fBackgroundcanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY, 10,10,10,1));
    
    TGHorizontalFrame *hframe = new TGHorizontalFrame(fPlotBkgdFrame,1000,40);
    
    TGTextButton *PlotAccidentalB = new TGTextButton(fPlotBkgdFrame,"&Accidentals");
    PlotAccidentalB->Connect("Clicked()","FitterGui",this,"PlotAccidental()");
    hframe->AddFrame(PlotAccidentalB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    TGTextButton *PlotLiHeB = new TGTextButton(fPlotBkgdFrame,"&Li/He");
    PlotLiHeB->Connect("Clicked()","FitterGui",this,"PlotLiHe()");
    hframe->AddFrame(PlotLiHeB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    TGTextButton *PlotFNB = new TGTextButton(fPlotBkgdFrame,"&Fast Neutrons");
    PlotFNB->Connect("Clicked()","FitterGui",this,"PlotFN()");
    hframe->AddFrame(PlotFNB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    TGTextButton *PlotAmCB = new TGTextButton(fPlotBkgdFrame,"&AmC");
    PlotAmCB->Connect("Clicked()","FitterGui",this,"PlotAmC()");
    hframe->AddFrame(PlotAmCB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    hframe->MoveResize(0,800,800,40);
    
    fPlotBkgdFrame->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    
    // Set a name to the main frame
    fPlotBkgdFrame->SetWindowName(" Background Menu ");
    
    // Map all subwindows of main frame
    fPlotBkgdFrame->MapSubwindows();
    
    // Initialize the layout algorithm
    fPlotBkgdFrame->Resize(fPlotBkgdFrame->GetDefaultSize());
    
    // Map main frame
    fPlotBkgdFrame->MapWindow();
    
}

void FitterGui::PlotAccidental()
{
    TH1D* AccHistograms[MaxDetectors];
    
    if(Analysis)//Hydrogen
    {
        TFile* HydrogenBackgroundsF = new TFile("./BackgroundSpectrum/HBackground/Backgrounds.root");
        
        for(Int_t i = 0; i < NADs; i++)
        {
            AccHistograms[i] = (TH1D*)HydrogenBackgroundsF->Get(Form("Accidentals_AD%d",i));
        }
        
        delete HydrogenBackgroundsF;
    }
    else//Gd
    {
        TFile* GdBackgroundsF = new TFile("./BackgroundSpectrum/GDBackground/Backgrounds.root");
        
        for(Int_t i = 0; i < NADs; i++)
        {
            AccHistograms[i] = (TH1D*)GdBackgroundsF->Get(Form("Accidentals_AD%d",i));
        }
        
        delete GdBackgroundsF;
    }
    // Draws function graphics in randomly choosen interval
    
    for(Int_t i = 0; i < NADs; i++)
    {
        fCanvas->cd(i+1);
        
        AccHistograms[i]->Draw();
    }
    
    fCanvas->Update();
}

void FitterGui::PlotFN()
{
    TH1D* FNHistograms[MaxDetectors];
    
    if(Analysis)//Hydrogen
    {
        TFile* HydrogenBackgroundsF = new TFile("./BackgroundSpectrum/HBackground/Backgrounds.root");
        
        for(Int_t i = 0; i < NADs; i++)
        {
            FNHistograms[i] = (TH1D*)HydrogenBackgroundsF->Get(Form("FN_AD%d",i));
        }
        
        delete HydrogenBackgroundsF;
    }
    else//Gd
    {
        TFile* GdBackgroundsF = new TFile("./BackgroundSpectrum/GDBackground/Backgrounds.root");
        
        for(Int_t i = 0; i < NADs; i++)
        {
            FNHistograms[i] = (TH1D*)GdBackgroundsF->Get(Form("FN_AD%d",i));
        }
        
        delete GdBackgroundsF;
    }
    // Draws function graphics in randomly choosen interval
    
    for(Int_t i = 0; i < NADs; i++)
    {
        fCanvas->cd(i+1);
        
        FNHistograms[i]->Draw();
    }
    
    fCanvas->Update();
}

void FitterGui::PlotLiHe()
{
    TH1D* LiHeHistograms[MaxDetectors];
    
    if(Analysis)//Hydrogen
    {
        TFile* HydrogenBackgroundsF = new TFile("./BackgroundSpectrum/HBackground/Backgrounds.root");
        
        for(Int_t i = 0; i < NADs; i++)
        {
            LiHeHistograms[i] = (TH1D*)HydrogenBackgroundsF->Get(Form("LiHe_AD%d",i));
        }
        
        delete HydrogenBackgroundsF;
    }
    else//Gd
    {
        TFile* GdBackgroundsF = new TFile("./BackgroundSpectrum/GDBackground/Backgrounds.root");
        
        for(Int_t i = 0; i < NADs; i++)
        {
            LiHeHistograms[i] = (TH1D*)GdBackgroundsF->Get(Form("LiHe_AD%d",i));
        }
        
        delete GdBackgroundsF;
    }
    // Draws function graphics in randomly choosen interval
    
    for(Int_t i = 0; i < NADs; i++)
    {
        fCanvas->cd(i+1);
        
        LiHeHistograms[i]->Draw();
    }
    
    fCanvas->Update();
}

void FitterGui::PlotAmC()
{
    TH1D* AmCHistograms[MaxDetectors];
    
    if(Analysis)//Hydrogen
    {
        TFile* HydrogenBackgroundsF = new TFile("./BackgroundSpectrum/HBackground/Backgrounds.root");
        
        for(Int_t i = 0; i < NADs; i++)
        {
            AmCHistograms[i] = (TH1D*)HydrogenBackgroundsF->Get(Form("AmC_AD%d",i));
        }
        
        delete HydrogenBackgroundsF;
    }
    else//Gd
    {
        TFile* GdBackgroundsF = new TFile("./BackgroundSpectrum/GDBackground/Backgrounds.root");
        
        for(Int_t i = 0; i < NADs; i++)
        {
            AmCHistograms[i] = (TH1D*)GdBackgroundsF->Get(Form("AmC_AD%d",i));
        }
        
        delete GdBackgroundsF;
    }
    // Draws function graphics in randomly choosen interval
    
    for(Int_t i = 0; i < NADs; i++)
    {
        fCanvas->cd(i+1);
        
        AmCHistograms[i]->Draw();
    }
    
    fCanvas->Update();
}

void FitterGui::PlotErrorBudget()
{
    TGMainFrame *fPlotErrorBudgetFrame = new TGMainFrame(gClient->GetRoot(),800,800,kMainFrame | kVerticalFrame);
    fPlotErrorBudgetFrame->SetName("fPlotErrorBudget");
    fPlotErrorBudgetFrame->SetWindowName(" Plot Error Budget ");
    //    fPlotADVisFrame->SetLayoutBroken(kTRUE);
    
    fBudgetCanvas = new TRootEmbeddedCanvas("BudgetCanvas",fPlotErrorBudgetFrame,700,700);
    BudgetCanvas = fBudgetCanvas->GetCanvas();
    
    TGHorizontalFrame *hframe = new TGHorizontalFrame(fPlotErrorBudgetFrame,1000,20);
    
    TGTextButton *PlotOnB = new TGTextButton(fPlotErrorBudgetFrame,"&Plot Turn On");
    PlotOnB->Connect("Clicked()","FitterGui",this,"PlotTurnOnBudget()");
    hframe->AddFrame(PlotOnB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    TGTextButton *PlotOffB = new TGTextButton(fPlotErrorBudgetFrame,"&Plot Turn Off");
    PlotOffB->Connect("Clicked()","FitterGui",this,"PlotTurnOffBudget()");
    hframe->AddFrame(PlotOffB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    TGTextButton *CalculateOnB = new TGTextButton(fPlotErrorBudgetFrame,"&Run Turn On");
    CalculateOnB->Connect("Clicked()","FitterGui",this,"RunTurnOnBudget()");
    hframe->AddFrame(CalculateOnB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    TGTextButton *CalculateOffB = new TGTextButton(fPlotErrorBudgetFrame,"&Run Turn Off");
    CalculateOffB->Connect("Clicked()","FitterGui",this,"RunTurnOffBudget()");
    hframe->AddFrame(CalculateOffB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    hframe->MoveResize(0,800,800,20);
    
    fPlotErrorBudgetFrame->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    
    fPlotErrorBudgetFrame->AddFrame(fBudgetCanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY, 10,10,10,1));
    
    // Map all subwindows of main frame
    fPlotErrorBudgetFrame->MapSubwindows();
    
    // Initialize the layout algorithm
    fPlotErrorBudgetFrame->Resize(fPlotErrorBudgetFrame->GetDefaultSize());
    
    // Map main frame
    fPlotErrorBudgetFrame->MapWindow();
}

void FitterGui:: PlotTurnOnBudget()
{
    char *dirname;
    if(Analysis)
    {
        if(FitSin22t13)
        {
            dirname= Form("./ChiSquare/Hydrogen/Combine%d/TurnOnBudget/S2/",CombineMode);
        }
        else
        {
            dirname = Form("./ChiSquare/Hydrogen/Combine%d/TurnOnBudget/DM/",CombineMode);
        }
    }
    else
    {
        if(FitSin22t13)
        {
            dirname = Form("./ChiSquare/Gadolinium/Combine%d/TurnOnBudget/S2/",CombineMode);
        }
        else
        {
            dirname = Form("./ChiSquare/Gadolinium/Combine%d/TurnOnBudget/DM/",CombineMode);
        }
    }
    
    char *ext= Form(".root");
    
    Int_t DataSet=2;
    
    NominalData* Data = new NominalData(Analysis,DataSet);
    
    Data->SetToyMCSamplesDirectory(ToyMCSampleDirectory);
    Data->SetPredictionDirectory(NominalPredictionsDirectory);
    Data->SetResponseDirectory(ResponseMatrixDirectory);
    Data->SetBkgCovDirectory(BkgCovarianceMatrixDirectory);
    Data->SetSysCovDirectory(SysCovarianceMatrixDirectory);
    
    Double_t sin22t13 = Data->GetSin22t13();
    Double_t deltaM = Data->GetDm2ee();
    delete Data;
    
    BudgetCanvas->Clear();
    BudgetCanvas->Resize();
    
    TH1D* Histo[50];
    TH1D* InterpolateHisto[50];
    Double_t minX=0;
    Double_t maxX=10;
    
    Int_t num = 0;
    TString fname[50];
    
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    
    if(files)
    {
        TSystemFile *file;
        TIter next(files);
        
        while ((file=(TSystemFile*)next()))
        {
            fname[num] = file->GetName();
            
            if (!file->IsDirectory() && fname[num].EndsWith(ext))
            {
                std::cout << fname[num] << std::endl;
                
                TFile* File = new TFile(dirname+fname[num]);
                
                Histo[num] = (TH1D*)File->Get(Form("Sin%f_Distribution_DeltaM%f_Period%d",sin22t13,deltaM,Period-1));
                
                delete File;
                
                Double_t Minimum = Histo[num]->GetMinimum();
                
                if(FitSin22t13)
                {
                    InterpolateHisto[num] = new TH1D(Form("#chi^{2} Sin%f distribution",sin22t13),Form("#chi^{2} Sin%f distribution",sin22t13),Histo[num]->GetXaxis()->GetNbins()*InterpolationFactor,0,0.2);
                }
                else
                {
                    InterpolateHisto[num] = new TH1D(Form("#chi^{2} DeltaM%f distribution",deltaM),Form("#chi^{2} DeltaM%f distribution",deltaM),Histo[num]->GetXaxis()->GetNbins()*InterpolationFactor,0.0015,0.0035);
                }
                
                for (Int_t i = 0; i<=((Histo[num]->GetXaxis()->GetNbins()-1)*InterpolationFactor); i++)
                {
                    Double_t width = (Histo[num]->GetXaxis()->GetXmax()-Histo[num]->GetXaxis()->GetXmin())/((Histo[num]->GetXaxis()->GetNbins()-1)*InterpolationFactor);
                    
                    Double_t xvalue =Histo[num]->GetXaxis()->GetXmin()+i*width;
                    
                    InterpolateHisto[num]->SetBinContent(i+1,Histo[num]->Interpolate(xvalue)-Minimum);
                    
                    if((TMath::Floor(InterpolateHisto[num]->GetBinContent(i+1))==0)&&minX<xvalue)//Look for minimum X where Δχ2 == 1
                        
                    {
                        minX = xvalue;
                    }
                    else if((TMath::Floor(InterpolateHisto[num]->GetBinContent(i+1))==0)&&maxX>xvalue)//Look for maximum X where Δχ2 == 1
                        
                    {
                        maxX = xvalue;
                    }
                }
                
                if(num<9)
                {
                    Histo[num]->SetLineColor(num+1);
                    InterpolateHisto[num]->SetLineColor(num+1);
                }
                else
                {
                    Histo[num]->SetLineColor(2*num+1);
                    InterpolateHisto[num]->SetLineColor(2*num+1);
                }
                
                num++;
            }
        }
        
        TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
        legend->SetTextFont(72);
        legend->SetTextSize(0.01);
        
        InterpolateHisto[0]->SetMaximum(2.5);
        InterpolateHisto[0]->SetStats(0);
        InterpolateHisto[0]->GetXaxis()->SetTitleOffset(1.3);
        InterpolateHisto[0]->GetXaxis()->SetTitleSize(0.03);
        InterpolateHisto[0]->GetXaxis()->SetLabelSize(0.025);
        
        if(FitSin22t13)
        {
            InterpolateHisto[0]->GetXaxis()->SetTitle("sin^{2}(2#theta_{13})");
        }
        else
        {
            InterpolateHisto[0]->GetXaxis()->SetTitle("#Delta^{2}m_{23}");
        }
        InterpolateHisto[0]->GetYaxis()->SetTitle("#Delta#chi^{2}");
        InterpolateHisto[0]->SetTitle("");
        InterpolateHisto[0]->Draw();
        
        TBox* B = new TBox(minX,0,maxX,1);
        B->SetFillStyle(0);
        B->SetLineStyle(3);
        
        for(Int_t i = 0; i < num; i++)
        {
            legend->AddEntry(Histo[i],fname[i],"lpe");
            
            InterpolateHisto[i]->Draw("same");
        }
        legend->Draw("same");
        B->Draw("samel");
        
        BudgetCanvas->Update();
        
        if(FitSin22t13)
        {
            BudgetCanvas->Print("./Images/ErrorBudget/TurnOnS2.eps");
        }
        else
        {
            BudgetCanvas->Print("./Images/ErrorBudget/TurnOnDM.eps");
        }
    }
}

void FitterGui:: PlotTurnOffBudget()
{
    char *dirname;
    
    if(Analysis)
    {
        if(FitSin22t13)
        {
            dirname= Form("./ChiSquare/Hydrogen/Combine%d/TurnOffBudget/S2/",CombineMode);
        }
        else
        {
            dirname = Form("./ChiSquare/Hydrogen/Combine%d/TurnOffBudget/DM/",CombineMode);
        }
    }
    else
    {
        if(FitSin22t13)
        {
            dirname = Form("./ChiSquare/Gadolinium/Combine%d/TurnOffBudget/S2/",CombineMode);
        }
        else
        {
            dirname = Form("./ChiSquare/Gadolinium/Combine%d/TurnOffBudget/DM/",CombineMode);
        }
    }
    
    std::cout << dirname << std::endl;
    char *ext= (Form(".root"));
    
    BudgetCanvas->Clear();
    BudgetCanvas->Resize();
    
    Int_t DataSet=2;
    
    NominalData* Data = new NominalData(Analysis,DataSet);
    
    Data->SetToyMCSamplesDirectory(ToyMCSampleDirectory);
    Data->SetPredictionDirectory(NominalPredictionsDirectory);
    Data->SetResponseDirectory(ResponseMatrixDirectory);
    Data->SetBkgCovDirectory(BkgCovarianceMatrixDirectory);
    Data->SetSysCovDirectory(SysCovarianceMatrixDirectory);
    
    Double_t sin22t13 = Data->GetSin22t13();
    Double_t deltaM = Data->GetDm2ee();
    delete Data;
    
    TH1D* Histo[50];
    TH1D* InterpolateHisto[50];
    
    Double_t minX=0;
    Double_t maxX=10;
    Int_t num = 0;
    TString fname[30];
    
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    
    if(files)
    {
        TSystemFile *file;
        TIter next(files);
        
        while ((file=(TSystemFile*)next()))
        {
            fname[num] = file->GetName();
            
            if (!file->IsDirectory() && fname[num].EndsWith(ext))
            {
                std::cout << fname[num] << std::endl;
                
                TFile* FileF = TFile::Open(dirname+fname[num]);
                
                Histo[num] = (TH1D*)FileF->Get(Form("Sin%f_Distribution_DeltaM%f_Period%d",sin22t13,deltaM,Period-1));
                
                Double_t Minimum = Histo[num]->GetMinimum();
                
                FileF->Close();
                
                if(FitSin22t13)
                {
                    InterpolateHisto[num] = new TH1D(Form("#chi^{2} Sin%f distribution",sin22t13),Form("#chi^{2} Sin%f distribution",sin22t13),Histo[num]->GetXaxis()->GetNbins()*InterpolationFactor,0,0.2);
                }
                else
                {
                    InterpolateHisto[num] = new TH1D(Form("#chi^{2} DeltaM%f distribution",deltaM),Form("#chi^{2} DeltaM%f distribution",deltaM),Histo[num]->GetXaxis()->GetNbins()*InterpolationFactor,0.0015,0.0035);
                }
                
                for (Int_t i = 0; i<=((Histo[num]->GetXaxis()->GetNbins()-1)*InterpolationFactor); i++)
                {
                    Double_t width = (Histo[num]->GetXaxis()->GetXmax()-Histo[num]->GetXaxis()->GetXmin())/((Histo[num]->GetXaxis()->GetNbins()-1)*InterpolationFactor);
                    
                    Double_t xvalue = Histo[num]->GetXaxis()->GetXmin()+i*width;
                    
                    InterpolateHisto[num]->SetBinContent(i+1,Histo[num]->Interpolate(xvalue)-Minimum);
                    
                    if((TMath::Floor(InterpolateHisto[num]->GetBinContent(i+1))==0)&&minX<xvalue)//Look for minimum X where Δχ2 == 1
                    {
                        minX = xvalue;
                        
                        std::cout << minX << std::endl;
                        
                    }
                    else if((TMath::Floor(InterpolateHisto[num]->GetBinContent(i+1))==0)&&maxX>xvalue)//Look for maximum X where Δχ2 == 1
                        
                    {
                        maxX = xvalue;
                    }
                }
                
                if(num<9)
                {
                    Histo[num]->SetLineColor(num+1);
                    InterpolateHisto[num]->SetLineColor(num+1);
                }
                else
                {
                    Histo[num]->SetLineColor(2*num+1);
                    InterpolateHisto[num]->SetLineColor(2*num+1);
                }
                num++;
                
            }
        }
        
        TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
        legend->SetTextFont(72);
        legend->SetTextSize(0.01);
        
        InterpolateHisto[0]->SetStats(0);
        InterpolateHisto[0]->SetTitle("");
        InterpolateHisto[0]->GetXaxis()->SetLabelSize(0.025);
        InterpolateHisto[0]->GetXaxis()->SetTitleOffset(1.3);
        InterpolateHisto[0]->GetXaxis()->SetTitleSize(0.03);
        InterpolateHisto[0]->GetYaxis()->SetTitle("#Delta#chi^{2}");
        
        if(FitSin22t13)
        {
            InterpolateHisto[0]->GetXaxis()->SetTitle("sin^{2}(2#theta_{13})");
        }
        else
        {
            InterpolateHisto[0]->GetXaxis()->SetTitle("#Delta^{2}m_{23}");
        }
        
        InterpolateHisto[0]->SetMaximum(2.5);
        InterpolateHisto[0]->Draw();
        TBox* B = new TBox(minX,0,maxX,1);
        B->SetFillStyle(0);
        B->SetLineStyle(3);
        
        for(Int_t i = 0; i < num; i++)
        {
            legend->AddEntry(Histo[i],fname[i],"lpe");
            InterpolateHisto[i]->Draw("same");
        }
        legend->Draw("same");
        B->Draw("samel");
        BudgetCanvas->Update();
        
        if(FitSin22t13)
        {
            BudgetCanvas->Print("./Images/ErrorBudget/TurnOffS2.eps");
        }
        else
        {
            BudgetCanvas->Print("./Images/ErrorBudget/TurnOffDM.eps");
        }
    }
}


void FitterGui:: RunTurnOnBudget()
{
    AutomaticBudget = 1; //  To loop the fitter over all possible systematics
    TurnOnBudget=1; // Generates Turn On Error Budget
    TurnOffBudget=0; // Generates Turn Off Error Budget
    
    RunFitter();
    
    TurnOnBudget=0;
    TurnOffBudget=0;
}

void FitterGui :: RunTurnOffBudget()
{
    AutomaticBudget = 1; //  To loop the fitter over all possible systematics
    TurnOnBudget=0; // Generates Turn On Error Budget
    TurnOffBudget=1; // Generates Turn Off Error Budget
    
    RunFitter();
    
    TurnOnBudget=0;
    TurnOffBudget=0;
}

void FitterGui::PlotADVis()
{
    TGMainFrame *fPlotADVisFrame = new TGMainFrame(gClient->GetRoot(),800,800,kMainFrame | kVerticalFrame);
    fPlotADVisFrame->SetName("fPlotADVis");
    fPlotADVisFrame->SetWindowName(" Plot AD Vis ");
    //    fPlotADVisFrame->SetLayoutBroken(kTRUE);
    
    fADVisCanvas = new TRootEmbeddedCanvas("ADViscanvas",fPlotADVisFrame,700,700);
    ADVisCanvas = fADVisCanvas->GetCanvas();
    
    TGHorizontalFrame *hframe = new TGHorizontalFrame(fPlotADVisFrame,1000,20);
    
    TGTextButton *PlotADVisB = new TGTextButton(fPlotADVisFrame,"&AD Visible");
    PlotADVisB->Connect("Clicked()","FitterGui",this,"PlotAllADVis()");
    hframe->AddFrame(PlotADVisB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    TGTextButton *PlotCombineB = new TGTextButton(fPlotADVisFrame,"&Combined");
    PlotCombineB->Connect("Clicked()","FitterGui",this,"PlotCombine()");
    hframe->AddFrame(PlotCombineB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    NominalBox = new TGButtonGroup(hframe,"Nominal/Random",kHorizontalFrame,TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),0xffffff);
    NominalBox->SetName("NominalBox");
    fN[0] = new TGRadioButton(NominalBox,"Nominal Spectra",0);
    fN[1] = new TGRadioButton(NominalBox,"Random Spectra",1);
    
    fN[0]->SetState(kButtonDown);
    
    for (Int_t i = 0; i < 2; ++i)
    {
        fN[i]->Connect("Pressed()", "FitterGui",  this, "DoNominal()");
    }
    NominalBox->Show();
    hframe->AddFrame(NominalBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    NominalBox->MoveResize(400,800,400,15);
    
    hframe->MoveResize(0,800,800,20);
    
    fPlotADVisFrame->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    
    fPlotADVisFrame->AddFrame(fADVisCanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY, 10,10,10,1));
    
    // Map all subwindows of main frame
    fPlotADVisFrame->MapSubwindows();
    
    // Initialize the layout algorithm
    fPlotADVisFrame->Resize(fPlotADVisFrame->GetDefaultSize());
    
    // Map main frame
    fPlotADVisFrame->MapWindow();
}

void FitterGui::PlotAllADVis()
{
    flagCombine = 0;
    ADVisCanvas->Clear();
    ADVisCanvas->SetCanvasSize(700,700);
    ADVisCanvas->Divide(3,3);
    
    const Int_t NearADs = 3;
    const Int_t FarADs = 3;
    
    TH1D* ADVisH[NearADs][FarADs];
    
    
    TFile* ADVisF = new TFile(("./RootOutputs/"+AnalysisString+Form("/Spectra/Combine%d/",CombineMode)+NominalString+"PredictedSpectrum.root").c_str());
    
    for(Int_t i = 0; i < NearADs; i++)
    {
        for(Int_t j = 0; j < FarADs; j++)
        {
            ADVisH[i][j] = (TH1D*)ADVisF->Get(Form("Vis Prediction AD%d from AD%d",j+1,i+1));
        }
    }
    
    delete ADVisF;
    
    // Draws function graphics in randomly choosen interval
    
    for(Int_t i = 0; i < NearADs; i++)
    {
        for(Int_t j = 0; j < FarADs; j++)
        {
            ADVisCanvas->cd(j+i*FarADs+1);
            
            ADVisH[i][j]->Draw();
        }
    }
    
    ADVisCanvas->Update();
    
    for(Int_t i = 0; i < NearADs; i++)
    {
        for(Int_t j = 0; j < FarADs; j++)
        {
            delete  ADVisH[i][j];
        }
    }
}

void FitterGui::PlotCombine()
{
    flagCombine = 1;
    ADVisCanvas->Clear();
    ADVisCanvas->Update();
    
    const Int_t NearADs = 3;
    const Int_t FarADs = 3;
    
    Int_t MaxCombineNearADs= 3;
    Int_t MaxCombineFarADs= 3;
    
    TH1D* ADCombine[NearADs][FarADs];
    
    if(CombineMode==1)
    {
        MaxCombineNearADs = 1;
        MaxCombineFarADs = 1;
    }
    else if(CombineMode==2)
    {
        ADVisCanvas->SetCanvasSize(800, 400);
        ADVisCanvas->Divide(2,1);
        
        MaxCombineNearADs = 2;
        MaxCombineFarADs = 1;
    }
    else
    {
        ADVisCanvas->Divide(3,3);
    }
    
    TFile* ADVisF = new TFile(("./RootOutputs/"+AnalysisString+Form("/Spectra/Combine%d/",CombineMode)+NominalString+"PredictedSpectrum.root").c_str());
    
    for(Int_t i = 0; i < MaxCombineNearADs; i++)
    {
        for(Int_t j = 0; j < MaxCombineFarADs; j++)
        {
            ADCombine[i][j] = (TH1D*)ADVisF->Get(Form("Combined Prediction AD%d from AD%d",j+1,i+1));
        }
    }
    
    delete ADVisF;
    
    // Draws function graphics in randomly choosen interval
    
    for(Int_t i = 0; i < MaxCombineNearADs; i++)
    {
        for(Int_t j = 0; j < MaxCombineFarADs; j++)
        {
            ADVisCanvas->cd(j+i*MaxCombineFarADs+1);
            
            ADCombine[i][j]->Draw();
        }
    }
    
    ADVisCanvas->Update();
    
    for(Int_t i = 0; i < MaxCombineNearADs; i++)
    {
        for(Int_t j = 0; j < MaxCombineFarADs; j++)
        {
            delete  ADCombine[i][j];
        }
    }
}

void FitterGui::PlotADTrue()
{
    //I can sum all true energy histograms and all reactors to produce this histogram.
    
    TGMainFrame *fPlotADTrueFrame = new TGMainFrame(gClient->GetRoot(),800,800,kMainFrame | kVerticalFrame);
    fPlotADTrueFrame->SetName("fPlotADTrue");
    fPlotADTrueFrame->SetWindowName(" Plot AD True ");
    //    fPlotADTrueFrame->SetLayoutBroken(kTRUE);
    
    fADTrueCanvas = new TRootEmbeddedCanvas("ADTruecanvas",fPlotADTrueFrame,700,700);
    ADTrueCanvas = fADTrueCanvas->GetCanvas();
    
    TGHorizontalFrame *hframe = new TGHorizontalFrame(fPlotADTrueFrame,1000,20);
    
    TGTextButton *PlotADFarTrueB = new TGTextButton(fPlotADTrueFrame,"&True AD Far");
    PlotADFarTrueB->Connect("Clicked()","FitterGui",this,"PlotFarADTrue()");
    hframe->AddFrame(PlotADFarTrueB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    TGTextButton *PlotADNearTrueB = new TGTextButton(fPlotADTrueFrame,"&True AD Near");
    PlotADNearTrueB->Connect("Clicked()","FitterGui",this,"PlotNearADTrue()");
    hframe->AddFrame(PlotADNearTrueB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    IndexBox = new TGButtonGroup(hframe,"Show far prediction from Near AD",kHorizontalFrame,TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),0xffffff);
    IndexBox->SetName("IndexBox");
    
    fI[0] = new TGRadioButton(IndexBox,"Near AD1 ",0);
    fI[1] = new TGRadioButton(IndexBox,"Near AD2 ",1);
    fI[2] = new TGRadioButton(IndexBox,"Near AD3 ",2);
    
    fI[0]->SetState(kButtonDown);
    
    for (Int_t i = 0; i < 3; ++i)
    {
        fI[i]->Connect("Pressed()", "FitterGui",  this, "DoNear()");
    }
    
    IndexBox->Show();
    
    hframe->AddFrame(IndexBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    IndexBox->MoveResize(400,800,400,15);
    
    hframe->MoveResize(0,800,800,20);
    
    fPlotADTrueFrame->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    
    fPlotADTrueFrame->AddFrame(fADTrueCanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY, 10,10,10,1));
    
    // Map all subwindows of main frame
    fPlotADTrueFrame->MapSubwindows();
    
    // Initialize the layout algorithm
    fPlotADTrueFrame->Resize(fPlotADTrueFrame->GetDefaultSize());
    
    // Map main frame
    fPlotADTrueFrame->MapWindow();
}

void FitterGui :: PlotSpectrumFraction()
{
    TGMainFrame *fPlotFractionFrame = new TGMainFrame(gClient->GetRoot(),800,800,kMainFrame | kVerticalFrame);
    fPlotFractionFrame->SetName("fPlotADVis");
    fPlotFractionFrame->SetWindowName(" Plot Spectrum Reactor Fraction ");
    //    fPlotADVisFrame->SetLayoutBroken(kTRUE);
    
    fFractionCanvas = new TRootEmbeddedCanvas("ADViscanvas",fPlotFractionFrame,800,800);
    FractionCanvas = fFractionCanvas->GetCanvas();
    
    TGHorizontalFrame *hframe = new TGHorizontalFrame(fPlotFractionFrame,1000,40);
    
    TGTextButton *PlotNearB = new TGTextButton(fPlotFractionFrame,"&Near Fraction");
    PlotNearB->Connect("Clicked()","FitterGui",this,"PlotNear()");
    hframe->AddFrame(PlotNearB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    TGTextButton *PlotFarB = new TGTextButton(fPlotFractionFrame,"&Far Fraction");
    PlotFarB->Connect("Clicked()","FitterGui",this,"PlotFar()");
    hframe->AddFrame(PlotFarB, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    hslider = new TGHSlider(fPlotFractionFrame,390);
    
    hslider->Connect("PositionChanged(Int_t)", "FitterGui", this, "DoChangeBin()");
    
    Int_t Max;
    
    if(Binning)
    {
        Max = 51;
    }
    else
    {
        Max = 37;//Max LBNL Vis binning
    }
    
    hslider->SetRange(0,Max-1);
    hslider->SetPosition(0);
    
    hframe->AddFrame(hslider, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    hframe->MoveResize(0,800,800,30);
    
    fPlotFractionFrame->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    fPlotFractionFrame->AddFrame(fFractionCanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY, 10,10,10,1));
    
    // Map all subwindows of main frame
    fPlotFractionFrame->MapSubwindows();
    
    // Initialize the layout algorithm
    fPlotFractionFrame->Resize(fPlotFractionFrame->GetDefaultSize());
    
    // Map main frame
    fPlotFractionFrame->MapWindow();
}
void FitterGui::PlotNear()
{
    flagNear=1;
    
    Int_t VNearADs;
    
    if(NADs==6)
    {
        VNearADs  = 3;
    }
    else
    {
        VNearADs = 4;
    }
    
    const Int_t NearADs = VNearADs;
    const Int_t Reactors = 6;
    
    FractionCanvas->Clear();
    FractionCanvas->Resize();
    FractionCanvas->Divide(Reactors,NearADs);
    
    TH1D* NearH[NearADs][Reactors];
    
    TFile* NearF = new TFile(("./RootOutputs/"+AnalysisString+"/Spectra/NearSpectrumFraction.root").c_str());
    
    for(Int_t i = 0; i < NearADs; i++)//6 = reactors
    {
        for(Int_t r = 0; r < Reactors; r++)//6 = reactors
        {
            NearH[i][r] = (TH1D*)NearF->Get(Form("AD%d Near Spectrum fraction from Reactor%d Vis%d",i+1,r+1,PlotBin+1));
        }
    }
    
    delete NearF;
    
    for(Int_t j = 0; j < Reactors; j++)
    {
        for(Int_t i = 0; i < NearADs; i++)
        {
            FractionCanvas->cd(j+i*Reactors+1);
            
            NearH[i][j]->Draw();
        }
    }
    
    FractionCanvas->Update();
    
    for(Int_t i = 0; i < NearADs; i++)//6 = reactors
    {
        for(Int_t r = 0; r < Reactors; r++)//6 = reactors
        {
            delete NearH[i][r];
        }
    }
}

void FitterGui::PlotFar()
{
    flagNear=0;
    
    const Int_t FarADs = 3;
    const Int_t NearADs = 3;
    
    FractionCanvas->Clear();
    FractionCanvas->Resize();
    FractionCanvas->Divide(NearADs,FarADs);
    
    TH1D* FarH[NearADs][FarADs];
    
    TFile* FarF = new TFile(("./RootOutputs/"+AnalysisString+"/Spectra/FarSpectrumFraction.root").c_str());
    
    for(Int_t i = 0; i < NearADs; i++)
    {
        for(Int_t j = 0; j < FarADs; j++)
        {
            FarH[i][j] = (TH1D*)FarF->Get(Form("AD%i Far Spectrum prediction from near AD%i Vis %d",j+1,i+1,PlotBin+1));
        }
    }
    
    delete FarF;
    
    for(Int_t i = 0; i < NearADs; i++)
    {
        for(Int_t j = 0; j < FarADs; j++)
        {
            FractionCanvas->cd(j+i*FarADs+1);
            
            FarH[i][j]->Draw();
        }
    }
    
    FractionCanvas->Update();
    
    for(Int_t i = 0; i < NearADs; i++)//6 = reactors
    {
        for(Int_t j = 0; j < FarADs; j++)//6 = reactors
        {
            delete FarH[i][j];
        }
    }
    
}

void FitterGui::PlotNearADTrue()
{
    flagFar = 0;
    
    Int_t VNearADs;
    
    if(NADs==6)
    {
        VNearADs = 3;
    }
    else
    {
        VNearADs = 4;
    }
    
    const Int_t NearADs = VNearADs;
    const Int_t Reactors = 6;
    
    ADTrueCanvas->Clear();
    ADTrueCanvas->Resize();
    ADTrueCanvas->Divide(Reactors,NearADs);
    
    TH1D* NearFractionH[NearADs][Reactors];
    TH1D* NearTrueH[NearADs][Reactors];
    
    TFile* NearF = new TFile(("./RootOutputs/"+AnalysisString+"/Spectra/NearSpectrumFraction.root").c_str());
    
    TH1D* NbinsH = (TH1D*)NearF->Get(Form("AD%d Near Spectrum fraction from Reactor%d Vis%d",1,1,1));
    
    for(Int_t i = 0; i < NearADs; i++)
    {
        for(Int_t r = 0; r < Reactors; r++)
        {
            for(Int_t Bin = 0; Bin < NbinsH->GetXaxis()->GetNbins(); Bin++)
            {
                NearFractionH[i][r] = (TH1D*)NearF->Get(Form("AD%d Near Spectrum fraction from Reactor%d Vis%d",i+1,r+1,Bin+1));
                if(Bin==0)
                {
                    NearTrueH[i][r] = (TH1D*)NearFractionH[i][r]->Clone();
                    NearTrueH[i][r]->SetTitle(Form("AD%d Near Spectrum fraction from Reactor%d",i+1,r+1));
                }
                else
                {
                    NearTrueH[i][r]->Add(NearFractionH[i][r]);
                }
            }
        }
    }
    
    delete NbinsH;
    
    delete NearF;
    
    
    for(Int_t i = 0; i < NearADs; i++)
    {
        for(Int_t r = 0; r < Reactors; r++)//6 = reactors
        {
            ADTrueCanvas->cd(r+Reactors*i+1);
            NearTrueH[i][r]->SetStats(kTRUE);
            NearTrueH[i][r]->Draw();
        }
    }
    
    ADTrueCanvas->Update();
    
    for(Int_t i = 0; i < NearADs; i++)//6 = reactors
    {
        for(Int_t r = 0; r < Reactors; r++)//6 = reactors
        {
            delete NearTrueH[i][r];
            
            delete NearFractionH[i][r];
        }
    }
}


void FitterGui::PlotFarADTrue()
{
    flagFar = 1;
    Int_t VNearADs,VFarADs;
    if(NADs==6)
    {
        VNearADs = 3;
        VFarADs = 3;
    }
    else
    {
        VNearADs = 4;
        VFarADs = 4;
    }
    const Int_t NearADs = VNearADs;
    const Int_t FarADs = VFarADs;
    const Int_t Reactors = 6;
    
    ADTrueCanvas->Clear();
    ADTrueCanvas->Resize();
    ADTrueCanvas->Divide(Reactors,FarADs);
    
    TH1D* FarH[FarADs][Reactors];
    TH1D* FarTrueH[FarADs][Reactors];
    
    TFile* FarF = new TFile(("./RootOutputs/"+AnalysisString+"/Spectra/FarSpectrumFraction.root").c_str());
    
    TH1D* NbinsH = (TH1D*)FarF->Get(Form("AD%i Far Spectrum fraction from Reactor%i and near AD%i, Week%i Vis%i",NearADs+1,1,1,1,1));
    
    for(Int_t j = 0; j < NearADs; j++)
    {
        for(Int_t r = 0; r < Reactors; r++)//6 = reactors
        {
            for(Int_t Bin = 0; Bin < NbinsH->GetXaxis()->GetNbins(); Bin++)
            {
                FarH[j][r] = (TH1D*)FarF->Get(Form("AD%i Far Spectrum fraction from Reactor%i and near AD%i, Week%i Vis%i", j+(NearADs)+1, r+1, NearTrueIndex+1, Period,Bin+1));
                
                if(Bin==0)
                {
                    FarTrueH[j][r] = (TH1D*)FarH[j][r]->Clone();
                    FarTrueH[j][r]->SetTitle(Form("Far AD%d Prediction from near AD%i Reactor%d fraction",j+1,NearTrueIndex+1,r+1));
                }
                else
                {
                    FarTrueH[j][r]->Add(FarH[j][r]);
                }
            }
        }
    }
    
    delete NbinsH;
    
    delete FarF;
    
    for(Int_t j = 0; j < FarADs; j++)
    {
        for(Int_t r = 0; r < Reactors; r++)//6 = reactors
        {
            ADTrueCanvas->cd(r+j*Reactors+1);
            
            FarTrueH[j][r]->Draw();
        }
    }
    
    ADTrueCanvas->Update();
    
    for(Int_t r = 0; r < Reactors; r++)//6 = reactors
    {
        for(Int_t j = 0; j < FarADs; j++)
        {
            delete FarH[j][r];
        }
    }
    
}

void FitterGui::PlotCov()
{
    TGMainFrame *fPlotCovFrame = new TGMainFrame(gClient->GetRoot(),800,800,kMainFrame | kVerticalFrame);
    fPlotCovFrame->SetName("fPlotCov");
    fPlotCovFrame->SetWindowName(" Plot Covariance Matrix ");
    //    fPlotADVisFrame->SetLayoutBroken(kTRUE);
    
    fCovCanvas = new TRootEmbeddedCanvas("CovCanvas",fPlotCovFrame,800,800);
    CovCanvas = fCovCanvas->GetCanvas();
    
    TGHorizontalFrame *hframe = new TGHorizontalFrame(fPlotCovFrame,1000,40);
    
    CovariancePlotBox = new TGComboBox(fPlotCovFrame,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    CovariancePlotBox->SetName("CovarianceMatrixBox");
    
    CovariancePlotBox->AddEntry("Vary Accidental Matrix",1);
    CovariancePlotBox->AddEntry("Vary LiHe Matrix",2);
    CovariancePlotBox->AddEntry("Vary FN Matrix",3);
    CovariancePlotBox->AddEntry("Vary AmC Matrix",4);
    CovariancePlotBox->AddEntry("Distort LiHe Matrix",5);
    CovariancePlotBox->AddEntry("Distort FN Matrix",6);
    CovariancePlotBox->AddEntry("Distort AmC Matrix",7);
    CovariancePlotBox->AddEntry("Reactor Spectrum Matrix",8);
    CovariancePlotBox->AddEntry("Reactor Power Matrix",9);
    CovariancePlotBox->AddEntry("IAV Matrix",10);
    CovariancePlotBox->AddEntry("NL Matrix",11);
    CovariancePlotBox->AddEntry("Resolution Matrix",12);
    CovariancePlotBox->AddEntry("Efficiency Matrix",13);
    CovariancePlotBox->AddEntry("Sin^{2}(2#theta_{12}) Matrix",14);
    CovariancePlotBox->AddEntry("Relative Energy Scale Matrix",15);
    CovariancePlotBox->AddEntry("Statistical Covariance Matrix",16);
    CovariancePlotBox->AddEntry("Background Covariance Matrix",17);
    CovariancePlotBox->AddEntry("Systematic Covariance Matrix",18);
    CovariancePlotBox->AddEntry("Total Covariance Matrix",19);
    
    CovariancePlotBox->Connect("Selected(Int_t)", "FitterGui", this, "ChoosePlotCov()");
    CovariancePlotBox->Select(19);//Run all covariance matrices as default
    CovariancePlotBox->MoveResize(0,0,200,22);
    fPlotCovFrame->AddFrame(CovariancePlotBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    
    hframe->MoveResize(0,800,800,30);
    
    CorrelationBox = new TGButtonGroup(hframe,"Plot Covariance/Correlation",kVerticalFrame,TGLabel::GetDefaultGC()(),TGLabel::GetDefaultFontStruct(),0xffffff);
    CorrelationBox->SetName("CorrelationBox");
    fC[0] = new TGRadioButton(CorrelationBox,"Plot Covariance Matrix",0);
    fC[1] = new TGRadioButton(CorrelationBox,"Plot Correlation Matrix",1);
    
    fC[0]->ChangeBackground(0xffffff);
    fC[1]->ChangeBackground(0xffffff);
    
    fC[0]->SetState(kButtonDown);
    
    for (Int_t i = 0; i < 2; ++i)
    {
        fC[i]->Connect("Pressed()", "FitterGui",  this, "ChoosePlotCovariance()");
    }
    
    CorrelationBox->Show();
    
    hframe->AddFrame(CorrelationBox, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    CorrelationBox->MoveResize(400,800,400,15);
    
    fPlotCovFrame->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    fPlotCovFrame->AddFrame(fCovCanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY, 10,10,10,1));
    
    // Map all subwindows of main frame
    fPlotCovFrame->MapSubwindows();
    
    // Initialize the layout algorithm
    fPlotCovFrame->Resize(fPlotCovFrame->GetDefaultSize());
    
    // Map main frame
    fPlotCovFrame->MapWindow();
    
}

void FitterGui:: ChoosePlotVariations()
{
    TGButton *btn = (TGButton *) gTQSender;
    Int_t id = btn->WidgetId();
    
    if (id == 0)
    {
        NominalString = "Ratio";
        fV[id]->SetState(kButtonUp);
    }
    else
    {
        NominalString = "Spectrum";
        fV[id]->SetState(kButtonUp);
    }
    
    if(id == 0)//Update the current active window
    {
        ChoosePlotRatioVariations();
    }
    else
    {
        ChoosePlotSpectrumVariations();
    }
    
    std::cout << "Plot Ratio/Spectrum : " << NominalString  << std::endl;
    
}

void FitterGui:: ChoosePlotCovariance()
{
    TGButton *btn = (TGButton *) gTQSender;
    Int_t id = btn->WidgetId();
    
    if (id == 0)
    {
        CorrelationString = "Covariance";
        fC[id]->SetState(kButtonUp);
    }
    else
    {
        CorrelationString = "Correlation";
        fC[id]->SetState(kButtonUp);
    }
    
    if(id == 0)//Update the current active window
    {
        PlotCovariance = 1;
    }
    else
    {
        PlotCovariance = 0;
    }
    
    ChoosePlotCov();
    
    std::cout << "Plot Covariance/Correlation: " << CorrelationString  << std::endl;
    
}

void FitterGui::ChooseVariations()
{
    switch (VariationsPlotBox->GetSelected())
    {
        case 0:
            break;
        case 1:
            VarString = "VaryAccidentals";
            std::cout << " Vary Accidental Matrix " << std::endl;
            break;
        case 2:
            VarString = "VaryLiHe";
            std::cout << " Vary LiHe Matrix " << std::endl;
            break;
        case 3:
            VarString = "VaryFastNeutrons";
            std::cout << " Vary Fast Neutrons Matrix " << std::endl;
            break;
        case 4:
            VarString = "VaryAmC";
            std::cout << " Vary AmC Matrix " << std::endl;
            break;
        case 5:
            VarString = "DistortLiHe";
            std::cout << " Distort LiHe Matrix " << std::endl;
            break;
        case 6:
            VarString = "DistortFastNeutrons";
            std::cout << " Distort Fast Neutrons Matrix " << std::endl;
            break;
        case 7:
            VarString = "DistortAmC";
            std::cout << " Distort AmC Matrix " << std::endl;
            break;
        case 8:
            VarString = "Isotope";
            std::cout << " Reactor Spectrum Matrix " << std::endl;
            break;
        case 9:
            VarString = "Power";
            std::cout << " Reactor Power Matrix " << std::endl;
            break;
        case 10:
            VarString = "IAV";
            std::cout << " IAV Matrix " << std::endl;
            break;
        case 11:
            VarString = "NL";
            std::cout << " NL Matrix " << std::endl;
            break;
        case 12:
            VarString = "Resolution";
            std::cout << " Resolution Matrix " << std::endl;
            break;
        case 13:
            VarString = "Efficiency";
            std::cout << " Efficiency Matrix " << std::endl;
            break;
        case 14:
            VarString = "Sin22t12";
            std::cout << " Sin22t12 Matrix " << std::endl;
            break;
        case 15:
            VarString = "RelativeEnergyScale";
            std::cout << " Energy Scale Matrix " << std::endl;
            break;
        default:
            std::cout << " Plotting default: NL" << std::endl;
            break;
    }
}

void FitterGui:: ChoosePlotRatioVariations()
{
    if(CombineMode==0)
    {
        MaxFar=3;
    }
    else if(CombineMode ==1)
    {
        MaxFar=1;
    }
    else if(CombineMode ==2)
    {
        MaxFar=1;
    }
    
    if(deleteFlag == 1)
    {
        for(Int_t far = 0; far< MaxFar; far++)
        {
            for(Int_t near = 0; near< MaxNear; near++)
            {
                for(Int_t sample = 0; sample< NSamples; sample++)
                {
                    delete RatioH[sample][far][near];
                }
            }
        }
        
        deleteFlag = 0;
    }
    
    TFile* VarF = new TFile(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/Spectrum/",CombineMode)+VarString+".root").c_str());
    
    VCanvas->Clear();
    
    if(VariationsPlotBox->GetSelected()<8)
    {
        VCanvas->SetCanvasSize(400,400*MaxFar);
        VCanvas->Divide(1,MaxFar);
        
        for(Int_t far = 0; far< MaxFar; far++)
        {
            for(Int_t near = 0; near< MaxNear; far++)
            {
                for(Int_t sample = 0; sample< NSamples; sample++)
                {
                    RatioH[sample][far][near] = (TH1D*)VarF->Get(Form("Relative error AD%d Spectra Sample%i Period%d",far, sample,Period-1));
                }
            }
        }
    }
    else
    {
        VCanvas->SetCanvasSize(400*MaxNear,400*MaxFar);
        VCanvas->Divide(MaxNear,MaxFar);
        
        for(Int_t far = 0; far< MaxFar; far++)
        {
            for(Int_t near = 0; near< MaxNear; near++)
            {
                for(Int_t sample = 0; sample< NSamples; sample++)
                {
                    RatioH[sample][far][near] = (TH1D*)VarF->Get(Form("Relative error Far AD%d from Near AD%d Spectra Sample%i Period%d",far,near, sample,Period-1));
                }
            }
        }
    }
    
    gStyle->SetOptStat(SetStats);
    
    Int_t MaxPlotNear;
    
    if(VariationsPlotBox->GetSelected()<8)
    {
        MaxPlotNear=1;
    }
    else
    {
        MaxPlotNear=MaxNear;
    }
    
    for(Int_t far = 0; far< MaxFar; far++)
    {
        for(Int_t near = 0; near < MaxPlotNear; near++)
        {
            VCanvas->cd(near+far*MaxPlotNear+1);
            
            RatioH[0][far][near]->SetTitle((VarString+" Ratio").c_str());
            RatioH[0][far][near]->GetYaxis()->SetRangeUser(-0.5,0.5);
            
            for(Int_t sample = 0; sample< NSamples; sample++)
            {
                RatioH[sample][far][near]->Draw("HIST same");
            }
        }
    }
    
    delete VarF;
    
    RatioH[0][0][0]->Reset();
    RatioH[0][0][0]->SetLineColor(kRed);
    RatioH[0][0][0]->Draw("HIST same");
    
    
    VCanvas->Update();
    
    VCanvas->Print(("./Images/"+AnalysisString+Form("/Variations/RatioCombine%d",CombineMode)+VarString+".eps").c_str(),".eps");
    deleteFlag = 1;
}

void FitterGui:: ChoosePlotSpectrumVariations()
{
    if(CombineMode==0)
    {
        MaxFar=3;
        MaxNear=3;
    }
    else if(CombineMode ==1)
    {
        MaxFar=1;
        MaxNear=1;
    }
    else if(CombineMode ==2)
    {
        MaxFar=1;
        MaxNear=2;
    }
    
    if(VariationsPlotBox->GetSelected()<8)
    {
        MaxNear = 1;
    }
    
    if(deleteFlagSpec == 1)
    {
        for(Int_t far = 0; far< MaxFar; far++)
        {
            for(Int_t near = 0; near< MaxNear; near++)
            {
                delete NominalSpectrumH[far][near];
                
                for(Int_t sample = 0; sample< NSamples; sample++)
                {
                    delete SpectrumH[sample][far][near];
                }
            }
        }
        
        deleteFlagSpec = 0;
    }
    
    TFile* VarF = new TFile(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/Spectrum/",CombineMode)+VarString+".root").c_str());
    
    VCanvas->Clear();
    
    if(VariationsPlotBox->GetSelected()<8)
    {
        VCanvas->SetCanvasSize(400,400*MaxFar);
        VCanvas->Divide(1,MaxFar);
        
        for(Int_t far = 0; far< MaxFar; far++)
        {
            for(Int_t sample = 0; sample< NSamples; sample++)
            {
                SpectrumH[sample][far][0] = (TH1D*)VarF->Get(Form("Difference AD%d Spectra Sample%i Period%d",far,sample, Period-1));
                
                if(SpectrumH[sample][far][0]->GetBinContent(0)<0)
                {
                    SpectrumH[sample][far][0]->Scale(-1);
                }
            }
            
            NominalSpectrumH[far][0] = (TH1D*)VarF->Get(Form("Far AD%d Varied Background Spectra Sample0 Period%d",far, Period-1));
            NominalSpectrumH[far][0]->Reset();
        }
    }
    else
    {
        VCanvas->SetCanvasSize(400*MaxNear,400*MaxFar);
        VCanvas->Divide(MaxNear,MaxFar);
        
        for(Int_t far = 0; far< MaxFar; far++)
        {
            for(Int_t near = 0; near< MaxNear; near++)
            {
                NominalSpectrumH[far][near] = (TH1D*)VarF->Get(Form("Nominal Far AD%d prediction from Near AD%d VisH Sample0 Period%d",far, near, Period-1));
                
                for(Int_t sample = 0; sample< NSamples; sample++)
                {
                    SpectrumH[sample][far][near] = (TH1D*)VarF->Get(Form("Varied Far AD%d prediction from Near AD%d VisH Sample%i Period%d",far, near,sample, Period-1));
                }
            }
        }
    }
    
    delete VarF;
    
    SpectrumH[0][0][0]->SetTitle((VarString+" Spectrum").c_str());
    gStyle->SetOptStat(SetStats);
    Int_t MaxPlotNear;
    
    if(VariationsPlotBox->GetSelected()<8)
    {
        SpectrumH[0][0][0]->GetYaxis()->SetRangeUser(-0.5,0.5);
        MaxPlotNear=1;
    }
    else
    {
        MaxPlotNear=MaxNear;
    }
    for(Int_t far = 0; far< MaxFar; far++)
    {
        for(Int_t near = 0; near< MaxPlotNear; near++)
        {
            VCanvas->cd(near+far*MaxPlotNear+1);
            
            for(Int_t sample = 0; sample< NSamples; sample++)
            {
                SpectrumH[sample][far][near]->Draw("HIST same");
            }
            
            NominalSpectrumH[far][near]->SetLineColor(kRed);
            NominalSpectrumH[far][near]->Draw("HIST same");
        }
    }
    
    
    
    VCanvas->Update();
    VCanvas->Print(("./Images/"+AnalysisString+Form("/Variations/SpectrumCombine%d",CombineMode)+VarString+".eps").c_str(),".eps");
    
    deleteFlagSpec = 1;
}

void FitterGui::ChoosePlotCov()
{
    TH2D* CovMatrix2H;
    TFile* CovF;
    
    Int_t DataSet =2;
    NominalData* Data = new NominalData(Analysis,DataSet);
    
    Double_t sen22t13 = Data->GetSin22t13();
    Double_t dm2_ee = Data->GetDm2ee();
    delete Data;
    string FitterCovS;
    
    switch (CovariancePlotBox->GetSelected())
    {
        case 0:
            std::cout << " No Matrix Selected " << std::endl;
            break;
        case 1:
            CovString = "VaryAccidental";
            FitterCovS = "Vacc Matrix";
            std::cout << " Vary Accidental Matrix " << std::endl;
            break;
        case 2:
            CovString = "VaryLiHe";
            FitterCovS = "VLiHe Matrix";
            std::cout << " Vary LiHe Matrix " << std::endl;
            break;
        case 3:
            CovString = "VaryFastNeutrons";
            FitterCovS = "VFN Matrix";
            std::cout << " Vary Fast Neutrons Matrix " << std::endl;
            break;
        case 4:
            CovString = "VaryAmC";
            FitterCovS = "VAmC Matrix";
            std::cout << " Vary AmC Matrix " << std::endl;
            break;
        case 5:
            CovString = "DistortLiHe";
            FitterCovS = "DLiHe Matrix";
            std::cout << " Distort LiHe Matrix " << std::endl;
            break;
        case 6:
            CovString = "DistortFastNeutrons";
            FitterCovS = "DFN Matrix";
            std::cout << " Distort Fast Neutrons Matrix " << std::endl;
            break;
        case 7:
            CovString = "DistortAmC";
            FitterCovS = "DAmC Matrix";
            std::cout << " Distort AmC Matrix " << std::endl;
            break;
        case 8:
            CovString = "Isotope";
            FitterCovS = "Isotope Matrix";
            std::cout << " Reactor Spectrum Matrix " << std::endl;
            break;
        case 9:
            CovString = "ReactorPower";
            FitterCovS = "Power Matrix";
            std::cout << " Reactor Power Matrix " << std::endl;
            break;
        case 10:
            CovString = "IAV";
            FitterCovS = "IAV Matrix";
            std::cout << " IAV Matrix " << std::endl;
            break;
        case 11:
            CovString = "NL";
            FitterCovS = "NL Matrix";
            std::cout << " NL Matrix " << std::endl;
            break;
        case 12:
            CovString = "Resolution";
            FitterCovS = "Reso Matrix";
            std::cout << " Resolution Matrix " << std::endl;
            break;
        case 13:
            CovString = "Efficiency";
            FitterCovS = "Efficiency Matrix";
            std::cout << " Efficiency Matrix " << std::endl;
            break;
        case 14:
            CovString = "Sin22t12";
            FitterCovS = "Sin22t12 Matrix";
            std::cout << " Sin22t12 Matrix " << std::endl;
            break;
        case 15:
            CovString = "RelativeEnergyScale";
            FitterCovS = "Relative Scale Matrix";
            std::cout << " Energy Scale Matrix " << std::endl;
            break;
        case 16:
            CovString = "StatisticalCovarianceMatrix";
            FitterCovS = "Statistical Covariance Matrix";
            CovF = new TFile(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/StatisticalCovarianceMatrix.root",CombineMode)).c_str());
            CovMatrix2H = (TH2D*)CovF->Get(Form("Statistical Covariance Matrix for sin22t13 %f and Δm2ee %f period%d",sen22t13,dm2_ee,Period-1));
            std::cout << " Statistical Covariance Matrix " << std::endl;
            break;
        case 17:
            CovString = "BackgroundCovarianceMatrix";
            FitterCovS = "Background Covariance Matrix";
            CovF = new TFile(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/BackgroundCovarianceMatrices.root",CombineMode)).c_str());
            CovMatrix2H = (TH2D*)CovF->Get("Background Covariance Matrix");
            std::cout << " Background Covariance Matrix " << std::endl;
            break;
        case 18:
            CovString = "SystematicCovarianceMatrix";
            FitterCovS = "Systematic Covariance Matrix";
            CovF = new TFile(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/SystematicCovarianceMatrices.root",CombineMode)).c_str());
            CovMatrix2H = (TH2D*)CovF->Get("Systematic Covariance Matrix");
            std::cout << " Systematic Covariance Matrix " << std::endl;
            break;
        case 19:
            CovString = "TotalCovarianceMatrix";
            FitterCovS = "Total Covariance Matrix";
            CovF = new TFile(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/TotalCovarianceMatrix.root",CombineMode)).c_str());
            CovMatrix2H = (TH2D*)CovF->Get("Total Covariance Matrix");
            std::cout << " Total Covariance Matrix " << std::endl;
            break;
        default:
            std::cout << " No Matrix Selected " << std::endl;
            break;
    }
    
    if(CovariancePlotBox->GetSelected()<16)
    {
        CovF = new TFile(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesRoot/",CombineMode)+CovString+".root").c_str());
        Char_t File[30];
        if(CovariancePlotBox->GetSelected()<8)//Background
        {
            sprintf(File,"Covariance Matrix%d",Period-1);
        }
        else//Systematic, need rescaling
        {
            sprintf(File,"After Covariance Matrix%d",Period-1);
        }
        
        CovMatrix2H = (TH2D*)CovF->Get(File);
        
        //        if(Print)
        //        {
        //            TH1D* DiagonalCov = new TH1D(("sigma "+CovString).c_str(),("sigma "+CovString).c_str(),CovMatrix2H->GetXaxis()->GetNbins(),0,CovMatrix2H->GetXaxis()->GetNbins());
        //
        //            for(Int_t i = 0; i < CovMatrix2H->GetXaxis()->GetNbins(); i++)
        //            {
        //                DiagonalCov->SetBinContent(i+1,CovMatrix2H->GetBinContent(i+1,i+1));//Diagonal of the covariance matrix = sigma
        //            }
        //
        //            TCanvas* sigmaC = new TCanvas("sigma","sigma");
        //
        //            DiagonalCov->Draw();
        //
        //            sigmaC->Print((Form("./Images/Sigma_Combine%d_",CombineMode)+CovString+".eps").c_str(),".eps");
        //
        //            delete sigmaC;
        //
        //            delete DiagonalCov;
        //        }
        
        delete CovF;
        
        delete CovMatrix2H;
    }
    
    if(PlotCovariance)
    {
        CovF = new TFile(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CovarianceMatricesFitterPeriod%d.root",CombineMode,Period-1)).c_str());
    }
    else
    {
        CovF = new TFile(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CorrelationMatrices.root",CombineMode)).c_str());
    }
    
    CovMatrix2H = (TH2D*)CovF->Get((FitterCovS).c_str());
    
    TPad* centerpad = new TPad("upperPad", "upperPad", 0,0,1,1);
    centerpad->Draw();
    centerpad->cd();
    
    CovMatrix2H->SetStats(kFALSE);
    if(PlotCovariance&&CombineMode!=0)
    {
        CovMatrix2H->GetYaxis()->SetRange(2,36);
        CovMatrix2H->GetXaxis()->SetRange(2,36);
    }
    
    CovMatrix2H->Draw("colz");
    CovMatrix2H->Draw("colz");
    CovMatrix2H->Draw("colz");
    CovMatrix2H->SetTitle("");
    
    centerpad->Update();
    
    TPaletteAxis *palette = (TPaletteAxis*)CovMatrix2H->GetListOfFunctions()->FindObject("palette");
    
    palette->SetLabelSize(.025);
    palette->SetX1NDC(0.902);
    palette->SetX2NDC(0.925);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
    
    centerpad->Modified();
    centerpad->Update();
    
    CovCanvas->cd();
    CovCanvas->Update();
    
    if(PlotCovariance)
    {
        CovCanvas->Print(("./Images/"+AnalysisString+Form("/CovarianceMatrices/Combine%d",CombineMode)+CovString+".eps").c_str(),".eps");
    }
    else
    {
        CovCanvas->Print(("./Images/"+AnalysisString+Form("/CovarianceMatrices/Combine%dCorrelation",CombineMode)+CovString+".eps").c_str());
    }
    
    delete CovF;
    
}
//
//void FitterGui::ChoosePlotCorr()
//{
//    if(Print)
//    {
//        TFile *f1 = new TFile(("./CovarianceMatrices/"+AnalysisString+Form("/Combine%d/CorrelationMatrices.root",CombineMode)).c_str());
//
//        TIter next(f1->GetListOfKeys());
//        TKey *key;
//
//        while ((key = (TKey*)next()))
//        {
//            TCanvas* CorrelationC = new TCanvas("");
//
//            TH1D* h = (TH1D*)key->ReadObj();
//            string TCorName = key->GetTitle();
//
//            Int_t size = (Int_t)TCorName.length();
//
//            for(Int_t j = 0; j<=size; j++)
//            {
//                for(Int_t i = 0; i <=j; i++)
//                {
//                    delSpaces(TCorName);
//                }
//            }
//            h->SetStats(0);
//            h->SetTitle("");
//            h->Draw("colz");
//
//            CorrelationC->Print(("./Images/"+AnalysisString+Form("/CovarianceMatrices/Combine%dCorrelation",CombineMode)+TCorName+".eps").c_str());
//            delete h;
//            delete CorrelationC;
//        }
//    }
//}

string FitterGui :: delSpaces(string &str)
{
    str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
    return str;
}

void FitterGui::PlotResponseMatrix()
{
    RunResponseMatrix();
    
    TGMainFrame *fPlotResponseFrame = new TGMainFrame(gClient->GetRoot(),800,800,kMainFrame | kVerticalFrame);
    fPlotResponseFrame->SetName("fResponse");
    fPlotResponseFrame->SetWindowName(" Plot Response Matrix");
    //    fPlotADVisFrame->SetLayoutBroken(kTRUE);
    
    fResponseCanvas = new TRootEmbeddedCanvas("ResponseCanvas",fPlotResponseFrame,800,800);
    ResponseCanvas = fResponseCanvas->GetCanvas();
    
    TGHorizontalFrame *hframe = new TGHorizontalFrame(fPlotResponseFrame,800,40);
    
    ResponseBox = new TGComboBox(fPlotResponseFrame,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    ResponseBox->SetName("ResponseBox");
    ResponseBox->AddEntry("Evis-Enu Matrix",0);
    ResponseBox->AddEntry("Enu-Evis Matrix",1);
    ResponseBox->AddEntry("Evis-Enu Positron Matrix",2);
    ResponseBox->AddEntry("Evis-Enu IAV Matrix",3);
    ResponseBox->AddEntry("Evis-Enu NL Matrix ",4);
    ResponseBox->AddEntry("Evis-Enu Resolution Matrix",5);
    ResponseBox->Resize(180,20);
    ResponseBox->Connect("Selected(Int_t)", "FitterGui", this, "ChooseResponseMatrix()");
    ResponseBox->Select(1);//Enu-Evis Default right now.
    
    hframe->AddFrame(ResponseBox, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    
    hframe->MoveResize(0,800,800,40);
    
    fPlotResponseFrame->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    fPlotResponseFrame->AddFrame(fResponseCanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY, 10,10,10,1));
    
    // Map all subwindows of main frame
    fPlotResponseFrame->MapSubwindows();
    
    // Initialize the layout algorithm
    fPlotResponseFrame->Resize(fPlotResponseFrame->GetDefaultSize());
    
    // Map main frame
    fPlotResponseFrame->MapWindow();
}

void FitterGui :: ChooseResponseMatrix()
{
    TH2D* EnergyMatrix;
    
    TFile* TransEnergyMatrixDataF = TFile::Open(("./ResponseMatrices/"+AnalysisString+"/NominalResponseMatrix.root").c_str());
    
    switch (ResponseBox->GetSelected())
    {
        case 0:
            EnergyMatrix = (TH2D*)TransEnergyMatrixDataF->Get("EvisEnu");
            std::cout << " Evis - Enu Matrix " << std::endl;
            break;
        case 1:
            EnergyMatrix = (TH2D*)TransEnergyMatrixDataF->Get("EnuEvis");
            std::cout << " Enu - Evis Matrix " << std::endl;
            break;
        case 2:
            EnergyMatrix = (TH2D*)TransEnergyMatrixDataF->Get("EvisEnuPos");
            std::cout << " Evis - Enu Positron Matrix " << std::endl;
            break;
        case 3:
            EnergyMatrix = (TH2D*)TransEnergyMatrixDataF->Get("EvisEnuIAV");
            std::cout << " Evis - Enu IAV Matrix " << std::endl;
            break;
        case 4:
            EnergyMatrix = (TH2D*)TransEnergyMatrixDataF->Get("EvisEnuNL");
            std::cout << " Evis - Enu NL Matrix " << std::endl;
            break;
        case 5:
            EnergyMatrix = (TH2D*)TransEnergyMatrixDataF->Get("EvisEnuReso");
            std::cout << " Evis - Enu Resolution Matrix " << std::endl;
            break;
        default:
            std::cout << " No Matrix Selected " << std::endl;
            break;
    }
    
    TransEnergyMatrixDataF->Close();
    gStyle->SetOptStat(SetStats);
    EnergyMatrix->Draw("colz");
    
    ResponseCanvas->Update();
}

void FitterGui :: ChooseCovarianceMatrix()
{
    Automatic = 0;
    VaryAccidentalMatrix=0;
    VaryLiHeMatrix=0;
    VaryFastNeutronsMatrix=0;
    VaryAmCMatrix=0;
    DistortLiHeMatrix=0;
    DistortFastNeutronsMatrix=0;
    DistortAmCMatrix=0;
    //Systematics
    IsotopeMatrix=0;
    ReactorPowerMatrix=0;
    EnergyScaleMatrix=0;
    //         EnergyOffsetMatrix=0;
    //         AbsoluteScaleMatrix=0;
    //         AbsoluteOffsetMatrix=0;
    IAVMatrix=0;
    NLMatrix=0;
    ResolutionMatrix=0;
    Sin22t12Matrix=0;
    EfficiencyMatrix=0;
    
    switch (CovarianceMatrixBox->GetSelected())
    {
        case 0:
            Automatic = 1;
            std::cout << " Run all Matrices " << std::endl;
            break;
        case 1:
            VaryAccidentalMatrix=1;
            std::cout << " Vary Accidental Matrix " << std::endl;
            break;
        case 2:
            VaryLiHeMatrix=1;
            std::cout << " Vary LiHe Matrix " << std::endl;
            break;
        case 3:
            VaryFastNeutronsMatrix=1;
            std::cout << " Vary Fast Neutrons Matrix " << std::endl;
            break;
        case 4:
            VaryAmCMatrix=1;
            std::cout << " Vary AmC Matrix " << std::endl;
            break;
        case 5:
            DistortLiHeMatrix=1;
            std::cout << " Distort LiHe Matrix " << std::endl;
            break;
        case 6:
            DistortFastNeutronsMatrix=1;
            std::cout << " Distort Fast Neutrons Matrix " << std::endl;
            break;
        case 7:
            DistortAmCMatrix=1;
            std::cout << " Distort AmC Matrix " << std::endl;
            break;
        case 8:
            IsotopeMatrix=1;
            std::cout << " Reactor Spectrum Matrix " << std::endl;
            break;
        case 9:
            ReactorPowerMatrix=1;
            std::cout << " Reactor Power Matrix " << std::endl;
            break;
        case 10:
            IAVMatrix=1;
            std::cout << " IAV Matrix " << std::endl;
            break;
        case 11:
            NLMatrix=1;
            std::cout << " NL Matrix " << std::endl;
            break;
        case 12:
            ResolutionMatrix=1;
            std::cout << " Resolution Matrix " << std::endl;
            break;
        case 13:
            EfficiencyMatrix=1;
            std::cout << " Efficiency Matrix " << std::endl;
            break;
        case 14:
            Sin22t12Matrix=1;
            std::cout << " Sin22t12 Matrix " << std::endl;
            break;
        case 15:
            EnergyScaleMatrix=1;
            std::cout << " Energy Scale Matrix " << std::endl;
            break;
        default:
            std::cout << " No Matrix Selected " << std::endl;
            break;
    }
}

void FitterGui :: DoNADs()
{
    NADs = NADsBox->GetSelected();
    
    std::cout << "NUMBER OF DETECTORS : " << NADs << std::endl;
}

void FitterGui :: DoPeriod()
{
    Period =(Int_t)PeriodBox->GetNumberEntry()->GetIntNumber();
    std::cout << "NUMBER OF Period : " << Period << std::endl;
}

void FitterGui :: ChooseToyMC()
{
    ToyMC = ToyMCBox->GetSelected();
    std::cout << "TOY MC INPUT IN FITTER : " << ToyMC << std::endl;
}

void FitterGui :: ChooseMinuit()
{
    Minuit = MinuitBox->GetSelected();
    std::cout << "MINUIT : " << Minuit << std::endl;
}

void FitterGui::RunResponseMatrix()
{
    Int_t DataSet=2;
    
    NominalData* Data = new NominalData(Analysis,DataSet);
    Double_t sin22t13 = Data->GetSin22t13();
    Double_t deltaM = Data->GetDm2ee();
    delete Data;
    
    CreateEnergyMatrix* GenMatrix = new CreateEnergyMatrix(Data);
    
    GenMatrix->GenerateEnergyMatrix(sin22t13,deltaM,Period);//GenerateLBNLEnergyMatrix to try their code
    
    delete GenMatrix;
    
    std::cout << " Finished Generating Response Matrix for sin22t13 :" << sin22t13 << " and deltaM : " << deltaM << std::endl;
}

void FitterGui::PlotChi()
{
    TGMainFrame *fPlotChi = new TGMainFrame(gClient->GetRoot(),800,800,kMainFrame | kVerticalFrame);
    fPlotChi->SetName("fPlotChi");
    fPlotChi->SetWindowName(" Plot ChiSquare Distribution ");
    //    fPlotADVisFrame->SetLayoutBroken(kTRUE);
    
    fChiCanvas = new TRootEmbeddedCanvas("ChiCanvas",fPlotChi,800,800);
    ChiCanvas = fChiCanvas->GetCanvas();
    
    Char_t FileName[100];
    Char_t ObjectName[100];
    
    Int_t DataSet=2;
    
    NominalData* Data = new NominalData(Analysis,DataSet);
    Double_t Sin22t13 = Data->GetSin22t13();
    Double_t Dm2_ee = Data->GetDm2ee();
    delete Data;
    
    sprintf(ObjectName,"Sin%f_Distribution_DeltaM%f_Period%d",Sin22t13,Dm2_ee,Period-1);
    
    if(Fit2D)
    {
        sprintf(FileName,("./ChiSquare/"+AnalysisString+Form("/Combine%d/2DChiSquare.root",CombineMode)).c_str());
        sprintf(ObjectName,"2DChiDistribution");
        
        TGHorizontalFrame *hframe = new TGHorizontalFrame(fPlotChi,1000,40);
        
        Sslider = new TGHSlider(fPlotChi,390);
        
        Sslider->Connect("PositionChanged(Int_t)", "FitterGui", this, "DoChangeSinSquareOscillationParameter()");
        
        Dslider = new TGHSlider(fPlotChi,390);
        
        Dslider->Connect("PositionChanged(Int_t)", "FitterGui", this, "DoChangeDeltaSquareOscillationParameter()");
        
        
        hframe->MoveResize(0,800,800,30);
        
        Sslider->SetRange(0,NFits-1);
        Sslider->SetPosition(0);
        
        Dslider->SetRange(0,NFits-1);
        Dslider->SetPosition(0);
        
        hframe->AddFrame(Sslider, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
        hframe->AddFrame(Dslider, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
        
        hframe->MoveResize(0,800,800,30);
        
        fPlotChi->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
        
    }
    else if(FitSin22t13)
    {
        sprintf(FileName,("./ChiSquare/"+AnalysisString+"/Combine%d/1DChiSquareFitSin22t13.root").c_str(),CombineMode);
    }
    else
    {
        sprintf(FileName,("./ChiSquare/"+AnalysisString+"/Combine%d/1DChiSquareFitDm.root").c_str(),CombineMode);
    }
    
    TFile* ChiFile = new TFile(FileName);
    ChiH = ChiFile->Get(ObjectName);
    gStyle->SetOptStat(SetStats);
    
    if(Fit2D)
    {
        ChiH->Draw("colz");
        ChiCanvas->Print(("./Images/"+AnalysisString+"/ChiSquareResults/2DFit.eps").c_str(),".eps");
    }
    else
    {
        ChiH->Draw();
        ChiCanvas->Print(("./Images/"+AnalysisString+"/ChiSquareResults/1DFit.eps").c_str(),".eps");
    }
    
    delete ChiFile;
    
    fPlotChi->AddFrame(fChiCanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY, 10,10,10,1));
    
    // Map all subwindows of main frame
    fPlotChi->MapSubwindows();
    
    // Initialize the layout algorithm
    fPlotChi->Resize(fPlotChi->GetDefaultSize());
    
    // Map main frame
    fPlotChi->MapWindow();
}

void FitterGui::DoChangeSinSquareOscillationParameter()
{
    TH2D* Histo = (TH2D*)ChiH->Clone();
    TH1D* Histo1D;
    
    Sin22t13SliderValue = Sslider->GetPosition();
    
    Histo1D = Histo->ProjectionX("SinSlice",Sin22t13SliderValue+1,Sin22t13SliderValue+1);
    
    Histo1D->SetStats(kFALSE);
    Histo1D->Draw();
    
    ChiCanvas->Update();
    
    ChiCanvas->Print(("./Images/"+AnalysisString+Form("/ChiSquareResults/SinChi2_%f.eps",Sin22t13SliderValue)).c_str(),".eps");
    
    std::cout << "Fit Distribution for Sin22t13: " << Sin22t13SliderValue << std::endl;
}

void FitterGui::DoChangeDeltaSquareOscillationParameter()
{
    TH2D* Histo = (TH2D*)ChiH->Clone();
    TH1D* Histo1D;
    
    DeltaMSilderValue = Dslider->GetPosition();
    
    Histo1D = Histo->ProjectionY("DeltaSlice",DeltaMSilderValue+1,DeltaMSilderValue+1);
    
    Histo1D->SetStats(kFALSE);
    Histo1D->Draw();
    
    ChiCanvas->Update();
    
    ChiCanvas->Print(("./Images/"+AnalysisString+Form("/ChiSquareResults/DMChi2_%f.eps",DeltaMSilderValue)).c_str(),".eps");
    
    std::cout << "Fit Distribution for DeltaM213: " << Sin22t13SliderValue << std::endl;
}

void FitterGui::GenerateToyMC()
{
    gBenchmark->Start("GenerateTree");
    std::cout << "Running Toy MC Tree Generation" << std::endl;

    if(!PrintOnConsole)
    {
        coutstream = cout.rdbuf(0);//change stream of cout
    }
    //To draw using a better palette:
    const Int_t NRGBs = 5;
    const Int_t NCont = 999;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    gStyle->SetTitleOffset(1.3);
    
    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    
    Int_t n_evis_bins;
    Int_t n_etrue_bins;
    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.
    
    Int_t DataSet=2;
    
    NominalData* NData = new NominalData(Analysis,DataSet);
    
    NData->SetToyMCSamplesDirectory(ToyMCSampleDirectory);
    NData->SetPredictionDirectory(NominalPredictionsDirectory);
    NData->SetResponseDirectory(ResponseMatrixDirectory);
    NData->SetBkgCovDirectory(BkgCovarianceMatrixDirectory);
    NData->SetSysCovDirectory(SysCovarianceMatrixDirectory);
    
    Double_t s22t13start= NData->GetSinStart();
    Double_t s22t13end=NData->GetSinEnd();
    Double_t dm2_eestart=NData->GetDmeeStart();
    Double_t dm2_eeend=NData->GetDmeeEnd();
    
    NominalData* VData = new NominalData(Analysis,DataSet);
    
    VData->SetToyMCSamplesDirectory(ToyMCSampleDirectory);
    VData->SetPredictionDirectory(NominalPredictionsDirectory);
    VData->SetResponseDirectory(ResponseMatrixDirectory);
    VData->SetBkgCovDirectory(BkgCovarianceMatrixDirectory);
    VData->SetSysCovDirectory(SysCovarianceMatrixDirectory);
    
    InitialEnergy = NData->GetEmin();
    FinalEnergy = NData->GetEmax();
    InitialVisibleEnergy = NData->GetEVisMin();
    FinalVisibleEnergy = NData->GetEVisMax();
    
    InitialEnergy = VData->GetEmin();
    FinalEnergy = VData->GetEmax();
    InitialVisibleEnergy = VData->GetEVisMin();
    FinalVisibleEnergy = VData->GetEVisMax();
    //  Linear binning
    if(Binning)
    {
        n_evis_bins = NData->GetNbins();
        n_etrue_bins = NData->GetNbins();
        
        n_evis_bins = VData->GetNbins();
        n_etrue_bins = VData->GetNbins();
        
        for (Int_t i = 0; i <= n_evis_bins; i++)
        {
            evis_bins[i] = 0.2 * i + 0.7;
            enu_bins[i] = 0.2 * i + InitialEnergy;
        }
    }
    //  Non-linear binning
    else
    {
        n_evis_bins=37;
        n_etrue_bins=39;
        
        for (Int_t i = 0; i <= n_etrue_bins; i++)
        {
            enu_bins[i] = 0.2 * i + InitialEnergy;
        }
        
        evis_bins[0] = 0.7;
        for (Int_t i = 0; i < n_evis_bins-1; i++)
        {
            evis_bins[i+1] = 0.2 * i + 1.0;
        }
        evis_bins[n_evis_bins] = FinalVisibleEnergy;
    }
    
    Int_t MaxNearCombine,MaxFarCombine;
    
    if(CombineMode==0)
    {
        MaxNearCombine =3;
        MaxFarCombine =3;
    }
    else if(CombineMode==1)
    {
        MaxNearCombine =1;
        MaxFarCombine =1;
    }
    else
    {
        MaxNearCombine =2;
        MaxFarCombine =1;
    }
    
    Double_t sin22t13;
    Double_t dm2_ee;
    
    NData->SetCombineMode(CombineMode); //0 is 9x9, 1 is 1x1 and 2 is 2x2
    //    Data->SetSin22t12(0);
    //Parameters of the model
    NData->SetAnalysis(Analysis);//  Gd or H data
    NData->SetBinning(Binning);//  0 for LBNL binning or 1 for Linear binning
    NData->SetNSteps(NFits);
    NData->SetWeeks(Period);
    NData->SetNReactorPeriods(NReactorPeriods);
    NData->SetBCWModel(NL[0]);
    NData->SetLBNLModel(NL[1]);
    NData->SetUnifiedModel(NL[2]);
    
    NData->SetToyMC(ToyMC);
    
    NData->SetAllRandomSystematics(0);
    NData->SetStatisticalFluctuation(0);
    
    VData->SetAllRandomSystematics(1);
    VData->SetStatisticalFluctuation(1);
    
    VData->SetCombineMode(CombineMode); //0 is 9x9, 1 is 1x1 and 2 is 2x2
    //    Data->SetSin22t12(0);
    //Parameters of the model
    VData->SetAnalysis(Analysis);//  Gd or H data
    VData->SetBinning(Binning);//  0 for LBNL binning or 1 for Linear binning
    VData->SetNSteps(NFits);
    VData->SetWeeks(Period); //   1 or 20 so far.
    VData->SetNReactorPeriods(NReactorPeriods);
    VData->SetBCWModel(NL[0]);
    VData->SetLBNLModel(NL[1]);
    VData->SetUnifiedModel(NL[2]);
    
    VData->SetToyMC(ToyMC);
    
    Double_t DeltaWidth;
    Double_t SinWidth;
    
    Prediction* DataPred = new Prediction(NData);//Using data to produce the prediction
    Prediction* NomPred = new Prediction(NData);//Using reactor model to produce the prediction
    
    TH1D* DataHisto = new TH1D("DataHisto","DataHisto",n_evis_bins,evis_bins);
    TH1D* DataHisto1 = new TH1D("DataHisto1","DataHisto1",n_evis_bins,evis_bins);
    
    TH1D* NominalHisto = new TH1D("NomHisto","NomHisto",n_evis_bins,evis_bins);
    TH1D* NominalHisto1 = new TH1D("NomHisto1","NomHisto1",n_evis_bins,evis_bins);
    
    TH1D* VariationHisto[MaxSystematics];
    TH1D* VariationHisto1[MaxSystematics];
    
    for(Int_t SystematicI = 8; SystematicI<MaxSystematics; SystematicI++)
    {
        VariationHisto[SystematicI] = new TH1D(Form("VariationHisto_%d",SystematicI),Form("VariationHisto_%d",SystematicI),n_evis_bins,evis_bins);
        VariationHisto1[SystematicI] = new TH1D(Form("VariationHisto1_%d",SystematicI),Form("VariationHisto1_%d",SystematicI),n_evis_bins,evis_bins);
    }
    
    std::cout << " GENERATING NOMINAL AND FLUCTUATED TOY MC AND DATA SAMPLES IN A 101 x 101 GRID " << std::endl;
    
    if(!StatisticalFluctuation)
    {
        TTree *TNom= new TTree("TNom","ToyMCHistograms");
        
        TNom->Branch("sin22t13",&sin22t13,"sin22t13/D");
        TNom->Branch("dm2_ee",&dm2_ee,"dm2_ee/D");
        if(CombineMode==0)
        {
            std::cout << "This has not been coded, since it is not possible to fit 9x9 predictions" << std::endl;
            exit(EXIT_FAILURE);
        }
        else if(CombineMode==1)
        {
            TNom->Branch("NominalHisto","TH1D",&NominalHisto);
            TNom->Branch("DataHisto","TH1D",&DataHisto);
        }
        else if(CombineMode==2)
        {
            TNom->Branch("NominalHistoDayaBay","TH1D",&NominalHisto);
            TNom->Branch("DataHistoDayaBay","TH1D",&DataHisto);
            
            TNom->Branch("NominalHistoLingAo","TH1D",&NominalHisto1);
            TNom->Branch("DataHistoLingAo","TH1D",&DataHisto1);
        }
        
        DeltaWidth = (dm2_eeend-dm2_eestart)*1./(NFits-1);
        SinWidth = (s22t13end-s22t13start)*1./(NFits-1);
        
        for(Int_t step=0;step < NFits;step++)
        {
            sin22t13 = SinWidth*step + s22t13start;
            
            for(Int_t stepM = 0; stepM < NFits ; stepM++)
            {
                dm2_ee = DeltaWidth*stepM + dm2_eestart;
                
                std::cout << "DATA PREDICTION" << std::endl;
                DataPred->MakePrediction(sin22t13, dm2_ee, 0, Period-1, 0,0); //P12E Data
                
                std::cout << "NOMINAL PREDICTION" << std::endl;
                NomPred->MakePrediction(sin22t13, dm2_ee, 0, Period-1, 1,0); //ToyMC nominal
                
                for (Int_t far = 0; far<MaxFarCombine; far++)
                {
                    for (Int_t near = 0; near<MaxNearCombine; near++)
                    {
                        if(CombineMode==1)
                        {
                            NominalHisto = NomPred->GetPrediction(far,near);
                            DataHisto = DataPred->GetPrediction(far,near);
                            
                            NominalHisto->SetName(Form("Combined%d Nominal Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                            NominalHisto->SetTitle(Form("Combined%d Nominal Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                            
                            DataHisto->SetName(Form("Combined%d Data Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                            DataHisto->SetTitle(Form("Combined%d Data Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                        }
                        else if(CombineMode==2)
                        {
                            if(near==0)
                            {
                                NominalHisto = NomPred->GetPrediction(far,near);
                                DataHisto = DataPred->GetPrediction(far,near);
                                
                                NominalHisto->SetName(Form("Daya Bay Combined%d Nominal Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                NominalHisto->SetTitle(Form("Daya Bay Combined%d Nominal Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                
                                DataHisto->SetName(Form("Daya Bay Combined%d Data Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                DataHisto->SetTitle(Form("Daya Bay Combined%d Data Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                
                            }
                            else
                            {
                                NominalHisto1 = NomPred->GetPrediction(far,near);
                                DataHisto1 = DataPred->GetPrediction(far,near);
                                
                                NominalHisto1->SetName(Form("Ling Ao Combined%d Nominal Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                NominalHisto1->SetTitle(Form("Ling Ao Combined%d Nominal Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                
                                DataHisto1->SetName(Form("Ling Ao Combined%d Data Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                DataHisto1->SetTitle(Form("Ling Ao Combined%d Data Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                            }
                        }
                    }
                }
                
                TNom->Fill();
                
                for (Int_t far = 0; far<3; far++)
                {
                    for (Int_t near = 0; near<3; near++)
                    {
                        //                            VarPred->DeletePrediction(far, near);
                        NomPred->DeletePrediction(far, near);
                        DataPred->DeletePrediction(far, near);
                    }
                }
            }
        }
        
        TFile* TreeF1 = new TFile(Form("./ToyMCTrees/ToyMCTreeCombined%d.root",CombineMode),"recreate");
        TNom->Write();
        delete TreeF1;
        delete TNom;
    }
    else
    {
        if(CholeskyVariations)
        {
            //Produce Variation Tree using Cholesky decomposition in each of the covariance matrices produced with the Toy MC to generate random spectra. This idea is to accelerate the production of prediction with variations, if it worked the production would be much faster.
            
            // NOTE:
            
            // This doesn't seem to work yet, not all the individual systematic variations are invertible although the total one (SystematicI = 8) is invertible. What about the statistical term? should it be added to the matrices or not? If so, maybe all matrices would be invertible but it wouldn't be possible to isolate the systematic error from the statistical error, that's why I didn't take that approach.
            
            Int_t Experiment = 0;
            
            TH1D* NominalPrediction[MaxNearCombine*MaxFarCombine];
            //            for (Int_t far = 0; far<MaxFarCombine; far++)
            //            {
            //                for (Int_t near = 0; near<MaxNearCombine; near++)
            //                {
            //                    NominalPrediction[near+MaxNearCombine*far] = new TH1D(Form("NominalPred%d_%d",far,near),Form("NominalPred%d_%d",far,near),n_evis_bins,evis_bins);
            //                }
            //            }
            DeltaWidth = (dm2_eeend-dm2_eestart)*1./(NFits-1);
            SinWidth = (s22t13end-s22t13start)*1./(NFits-1);
            
            NomPred->LoadRootCovarianceMatrices(Period-1);//Load covariance matrices (without prediction included)
            
            TTree *TVar= new TTree("TVar","ToyMCHistograms");
            
            TVar->Branch("sin22t13",&sin22t13,"sin22t13/D");
            TVar->Branch("dm2_ee",&dm2_ee,"dm2_ee/D");
            TVar->Branch("Experiment",&Experiment,"Experiment/I");
            
            for(Int_t SystematicI = 8; SystematicI < MaxSystematics; SystematicI++)
            {
                if(CombineMode==0)
                {
                    std::cout << "This has not been coded, since it is not possible to fit 9x9 predictions" << std::endl;
                    exit(EXIT_FAILURE);
                }
                else if(CombineMode==1)
                {
                    TVar->Branch(Form("VariationHisto_%d",SystematicI),"TH1D",&VariationHisto[SystematicI]);
                }
                else if(CombineMode==2)
                {
                    TVar->Branch(Form("VariationHistoDayaBay_%d",SystematicI),"TH1D",&VariationHisto[SystematicI]);
                    TVar->Branch(Form("VariationHistoLingAo_%d",SystematicI),"TH1D",&VariationHisto1[SystematicI]);
                }
            }
            
            for(Int_t step=0;step < NFits;step++)
            {
                sin22t13 = SinWidth*step + s22t13start;
                
                for(Int_t stepM = 0; stepM < NFits ; stepM++)
                {
                    dm2_ee = DeltaWidth*stepM + dm2_eestart;
                    
                    std::cout << "VARIED PREDICTION" << std::endl;
                    NomPred->MakePrediction(sin22t13, dm2_ee, 0, Period-1, 1,0); //Nominal ToyMC with no variations
                    //Load normalize covariance matrices, multiply by nomprediction in each sin22t3 and dm2_ee point
                    
                    for (Int_t far = 0; far<MaxFarCombine; far++)
                    {
                        for (Int_t near = 0; near<MaxNearCombine; near++)
                        {
                            if(CombineMode==1)
                            {
                                NominalPrediction[near+MaxNearCombine*far] = NomPred->GetPrediction(far,near);
                                
                                NominalPrediction[near+MaxNearCombine*far]->SetName(Form("Combined%d Varied Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                NominalPrediction[near+MaxNearCombine*far]->SetTitle(Form("Combined%d Varied Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                            }
                            else if(CombineMode==2)
                            {
                                if(near==0)
                                {
                                    NominalPrediction[near+MaxNearCombine*far] = NomPred->GetPrediction(far,near);
                                    
                                    NominalPrediction[near+MaxNearCombine*far]->SetName(Form("Daya Bay Combined%d Varied Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                    NominalPrediction[near+MaxNearCombine*far]->SetTitle(Form("Daya Bay Combined%d Varied Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                }
                                else
                                {
                                    NominalPrediction[near+MaxNearCombine*far] = NomPred->GetPrediction(far,near);
                                    
                                    NominalPrediction[near+MaxNearCombine*far]->SetName(Form("Ling Ao Combined%d Varied Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                    NominalPrediction[near+MaxNearCombine*far]->SetTitle(Form("Ling Ao Combined%d Varied Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                }
                            }
                        }
                    }
                    // Apply Cholesky decomposition
                    NomPred->ProduceCovToyMCSample(Period-1,&NominalPrediction[0]);
                    
                    for(Int_t SystematicI = 8; SystematicI < MaxSystematics; SystematicI++)
                    {
                        //Get different Prediction = L*rand[x] as many times as experiments
                        for(Int_t experiment = 0; experiment < MaxExperiments ; experiment++)
                        {
                            TH1D* LongPrediction = NomPred->GetToyMCSample(SystematicI);//returns the whole vector (dimensions = 9*,2* or 1*n_evis_bins)
                            
                            for(Int_t i =0;i<MaxNearCombine*n_evis_bins;i++)
                            {
                                if(i<n_evis_bins)
                                {
                                    VariationHisto[SystematicI]->SetBinContent(i+1,LongPrediction->GetBinContent(i+1));
                                }
                                else
                                {
                                    VariationHisto1[SystematicI]->SetBinContent(i-n_evis_bins+1,LongPrediction->GetBinContent(i+1));
                                }
                            }
                            
                            TVar->Fill();
                            
                            Experiment++;
                        }
                    }
                    
                    for (Int_t far = 0; far<NADs/2; far++)
                    {
                        for (Int_t near = 0; near<NADs/2; near++)
                        {
                            NomPred->DeletePrediction(far, near);
                            //                            delete NominalPrediction[near+(NADs/2*far)];
                        }
                    }
                }
            }
            
            TFile* TreeF2 = new TFile(Form("./ToyMCTrees/VariationsToyMCTreeCombined%d.root",CombineMode),"recreate");
            TVar->Write();
            delete TreeF2;
            delete TVar;
        }
        else
        {
            Prediction* VarPred = new Prediction(VData);
            
            if(Fake_Experiments)
            {
                Int_t Experiment = 0;
                Int_t experiments = 100;

                std::cout << " GENERATING 100 FLUCTUATED TOY MC SAMPLES IN A 21 x 21 GRID" << std::endl << std::endl;
                
                TTree *TFake= new TTree("TFake","ToyMCHistogramsExperiments");
                
                for(Int_t SystematicI = 8; SystematicI < MaxSystematics; SystematicI++)
                {
                    if(CombineMode==0)
                    {
                        std::cout << "This has not been coded, since it is not possible to fit 9x9 predictions" << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    else if(CombineMode==1)
                    {
                        TFake->Branch("VariationHisto","TH1D",&VariationHisto[SystematicI]);
                    }
                    else if(CombineMode==2)
                    {
                        TFake->Branch("VariationHistoDayaBay","TH1D",&VariationHisto[SystematicI]);
                        TFake->Branch("VariationHistoLingAo","TH1D",&VariationHisto1[SystematicI]);
                    }
                }
                TFake->Branch("sin22t13",&sin22t13,"sin22t13/D");
                TFake->Branch("dm2_ee",&dm2_ee,"dm2_ee/D");
                TFake->Branch("Experiment",&Experiment,"Experiment/I");

                Int_t steps = 21;

                DeltaWidth = (dm2_eeend-dm2_eestart)*1./(steps-1);
                SinWidth = (s22t13end-s22t13start)*1./(steps-1);

                for(Int_t SystematicI = 8; SystematicI < MaxSystematics; SystematicI++)
                {
                    for(Int_t step=0;step < steps;step++)
                    {
                        sin22t13 = SinWidth*step + s22t13start;
                        
                        for(Int_t stepM = 0; stepM < steps ; stepM++)
                        {
                            dm2_ee = DeltaWidth*stepM + dm2_eestart;
                            
                            for(Int_t i =0;i<experiments;i++)
                            {
                                Experiment = i;

                                VarPred->MakePrediction(sin22t13, dm2_ee, 1, Period-1, 1,0); //ToyMC with variations
                                
                                for (Int_t far = 0; far<MaxFarCombine; far++)
                                {
                                    for (Int_t near = 0; near<MaxNearCombine; near++)
                                    {
                                        if(CombineMode==1)
                                        {
                                            VariationHisto[SystematicI] = VarPred->GetPrediction(far,near);
                                            VariationHisto[SystematicI]->SetName(Form("Combined%d Prediction+variations_experiment%d,sin_%f,DM_%f week%d",i,CombineMode,sin22t13,dm2_ee,Period));
                                            VariationHisto[SystematicI]->SetTitle(Form("Combined%d Prediction+variations_experiment%d,sin_%f,DM_%f week%d",i,CombineMode,sin22t13,dm2_ee,Period));
                                        }
                                        else if(CombineMode==2)
                                        {
                                            if(near==0)
                                            {
                                                VariationHisto[SystematicI] = VarPred->GetPrediction(far,near);
                                                VariationHisto[SystematicI]->SetName(Form("Daya Bay Combined%d Prediction+variations_experiment%d,sin_%f,DM_%f week%d",i,CombineMode,sin22t13,dm2_ee,Period));
                                                VariationHisto[SystematicI]->SetTitle(Form("Daya Bay Combined%d Prediction+variations_experiment%d,sin_%f,DM_%f week%d",i,CombineMode,sin22t13,dm2_ee,Period));
                                            }
                                            else
                                            {
                                                
                                                VariationHisto1[SystematicI] = VarPred->GetPrediction(far,near);
                                                VariationHisto1[SystematicI]->SetName(Form("Ling Ao Combined%d Prediction+variations_experiment%d,sin_%f,DM_%f week%d",i,CombineMode,sin22t13,dm2_ee,Period));
                                                VariationHisto1[SystematicI]->SetTitle(Form("Ling Ao Combined%d Prediction+variations_experiment%d,sin_%f,DM_%f week%d",i,CombineMode,sin22t13,dm2_ee,Period));
                                            }
                                        }
                                    }
                                    
                                    TFake->Fill();
                                    
                                }
                                for (Int_t far = 0; far<NADs/2; far++)
                                {
                                    for (Int_t near = 0; near<NADs/2; near++)
                                    {
                                        VarPred->DeletePrediction(far, near);
                                    }
                                }
                            }
                        }
                    }
                    
                    TFile* FTreeF = new TFile(Form("./ToyMCTrees/FakeExperiments%d_Combined_%d.root",SystematicI,CombineMode),"recreate");
                        TFake->Write();
                    delete FTreeF;
                    delete TFake;
                }
            }
            else
            {
                
                TTree *TVar= new TTree("TVar","ToyMCHistograms");
                
                TVar->Branch("sin22t13",&sin22t13,"sin22t13/D");
                TVar->Branch("dm2_ee",&dm2_ee,"dm2_ee/D");
                
                DeltaWidth = (dm2_eeend-dm2_eestart)*1./(NFits-1);
                SinWidth = (s22t13end-s22t13start)*1./(NFits-1);
                
                for(Int_t SystematicI = 8; SystematicI < MaxSystematics; SystematicI++)
                {
                    if(CombineMode==0)
                    {
                        std::cout << "This has not been coded, since it is not possible to fit 9x9 predictions" << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    else if(CombineMode==1)
                    {
                        TVar->Branch("VariationHisto","TH1D",&VariationHisto[SystematicI]);
                    }
                    else if(CombineMode==2)
                    {
                        TVar->Branch("VariationHistoDayaBay","TH1D",&VariationHisto[SystematicI]);
                        TVar->Branch("VariationHistoLingAo","TH1D",&VariationHisto1[SystematicI]);
                    }
                    
                    for(Int_t step=0;step < NFits;step++)
                    {
                        sin22t13 = SinWidth*step + s22t13start;
                        
                        for(Int_t stepM = 0; stepM < NFits ; stepM++)
                        {
                            dm2_ee = DeltaWidth*stepM + dm2_eestart;
                            
                            std::cout << "VARIED PREDICTION" << std::endl;
                            VarPred->MakePrediction(sin22t13, dm2_ee, 1, Period-1, 1,0); //ToyMC with variations
                            
                            for (Int_t far = 0; far<MaxFarCombine; far++)
                            {
                                for (Int_t near = 0; near<MaxNearCombine; near++)
                                {
                                    if(CombineMode==1)
                                    {
                                        VariationHisto[SystematicI] = VarPred->GetPrediction(far,near);
                                        
                                        VariationHisto[SystematicI]->SetName(Form("Combined%d Varied Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                        VariationHisto[SystematicI]->SetTitle(Form("Combined%d Varied Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                    }
                                    else if(CombineMode==2)
                                    {
                                        if(near==0)
                                        {
                                            VariationHisto[SystematicI] = VarPred->GetPrediction(far,near);
                                            
                                            VariationHisto[SystematicI]->SetName(Form("Daya Bay Combined%d Varied Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                            VariationHisto[SystematicI]->SetTitle(Form("Daya Bay Combined%d Varied Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                        }
                                        else
                                        {
                                            VariationHisto1[SystematicI] = VarPred->GetPrediction(far,near);
                                            
                                            VariationHisto1[SystematicI]->SetName(Form("Ling Ao Combined%d Varied Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                            VariationHisto1[SystematicI]->SetTitle(Form("Ling Ao Combined%d Varied Prediction,sin_%f,DM_%f week%d",CombineMode,sin22t13,dm2_ee,Period));
                                        }
                                    }
                                }
                            }
                            
                            TVar->Fill();
                            
                            for (Int_t far = 0; far<NADs/2; far++)
                            {
                                for (Int_t near = 0; near<NADs/2; near++)
                                {
                                    VarPred->DeletePrediction(far, near);
                                }
                            }
                        }
                    }
                    
                    TFile* TreeF2 = new TFile(Form("./ToyMCTrees/Variations%d_ToyMCTreeCombined%d.root",SystematicI,CombineMode),"recreate");
                    TVar->Write();
                    delete TreeF2;
                    delete TVar;
                }
            }
            
            delete VarPred;
        }
    }

    delete DataPred;
    delete NomPred;
    
    cout.rdbuf(coutstream);
    gBenchmark->Show("GenerateTree");
}

void FitterGui :: RunFitterTests()
{
    gBenchmark->Start("TestFitter");
    std::cout << "Running Fitter TESTs" << std::endl;

    if(!PrintOnConsole)
    {
        coutstream = cout.rdbuf(0);//change stream of cout
    }
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
    
    Int_t DataSet=2;
    NominalData* Data = new NominalData(Analysis,DataSet);
    
    Data->SetToyMCSamplesDirectory(ToyMCSampleDirectory);
    Data->SetPredictionDirectory(NominalPredictionsDirectory);
    Data->SetResponseDirectory(ResponseMatrixDirectory);
    Data->SetBkgCovDirectory(BkgCovarianceMatrixDirectory);
    Data->SetSysCovDirectory(SysCovarianceMatrixDirectory);
    
    //  No variations in the ToyMC Prediction, if desired set to 1:
    Data->SetAllRandomSystematics(0);
    Data->SetStatisticalFluctuation(StatisticalFluctuation);
    
    Data->SetCombineMode(CombineMode); //0 is 9x9, 1 is 1x1 and 2 is 2x2
    Data->SetUseToyMCTree(UseToyMCTree);
    //    Data->SetSin22t12(0);
    //    Data->SetSin22t13(0);
    
    //Parameters of the model
    Data->SetAnalysis(Analysis);//  Gd or H data
    Data->SetBinning(Binning);//  0 for LBNL binning or 1 for Linear binning
    Data->SetNSteps(NFits);// 101 in the final version.
    Data->SetWeeks(Period);
    Data->SetNReactorPeriods(NReactorPeriods);
    Data->SetBCWModel(NL[0]);
    Data->SetLBNLModel(NL[1]);
    Data->SetUnifiedModel(NL[2]);
    
    Data->SetToyMC(1);//Always TOYMC fitter tests.
    
    Fitter* Fit = new Fitter(Data);
    
    Int_t GridSamples = 21;
    //Use it with Combine = 2, stats on, and use toymctree.
    
    Fit->TestRandomExperimentsChi2(Minuit,Fit2D,FitSin22t13, this , MaxExperiments, GridSamples);
    
    delete Fit;
    
    cout.rdbuf(coutstream);
    gBenchmark->Show("TestFitter");
    
}

