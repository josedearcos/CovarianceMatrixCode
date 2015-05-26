#pragma once
#include <TQObject.h>
#include <RQ_OBJECT.h>
#include "Riostream.h"

#ifndef ROOT_TGDockableFrame
#include "TGDockableFrame.h"
#endif
#ifndef ROOT_TGMenu
#include "TGMenu.h"
#endif
#ifndef ROOT_TGMdiDecorFrame
#include "TGMdiDecorFrame.h"
#endif
#ifndef ROOT_TG3DLine
#include "TG3DLine.h"
#endif
#ifndef ROOT_TGMdiFrame
#include "TGMdiFrame.h"
#endif
#ifndef ROOT_TGMdiMainFrame
#include "TGMdiMainFrame.h"
#endif
#ifndef ROOT_TGMdiMenu
#include "TGMdiMenu.h"
#endif
#ifndef ROOT_TGColorDialog
#include "TGColorDialog.h"
#endif
#ifndef ROOT_TGListBox
#include "TGListBox.h"
#endif
#ifndef ROOT_TGNumberEntry
#include "TGNumberEntry.h"
#endif
#ifndef ROOT_TGScrollBar
#include "TGScrollBar.h"
#endif
#ifndef ROOT_TGComboBox
#include "TGComboBox.h"
#endif
#ifndef ROOT_TGuiBldHintsEditor
#include "TGuiBldHintsEditor.h"
#endif
#ifndef ROOT_TGuiBldNameFrame
#include "TGuiBldNameFrame.h"
#endif
#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif
#ifndef ROOT_TGFileDialog
#include "TGFileDialog.h"
#endif
#ifndef ROOT_TGShutter
#include "TGShutter.h"
#endif
#ifndef ROOT_TGButtonGroup
#include "TGButtonGroup.h"
#endif
#ifndef ROOT_TGCanvas
#include "TGCanvas.h"
#endif
#ifndef ROOT_TGFSContainer
#include "TGFSContainer.h"
#endif
#ifndef ROOT_TGFontDialog
#include "TGFontDialog.h"
#endif
#ifndef ROOT_TGuiBldEditor
#include "TGuiBldEditor.h"
#endif
#ifndef ROOT_TGColorSelect
#include "TGColorSelect.h"
#endif
#ifndef ROOT_TGButton
#include "TGButton.h"
#endif
#ifndef ROOT_TGFSComboBox
#include "TGFSComboBox.h"
#endif
#ifndef ROOT_TGLabel
#include "TGLabel.h"
#endif
#ifndef ROOT_TGMsgBox
#include "TGMsgBox.h"
#endif
#ifndef ROOT_TRootGuiBuilder
#include "TRootGuiBuilder.h"
#endif
#ifndef ROOT_TGTab
#include "TGTab.h"
#endif
#ifndef ROOT_TGListView
#include "TGListView.h"
#endif
#ifndef ROOT_TGSplitter
#include "TGSplitter.h"
#endif
#ifndef ROOT_TGStatusBar
#include "TGStatusBar.h"
#endif
#ifndef ROOT_TGListTree
#include "TGListTree.h"
#endif
#ifndef ROOT_TGuiBldGeometryFrame
#include "TGuiBldGeometryFrame.h"
#endif
#ifndef ROOT_TGToolTip
#include "TGToolTip.h"
#endif
#ifndef ROOT_TGToolBar
#include "TGToolBar.h"
#endif
#ifndef ROOT_TRootEmbeddedCanvas
#include "TRootEmbeddedCanvas.h"
#endif
#ifndef ROOT_TCanvas
#include "TCanvas.h"
#endif
#ifndef ROOT_TGuiBldDragManager
#include "TGuiBldDragManager.h"
#endif

#include "TGProgressBar.h"
#include "TGSlider.h"

const Int_t MaximumSamples = 10000;

class TGWindow;
class TGMainFrame;
class TRootEmbeddedCanvas;

class FitterGui
{
    RQ_OBJECT("FitterGui")
private:
    
    TObject* ChiH;
    
    string AnalysisString;
    string NominalString;
    string CorrelationString;
    string CovString;
    string VarString;
    TGMainFrame *fMain;
    
    TRootEmbeddedCanvas *fMaincanvas;
    TRootEmbeddedCanvas *fBackgroundcanvas;
    TRootEmbeddedCanvas *fADVisCanvas;
    TRootEmbeddedCanvas *fADTrueCanvas;
    TRootEmbeddedCanvas *fFractionCanvas;
    TRootEmbeddedCanvas *fResponseCanvas;
    TRootEmbeddedCanvas *fCovCanvas;
    TRootEmbeddedCanvas *fChiCanvas;
    TRootEmbeddedCanvas *fVariationsCanvas;
    TRootEmbeddedCanvas *fBudgetCanvas;
    TGHSlider* hslider;
    TGHSlider* Sslider;
    TGHSlider* Dslider;
    
    TCanvas* fCanvas;
    TCanvas* VCanvas;
    TCanvas* ADVisCanvas;
    TCanvas* ADTrueCanvas;
    TCanvas* ChiCanvas;
    TCanvas* FractionCanvas;
    TCanvas* ResponseCanvas;
    TCanvas* CovCanvas;
    TCanvas* BudgetCanvas;

    TGHProgressBar *ToyMCProgressBar;
    TGHProgressBar *FitterProgressBar;
    TGCompositeFrame *fToyMCFrame2;
    TGCompositeFrame *fLoadInputs2;
    TGCompositeFrame* fFitterFrame2;
    
    Int_t visible_bins;
    Int_t true_bins;
    Int_t Period;
    Int_t NReactorPeriods;
    Int_t DataPeriods;
    Int_t NADs;
    Int_t NSamples;
    Int_t CombineMode;
    Int_t NFits;
    Int_t NCovarianceMatrices;
    bool Analysis;
    bool flagCombine;
    bool flagNear;
    bool flagFar;
    bool Binning;
    Int_t PlotBin;
    bool Automatic;
    bool ToyMC;
    bool Minuit;
    bool deleteFlag;
    bool deleteFlagSpec;
    Int_t MaxFar;
    Int_t MaxNear;
    Int_t NearTrueIndex;
    
    Double_t Sin22t13SliderValue;
    Double_t DeltaMSilderValue;
    
    bool UseToyMCTree;

    bool Fit2D;
    bool FitSin22t13;
    //Manual Control:
    //To calculate Covariance Matrices set to 1. Activate only one at a time.
    //Backgrounds
    bool VaryAccidentalMatrix;
    bool VaryLiHeMatrix;
    bool VaryFastNeutronsMatrix;
    bool VaryAmCMatrix;
    bool DistortLiHeMatrix;
    bool DistortFastNeutronsMatrix;
    bool DistortAmCMatrix;
    //Systematics
    bool IsotopeMatrix;
    bool ReactorPowerMatrix;
    bool EnergyScaleMatrix;
    //        bool EnergyOffsetMatrix;
    //        bool AbsoluteScaleMatrix;
    //        bool AbsoluteOffsetMatrix;
    bool IAVMatrix;
    bool NLMatrix;
    bool ResolutionMatrix;
    bool Sin22t12Matrix;
    bool EfficiencyMatrix;
    
    bool PlotCovariance;
    string delSpaces(string &str);

    TGNumberEntry* NSamplesBox;
    TGNumberEntry* PeriodBox;
    TGNumberEntry* NReactorPeriodBox;
    TGNumberEntry* CombineBox;
    TGNumberEntry* NFitsBox;
    TGButtonGroup* VariationsBox;
    TGButtonGroup* CorrelationBox;
    TGButtonGroup* NLBox;
    TGButtonGroup* NominalBox;
    TGButtonGroup* IndexBox;
    TGButtonGroup* FitterBox;
    TGButtonGroup* Fitter1DBox;
    TGRadioButton* fR[3];
    TGRadioButton* fN[2];
    TGRadioButton* fI[3];
    TGRadioButton* fF[4];
    TGRadioButton* fV[2];
    TGRadioButton* fC[2];
    TGRadioButton* f1DF[2];
    
    TGComboBox* HGdBox;
    TGComboBox* BinningBox;
    TGComboBox* FluctuationBox;
    TGComboBox* ToyMCTreeBox;
    TGComboBox* NADsBox;
    TGComboBox* ToyMCBox;
    TGComboBox* MinuitBox;
    TGComboBox* ResponseBox;
    TGComboBox* CovarianceMatrixBox;
    TGComboBox* CovariancePlotBox;
    TGComboBox* VariationsPlotBox;
    bool NL[3];
    bool FitterMode[2];
    bool Fitter1DMode[2];
    
    bool AutomaticBudget; //  To loop the fitter over all possible systematics
    bool TurnOnBudget; // Generates Turn On Error Budget
    bool TurnOffBudget;// Generates Turn Off Error Budget
    bool StatisticalFluctuation;
    
    string ResponseMatrixDirectory;
    string ToyMCSampleDirectory;
    string NominalPredictionsDirectory;
    string SysCovarianceMatrixDirectory;
    string BkgCovarianceMatrixDirectory;
public:
    FitterGui(const TGWindow *p,UInt_t w,UInt_t h);
    virtual ~FitterGui();
    void DoToyMC();
    void ChooseToyMC();
    void ChooseMinuit();
    void DoReadInputs();
    void DoFitter();
    void RunFitter();
    void RunToyMC();
    void RunFlux();
    void RunFitterTests();
    void LoadResponseMatrix();//To convert data to true energy
    void LoadNominalPredictions();//To fit and to create covariance matrices
    void LoadToyMCSamples();//To produce covariance matrices
    void Update(Int_t);
    void UpdateFitter(Int_t);
    void DoSamples();
    void DoCombine();
    void DoBinning();
    void CalculateBinning();
    
    void DoNLModel();
    void DoFitterMode();
    void Do1DFitterMode();
    void DoAnalysis();
    void DoNADs();
    void DoPeriod();
    void DoNReactorPeriods();
    void PlotBkgd();
    
    void PlotADVis();
    void PlotADTrue();
    
    void PlotSpectrumFraction();
    void PlotFar();
    void PlotNear();
    void PlotFarADTrue();
    void PlotNearADTrue();
    void PlotVariations();
    void PlotCov();
    void ChooseCovarianceMatrix();
    void ChoosePlotCov();
    void ChoosePlotCovariance();
    void ChoosePlotVariations();
    void ChoosePlotRatioVariations();
    void ChoosePlotSpectrumVariations();
    void PlotResponseMatrix();
    void ChooseResponseMatrix();
    void ChooseVariations();
    void PlotAccidental();
    void PlotLiHe();
    void PlotFN();
    void PlotAmC();

    void DoNear();
    void PlotCombine();
    void PlotAllADVis();
    void DoNominal();
    void DoChangeBin();
    void PlotErrorBudget();
    void PlotTurnOnBudget();
    void PlotTurnOffBudget();
    void RunTurnOnBudget();
    void RunTurnOffBudget();

    void RunResponseMatrix();
    void DoFits();

    void PlotChi();
    void GenerateToyMC();
    void DoChangeSinSquareOscillationParameter();
    void DoChangeDeltaSquareOscillationParameter();
    
    void DoStatisticalFluctuation();
    void DoUseToyMCSamples();
};

