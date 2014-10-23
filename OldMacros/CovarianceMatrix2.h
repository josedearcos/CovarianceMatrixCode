#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TMath.h"
#include <vector>
#include "Prediction.h"
#include "ReactorSpectrumMultiple.h"
#include "NominalData.h"
#include <math.h>

const Int_t NSamples = 1;//Number of samples used in the generation of covariance matrices

//const bool PreRebin = 1; //0 to rebin the Energy Matrix after its calculation, 1 to rebin the IAV matrix before the Energy Matrix calculation.

const Int_t MaxPeriods = 1;
const Int_t MaxDetectors = 8;
const Int_t MaxNearDetectors =4;
const Int_t MaxFarDetectors =4;
const Int_t MaxNbins=51;
const Int_t Halls=3;
const bool AddBackgrounds=0;//To add backgrounds
const bool VariateRate = 1; // To vary backgrounds rate
const bool DistortBackgrounds=1;//To vary backgrounds shape
const bool WriteOutput=0;//To save the covariance matrices in a .txt file.
const Int_t MatrixBins=240;//Same than IAV matrix for Gd. When the IAV H matrix is produced this can be selected through NominalData.h
const Int_t n_bcw_positron_nl = 1000;

const Double_t Me = 0.510999; // MeV
const Double_t Mn = 939.565; // MeV
const Double_t Mp = 938.272; // MeV
const Int_t n_unified_nl_points = 500;

class CovarianceMatrix2
{
private:
    TRandom3* rand;
    Prediction Pred;

    enum Systematic{ReactorE, EnergyE, IAVE, NLE, ResolutionE};
    Systematic SystematicE;
    enum Background{AccidentalE, LiHeE, FastNeutronsE, AmCE};
    Background BackgroundE;
    enum NLModel{BCWE, LBNLE, IHEPE};//Add here new NL models
    NLModel NLModelE;
    
    //AD configuration parameters:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    Int_t hall;
    
    //Binning parameters:
    Int_t Nbins;
    Int_t InitialBins;
    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    Int_t Nweeks;
    Int_t RebinFactor;
    Double_t BinWidth;
    Int_t TotalBins;
    char* OutputFileName;
    
    //Energy Matrix:
    Double_t PositronEnergy;
    Double_t MaxPositronEnergy;
    Int_t MaxPositronEnergyIndex;
    Int_t PositronEnergyIndex;
    Double_t Norm[MaxDetectors][MatrixBins];
    TH1F* PositronTrueSpectrumH[MaxDetectors][MatrixBins];
    TH1F* PositronIAVSpectrumH[MaxDetectors][MatrixBins];
    TH1F* PositronNLSpectrumH[MaxDetectors][MatrixBins];
    TH1F* PositronVisibleSpectrumH[MaxDetectors][MatrixBins];
    TH1F* EnergySlice[MaxDetectors][MatrixBins];
    
    //IAV:
    Double_t IAVError; // relative uncertainty of the IAV thickness
    Double_t IAVError[MaxDetectors]; // relative uncertainty of the IAV thickness
    Double_t IAVMatrix[MatrixBins][MatrixBins];
    
    //Non linearity:
    // Detector response non-linearlity function
    TF1 * nl_func;
    
    //Interpolation vectors
    vector<Int_t> sign;
    vector<Int_t> IAVEnergyIdx;
    vector<Double_t> IAVEnergy;
    vector<Double_t>  dEtrue;
    vector<Double_t>  binScaling;
    vector<Double_t>  dNdE;
    
    //BCW Model
    Double_t m_bcw_elec_nl_par[5];
    Double_t m_bcw_elec_nl_par_nominal[5];
    Double_t m_bcw_elec_nl_par_error[5];
    
    Double_t m_bcw_positron_nl_e[n_bcw_positron_nl];
    Double_t m_bcw_positron_nl_fac[n_bcw_positron_nl];
    
    Double_t m_rel_escale[MaxDetectors];
    Double_t m_rel_escale_nominal[MaxDetectors];
    Double_t m_rel_escale_error[MaxDetectors];
    Double_t m_rel_eoffset[MaxDetectors];
    TGraph* g_bcw_elec_nl_error[2];
    
    //LBNL Model
    Double_t total_err;
    Double_t m_lbnl_nl_par[3];
    Double_t m_lbnl_nl_par_nominal[3];
    Double_t m_lbnl_nl_par_error[3];
    Double_t m_lbnl_positron_nl_e[300];
    Double_t m_lbnl_positron_nl_fac[300];
    Double_t m_lbnl_positron_nl_err[3][300];
    
    //IHEP Model
    Double_t m_ihep_nl_par[5];
    Double_t m_ihep_nl_par_nominal[5];
    Double_t m_ihep_nl_par_error[5];
    
    //Unified Model
    Int_t m_num_unified_nl_pars;
    
    TGraph* g_unified_positron_nl;
    TGraph* g_unified_positron_nl_pulls[10];
    
    Double_t m_unified_positron_nl_e[n_unified_nl_points];
    Double_t m_unified_positron_nl_fac[n_unified_nl_points];
    Double_t m_unified_positron_nl_err[10][n_unified_nl_points];
    
    Double_t m_unified_nl_par_covmatrix[4][4];
    Double_t m_unified_nl_par_covmatrix_l[4][4];
    
    //Resolution
    // Detector resolution function
    TF1 * reso_func;
    
    Double_t resolutionRange; //Range of resolution
    Double_t m_detectorResolution; // Detector energy resolution parameter
    Double_t m_detectorResolution_nominal; // Detector energy resolution parameter
    Double_t m_detectorResolution_error; // Detector energy resolution parameter
    Double_t m_detectorResolution_error_uncorr; // Detector energy resolution parameter
    Double_t m_detectorResolution_bias[MaxDetectors]; // Random biases of the detector resolution.
    
    Double_t m_abs_escale;
    Double_t m_abs_escale_nominal;
    Double_t m_abs_escale_error;
    
    Double_t m_abs_eoffset;
    Double_t m_abs_eoffset_error;
    
    Double_t m_rel_eoffset_error;    //Relative errors:
    Double_t HAccidentalError[MaxDetectors];
    Double_t HLiHeError[Halls];
    Double_t HFastNeutronsError[Halls];
    Double_t HAmCError[Halls];
    Double_t DistortAcc;
    Double_t DistortLiHe;
    Double_t DistortFN;
    Double_t DistortAmC;

    Double_t ScaleFactorAccidental[MaxDetectors][MaxPeriods];
    Double_t ScaleFactorLiHe[MaxDetectors][MaxPeriods];
    Double_t ScaleFactorFastNeutrons[MaxDetectors][MaxPeriods];
    Double_t ScaleFactorAmC[MaxDetectors][MaxPeriods];

    //Choose matrix
    bool AccidentalMatrix;
    bool LiHeMatrix;
    bool FastNeutronsMatrix;
    bool AmCMatrix;
    
    char* OutputFileName;
    
    TH1F* PredictionH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    TH1F* RandomPredictionH[MaxNearDetectors][MaxNearDetectors][MaxPeriods];

//    TH1F* RandomPredictionH1[MaxNearDetectors][MaxNearDetectors][MaxPeriods];
  
    TH1F* NearHallSpectrumH[MaxNearDetectors][MaxPeriods];
    TH1F* RandomNearHallSpectrumH[MaxNearDetectors][MaxPeriods];
    
    TH1F* BackgroundSpectrumH[MaxDetectors][MaxPeriods];
    TH1F* AccidentalsH[MaxDetectors][MaxPeriods];
    TH1F* LiHeH[MaxDetectors][MaxPeriods];
    TH1F* FastNeutronsH[MaxDetectors][MaxPeriods];
    TH1F* AmCH[MaxDetectors][MaxPeriods];
    TH1F* RandomAccidentalsH[MaxDetectors][MaxPeriods];
    TH1F* RandomLiHeH[MaxDetectors][MaxPeriods];
    TH1F* RandomFastNeutronsH[MaxDetectors][MaxPeriods];
    TH1F* RandomAmCH[MaxDetectors][MaxPeriods];
    
    TH2F* TrueEnergyToVisibleEnergy;
    TH2F* EnergyMatrixH[MaxDetectors];
    TH2F* CovMatrix2H;
    TH2F* StatCov2H;
    TH2F* BkgCov2H;
//    TH1F* TestGaussianH;

    Double_t Sigma_Near[MaxFarDetectors][MaxNearDetectors][MaxPeriods][MaxNbins];
    Double_t Sigma_Far[MaxFarDetectors][MaxNearDetectors][MaxPeriods][MaxNbins];
    Double_t CovStat[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t NormCovStat[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t CovBkg[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t NormCovBkg[9*MaxNbins][9*MaxNbins][MaxPeriods];

    //Main Data loaded from Theta13-inputs_20week.txt if 20 weeks is used. Otherwise change input file.
    Int_t ObsEvts[MaxDetectors][MaxPeriods];
    Double_t Livetime[MaxDetectors][MaxPeriods];
    Double_t MuonVetoEff[MaxDetectors][MaxPeriods];
    Double_t DMCEff[MaxDetectors][MaxPeriods];
    Double_t TargetMass[MaxDetectors][MaxPeriods];
    Double_t BgEvts[MaxDetectors][MaxPeriods];
    Double_t AccEvts[MaxDetectors][MaxPeriods];
    Double_t AccErr[MaxDetectors][MaxPeriods];
    Double_t Li9Evts[MaxDetectors][MaxPeriods];
    Double_t Li9Err[MaxDetectors][MaxPeriods];
    Double_t FnEvts[MaxDetectors][MaxPeriods];
    Double_t FnErr[MaxDetectors][MaxPeriods];
    Double_t AmcEvts[MaxDetectors][MaxPeriods];
    Double_t AmcErr[MaxDetectors][MaxPeriods];
    Double_t AlnEvts[MaxDetectors][MaxPeriods];
    Double_t AlnErr[MaxDetectors][MaxPeriods];
    Double_t ErrEvts[MaxDetectors][MaxPeriods];
    
    //Load functions
    void LoadPredictions();
    void LoadBackgrounds();
    void LoadIavCorrection();
    void LoadNearHall();

    //Functions to vary background shapes. So far the same ones than LBNL.
    void FluctuateBackgrounds();
    TF1* GetDistortionFunction(Double_t amount);
    TF1* GetFastNeutronsDistortionFunction(Double_t amount);
    
    //Functions to calculate the energy matrix
    TF1* GetNeutrinoToVisibleFunction();
    void GetEnergyShift(Int_t TrueEnergyIndex);
    void GetIAVShift(Int_t TrueEnergyIndex);
    void GetNLShift(Int_t TrueEnergyIndex);
    void GetResolutionShift(Int_t TrueEnergyIndex);
    void CreateEnergyMatrix();
    void NormalizeEnergyMatrix(Int_t TrueEnergyIndex);
    void NLInterpolation();

    //Functions to generate the respective covariance matrices
    void GenerateStatisticalCovarianceMatrix();
    void GenerateCovarianceMatrix();
    
    Double_t VisibleEnergy0F(Double_t * x, Double_t * par);// Zeroth order true to visible function.
    Double_t VisibleEnergy1F(Double_t * x, Double_t * par);// First order true to visible function.
    
    Double_t nl_func(Double_t * x, Double_t * par);
    Double_t nl_func_bcw(Double_t* x, Double_t* par);
    Double_t nl_func_lbnl(Double_t* x, Double_t* par);
    
    Double_t reso_func_bcw(Double_t * x, Double_t * par);
    
public:
    CovarianceMatrix2();
    CovarianceMatrix2(NominalData Data);

    void CovarianceMatrixMain();

    void SetAccidentalMatrix(bool AccMatrix);
    void SetLiHeMatrix(bool LiMatrix);
    void SetFastNeutronsMatrix(bool FNmatrix);
    void SetAmCMatrix(bool AmCmatrix);
    
    void SetIsotopeMatrix(bool IsotopeMatrix);
    void SetReactorPowerMatrix(bool ReactorPowerMatrix);
    void SetIAVMatrix(bool IAVMatrix);
    void SetNLMatrix(bool NLMatrix);
    void SetResolutionMatrix(bool ResolutionMatrix);
    
    void RandomEnergyScaleMatrix();
    void RandomIAVMatrix();
    void RandomNLMatrix();
    void RandomResolutionMatrix();
    
    void SetBCWModel(bool BCW);
    void SetLBNLModel(bool LBNL);
    void SetIHEPModel(bool IHEP);
    
    void LoadNLParameters();

};
CovarianceMatrix2 :: CovarianceMatrix2()
{
    Nom = new NominalData();
    rand = new TRandom3();
    Pred = new Prediction();

    Nbins = Nom->GetNbins();
    InitialEnergy = Nom->GetEmin();
    FinalEnergy = Nom->GetEmax();
    InitialBins = Nom->GetInitialBins();
    RebinFactor = MatrixBins/(Nbins + InitialBins);
    MaxPositronEnergy=0;
    
    TotalBins = MatrixBins;
    BinWidth=(FinalVisibleEnergy-InitialVisibleEnergy)/TotalBins;

//    if(PreRebin)
//    {
//        TotalBins = Nbins;
//        BinWidth=(FinalEnergy-InitialEnergy)/TotalBins;
//    }
    
    Nweeks = Nom->GetWeeks();
    
    NADs = Nom->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    if(NADs == 8) //ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
    
    for(Int_t i=0;i<NADs;i++)
    {
        HAccidentalError[i]=Nom->GetHAccidentalError(i);
    }
    
    for(Int_t i=0;i<3;i++)
    {
        HLiHeError[i]=Nom->GetHLiHeError(i);
    }
    
    for(Int_t i=0;i<3;i++)
    {
        HFastNeutronsError[i]=Nom->GetHFNError(i);
    }
    
    for(Int_t i=0;i<3;i++)
    {
        HAmCError[i]=Nom->GetHAmCError(i);
    }
    
    IAVError=Nom->GetIAVError();
    
    m_abs_escale =1;
    m_abs_escale_nominal =1;
    m_abs_escale_error = 0.01;
    
    m_abs_eoffset = 0.0;
    m_abs_eoffset_error = 0.08;  //(MeV)
    
    m_rel_eoffset_error = 0.013; //(MeV)
    
    for(int idet=0; idet<3; idet++)
    {
        m_rel_escale[idet] = 1.0;
        m_rel_escale_error[idet] = 0.0035; // 0.35%
        m_rel_escale_nominal[idet] = m_rel_escale[idet];
        m_rel_eoffset[idet] = 0.0;
    }
       
    BackgroundE=4;//A number greater than the amount of different backgrounds included in the model (Right now Acc, LiHe, FN and AmC are indexes 0,1,2,3 respectively, if C(α,n) is added use 5);
    DistortAcc = 0;
    DistortLiHe= 0;
    DistortFN  = 0;
    DistortAmC = 0;
}

CovarianceMatrix2 :: CovarianceMatrix2(NominalData Data)
{
    rand = new TRandom3();
    Pred = new Prediction(Data);
    
    Nbins = Data.GetNbins();
    InitialEnergy = Data.GetEmin();
    FinalEnergy = Data.GetEmax();
    InitialVisibleEnergy =Data.GetEVisMin();
    FinalVisibleEnergy = Data.GetEVisMax();
    InitialBins = Data.GetInitialBins();

    RebinFactor = MatrixBins/(Nbins + InitialBins);
    MaxPositronEnergy=0;
    
    TotalBins = MatrixBins;
    BinWidth=(FinalVisibleEnergy-InitialVisibleEnergy)/TotalBins;
    
//    if(PreRebin)
//    {
//        TotalBins = Nbins+InitialBins;
//        BinWidth=(FinalEnergy-InitialEnergy)/TotalBins;
//    }
//    
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
    
    for(Int_t i=0;i<NADs;i++)
    {
        HAccidentalError[i]=Data.GetHAccidentalError(i);
    }
    
    for(Int_t i=0;i<3;i++)
    {
        HLiHeError[i]=Data.GetHLiHeError(i);
    }
    
    for(Int_t i=0;i<3;i++)
    {
        HFastNeutronsError[i]=Data.GetHFNError(i);
    }

    for(Int_t i=0;i<3;i++)
    {
        HAmCError[i]=Data.GetHAmCError(i);
    }
    
    IAVError=Data.GetIAVError();
  
    m_abs_escale =1;
    m_abs_escale_nominal =1;
    m_abs_escale_error = 0.01;
    
    m_abs_eoffset = 0.0;
    m_abs_eoffset_error = 0.08;  //(MeV)
    
    m_rel_eoffset_error = 0.013; //  (MeV)
    
    for(int idet=0; idet<3; idet++)
    {
        m_rel_escale[idet] = 1.0;
        m_rel_escale_error[idet] = 0.0035; // 0.35%
        m_rel_escale_nominal[idet] = m_rel_escale[idet];
        m_rel_eoffset[idet] = 0.0;
    }
    
    BackgroundE=4;//A number greater than the amount of different backgrounds included in the model (Right now Acc, LiHe, FN and AmC are indexes 0,1,2,3 respectively, if C(α,n) is added use 5);
    DistortAcc = 0;
    DistortLiHe= 0;
    DistortFN  = 0;
    DistortAmC = 0;
}

void CovarianceMatrix2 :: CovarianceMatrixMain()
{
    LoadBackgrounds();
    LoadPredictions();
    LoadNearHall();
    LoadIavCorrection();
    GenerateStatisticalCovarianceMatrix();

    CovMatrix2H=(TH2F*)StatCov2H->Clone("Covariance Matrix");
    CovMatrix2H->Reset();
    
//    TestGaussianH = new TH1F("Test Gaussian","Test Gaussian",100,-50,50);
//    TestGaussianH->Reset();
    cout << "Randomize Systematic #" << SystematicE << endl;

    CovarianceMatrix2* NLObject = new CovarianceMatrix2();//Finally the standard constructor is useful
    switch (NLModelE)
    {
        case 0://BCW NL Model
            NLObject->SetBCWModel(1);
            break;
        case 1://LBNL NL Model
            NLObject->SetLBNLModel(1);
            break;
        case 2://IHEP NL Model
            NLObject->SetIHEPModel(1);
            break;
        default:
    }
    NLObject->LoadNLParameters();
    for (Int_t samples = 0; samples<NSamples; samples++)
    {
        switch (SystematicE)
        {
            case 0://Vary Reactor
                Pred.MakePrediction();
                break;
            case 1://Vary Energy Scale
                RandomEnergyScaleMatrix();
                break;
            case 2://Vary IAV
                RandomIAVMatrix();
                break;
            case 3://Vary NL
                NLObject->LoadNLParameters();//Reset
                RandomNLMatrix();
                break;
            case 4://Vary Resolution
                RandomResolutionMatrix();
                break;
            default://Add nominal systematics
        }
        
        switch (NLModelE)
        {
            case 0://BCW NL Model
                cout << "Using BCW NL Model"<< endl;
                nl_func = new TF1("nl_func",NLObject,&CovarianceMatrix2::nl_func_bcw,InitialVisibleEnergy,FinalVisibleEnergy,7,"CovarianceMatrix2","nl_func_bcw");
                for (Int_t idet = 0; idet < ADsEH3; idet++)
                {
                    //NL function set up
                    nl_func->SetParameters(m_bcw_elec_nl_par[0],m_bcw_elec_nl_par[1],m_bcw_elec_nl_par[2], m_bcw_elec_nl_par[3], m_bcw_elec_nl_par[4], m_abs_escale * m_rel_escale[idet], m_abs_eoffset+m_rel_eoffset[idet]);
                }
                break;
            case 1://LBNL NL Model
                cout << "Using LBNL NL Model"<< endl;
                nl_func = new TF1("nl_func",NLObject,&CovarianceMatrix2::nl_func_lbnl,InitialVisibleEnergy,FinalVisibleEnergy,5,"CovarianceMatrix2","nl_func_lbnl");
                for (Int_t idet = 0; idet < ADsEH3; idet++)
                {
                    nl_func->SetParameters(m_lbnl_nl_par[0],m_lbnl_nl_par[1],m_lbnl_nl_par[2],m_abs_escale * m_rel_escale[idet],m_abs_eoffset + m_rel_eoffset[idet]);
                }
                
                break;
            case 2://IHEP NL Model
                break;
            default:
        }

        //Resolution function set up
        CovarianceMatrix2* ResoObject = new CovarianceMatrix2();//Finally the standard constructor is useful
        
        reso_func = new TF1("reso_func",ResoObject,&CovarianceMatrix2::reso_func_bcw,0,20,3,"CovarianceMatrix2","reso_func_bcw");
        reso_func->SetParameters(0.022,0.077,0.018); // based on Bryce's TN
        //BCW values http://dayabay.ihep.ac.cn/DocDB/0087/008768/013/6AdAnalysis-BCW.pdf are different! 0.13226034931, 0.32604048828, 0.26196488314
        
        //Resolution
        resolutionRange = 8; // Why 8σ? Seems chosen trivially but in my opinion it's the way to limit the range of the convolution, for a Normal distribution this covers up to 99.99999... of the area

        for(Int_t i=0;i<NADs;i++)
        {
            m_detectorResolution_bias[i] = 0;//this will be used to distort it
        }
        
        //Energy shift function set up
        VisibleF=GetNeutrinoToVisibleFunction(0);//0 for 0th order, 1 for 1st order
        NLInterpolation();
        for(Int_t TrueEnergyIndex=0; TrueEnergyIndex<TotalBins; TrueEnergyIndex++)
        {
            GetEnergyShift(TrueEnergyIndex);
            GetIAVShift(TrueEnergyIndex);
            GetNLShift(TrueEnergyIndex);
            GetResolutionShift(TrueEnergyIndex);
            NormalizeEnergyMatrix(TrueEnergyIndex);
        }
        CreateEnergyMatrix();

        FluctuateBackgrounds();
        GenerateCovarianceMatrix();
        CovMatrix2H->Add(BkgCov2H);
        
        delete nl_func;
        delete reso_func;
        delete VisibleF;
        ResoObject->~CovarianceMatrix2();
        NLObject->~CovarianceMatrix2();

        for (Int_t idet = 0; idet<ADsEH3; idet++)
        {
            delete EnergyMatrixH[idet];
            for (Int_t TrueEnergyIndex= 0; TrueEnergyIndex<TotalBins; TrueEnergyIndex++)
            {
                delete PositronTrueSpectrumH[idet][TrueEnergyIndex];
                delete PositronIAVSpectrumH[idet][TrueEnergyIndex];
                delete PositronNLSpectrumH[idet][TrueEnergyIndex];
                delete PositronVisibleSpectrumH[idet][TrueEnergyIndex];
                delete EnergySlice[idet][TrueEnergyIndex];
            }
        }
    }
    
    CovMatrix2H->Scale(1./(NSamples));
   // CovMatrix2H->Draw("colz");
    SaveTotalSpectrum();
    SaveCovarianceMatrix();
}

void CovarianceMatrix2 :: GenerateCovarianceMatrix()
{    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        //Add nominal backgrounds
        for(Int_t AD=0; AD<NADs; AD++)
        {
            BackgroundSpectrumH[AD][week]=(TH1F*)AccidentalsH[0][0]->Clone();
            BackgroundSpectrumH[AD][week]->Reset();

            if(AddBackgrounds)//Nominal backgrounds for Fpred
            {
                BackgroundSpectrumH[AD][week]->Add(AccidentalsH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(LiHeH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(FastNeutronsH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(AmCH[AD][week]);
            }
        }
        for (Int_t near=0; near<(ADsEH1+ADsEH2); near++)
        {
            RandomNearHallSpectrumH[near][week]=(TH1F*)NearHallSpectrumH[near][week]->Clone(Form("Corrected Random Near AD%d",near));
            RandomNearHallSpectrumH[near][week]->Reset();
            
//Add detector effects //Multiply Matrix and NearHallSpectrum
            for(Int_t i=0;i<Nbins;i++)
            {
                Int_t Sum = 0;

                for(Int_t j=0; j<Nbins; j++)
                {
                    Sum = Sum+(EnergyMatrixH[near]->GetBinContent(i+InitialBins+1,j+InitialBins+1)*NearHallSpectrumH[near][week]->GetBinContent(i+1));
//                    cout << "SUM " << Sum << endl;
//                        cout << " i: " << i <<  "; j: " << j << endl;
                }
                RandomNearHallSpectrumH[near][week]->SetBinContent(i+1,Sum);
            }
            
            if(AddBackgrounds)//Add fluctuated backgrounds for Fobs
            {
                RandomNearHallSpectrumH[near][week]->Add(RandomAccidentalsH[near][week]);
                RandomNearHallSpectrumH[near][week]->Add(RandomLiHeH[near][week]);
                RandomNearHallSpectrumH[near][week]->Add(RandomFastNeutronsH[near][week]);
                RandomNearHallSpectrumH[near][week]->Add(RandomAmCH[near][week]);
            }
//Add detector effects //Multiply Matrix and Prediction (Far Spectrum)
            for (Int_t far=0; far<ADsEH3; far++)//ADsEH3
            {
                RandomPredictionH[far][near][week]=(TH1F*)PredictionH[far][near][week]->Clone(Form("Corrected Random Prediction %d from Near AD%d",far,near));
                RandomPredictionH[far][near][week]->Reset();

//                RandomPredictionH1[far][near][week]=(TH1F*)PredictionH[far][near][week]->Clone();
                
                for(Int_t i=0;i<Nbins;i++)
                {
                    Int_t Sum = 0;
                    
                    for(Int_t j=0; j<Nbins; j++)
                    {
                          Sum = Sum+(EnergyMatrixH[far]->GetBinContent(i+InitialBins+1,j+InitialBins+1)*PredictionH[far][near][week]->GetBinContent(i+1));
//                        cout << "SUM " << Sum << endl;
//                        cout << " i: " << i <<  "; j: " << j << endl;
                    }
                    if (Sum<0)
                    {
                        Sum=0;//Fix temporary negative sum in the last bins
                    }
                    RandomPredictionH[far][near][week]->SetBinContent(i+1,Sum);
                }
                
                if(AddBackgrounds)
                {
                    RandomPredictionH[far][near][week]->Add(RandomAccidentalsH[ADsEH1+ADsEH2+far][week]);
                    RandomPredictionH[far][near][week]->Add(RandomLiHeH[ADsEH1+ADsEH2+far][week]);
                    RandomPredictionH[far][near][week]->Add(RandomFastNeutronsH[ADsEH1+ADsEH2+far][week]);
                    RandomPredictionH[far][near][week]->Add(RandomAmCH[ADsEH1+ADsEH2+far][week]);

                    FluctuateBackgrounds();

//                    RandomPredictionH1[far][near][week]->Add(RandomAccidentalsH[ADsEH1+ADsEH2+far][week]);
//                    RandomPredictionH1[far][near][week]->Add(RandomLiHeH[ADsEH1+ADsEH2+far][week]);
//                    RandomPredictionH1[far][near][week]->Add(RandomFastNeutronsH[ADsEH1+ADsEH2+far][week]);
//                    RandomPredictionH1[far][near][week]->Add(RandomAmCH[ADsEH1+ADsEH2+far][week]);
                }
            }
        }
        
        TFile* CorrectedPredictionsF = TFile::Open("./RootOutputs/DetectorCorrectedSpectrum.root","recreate");
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
            {
                for (Int_t far =0; far<ADsEH3; far++)
                {
                    RandomPredictionH[far][near][week]->Write();
                    RandomNearHallSpectrumH[near][week]->Write();
                }
            }
        }
        CorrectedPredictionsF->Close();
    }

    TH2F* ABkgCov2H = new TH2F("Background Cov Matrix","Background Covariance Matrix",Nbins*9,0,Nbins*9,Nbins*9,0,Nbins*9);
    BkgCov2H=(TH2F*)ABkgCov2H->Clone("Bkg covariance Matrix");
    delete ABkgCov2H;

    Int_t x =0;
    Int_t y =0;
//    Int_t icounter = 1;
//    Int_t jcounter = 1;

    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<ADsEH3; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
//                cout <<"\n"<< "Ni" << Ni1 << Ni2 << Ni3 << Ni4 << "\n";
//                cout << "Fi"<< Fi1 << Fi2 << Fi3 << Fi4 << "\n";
//                cout << "====================================" << icounter++ <<"\n";
//                cout << Ni1*Fi1 << Ni2*Fi2 << Ni3*Fi3 << Ni4*Fi4;
                
                for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}

                    for (Int_t farj=0; farj<ADsEH3; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
//                        cout <<"\n"<< "Nj"<< Nj1 << Nj2 << Nj3 << Nj4 << "\n";
//                        cout << "Fj"<< Fj1 << Fj2 << Fj3 << Fj4 << "\n";
//                        cout << "===================================="<< jcounter++ <<"\n";
//                        cout << Nj1*Fj1 << Nj2*Fj2 << Nj3*Fj3 << Nj4*Fj4 <<"\n";
//                        cout <<"\n";
                        
                        for (Int_t i = 0; i<Nbins; i++)
                        {//columns
                            //cout << "x = " << i << "+ ("<< Ni1*Fi1<< "+"<< Ni2*Fi2 << "+"<< Ni3*Fi3 << "+"<< Ni4*Fi4 << "-1)*41"<< " that is "<< i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*TotalBins << "\n" ;

                            for (Int_t j = 0; j<Nbins; j++)
                            {//rows
                                x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*Nbins;
                               
                                CovBkg[x][y][week]=(RandomPredictionH[fari][neari][week]->GetBinContent(i+1)-PredictionH[fari][neari][week]->GetBinContent(i+1)-BackgroundSpectrumH[fari+ADsEH1+ADsEH2][week]->GetBinContent(i+1))*(RandomPredictionH[farj][nearj][week]->GetBinContent(j+1)-PredictionH[farj][nearj][week]->GetBinContent(j+1)-BackgroundSpectrumH[farj+ADsEH1+ADsEH2][week]->GetBinContent(j+1));
                                
                              //  cout << CovBkg[x][y][week] << "\n";
                            }
                          //  TestGaussianH->Fill(RandomPredictionH[fari][neari][week]->GetBinContent(i+1)-RandomPredictionH1[fari][neari][week]->GetBinContent(i+1));
                          //  TestGaussianH->Fill(RandomPredictionH[fari][neari][week]->GetBinContent(i+1)-PredictionH[fari][neari][week]->GetBinContent(i+1)-BackgroundSpectrumH[fari+ADsEH1+ADsEH2][week]->GetBinContent(i+1));

                        }
                    }
                }
            }
        }
    }
   // TestGaussianH->Draw();
    Int_t x =0;
    Int_t y =0;
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<ADsEH3; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<ADsEH3; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
                        for (Int_t i = 0; i<Nbins; i++)
                        {//columns
                            for (Int_t j = 0; j<Nbins; j++)
                            {//rows
                                x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*Nbins;
                                
                                NormCovBkg[x][y][week]=(CovBkg[x][y][week])/(sqrt(CovBkg[x][x][week]*CovBkg[y][y][week]));
                                if(CovBkg[x][y][week]==0&&((CovBkg[x][x][week]||CovBkg[y][y][week])==0))
                                {
                                    NormCovBkg[x][y][week]=0;//To avoid (0/0) nans when bin contents are practically the same; (this happens when backgrounds are not varied)
                                }
                                BkgCov2H->SetBinContent(x+1,y+1,NormCovBkg[x][y][week]);
                              //  cout <<NormCovBkg[x][y][week];
                               // cout << " " << CovBkg[x][y][week] << " " <<  CovBkg[x][x][week] << " " << CovBkg[y][y][week] << " " << sqrt(CovBkg[x][x][week]*CovBkg[y][y][week]) << "\n";

                            }
                        }
                    }
                }
            }
        }
    }
    
//    BkgCov2H->Draw("colz");
}

////////////////////////////////////////////////////////////////////////////////////
//// I need to correct backgrounds for efficiencies, and maybe events too (Check)
////////////////////////////////////////////////////////////////////////////////////
void CovarianceMatrix2 :: GenerateStatisticalCovarianceMatrix()
{
    for (Int_t week = 0; week<Nweeks; week++)
    {
        //Add nominal backgrounds
        for(Int_t AD=0; AD<NADs; AD++)
        {
            BackgroundSpectrumH[AD][week]=(TH1F*)AccidentalsH[0][0]->Clone();
            BackgroundSpectrumH[AD][week]->Reset();

            if(AddBackgrounds)
            {
                BackgroundSpectrumH[AD][week]->Add(AccidentalsH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(LiHeH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(FastNeutronsH[AD][week]);
                BackgroundSpectrumH[AD][week]->Add(AmCH[AD][week]);
            }
        }
        for (Int_t far=0; far<ADsEH3; far++)
        {
            for (Int_t near=0; near<(ADsEH1+ADsEH2); near++)
            {
                for (Int_t pts = 0; pts < Nbins; pts++)
                {
                    Sigma_Far[far][near][week][pts]=sqrt(PredictionH[far][near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[ADsEH1+ADsEH2+far][week]->GetBinContent(pts+1));
                    
                    Sigma_Near[far][near][week][pts]=(PredictionH[far][near][week]->GetBinContent(pts+1)/NearHallSpectrumH[near][week]->GetBinContent(pts+1))*sqrt(NearHallSpectrumH[near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[near][week]->GetBinContent(pts+1));
                }
            }
        }
    }
    
    TH2F* AStatCov2H = new TH2F("Statistical Cov Matrix","Statistical Covariance Matrix",Nbins*9,0,Nbins*9,Nbins*9,0,Nbins*9);
    StatCov2H=(TH2F*)AStatCov2H->Clone("Statistical Covariance Matrix");
    delete AStatCov2H;
    Int_t x =0;
    Int_t y =0;
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<ADsEH3; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<ADsEH3; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
                        for (Int_t i = 0; i<Nbins; i++)
                        {//columns
                            for (Int_t j = 0; j<Nbins; j++)
                            {//rows
                                x= i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*Nbins;
                                
                                //Near component correlated
                                if(neari==nearj && fari!=farj)
                                {
                                    CovStat[x][y][week]=Sigma_Near[fari][neari][week][i]*Sigma_Near[farj][nearj][week][j];
                                }
                                //Far component correlated
                                if(fari==farj && neari!=nearj)
                                {
                                    CovStat[x][y][week]=Sigma_Far[fari][neari][week][i]*Sigma_Far[farj][nearj][week][j];
                                }
                                if(neari==nearj && fari==farj)
                                {
                                //General covariance
                                    CovStat[x][y][week]=(Sigma_Near[fari][neari][week][i]*Sigma_Near[fari][neari][week][j])+(Sigma_Far[farj][nearj][week][i]*Sigma_Far[farj][nearj][week][j]);
                                }
                                //Uncorrelated terms
                                if(neari!=nearj && fari!=farj)
                                {
                                    CovStat[x][y][week]=0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Int_t x =0;
    Int_t y =0;
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
        {
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<ADsEH3; fari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<ADsEH3; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
                        for (Int_t i = 0; i<Nbins; i++)
                        {//columns
                            for (Int_t j = 0; j<Nbins; j++)
                            {//rows
                                x= i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*Nbins;
                                
                                NormCovStat[x][y][week]=(CovStat[x][y][week])/(sqrt(CovStat[x][x][week]*CovStat[y][y][week]));
                                StatCov2H->SetBinContent(x+1,y+1,NormCovStat[x][y][week]);
                            }
                        }
                    }
                }
            }
        }
    }
    
  // StatCov2H->Draw("colz");

    if(WriteOutput)
    {
        ofstream outf("CovarianceMatrices/StatisticalCovarianceMatrix.txt");
        Int_t x =0;
        Int_t y =0;
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
                if(neari==1){Ni2++;}
                if(neari==2){Ni3++;}
                if(neari==3){Ni4++;}
                
                for (Int_t fari=0; fari<ADsEH3; fari++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
                    if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
                    if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
                    if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
                    
                    for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
                        if(nearj==1){Nj2++;}
                        if(nearj==2){Nj3++;}
                        if(nearj==3){Nj4++;}
                        
                        for (Int_t farj=0; farj<ADsEH3; farj++)
                        {
                            //Logic for the 2D matrix index done up to 8 ADs
                            if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
                            if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
                            if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
                            if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
                            for (Int_t i = 0; i<Nbins; i++)
                            {//columns
                                for (Int_t j = 0; j<Nbins; j++)
                                {//rows
                                    x= i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins;
                                    y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*Nbins;
                                    
                                    outf << NormCovStat[x][y][week] << " ";
                                }
                            }
                        }
                    }
                }
            }
            outf << endl;
        }
        outf.close();
    }
}

void CovarianceMatrix2 :: RandomEnergyScaleMatrix()
{
    
}
void CovarianceMatrix2 :: RandomIAVMatrix()
{
    for (Int_t AD=0; AD<NADs; AD++)
    {
        IAVError[AD]= (1+IAVError*rand->Gaus(0,1)); //Each AD is varied individually.
    }
    
    for(Int_t i=0; i<MatrixBins; i++)
    {
        for(Int_t j=0; i<MatrixBins; i++)
        {
//            cout<<"IAV ORIGINAL" << IAVMatrix[i][j] << endl;
            
            IAVMatrix[i][j]=IAVError[AD]*IAVMatrix[i][j];
            
//            cout<< "IAV ALTERED" << IAVMatrix[i][j] << endl;
        }
    }
}
void CovarianceMatrix2 :: RandomNLMatrix()
{
    switch (NLModelE)
    {
        case 0://BCW NL Model
            //Randomize BCW nonlinear model parameters
            for (Int_t i = 0; i < 5; i++)
            {
                m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i] + m_bcw_elec_nl_par_error[i] * rand->Gaus(0,1);
                cout << " IF ANY OF THIS NUMBERS IS NEGATIVE I HAVE TO CHANGE THIS PART OF THE CODE, CHECK THAT NO PARAMETER IS NEGATIVE" << m_bcw_elec_nl_par[i] << endl;
            }
            break;
        case 1:
            //Randomize LBNL nonlinear model parameters
            for (Int_t i = 0; i < 3; i++)
            {
                m_lbnl_nl_par[i] = m_lbnl_nl_par_nominal[i] + m_lbnl_nl_par_error[i] * rand->Gaus(0,1);
                cout << " IF ANY OF THIS NUMBERS IS NEGATIVE I HAVE TO CHANGE THIS PART OF THE CODE, CHECK THAT NO PARAMETER IS NEGATIVE" <<  m_lbnl_nl_par[i]  << endl;
            }
            break;

        case 2:
            //Randomize IHEP nonlinear model parameters
            break;

        default:
    }
}
void CovarianceMatrix2 :: RandomResolutionMatrix()
{
    
}

void CovarianceMatrix2 :: NLInterpolation()
{
    for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
    {
        IAVEnergy.push_back(nl_func->GetX(VisibleEnergyIndex*BinWidth));

        IAVEnergyIdx.push_back((Int_t)(IAVEnergy[VisibleEnergyIndex]/BinWidth));
//        cout << "IAV Energy " << IAVEnergy[VisibleEnergyIndex] << endl;
//        cout << "IAV Energy Idx " << IAVEnergyIdx[VisibleEnergyIndex] << endl;
//        cout << " VisibleEnergyIndex" << VisibleEnergyIndex << endl;
        
        binScaling.push_back(nl_func->Derivative(IAVEnergy[VisibleEnergyIndex]));
        
//        cout << "BIN SCALING" << binScaling[VisibleEnergyIndex] << endl;
        
        sign.push_back(1);
        dEtrue.push_back(IAVEnergy[VisibleEnergyIndex] - IAVEnergyIdx[VisibleEnergyIndex]*BinWidth);
//        cout << "dEtrue" << dEtrue[VisibleEnergyIndex] << endl;
        if (dEtrue[VisibleEnergyIndex] < 0)
        {
            sign.push_back(-1);
        }
//        cout << "sign" << sign[VisibleEnergyIndex] << endl;
    }
}

void CovarianceMatrix2 :: CreateEnergyMatrix()
{
//    if(PreRebin)
//    {
//        TFile* SaveSpectrumDataF = TFile::Open("./RootOutputs/EnergyMatrixComparison/Energy51.root","recreate");
//    }
//    else
//    {
        TFile* SaveSpectrumDataF = TFile::Open("./RootOutputs/EnergyMatrixComparison/Energy240.root","recreate");
//    }
    // This method should be used with the relative oscillation method
    for(int idet=0; idet<ADsEH3; idet++)
    {
        EnergyMatrixH[idet] = new TH2F(Form("EvisEnu%d",idet),Form("EvisEnu%d",idet),TotalBins,InitialVisibleEnergy,FinalVisibleEnergy,TotalBins,InitialVisibleEnergy,FinalVisibleEnergy);
        
        for(Int_t TrueEnergyIndex=1; TrueEnergyIndex<=TotalBins; TrueEnergyIndex++)
        {
            for(Int_t VisibleBin=1;VisibleBin<=TotalBins;VisibleBin++)
            {
                EnergyMatrixH[idet]->SetBinContent(TrueEnergyIndex,VisibleBin,EnergySlice[idet][TrueEnergyIndex-1]->GetBinContent(VisibleBin));

            }
            // UNCOMMENT FOLLOWING LINES TO SAVE ENERGY MATRIX SLICES IN EACH PRODUCTION STEP
//                        PositronTrueSpectrumH[idet][TrueEnergyIndex-1]->Write();
//                        PositronIAVSpectrumH[idet][TrueEnergyIndex-1]->Write();
//                        PositronNLSpectrumH[idet][TrueEnergyIndex-1]->Write();
//                        PositronVisibleSpectrumH[idet][TrueEnergyIndex-1]->Write();
//                        EnergySlice[idet][TrueEnergyIndex-1]->Write();
        }
        
//        if(!PreRebin)
//        {
//            EnergyMatrixH[idet]->Write();//Before the rebin
            EnergyMatrixH[idet]->Rebin2D(RebinFactor,RebinFactor);
            EnergyMatrixH[idet]->Scale(1./RebinFactor);
//        }
        EnergyMatrixH[idet]->Write();//After the rebin
    }
    SaveSpectrumDataF->Close();
}

void CovarianceMatrix2 :: GetEnergyShift(Int_t TrueEnergyIndex)
{
    Double_t NeutrinoEnergy = (TrueEnergyIndex)*BinWidth;
    PositronEnergy = VisibleF->Eval(NeutrinoEnergy);
    
    if(MaxPositronEnergy==0)
    {
        MaxPositronEnergy = VisibleF->Eval((TotalBins-1)*BinWidth);
        MaxPositronEnergyIndex = Int_t((MaxPositronEnergy-BinWidth)/BinWidth);
    }
    
    PositronEnergyIndex=(Int_t)(PositronEnergy/BinWidth);

    for(int idet=0; idet<ADsEH3; idet++)
    {
        //Calculate Energy Shift        
        PositronTrueSpectrumH[idet][TrueEnergyIndex]= new TH1F(Form("PositronTrueSpectrum%d,%d",idet,TrueEnergyIndex),Form("PositronTrueSpectrum%d,%d",idet,TrueEnergyIndex), TotalBins,InitialVisibleEnergy,FinalVisibleEnergy);
        
        //Reset values to 0
        if (PositronEnergy < 1.022||TrueEnergyIndex>=MaxPositronEnergyIndex)//is about the minimum visible energy (1.017999 for 0th order, 1.01574 for 1st order calculation)
        {
            for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
            {
                PositronTrueSpectrumH[idet][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, 0);
            }
        }
        else
        {
            PositronTrueSpectrumH[idet][TrueEnergyIndex]->SetBinContent(PositronEnergyIndex+1, 1); //Change this to include ADs by setting PredictionH[idet][0][0]->GetBinContent(PositronEnergyIndex+1)
        }
    }
}

void CovarianceMatrix2 :: GetIAVShift(Int_t TrueEnergyIndex)
{
    for(int idet=0; idet<ADsEH3; idet++)
    {
        //IAV
        PositronIAVSpectrumH[idet][TrueEnergyIndex]=(TH1F*)PositronTrueSpectrumH[idet][TrueEnergyIndex]->Clone(Form("IAVSpectrum%d,%d",idet,TrueEnergyIndex));
        PositronIAVSpectrumH[idet][TrueEnergyIndex]->SetTitle("IAV Spectrum");
        //Calculate IAV Shift
        if (PositronEnergy >= 1.022)//1.8-0.78 is about the minimum visible energy (1.017999 for 0th order, 1.01574 for 1st order calculation)
        {
            //Reset values to 0
            for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
            {
               PositronIAVSpectrumH[idet][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, 0);
            }
            
            for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TrueEnergyIndex+1; VisibleEnergyIndex++)
            {
                PositronIAVSpectrumH[idet][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, PositronIAVSpectrumH[idet][TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1) + IAVMatrix[PositronEnergyIndex][VisibleEnergyIndex] * PositronTrueSpectrumH[idet][TrueEnergyIndex]->GetBinContent(PositronEnergyIndex+1));
            }
        }
    }
}

void CovarianceMatrix2 :: GetNLShift(Int_t TrueEnergyIndex)
{
    for(Int_t idet=0; idet<ADsEH3; idet++)
    {
        //NL
        PositronNLSpectrumH[idet][TrueEnergyIndex] = (TH1F*)PositronIAVSpectrumH[idet][TrueEnergyIndex]->Clone(Form("NLSpectrum%d,%d",idet,TrueEnergyIndex));
        PositronNLSpectrumH[idet][TrueEnergyIndex]->SetTitle("NL Spectrum");
        
        //Calculate NL Shift
        for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
        {
            if(IAVEnergyIdx[VisibleEnergyIndex]==TotalBins-1)
            {
                dNdE.push_back(PositronIAVSpectrumH[idet][TrueEnergyIndex]->GetBinContent(IAVEnergyIdx[VisibleEnergyIndex]+1));
            }
            else
            {
                dNdE.push_back((TMath::Abs(dEtrue[VisibleEnergyIndex])/BinWidth)*PositronIAVSpectrumH[idet][TrueEnergyIndex]->GetBinContent(IAVEnergyIdx[VisibleEnergyIndex]+sign[VisibleEnergyIndex]+1)+(1 - TMath::Abs(dEtrue[VisibleEnergyIndex])/BinWidth)* PositronIAVSpectrumH[idet][TrueEnergyIndex]->GetBinContent(IAVEnergyIdx[VisibleEnergyIndex]+1));
            }
        }
        for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
        {
        PositronNLSpectrumH[idet][TrueEnergyIndex]->SetBinContent(VisibleEnergyIndex+1, dNdE[VisibleEnergyIndex]/binScaling[VisibleEnergyIndex]);
        dNdE.pop_back();//Otherwise I keep pushing them down the stack
        }
    }
}

void CovarianceMatrix2 :: GetResolutionShift(Int_t TrueEnergyIndex)
{
    for(Int_t idet=0; idet<ADsEH3; idet++)
    {        
        //Resolution
        PositronVisibleSpectrumH[idet][TrueEnergyIndex] = (TH1F*)PositronNLSpectrumH[idet][TrueEnergyIndex]->Clone(Form("VisibleSpectrum%d,%d",idet,TrueEnergyIndex));
        PositronVisibleSpectrumH[idet][TrueEnergyIndex]->SetTitle("Visible Spectrum");

        //Calculate Resolution effect

        for(Int_t VisibleEnergyIndex=0; VisibleEnergyIndex<TotalBins; VisibleEnergyIndex++)
        {
            Double_t sigma = (reso_func->Eval(VisibleEnergyIndex*BinWidth) + m_detectorResolution_bias[idet]) * VisibleEnergyIndex*BinWidth;//Sigma = FReso(Energy) * Energy
            //                cout << "SIGMA " << sigma<<endl;
            Double_t minDetE = VisibleEnergyIndex*BinWidth - resolutionRange*sigma;
            Double_t maxDetE = VisibleEnergyIndex*BinWidth + resolutionRange*sigma;
            Int_t minDetEIdx = (Int_t)((minDetE)/BinWidth);
            Int_t maxDetEIdx = (Int_t)((maxDetE)/BinWidth);
            
            if(minDetEIdx < 0){minDetEIdx = 0;}
            if(maxDetEIdx >= TotalBins){maxDetEIdx = TotalBins-1;}
            
            for(Int_t detIdx=minDetEIdx; detIdx<=maxDetEIdx; detIdx++)
            {
                if(detIdx==0)
                {
                    continue;
                }
                Double_t gausFactor = TMath::Gaus((VisibleEnergyIndex-detIdx)*BinWidth,0,sigma,true);
                
                PositronVisibleSpectrumH[idet][TrueEnergyIndex]->SetBinContent(detIdx+1, (PositronVisibleSpectrumH[idet][TrueEnergyIndex]->GetBinContent(detIdx+1)+PositronNLSpectrumH[idet][TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1)*gausFactor));
                
               // cout << (PositronVisibleSpectrumH[idet][TrueEnergyIndex]->GetBinContent(detIdx+1)+PositronNLSpectrumH[idet][TrueEnergyIndex]->GetBinContent(VisibleEnergyIndex+1)*gausFactor) << "CHECK THIS" << endl;
            }
        }
        Norm[idet][TrueEnergyIndex] = PositronVisibleSpectrumH[idet][TrueEnergyIndex]->Integral();
    }
}

void CovarianceMatrix2 :: NormalizeEnergyMatrix(Int_t TrueEnergyIndex)
{
    for(Int_t idet=0; idet<ADsEH3; idet++)
    {
        EnergySlice[idet][TrueEnergyIndex] = PositronVisibleSpectrumH[idet][TrueEnergyIndex]->Clone(Form("Energy Slice%d,%d",idet,TrueEnergyIndex));
        if(Norm[idet][TrueEnergyIndex]!=0)
        {
            EnergySlice[idet][TrueEnergyIndex]->Scale(1./Norm[idet][TrueEnergyIndex]);
        }
        else
        {
            EnergySlice[idet][TrueEnergyIndex]->Reset();
        }
    }
}

TF1* CovarianceMatrix2 :: GetNeutrinoToVisibleFunction(Int_t order)
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

void CovarianceMatrix2 :: LoadNLParameters()
{
    switch (NLModelE)
    {
        case 0://BCW NL Model
            cout << "LOAD BCW NL PARAMETERS" << endl;
            ifstream bcw_positron_data("bcw_nl_data/positron.dat");
            for (Int_t i = 0; i < n_bcw_positron_nl; i++)
            {
                bcw_positron_data >> m_bcw_positron_nl_e[i] >> m_bcw_positron_nl_fac[i];
            }
            bcw_positron_data.close();
            
            ifstream bcw_elec_data("bcw_nl_data/par.dat");
            for (Int_t i = 0; i < 3; i++)
            {
                bcw_elec_data >> m_bcw_elec_nl_par_nominal[i]  >> m_bcw_elec_nl_par_error[i];
                m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];
            }
            bcw_elec_data.close();
            
            for (Int_t i = 3; i < 5; i++)
            { // those describe additional uncertainty for bcw model that are described in Doc-XXXX
                m_bcw_elec_nl_par[i] = 0;
                m_bcw_elec_nl_par[i] = m_bcw_elec_nl_par_nominal[i];
                m_bcw_elec_nl_par_error[i] = 1;
            }
            
            TFile *bcw_ele_err_file = new TFile("bcw_nl_data/ele_err.root");
            g_bcw_elec_nl_error[0] = (TGraph*)bcw_ele_err_file->Get("g_up")->Clone();
            g_bcw_elec_nl_error[1] = (TGraph*)bcw_ele_err_file->Get("g_down")->Clone();
            bcw_ele_err_file->Close();
            break;
        case 1://LBNL NL Model
            cout << "LOAD LBNL NL PARAMETERS" << endl;
            ifstream lbnl_positron_data("lbnl_nl_data/lbnl_positron_nl.txt");
            if (!lbnl_positron_data.is_open())
            {
                cout << "Error: cannot find LBNL non-linearity curve!!!" << endl;
                exit(0);
            }
            for (Int_t i = 0; i < 300; i++)
            {
                lbnl_positron_data >> m_lbnl_positron_nl_e[i] >> m_lbnl_positron_nl_fac[i] >> total_err >> m_lbnl_positron_nl_err[0][i] >> m_lbnl_positron_nl_err[1][i] >> m_lbnl_positron_nl_err[2][i];
            }
            for (Int_t i = 0; i < 3; i++)
            {
                m_lbnl_nl_par[i] = 0;
                m_lbnl_nl_par_nominal[i] = 0;
                m_lbnl_nl_par_error[i] = 1.0;
            }
            break;
        case 2://IHEP NL Model
            break;
        default:
    }
}

void CovarianceMatrix2 :: LoadBackgrounds()
{
    TFile* BackgroundsF = TFile::Open("./BackgroundSpectrum/Backgrounds.root");
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t AD =0; AD<NADs; AD++)
        {
            AccidentalsH[AD][week]= (TH1F*)gDirectory->Get(Form("Accidentals_AD%i",AD+1));
            LiHeH[AD][week]=(TH1F*)gDirectory->Get("hist_Bkg_StrongAmC");//Missing LiHe inputs so far
            FastNeutronsH[AD][week]=(TH1F*)gDirectory->Get("FN");
            AmCH[AD][week]=(TH1F*)gDirectory->Get("hist_Bkg_StrongAmC");
        }
    }
      BackgroundsF->Close();
}

//Randomize the events (rate) inside errors taking into account the corresponding correlations between backgrounds and ADs. Also random shape variations are included using distortion functions.
void CovarianceMatrix2 :: FluctuateBackgrounds()
{
    for (Int_t week = 0; week<Nweeks; week++)
    {
            for (Int_t AD =0; AD<NADs; AD++)
            {
                ScaleFactorAccidental[AD][week]=1;
                ScaleFactorLiHe[AD][week]=1;
                ScaleFactorFastNeutrons[AD][week]=1;
                ScaleFactorAmC[AD][week]=1;
                
                hall=2;
                if(AD<ADsEH1)
                {
                    hall=1;
                }
                if(AD>=ADsEH1+ADsEH2)
                {
                    hall=3;
                }
                
                if(AccidentalMatrix)
                {
                    rand->SetSeed(0);
                    ScaleFactorAccidental[AD][week]=(1+HAccidentalError[AD]*rand->Gaus(0,1));
                }
                if(LiHeMatrix)
                {
                    rand->SetSeed(0);
                    ScaleFactorLiHe[AD][week]=(1+HLiHeError[hall]*rand->Gaus(0,1));
                }
                if(FastNeutronsMatrix)
                {
                    rand->SetSeed(0);
                    ScaleFactorFastNeutrons[AD][week]=(1+HFastNeutronsError[hall]*rand->Gaus(0,1));
                }
                if(AmCMatrix)
                {
                    rand->SetSeed(0);
                    ScaleFactorAmC[AD][week]=(1+HAmCError[hall]*rand->Gaus(0,1));
                }
            }
            
            //Accidentals uncorrelated
            //LiHe and Fast neutrons correlated by hall
            
            ScaleFactorLiHe[1][week]=ScaleFactorLiHe[0][week];
            ScaleFactorLiHe[4][week]=ScaleFactorLiHe[3][week];
            ScaleFactorLiHe[5][week]=ScaleFactorLiHe[3][week];
            ScaleFactorFastNeutrons[1][week]=ScaleFactorFastNeutrons[0][week];
            ScaleFactorFastNeutrons[4][week]=ScaleFactorFastNeutrons[3][week];
            ScaleFactorFastNeutrons[5][week]=ScaleFactorFastNeutrons[3][week];
            
            //AmC all correlated
            
            ScaleFactorAmC[1][week]=ScaleFactorAmC[0][week];
            ScaleFactorAmC[2][week]=ScaleFactorAmC[0][week];
            ScaleFactorAmC[3][week]=ScaleFactorAmC[0][week];
            ScaleFactorAmC[4][week]=ScaleFactorAmC[0][week];
            ScaleFactorAmC[5][week]=ScaleFactorAmC[0][week];
            
            for (Int_t AD =0; AD<NADs; AD++)
            {
                RandomAccidentalsH[AD][week]=(TH1F*)AccidentalsH[AD][week]->Clone();
                RandomLiHeH[AD][week]=(TH1F*)LiHeH[AD][week]->Clone();
                RandomFastNeutronsH[AD][week]=(TH1F*)FastNeutronsH[AD][week]->Clone();
                RandomAmCH[AD][week]=(TH1F*)AmCH[AD][week]->Clone();
                
                if(VariateRate)
                {
                    RandomAccidentalsH[AD][week]->Scale(ScaleFactorAccidental[AD][week]);
                    RandomLiHeH[AD][week]->Scale(ScaleFactorLiHe[AD][week]);
                    RandomFastNeutronsH[AD][week]->Scale(ScaleFactorFastNeutrons[AD][week]);
                    RandomAmCH[AD][week]->Scale(ScaleFactorAmC[AD][week]);
                   // RandomAmCH[AD][week]->Draw("same");
                    //            RandomLiHeH[AD][week]->Draw("same");
                     //           RandomFastNeutronsH[AD][week]->Draw("same");
                    // RandomAccidentalsH[AD][week]->Draw("same");

                }

            }

            if(DistortBackgrounds)
            {
                if(AccidentalMatrix)
                {
                    TF1* func_acc=GetDistortionFunction(DistortAcc);
                    
                    for(Int_t iAD=0;iAD<NADs;iAD++)
                    {
                        RandomAccidentalsH[iAD][week]->Multiply(func_acc);
                        RandomAccidentalsH[iAD][week]->Scale(RandomAccidentalsH[iAD][week]->Integral()/AccidentalsH[iAD][week]->Integral());
                        //      cout<<RandomAccidentalsH[AD][week]->Integral()/AccidentalsH[AD][week]->Integral()<<"\n";
                        //      func_acc->Draw();
                    }
//                    RandomAccidentalsH[0][week]->Draw("same");

                    delete func_acc;
                }
                if(LiHeMatrix)
                {
                    TF1* func_LiHe=GetDistortionFunction(DistortLiHe);
                    for(Int_t iAD=0;iAD<NADs;iAD++)
                    {
                        RandomLiHeH[iAD][week]->Multiply(func_LiHe);
                        RandomLiHeH[iAD][week]->Scale(RandomLiHeH[iAD][week]->Integral()/LiHeH[iAD][week]->Integral());
                    //                    cout<<RandomLiHeH[AD][week]->Integral()/LiHeH[AD][week]->Integral()<<"\n";
                    //      func_LiHe->Draw("same");
                    }
                    delete func_LiHe;
                }
                if(FastNeutronsMatrix)
                {

//FN shape distortions are applied taking into account correlations in each EH
                    TF1* func_FN=GetFastNeutronsDistortionFunction(DistortFN);

                        for(Int_t iAD=0;iAD<ADsEH1;iAD++)
                        {
                            RandomFastNeutronsH[iAD][week]->Multiply(func_FN);
                            RandomFastNeutronsH[iAD][week]->Scale(RandomFastNeutronsH[iAD][week]->Integral()/FastNeutronsH[iAD][week]->Integral());
                        }
                    func_FN=GetFastNeutronsDistortionFunction(DistortFN);

                        for(Int_t iAD=ADsEH1;iAD<ADsEH2+ADsEH1;iAD++)
                        {
                            RandomFastNeutronsH[iAD][week]->Multiply(func_FN);
                            RandomFastNeutronsH[iAD][week]->Scale(RandomFastNeutronsH[iAD][week]->Integral()/FastNeutronsH[iAD][week]->Integral());
                        }
                    func_FN=GetFastNeutronsDistortionFunction(DistortFN);

                        for(Int_t iAD=ADsEH1+ADsEH2;iAD<ADsEH1+ADsEH2+ADsEH3;iAD++)
                        {
                            RandomFastNeutronsH[iAD][week]->Multiply(func_FN);
                            RandomFastNeutronsH[iAD][week]->Scale(RandomFastNeutronsH[iAD][week]->Integral()/FastNeutronsH[iAD][week]->Integral());
                        }
                    
//                    for(Int_t iAD=0;iAD<NADs;iAD++)
//                    {
//                       RandomFastNeutronsH[iAD][week]->Draw("same");
//                    }
                    //  func_FN->Draw("same");
                    //   RandomFastNeutronsH[0][week]->Draw("same");

                    delete func_FN;
                }
                if(AmCMatrix)
                {
                    TF1* func_AmC=GetDistortionFunction(DistortAmC);
                    for(Int_t iAD=0;iAD<NADs;iAD++)
                    {
                        RandomAmCH[iAD][week]->Multiply(func_AmC);
                        RandomAmCH[iAD][week]->Scale(RandomAmCH[iAD][week]->Integral()/AmCH[iAD][week]->Integral());
                       //RandomAmCH[iAD][week]->Draw("same");

                        //     func_AmC->Draw("same");
                    }

                    delete func_AmC;
                }
            }
    }
}

void CovarianceMatrix2 :: LoadPredictions()
{
TFile* FarHallPredictionsF = TFile::Open("./RootOutputs/FarSpectrumFraction.root");
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            for (Int_t far =0; far<ADsEH3; far++)
            {
                PredictionH[far][near][week] = (TH1F*)gDirectory->Get(Form("AD%i Far Spectrum prediction from near AD%i",far+1,near+1));
            }
        }
    }
    FarHallPredictionsF->Close();
}

void CovarianceMatrix2 :: LoadNearHall()
{
    Char_t filenameNear[1024];
    
    TFile* NearHallDataF = TFile::Open("./RootOutputs/NominalOutputs/Oscillation.root"); //This should be real data. It has to be fixed.
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            sprintf(filenameNear,"Total spectrum after oscillation at AD%i",near+1);
            NearHallDataF->cd("Total AD Spectra after oscillation");
            NearHallSpectrumH[near][week] = (TH1F*)gDirectory->Get(filenameNear);
        }
    }
    NearHallDataF->Close();
    
    TFile* PaddedPredictionsF = TFile::Open("./RootOutputs/PaddedSpectrum.root","recreate");
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            for (Int_t far =0; far<ADsEH3; far++)
            {
                PredictionH[far][near][week]->Write();
                NearHallSpectrumH[near][week]->Write();
            }
        }
    }
    PaddedPredictionsF->Close();
}


void CovarianceMatrix2 :: LoadIavCorrection()//From Bryce Littlejohn's results. //I don't use different IAV matrix for each AD since the analysis it's been done for just 1 of them, assume identical ADs.
{
    TFile * f = new TFile("./IavDistortion/IAVDistortion.root");
    TH2F * Correction = (TH2F*)f->Get("Correction");
    
    cout << "Reading IAV correction file" << endl;
    
//    if(PreRebin)
//    {
//        Correction->Rebin2D(RebinFactor,RebinFactor);
//    }
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

void CovarianceMatrix2 :: SaveTotalSpectrum()
{
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            RandomNearHallSpectrumH[near][week]->Add(BackgroundSpectrumH[near][week]);

            for (Int_t far =0; far<ADsEH3; far++)
            {
                RandomPredictionH[far][near][week]->Add(BackgroundSpectrumH[ADsEH1+ADsEH2+far][week]);
            }
        }
     }

    TFile* SaveSpectrumDataF = TFile::Open("./RootOutputs/SpectrumWithBackgrounds.root","recreate");
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            RandomNearHallSpectrumH[near][week]->Write();
            
            for (Int_t far =0; far<ADsEH3; far++)
            {
                RandomPredictionH[far][near][week]->Write();
            }
        }
    }
   SaveSpectrumDataF->Close();
}

void CovarianceMatrix2 :: SaveCovarianceMatrix()
{
    TFile* SaveCovarianceMatrixF = TFile::Open("./CovarianceMatrices/CovarianceMatrixMatrices.root","recreate");
    StatCov2H->Write();
    switch (BackgroundE)
    {
        case 0://Vary Accidentals
            CovMatrix2H->SetName("Accidental Covariance Matrix");
            CovMatrix2H->SetTitle("Accidental Covariance Matrix");
            break;
        case 1://Vary LiHe
            CovMatrix2H->SetName("LiHe Covariance Matrix");
            CovMatrix2H->SetTitle("LiHe Covariance Matrix");
            break;
        case 2://Vary Fast Neutrons
            CovMatrix2H->SetName("FN Covariance Matrix");
            CovMatrix2H->SetTitle("FN Covariance Matrix");
            break;
        case 3://Vary AmC
            CovMatrix2H->SetName("AmC Covariance Matrix");
            CovMatrix2H->SetTitle("AmC Covariance Matrix");
            break;
        default://Add nominal backgrounds
            CovMatrix2H->SetName("Nominal Covariance Matrix");
            CovMatrix2H->SetTitle("Nominal Covariance Matrix");
    }
   // CovMatrix2H->Draw("colz");
    CovMatrix2H->Write();
    SaveCovarianceMatrixF->Close();
    
    if(WriteOutput)
    {
        ofstream outf("CovarianceMatrices/BackgroundCovarianceMatrix.txt");
        Int_t x =0;
        Int_t y =0;
            for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
            {
                //Logic for the 2D matrix index done up to 8 ADs
                if(neari==0){Int_t Ni1=1;Int_t Ni2=0;Int_t Ni3=0;Int_t Ni4=0;}
                if(neari==1){Ni2++;}
                if(neari==2){Ni3++;}
                if(neari==3){Ni4++;}
                
                for (Int_t fari=0; fari<ADsEH3; fari++)
                {
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(Ni1!=Ni2){Int_t Fi1=fari+1;Int_t Fi2 = 0;Int_t Fi3 = 0;Int_t Fi4 = 0;}
                    if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Int_t Fi2=fari+1;}
                    if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Int_t Fi3=fari+1;}
                    if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Int_t Fi4=fari+1;}
                    
                    for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(nearj==0){Int_t Nj1=1;Int_t Nj2=0;Int_t Nj3=0;Int_t Nj4=0;}
                        if(nearj==1){Nj2++;}
                        if(nearj==2){Nj3++;}
                        if(nearj==3){Nj4++;}
                        
                        for (Int_t farj=0; farj<ADsEH3; farj++)
                        {
                            //Logic for the 2D matrix index done up to 8 ADs
                            if(Nj1!=Nj2){Int_t Fj1=farj+1;Int_t Fj2 = 0;Int_t Fj3 = 0;Int_t Fj4 = 0;}
                            if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Int_t Fj2=farj+1;}
                            if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Int_t Fj3=farj+1;}
                            if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Int_t Fj4=farj+1;}
                            
                            for (Int_t i = 0; i<Nbins; i++)
                            {//columns
                                for (Int_t j = 0; j<Nbins; j++)
                                {//rows
                                    x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*Nbins;
                                    y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*Nbins;
                                    
                                    outf << CovMatrix2H->GetBinContent(x+1,y+1)<< " ";
                                }
                            }
                        }
                    }
                }
            }
            outf << endl;
            outf.close();
    }
}

void CovarianceMatrix2 :: SetAccidentalMatrix(bool AccMatrix)
{
    AccidentalMatrix=AccMatrix;
    if(AccMatrix)
    {
        BackgroundE = AccidentalE;
        DistortAcc=0.2;
    }
}

void CovarianceMatrix2 :: SetLiHeMatrix(bool LiMatrix)
{
    LiHeMatrix=LiMatrix;
    if(LiMatrix)
    {
        BackgroundE = LiHeE;
        DistortLiHe=0.2;
    }
}

void CovarianceMatrix2 :: SetFastNeutronsMatrix(bool FNmatrix)
{
    FastNeutronsMatrix=FNmatrix;
    if(FNmatrix)
    {
        BackgroundE = FastNeutronsE;
        DistortFN=0.2;
    }
}

void CovarianceMatrix2 :: SetAmCMatrix(bool AmCmatrix)
{
    AmCMatrix=AmCmatrix;
    if(AmCMatrix)
    {
        BackgroundE = AmCE;
        DistortAmC=0.2;
    }
}

void CovarianceMatrix2 :: SetIsotopeMatrix(bool IsotopeMatrix)
{
    if(IsotopeMatrix)
    {
        Pred.SetRandomIsotopeFraction(IsotopeMatrix);
        SystematicE = ReactorE;
    }
}

void CovarianceMatrix2 :: SetReactorPowerMatrix(bool ReactorPowerMatrix)
{
    if(ReactorPowerMatrix)
    {
        Pred.SetRandomReactorPower(ReactorPowerMatrix);
        SystematicE = ReactorE;
    }
}

void CovarianceMatrix2 :: SetEnergyScaleMatrix(bool EnergyScaleMatrix)
{
    if(EnergyScaleMatrix)
    {
        RandomEnergyScaleMatrix();
        SystematicE = EnergyE;
    }
}

void CovarianceMatrix2 :: SetIAVMatrix(bool IAVMatrix)
{
    if(IAVMatrix)
    {
        RandomIAVMatrix();
        SystematicE = IAVE;
    }
}

void CovarianceMatrix2 :: SetNLMatrix(bool NLMatrix)
{
    if(NLMatrix)
    {
        RandomNLMatrix();
        SystematicE = NLE;
    }
}

void CovarianceMatrix2 :: SetResolutionMatrix(bool ResolutionMatrix)
{
    if(ResolutionMatrix)
    {
        RandomResolutionMatrix();
        SystematicE = ResolutionE;
    }
}

void CovarianceMatrix2 :: SetBCWModel(bool bcw)
{
    if(bcw)
    {
        NLModelE=BCWE;
    }
}
    
void CovarianceMatrix2 :: SetLBNLModel(bool lbnl)
{
    if(lbnl)
    {
        NLModelE=LBNLE;
    }
}
    
void CovarianceMatrix2 :: SetIHEPModel(bool ihep)
{
    if(ihep)
    {
        NLModelE=IHEPE;
    }
}

TF1* CovarianceMatrix2 :: GetNeutrinoToVisibleFunction(Int_t order)
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

TF1* CovarianceMatrix2 :: GetDistortionFunction(Double_t amount)
{
    TF1 *func = new TF1("func","TMath::Abs([0]+[1]*x)",InitialEnergy,FinalEnergy);
    rand->SetSeed(0);
    Double_t slope=amount*rand->Gaus(0,1);
    Double_t anchor_point=3.5;
    //want offset to be set by requiring func at anchor point to be 1
    Double_t offset=(1-slope*anchor_point);
    func->SetParameter(0,offset);
    func->SetParameter(1,slope);
    
    return func;
}

TF1* CovarianceMatrix2 :: GetFastNeutronsDistortionFunction(Double_t amount)
{
    TF1 *func = new TF1("func","[0]/(0.2*pow(x,0.1))+[1]",InitialEnergy,FinalEnergy);
    rand->SetSeed(0);
    Double_t scaling =amount*rand->Gaus(0,1);
    func->SetParameter(0,scaling);
    // func->SetParameter(1,0); // see if it works
    //set offset so that func(FinalEnergy)=1;
    func->SetParameter(1,1-1*func->Eval(10));//I'm using 10 instead of 12 because when I fit the FN I use the range 10-100 MeV.
    
    return func;
}

// First order true to visible function.
//(x-(Mn-Mp))*(1 - x/Mn*(1.0 - [0]*sqrt(1 - Me*Me/(x-(Mn-Mp))/(x-(Mn-Mp)))))- ((Mn-Mp)*(Mn-Mp) - Me*Me)/(2*Mn)-Me]
Double_t VisibleEnergy1F(Double_t* x, Double_t* par)
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

Double_t VisibleEnergy0F(Double_t* x, Double_t* par)
{
    return x[0]-par[0];
}

Double_t CovarianceMatrix2 :: reso_func_bcw(Double_t* x, Double_t* par)
{
    Double_t e_orig = x[0];
    Double_t e_sigma = 1.0;
    
    if (e_orig > 0)//To avoid NaNs.
    {
        e_sigma = TMath::Sqrt(par[0]*par[0] + par[1]*par[1]/e_orig + par[2]*par[2]/e_orig/e_orig);
    }
    
    return e_sigma;
}

Double_t CovarianceMatrix2 :: nl_func_bcw(Double_t* x, Double_t* par)
{
    //Input error
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

// LBNL non-linearityr function, based on beta-gamma spectra
Double_t CovarianceMatrix2 :: nl_func_lbnl(Double_t * x, Double_t * par)
{
    
    //par[0]: flat energy scale shift for electron
    //par[1]: size of energy scale shift propotional to exp(-1.5*eVis) for electron
    //par[2]: size of energy scale shift due to Ge68 calibration point
    
    Double_t e_positron_true = x[0];
    Double_t escale_par = par[3]; // Add flat energy scale parameter
    Double_t escale_offset = par[4]; // Add fixed energy offset
    
    Double_t scinti_nl_fac = 1;
    Double_t err[3];
    
    if (e_positron_true < m_lbnl_positron_nl_e[0])
    {
        scinti_nl_fac = m_lbnl_positron_nl_fac[0];
        for (Int_t ierr = 0; ierr < 3; ierr++)
        {
            err[ierr] = m_lbnl_positron_nl_err[ierr][0];
        }
    }
    else if (e_positron_true > m_lbnl_positron_nl_e[299])
    {
        scinti_nl_fac = m_lbnl_positron_nl_fac[299];
        for (Int_t ierr = 0; ierr < 3; ierr++)
        {
            err[ierr] = m_lbnl_positron_nl_err[ierr][299];
        }
    }
    else
    {
        for (Int_t i = 0; i < 299; i++)
        {
            if (e_positron_true >= m_lbnl_positron_nl_e[i] && e_positron_true < m_lbnl_positron_nl_e[i+1])
            {
                scinti_nl_fac = ((m_lbnl_positron_nl_e[i+1] - e_positron_true)*m_lbnl_positron_nl_fac[i] + (e_positron_true - m_lbnl_positron_nl_e[i])*m_lbnl_positron_nl_fac[i+1]) / (m_lbnl_positron_nl_e[i+1] - m_lbnl_positron_nl_e[i]);
                for (Int_t ierr = 0; ierr < 3; ierr++)
                {
                    err[ierr] = ((m_lbnl_positron_nl_e[i+1] - e_positron_true)*m_lbnl_positron_nl_err[ierr][i] + (e_positron_true - m_lbnl_positron_nl_e[i])*m_lbnl_positron_nl_err[ierr][i+1]) / (m_lbnl_positron_nl_e[i+1] - m_lbnl_positron_nl_e[i]);
                }
                break;
            }
        }
        
    }
    double random_nl_fac = scinti_nl_fac;
    
    
    for (Int_t ierr = 0; ierr < 3; ierr++)
    {
        random_nl_fac += par[ierr]*err[ierr];
    }
    
    double visibleE = random_nl_fac * e_positron_true;
    //  cout << "NL " << e_positron_true << " " << visibleE/e_positron_true << endl;
    
    return visibleE  * escale_par + escale_offset;
}
