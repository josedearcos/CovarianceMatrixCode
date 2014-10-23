class StatisticalCovarianceMatrix
{
private:
    NominalData* Nom;
    
    bool SimpleReactorModel;
    
    //AD configuration parameters:
    Int_t NADs;
    Int_t ADsEH1;
    Int_t ADsEH2;
    Int_t ADsEH3;
    
    //Binning parameters:
    bool LinearBinning;
    Int_t n_evis_bins;
    Int_t n_etrue_bins;

    Double_t InitialEnergy;
    Double_t FinalEnergy;
    Double_t InitialVisibleEnergy;
    Double_t FinalVisibleEnergy;
    Int_t Nweeks;
    
    TH1F* OriginalPredictionH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    TH1F* PredictionTrueH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    TH1F* PredictionVisH[MaxFarDetectors][MaxNearDetectors][MaxPeriods];
    
    TH1F* OriginalNearHallSpectrumH[MaxNearDetectors][MaxPeriods];
    TH1F* NearHallSpectrumTrueH[MaxNearDetectors][MaxPeriods];
    TH1F* NearHallSpectrumVisH[MaxNearDetectors][MaxPeriods];
    
    TH1F* BackgroundSpectrumH[MaxDetectors][MaxPeriods];
    TH1F* AccidentalsH[MaxDetectors][MaxPeriods];
    TH1F* LiHeH[MaxDetectors][MaxPeriods];
    TH1F* FastNeutronsH[MaxDetectors][MaxPeriods];
    TH1F* AmCH[MaxDetectors][MaxPeriods];
    
    TH2F* NominalResponseMatrixH[MaxDetectors];

    Double_t evis_bins[MaxNbins+1]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins +1 for the 12MeV limit.
    Double_t enu_bins[MaxNbins+1]; // 39 bins between 1.8 and 9.6 MeV +1 for the 9.6 limit.

    Double_t Sigma_Near[MaxFarDetectors][MaxNearDetectors][MaxPeriods][MaxNbins];
    Double_t Sigma_Far[MaxFarDetectors][MaxNearDetectors][MaxPeriods][MaxNbins];
  
    Double_t CovStat[9*MaxNbins][9*MaxNbins][MaxPeriods];
    Double_t NormCovStat[9*MaxNbins][9*MaxNbins][MaxPeriods];
    TH2F* StatCov2H;

    void GenerateStatisticalCovarianceMatrix();
    void LoadPredictions();
    void LoadBackgrounds();
    void LoadNearHall();
    void LoadResponseMatrix();
    void SaveStatisticalCovarianceMatrix();
    
public:
    StatisticalCovarianceMatrix();
    StatisticalCovarianceMatrix(NominalData*);
    
    void StatisticalCovarianceMatrixMain();
};

StatisticalCovarianceMatrix :: StatisticalCovarianceMatrix()
{

}
StatisticalCovarianceMatrix :: StatisticalCovarianceMatrix(NominalData* Data)
{
    SimpleReactorModel = Data->GetSimpleReactorModel();

    n_evis_bins = Data->GetNbins();
    n_etrue_bins = Data->GetNbins();

    InitialEnergy = Data->GetEmin();
    FinalEnergy = Data->GetEmax();
    InitialVisibleEnergy =Data->GetEVisMin();
    FinalVisibleEnergy = Data->GetEVisMax();
    Nweeks = Data->GetWeeks();

    LinearBinning = Data->GetBinning();
    
    //Linear binning
    if(LinearBinning)
    {
        for (Int_t i = 0; i <= n_evis_bins; i++)
        {
            evis_bins[i] = 0.2 * i + InitialEnergy;
            enu_bins[i] = 0.2 * i + InitialEnergy;
        }
    }
    //Non-linear binning
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

    NADs = Data->GetADs();
    ADsEH1 = 2;
    ADsEH2 = 1;
    ADsEH3 = 3;
    
    if(NADs == 8) //ADs can only be 6 or 8
    {
        ADsEH2 = 2;
        ADsEH3 = 4;
    }
}
void StatisticalCovarianceMatrix :: StatisticalCovarianceMatrixMain()
{
    LoadPredictions();
    LoadNearHall();//Order matters, this have to be after LoadPredictions()
    LoadBackgrounds();
    LoadResponseMatrix();
    
    StatCov2H = new TH2F("Statistical Covariance Matrix","Statistical Covariance Matrix",n_evis_bins*9,0,n_evis_bins*9,n_evis_bins*9,0,n_evis_bins*9);
    GenerateStatisticalCovarianceMatrix();
    SaveStatisticalCovarianceMatrix();
    delete StatCov2H;
}
////////////////////////////////////////////////////////////////////////////////////
//// I need to correct backgrounds for efficiencies, and maybe events too (Check)
////////////////////////////////////////////////////////////////////////////////////
void StatisticalCovarianceMatrix :: GenerateStatisticalCovarianceMatrix()
{
    TFile* Check = TFile::Open("./StatisticalSpectrum.root","recreate");

    //Nominal in this one:
    for (Int_t week = 0; week<Nweeks; week++)
    {
        Double_t SumNear,SumFar;
        for (Int_t near=0; near<(ADsEH1+ADsEH2); near++)
        {
            NearHallSpectrumVisH[near][week] = new TH1F(Form("Near Vis Nominal AD%d",near),Form("Near Vis Nominal AD%d",near),n_evis_bins,evis_bins);
            
            for(Int_t j=1; j<=n_evis_bins; j++)
            {
                SumNear = 0;
                
                for(Int_t i=1;i<=n_etrue_bins;i++)
                {
                    SumNear = SumNear+(NominalResponseMatrixH[0]->GetBinContent(i,j)*NearHallSpectrumTrueH[near][week]->GetBinContent(i));
                }
                NearHallSpectrumVisH[near][week]->SetBinContent(j,SumNear);
            }
            NearHallSpectrumVisH[near][week]->Write();

            for (Int_t far=0; far<ADsEH3; far++)
            {
                PredictionVisH[far][near][week] = new TH1F(Form("Far Nominal Visible Spectrum at AD%d from AD%d",far+1,near+1),Form("Far Nominal Visible Spectrum at AD%d from AD%d",far+1,near+1),n_evis_bins,evis_bins);
                for(Int_t j=1; j<=n_evis_bins; j++)
                {
                    SumFar = 0;
                    for(Int_t i=1;i<=n_etrue_bins;i++)
                    {
                        SumFar = SumFar+(NominalResponseMatrixH[0]->GetBinContent(i,j)*PredictionTrueH[far][near][week]->GetBinContent(i));
                    }
                    PredictionVisH[far][near][week]->SetBinContent(j,SumFar);
                }
                PredictionVisH[far][near][week]->Write();
//                CopyPredictionVisH[far][near][week]=(TH1F*)PredictionVisH[far][near][week]->Clone();
            }
        }
        Check->Close();
        //Add nominal backgrounds
        for(Int_t AD=0; AD<NADs; AD++)
        {
            BackgroundSpectrumH[AD][week]=(TH1F*)AccidentalsH[0][0]->Clone();
            BackgroundSpectrumH[AD][week]->Reset();
            
            BackgroundSpectrumH[AD][week]->Add(AccidentalsH[AD][week]);
            BackgroundSpectrumH[AD][week]->Add(LiHeH[AD][week]);
            BackgroundSpectrumH[AD][week]->Add(FastNeutronsH[AD][week]);
            BackgroundSpectrumH[AD][week]->Add(AmCH[AD][week]);
        }
        for (Int_t far=0; far<ADsEH3; far++)
        {
            for (Int_t near=0; near<(ADsEH1+ADsEH2); near++)
            {
                for (Int_t pts = 0; pts < n_evis_bins; pts++)
                {
                    Sigma_Far[far][near][week][pts]=sqrt(PredictionVisH[far][near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[far][week]->GetBinContent(pts+1));

                    //N OBSERVED, F PREDICTED FROM N OBSERVED, BACKGROUNDS PREDICTED (WITHOUT ANY FLUCTUATION)
                    Sigma_Near[far][near][week][pts]=(PredictionVisH[far][near][week]->GetBinContent(pts+1)/NearHallSpectrumVisH[near][week]->GetBinContent(pts+1))*sqrt(NearHallSpectrumVisH[near][week]->GetBinContent(pts+1)+BackgroundSpectrumH[near][week]->GetBinContent(pts+1));
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
            Int_t Ni1,Ni2,Ni3,Ni4;
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Ni1=1; Ni2=0; Ni3=0; Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<ADsEH3; fari++)
            {
                Int_t Fi1,Fi2,Fi3,Fi4;
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Fi1=fari+1;Fi2 = 0;Fi3 = 0;Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3; Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3; Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3; Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                {
                    Int_t Nj1,Nj2,Nj3,Nj4;
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<ADsEH3; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        Int_t Fj1,Fj2,Fj3,Fj4;
                        if(Nj1!=Nj2){Fj1=farj+1;Fj2 = 0;Fj3 = 0;Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3; Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3; Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3; Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            
                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                
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
    x = 0;
    y = 0;
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
        {
            Int_t Ni1,Ni2,Ni3,Ni4;
            //Logic for the 2D matrix index done up to 8 ADs
            if(neari==0){Ni1=1; Ni2=0; Ni3=0; Ni4=0;}
            if(neari==1){Ni2++;}
            if(neari==2){Ni3++;}
            if(neari==3){Ni4++;}
            
            for (Int_t fari=0; fari<ADsEH3; fari++)
            {
                Int_t Fi1,Fi2,Fi3,Fi4;
                //Logic for the 2D matrix index done up to 8 ADs
                if(Ni1!=Ni2){Fi1=fari+1;Fi2 = 0;Fi3 = 0;Fi4 = 0;}
                if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Fi2=fari+1;}
                if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Fi3=fari+1;}
                if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Fi4=fari+1;}
                
                for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                {
                    Int_t Nj1,Nj2,Nj3,Nj4;
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                    if(nearj==1){Nj2++;}
                    if(nearj==2){Nj3++;}
                    if(nearj==3){Nj4++;}
                    
                    for (Int_t farj=0; farj<ADsEH3; farj++)
                    {
                        //Logic for the 2D matrix index done up to 8 ADs
                        Int_t Fj1,Fj2,Fj3,Fj4;
                        if(Nj1!=Nj2){Fj1=farj+1;Fj2 = 0;Fj3 = 0;Fj4 = 0;}
                        if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3;Fj2=farj+1;}
                        if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3;Fj3=farj+1;}
                        if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3;Fj4=farj+1;}
                        
                        for (Int_t i = 0; i<n_evis_bins; i++)
                        {//columns
                            
                            for (Int_t j = 0; j<n_evis_bins; j++)
                            {//rows
                                x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                                y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                
                                if ((CovStat[x][x][week]*CovStat[y][y][week])==0)
                                {
                                    NormCovStat[x][y][week]=0;
                                }
                                else
                                {
                                    NormCovStat[x][y][week]=(CovStat[x][y][week])/(sqrt(CovStat[x][x][week]*CovStat[y][y][week]));
                                }
                                
                                StatCov2H->SetBinContent(x+1,y+1,NormCovStat[x][y][week]);
                            }
                        }
                    }
                }
            }
        }
    }
    // StatCov2H->Draw("colz");
}

void StatisticalCovarianceMatrix :: SaveStatisticalCovarianceMatrix()
{
    //Save statistical matrix
    TFile* SaveStatisticalCovMatrixF = TFile::Open("./CovarianceMatrices/StatisticalCovarianceMatrix.root","recreate");
    StatCov2H->Write();
    SaveStatisticalCovMatrixF->Close();
    
    //Save in a txt file
    if(WriteOutput)
    {
        ofstream statf("CovarianceMatrices/StatisticalCovarianceMatrix.txt");
        Int_t x =0;
        Int_t y =0;
        for (Int_t week = 0; week<Nweeks; week++)
        {
            for (Int_t neari=0; neari<(ADsEH1+ADsEH2); neari++)
            {
                Int_t Ni1,Ni2,Ni3,Ni4;
                //Logic for the 2D matrix index done up to 8 ADs
                if(neari==0){Ni1=1; Ni2=0; Ni3=0; Ni4=0;}
                if(neari==1){Ni2++;}
                if(neari==2){Ni3++;}
                if(neari==3){Ni4++;}
                
                for (Int_t fari=0; fari<ADsEH3; fari++)
                {
                    Int_t Fi1,Fi2,Fi3,Fi4;
                    //Logic for the 2D matrix index done up to 8 ADs
                    if(Ni1!=Ni2){Fi1=fari+1;Fi2 = 0;Fi3 = 0;Fi4 = 0;}
                    if(Ni1==Ni2&&Ni2!=Ni3){Fi1=ADsEH3;Fi2=fari+1;}
                    if(Ni2==Ni3&&Ni3!=Ni4){Fi2=ADsEH3;Fi3=fari+1;}
                    if(Ni3==Ni4&&Ni4==1){Fi3=ADsEH3;Fi4=fari+1;}
                    
                    for (Int_t nearj=0; nearj<(ADsEH1+ADsEH2); nearj++)
                    {
                        Int_t Nj1,Nj2,Nj3,Nj4;
                        //Logic for the 2D matrix index done up to 8 ADs
                        if(nearj==0){Nj1=1;Nj2=0;Nj3=0;Nj4=0;}
                        if(nearj==1){Nj2++;}
                        if(nearj==2){Nj3++;}
                        if(nearj==3){Nj4++;}
                        
                        for (Int_t farj=0; farj<ADsEH3; farj++)
                        {
                            //Logic for the 2D matrix index done up to 8 ADs
                            Int_t Fj1,Fj2,Fj3,Fj4;
                            if(Nj1!=Nj2){Fj1=farj+1;Fj2 = 0;Fj3 = 0;Fj4 = 0;}
                            if(Nj1==Nj2&&Nj2!=Nj3){Fj1=ADsEH3; Fj2=farj+1;}
                            if(Nj2==Nj3&&Nj3!=Nj4){Fj2=ADsEH3; Fj3=farj+1;}
                            if(Nj3==Nj4&&Nj4==1){Fj3=ADsEH3; Fj4=farj+1;}
                            
                            for (Int_t i = 0; i<n_evis_bins; i++)
                            {//columns
                                
                                for (Int_t j = 0; j<n_evis_bins; j++)
                                {//rows
                                    x = i+(Ni1*Fi1+Ni2*Fi2+Ni3*Fi3+Ni4*Fi4-1)*n_evis_bins;
                                    y = j+(Nj1*Fj1+Nj2*Fj2+Nj3*Fj3+Nj4*Fj4-1)*n_evis_bins;
                                    
                                    statf << NormCovStat[x][y][week] << " ";
                                }
                            }
                        }
                    }
                }
            }
            statf << std::endl;
        }
        statf.close();
    }
}

void StatisticalCovarianceMatrix :: LoadPredictions()
{
    TFile* FarHallPredictionsF = TFile::Open("./RootOutputs/FarSpectrumFraction.root");
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            for (Int_t far =0; far<ADsEH3; far++)
            {
                OriginalPredictionH[far][near][week] = (TH1F*)gDirectory->Get(Form("AD%i Far Spectrum prediction from near AD%i",far+1,near+1));
                PredictionTrueH[far][near][week]= (TH1F*)OriginalPredictionH[far][near][week]->Clone();
                if(SimpleReactorModel)
                {
                    PredictionTrueH[far][near][week]=(TH1F*)OriginalPredictionH[far][near][week]->Rebin(n_etrue_bins,Form("AD%i Far Spectrum prediction from near AD%i",far+1,near+1),enu_bins);
                }
                else
                {
                    PredictionTrueH[far][near][week]=(TH1F*)OriginalPredictionH[far][near][week]->Clone();
                }
            }
        }
    }
    FarHallPredictionsF->Close();
}

void StatisticalCovarianceMatrix :: LoadNearHall()
{
    Char_t filenameNear[100];
    
    TFile* NearHallDataF = TFile::Open("./RootOutputs/NominalOutputs/Oscillation.root"); //This should be real data. It has to be fixed.
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            sprintf(filenameNear,"Total spectrum after oscillation at AD%i",near+1);
            NearHallDataF->cd("Total AD Spectra after oscillation");
            OriginalNearHallSpectrumH[near][week] = (TH1F*)gDirectory->Get(filenameNear);
            NearHallSpectrumTrueH[near][week]= (TH1F*)OriginalNearHallSpectrumH[near][week]->Clone();
            if(SimpleReactorModel)
            {
                NearHallSpectrumTrueH[near][week]=(TH1F*)OriginalNearHallSpectrumH[near][week]->Rebin(n_etrue_bins,Form("AD%i Near Spectrum",near+1),enu_bins);
            }
            else
            {
                NearHallSpectrumTrueH[near][week]=(TH1F*)OriginalNearHallSpectrumH[near][week]->Clone();
            }
        }
    }
    NearHallDataF->Close();
    
    TFile* SpectrumF = TFile::Open("./RootOutputs/TrueSpectrum.root","recreate");
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t near =0; near<(ADsEH1+ADsEH2); near++)
        {
            NearHallSpectrumTrueH[near][week]->Write();
            
            for (Int_t far =0; far<ADsEH3; far++)
            {
                PredictionTrueH[far][near][week]->Write();
            }
        }
    }
    SpectrumF->Close();
}

void StatisticalCovarianceMatrix :: LoadResponseMatrix()
{
    Char_t filenameResponse[100];
    
    Int_t NLModelE = 1;//FIX THIS PROBABLY WITH  Polymorphism
    
    if(LinearBinning==0)
    {
        switch (NLModelE)
        {
            case 0://BCW NL Model
                sprintf(filenameResponse,"./NominalResponseMatrices/LBNLBinning/NominalResponseBCWModel.root");
                break;
            case 1://LBNL NL Model
                sprintf(filenameResponse,"./NominalResponseMatrices/LBNLBinning/NominalResponseLBNLModel.root");
                break;
            case 2://IHEP NL Model
                sprintf(filenameResponse,"./NominalResponseMatrices/LBNLBinning/NominalResponseIHEPModel.root");
                break;
            case 3://Unified NL Model
                sprintf(filenameResponse,"./NominalResponseMatrices/LBNLBinning/NominalResponseUnifiedModel.root");
                break;
        }
    }
    else
    {
        switch (NLModelE)
        {
            case 0://BCW NL Model
                sprintf(filenameResponse,"./NominalResponseMatrices/LinearBinning/NominalResponseBCWModel.root");
                break;
            case 1://LBNL NL Model
                sprintf(filenameResponse,"./NominalResponseMatrices/LinearBinning/NominalResponseLBNLModel.root");
                break;
            case 2://IHEP NL Model
                sprintf(filenameResponse,"./NominalResponseMatrices/LinearBinning/NominalResponseIHEPModel.root");
                break;
            case 3://Unified NL Model
                sprintf(filenameResponse,"./NominalResponseMatrices/LinearBinning/NominalResponseUnifiedModel.root");
                break;
        }
    }

    
    TFile* ResponseF = TFile::Open(filenameResponse);
    
    for (Int_t AD =0; AD<ADsEH3; AD++)
    {
        NominalResponseMatrixH[AD] = (TH2F*)gDirectory->Get(Form("EvisEnu%i;2",AD));
    }
    
    ResponseF->Close();
}

void StatisticalCovarianceMatrix :: LoadBackgrounds()
{
    TFile* BackgroundsF = TFile::Open("./BackgroundSpectrum/Backgrounds.root");
    
    for (Int_t week = 0; week<Nweeks; week++)
    {
        for (Int_t AD =0; AD<NADs; AD++)
        {
            AccidentalsH[AD][week]= (TH1F*)gDirectory->Get(Form("Accidentals_AD%i",AD+1));
            LiHeH[AD][week]=(TH1F*)gDirectory->Get("LiHe");//Missing LiHe inputs so far in Hydrogen Analysis
            FastNeutronsH[AD][week]=(TH1F*)gDirectory->Get("FN");
            AmCH[AD][week]=(TH1F*)gDirectory->Get("AmC");
        }
    }
    BackgroundsF->Close();
}