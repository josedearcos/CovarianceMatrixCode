void plot_statistical_matrix(){
    
    TFile* MatrixF = TFile::Open("./CovarianceMatrices/CovarianceMatrixMatrices.root");
        
        TH1F* MatrixH=(TH2D*)gDirectory->Get("Statistical Covariance Matrix");

        MatrixH->SetDirectory(0);
        MatrixH->Draw("colz");

    MatrixF->Close();
}