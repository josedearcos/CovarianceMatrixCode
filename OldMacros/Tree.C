

void Tree(){
    
    TFile f("./RootOutputs/FarSpectrumFraction.root", "update");

    TFile f("./RootOutputs/toySpectra.root", "recreate");

    Char_t branchname[1024];

    TTree *tree = new TTree("tree", "Toy Spectra");
    TH1F* PredH = 0;
    
    for(int idet=0;idet<Ndetectors;++idet)
    {
        sprintf(branchname,"SpectrumAD%i",idet+1);
        tree->Branch(branchname,"TH1F", h_ad[idet]);
    }
    
    for (Int_t i = 0; i < nentries; i++){
        new_v= gRandom->Gaus(0, 1);
        newBranch->Fill();
    }

}