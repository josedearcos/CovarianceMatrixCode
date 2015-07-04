void RunSmoothMatrix()
{
    std::cout << "Starting smoothing" << std::endl;
    
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);
    
    const Int_t NADs = 6;
    const Int_t volumesX = 2;//or 10
    const Int_t volumesY = 1;//or 10
    const Int_t TrueBins = 39;
    const Int_t VisBins = 34;
    
    //Show all the outputs from each respective systematic:
    
    TH2D* Histo[100][NADs][volumesX][volumesY];
    TH2D* SmoothedHisto[100][NADs][volumesX][volumesY];
    TH2D* ComparedHisto[100][NADs][volumesX][volumesY];

    Int_t num = 1;
    
    TString fname[100];
    
    const char *dirname;
    
    dirname= "../";
    char *ext= (Form("34_39.root"));
    
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
                for(Int_t AD=0; AD<NADs; AD++)
                {
                    for(Int_t x=0; x<volumesX; x++)
                    {
                        for(Int_t y=0; y<volumesY; y++)
                        {
                            Histo[num][AD][x][y] = (TH2D*)FileF->Get(Form("FineEvisEnu%i,Cell%d,%d",AD+1,x,y));
                        }
                    }
                }
                FileF->Close();
                num++;
            }
        }
    }
    
    num--;
    
    std::cout << num << " Files loaded and smoothed" << std::endl;

    //Save:
    
    TString newname;
    
    for(Int_t index = 1; index < num; index++)
    {
        newname = "Smoothed";

        newname.Append(fname[index]);
        
        std::cout << newname << std::endl;
        
        TFile* SaveFile = new TFile(newname,"recreate");
        
        for(Int_t AD=0; AD<NADs; AD++)
        {
            for(Int_t x=0; x<volumesX; x++)
            {
                for(Int_t y=0; y<volumesY; y++)
                {
                    SmoothedHisto[index][AD][x][y] = (TH2D*)Histo[index][AD][x][y]->Clone(Form("SmoothedEvisEnu%i,Cell%d,%d",AD+1,x,y));
                    SmoothedHisto[index][AD][x][y]->Smooth(1);
                    
                    SmoothedHisto[index][AD][x][y]->Write();
                }
            }
        }
        delete SaveFile;
    }
    
    std::cout << "Files saved" << std::endl;
    
    //Compare original and smoothed

    for(Int_t index = 1; index < num; index++)
    {
        newname = "Compared";
        
        newname.Append(fname[index]);
        
        std::cout << newname << std::endl;
        
        TFile* SaveFile = new TFile(newname,"recreate");
        
        for(Int_t AD=0; AD<NADs; AD++)
        {
            for(Int_t x=0; x<volumesX; x++)
            {
                for(Int_t y=0; y<volumesY; y++)
                {
                    ComparedHisto[index][AD][x][y] = (TH2D*)Histo[index][AD][x][y]->Clone(Form("ComparedEvisEnu%i,Cell%d,%d",AD+1,x,y));
                    ComparedHisto[index][AD][x][y]->Add(SmoothedHisto[index][AD][x][y],-1);

                    ComparedHisto[index][AD][x][y]->Write();
                }
            }
        }
        delete SaveFile;
    }

}