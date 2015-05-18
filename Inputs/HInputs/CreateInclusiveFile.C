#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include <sstream>
#include <string>
#include "TString.h"
using namespace std;

int CreateInclusiveFile(){
    
    Int_t VolumeID = 0;
    Int_t index = 0;

    string FileS = "Empty";
    
    while(index<2)
    {
        if(index==0)
        {
            FileS = "GdLS";
            VolumeID = 0;
        }
        else
        {
            FileS = "LS";
            VolumeID = 1;
        }
        
        std::cout << " Creating " << FileS << " Inclusive File " << std::endl;
        
        index++;
        
        const Int_t NADs =8;
        const Int_t DataStored = 17;
        Int_t NewPeriods = 101;
        //Read Inclusive File
        std::string line;
        Int_t linenum=0;//<---caution: only increments for lines that do not begin with #
        Int_t week=0;
        Double_t readvals[DataStored][NADs]={0};
        Int_t Nweeks;
        
        ifstream mainfile(("./nH_"+FileS+"_table.txt").c_str());
        while(!mainfile.eof())
        {
            std::getline(mainfile,line);
            std::string firstchar = line.substr(0,1);
            
            if(firstchar=="#") continue;//<-- ignore lines with comments
            
            //Special numbers
            if(linenum == 0)
            {
                Nweeks=atoi(line.c_str());
                //std::cout << "The number of periods is " << Nweeks << std::endl;
            }
            
            if(linenum == 1){}
            //        {
            //            if(atoi(firstchar.c_str())==0) isMC=true;
            //            if(atoi(firstchar.c_str())==1) isMC=false;
            //            cout << "Simflag: " << isMC << " (1-->MC, 0-->Data)" << endl;
            //        }
            
            if(linenum > 2)
            {
                //std::cout << "reading " << line << std::endl;
                std::istringstream iss(line);
                Int_t row=0;
                Int_t column=0;
                while(iss)
                {
                    std::string sub; iss >> sub;
                    if(column==0) week=atoi(sub.c_str());
                    
                    if(column==1) row=atoi(sub.c_str());
                    
                    //                if(data==8||data==9||data==11||data==13||data==15||data==17)
                    //                {
                    //                    if(column>1 && sub!="")
                    //                    {
                    //                        readvals[row][column-2]+=(atof(sub.c_str())*(atof(sub.c_str()));//quadratic sum of errors?
                    //                    }
                    //                }
                    //                else
                    //                {
                    if(column>1 && sub!="") readvals[row][column-2]+=atof(sub.c_str());
                    //                }
                    column+=1;
                }//looping over columns
            }
            
            linenum++;//only lines >2
        }
        
        
        //Create a file with evenly distributed data, it is fake, need real data:
        
        FILE *f = fopen(("./nH_"+FileS+"_table_Inclusive.txt").c_str(), "w");
        
        if (f == NULL)
        {
            printf("Error opening file!\n");
            exit(1);
        }
        
        fprintf(f,"# Total number of time periods\n\
1\n\
# Volume 0=GD-LS, 1=LS Volume\n\
%i\n\
# Delta M^2_{32} and uncertainty  (in eV^2)\n\
0.00243  0.00013\n\
# ===============================================================\n\
# Note, AD1/2 for EH1AD1/2, AD3/4 for EH2AD1/2, AD5/6/7/8 for EH3AD1/2/3/4\n\
# First column is time period number, second column is row number\n\
# The nominal period for each set of entries is one week\n\
# Row 0 ==>  Start UTC  |  End UTC | Start date and time\n\
# Row 1 ==>  Observed number of events in AD1 to AD8\n\
# Row 2 ==>  Live time in days for AD1 to AD8\n\
# Row 3 ==>  Muon veto efficiency for AD1 to AD8\n\
# Row 4 ==>  Multiplicity cut efficiency for AD1 to AD8\n\
# Row 5 ==>  Target Mass (kg) for AD1 to AD8\n\
# Row 6 ==>  Hydrogen capture fraction (%) for AD1 to AD8\n\
# Row 7 ==>  Efficiency of Prompt energy, Delayed energy, time, distance\n\
# Row 8 ==>  Uncorrelated detector efficiency uncertainty  (relative error in %)\n\
# Row 9 ==>  Uncorrelated reactor power uncertainty (relative error in %)\n\
#### Note: The background estimates are corrected for muon veto efficiency and multiplicity efficiency\n\
# Row 10 ==>  Expected number of accidental events per live day (AD1 to AD8)\n\
# Row 11 ==>  Absolute uncertainty on the accidental bkg per live day\n\
# Row 12 ==>  Expected number of fast-n events per live day (AD1 to AD8)\n\
# Row 13 ==>  Absolute uncertainty on the fast-n bkg per live day\n\
# Row 14 ==>  Expected number of li9/he8 events per live day (AD1 to AD8)\n\
# Row 15 ==>  Absolute uncertainty on the li9/he8 bkg per live day\n\
# Row 16 ==>  Expected number of AmC-corr events per live day (AD1 to AD8)\n\
# Row 17 ==>  Absolute uncertainty on the AmC-corr bkg per live day\n\
# ===============================================================\n"
,VolumeID);
        
        
        for(Int_t data = 0; data <= DataStored; data++)
        {
            fprintf(f, "%d %d ", 1, data);
            
            for(Int_t AD = 0; AD < NADs; AD++)
            {
                if(data==0)
                {
                    if(AD<3)//To avoid date wrong input.
                    {
                        fprintf(f,"%f ",readvals[data][AD]);
                    }
                }
                else if(data==3||data==4||data==5||data==6||data==7||data==8||data==9||data==11||data==13||data==15||data==17)
                {
                    fprintf(f,"%f ",readvals[data][AD]/double(NewPeriods));//Divide Efficiencies, sum rates
                }
                else
                {
                    fprintf(f,"%f ",readvals[data][AD]);
                }
            }
            
            fprintf(f,"\n");
        }
        
        
        fclose(f);
    }
    return 0;
}