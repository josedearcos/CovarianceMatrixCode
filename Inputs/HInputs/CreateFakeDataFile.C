#include <iostream>
#include <fstream>
#include "stdio.h"
#include <stdlib.h>
#include <sstream>
#include <string>
#include "TString.h"
using namespace std;

int CreateFakeDataFile(){
    
    const Int_t NADs =6;
    const Int_t DataStored = 20;
    Int_t NewPeriods = 101;
    //Read Inclusive File
    std::string line;
    Int_t linenum=0;//<---caution: only increments for lines that do not begin with #
    Int_t week=0;
    Double_t readvals[20][6]={0};
    Int_t Nweeks;
    
    ifstream mainfile("./P12E_Inclusive.txt");
    while(!mainfile.eof())
    {
        std::getline(mainfile,line);
        std::string firstchar = line.substr(0,1);
        
        if(firstchar=="#") continue;//<-- ignore lines with comments
        
        //Special numbers
        if(linenum == 0)
        {
            Nweeks=atoi(line.c_str());
            std::cout << "The number of periods is " << Nweeks << std::endl;
        }
        
        if(linenum == 1){}
        //        {
        //            if(atoi(firstchar.c_str())==0) isMC=true;
        //            if(atoi(firstchar.c_str())==1) isMC=false;
        //            cout << "Simflag: " << isMC << " (1-->MC, 0-->Data)" << endl;
        //        }
        
        if(linenum > 2)
        {
            std::cout << "reading " << line << std::endl;
            std::istringstream iss(line);
            Int_t row=0;
            Int_t column=0;
            while(iss)
            {
                std::string sub; iss >> sub;
                if(column==0) week=atoi(sub.c_str());
                
                if(column==1) row=atoi(sub.c_str());
                
                if(column>1 && sub!="") readvals[row][column-2]=atof(sub.c_str());
                column+=1;
            }//looping over columns
        }
        
        linenum++;//only lines >2
    }
    
    //Create a file with evenly distributed data, it is fake, need real data:
    
    FILE *f = fopen(Form("P12E_%d.txt",NewPeriods), "w");

    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    
    fprintf(f,"# Total number of time periods, FAKE FILE, NEED A REAL ONE\n");
    fprintf(f, "%d\n", NewPeriods);
    fprintf(f,"# DataFlag 1=Data, 0=MC\n");
    fprintf(f, "1\n");
    fprintf(f,"# Delta M^2_{32} and uncertainty  (in eV^2)\n");
    fprintf(f, "0.00243 0.00013 \n");
    fprintf(f, "# ===============================================================\n\
# First column is time period number, second column is row number\n\
# The nominal period for each set of entries is one week\n\
#\n\
# Row 0 ==>  Start UTC  |  End UTC | Start date and time\n\
# Row 1 ==>  Observed number of events in AD1 to AD6\n\
# Row 2 ==>  Live time in days for AD1 to AD6\n\
# Row 3 ==>  Muon veto efficiency for AD1 to AD6\n\
# Row 4 ==>  Multiplicity cut efficiency for AD1 to AD6\n\
# Row 5 ==>  6 MeV neutron cut efficiency for AD1 to AD6\n\
# Row 6 ==>  Uncorrelated reactor power uncertainty (relative error in %)\n\
# Row 7 ==>  Total IBD cut efficiency uncertainty (relative error in %)\n\
# Row 8 ==>  Total target mass for AD1 to AD6 (kg)\n\
#### Note: The background estimates are *not* corrected for efficiencies anymore\n\
# Row 9  ==>  Expected number of bkg events per live day (AD1 to AD6)\n\
# Row 10 ==>  Absolute uncertainty on the bkg estimates per live day\n\
# Row 11 ==>  Expected number of accidental events per live day (AD1 to AD6)\n\
# Row 12 ==>  Absolute uncertainty on the accidental bkg per live day\n\
# Row 13 ==>  Expected number of li9/he8 events per live day (AD1 to AD6)\n\
# Row 14 ==>  Absolute uncertainty on the li9/he8 bkg per live day\n\
# Row 15 ==>  Expected number of fast-n events per live day (AD1 to AD6)\n\
# Row 16 ==>  Absolute uncertainty on the fast-n bkg per live day\n\
# Row 17 ==>  Expected number of AmC-corr events per live day (AD1 to AD6)\n\
# Row 18 ==>  Absolute uncertainty on the AmC-corr bkg per live day\n\
# Row 19 ==>  Expected number of alpha-n events per live day (AD1 to AD6)\n\
# Row 20 ==>  Absolute uncertainty on the alpha-n bkg per live day)\n\
# ===============================================================\n"
    );
    
    for(Int_t period = 1; period <=NewPeriods; period++)
    {
        for(Int_t data = 0; data <= DataStored; data++)
        {
            fprintf(f, "%d %d ", period, data);//An average of the inclusive data for the expected number of events

            for(Int_t AD = 0; AD < NADs; AD++)
            {`
                if(data==0)
                {
                    if(AD<3)//To avoid date wrong input.
                    {
                        fprintf(f,"%f ",readvals[data][AD]);
                    }
                }
                else if(data==1||data==2||data==9||data==11||data==13||data==15||data==17||data==19)
                {
                    fprintf(f,"%f ",readvals[data][AD]/double(NewPeriods));
                }
                else
                {
                    fprintf(f,"%f ",readvals[data][AD]);
                }
            }
            
            fprintf(f,"\n");
        }
    }
    
    fclose(f);
    
    return 0;
}