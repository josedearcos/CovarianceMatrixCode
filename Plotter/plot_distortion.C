#include "TRandom3.h"

void plot_distortion(){
    
    TRandom3* rand = new TRandom3();

//    TF1 *func = new TF1("func","TMath::Abs([0]+[1]*x)",0,12);
//
//    rand->SetSeed(0);
//    Double_t slope=0.2*rand->Gaus(0,1);
//    Double_t anchor_point=3.5;
//    //want offset to be set by requiring func at anchor point to be 1
//    Double_t offset=(1-slope*anchor_point);
//    func->SetParameter(0,offset);
//    func->SetParameter(1,slope);
    
    //FN distortion
    TF1 *func = new TF1("func","[0]/(0.2*pow(x,0.1))+[1]",0,12);
    rand->SetSeed(0);
    Double_t scaling =0.2*rand->Gaus(0,1);
    func->SetParameter(0,scaling);
    func->SetParameter(1,0);
    //set offset so that func(FinalEnergy)=1;
    func->SetParameter(1,1-1*func->Eval(12));
    
    func->Draw();
    
    
   }
