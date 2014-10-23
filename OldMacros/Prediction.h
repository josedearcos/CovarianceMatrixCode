#pragma once
#include <stdio.h>
#include "CrossSection.h"
#include "ReactorSpectrumMultiple.h"
#include "AntineutrinoSpectrum.h"
#include "OscillationReactor.h"
#include "Oscillation.h"
#include "NominalData.h"

class Prediction
{
private:
    NominalData* Data;
    bool RandomSin22t12;
    bool RandomIsotopeFraction;
    bool RandomReactorPower;
    
    Double_t Sin22t13;

public:
    Prediction();
    Prediction(NominalData*);

    void MakePrediction();
    void SetRandomS22t12(bool);
    void SetRandomIsotopeFraction(bool);
    void SetRandomReactorPower(bool);
};

Prediction :: Prediction()
{
    Data = new NominalData();
    //Randomize parameters: 1 to active, 0 to desactivate
    RandomSin22t12 = 0;
    RandomIsotopeFraction = 0;
    RandomReactorPower = 0;
}
Prediction :: Prediction(NominalData* data)
{
    //Randomize parameters: 1 to active, 0 to desactivate
    RandomSin22t12 = 0;
    RandomIsotopeFraction = 0;
    RandomReactorPower = 0;
    Data = (NominalData*)data;
}

//Run one time with random parameters = 0 to generate the Expected Oscillation file.
void Prediction :: MakePrediction()
{
        std::cout << "Calculating Prediction" << std::endl;
    
        ReactorSpectrumMultiple* Reactor = new ReactorSpectrumMultiple(Data);
        AntineutrinoSpectrum* Antineutrino= new AntineutrinoSpectrum(Data);
        OscillationReactor* OscRea= new OscillationReactor(Data);
        Oscillation* Osc= new Oscillation(Data);
    
        Reactor->MultipleReactorSpectrumMain(RandomIsotopeFraction,RandomReactorPower);
        Antineutrino->AntineutrinoSpectrumMain(RandomIsotopeFraction,RandomReactorPower);
        OscRea->OscillationFromReactorData(RandomSin22t12,RandomIsotopeFraction,RandomReactorPower);
        Osc->OscillationFromNearHallData(RandomSin22t12,RandomIsotopeFraction,RandomReactorPower);
        //Here I should add the distortion made by systematics
    
        OscRea->~OscillationReactor();
        Osc->~Oscillation();
        Antineutrino->~AntineutrinoSpectrum();
        Reactor->~ReactorSpectrumMultiple();
}

void Prediction :: SetRandomS22t12(bool randomSin22t12)
{
     RandomSin22t12 = randomSin22t12;
}

void Prediction :: SetRandomReactorPower(bool random_reactor_power)
{
    RandomReactorPower=random_reactor_power;
}

void Prediction :: SetRandomIsotopeFraction(bool random_isotope_fraction)
{
    RandomIsotopeFraction=random_isotope_fraction;
}