#include <stdio.h>
#include "CrossSection.h"
#include "ReactorSpectrumMultiple.h"
#include "AntineutrinoSpectrum.h"
#include "OscillationReactor.h"
#include "Oscillation.h"
#include "NominalData.h"

int main(void)
{
    TH1::AddDirectory(kFALSE);

        NominalData* Data;
        bool RandomSin22t13;
        bool RandomIsotopeFraction;
        bool RandomReactorPower;
        
        void MakePrediction();
        void SetRandomIsotopeFraction(bool);
        void SetRandomReactorPower(bool);
    
        //Randomize parameters: 1 to active, 0 to desactivate
        RandomSin22t13 = 0;
        RandomIsotopeFraction = 0;
        RandomReactorPower = 0;
    
    //Run one time with random parameters = 0 to generate the Expected Oscillation file.

    
        RandomReactorPower=0;
    
        RandomIsotopeFraction=1;
        std::cout << "Calculating Prediction" << std::endl;

        Data = new NominalData();
        ReactorSpectrumMultiple Reactor(Data);
        AntineutrinoSpectrum Antineutrino(Data); // Antineutrino Spectrum class, will multiply the cross section and the reactor spectrum.
        OscillationReactor OscRea(Data);
        Oscillation Osc(Data);
        
        Reactor.MultipleReactorSpectrumMain(RandomIsotopeFraction,RandomReactorPower,Data);
        Antineutrino.AntineutrinoSpectrumMain((RandomIsotopeFraction||RandomReactorPower));
        OscRea.OscillationFromReactorData((RandomSin22t13||RandomIsotopeFraction||RandomReactorPower));
        Osc.OscillationFromNearHallData((RandomSin22t13||RandomIsotopeFraction||RandomReactorPower));
        //Here I should add the distortion made by systematics
        
        OscRea.~OscillationReactor();
        Osc.~Oscillation();
        Antineutrino.~AntineutrinoSpectrum();
        Reactor.~ReactorSpectrumMultiple();

}