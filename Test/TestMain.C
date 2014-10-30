#include "Test.h"


int main()
{
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);
    
    std::cout << " TEST MAIN " << std::endl;
    gStyle->SetOptStat(1111111);
    Test* TestObject = new Test();
    TestObject->TestCovarianceMatrixSamples();
    delete TestObject;
    return 0;
}
