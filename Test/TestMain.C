#include "Test.h"


int main()
{
 
    std::cout << " TEST MAIN " << std::endl;
    gStyle->SetOptStat(1111111);
    Test* TestObject = new Test();
    TestObject->LoganCrossCheck();
    delete TestObject;
    return 0;
}
