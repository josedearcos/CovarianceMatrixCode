#include "nHToyMcTest.h"

int main()
{
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);

    nHToyMcTest* Test = new nHToyMcTest();
    Test->nHToyMcTestMain();
    return 0;
}