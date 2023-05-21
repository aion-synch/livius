#include <iostream>
#include "colors.h"
#include "Tests/test_cases.h"

using namespace std;

int main(int argc, char* argv[])
{
    corenc::color::color_output("A finite element method library", corenc::color::CYAN);
    corenc::test_cases tests;
    tests.perform();
    return 0;
}
