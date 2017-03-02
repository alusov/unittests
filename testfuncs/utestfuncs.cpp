#include <limits.h>
#include "gtest/gtest.h"
#define FUNCDESCR "../funcdesc.json"
#include "utestfuncs.hpp"
#include <iostream>


 
int main(int argc, char **argv) 
{
 ::testing::InitGoogleTest(&argc, argv);
 return RUN_ALL_TESTS();
}
