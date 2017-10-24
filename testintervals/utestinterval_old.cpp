#include <limits.h>
#include "gtest/gtest.h"
#include "utestinterval.hpp"
 
int main(int argc, char **argv) 
{
 if(argc==1){
   std::cout << "usage: " << argv[0] << ' ' << "'path to json file'" << '\n';
   return 1;
 }
 JSONPATH = argv[1];
 ::testing::InitGoogleTest(&argc, argv);
 return RUN_ALL_TESTS();
}
