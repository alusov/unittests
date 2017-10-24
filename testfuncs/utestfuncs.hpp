#ifndef UTESTFUNCS_HPP
#define UTESTFUNCS_HPP

#include <limits.h>
#include "gtest/gtest.h"
#include "testfuncs/testfuncs.hpp"
#include "testfuncs/benchmarks.hpp"
#include "expression/expr.hpp"
#include "expression/algorithm.hpp"



#define EPSILON 0.001
#define OUTEXPR false


using namespace snowgoose::expression;
using namespace OPTITEST;

class FuncsTest : public ::testing::Test {
 protected:
	FuncsTest()
	{
	}
	void Test(const Benchmark<double> &bm, double eps = EPSILON)
	{
		std::vector<double> globMinX;
		for(int i=0; i < bm.getDim(); i++) 
			globMinX.push_back(bm.getGlobMinX(i));
		double globMinY = bm.calcFunc(globMinX);
		double expected = bm.getGlobMinY();
		double epsilon = eps;
		ASSERT_NEAR(expected, globMinY, epsilon);
	}
};
	 
TEST_F(FuncsTest, TestAckley1)
{
	int N = 3;
	Test(Ackley1Benchmark<double>(N));
}


TEST_F(FuncsTest, TestAckley3)
{
	Test(Ackley3Benchmark<double>());
}


TEST_F(FuncsTest, TestAckley4)
{
	Test(Ackley4Benchmark<double>());
}


TEST_F(FuncsTest, TestAlpine2)
{
    	Test(Alpine2Benchmark<double>());
}


#endif






