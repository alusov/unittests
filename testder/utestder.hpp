#ifndef UTESTFUNCS_HPP
#define UTESTFUNCS_HPP

#include <limits.h>
#include "gtest/gtest.h"
#include "testfuncs/testfuncs.hpp"
#include "testfuncs/benchmarks.hpp"
#include "expression/expr.hpp"
#include "expression/algder.hpp"
#include "derivatives/valder.hpp"
#include "pointgen/randpointgen.hpp"
#include "box/box.hpp"


#define EPSILON 0.1
#define PERCENT 1
#define DELTA 0.0001

using namespace snowgoose;
using namespace snowgoose::expression;
using namespace snowgoose::derivative;
using namespace OPTITEST;


class DerivativeTest : public ::testing::Test {
 protected:
	DerivativeTest()
	{
	}
    
	std::vector<double> getRandomPoint(const Benchmark<double> &bm)
	{
		const int dim = bm.getDim();
		Box<double> box(dim);
	
		for (int i = 0; i < dim; i++)
		{
			box.mA[i] = bm.getBounds(i).first;
			box.mB[i] = bm.getBounds(i).second;
		}
		RandomPointGenerator<double> rg(box);
		std::vector<double> point(dim, 0.0);
		rg.getPoint(point.data());
		return point;
	}
        
	void TestDerivative(const Benchmark<double> &bm)
	{
		std::vector<double> point = getRandomPoint(bm);
		auto der = bm.calcGrad(point);
		double func_val = bm.calcFunc(point);
		double func_val_by_der = der.value();
		double epsilon = EPSILON;
		ASSERT_NEAR(func_val, func_val_by_der, epsilon);
		
		auto grad = der.grad();
			
		for(int i=0; i < point.size(); ++i)
		{
		    auto new_point = std::vector<double>(point);
		    new_point[i]+=DELTA;
		    double new_func_val = bm.calcFunc(new_point);
		    double partial_derivative = (new_func_val - func_val)/DELTA;    
		    std::cout << "func_val=" << func_val << " new_func_val=" << new_func_val << '\n';
		    double expected = 100;
		    double derivative_relative_value = (partial_derivative/grad[i])*100;
		    double persent = PERCENT;
		    ASSERT_NEAR(expected, derivative_relative_value, persent);
		}
	}
};



TEST_F(DerivativeTest, TestDerivativeAckley1)
{
	int N = 3;
	TestDerivative(Ackley1Benchmark<double>(N));
}

TEST_F(DerivativeTest, TestDerivativeAckley3)
{
	TestDerivative(Ackley3Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeAckley4)
{
	TestDerivative(Ackley4Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeAlpine2)
{
	TestDerivative(Alpine2Benchmark<double>());
}


#endif






