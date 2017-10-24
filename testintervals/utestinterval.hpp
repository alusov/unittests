#ifndef UTESTINTERVAL_HPP
#define UTESTINTERVAL_HPP

#include <limits.h>
#include <algorithm>
#include "gtest/gtest.h"
#include "testfuncs/testfuncs.hpp"
#include "testfuncs/benchmarks.hpp"
#include "expression/expr.hpp"
#include "expression/algorithm.hpp"
#include "pointgen/randpointgen.hpp"
#include "box/box.hpp"

#define OFFSET 0.001



using namespace snowgoose::expression;
using namespace snowgoose;
using namespace OPTITEST;

class IntervalTest : public ::testing::Test {
 protected:

	IntervalTest() 
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

	std::vector<Interval<double>> getIntervals(const std::vector<double> &point)
	{
		std::vector<Interval<double>> intervals;
		std::for_each(point.begin(), point.end(), [&](double x) { intervals.push_back({ x - OFFSET, x + OFFSET }); });
		return intervals;
	}

	void TestInterval(const Benchmark<double> &bm)
	{
		auto point = getRandomPoint(bm);
		double funcValue = bm.calcFunc(point);
		auto intervals = getIntervals(point);
		auto interval = bm.calcInterval(intervals);

		double lowBound = interval.lb();
		double upperBound = interval.rb();

		ASSERT_GE(funcValue, lowBound);
		ASSERT_LE(funcValue, upperBound);
	}

};

TEST_F(IntervalTest, TestIntervalAckley1)
{
	int N = 3;
	TestInterval(Ackley1Benchmark<double>(N));
}

TEST_F(IntervalTest, TestIntervalAckley3)
{
	TestInterval(Ackley3Benchmark<double>());
}

TEST_F(IntervalTest, TestIntervalAckley4)
{
	TestInterval(Ackley4Benchmark<double>());
}

TEST_F(IntervalTest, TestIntervalAlpine2)
{
	TestInterval(Alpine2Benchmark<double>());
}


#endif
	 






