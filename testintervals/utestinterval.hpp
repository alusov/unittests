#ifndef UTESTINTERVAL_HPP
#define UTESTINTERVAL_HPP

#include <limits.h>
#include <algorithm>
#include "gtest/gtest.h"
#include "expression/testfuncs.hpp"
#include "expression/expr.hpp"
#include "expression/algorithm.hpp"
#include "descfunc/descfunc.hpp"
#include "descfunc/keys.hpp"
#include "pointgen/randpointgen.hpp"
#include "box/box.hpp"

#define OFFSET 0.001

const char* JSONPATH;

using namespace snowgoose::expression;
using namespace snowgoose;

class IntervalTest : public ::testing::Test {
 protected:

  IntervalTest() : dfr(JSONPATH)
  {
  }

  virtual void SetUp() 
  {
		
  }
 
  virtual void TearDown() 
  {
  }


	std::vector<double> getRandomPoint(const std::string& key, int customDim = 0)
	{
		auto desc = dfr.getdesr(key);
		const int dim = desc.anyDim ? customDim : desc.dim;
		Box<double> box(dim);
		
		for (int i = 0; i < dim; i++)
		{
			int boundIndex = desc.anyDim ? 0 : i;
			box.mA[i] = desc.bounds[boundIndex].first;
			box.mB[i] = desc.bounds[boundIndex].second;
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

	void TestInterval(const std::string& key,const Expr<double> &exprFunc,const Expr<Interval<double>> &exprInterval, int customDim = 0)
	{
		auto point = getRandomPoint(key, customDim);
		double funcValue = exprFunc.calc(point, FuncAlg<double>());
		auto intervals = getIntervals(point);
		auto interval = exprInterval.calc(intervals, InterEvalAlg<double>());

		double lowBound = interval.lb();
		double upperBound = interval.rb();

		ASSERT_GE(funcValue, lowBound);
		ASSERT_LE(funcValue, upperBound);
	}

	DescFuncReader dfr;
};

TEST_F(IntervalTest, TestIntervalAckley1)
{
	int N = 3;
	TestInterval(K.Ackley1, Ackley1<double>(N), Ackley1<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalAckley2)
{
	int N = 4;
	TestInterval(K.Ackley2, Ackley2<double>(N), Ackley2<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalAckley3)
{
	TestInterval(K.Ackley3, Ackley3<double>(), Ackley3<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalAckley4)
{
	TestInterval(K.Ackley4, Ackley4<double>(), Ackley4<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalAdjiman)
{
	TestInterval(K.Adjiman, Adjiman<double>(), Adjiman<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalAlpine1)
{
	int N = 3;
	TestInterval(K.Alpine1, Alpine1<double>(N), Alpine1<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalAlpine2)
{
	TestInterval(K.Alpine2, Alpine2<double>(), Alpine2<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBrad)
{
	TestInterval(K.Brad, Brad<double>(), Brad<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBartelsConn)
{
	TestInterval(K.BartelsConn, BartelsConn<double>(), BartelsConn<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBeale)
{
	TestInterval(K.Beale, Beale<double>(), Beale<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBiggsExpr2)
{
	TestInterval(K.BiggsEXP2, BiggsExpr2<double>(), BiggsExpr2<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBiggsExpr3)
{
	TestInterval(K.BiggsEXP3, BiggsExpr3<double>(), BiggsExpr3<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBiggsExpr4)
{
	TestInterval(K.BiggsEXP4, BiggsExpr4<double>(), BiggsExpr4<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBiggsExpr5)
{
	TestInterval(K.BiggsEXP5, BiggsExpr5<double>(), BiggsExpr5<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBiggsExpr6)
{
	TestInterval(K.BiggsEXP6, BiggsExpr6<double>(), BiggsExpr6<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBird)
{
	TestInterval(K.Bird, Bird<double>(), Bird<Interval<double>>());
}

#endif
	 






