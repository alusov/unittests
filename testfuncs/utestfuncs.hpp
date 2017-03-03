#ifndef UTESTFUNCS_HPP
#define UTESTFUNCS_HPP

#include <limits.h>
#include "gtest/gtest.h"
#include "expression/testfuncs.hpp"
#include "expression/expr.hpp"
#include "expression/algorithm.hpp"
#include "descfunc/descfunc.hpp"
#include "descfunc/keys.hpp"


#define EPSILON 0.0001
#define OUTEXPR true

const char* JSONPATH;

using namespace snowgoose::expression;

class FuncsTest : public ::testing::Test {
 protected:

  FuncsTest() : dfr(JSONPATH)
  {
  }
  
  virtual void SetUp() 
  {
		
  }
 
  virtual void TearDown() 
  {
  }

	void TestCustomDim(const std::string& key,const Expr<double>& expr, int dim, int indexGlobMin =0)
	{
		if (OUTEXPR) std::cout << '\n' << expr << '\n';
		auto desc = dfr.getdesr(key);
		std::vector<double> globMinX = std::vector<double>(dim, desc.globMinX[indexGlobMin][0]);
		double globMinY = expr.calc(globMinX, FuncAlg<double>());		
		double expected = desc.globMinY;
		double epsilon = EPSILON;
		ASSERT_NEAR(expected, globMinY, epsilon);
	}

	void Test(const std::string& key,const Expr<double>& expr, int indexGlobMin = 0)
	{
		if (OUTEXPR) std::cout << '\n' << expr << '\n';
		auto desc = dfr.getdesr(key);
		int dim = desc.dim;
		std::vector<double> globMinX = desc.globMinX[indexGlobMin];
		double globMinY = expr.calc(globMinX, FuncAlg<double>());
		double expected = desc.globMinY;
		double epsilon = EPSILON;
		ASSERT_NEAR(expected, globMinY, epsilon);
	}

	DescFuncReader dfr;
};
	 
TEST_F(FuncsTest, TestAckley1)
{
	int N = 3;
	auto expr = Ackley1<double>(N);
	TestCustomDim(K.Ackley1, expr, N);
}

TEST_F(FuncsTest, TestAckley2)
{
	int N = 4;
	auto expr = Ackley2<double>(N);
	TestCustomDim(K.Ackley2, expr, N);
}

TEST_F(FuncsTest, TestAckley3)
{
	Test(K.Ackley3, Ackley3<double>());
}

/*this test is for full interval library version only*/
/*
TEST_F(FuncsTest, TestAckley4)
{
	Test(K.Ackley4, Ackley4<double>());
}
*/

TEST_F(FuncsTest, TestAdjiman)
{
	Test(K.Adjiman, Adjiman<double>());
}

TEST_F(FuncsTest, TestAlpine1)
{
	int N = 3;
	auto expr = Alpine1<double>(N);
	TestCustomDim(K.Alpine1, expr, N);
}

TEST_F(FuncsTest, TestAlpine2)
{
    Test(K.Alpine2, Alpine2<double>());
}

TEST_F(FuncsTest, TestBrad)
{
	Test(K.Brad, Brad<double>());
}

TEST_F(FuncsTest, TestBartelsConn)
{
	Test(K.BartelsConn, BartelsConn<double>());
}

TEST_F(FuncsTest, TestBeale)
{
	Test(K.Beale, Beale<double>());
}

TEST_F(FuncsTest, TestBiggsExpr2)
{
	Test(K.BiggsEXP2, BiggsExpr2<double>());
}

TEST_F(FuncsTest, TestBiggsExpr3)
{
	Test(K.BiggsEXP3, BiggsExpr3<double>());
}

TEST_F(FuncsTest, TestBiggsExpr4)
{
	Test(K.BiggsEXP4, BiggsExpr4<double>());
}

TEST_F(FuncsTest, TestBiggsExpr5)
{
	Test(K.BiggsEXP5, BiggsExpr5<double>());
}

TEST_F(FuncsTest, TestBiggsExpr6)
{
	Test(K.BiggsEXP6, BiggsExpr6<double>());
}

TEST_F(FuncsTest, TestBird)
{
	Test(K.Bird, Bird<double>());
}

#endif






