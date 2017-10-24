#ifndef UTESTINTERVAL_HPP
#define UTESTINTERVAL_HPP

#include <limits.h>
#include <algorithm>
#include "gtest/gtest.h"
#include "testfuncs/testfuncs.hpp"
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
using namespace OPTITEST;

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
		double funcValue = exprFunc.calc(FuncAlg<double>(point));
		auto intervals = getIntervals(point);
		auto interval = exprInterval.calc(InterEvalAlg<double>(intervals));

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

//*****************************************



TEST_F(IntervalTest, TestIntervalBohachevsky1)
{
        TestInterval(K.Bohachevsky1, Bohachevsky1<double>(), Bohachevsky1<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBohachevsky2)
{
        TestInterval(K.Bohachevsky2, Bohachevsky2<double>(), Bohachevsky2<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBohachevsky3)
{
        TestInterval(K.Bohachevsky3, Bohachevsky3<double>(), Bohachevsky3<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBooth)
{
        TestInterval(K.Booth, Booth<double>(), Booth<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBoxBettsQuadraticSum)
{
        TestInterval(K.BoxBettsQuadraticSum, BoxBettsQuadraticSum<double>(), BoxBettsQuadraticSum<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBraninRCOS)
{
        TestInterval(K.BraninRCOS, BraninRCOS<double>(), BraninRCOS<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBraninRCOS2)
{
        TestInterval(K.BraninRCOS2, BraninRCOS2<double>(), BraninRCOS2<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBrent)
{
        TestInterval(K.Brent, Brent<double>(), Brent<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBrown)
{
        int N = 3;
        TestInterval(K.Brown, Brown<double>(N), Brown<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalBukin2)
{
        TestInterval(K.Bukin2, Bukin2<double>(), Bukin2<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBukin4)
{
        TestInterval(K.Bukin4, Bukin4<double>(), Bukin4<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalBukin6)
{
        TestInterval(K.Bukin6, Bukin6<double>(), Bukin6<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalCamelSixHump)
{
        TestInterval(K.CamelSixHump, CamelSixHump<double>(), CamelSixHump<Interval<double>>());
}


TEST_F(IntervalTest, TestIntervalCamelThreeHump)
{
        TestInterval(K.CamelThreeHump, CamelThreeHump<double>(), CamelThreeHump<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalChichinadze)
{
        TestInterval(K.Chichinadze, Chichinadze<double>(), Chichinadze<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalChungReynolds)
{
        int N = 3;
        TestInterval(K.ChungReynolds, ChungReynolds<double>(N), ChungReynolds<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalColville)
{
        TestInterval(K.Colville, Colville<double>(), Colville<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalComplex)
{
        TestInterval(K.Complex, Complex<double>(), Complex<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalCosineMixture)
{
	TestInterval(K.CosineMixture, CosineMixture<double>(), CosineMixture<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalCrossInTray)
{
	TestInterval(K.CrossInTray, CrossInTray<double>(), CrossInTray<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalCrossLeg)
{
	TestInterval(K.CrossLeg, CrossLeg<double>(), CrossLeg<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalCube)
{
        TestInterval(K.Cube, Cube<double>(), Cube<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalDavis)
{
        TestInterval(K.Davis, Davis<double>(), Davis<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalDeb1)
{
        int N = 3;
        TestInterval(K.Deb1, Deb1<double>(N), Deb1<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalDeckkersAarts)
{
        TestInterval(K.DeckkersAarts, DeckkersAarts<double>(), DeckkersAarts<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalDixonPrice)
{
        TestInterval(K.DixonPrice, DixonPrice<double>(), DixonPrice<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalDolan)
{
        TestInterval(K.Dolan, Dolan<double>(), Dolan<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalDropWave)
{
        TestInterval(K.DropWave, DropWave<double>(), DropWave<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalEasom)
{
        TestInterval(K.Easom, Easom<double>(), Easom<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalEggCrate)
{
        TestInterval(K.EggCrate, EggCrate<double>(), EggCrate<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalEggHolder)
{
        TestInterval(K.EggHolder, EggHolder<double>(), EggHolder<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalElAttarVidyasagarDutt)
{
        TestInterval(K.ElAttarVidyasagarDutt, ElAttarVidyasagarDutt<double>(), ElAttarVidyasagarDutt<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalEngvall)
{
        TestInterval(K.Engvall, Engvall<double>(), Engvall<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalExp2)
{
        TestInterval(K.Exp2, Exp2<double>(), Exp2<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalExponential)
{
        int N = 3;
        TestInterval(K.Exponential, Exponential<double>(N), Exponential<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalFreudensteinRoth)
{
        TestInterval(K.FreudensteinRoth, FreudensteinRoth<double>(), FreudensteinRoth<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalGoldsteinPrice)
{
        TestInterval(K.GoldsteinPrice, GoldsteinPrice<double>(), GoldsteinPrice<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalGramacyLee2)
{
        TestInterval(K.GramacyLee2, GramacyLee2<double>(), GramacyLee2<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalGramacyLee3)
{
        TestInterval(K.GramacyLee3, GramacyLee3<double>(), GramacyLee3<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalGriewank)
{
        int N = 3;
        TestInterval(K.Griewank, Griewank<double>(N), Griewank<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalHansen)
{
        TestInterval(K.Hansen, Hansen<double>(), Hansen<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalHartman3)
{
        TestInterval(K.Hartman3, Hartman3<double>(), Hartman3<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalHartman6)
{
        TestInterval(K.Hartman6, Hartman6<double>(), Hartman6<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalHelicalValley)
{
        TestInterval(K.HelicalValley, HelicalValley<double>(), HelicalValley<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalHimmelblau)
{
        TestInterval(K.Himmelblau, Himmelblau<double>(), Himmelblau<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalHosaki)
{
        TestInterval(K.Hosaki, Hosaki<double>(), Hosaki<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalJennrichSampson)
{
        TestInterval(K.JennrichSampson, JennrichSampson<double>(), JennrichSampson<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalKeane)
{
        TestInterval(K.Keane, Keane<double>(), Keane<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalLangerman5)
{
        TestInterval(K.Langerman5, Langerman5<double>(), Langerman5<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalLeon)
{
        TestInterval(K.Leon, Leon<double>(), Leon<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalMatyas)
{
        TestInterval(K.Matyas, Matyas<double>(), Matyas<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalMcCormick)
{
        TestInterval(K.McCormick, McCormick<double>(), McCormick<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalMieleCantrell)
{
        TestInterval(K.MieleCantrell, MieleCantrell<double>(), MieleCantrell<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalMishra3)
{
        TestInterval(K.Mishra3, Mishra3<double>(), Mishra3<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalMishra4)
{
        TestInterval(K.Mishra4, Mishra4<double>(), Mishra4<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalMishra5)
{
        TestInterval(K.Mishra5, Mishra5<double>(), Mishra5<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalMishra6)
{
        TestInterval(K.Mishra6, Mishra6<double>(), Mishra6<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalMishra7)
{
        TestInterval(K.Mishra7, Mishra7<double>(), Mishra7<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalMishra8)
{
        TestInterval(K.Mishra8, Mishra8<double>(), Mishra8<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalMishra9)
{
        TestInterval(K.Mishra9, Mishra9<double>(), Mishra9<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalParsopoulos)
{
        TestInterval(K.Parsopoulos, Parsopoulos<double>(), Parsopoulos<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalPathological)
{
        int N = 3;
        TestInterval(K.Pathological, Pathological<double>(N), Pathological<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalPeriodic)
{
        TestInterval(K.Periodic, Periodic<double>(), Periodic<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalPinter)
{
        int N = 3;
        TestInterval(K.Pinter, Pinter<double>(N), Pinter<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalPowellSingular2)
{
        int N = 8;
        TestInterval(K.PowellSingular2, PowellSingular2<double>(N), PowellSingular2<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalPowellSum)
{
        int N = 3;
        TestInterval(K.PowellSum, PowellSum<double>(N), PowellSum<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalPrice1)
{
        TestInterval(K.Price1, Price1<double>(), Price1<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalPrice2)
{
        TestInterval(K.Price2, Price2<double>(), Price2<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalPrice3)
{
        TestInterval(K.Price3, Price3<double>(), Price3<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalPrice4)
{
        TestInterval(K.Price4, Price4<double>(), Price4<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalProblem02)
{
        TestInterval(K.Problem02, Problem02<double>(), Problem02<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalProblem04)
{
        TestInterval(K.Problem04, Problem04<double>(), Problem04<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalProblem05)
{
        TestInterval(K.Problem05, Problem05<double>(), Problem05<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalProblem06)
{
        TestInterval(K.Problem06, Problem06<double>(), Problem06<Interval<double>>());
}


TEST_F(IntervalTest, TestIntervalQing)
{
        TestInterval(K.Qing, Qing<double>(), Qing<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalQuadratic)
{
        TestInterval(K.Quadratic, Quadratic<double>(), Quadratic<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalQuintic)
{
        int N = 3;
        TestInterval(K.Quintic, Quintic<double>(N), Quintic<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalRosenbrock)
{
        int N = 3;
        TestInterval(K.Rosenbrock, Rosenbrock<double>(N), Rosenbrock<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalRosenbrockModified)
{
        TestInterval(K.RosenbrockModified, RosenbrockModified<double>(), RosenbrockModified<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalRotatedEllipse)
{
        TestInterval(K.RotatedEllipse, RotatedEllipse<double>(), RotatedEllipse<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalRotatedEllipse2)
{
        TestInterval(K.RotatedEllipse2, RotatedEllipse2<double>(), RotatedEllipse2<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalScahffer1)
{
        TestInterval(K.Scahffer1, Scahffer1<double>(), Scahffer1<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalScahffer3)
{
        TestInterval(K.Scahffer3, Scahffer3<double>(), Scahffer3<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalScahffer4)
{
        TestInterval(K.Scahffer4, Scahffer4<double>(), Scahffer4<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalScahffer2_6)
{
        TestInterval(K.Scahffer2_6, Scahffer2_6<double>(), Scahffer2_6<Interval<double>>());
}


TEST_F(IntervalTest, TestIntervalSchafferF6)
{
        int N = 3;
        TestInterval(K.SchafferF6, SchafferF6<double>(N), SchafferF6<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalSchmidtVetters)
{
        TestInterval(K.SchmidtVetters, SchmidtVetters<double>(), SchmidtVetters<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalSchumerSteiglitz)
{
        int N = 3;
        TestInterval(K.SchumerSteiglitz, SchumerSteiglitz<double>(N), SchumerSteiglitz<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalSchwefel)
{
        int N = 3;
        TestInterval(K.Schwefel, Schwefel<double>(N, 1.0), Schwefel<Interval<double>>(N, 1.0), N);
}

TEST_F(IntervalTest, TestIntervalSchwefel1_2)
                      
{
        int N = 3;
        TestInterval(K.Schwefel1_2, Schwefel1_2<double>(N), Schwefel1_2<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalSchwefel2_20)
{
        int N = 3;
        TestInterval(K.Schwefel2_20, Schwefel2_20<double>(N), Schwefel2_20<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalSchwefel2_22)
{
        int N = 3;
        TestInterval(K.Schwefel2_22, Schwefel2_22<double>(N), Schwefel2_22<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalSchwefel2_23)
{
        int N = 3;
        TestInterval(K.Schwefel2_23, Schwefel2_23<double>(N), Schwefel2_23<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalSchwefel2_26)
{
        TestInterval(K.Schwefel2_26, Schwefel2_26<double>(), Schwefel2_26<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalSchwefel2_36)
{
	TestInterval(K.Schwefel2_36, Schwefel2_36<double>(), Schwefel2_36<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalSchwefel2_4)
{
        int N = 3;
        TestInterval(K.Schwefel2_4, Schwefel2_4<double>(N), Schwefel2_4<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalShekel10)
{
        TestInterval(K.Shekel10, Shekel10<double>(), Shekel10<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalShekel5)
{
        TestInterval(K.Shekel5, Shekel5<double>(), Shekel5<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalShekel7)
{
        TestInterval(K.Shekel7, Shekel7<double>(), Shekel7<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalShubert)
{
        TestInterval(K.Shubert, Shubert<double>(), Shubert<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalShubert2)
{
        TestInterval(K.Shubert2, Shubert2<double>(), Shubert2<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalShubert3)
{
        TestInterval(K.Shubert3, Shubert3<double>(), Shubert3<Interval<double>>());
}


TEST_F(IntervalTest, TestIntervalSolomon)
{
        TestInterval(K.Solomon, Solomon<double>(), Solomon<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalSphere)
{
        int N = 3;
        TestInterval(K.Sphere, Sphere<double>(N), Sphere<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalStrechedVSineWave)
{
        int N = 3;
        TestInterval(K.StrechedVSineWave, StrechedVSineWave<double>(N), StrechedVSineWave<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalStyblinskiTang)
{
        TestInterval(K.StyblinskiTang, StyblinskiTang<double>(), StyblinskiTang<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalSumSquares)
{
        int N = 3;
        TestInterval(K.SumSquares, SumSquares<double>(N), SumSquares<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalTable1HolderTable1)
{
        TestInterval(K.Table1HolderTable1, Table1HolderTable1<double>(), Table1HolderTable1<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalTable2HolderTable2)
{
        TestInterval(K.Table2HolderTable2, Table2HolderTable2<double>(), Table2HolderTable2<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalTable3Carrom)
{
        TestInterval(K.Table3Carrom, Table3Carrom<double>(), Table3Carrom<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalTesttubeHolder)
{
        TestInterval(K.TesttubeHolder, TesttubeHolder<double>(), TesttubeHolder<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalTrecanni)
{
        TestInterval(K.Trecanni, Trecanni<double>(), Trecanni<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalTrefethen)
{
        TestInterval(K.Trefethen, Trefethen<double>(), Trefethen<Interval<double>>());
}


TEST_F(IntervalTest, TestIntervalTrid10)
{
	TestInterval(K.Trid10, Trid10<double>(), Trid10<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalTrid6)
{
	TestInterval(K.Trid6, Trid6<double>(), Trid6<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalTrigonometric1)
{
        int N = 3;
        TestInterval(K.Trigonometric1, Trigonometric1<double>(N), Trigonometric1<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalTrigonometric2)
{
        int N = 3;
        TestInterval(K.Trigonometric2, Trigonometric2<double>(N), Trigonometric2<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalTripod)
{
        TestInterval(K.Tripod, Tripod<double>(), Tripod<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalUrsem1)
{
        TestInterval(K.Ursem1, Ursem1<double>(), Ursem1<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalUrsem3)
{
        TestInterval(K.Ursem3, Ursem3<double>(), Ursem3<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalUrsem4)
{
        TestInterval(K.Ursem4, Ursem4<double>(), Ursem4<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalUrsemWaves)
{
        TestInterval(K.UrsemWaves, UrsemWaves<double>(), UrsemWaves<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalVenterSobiezcczanskiSobieski)
{
        TestInterval(K.VenterSobiezcczanskiSobieski, VenterSobiezcczanskiSobieski<double>(), VenterSobiezcczanskiSobieski<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalWWavy)
{
        int N = 3;
        TestInterval(K.WWavy, WWavy<double>(N), WWavy<Interval<double>>(N), N);
}



TEST_F(IntervalTest, TestIntervalWatson)
{
        TestInterval(K.Watson, Watson<double>(), Watson<Interval<double>>());
}


TEST_F(IntervalTest, TestIntervalWayburnSeader1)
{
        TestInterval(K.WayburnSeader1, WayburnSeader1<double>(), WayburnSeader1<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalWayburnSeader2)
{
        TestInterval(K.WayburnSeader2, WayburnSeader2<double>(), WayburnSeader2<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalWayburnSeader3)
{
        TestInterval(K.WayburnSeader3, WayburnSeader3<double>(), WayburnSeader3<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalWeierstrass)
{
        int N = 3;
        TestInterval(K.Weierstrass, Weierstrass<double>(N), Weierstrass<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalWhitley)
{
        int N = 3;
        TestInterval(K.Whitley, Whitley<double>(N), Whitley<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalWolfe)
{
	TestInterval(K.Wolfe, Wolfe<double>(), Wolfe<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalXinSheYang2)
{
        int N = 3;
        TestInterval(K.XinSheYang2, XinSheYang2<double>(N), XinSheYang2<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalXinSheYang3)
{
        int N = 3;
        TestInterval(K.XinSheYang3, XinSheYang3<double>(N), XinSheYang3<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalXinSheYang4)
{
        int N = 3;
        TestInterval(K.XinSheYang4, XinSheYang4<double>(N), XinSheYang4<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalZakharov)
{
        int N = 3;
        TestInterval(K.Zakharov, Zakharov<double>(N), Zakharov<Interval<double>>(N), N);
}

TEST_F(IntervalTest, TestIntervalZettl)
{
        TestInterval(K.Zettl, Zettl<double>(), Zettl<Interval<double>>());
}

TEST_F(IntervalTest, TestIntervalZirilli)
{
        TestInterval(K.Zirilli, Zirilli<double>(), Zirilli<Interval<double>>());
}

#endif
	 






