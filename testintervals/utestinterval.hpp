#ifndef UTESTINTERVAL_HPP
#define UTESTINTERVAL_HPP

#include <limits.h>
#include <algorithm>
#include "gtest/gtest.h"
#include "testfuncs/manydim/testfuncs.hpp"
#include "testfuncs/manydim/benchmarks.hpp"
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
			box.mA[i] = bm.getBounds()[i].first;
			box.mB[i] = bm.getBounds()[i].second;
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

TEST_F(IntervalTest, IntervalTestAckley1)
{
	TestInterval(Ackley1Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestAckley2)
{
	TestInterval(Ackley2Benchmark<double>(4));
}

TEST_F(IntervalTest, IntervalTestAckley3)
{
	TestInterval(Ackley3Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestAdjiman)
{
	TestInterval(AdjimanBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestAlpine1)
{
	TestInterval(Alpine1Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestAlpine2)
{
    TestInterval(Alpine2Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBrad)
{
	TestInterval(BradBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBartelsConn)
{
	TestInterval(BartelsConnBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBeale)
{
	TestInterval(BealeBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBiggsbenchmark2)
{
	TestInterval(BiggsEXP2Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBiggsbenchmark3)
{
	TestInterval(BiggsEXP3Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBiggsbenchmark4)
{
	TestInterval(BiggsEXP4Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBiggsbenchmark5)
{
	TestInterval(BiggsEXP5Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBiggsbenchmark6)
{
	TestInterval(BiggsEXP6Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBird)
{
	TestInterval(BirdBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBohachevsky1)
{
        TestInterval(Bohachevsky1Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBohachevsky2)
{
        TestInterval(Bohachevsky2Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBohachevsky3)
{
        TestInterval(Bohachevsky3Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBooth)
{
        TestInterval(BoothBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBoxBettsQuadraticSum)
{
        TestInterval(BoxBettsQuadraticSumBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBraninRCOS)
{
        TestInterval(BraninRCOSBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBraninRCOS2)
{
        TestInterval(BraninRCOS2Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBrent)
{
        TestInterval(BrentBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBrown)
{
        
        TestInterval(BrownBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestBukin2)
{
        TestInterval(Bukin2Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBukin4)
{
        TestInterval(Bukin4Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestBukin6)
{
        TestInterval(Bukin6Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestCamelSixHump)
{
        TestInterval(CamelSixHumpBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestCamelThreeHump)
{
        TestInterval(CamelThreeHumpBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestChichinadze)
{
        TestInterval(ChichinadzeBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestChungReynolds)
{
        
        TestInterval(ChungReynoldsBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestColville)
{
        TestInterval(ColvilleBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestComplex)
{
        TestInterval(ComplexBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestCosineMixture)
{
	TestInterval(CosineMixtureBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestCrossInTray)
{
	TestInterval(CrossInTrayBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestCrossLeg)
{
	TestInterval(CrossLegBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestCube)
{
        TestInterval(CubeBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestDavis)
{
        TestInterval(DavisBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestDeb1)
{
        TestInterval(Deb1Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestDeckkersAarts)
{
        TestInterval(DeckkersAartsBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestDixonPrice)
{
        TestInterval(DixonPriceBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestDolan)
{
        TestInterval(DolanBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestDropWave)
{
        TestInterval(DropWaveBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestEasom)
{
        TestInterval(EasomBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestEggCrate)
{
        TestInterval(EggCrateBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestEggHolder)
{
        TestInterval(EggHolderBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestElAttarVidyasagarDutt)
{
        TestInterval(ElAttarVidyasagarDuttBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestEngvall)
{
        TestInterval(EngvallBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestExp2)
{
        TestInterval(Exp2Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestExponential)
{
        TestInterval(ExponentialBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestFreudensteinRoth)
{
        TestInterval(FreudensteinRothBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestGoldsteinPrice)
{
        TestInterval(GoldsteinPriceBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestGramacyLee2)
{
        TestInterval(GramacyLee2Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestGramacyLee3)
{
        TestInterval(GramacyLee3Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestGriewank)
{
        
        TestInterval(GriewankBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestHansen)
{
        TestInterval(HansenBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestHartman3)
{
        TestInterval(Hartman3Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestHartman6)
{
        TestInterval(Hartman6Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestHelicalValley)
{
        TestInterval(HelicalValleyBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestHimmelblau)
{
        TestInterval(HimmelblauBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestHosaki)
{
        TestInterval(HosakiBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestJennrichSampson)
{
        TestInterval(JennrichSampsonBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestKeane)
{
        TestInterval(KeaneBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestLangerman5)
{
        TestInterval(Langerman5Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestLeon)
{
        TestInterval(LeonBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestMatyas)
{
        TestInterval(MatyasBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestMcCormick)
{
        TestInterval(McCormickBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestMieleCantrell)
{
        TestInterval(MieleCantrellBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestMishra3)
{
        TestInterval(Mishra3Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestMishra4)
{
        TestInterval(Mishra4Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestMishra5)
{
        TestInterval(Mishra5Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestMishra6)
{
        TestInterval(Mishra6Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestMishra7)
{
        TestInterval(Mishra7Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestMishra8)
{
        TestInterval(Mishra8Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestMishra9)
{
        TestInterval(Mishra9Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestParsopoulos)
{
        TestInterval(ParsopoulosBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestPathological)
{       
        TestInterval(PathologicalBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestPeriodic)
{
        TestInterval(PeriodicBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestPinter)
{
        TestInterval(PinterBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestPowellSingular2)
{
        TestInterval(PowellSingular2Benchmark<double>(8));
}

TEST_F(IntervalTest, IntervalTestPowellSum)
{
        TestInterval(PowellSingular2Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestPrice1)
{
        TestInterval(Price1Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestPrice2)
{
        TestInterval(Price2Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestPrice3)
{
        TestInterval(Price3Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestPrice4)
{
        TestInterval(Price4Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestProblem02)
{
        TestInterval(Problem02Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestProblem04)
{
        TestInterval(Problem04Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestProblem05)
{
        TestInterval(Problem05Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestProblem06)
{
        TestInterval(Problem06Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestQing)
{
        TestInterval(QingBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestQuadratic)
{
        TestInterval(QuadraticBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestQuintic)
{
        TestInterval(QuinticBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestRosenbrock)
{
        TestInterval(RosenbrockBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestRosenbrockModified)
{
        TestInterval(RosenbrockModifiedBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestRotatedEllipse)
{
        TestInterval(RotatedEllipseBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestRotatedEllipse2)
{
        TestInterval(RotatedEllipse2Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestScahffer1)
{
        TestInterval(Scahffer1Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestScahffer3)
{
        TestInterval(Scahffer3Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestScahffer4)
{
        TestInterval(Scahffer4Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestScahffer2_6)
{
        TestInterval(Scahffer2_6Benchmark<double>());
}


TEST_F(IntervalTest, IntervalTestSchafferF6)
{
        
        TestInterval(SchafferF6Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestSchmidtVetters)
{
        TestInterval(SchmidtVettersBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestSchumerSteiglitz)
{
        TestInterval(SchumerSteiglitzBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestSchwefel)
{
        TestInterval(SchwefelBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestSchwefel1_2)                   
{

        TestInterval(Schwefel1_2Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestSchwefel2_20)
{
        TestInterval(Schwefel2_20Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestSchwefel2_22)
{
        TestInterval(Schwefel2_20Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestSchwefel2_23)
{
        TestInterval(Schwefel2_23Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestSchwefel2_26)
{
        TestInterval(Schwefel2_26Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestSchwefel2_36)
{
	TestInterval(Schwefel2_36Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestSchwefel2_4)
{
        TestInterval(Schwefel2_4Benchmark<double>(3));
}


TEST_F(IntervalTest, IntervalTestShekel10)
{
        TestInterval(Shekel10Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestShekel5)
{
        TestInterval(Shekel5Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestShekel7)
{
        TestInterval(Shekel7Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestShubert)
{
        TestInterval(ShubertBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestShubert2)
{
        TestInterval(Shubert2Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestShubert3)
{
        TestInterval(Shubert3Benchmark<double>());
}


TEST_F(IntervalTest, IntervalTestSolomon)
{
        TestInterval(SolomonBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestSphere)
{
        TestInterval(SphereBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestStrechedVSineWave)
{
        TestInterval(StrechedVSineWaveBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestStyblinskiTang)
{
        TestInterval(StyblinskiTangBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestSumSquares)
{
        TestInterval(SumSquaresBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestTable1HolderTable1)
{
        TestInterval(Table1HolderTable1Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestTable2HolderTable2)
{
        TestInterval(Table2HolderTable2Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestTable3Carrom)
{
        TestInterval(Table3CarromBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestTesttubeHolder)
{
        TestInterval(TesttubeHolderBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestTrecanni)
{
        TestInterval(TrecanniBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestTrefethen)
{
        TestInterval(TrefethenBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestTrid10)
{
	TestInterval(Trid10Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestTrid6)
{
	TestInterval(Trid6Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestTrigonometric1)
{
        TestInterval(Trigonometric1Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestTrigonometric2)
{
        TestInterval(Trigonometric2Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestTripod)
{
        TestInterval(TripodBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestUrsem1)
{
        TestInterval(Ursem1Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestUrsem3)
{
        TestInterval(Ursem3Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestUrsem4)
{
        TestInterval(Ursem4Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestUrsemWaves)
{
        TestInterval(UrsemWavesBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestVenterSobiezcczanskiSobieski)
{
        TestInterval(VenterSobiezcczanskiSobieskiBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestWWavy)
{
        
        TestInterval(WWavyBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestWayburnSeader1)
{
        TestInterval(WayburnSeader1Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestWayburnSeader2)
{
        TestInterval(WayburnSeader2Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestWayburnSeader3)
{
        TestInterval(WayburnSeader3Benchmark<double>());
}

TEST_F(IntervalTest, IntervalTestWeierstrass)
{
        
        TestInterval(WeierstrassBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestWhitley)
{
        
        TestInterval(WhitleyBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestWolfe)
{
	TestInterval(WolfeBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestXinSheYang2)
{
        
        TestInterval(XinSheYang2Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestXinSheYang3)
{
        
        TestInterval( XinSheYang3Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestXinSheYang4)
{
        
        TestInterval(XinSheYang4Benchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestZakharov)
{
        
        TestInterval(ZakharovBenchmark<double>(3));
}

TEST_F(IntervalTest, IntervalTestZettl)
{
        TestInterval(ZettlBenchmark<double>());
}

TEST_F(IntervalTest, IntervalTestZirilli)
{
        TestInterval(ZirilliBenchmark<double>());
}


#endif
	 






