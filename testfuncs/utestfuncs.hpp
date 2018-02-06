#ifndef UTESTFUNCS_HPP
#define UTESTFUNCS_HPP

#include "gtest/gtest.h"
#include "testfuncs/manydim/benchmarks.hpp"
#define EPSILON 0.001

class FuncsTest : public ::testing::Test {
 protected:
	FuncsTest()
	{
	}

	void Test(const Benchmark<double> &bm, double eps = EPSILON)
	{
		double globMinY = bm.calcFunc(bm.getGlobMinX());
		double expected = bm.getGlobMinY();
		double epsilon = eps;
		ASSERT_NEAR(expected, globMinY, epsilon);
	}

};

	 
TEST_F(FuncsTest, TestAckley1)
{
	Test(Ackley1Benchmark<double>(3));
}

TEST_F(FuncsTest, TestAckley2)
{
	Test(Ackley2Benchmark<double>(4));
}

TEST_F(FuncsTest, TestAckley3)
{
	Test(Ackley3Benchmark<double>());
}

TEST_F(FuncsTest, TestAdjiman)
{
	Test(AdjimanBenchmark<double>());
}

TEST_F(FuncsTest, TestAlpine1)
{
	Test(Alpine1Benchmark<double>(3));
}

TEST_F(FuncsTest, TestAlpine2)
{
    Test(Alpine2Benchmark<double>());
}

TEST_F(FuncsTest, TestBrad)
{
	Test(BradBenchmark<double>());
}

TEST_F(FuncsTest, TestBartelsConn)
{
	Test(BartelsConnBenchmark<double>());
}

TEST_F(FuncsTest, TestBeale)
{
	Test(BealeBenchmark<double>());
}

TEST_F(FuncsTest, TestBiggsbenchmark2)
{
	Test(BiggsEXP2Benchmark<double>());
}

TEST_F(FuncsTest, TestBiggsbenchmark3)
{
	Test(BiggsEXP3Benchmark<double>());
}

TEST_F(FuncsTest, TestBiggsbenchmark4)
{
	Test(BiggsEXP4Benchmark<double>());
}

TEST_F(FuncsTest, TestBiggsbenchmark5)
{
	Test(BiggsEXP5Benchmark<double>());
}

TEST_F(FuncsTest, TestBiggsbenchmark6)
{
	Test(BiggsEXP6Benchmark<double>());
}

TEST_F(FuncsTest, TestBird)
{
	Test(BirdBenchmark<double>());
}

TEST_F(FuncsTest, TestBohachevsky1)
{
        Test(Bohachevsky1Benchmark<double>());
}

TEST_F(FuncsTest, TestBohachevsky2)
{
        Test(Bohachevsky2Benchmark<double>());
}

TEST_F(FuncsTest, TestBohachevsky3)
{
        Test(Bohachevsky3Benchmark<double>());
}

TEST_F(FuncsTest, TestBooth)
{
        Test(BoothBenchmark<double>());
}

TEST_F(FuncsTest, TestBoxBettsQuadraticSum)
{
        Test(BoxBettsQuadraticSumBenchmark<double>());
}

TEST_F(FuncsTest, TestBraninRCOS)
{
        Test(BraninRCOSBenchmark<double>());
}

TEST_F(FuncsTest, TestBraninRCOS2)
{
        Test(BraninRCOS2Benchmark<double>());
}

TEST_F(FuncsTest, TestBrent)
{
        Test(BrentBenchmark<double>());
}

TEST_F(FuncsTest, TestBrown)
{
        
        Test(BrownBenchmark<double>(3));
}

TEST_F(FuncsTest, TestBukin2)
{
        Test(Bukin2Benchmark<double>());
}

TEST_F(FuncsTest, TestBukin4)
{
        Test(Bukin4Benchmark<double>());
}

TEST_F(FuncsTest, TestBukin6)
{
        Test(Bukin6Benchmark<double>());
}

TEST_F(FuncsTest, TestCamelSixHump)
{
        Test(CamelSixHumpBenchmark<double>());
}

TEST_F(FuncsTest, TestCamelThreeHump)
{
        Test(CamelThreeHumpBenchmark<double>());
}

TEST_F(FuncsTest, TestChichinadze)
{
        Test(ChichinadzeBenchmark<double>());
}

TEST_F(FuncsTest, TestChungReynolds)
{
        
        Test(ChungReynoldsBenchmark<double>(3));
}

TEST_F(FuncsTest, TestColville)
{
        Test(ColvilleBenchmark<double>());
}

TEST_F(FuncsTest, TestComplex)
{
        Test(ComplexBenchmark<double>());
}

TEST_F(FuncsTest, TestCosineMixture)
{
	Test(CosineMixtureBenchmark<double>());
}

TEST_F(FuncsTest, TestCrossInTray)
{
	Test(CrossInTrayBenchmark<double>());
}

TEST_F(FuncsTest, TestCrossLeg)
{
	Test(CrossLegBenchmark<double>());
}

TEST_F(FuncsTest, TestCube)
{
        Test(CubeBenchmark<double>());
}

TEST_F(FuncsTest, TestDavis)
{
        Test(DavisBenchmark<double>());
}

TEST_F(FuncsTest, TestDeb1)
{
        Test(Deb1Benchmark<double>(3));
}

TEST_F(FuncsTest, TestDeckkersAarts)
{
        Test(DeckkersAartsBenchmark<double>(), 0.1);
}

TEST_F(FuncsTest, TestDixonPrice)
{
        Test(DixonPriceBenchmark<double>());
}

TEST_F(FuncsTest, TestDolan)
{
        Test(DolanBenchmark<double>());
}

TEST_F(FuncsTest, TestDropWave)
{
        Test(DropWaveBenchmark<double>());
}

TEST_F(FuncsTest, TestEasom)
{
        Test(EasomBenchmark<double>());
}

TEST_F(FuncsTest, TestEggCrate)
{
        Test(EggCrateBenchmark<double>());
}

TEST_F(FuncsTest, TestEggHolder)
{
        Test(EggHolderBenchmark<double>());
}

TEST_F(FuncsTest, TestElAttarVidyasagarDutt)
{
        Test(ElAttarVidyasagarDuttBenchmark<double>());
}

TEST_F(FuncsTest, TestEngvall)
{
        Test(EngvallBenchmark<double>());
}

TEST_F(FuncsTest, TestExp2)
{
        Test(Exp2Benchmark<double>());
}

TEST_F(FuncsTest, TestExponential)
{
        Test(ExponentialBenchmark<double>(3));
}

TEST_F(FuncsTest, TestFreudensteinRoth)
{
        Test(FreudensteinRothBenchmark<double>());
}

TEST_F(FuncsTest, TestGoldsteinPrice)
{
        Test(GoldsteinPriceBenchmark<double>());
}

TEST_F(FuncsTest, TestGramacyLee2)
{
        Test(GramacyLee2Benchmark<double>());
}

TEST_F(FuncsTest, TestGramacyLee3)
{
        Test(GramacyLee3Benchmark<double>());
}

TEST_F(FuncsTest, TestGriewank)
{
        
        Test(GriewankBenchmark<double>(3));
}

TEST_F(FuncsTest, TestHansen)
{
        Test(HansenBenchmark<double>());
}

TEST_F(FuncsTest, TestHartman3)
{
        Test(Hartman3Benchmark<double>());
}

TEST_F(FuncsTest, TestHartman6)
{
        Test(Hartman6Benchmark<double>());
}

TEST_F(FuncsTest, TestHelicalValley)
{
        Test(HelicalValleyBenchmark<double>());
}

TEST_F(FuncsTest, TestHimmelblau)
{
        Test(HimmelblauBenchmark<double>());
}

TEST_F(FuncsTest, TestHosaki)
{
        Test(HosakiBenchmark<double>());
}

TEST_F(FuncsTest, TestJennrichSampson)
{
        Test(JennrichSampsonBenchmark<double>(), 0.01);
}

TEST_F(FuncsTest, TestKeane)
{
        Test(KeaneBenchmark<double>());
}

TEST_F(FuncsTest, TestLangerman5)
{
        Test(Langerman5Benchmark<double>());
}

TEST_F(FuncsTest, TestLeon)
{
        Test(LeonBenchmark<double>());
}

TEST_F(FuncsTest, TestMatyas)
{
        Test(MatyasBenchmark<double>());
}

TEST_F(FuncsTest, TestMcCormick)
{
        Test(McCormickBenchmark<double>());
}

TEST_F(FuncsTest, TestMieleCantrell)
{
        Test(MieleCantrellBenchmark<double>());
}

TEST_F(FuncsTest, TestMishra3)
{
        Test(Mishra3Benchmark<double>(), 001);
}

TEST_F(FuncsTest, TestMishra4)
{
        Test(Mishra4Benchmark<double>(), 0.01);
}

TEST_F(FuncsTest, TestMishra5)
{
        Test(Mishra5Benchmark<double>());
}

TEST_F(FuncsTest, TestMishra6)
{
        Test(Mishra6Benchmark<double>());
}

TEST_F(FuncsTest, TestMishra7)
{
        Test(Mishra7Benchmark<double>());
}

TEST_F(FuncsTest, TestMishra8)
{
        Test(Mishra8Benchmark<double>());
}

TEST_F(FuncsTest, TestMishra9)
{
        Test(Mishra9Benchmark<double>());
}

TEST_F(FuncsTest, TestParsopoulos)
{
        Test(ParsopoulosBenchmark<double>());
}

TEST_F(FuncsTest, TestPathological)
{       
        Test(PathologicalBenchmark<double>(3));
}

TEST_F(FuncsTest, TestPeriodic)
{
        Test(PeriodicBenchmark<double>());
}

TEST_F(FuncsTest, TestPinter)
{
        Test(PinterBenchmark<double>(3));
}

TEST_F(FuncsTest, TestPowellSingular2)
{
        Test(PowellSingular2Benchmark<double>(8));
}

TEST_F(FuncsTest, TestPowellSum)
{
        Test(PowellSingular2Benchmark<double>(3));
}

TEST_F(FuncsTest, TestPrice1)
{
        Test(Price1Benchmark<double>());
}

TEST_F(FuncsTest, TestPrice2)
{
        Test(Price2Benchmark<double>());
}

TEST_F(FuncsTest, TestPrice3)
{
        Test(Price3Benchmark<double>());
}

TEST_F(FuncsTest, TestPrice4)
{
        Test(Price4Benchmark<double>());
}

TEST_F(FuncsTest, TestProblem02)
{
        Test(Problem02Benchmark<double>());
}

TEST_F(FuncsTest, TestProblem04)
{
        Test(Problem04Benchmark<double>());
}

TEST_F(FuncsTest, TestProblem05)
{
        Test(Problem05Benchmark<double>());
}

TEST_F(FuncsTest, TestProblem06)
{
        Test(Problem06Benchmark<double>());
}

TEST_F(FuncsTest, TestQing)
{
        Test(QingBenchmark<double>());
}

TEST_F(FuncsTest, TestQuadratic)
{
        Test(QuadraticBenchmark<double>(), 0.1);
}

TEST_F(FuncsTest, TestQuintic)
{
        Test(QuinticBenchmark<double>(3));
}

TEST_F(FuncsTest, TestRosenbrock)
{
        Test(RosenbrockBenchmark<double>(3));
}

TEST_F(FuncsTest, TestRosenbrockModified)
{
        Test(RosenbrockModifiedBenchmark<double>());
}

TEST_F(FuncsTest, TestRotatedEllipse)
{
        Test(RotatedEllipseBenchmark<double>());
}

TEST_F(FuncsTest, TestRotatedEllipse2)
{
        Test(RotatedEllipse2Benchmark<double>());
}

TEST_F(FuncsTest, TestScahffer1)
{
        Test(Scahffer1Benchmark<double>());
}

TEST_F(FuncsTest, TestScahffer3)
{
        Test(Scahffer3Benchmark<double>());
}

TEST_F(FuncsTest, TestScahffer4)
{
        Test(Scahffer4Benchmark<double>());
}

TEST_F(FuncsTest, TestScahffer2_6)
{
        Test(Scahffer2_6Benchmark<double>());
}


TEST_F(FuncsTest, TestSchafferF6)
{
        
        Test(SchafferF6Benchmark<double>(3));
}

TEST_F(FuncsTest, TestSchmidtVetters)
{
        Test(SchmidtVettersBenchmark<double>(), 0.01);
}

TEST_F(FuncsTest, TestSchumerSteiglitz)
{
        Test(SchumerSteiglitzBenchmark<double>(3));
}

TEST_F(FuncsTest, TestSchwefel)
{
        Test(SchwefelBenchmark<double>(3));
}

TEST_F(FuncsTest, TestSchwefel1_2)                   
{

        Test(Schwefel1_2Benchmark<double>(3));
}

TEST_F(FuncsTest, TestSchwefel2_20)
{
        Test(Schwefel2_20Benchmark<double>(3));
}

TEST_F(FuncsTest, TestSchwefel2_22)
{
        Test(Schwefel2_20Benchmark<double>(3));
}

TEST_F(FuncsTest, TestSchwefel2_23)
{
        Test(Schwefel2_23Benchmark<double>(3));
}

TEST_F(FuncsTest, TestSchwefel2_26)
{
        Test(Schwefel2_26Benchmark<double>());
}

TEST_F(FuncsTest, TestSchwefel2_36)
{
	Test(Schwefel2_36Benchmark<double>());
}

TEST_F(FuncsTest, TestSchwefel2_4)
{
        Test(Schwefel2_4Benchmark<double>(3));
}


TEST_F(FuncsTest, TestShekel10)
{
        Test(Shekel10Benchmark<double>());
}

TEST_F(FuncsTest, TestShekel5)
{
        Test(Shekel5Benchmark<double>());
}

TEST_F(FuncsTest, TestShekel7)
{
        Test(Shekel7Benchmark<double>());
}

TEST_F(FuncsTest, TestShubert)
{
        Test(ShubertBenchmark<double>());
}

TEST_F(FuncsTest, TestShubert2)
{
        Test(Shubert2Benchmark<double>());
}

TEST_F(FuncsTest, TestShubert3)
{
        Test(Shubert3Benchmark<double>());
}


TEST_F(FuncsTest, TestSolomon)
{
        Test(SolomonBenchmark<double>());
}

TEST_F(FuncsTest, TestSphere)
{
        Test(SphereBenchmark<double>(3));
}

TEST_F(FuncsTest, TestStrechedVSineWave)
{
        Test(StrechedVSineWaveBenchmark<double>(3));
}

TEST_F(FuncsTest, TestStyblinskiTang)
{
        Test(StyblinskiTangBenchmark<double>());
}

TEST_F(FuncsTest, TestSumSquares)
{
        Test(SumSquaresBenchmark<double>(3));
}

TEST_F(FuncsTest, TestTable1HolderTable1)
{
        Test(Table1HolderTable1Benchmark<double>());
}

TEST_F(FuncsTest, TestTable2HolderTable2)
{
        Test(Table2HolderTable2Benchmark<double>());
}

TEST_F(FuncsTest, TestTable3Carrom)
{
        Test(Table3CarromBenchmark<double>());
}

TEST_F(FuncsTest, TestTesttubeHolder)
{
        Test(TesttubeHolderBenchmark<double>());
}

TEST_F(FuncsTest, TestTrecanni)
{
        Test(TrecanniBenchmark<double>());
}

TEST_F(FuncsTest, TestTrefethen)
{
        Test(TrefethenBenchmark<double>());
}

TEST_F(FuncsTest, TestTrid10)
{
	Test(Trid10Benchmark<double>());
}

TEST_F(FuncsTest, TestTrid6)
{
	Test(Trid6Benchmark<double>());
}

TEST_F(FuncsTest, TestTrigonometric1)
{
        Test(Trigonometric1Benchmark<double>(3));
}

TEST_F(FuncsTest, TestTrigonometric2)
{
        Test(Trigonometric2Benchmark<double>(3));
}

TEST_F(FuncsTest, TestTripod)
{
        Test(TripodBenchmark<double>());
}

TEST_F(FuncsTest, TestUrsem1)
{
        Test(Ursem1Benchmark<double>());
}

TEST_F(FuncsTest, TestUrsem3)
{
        Test(Ursem3Benchmark<double>());
}

TEST_F(FuncsTest, TestUrsem4)
{
        Test(Ursem4Benchmark<double>());
}

TEST_F(FuncsTest, TestUrsemWaves)
{
        Test(UrsemWavesBenchmark<double>());
}

TEST_F(FuncsTest, TestVenterSobiezcczanskiSobieski)
{
        Test(VenterSobiezcczanskiSobieskiBenchmark<double>());
}

TEST_F(FuncsTest, TestWWavy)
{
        
        Test(WWavyBenchmark<double>(3));
}

TEST_F(FuncsTest, TestWayburnSeader1)
{
        Test(WayburnSeader1Benchmark<double>());
}

TEST_F(FuncsTest, TestWayburnSeader2)
{
        Test(WayburnSeader2Benchmark<double>());
}

TEST_F(FuncsTest, TestWayburnSeader3)
{
        Test(WayburnSeader3Benchmark<double>());
}

TEST_F(FuncsTest, TestWeierstrass)
{
        
        Test(WeierstrassBenchmark<double>(3));
}

TEST_F(FuncsTest, TestWhitley)
{
        
        Test(WhitleyBenchmark<double>(3));
}

TEST_F(FuncsTest, TestWolfe)
{
	Test(WolfeBenchmark<double>());
}

TEST_F(FuncsTest, TestXinSheYang2)
{
        
        Test(XinSheYang2Benchmark<double>(3));
}

TEST_F(FuncsTest, TestXinSheYang3)
{
        
        Test( XinSheYang3Benchmark<double>(3));
}

TEST_F(FuncsTest, TestXinSheYang4)
{
        
        Test(XinSheYang4Benchmark<double>(3));
}

TEST_F(FuncsTest, TestZakharov)
{
        
        Test(ZakharovBenchmark<double>(3));
}

TEST_F(FuncsTest, TestZettl)
{
        Test(ZettlBenchmark<double>());
}

TEST_F(FuncsTest, TestZirilli)
{
        Test(ZirilliBenchmark<double>());
}


#endif






