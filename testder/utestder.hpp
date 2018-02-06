#ifndef UTESTFUNCS_HPP
#define UTESTFUNCS_HPP

#include <limits.h>
#include "gtest/gtest.h"
#include "testfuncs/manydim/testfuncs.hpp"
#include "testfuncs/manydim/benchmarks.hpp"
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
			box.mA[i] = bm.getBounds()[i].first;
			box.mB[i] = bm.getBounds()[i].second;
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
	TestDerivative(Ackley1Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeAckley2)
{
	TestDerivative(Ackley2Benchmark<double>(4));
}

TEST_F(DerivativeTest, TestDerivativeAckley3)
{
	TestDerivative(Ackley3Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeAdjiman)
{
	TestDerivative(AdjimanBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeAlpine1)
{
	TestDerivative(Alpine1Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeAlpine2)
{
    TestDerivative(Alpine2Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBrad)
{
	TestDerivative(BradBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBartelsConn)
{
	TestDerivative(BartelsConnBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBeale)
{
	TestDerivative(BealeBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBiggsbenchmark2)
{
	TestDerivative(BiggsEXP2Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBiggsbenchmark3)
{
	TestDerivative(BiggsEXP3Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBiggsbenchmark4)
{
	TestDerivative(BiggsEXP4Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBiggsbenchmark5)
{
	TestDerivative(BiggsEXP5Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBiggsbenchmark6)
{
	TestDerivative(BiggsEXP6Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBird)
{
	TestDerivative(BirdBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBohachevsky1)
{
        TestDerivative(Bohachevsky1Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBohachevsky2)
{
        TestDerivative(Bohachevsky2Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBohachevsky3)
{
        TestDerivative(Bohachevsky3Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBooth)
{
        TestDerivative(BoothBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBoxBettsQuadraticSum)
{
        TestDerivative(BoxBettsQuadraticSumBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBraninRCOS)
{
        TestDerivative(BraninRCOSBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBraninRCOS2)
{
        TestDerivative(BraninRCOS2Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBrent)
{
        TestDerivative(BrentBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBrown)
{
        
        TestDerivative(BrownBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeBukin2)
{
        TestDerivative(Bukin2Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBukin4)
{
        TestDerivative(Bukin4Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeBukin6)
{
        TestDerivative(Bukin6Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeCamelSixHump)
{
        TestDerivative(CamelSixHumpBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeCamelThreeHump)
{
        TestDerivative(CamelThreeHumpBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeChichinadze)
{
        TestDerivative(ChichinadzeBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeChungReynolds)
{
        
        TestDerivative(ChungReynoldsBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeColville)
{
        TestDerivative(ColvilleBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeComplex)
{
        TestDerivative(ComplexBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeCosineMixture)
{
	TestDerivative(CosineMixtureBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeCrossInTray)
{
	TestDerivative(CrossInTrayBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeCrossLeg)
{
	TestDerivative(CrossLegBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeCube)
{
        TestDerivative(CubeBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeDavis)
{
        TestDerivative(DavisBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeDeb1)
{
        TestDerivative(Deb1Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeDeckkersAarts)
{
        TestDerivative(DeckkersAartsBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeDixonPrice)
{
        TestDerivative(DixonPriceBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeDolan)
{
        TestDerivative(DolanBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeDropWave)
{
        TestDerivative(DropWaveBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeEasom)
{
        TestDerivative(EasomBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeEggCrate)
{
        TestDerivative(EggCrateBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeEggHolder)
{
        TestDerivative(EggHolderBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeElAttarVidyasagarDutt)
{
        TestDerivative(ElAttarVidyasagarDuttBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeEngvall)
{
        TestDerivative(EngvallBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeExp2)
{
        TestDerivative(Exp2Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeExponential)
{
        TestDerivative(ExponentialBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeFreudensteinRoth)
{
        TestDerivative(FreudensteinRothBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeGoldsteinPrice)
{
        TestDerivative(GoldsteinPriceBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeGramacyLee2)
{
        TestDerivative(GramacyLee2Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeGramacyLee3)
{
        TestDerivative(GramacyLee3Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeGriewank)
{
        
        TestDerivative(GriewankBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeHansen)
{
        TestDerivative(HansenBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeHartman3)
{
        TestDerivative(Hartman3Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeHartman6)
{
        TestDerivative(Hartman6Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeHelicalValley)
{
        TestDerivative(HelicalValleyBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeHimmelblau)
{
        TestDerivative(HimmelblauBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeHosaki)
{
        TestDerivative(HosakiBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeJennrichSampson)
{
        TestDerivative(JennrichSampsonBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeKeane)
{
        TestDerivative(KeaneBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeLangerman5)
{
        TestDerivative(Langerman5Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeLeon)
{
        TestDerivative(LeonBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeMatyas)
{
        TestDerivative(MatyasBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeMcCormick)
{
        TestDerivative(McCormickBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeMieleCantrell)
{
        TestDerivative(MieleCantrellBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeMishra3)
{
        TestDerivative(Mishra3Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeMishra4)
{
        TestDerivative(Mishra4Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeMishra5)
{
        TestDerivative(Mishra5Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeMishra6)
{
        TestDerivative(Mishra6Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeMishra7)
{
        TestDerivative(Mishra7Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeMishra8)
{
        TestDerivative(Mishra8Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeMishra9)
{
        TestDerivative(Mishra9Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeParsopoulos)
{
        TestDerivative(ParsopoulosBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativePathological)
{       
        TestDerivative(PathologicalBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativePeriodic)
{
        TestDerivative(PeriodicBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativePinter)
{
        TestDerivative(PinterBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativePowellSingular2)
{
        TestDerivative(PowellSingular2Benchmark<double>(8));
}

TEST_F(DerivativeTest, TestDerivativePowellSum)
{
        TestDerivative(PowellSingular2Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativePrice1)
{
        TestDerivative(Price1Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativePrice2)
{
        TestDerivative(Price2Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativePrice3)
{
        TestDerivative(Price3Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativePrice4)
{
        TestDerivative(Price4Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeProblem02)
{
        TestDerivative(Problem02Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeProblem04)
{
        TestDerivative(Problem04Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeProblem05)
{
        TestDerivative(Problem05Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeProblem06)
{
        TestDerivative(Problem06Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeQing)
{
        TestDerivative(QingBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeQuadratic)
{
        TestDerivative(QuadraticBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeQuintic)
{
        TestDerivative(QuinticBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeRosenbrock)
{
        TestDerivative(RosenbrockBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeRosenbrockModified)
{
        TestDerivative(RosenbrockModifiedBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeRotatedEllipse)
{
        TestDerivative(RotatedEllipseBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeRotatedEllipse2)
{
        TestDerivative(RotatedEllipse2Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeScahffer1)
{
        TestDerivative(Scahffer1Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeScahffer3)
{
        TestDerivative(Scahffer3Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeScahffer4)
{
        TestDerivative(Scahffer4Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeScahffer2_6)
{
        TestDerivative(Scahffer2_6Benchmark<double>());
}


TEST_F(DerivativeTest, TestDerivativeSchafferF6)
{
        
        TestDerivative(SchafferF6Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeSchmidtVetters)
{
        TestDerivative(SchmidtVettersBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeSchumerSteiglitz)
{
        TestDerivative(SchumerSteiglitzBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeSchwefel)
{
        TestDerivative(SchwefelBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeSchwefel1_2)                   
{

        TestDerivative(Schwefel1_2Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeSchwefel2_20)
{
        TestDerivative(Schwefel2_20Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeSchwefel2_22)
{
        TestDerivative(Schwefel2_20Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeSchwefel2_23)
{
        TestDerivative(Schwefel2_23Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeSchwefel2_26)
{
        TestDerivative(Schwefel2_26Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeSchwefel2_36)
{
	TestDerivative(Schwefel2_36Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeSchwefel2_4)
{
        TestDerivative(Schwefel2_4Benchmark<double>(3));
}


TEST_F(DerivativeTest, TestDerivativeShekel10)
{
        TestDerivative(Shekel10Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeShekel5)
{
        TestDerivative(Shekel5Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeShekel7)
{
        TestDerivative(Shekel7Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeShubert)
{
        TestDerivative(ShubertBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeShubert2)
{
        TestDerivative(Shubert2Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeShubert3)
{
        TestDerivative(Shubert3Benchmark<double>());
}


TEST_F(DerivativeTest, TestDerivativeSolomon)
{
        TestDerivative(SolomonBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeSphere)
{
        TestDerivative(SphereBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeStrechedVSineWave)
{
        TestDerivative(StrechedVSineWaveBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeStyblinskiTang)
{
        TestDerivative(StyblinskiTangBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeSumSquares)
{
        TestDerivative(SumSquaresBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeTable1HolderTable1)
{
        TestDerivative(Table1HolderTable1Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeTable2HolderTable2)
{
        TestDerivative(Table2HolderTable2Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeTable3Carrom)
{
        TestDerivative(Table3CarromBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeTesttubeHolder)
{
        TestDerivative(TesttubeHolderBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeTrecanni)
{
        TestDerivative(TrecanniBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeTrefethen)
{
        TestDerivative(TrefethenBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeTrid10)
{
	TestDerivative(Trid10Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeTrid6)
{
	TestDerivative(Trid6Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeTrigonometric1)
{
        TestDerivative(Trigonometric1Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeTrigonometric2)
{
        TestDerivative(Trigonometric2Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeTripod)
{
        TestDerivative(TripodBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeUrsem1)
{
        TestDerivative(Ursem1Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeUrsem3)
{
        TestDerivative(Ursem3Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeUrsem4)
{
        TestDerivative(Ursem4Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeUrsemWaves)
{
        TestDerivative(UrsemWavesBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeVenterSobiezcczanskiSobieski)
{
        TestDerivative(VenterSobiezcczanskiSobieskiBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeWWavy)
{
        
        TestDerivative(WWavyBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeWayburnSeader1)
{
        TestDerivative(WayburnSeader1Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeWayburnSeader2)
{
        TestDerivative(WayburnSeader2Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeWayburnSeader3)
{
        TestDerivative(WayburnSeader3Benchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeWeierstrass)
{
        
        TestDerivative(WeierstrassBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeWhitley)
{
        
        TestDerivative(WhitleyBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeWolfe)
{
	TestDerivative(WolfeBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeXinSheYang2)
{
        
        TestDerivative(XinSheYang2Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeXinSheYang3)
{
        
        TestDerivative( XinSheYang3Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeXinSheYang4)
{
        
        TestDerivative(XinSheYang4Benchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeZakharov)
{
        
        TestDerivative(ZakharovBenchmark<double>(3));
}

TEST_F(DerivativeTest, TestDerivativeZettl)
{
        TestDerivative(ZettlBenchmark<double>());
}

TEST_F(DerivativeTest, TestDerivativeZirilli)
{
        TestDerivative(ZirilliBenchmark<double>());
}



#endif






