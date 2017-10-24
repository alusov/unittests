#ifndef UTESTFUNCS_HPP
#define UTESTFUNCS_HPP

#include <limits.h>
#include "gtest/gtest.h"
#include "testfuncs/testfuncs.hpp"
#include "expression/expr.hpp"
#include "expression/algder.hpp"
#include "derivatives/valder.hpp"
#include "descfunc/descfunc.hpp"
#include "descfunc/keys.hpp"
#include "pointgen/randpointgen.hpp"
#include "box/box.hpp"


#define EPSILON 0.1
#define PERCENT 1
#define DELTA 0.0001

const char* JSONPATH;

using namespace snowgoose;
using namespace snowgoose::expression;
using namespace snowgoose::derivative;
using namespace OPTITEST;


class DerivativeTest : public ::testing::Test {
 protected:
	DerivativeTest() : dfr(JSONPATH)
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
        
	void TestDerivative(const std::string& key, const Expr<double>& exprFunc, const Expr<ValDer<double>>& exprDer, int customDim=1)
	{
        std::vector<double> point = getRandomPoint(key, customDim);
		auto der = exprDer.calc(ValDerAlg<double>(point));
        double func_val = exprFunc.calc(FuncAlg<double>(point));
        double func_val_by_der = der.value();
        double epsilon = EPSILON;
        ASSERT_NEAR(func_val, func_val_by_der, epsilon);
        
        auto grad = der.grad();
        	
        for(int i=0; i < point.size(); ++i)
        {
            auto new_point = std::vector<double>(point);
            new_point[i]+=DELTA;
            double new_func_val = exprFunc.calc(FuncAlg<double>(new_point));;
            double partial_derivative = (new_func_val - func_val)/DELTA;    
            std::cout << "func_val=" << func_val << " new_func_val=" << new_func_val << '\n';
            double expected = 100;
            double derivative_relative_value = (partial_derivative/grad[i])*100;
            double persent = PERCENT;
            ASSERT_NEAR(expected, derivative_relative_value, persent);
        }
	}
	DescFuncReader dfr;
};



TEST_F(DerivativeTest, TestDerivativeAckley1)
{
	int N = 3;
	TestDerivative(K.Ackley1, Ackley1<double>(N), Ackley1<ValDer<double>>(N), N);
}



TEST_F(DerivativeTest, TestDerivativeAckley2)
{
	int N = 4;
	TestDerivative(K.Ackley2, Ackley2<double>(N), Ackley2<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeAckley3)
{
	TestDerivative(K.Ackley3, Ackley3<double>(), Ackley3<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeAckley4)
{
	TestDerivative(K.Ackley4, Ackley4<double>(), Ackley4<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeAdjiman)
{
	TestDerivative(K.Adjiman, Adjiman<double>(), Adjiman<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeAlpine1)
{
	int N = 3;
	TestDerivative(K.Alpine1, Alpine1<double>(N), Alpine1<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeAlpine2)
{
	TestDerivative(K.Alpine2, Alpine2<double>(), Alpine2<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBrad)
{
	TestDerivative(K.Brad, Brad<double>(), Brad<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBartelsConn)
{
	TestDerivative(K.BartelsConn, BartelsConn<double>(), BartelsConn<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBeale)
{
	TestDerivative(K.Beale, Beale<double>(), Beale<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBiggsExpr2)
{
	TestDerivative(K.BiggsEXP2, BiggsExpr2<double>(), BiggsExpr2<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBiggsExpr3)
{
	TestDerivative(K.BiggsEXP3, BiggsExpr3<double>(), BiggsExpr3<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBiggsExpr4)
{
	TestDerivative(K.BiggsEXP4, BiggsExpr4<double>(), BiggsExpr4<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBiggsExpr5)
{
	TestDerivative(K.BiggsEXP5, BiggsExpr5<double>(), BiggsExpr5<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBiggsExpr6)
{
	TestDerivative(K.BiggsEXP6, BiggsExpr6<double>(), BiggsExpr6<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBird)
{
	TestDerivative(K.Bird, Bird<double>(), Bird<ValDer<double>>());
}


TEST_F(DerivativeTest, TestDerivativeBohachevsky1)
{
        TestDerivative(K.Bohachevsky1, Bohachevsky1<double>(), Bohachevsky1<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBohachevsky2)
{
        TestDerivative(K.Bohachevsky2, Bohachevsky2<double>(), Bohachevsky2<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBohachevsky3)
{
        TestDerivative(K.Bohachevsky3, Bohachevsky3<double>(), Bohachevsky3<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBooth)
{
        TestDerivative(K.Booth, Booth<double>(), Booth<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBoxBettsQuadraticSum)
{
        TestDerivative(K.BoxBettsQuadraticSum, BoxBettsQuadraticSum<double>(), BoxBettsQuadraticSum<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBraninRCOS)
{
        TestDerivative(K.BraninRCOS, BraninRCOS<double>(), BraninRCOS<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBraninRCOS2)
{
        TestDerivative(K.BraninRCOS2, BraninRCOS2<double>(), BraninRCOS2<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBrent)
{
        TestDerivative(K.Brent, Brent<double>(), Brent<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBukin2)
{
        TestDerivative(K.Bukin2, Bukin2<double>(), Bukin2<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBukin4)
{
        TestDerivative(K.Bukin4, Bukin4<double>(), Bukin4<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeBukin6)
{
        TestDerivative(K.Bukin6, Bukin6<double>(), Bukin6<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeCamelSixHump)
{
        TestDerivative(K.CamelSixHump, CamelSixHump<double>(), CamelSixHump<ValDer<double>>());
}


TEST_F(DerivativeTest, TestDerivativeCamelThreeHump)
{
        TestDerivative(K.CamelThreeHump, CamelThreeHump<double>(), CamelThreeHump<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeChichinadze)
{
        TestDerivative(K.Chichinadze, Chichinadze<double>(), Chichinadze<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeChungReynolds)
{
        int N = 3;
        TestDerivative(K.ChungReynolds, ChungReynolds<double>(N), ChungReynolds<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeColville)
{
        TestDerivative(K.Colville, Colville<double>(), Colville<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeComplex)
{
        TestDerivative(K.Complex, Complex<double>(), Complex<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeCosineMixture)
{
	TestDerivative(K.CosineMixture, CosineMixture<double>(), CosineMixture<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeCrossInTray)
{
	TestDerivative(K.CrossInTray, CrossInTray<double>(), CrossInTray<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeCrossLeg)
{
	TestDerivative(K.CrossLeg, CrossLeg<double>(), CrossLeg<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeCube)
{
        TestDerivative(K.Cube, Cube<double>(), Cube<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeDavis)
{
        TestDerivative(K.Davis, Davis<double>(), Davis<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeDeb1)
{
        int N = 3;
        TestDerivative(K.Deb1, Deb1<double>(N), Deb1<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeDeckkersAarts)
{
        TestDerivative(K.DeckkersAarts, DeckkersAarts<double>(), DeckkersAarts<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeDixonPrice)
{
        TestDerivative(K.DixonPrice, DixonPrice<double>(), DixonPrice<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeDolan)
{
        TestDerivative(K.Dolan, Dolan<double>(), Dolan<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeDropWave)
{
        TestDerivative(K.DropWave, DropWave<double>(), DropWave<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeEasom)
{
        TestDerivative(K.Easom, Easom<double>(), Easom<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeEggCrate)
{
        TestDerivative(K.EggCrate, EggCrate<double>(), EggCrate<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeEggHolder)
{
        TestDerivative(K.EggHolder, EggHolder<double>(), EggHolder<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeElAttarVidyasagarDutt)
{
        TestDerivative(K.ElAttarVidyasagarDutt, ElAttarVidyasagarDutt<double>(), ElAttarVidyasagarDutt<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeEngvall)
{
        TestDerivative(K.Engvall, Engvall<double>(), Engvall<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeExp2)
{
        TestDerivative(K.Exp2, Exp2<double>(), Exp2<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeExponential)
{
        int N = 3;
        TestDerivative(K.Exponential, Exponential<double>(N), Exponential<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeFreudensteinRoth)
{
        TestDerivative(K.FreudensteinRoth, FreudensteinRoth<double>(), FreudensteinRoth<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeGoldsteinPrice)
{
        TestDerivative(K.GoldsteinPrice, GoldsteinPrice<double>(), GoldsteinPrice<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeGramacyLee2)
{
        TestDerivative(K.GramacyLee2, GramacyLee2<double>(), GramacyLee2<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeGramacyLee3)
{
        TestDerivative(K.GramacyLee3, GramacyLee3<double>(), GramacyLee3<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeGriewank)
{
        int N = 3;
        TestDerivative(K.Griewank, Griewank<double>(N), Griewank<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeHansen)
{
        TestDerivative(K.Hansen, Hansen<double>(), Hansen<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeHartman3)
{
        TestDerivative(K.Hartman3, Hartman3<double>(), Hartman3<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeHartman6)
{
        TestDerivative(K.Hartman6, Hartman6<double>(), Hartman6<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeHelicalValley)
{
        TestDerivative(K.HelicalValley, HelicalValley<double>(), HelicalValley<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeHimmelblau)
{
        TestDerivative(K.Himmelblau, Himmelblau<double>(), Himmelblau<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeHosaki)
{
        TestDerivative(K.Hosaki, Hosaki<double>(), Hosaki<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeJennrichSampson)
{
        TestDerivative(K.JennrichSampson, JennrichSampson<double>(), JennrichSampson<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeKeane)
{
        TestDerivative(K.Keane, Keane<double>(), Keane<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeLangerman5)
{
        TestDerivative(K.Langerman5, Langerman5<double>(), Langerman5<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeLeon)
{
        TestDerivative(K.Leon, Leon<double>(), Leon<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeMatyas)
{
        TestDerivative(K.Matyas, Matyas<double>(), Matyas<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeMcCormick)
{
        TestDerivative(K.McCormick, McCormick<double>(), McCormick<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeMieleCantrell)
{
        TestDerivative(K.MieleCantrell, MieleCantrell<double>(), MieleCantrell<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeMishra3)
{
        TestDerivative(K.Mishra3, Mishra3<double>(), Mishra3<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeMishra4)
{
        TestDerivative(K.Mishra4, Mishra4<double>(), Mishra4<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeMishra5)
{
        TestDerivative(K.Mishra5, Mishra5<double>(), Mishra5<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeMishra6)
{
        TestDerivative(K.Mishra6, Mishra6<double>(), Mishra6<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeMishra7)
{
        TestDerivative(K.Mishra7, Mishra7<double>(), Mishra7<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeMishra8)
{
        TestDerivative(K.Mishra8, Mishra8<double>(), Mishra8<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeMishra9)
{
        TestDerivative(K.Mishra9, Mishra9<double>(), Mishra9<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeParsopoulos)
{
        TestDerivative(K.Parsopoulos, Parsopoulos<double>(), Parsopoulos<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativePathological)
{
        int N = 3;
        TestDerivative(K.Pathological, Pathological<double>(N), Pathological<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativePeriodic)
{
        TestDerivative(K.Periodic, Periodic<double>(), Periodic<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativePinter)
{
        int N = 3;
        TestDerivative(K.Pinter, Pinter<double>(N), Pinter<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativePowellSingular2)
{
        int N = 8;
        TestDerivative(K.PowellSingular2, PowellSingular2<double>(N), PowellSingular2<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativePowellSum)
{
        int N = 3;
        TestDerivative(K.PowellSum, PowellSum<double>(N), PowellSum<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativePrice1)
{
        TestDerivative(K.Price1, Price1<double>(), Price1<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativePrice2)
{
        TestDerivative(K.Price2, Price2<double>(), Price2<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativePrice3)
{
        TestDerivative(K.Price3, Price3<double>(), Price3<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativePrice4)
{
        TestDerivative(K.Price4, Price4<double>(), Price4<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeProblem02)
{
        TestDerivative(K.Problem02, Problem02<double>(), Problem02<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeProblem04)
{
        TestDerivative(K.Problem04, Problem04<double>(), Problem04<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeProblem05)
{
        TestDerivative(K.Problem05, Problem05<double>(), Problem05<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeProblem06)
{
        TestDerivative(K.Problem06, Problem06<double>(), Problem06<ValDer<double>>());
}


TEST_F(DerivativeTest, TestDerivativeQing)
{
        TestDerivative(K.Qing, Qing<double>(), Qing<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeQuadratic)
{
        TestDerivative(K.Quadratic, Quadratic<double>(), Quadratic<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeQuintic)
{
        int N = 3;
        TestDerivative(K.Quintic, Quintic<double>(N), Quintic<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeRosenbrock)
{
        int N = 3;
        TestDerivative(K.Rosenbrock, Rosenbrock<double>(N), Rosenbrock<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeRosenbrockModified)
{
        TestDerivative(K.RosenbrockModified, RosenbrockModified<double>(), RosenbrockModified<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeRotatedEllipse)
{
        TestDerivative(K.RotatedEllipse, RotatedEllipse<double>(), RotatedEllipse<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeRotatedEllipse2)
{
        TestDerivative(K.RotatedEllipse2, RotatedEllipse2<double>(), RotatedEllipse2<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeScahffer1)
{
        TestDerivative(K.Scahffer1, Scahffer1<double>(), Scahffer1<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeScahffer3)
{
        TestDerivative(K.Scahffer3, Scahffer3<double>(), Scahffer3<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeScahffer4)
{
        TestDerivative(K.Scahffer4, Scahffer4<double>(), Scahffer4<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeScahffer2_6)
{
        TestDerivative(K.Scahffer2_6, Scahffer2_6<double>(), Scahffer2_6<ValDer<double>>());
}


TEST_F(DerivativeTest, TestDerivativeSchafferF6)
{
        int N = 3;
        TestDerivative(K.SchafferF6, SchafferF6<double>(N), SchafferF6<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeSchmidtVetters)
{
        TestDerivative(K.SchmidtVetters, SchmidtVetters<double>(), SchmidtVetters<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeSchumerSteiglitz)
{
        int N = 3;
        TestDerivative(K.SchumerSteiglitz, SchumerSteiglitz<double>(N), SchumerSteiglitz<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeSchwefel)
{
        int N = 3;
        TestDerivative(K.Schwefel, Schwefel<double>(N, 1.0), Schwefel<ValDer<double>>(N, 1.0), N);
}

TEST_F(DerivativeTest, TestDerivativeSchwefel1_2)
                      
{
        int N = 3;
        TestDerivative(K.Schwefel1_2, Schwefel1_2<double>(N), Schwefel1_2<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeSchwefel2_20)
{
        int N = 3;
        TestDerivative(K.Schwefel2_20, Schwefel2_20<double>(N), Schwefel2_20<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeSchwefel2_22)
{
        int N = 3;
        TestDerivative(K.Schwefel2_22, Schwefel2_22<double>(N), Schwefel2_22<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeSchwefel2_23)
{
        int N = 3;
        TestDerivative(K.Schwefel2_23, Schwefel2_23<double>(N), Schwefel2_23<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeSchwefel2_26)
{
        TestDerivative(K.Schwefel2_26, Schwefel2_26<double>(), Schwefel2_26<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeSchwefel2_36)
{
	TestDerivative(K.Schwefel2_36, Schwefel2_36<double>(), Schwefel2_36<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeSchwefel2_4)
{
        int N = 3;
        TestDerivative(K.Schwefel2_4, Schwefel2_4<double>(N), Schwefel2_4<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeShekel10)
{
        TestDerivative(K.Shekel10, Shekel10<double>(), Shekel10<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeShekel5)
{
        TestDerivative(K.Shekel5, Shekel5<double>(), Shekel5<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeShekel7)
{
        TestDerivative(K.Shekel7, Shekel7<double>(), Shekel7<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeShubert)
{
        TestDerivative(K.Shubert, Shubert<double>(), Shubert<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeShubert2)
{
        TestDerivative(K.Shubert2, Shubert2<double>(), Shubert2<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeShubert3)
{
        TestDerivative(K.Shubert3, Shubert3<double>(), Shubert3<ValDer<double>>());
}


TEST_F(DerivativeTest, TestDerivativeSolomon)
{
        TestDerivative(K.Solomon, Solomon<double>(), Solomon<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeSphere)
{
        int N = 3;
        TestDerivative(K.Sphere, Sphere<double>(N), Sphere<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeStrechedVSineWave)
{
        int N = 3;
        TestDerivative(K.StrechedVSineWave, StrechedVSineWave<double>(N), StrechedVSineWave<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeStyblinskiTang)
{
        TestDerivative(K.StyblinskiTang, StyblinskiTang<double>(), StyblinskiTang<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeSumSquares)
{
        int N = 3;
        TestDerivative(K.SumSquares, SumSquares<double>(N), SumSquares<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeTable1HolderTable1)
{
        TestDerivative(K.Table1HolderTable1, Table1HolderTable1<double>(), Table1HolderTable1<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeTable2HolderTable2)
{
        TestDerivative(K.Table2HolderTable2, Table2HolderTable2<double>(), Table2HolderTable2<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeTable3Carrom)
{
        TestDerivative(K.Table3Carrom, Table3Carrom<double>(), Table3Carrom<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeTesttubeHolder)
{
        TestDerivative(K.TesttubeHolder, TesttubeHolder<double>(), TesttubeHolder<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeTrecanni)
{
        TestDerivative(K.Trecanni, Trecanni<double>(), Trecanni<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeTrefethen)
{
        TestDerivative(K.Trefethen, Trefethen<double>(), Trefethen<ValDer<double>>());
}


TEST_F(DerivativeTest, TestDerivativeTrid10)
{
	TestDerivative(K.Trid10, Trid10<double>(), Trid10<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeTrid6)
{
	TestDerivative(K.Trid6, Trid6<double>(), Trid6<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeTrigonometric1)
{
        int N = 3;
        TestDerivative(K.Trigonometric1, Trigonometric1<double>(N), Trigonometric1<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeTrigonometric2)
{
        int N = 3;
        TestDerivative(K.Trigonometric2, Trigonometric2<double>(N), Trigonometric2<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeTripod)
{
        TestDerivative(K.Tripod, Tripod<double>(), Tripod<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeUrsem1)
{
        TestDerivative(K.Ursem1, Ursem1<double>(), Ursem1<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeUrsem3)
{
        TestDerivative(K.Ursem3, Ursem3<double>(), Ursem3<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeUrsem4)
{
        TestDerivative(K.Ursem4, Ursem4<double>(), Ursem4<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeUrsemWaves)
{
        TestDerivative(K.UrsemWaves, UrsemWaves<double>(), UrsemWaves<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeVenterSobiezcczanskiSobieski)
{
        TestDerivative(K.VenterSobiezcczanskiSobieski, VenterSobiezcczanskiSobieski<double>(), VenterSobiezcczanskiSobieski<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeWWavy)
{
        int N = 3;
        TestDerivative(K.WWavy, WWavy<double>(N), WWavy<ValDer<double>>(N), N);
}


TEST_F(DerivativeTest, TestDerivativeWayburnSeader1)
{
        TestDerivative(K.WayburnSeader1, WayburnSeader1<double>(), WayburnSeader1<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeWayburnSeader2)
{
        TestDerivative(K.WayburnSeader2, WayburnSeader2<double>(), WayburnSeader2<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeWayburnSeader3)
{
        TestDerivative(K.WayburnSeader3, WayburnSeader3<double>(), WayburnSeader3<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeWeierstrass)
{
        int N = 3;
        TestDerivative(K.Weierstrass, Weierstrass<double>(N), Weierstrass<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeWhitley)
{
        int N = 3;
        TestDerivative(K.Whitley, Whitley<double>(N), Whitley<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeWolfe)
{
	TestDerivative(K.Wolfe, Wolfe<double>(), Wolfe<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeXinSheYang2)
{
        int N = 3;
        TestDerivative(K.XinSheYang2, XinSheYang2<double>(N), XinSheYang2<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeXinSheYang3)
{
        int N = 3;
        TestDerivative(K.XinSheYang3, XinSheYang3<double>(N), XinSheYang3<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeXinSheYang4)
{
        int N = 3;
        TestDerivative(K.XinSheYang4, XinSheYang4<double>(N), XinSheYang4<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeZakharov)
{
        int N = 3;
        TestDerivative(K.Zakharov, Zakharov<double>(N), Zakharov<ValDer<double>>(N), N);
}

TEST_F(DerivativeTest, TestDerivativeZettl)
{
        TestDerivative(K.Zettl, Zettl<double>(), Zettl<ValDer<double>>());
}

TEST_F(DerivativeTest, TestDerivativeZirilli)
{
        TestDerivative(K.Zirilli, Zirilli<double>(), Zirilli<ValDer<double>>());
}



#endif






