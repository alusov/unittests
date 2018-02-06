/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   utestglobopt.hpp
 * Author: alusov
 *
 * Created on February 27, 2017, 11:00 AM
 */

#ifndef UTESTGLOBOPT_HPP
#define UTESTGLOBOPT_HPP


#include <limits>
#include "gtest/gtest.h"
#include <memory>

#include "exprbndsupp.hpp"
#include "testfuncs/manydim/testfuncs.hpp"
#include "testfuncs/manydim/benchmarks.hpp"
#include "expression/expr.hpp"
#include "expression/algorithm.hpp"
#include "oneobj/contboxconstr/benchmarkfunc.hpp"

#include <bags/bfsbag.hpp>
#include <cutfact/lbcutfact/recordsupp.hpp>
#include <cutfact/lbcutfact/lbcutfactory.hpp>
#include <applycut/basic/serialcutapp.hpp>
#include <decomp/bisectdecomp.hpp>
#include <solver/basesolver.hpp>
#include <box/boxutils.hpp>

#define EPSILON 0.001
#define MAX_COUNT 10000

using namespace snowgoose::expression;
using namespace OPTITEST;

class GlobOptTest : public ::testing::Test {
 protected:

  GlobOptTest()
  {
  }

  void testglobopt(const PtrBench<double> bm, double epsilon = EPSILON, int count = MAX_COUNT) 
  {
    int dim = bm->getDim();
    // Setup problem
    OPTITEST::BenchmarkProblemFactory fact(bm);
    COMPI::MPProblem<double> *mpp = fact.getProblem();

    //Setup bag of sub problems 
    NUC::Sub<double> sub(0, std::numeric_limits<double>::max(), *(mpp->mBox));
    NUC::BFSBag<double> bag;
    bag.putSub(sub);

    //Setup Cut Factory
    NUC::RecordSupplier<double> rs(std::numeric_limits<double>::max(), dim);
    auto pf = mpp->mObjectives.at(0);
    ExprBoundSupplier ibs(bm);
    NUC::LBCutFactory<double> cf(EPSILON, rs, ibs);

    // Setup decomposer
    NUC::BisectDecomp<double> bisdec;
    // Setup cut applicator 
    NUC::SerialCutApplicator<double> cutapp;
    // Setup solver
    NUC::BaseSolver<double> solver(bag, bisdec, cf, cutapp);

    int cnt = 0;
    auto st = [&](const NUC::BaseSolver<double>& bs){
        const int maxCnt = count;
        if (cnt++ > maxCnt) 
          return true;
        else 
          return false;
    };
    // Set stopper for solver
    solver.setStopper(st);

    double x[dim];
    //Setup sub evaluators
    auto sf = [&](NUC::Sub<double>& s) {
        snowgoose::BoxUtils::getCenter(s.mBox, x);
        double v = pf->func(x);
        rs.updateRv(v,x);
        s.mScore = ibs.getBound(s.mBox);
    };
    solver.addSubEval(sf);

    // Run solver
    solver.solve(); 
  
    double globMinY = rs.getBound(sub.mBox);		
    double expected = bm->getGlobMinY();
    
    std::cout << "Best value found : " << rs.getBound(sub.mBox) << '\n';
    std::cout << "Number of steps: " << cnt << '\n';
    double* record = rs.getPoint();
    if(record != nullptr){
        for(int i=0; i < dim; i++)
            std::cout << "x[" << i << "]=" << std::setprecision(15) << *(record+i) << ' ';
        std::cout << '\n';
    }
        
        
    ASSERT_NEAR(expected, globMinY, epsilon);
}


};


TEST_F(GlobOptTest, TestGlobOptAckley1)
{
	testglobopt(std::make_shared<Ackley1Benchmark<double>>(3));
}



TEST_F(GlobOptTest, TestGlobOptAckley2)
{
	testglobopt(std::make_shared<Ackley2Benchmark<double>>(4));
}

TEST_F(GlobOptTest, TestGlobOptAckley3)
{
	testglobopt(std::make_shared<Ackley3Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptAdjiman)
{
	testglobopt(std::make_shared<AdjimanBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptAlpine1)
{
	testglobopt(std::make_shared<Alpine1Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptAlpine2)
{
    testglobopt(std::make_shared<Alpine2Benchmark<double>>(), 0.01);
}

TEST_F(GlobOptTest, TestGlobOptBrad)
{
	testglobopt(std::make_shared<BradBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBartelsConn)
{
	testglobopt(std::make_shared<BartelsConnBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBeale)
{
	testglobopt(std::make_shared<BealeBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBiggsbenchmark2)
{
	testglobopt(std::make_shared<BiggsEXP2Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBiggsbenchmark3)
{
	testglobopt(std::make_shared<BiggsEXP3Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBiggsbenchmark4)
{
	testglobopt(std::make_shared<BiggsEXP4Benchmark<double>>(), 0.01);
}

TEST_F(GlobOptTest, TestGlobOptBiggsbenchmark5)
{
	testglobopt(std::make_shared<BiggsEXP5Benchmark<double>>(), 0.1);
}

TEST_F(GlobOptTest, TestGlobOptBiggsbenchmark6)
{
	testglobopt(std::make_shared<BiggsEXP6Benchmark<double>>(), 1);
}

TEST_F(GlobOptTest, TestGlobOptBird)
{
	testglobopt(std::make_shared<BirdBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBohachevsky1)
{
        testglobopt(std::make_shared<Bohachevsky1Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBohachevsky2)
{
        testglobopt(std::make_shared<Bohachevsky2Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBohachevsky3)
{
        testglobopt(std::make_shared<Bohachevsky3Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBooth)
{
        testglobopt(std::make_shared<BoothBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBoxBettsQuadraticSum)
{
        testglobopt(std::make_shared<BoxBettsQuadraticSumBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBraninRCOS)
{
        testglobopt(std::make_shared<BraninRCOSBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBraninRCOS2)
{
        testglobopt(std::make_shared<BraninRCOS2Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBrent)
{
        testglobopt(std::make_shared<BrentBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBrown)
{
        
        testglobopt(std::make_shared<BrownBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptBukin2)
{
        testglobopt(std::make_shared<Bukin2Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBukin4)
{
        testglobopt(std::make_shared<Bukin4Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBukin6)
{
        testglobopt(std::make_shared<Bukin6Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptCamelSixHump)
{
        testglobopt(std::make_shared<CamelSixHumpBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptCamelThreeHump)
{
        testglobopt(std::make_shared<CamelThreeHumpBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptChichinadze)
{
        testglobopt(std::make_shared<ChichinadzeBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptChungReynolds)
{
        
        testglobopt(std::make_shared<ChungReynoldsBenchmark<double>>(3), 0.3, 100000);
}

TEST_F(GlobOptTest, TestGlobOptColville)
{
        testglobopt(std::make_shared<ColvilleBenchmark<double>>(), 0.1, 100000);
}

TEST_F(GlobOptTest, TestGlobOptComplex)
{
        testglobopt(std::make_shared<ComplexBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptCosineMixture)
{
	testglobopt(std::make_shared<CosineMixtureBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptCrossInTray)
{
	testglobopt(std::make_shared<CrossInTrayBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptCrossLeg)
{
	testglobopt(std::make_shared<CrossLegBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptCube)
{
        testglobopt(std::make_shared<CubeBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptDavis)
{
        testglobopt(std::make_shared<DavisBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptDeb1)
{
        testglobopt(std::make_shared<Deb1Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptDeckkersAarts)
{
        testglobopt(std::make_shared<DeckkersAartsBenchmark<double>>(), 0.1);
}

TEST_F(GlobOptTest, TestGlobOptDixonPrice)
{
        testglobopt(std::make_shared<DixonPriceBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptDolan)
{
        testglobopt(std::make_shared<DolanBenchmark<double>>(), 0.001, 100000);
}

TEST_F(GlobOptTest, TestGlobOptDropWave)
{
        testglobopt(std::make_shared<DropWaveBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptEasom)
{
        testglobopt(std::make_shared<EasomBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptEggCrate)
{
        testglobopt(std::make_shared<EggCrateBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptEggHolder)
{
        testglobopt(std::make_shared<EggHolderBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptElAttarVidyasagarDutt)
{
        testglobopt(std::make_shared<ElAttarVidyasagarDuttBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptEngvall)
{
        testglobopt(std::make_shared<EngvallBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptExp2)
{
        testglobopt(std::make_shared<Exp2Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptExponential)
{
        testglobopt(std::make_shared<ExponentialBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptFreudensteinRoth)
{
        testglobopt(std::make_shared<FreudensteinRothBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptGoldsteinPrice)
{
        testglobopt(std::make_shared<GoldsteinPriceBenchmark<double>>(), 0.01, 100000);
}

TEST_F(GlobOptTest, TestGlobOptGramacyLee2)
{
        testglobopt(std::make_shared<GramacyLee2Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptGramacyLee3)
{
        testglobopt(std::make_shared<GramacyLee3Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptGriewank)
{
        
        testglobopt(std::make_shared<GriewankBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptHansen)
{
        testglobopt(std::make_shared<HansenBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptHartman3)
{
        testglobopt(std::make_shared<Hartman3Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptHartman6)
{
        testglobopt(std::make_shared<Hartman6Benchmark<double>>(), 0.001, 100000);
}

TEST_F(GlobOptTest, TestGlobOptHelicalValley)
{
        testglobopt(std::make_shared<HelicalValleyBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptHimmelblau)
{
        testglobopt(std::make_shared<HimmelblauBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptHosaki)
{
        testglobopt(std::make_shared<HosakiBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptJennrichSampson)
{
        testglobopt(std::make_shared<JennrichSampsonBenchmark<double>>(), 0.01);
}

TEST_F(GlobOptTest, TestGlobOptKeane)
{
        testglobopt(std::make_shared<KeaneBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptLangerman5)
{
        testglobopt(std::make_shared<Langerman5Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptLeon)
{
        testglobopt(std::make_shared<LeonBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMatyas)
{
        testglobopt(std::make_shared<MatyasBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMcCormick)
{
        testglobopt(std::make_shared<McCormickBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMieleCantrell)
{
        testglobopt(std::make_shared<MieleCantrellBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMishra3)
{
        testglobopt(std::make_shared<Mishra3Benchmark<double>>(), 001);
}

TEST_F(GlobOptTest, TestGlobOptMishra4)
{
        testglobopt(std::make_shared<Mishra4Benchmark<double>>(), 0.01);
}

TEST_F(GlobOptTest, TestGlobOptMishra5)
{
        testglobopt(std::make_shared<Mishra5Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMishra6)
{
        testglobopt(std::make_shared<Mishra6Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMishra7)
{
        testglobopt(std::make_shared<Mishra7Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMishra8)
{
        testglobopt(std::make_shared<Mishra8Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMishra9)
{
        testglobopt(std::make_shared<Mishra9Benchmark<double>>(), 0.1, 100000);
}

TEST_F(GlobOptTest, TestGlobOptParsopoulos)
{
        testglobopt(std::make_shared<ParsopoulosBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptPathological)
{       
        testglobopt(std::make_shared<PathologicalBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptPeriodic)
{
        testglobopt(std::make_shared<PeriodicBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptPinter)
{
        testglobopt(std::make_shared<PinterBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptPowellSingular2)
{
        testglobopt(std::make_shared<PowellSingular2Benchmark<double>>(8), EPSILON, 10 * MAX_COUNT);
}

TEST_F(GlobOptTest, TestGlobOptPowellSum)
{
        testglobopt(std::make_shared<PowellSingular2Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptPrice1)
{
        testglobopt(std::make_shared<Price1Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptPrice2)
{
        testglobopt(std::make_shared<Price2Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptPrice3)
{
        testglobopt(std::make_shared<Price3Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptPrice4)
{
        testglobopt(std::make_shared<Price4Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptProblem02)
{
        testglobopt(std::make_shared<Problem02Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptProblem04)
{
        testglobopt(std::make_shared<Problem04Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptProblem05)
{
        testglobopt(std::make_shared<Problem05Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptProblem06)
{
        testglobopt(std::make_shared<Problem06Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptQing)
{
        testglobopt(std::make_shared<QingBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptQuadratic)
{
        testglobopt(std::make_shared<QuadraticBenchmark<double>>(), 0.1);
}

TEST_F(GlobOptTest, TestGlobOptQuintic)
{
        testglobopt(std::make_shared<QuinticBenchmark<double>>(3), 0.001, 100000);
}

TEST_F(GlobOptTest, TestGlobOptRosenbrock)
{
        testglobopt(std::make_shared<RosenbrockBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptRosenbrockModified)
{
        testglobopt(std::make_shared<RosenbrockModifiedBenchmark<double>>(), 1, 100000);
}

TEST_F(GlobOptTest, TestGlobOptRotatedEllipse)
{
        testglobopt(std::make_shared<RotatedEllipseBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptRotatedEllipse2)
{
        testglobopt(std::make_shared<RotatedEllipse2Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptScahffer1)
{
        testglobopt(std::make_shared<Scahffer1Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptScahffer3)
{
        testglobopt(std::make_shared<Scahffer3Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptScahffer4)
{
        testglobopt(std::make_shared<Scahffer4Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptScahffer2_6)
{
        testglobopt(std::make_shared<Scahffer2_6Benchmark<double>>());
}


TEST_F(GlobOptTest, TestGlobOptSchafferF6)
{
        
        testglobopt(std::make_shared<SchafferF6Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptSchmidtVetters)
{
        testglobopt(std::make_shared<SchmidtVettersBenchmark<double>>(), 0.001, 100000);
}

TEST_F(GlobOptTest, TestGlobOptSchumerSteiglitz)
{
        testglobopt(std::make_shared<SchumerSteiglitzBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptSchwefel)
{
        testglobopt(std::make_shared<SchwefelBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptSchwefel1_2)                   
{

        testglobopt(std::make_shared<Schwefel1_2Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptSchwefel2_20)
{
        testglobopt(std::make_shared<Schwefel2_20Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptSchwefel2_22)
{
        testglobopt(std::make_shared<Schwefel2_20Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptSchwefel2_23)
{
        testglobopt(std::make_shared<Schwefel2_23Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptSchwefel2_26)
{
        testglobopt(std::make_shared<Schwefel2_26Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptSchwefel2_36)
{
	testglobopt(std::make_shared<Schwefel2_36Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptSchwefel2_4)
{
        testglobopt(std::make_shared<Schwefel2_4Benchmark<double>>(3));
}


TEST_F(GlobOptTest, TestGlobOptShekel10)
{
        testglobopt(std::make_shared<Shekel10Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptShekel5)
{
        testglobopt(std::make_shared<Shekel5Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptShekel7)
{
        testglobopt(std::make_shared<Shekel7Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptShubert)
{
        testglobopt(std::make_shared<ShubertBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptShubert2)
{
        testglobopt(std::make_shared<Shubert2Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptShubert3)
{
        testglobopt(std::make_shared<Shubert3Benchmark<double>>());
}


TEST_F(GlobOptTest, TestGlobOptSolomon)
{
        testglobopt(std::make_shared<SolomonBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptSphere)
{
        testglobopt(std::make_shared<SphereBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptStrechedVSineWave)
{
        testglobopt(std::make_shared<StrechedVSineWaveBenchmark<double>>(3), 0.01);
}

TEST_F(GlobOptTest, TestGlobOptStyblinskiTang)
{
        testglobopt(std::make_shared<StyblinskiTangBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptSumSquares)
{
        testglobopt(std::make_shared<SumSquaresBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptTable1HolderTable1)
{
        testglobopt(std::make_shared<Table1HolderTable1Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptTable2HolderTable2)
{
        testglobopt(std::make_shared<Table2HolderTable2Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptTable3Carrom)
{
        testglobopt(std::make_shared<Table3CarromBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptTesttubeHolder)
{
        testglobopt(std::make_shared<TesttubeHolderBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptTrecanni)
{
        testglobopt(std::make_shared<TrecanniBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptTrefethen)
{
        testglobopt(std::make_shared<TrefethenBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptTrid10)
{
	testglobopt(std::make_shared<Trid10Benchmark<double>>(), 1200);
}

TEST_F(GlobOptTest, TestGlobOptTrid6)
{
	testglobopt(std::make_shared<Trid6Benchmark<double>>(), 10, 400000);
}

TEST_F(GlobOptTest, TestGlobOptTrigonometric1)
{
        testglobopt(std::make_shared<Trigonometric1Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptTrigonometric2)
{
        testglobopt(std::make_shared<Trigonometric2Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptTripod)
{
        testglobopt(std::make_shared<TripodBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptUrsem1)
{
        testglobopt(std::make_shared<Ursem1Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptUrsem3)
{
        testglobopt(std::make_shared<Ursem3Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptUrsem4)
{
        testglobopt(std::make_shared<Ursem4Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptUrsemWaves)
{
        testglobopt(std::make_shared<UrsemWavesBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptVenterSobiezcczanskiSobieski)
{
        testglobopt(std::make_shared<VenterSobiezcczanskiSobieskiBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptWWavy)
{
        
        testglobopt(std::make_shared<WWavyBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptWatson)
{
        
        testglobopt(std::make_shared<WatsonBenchmark<double>>(), 7);
}

TEST_F(GlobOptTest, TestGlobOptWayburnSeader1)
{
        testglobopt(std::make_shared<WayburnSeader1Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptWayburnSeader2)
{
        testglobopt(std::make_shared<WayburnSeader2Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptWayburnSeader3)
{
        testglobopt(std::make_shared<WayburnSeader3Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptWeierstrass)
{
        
        testglobopt(std::make_shared<WeierstrassBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptWhitley)
{
        
        testglobopt(std::make_shared<WhitleyBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptWolfe)
{
	testglobopt(std::make_shared<WolfeBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptXinSheYang2)
{
        
        testglobopt(std::make_shared<XinSheYang2Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptXinSheYang3)
{
        
        testglobopt(std::make_shared< XinSheYang3Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptXinSheYang4)
{
        
        testglobopt(std::make_shared<XinSheYang4Benchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptZakharov)
{
        
        testglobopt(std::make_shared<ZakharovBenchmark<double>>(3));
}

TEST_F(GlobOptTest, TestGlobOptZettl)
{
        testglobopt(std::make_shared<ZettlBenchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptZirilli)
{
        testglobopt(std::make_shared<ZirilliBenchmark<double>>());
}



#endif /* UTESTGLOBOPT_HPP */

