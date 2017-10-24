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

#include "exprbndsupp.hpp"
#include "testfuncs/testfuncs.hpp"
#include "expression/expr.hpp"
#include "expression/algorithm.hpp"
#include "descfunc/descfunc.hpp"
#include "descfunc/keys.hpp"
#include "oneobj/contboxconstr/exprfunc.hpp"

#include <bags/bfsbag.hpp>
#include <cutfact/lbcutfact/recordsupp.hpp>
#include <cutfact/lbcutfact/lbcutfactory.hpp>
#include <applycut/basic/serialcutapp.hpp>
#include <decomp/bisectdecomp.hpp>
#include <solver/basesolver.hpp>
#include <box/boxutils.hpp>

#define EPSILON 0.001
#define MAX_COUNT 10000

const char* JSONPATH;

using namespace snowgoose::expression;
using namespace OPTITEST;
  


class GlobOptTest : public ::testing::Test {
 protected:

  GlobOptTest() : dfr(JSONPATH)
  {
  }
  
  virtual void SetUp() 
  {
		
  }
  
  virtual void TearDown() 
  {
      
  }

  void testglobopt(const std::string& key,const Expr<double>& expr, const Expr<Interval<double>> & exprInterval, int customDim = 0, double epsilon = EPSILON, int count = MAX_COUNT) 
  {
    auto desc = dfr.getdesr(key);
    int dim = desc.anyDim ? customDim : desc.dim;
    Bounds bounds = desc.anyDim ? Bounds(dim, desc.bounds[0]) : desc.bounds;

    // Setup problem
    OPTITEST::ExprProblemFactory fact(expr, bounds);
    COMPI::MPProblem<double> *mpp = fact.getProblem();

    //Setup bag of sub problems
    
    NUC::Sub<double> sub(0, std::numeric_limits<double>::max(), *(mpp->mBox));
    NUC::BFSBag<double> bag;
    bag.putSub(sub);

    //Setup Cut Factory
    NUC::RecordSupplier<double> rs(std::numeric_limits<double>::max(), dim);
    auto pf = mpp->mObjectives.at(0);
    ExprBoundSupplier ibs(dim, exprInterval);
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
    double expected = desc.globMinY;
    
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

DescFuncReader dfr;
};



TEST_F(GlobOptTest, TestGlobOptAckley1)
{
    int N = 3;
	auto expr = Ackley1<double>(N);
    auto exprInterval = Ackley1<Interval<double>>(N);
    testglobopt(K.Ackley1, expr, exprInterval, N);
}

TEST_F(GlobOptTest, TestGlobOptAckley2)
{
	int N = 4;
	auto expr = Ackley2<double>(N);
    auto exprInterval = Ackley2<Interval<double>>(N);
	testglobopt(K.Ackley2, expr, exprInterval, N);
}

TEST_F(GlobOptTest, TestGlobOptAckley3)
{
	auto expr = Ackley3<double>();
    auto exprInterval = Ackley3<Interval<double>>();
	testglobopt(K.Ackley3, expr, exprInterval);
}


TEST_F(GlobOptTest, TestGlobOptAdjiman)
{
    auto expr = Adjiman<double>();
    auto exprInterval = Adjiman<Interval<double>>();
	testglobopt(K.Adjiman, expr, exprInterval);
}

TEST_F(GlobOptTest, TestGlobOptAlpine1)
{  
    int N = 3;
    auto expr = Alpine1<double>(N);
    auto exprInterval = Alpine1<Interval<double>>(N);
    testglobopt(K.Alpine1, expr, exprInterval, N);
}

TEST_F(GlobOptTest, TestGlobOptAlpine2)
{
    auto expr = Alpine2<double>();
    auto exprInterval = Alpine2<Interval<double>>();
    testglobopt(K.Alpine2, expr, exprInterval, 0, 0.01);
}

TEST_F(GlobOptTest, TestGlobOptBrad)
{
    auto expr = Brad<double>();
    auto exprInterval = Brad<Interval<double>>();
    testglobopt(K.Brad, expr, exprInterval);
}

TEST_F(GlobOptTest, TestGlobOptBartelsConn)
{
    auto expr = BartelsConn<double>();
    auto exprInterval = BartelsConn<Interval<double>>();
    testglobopt(K.BartelsConn, expr, exprInterval);
}

TEST_F(GlobOptTest, TestGlobOptBeale)
{
    auto expr = Beale<double>();
    auto exprInterval = Beale<Interval<double>>();
    testglobopt(K.Beale, expr, exprInterval);
}

TEST_F(GlobOptTest, TestGlobOptBiggsExpr2)
{
    auto expr = BiggsExpr2<double>();
    auto exprInterval = BiggsExpr2<Interval<double>>();
    testglobopt(K.BiggsEXP2, expr, exprInterval);
}

TEST_F(GlobOptTest, TestGlobOptBiggsExpr3)
{
    auto expr = BiggsExpr3<double>();
    auto exprInterval = BiggsExpr3<Interval<double>>();
    testglobopt(K.BiggsEXP3, expr, exprInterval);
}

TEST_F(GlobOptTest, TestGlobOptBiggsExpr4)
{
 
    auto expr = BiggsExpr4<double>();
    auto exprInterval = BiggsExpr4<Interval<double>>();
    testglobopt(K.BiggsEXP4, expr, exprInterval, 0, 0.01);
}

TEST_F(GlobOptTest, TestGlobOptBiggsExpr5)
{
    auto expr = BiggsExpr5<double>();
    auto exprInterval = BiggsExpr5<Interval<double>>();
    testglobopt(K.BiggsEXP5, expr, exprInterval, 0, 0.1);
}

TEST_F(GlobOptTest, TestGlobOptBiggsExpr6)
{    
    auto expr = BiggsExpr6<double>();
    auto exprInterval = BiggsExpr6<Interval<double>>();
    testglobopt(K.BiggsEXP6, expr, exprInterval, 0, 1);
}

TEST_F(GlobOptTest, TestGlobOptBird)
{
    auto expr = Bird<double>();
    auto exprInterval = Bird<Interval<double>>();
    testglobopt(K.Bird, expr, exprInterval);
}

//****************************************************

TEST_F(GlobOptTest, TestGlobOptBohachevsky1)
{
        testglobopt(K.Bohachevsky1, Bohachevsky1<double>(), Bohachevsky1<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBohachevsky2)
{
        testglobopt(K.Bohachevsky2, Bohachevsky2<double>(), Bohachevsky2<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBohachevsky3)
{
        testglobopt(K.Bohachevsky3, Bohachevsky3<double>(), Bohachevsky3<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBooth)
{
        testglobopt(K.Booth, Booth<double>(), Booth<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBoxBettsQuadraticSum)
{
        testglobopt(K.BoxBettsQuadraticSum, BoxBettsQuadraticSum<double>(), BoxBettsQuadraticSum<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBraninRCOS)
{
        testglobopt(K.BraninRCOS, BraninRCOS<double>(), BraninRCOS<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBraninRCOS2)
{
        testglobopt(K.BraninRCOS2, BraninRCOS2<double>(), BraninRCOS2<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBrent)
{
        testglobopt(K.Brent, Brent<double>(), Brent<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBrown)
{
        int N = 3;
        testglobopt(K.Brown, Brown<double>(N), Brown<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptBukin2)
{
        testglobopt(K.Bukin2, Bukin2<double>(), Bukin2<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBukin4)
{
        testglobopt(K.Bukin4, Bukin4<double>(), Bukin4<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptBukin6)
{
        testglobopt(K.Bukin6, Bukin6<double>(), Bukin6<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptCamelSixHump)
{
        testglobopt(K.CamelSixHump, CamelSixHump<double>(), CamelSixHump<Interval<double>>());
}


TEST_F(GlobOptTest, TestGlobOptCamelThreeHump)
{
        testglobopt(K.CamelThreeHump, CamelThreeHump<double>(), CamelThreeHump<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptChichinadze)
{
        testglobopt(K.Chichinadze, Chichinadze<double>(), Chichinadze<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptChungReynolds)
{
        int N = 3;
        testglobopt(K.ChungReynolds, ChungReynolds<double>(N), ChungReynolds<Interval<double>>(N), N, 0.3, 100000);
}

TEST_F(GlobOptTest, TestGlobOptColville)
{
        testglobopt(K.Colville, Colville<double>(), Colville<Interval<double>>(), 0, 0.1, 100000);
}

TEST_F(GlobOptTest, TestGlobOptComplex)
{
        testglobopt(K.Complex, Complex<double>(), Complex<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptCosineMixture)
{
	testglobopt(K.CosineMixture, CosineMixture<double>(), CosineMixture<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptCrossInTray)
{
	testglobopt(K.CrossInTray, CrossInTray<double>(), CrossInTray<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptCrossLeg)
{
	testglobopt(K.CrossLeg, CrossLeg<double>(), CrossLeg<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptCube)
{
        testglobopt(K.Cube, Cube<double>(), Cube<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptDavis)
{
        testglobopt(K.Davis, Davis<double>(), Davis<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptDeb1)
{
        int N = 3;
        testglobopt(K.Deb1, Deb1<double>(N), Deb1<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptDeckkersAarts)
{
        testglobopt(K.DeckkersAarts, DeckkersAarts<double>(), DeckkersAarts<Interval<double>>(), 0, 0.1);
}

TEST_F(GlobOptTest, TestGlobOptDixonPrice)
{
        testglobopt(K.DixonPrice, DixonPrice<double>(), DixonPrice<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptDolan)
{
        testglobopt(K.Dolan, Dolan<double>(), Dolan<Interval<double>>(), 0, 0.001, 100000);
}

TEST_F(GlobOptTest, TestGlobOptDropWave)
{
        testglobopt(K.DropWave, DropWave<double>(), DropWave<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptEasom)
{
        testglobopt(K.Easom, Easom<double>(), Easom<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptEggCrate)
{
        testglobopt(K.EggCrate, EggCrate<double>(), EggCrate<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptEggHolder)
{
        testglobopt(K.EggHolder, EggHolder<double>(), EggHolder<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptElAttarVidyasagarDutt)
{
        testglobopt(K.ElAttarVidyasagarDutt, ElAttarVidyasagarDutt<double>(), ElAttarVidyasagarDutt<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptEngvall)
{
        testglobopt(K.Engvall, Engvall<double>(), Engvall<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptExp2)
{
        testglobopt(K.Exp2, Exp2<double>(), Exp2<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptExponential)
{
        int N = 3;
        testglobopt(K.Exponential, Exponential<double>(N), Exponential<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptFreudensteinRoth)
{
        testglobopt(K.FreudensteinRoth, FreudensteinRoth<double>(), FreudensteinRoth<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptGoldsteinPrice)
{
        testglobopt(K.GoldsteinPrice, GoldsteinPrice<double>(), GoldsteinPrice<Interval<double>>(), 0, 0.01, 100000);
}

TEST_F(GlobOptTest, TestGlobOptGramacyLee2)
{
        testglobopt(K.GramacyLee2, GramacyLee2<double>(), GramacyLee2<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptGramacyLee3)
{
        testglobopt(K.GramacyLee3, GramacyLee3<double>(), GramacyLee3<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptGriewank)
{
        int N = 3;
        testglobopt(K.Griewank, Griewank<double>(N), Griewank<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptHansen)
{
        testglobopt(K.Hansen, Hansen<double>(), Hansen<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptHartman3)
{
        testglobopt(K.Hartman3, Hartman3<double>(), Hartman3<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptHartman6)
{
        testglobopt(K.Hartman6, Hartman6<double>(), Hartman6<Interval<double>>(), 0, 0.001, 100000);
}

TEST_F(GlobOptTest, TestGlobOptHelicalValley)
{
        testglobopt(K.HelicalValley, HelicalValley<double>(), HelicalValley<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptHimmelblau)
{
        testglobopt(K.Himmelblau, Himmelblau<double>(), Himmelblau<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptHosaki)
{
        testglobopt(K.Hosaki, Hosaki<double>(), Hosaki<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptJennrichSampson)
{
        testglobopt(K.JennrichSampson, JennrichSampson<double>(), JennrichSampson<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptKeane)
{
        testglobopt(K.Keane, Keane<double>(), Keane<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptLangerman5)
{
        testglobopt(K.Langerman5, Langerman5<double>(), Langerman5<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptLeon)
{
        testglobopt(K.Leon, Leon<double>(), Leon<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMatyas)
{
        testglobopt(K.Matyas, Matyas<double>(), Matyas<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMcCormick)
{
        testglobopt(K.McCormick, McCormick<double>(), McCormick<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMieleCantrell)
{
        testglobopt(K.MieleCantrell, MieleCantrell<double>(), MieleCantrell<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMishra3)
{
        testglobopt(K.Mishra3, Mishra3<double>(), Mishra3<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMishra4)
{
        testglobopt(K.Mishra4, Mishra4<double>(), Mishra4<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMishra5)
{
        testglobopt(K.Mishra5, Mishra5<double>(), Mishra5<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMishra6)
{
        testglobopt(K.Mishra6, Mishra6<double>(), Mishra6<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMishra7)
{
        testglobopt(K.Mishra7, Mishra7<double>(), Mishra7<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMishra8)
{
        testglobopt(K.Mishra8, Mishra8<double>(), Mishra8<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptMishra9)
{
        testglobopt(K.Mishra9, Mishra9<double>(), Mishra9<Interval<double>>(),0,0.1, 100000);
}

TEST_F(GlobOptTest, TestGlobOptParsopoulos)
{
        testglobopt(K.Parsopoulos, Parsopoulos<double>(), Parsopoulos<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptPathological)
{
        int N = 3;
        testglobopt(K.Pathological, Pathological<double>(N), Pathological<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptPeriodic)
{
        testglobopt(K.Periodic, Periodic<double>(), Periodic<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptPinter)
{
        int N = 3;
        testglobopt(K.Pinter, Pinter<double>(N), Pinter<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptPowellSingular2)
{
        int N = 8;
        testglobopt(K.PowellSingular2, PowellSingular2<double>(N), PowellSingular2<Interval<double>>(N), N, EPSILON, 10 * MAX_COUNT);
}

TEST_F(GlobOptTest, TestGlobOptPowellSum)
{
        int N = 3;
        testglobopt(K.PowellSum, PowellSum<double>(N), PowellSum<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptPrice1)
{
        testglobopt(K.Price1, Price1<double>(), Price1<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptPrice2)
{
        testglobopt(K.Price2, Price2<double>(), Price2<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptPrice3)
{
        testglobopt(K.Price3, Price3<double>(), Price3<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptPrice4)
{
        testglobopt(K.Price4, Price4<double>(), Price4<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptProblem02)
{
        testglobopt(K.Problem02, Problem02<double>(), Problem02<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptProblem04)
{
        testglobopt(K.Problem04, Problem04<double>(), Problem04<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptProblem05)
{
        testglobopt(K.Problem05, Problem05<double>(), Problem05<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptProblem06)
{
        testglobopt(K.Problem06, Problem06<double>(), Problem06<Interval<double>>());
}


TEST_F(GlobOptTest, TestGlobOptQing)
{
        testglobopt(K.Qing, Qing<double>(), Qing<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptQuadratic)
{
        testglobopt(K.Quadratic, Quadratic<double>(), Quadratic<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptQuintic)
{
        int N = 3;
        testglobopt(K.Quintic, Quintic<double>(N), Quintic<Interval<double>>(N), N, 0.001, 100000);
}

TEST_F(GlobOptTest, TestGlobOptRosenbrock)
{
        int N = 3;
        testglobopt(K.Rosenbrock, Rosenbrock<double>(N), Rosenbrock<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptRosenbrockModified)
{
        testglobopt(K.RosenbrockModified, RosenbrockModified<double>(), RosenbrockModified<Interval<double>>(), 0, 0.5, 10000);
}

TEST_F(GlobOptTest, TestGlobOptRotatedEllipse)
{
        testglobopt(K.RotatedEllipse, RotatedEllipse<double>(), RotatedEllipse<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptRotatedEllipse2)
{
        testglobopt(K.RotatedEllipse2, RotatedEllipse2<double>(), RotatedEllipse2<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptScahffer1)
{
        testglobopt(K.Scahffer1, Scahffer1<double>(), Scahffer1<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptScahffer3)
{
        testglobopt(K.Scahffer3, Scahffer3<double>(), Scahffer3<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptScahffer4)
{
        testglobopt(K.Scahffer4, Scahffer4<double>(), Scahffer4<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptScahffer2_6)
{
        testglobopt(K.Scahffer2_6, Scahffer2_6<double>(), Scahffer2_6<Interval<double>>());
}


TEST_F(GlobOptTest, TestGlobOptSchafferF6)
{
        int N = 3;
        testglobopt(K.SchafferF6, SchafferF6<double>(N), SchafferF6<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptSchmidtVetters)
{
        testglobopt(K.SchmidtVetters, SchmidtVetters<double>(), SchmidtVetters<Interval<double>>(), 0, 0.001, 100000);
}

TEST_F(GlobOptTest, TestGlobOptSchumerSteiglitz)
{
        int N = 3;
        testglobopt(K.SchumerSteiglitz, SchumerSteiglitz<double>(N), SchumerSteiglitz<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptSchwefel)
{
        int N = 3;
        testglobopt(K.Schwefel, Schwefel<double>(N, 1.0), Schwefel<Interval<double>>(N, 1.0), N);
}

TEST_F(GlobOptTest, TestGlobOptSchwefel1_2)
                      
{
        int N = 3;
        testglobopt(K.Schwefel1_2, Schwefel1_2<double>(N), Schwefel1_2<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptSchwefel2_20)
{
        int N = 3;
        testglobopt(K.Schwefel2_20, Schwefel2_20<double>(N), Schwefel2_20<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptSchwefel2_22)
{
        int N = 3;
        testglobopt(K.Schwefel2_22, Schwefel2_22<double>(N), Schwefel2_22<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptSchwefel2_23)
{
        int N = 3;
        testglobopt(K.Schwefel2_23, Schwefel2_23<double>(N), Schwefel2_23<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptSchwefel2_26)
{
        testglobopt(K.Schwefel2_26, Schwefel2_26<double>(), Schwefel2_26<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptSchwefel2_36)
{
	testglobopt(K.Schwefel2_36, Schwefel2_36<double>(), Schwefel2_36<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptSchwefel2_4)
{
        int N = 3;
        testglobopt(K.Schwefel2_4, Schwefel2_4<double>(N), Schwefel2_4<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptShekel10)
{
        testglobopt(K.Shekel10, Shekel10<double>(), Shekel10<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptShekel5)
{
        testglobopt(K.Shekel5, Shekel5<double>(), Shekel5<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptShekel7)
{
        testglobopt(K.Shekel7, Shekel7<double>(), Shekel7<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptShubert)
{
        testglobopt(K.Shubert, Shubert<double>(), Shubert<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptShubert2)
{
        testglobopt(K.Shubert2, Shubert2<double>(), Shubert2<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptShubert3)
{
        testglobopt(K.Shubert3, Shubert3<double>(), Shubert3<Interval<double>>());
}


TEST_F(GlobOptTest, TestGlobOptSolomon)
{
        testglobopt(K.Solomon, Solomon<double>(), Solomon<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptSphere)
{
        int N = 3;
        testglobopt(K.Sphere, Sphere<double>(N), Sphere<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptStrechedVSineWave)
{
        int N = 3;
        testglobopt(K.StrechedVSineWave, StrechedVSineWave<double>(N), StrechedVSineWave<Interval<double>>(N), N, 0.01);
}

TEST_F(GlobOptTest, TestGlobOptStyblinskiTang)
{
        testglobopt(K.StyblinskiTang, StyblinskiTang<double>(), StyblinskiTang<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptSumSquares)
{
        int N = 3;
        testglobopt(K.SumSquares, SumSquares<double>(N), SumSquares<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptTable1HolderTable1)
{
        testglobopt(K.Table1HolderTable1, Table1HolderTable1<double>(), Table1HolderTable1<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptTable2HolderTable2)
{
        testglobopt(K.Table2HolderTable2, Table2HolderTable2<double>(), Table2HolderTable2<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptTable3Carrom)
{
        testglobopt(K.Table3Carrom, Table3Carrom<double>(), Table3Carrom<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptTesttubeHolder)
{
        testglobopt(K.TesttubeHolder, TesttubeHolder<double>(), TesttubeHolder<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptTrecanni)
{
        testglobopt(K.Trecanni, Trecanni<double>(), Trecanni<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptTrefethen)
{
        testglobopt(K.Trefethen, Trefethen<double>(), Trefethen<Interval<double>>());
}


TEST_F(GlobOptTest, TestGlobOptTrid10)
{
	testglobopt(K.Trid10, Trid10<double>(), Trid10<Interval<double>>(),0, 1200);
}

TEST_F(GlobOptTest, TestGlobOptTrid6)
{
	testglobopt(K.Trid6, Trid6<double>(), Trid6<Interval<double>>(),0, 10, 400000);
}


TEST_F(GlobOptTest, TestGlobOptTrigonometric1)
{
        int N = 3;
        testglobopt(K.Trigonometric1, Trigonometric1<double>(N), Trigonometric1<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptTrigonometric2)
{
        int N = 3;
        testglobopt(K.Trigonometric2, Trigonometric2<double>(N), Trigonometric2<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptTripod)
{
        testglobopt(K.Tripod, Tripod<double>(), Tripod<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptUrsem1)
{
        testglobopt(K.Ursem1, Ursem1<double>(), Ursem1<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptUrsem3)
{
        testglobopt(K.Ursem3, Ursem3<double>(), Ursem3<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptUrsem4)
{
        testglobopt(K.Ursem4, Ursem4<double>(), Ursem4<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptUrsemWaves)
{
        testglobopt(K.UrsemWaves, UrsemWaves<double>(), UrsemWaves<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptVenterSobiezcczanskiSobieski)
{
        testglobopt(K.VenterSobiezcczanskiSobieski, VenterSobiezcczanskiSobieski<double>(), VenterSobiezcczanskiSobieski<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptWWavy)
{
        int N = 3;
        testglobopt(K.WWavy, WWavy<double>(N), WWavy<Interval<double>>(N), N);
}



TEST_F(GlobOptTest, TestGlobOptWatson)
{
        testglobopt(K.Watson, Watson<double>(), Watson<Interval<double>>(),0,7);
}


TEST_F(GlobOptTest, TestGlobOptWayburnSeader1)
{
        testglobopt(K.WayburnSeader1, WayburnSeader1<double>(), WayburnSeader1<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptWayburnSeader2)
{
        testglobopt(K.WayburnSeader2, WayburnSeader2<double>(), WayburnSeader2<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptWayburnSeader3)
{
        testglobopt(K.WayburnSeader3, WayburnSeader3<double>(), WayburnSeader3<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptWeierstrass)
{
        int N = 3;
        testglobopt(K.Weierstrass, Weierstrass<double>(N), Weierstrass<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptWhitley)
{
        int N = 3;
        testglobopt(K.Whitley, Whitley<double>(N), Whitley<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptWolfe)
{
	testglobopt(K.Wolfe, Wolfe<double>(), Wolfe<Interval<double>>());
}

TEST_F(GlobOptTest, TestGlobOptXinSheYang2)
{
        int N = 3;
        testglobopt(K.XinSheYang2, XinSheYang2<double>(N), XinSheYang2<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptXinSheYang3)
{
        int N = 3;
        testglobopt(K.XinSheYang3, XinSheYang3<double>(N), XinSheYang3<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptXinSheYang4)
{
        int N = 3;
        testglobopt(K.XinSheYang4, XinSheYang4<double>(N), XinSheYang4<Interval<double>>(N), N);
}

TEST_F(GlobOptTest, TestGlobOptZakharov)
{
        int N = 3;
        testglobopt(K.Zakharov, Zakharov<double>(N), Zakharov<Interval<double>>(N), N);
}


#endif /* UTESTGLOBOPT_HPP */

