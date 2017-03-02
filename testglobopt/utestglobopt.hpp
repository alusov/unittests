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
#include "expression/testfuncs.hpp"
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

using namespace snowgoose::expression;
  


class GlobOptTest : public ::testing::Test {
 protected:

  GlobOptTest() : dfr(FUNCDESCR)
  {
  }
  
  virtual void SetUp() 
  {
		
  }
  
  virtual void TearDown() 
  {
      
  }

  void testglobopt(const std::string& key,const Expr<double>& expr, const Expr<Interval<double>> & exprInterval, int customDim = 0, double epsilon = EPSILON) 
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
    NUC::RecordSupplier<double> rs(std::numeric_limits<double>::max());
    COMPI::Functor<double>* pf = mpp->mObjectives.at(0);
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
        const int maxCnt = MAX_COUNT;
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
            std::cout << "x[" << i << "]=" << *(record+i) << ' ';
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

/*this test is for full interval library version only*/
/*
TEST_F(GlobOptTest, TestGlobOptAckley4)
{
	auto expr = Ackley4<double>();
    auto exprInterval = Ackley4<Interval<double>>();
	testglobopt(K.Ackley4, expr, exprInterval);
}*/

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


#endif /* UTESTGLOBOPT_HPP */

