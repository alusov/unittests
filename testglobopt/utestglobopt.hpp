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
#include "testfuncs/testfuncs.hpp"
#include "testfuncs/benchmarks.hpp"
#include "expression/expr.hpp"
#include "expression/algorithm.hpp"
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
    OPTITEST::ExprProblemFactory fact(bm);
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
    	int N = 3;
    	testglobopt(std::make_shared<Ackley1Benchmark<double>>(N));
}

TEST_F(GlobOptTest, TestGlobOptAckley3)
{
	testglobopt(std::make_shared<Ackley3Benchmark<double>>());
}

TEST_F(GlobOptTest, TestGlobOptAlpine2)
{
    	testglobopt(std::make_shared<Alpine2Benchmark<double>>(), 0.01);
}

#endif /* UTESTGLOBOPT_HPP */

