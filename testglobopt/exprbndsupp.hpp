/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   exprbndsupp.hpp
 * Author: alusov
 *
 * Created on February 27, 2017, 10:27 AM
 */

#ifndef EXPRBNDSUPP_HPP
#define EXPRBNDSUPP_HPP

#include <math.h>
#include <box/boxutils.hpp>
#include "expression/algorithm.hpp"
#include <cutfact/lbcutfact/boundsupp.hpp>
#include "expression/expr.hpp"


using namespace snowgoose::interval;
using namespace snowgoose::expression;

/**
 * Interval analysis lower bound supplier for any expression 
 */
class ExprBoundSupplier : public NUC::BoundSupplier <double>{
public:

    /**
     * Constructor
     * @param n problem dimension
     */
    ExprBoundSupplier(int n, Expr<Interval<double>> expr) : mExpr(expr), mN(n) {

    }

    /**
     * Retrieve 
     * @param box to find the bound
     * @return bound
     */
    double getBound(const snowgoose::Box<double>& box) {
        std::vector<Interval<double>> intervals;
        for(int i=0; i < mN; i++)
        {
            Interval<double> interval(box.mA[i], box.mB[i]);
            intervals.push_back(interval);
        }
        return mExpr.calc(intervals, InterEvalAlg<double>()).lb();
    }

private:
    Expr<Interval<double>> mExpr;
    int mN;
};


#endif /* EXPRBNDSUPP_HPP */

