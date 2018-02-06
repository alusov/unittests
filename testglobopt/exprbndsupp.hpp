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
#include "testfuncs/manydim/benchmarks.hpp"


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
    ExprBoundSupplier(const PtrBench<double> benchmark) : bm(benchmark) {

    }

    /**
     * Retrieve 
     * @param box to find the bound
     * @return low bound
     */
    double getBound(const snowgoose::Box<double>& box) {
        std::vector<Interval<double>> intervals;
        for(int i=0; i < bm->getDim(); i++)
        {
            Interval<double> interval(box.mA[i], box.mB[i]);
            intervals.push_back(interval);
        }
        return bm->calcInterval(intervals).lb();
    }

private:
    PtrBench<double> bm;
};


#endif /* EXPRBNDSUPP_HPP */

