//
// Created by Kave Salamatian on 09/06/2018.
//

#ifndef VALENTIN_OT_H
#define VALENTIN_OT_H
#include <armadillo>

using namespace arma;

class OTResult {
    public:
        int result; // 1 for success, 0 for fail
        double optCost;
        arma::mat optPlan;
};

OTResult sinkhorn_knopp(rowvec &, vec &, mat &, double , double , double , int );

#endif //VALENTIN_OT_H
