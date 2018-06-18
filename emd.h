//
// Created by lamure on 13/06/18.
//

#ifndef VALENTIN_EMD_H
#define VALENTIN_EMD_H


#include <iostream>
#include <vector>
#include "network_simplex_simple.h"


#include <iostream>
#include <vector>

#include "network_simplex_simple.h"

using namespace lemon;

typedef unsigned int node_id_type;

enum ProblemType {
    INFEASIBLE,
    OPTIMAL,
    MAX_ITER_REACHED,
    UNBOUNDED

};


int EMD_wrap(int n1,int n2, double *X, double *Y,double *D, double *G, double* alpha, double* beta, double *cost, int maxIter);




#endif //VALENTIN_EMD_H
