//
// Created by lamure on 08/06/18.
//

#ifndef VALENTIN_EMD_TEST_HPP
#define VALENTIN_EMD_TEST_HPP

#endif //VALENTIN_EMD_TEST_HPP

#include <iostream>
#include <vector>

#include "armadillo"



using namespace std;

using namespace arma;

typedef unsigned int node_id_type;

enum ProblemType {
    INFEASIBLE,
    OPTIMAL,
    UNBOUNDED,
    MAX_ITER_REACHED
};

int EMD(int n1,int n2, mat max);

#endif