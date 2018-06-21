//
// Created by lamure on 08/06/18.
//

/* definir histogramme 1 et 2 plus matrice de distance
 * le 1er histogramme est la liste des voisisn des landmark
 * la 2eme la liste des voisins du landmark qui ont changés
 */

#include <iostream>
#include <vector>

#include "network_simplex_simple.h"

#include "emd.h"

#include <chrono>

#include <stdio.h>


using namespace lemon;
using namespace std;

// all types should be signed
typedef int64_t arc_id_type; // {short, int, int64_t} ; Should be able to handle (n1*n2+n1+n2) with n1 and n2 the number of nodes (INT_MAX = 46340^2, I64_MAX = 3037000500^2)
typedef double supply_type; // {float, double, int, int64_t} ; Should be able to handle the sum of supplies and *should be signed* (a demand is a negative supply)
typedef double cost_type;  // {float, double, int, int64_t} ; Should be able to handle (number of arcs * maximum cost) and *should be signed*

struct TsFlow {
    int from, to;
    supply_type amount;
};

template<typename T>
cost_type sqrdistance2d(T* a, T* b) {
    return (a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]);
}

double EMD_wrap(int n1, int n2, double *X, double *Y, double *D, double *G,
             double* alpha, double* beta, int maxIter)  {
// beware M and C anre strored in row major C style!!!
    int n, m, i, cur;
    int cost =0;

    typedef FullBipartiteDigraph Digraph;
    DIGRAPH_TYPEDEFS(FullBipartiteDigraph);

    // Get the number of non zero coordinates for r and c
    n=0;
    for (int i=0; i<n1; i++) {
        double val=*(X+i);
        if (val>0) {
            n++;
        }else if(val<0){
            return INFEASIBLE;
        }
    }
    m=0;
    for (int i=0; i<n2; i++) {
        double val=*(Y+i);
        if (val>0) {
            m++;
        }else if(val<0){
            return INFEASIBLE;
        }
    }

    // Define the graph

    std::vector<int> indI(n), indJ(m);
    std::vector<double> weights1(n), weights2(m);
    Digraph di(n, m);
    NetworkSimplexSimple<Digraph,double,double, node_id_type> net(di, true, n+m, n*m, maxIter);

    // Set supply and demand, don't account for 0 values (faster)

    cur=0;
    for (int i=0; i<n1; i++) {
        double val=*(X+i);
        if (val>0) {
            weights1[ cur ] = val;
            indI[cur++]=i;
        }
    }

    // Demand is actually negative supply...

    cur=0;
    for (int i=0; i<n2; i++) {
        double val=*(Y+i);
        if (val>0) {
            weights2[ cur ] = -val;
            indJ[cur++]=i;
        }
    }


    net.supplyMap(&weights1[0], n, &weights2[0], m);

    // Set the cost of each edge
    for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
            double val=*(D+indI[i]*n2+indJ[j]);
            net.setCost(di.arcFromId(i*m+j), val);
        }
    }


    // Solve the problem with the network simplex algorithm

    int ret=net.run();

    double resultdist = net.totalCost();  // resultdist is the EMD

   /* if (ret==(int)net.OPTIMAL ) {
        cost = 0;
        Arc a; di.first(a);
        for (; a != INVALID; di.next(a)) {
            int i = di.source(a);
            int j = di.target(a);
            double flow = net.flow(a);
            cout << " i 1 : " << indI[i] << " ind J " << n2+indJ[j-n] << " D " << D << endl;
            cout << "multiplication  : " << D+indI[i]*n2+indJ[j-n] <<endl;
            cost = cost + flow * (*(D+indI[i]*n2+indJ[j-n]));
            *(G+indI[i]*n2+indJ[j-n]) = flow;

        }

    }*/


    return resultdist;
}