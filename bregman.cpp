//
// Created by Kave Salamatian on 07/06/2018.
//
#include <string>
#include <fstream>
#include <utility>
#include <iostream>
#include <math.h>
#include<algorithm>

#include <armadillo>
#include "ot.h"

using namespace std;
using namespace arma;




double find_theta(double x, double theta_max) {

    double theta = theta_max;
    double logx = log(x);
    double X, Y;
    for (int k = 1; k <= 20; k++) {
        X = x * (1 - pow(x, -theta)) - theta * logx;
        Y = logx * (pow(x, 1 - theta) - 1);
        theta = theta - X / Y;
    }
    return theta;
}



OTResult sinkhorn_knopp(rowvec  &a, vec &b, mat &C, double epsilon, double theta_max , double tolerance, int maxIter){
    // Compute N dual-Sinkhorn divergences with an overrelaxation version of Cuturi's algortihm
    //
    //---------------------------
    //Required Inputs:
    //---------------------------
    //a is a d_1 x N matrix, where each column vector is in the probability simplex
    //b is a d2 x N matrix of N vectors in the probability simplex
    //C is a d1 x d2 cost matrix of pairwise distances between bins described in a_i and b_j
    //epsilon is the regularization parameter
    //theta_max\in [1;2) is the maximum value of the overrelaxation parameter
    //tolerance >0
    //maxIter: maximal number of iterations.
    //---------------------------
    //Output
    //---------------------------
    // D : vector of N dual-sinkhorn divergences
    // u : d1 x N matrix of left scalings
    // v : d2 x N matrix of right scalings
    // The smoothed optimal transport between (a_i,b_i) can be recovered as
    // T_i = diag(u(:,i)) * exp(-C/epsilon) * diag(v(:,i));


    if ((theta_max<1) || (theta_max>=2)) {
        cout<<"The overrelaxation parameter theta_max should be in [1;2)";
        theta_max=max(1.0,min(1.99999,theta_max));
    }
    if (a.min()<=0) {
        cout<<"The input data a should contain non negative values"<<endl;
        a.transform( [](double val) { return max(1e-15, val); } );
    }
    if (b.min()<=0) {
        cout<<"The input datb a should contain non negative values"<<endl;
        b.transform( [](double val) { return max(1e-15, val); } );
    }
    if ((sum(a)-1) >1e-8) {
        cout<<"The input data a does not contain normalized histograms. Sum is "<<sum(a)<<endl;
        rowvec aa(a.n_cols);
        aa.fill(sum(a));
        a=a/aa;
    }
    if ((sum(b)-1) >1e-15) {
        cout<<"The input data a does not contain normalized histograms"<<endl;
        rowvec bb(a.n_cols);
        bb.fill(sum(b));
        b=b/bb;
    }
    if ((a.n_cols != C.n_rows) || (b.n_rows!=C.n_cols)) {
        cout<<"The dimension of the cost matrix C is not size(a,1) x size(b,1)"<<endl;
        OTResult otresult;
        otresult.result = 0;
        return otresult;
    }
//    for(int n=0; n<b.n_rows; n++){
//        cout<<C(0,n)<<":"<<C(1,n)<<":"<<C(2,n)<<" '";
//    }

    // To prevent from numerical errors
    mat K = mat(C).transform([&epsilon](double val) { return max(1e-64,exp(-val/epsilon)); });
    // Initialization of Left scaling Factors, N column vectors.

    rowvec u;
    u.ones(a.n_cols);
    vec v;
    v.ones(b.n_rows);
    float delta=0.001;
    int compt=0;
    rowvec u_tilde, uu;
    vec v_tilde, vv;
    double criterion;
    double theta;
    for (int compt=0; compt<maxIter; compt=compt+1) {
        //SK iteration
        u_tilde = a / (K * v).t();
        //find the theta ensuring the decrease of the Lyapunov function:
        uu=u/u_tilde;
        theta=min(theta_max,max(1.0,find_theta(uu.min(),theta_max)-delta));
        //overrelaxation
        u = u.transform([&theta](double val){return pow(val,1-theta);})%u_tilde.transform([&theta](double val){return pow(val,theta);});
        //SK iteration
        v_tilde = b / (u*K).t();
        //find the theta ensuring the decrease of the Lyapunov function:
        vv=v/v_tilde;
        theta=min(theta_max,max(1.0,find_theta(vv.min(),theta_max)-delta));
        //overrelaxation
        v = v.transform([&theta](double val){return pow(val,1-theta);})%v_tilde.transform([&theta](double val){return pow(val,theta);});

        if ((compt%20==1) || (compt==maxIter)){
            criterion=sum(abs((u*K).t()%v-b));
            if ((criterion<tolerance)){
                break;
            }
        }
    }

    u = a / (K * v).t();
    //Compute the value of the Sinkhorn divergences
    float D=sum(u.t()%((C%K)*v));
    OTResult otresult;
    otresult.result =1;
    otresult.optCost=D;
    otresult.optPlan= diagmat(u)*K*diagmat(v);
    return(otresult);
}






