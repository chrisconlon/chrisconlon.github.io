//
//  msllcpp.h
//  
//
//  Created by Christopher Conlon on 9/14/11.
//  Copyright 2011 Columbia Univeristy. All rights reserved.
//

#ifndef _msllcpp_h
#define _msllcpp_h

#include "armadillo"

using namespace arma;
using namespace std;

mat create_utils(vec theta);
mat probs_only(const mat utils, const vec theta);

double get_loglik(const double *theta);
double get_loglikg(const double *theta, double *grad);
double get_stderr(const double *theta, double *sevec);

vec parallel_likelihood(const mat utils, const vec theta, double &ll,int mode);
vec compute_se(const mat utils, const vec theta);

rowvec score_contribution(const mat pij, const vec qvec);
mat score_contribution2(const mat pij, const vec qvec, mat &A, vec &B);
rowvec gradient_contribution(const mat pij, const vec qvec);
vec compute_gradient(const mat scores);

mat get_probs(const double *theta);
double sum_FE(const double *theta);
void get_jacobian(double *jac);

#endif
