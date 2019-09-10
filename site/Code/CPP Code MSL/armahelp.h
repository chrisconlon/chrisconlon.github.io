#ifndef ARMAHELP_H
#define ARMAHELP_H

#include <iostream>
#include "armadillo"

using namespace arma;

mat pos_log(mat x);
mat element_divide(const mat A, const mat B);
mat sparse_outer(const vec x, const uvec nonzeros);
uvec counting_vector(const int a, const int offset);
vec accu_vec(const uvec	index, const vec x);

#endif
