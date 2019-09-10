#include <iostream>
#include "armadillo"
#include "armahelp.h"

mat pos_log(mat x){
    uvec y=find(x);
    x.elem(y) = log(x.elem(y));
    return x;
}

mat element_divide(const mat A, const mat B){
    mat C;
    if(A.n_rows != B.n_rows || A.n_cols != B.n_cols){
        cout << "Matrices must have same number of rows/columns" << endl;
        exit(1);
    }
    C.copy_size(A);
# pragma omp parallel
{
# pragma omp parallel for
    for(int i=0;i<A.n_elem; i++){
        C[i] =B[i] > 0 ? A[i] / B[i] : 0.0;
    }

}
    return C;
}

mat sparse_outer(const vec x, const uvec nonzeros){
    vec a= x.elem(nonzeros);
    mat B = a * trans(a);
}

uvec counting_vector(const int a, const int offset){
	uvec y = zeros<uvec>(a,1);
	for(int i=0; i < a; i++){
		y[i] = i + offset;
	}
        return y;
}


vec accu_vec(const uvec index, const vec x){
	int ylen = max(index)+1;
	vec y = zeros<vec>(ylen,1);

	for(int i=0; i<index.n_elem ; i++){
		y(index(i)) += as_scalar(x(i));	
	}
	return y;
}
