#include <iostream>
#include <mat.h>
#include "armadillo"
#include "armahelp.h"
#include "globals.h"
#include "msllcpp.h"
#include "mex.h" 

using namespace arma;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
    const mxArray *theta_array;
    double *output, ll=0.0;

    if (nrhs != 2) { 
	    mexErrMsgIdAndTxt( "MATLAB:choiceprob:invalidNumInputs",
                          "Two input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "MATLAB:choiceprob:maxlhs",
                          "Too many output arguments."); 
    }   
    
    VerifyDataset(prhs[0]);

    if(mxGetM(theta_array=prhs[1])!=N+L+J){
        mexErrMsgIdAndTxt( "MATLAB:choiceprob:theta_array",
                          "Theta is wrong dimension. Should be %d", N+L+J); 
    }

    mat pjt = get_probs(mxGetPr(theta_array));

    plhs[0] = mxCreateNumericMatrix(J+1,N, mxDOUBLE_CLASS, mxREAL); 
    memcpy(mxGetPr(plhs[0]), pjt.memptr(),N*(J+1)*sizeof(double));
    
    return;
} 

