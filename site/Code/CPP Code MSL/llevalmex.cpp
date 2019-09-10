#include <iostream>
#include <mat.h>
#include "armadillo"
#include "armahelp.h"
#include "globals.h"
#include "msllcpp.h"
#include "knitro.h"
#include "omp.h"
#include "mex.h"

using namespace arma;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
    const mxArray *theta_array, *grad_array;
    double *x, *output,*grad,*std_err;
    int n;

    if (nrhs != 2) { 
	    mexErrMsgIdAndTxt( "MATLAB:lleval:invalidNumInputs",
                          "Two input arguments required."); 
    } else if (nlhs > 3) {
	    mexErrMsgIdAndTxt( "MATLAB:lleval:maxlhs",
                          "Too many output arguments."); 
    }   
    
    VerifyDataset(prhs[0]);

    if(mxGetM(theta_array=prhs[1])!=NFE+L+J){
        mexErrMsgIdAndTxt( "MATLAB:lleval:theta_array",
                          "Theta is wrong dimension. Should be %d", NFE+L+J); 
    }

    n = N + J + L;
    cout << "Theta dimension:" << mxGetM(theta_array) << " and n=" << n << endl;

    
    x 	    = new double [n];
    grad    = new double [n];

    memcpy(x, mxGetPr(theta_array),sizeof(double)*n);
    plhs[0] = mxCreateNumericMatrix(1,1, mxDOUBLE_CLASS, mxREAL); 
    plhs[1] = mxCreateNumericMatrix(n,1, mxDOUBLE_CLASS, mxREAL); 
    plhs[2] = mxCreateNumericMatrix(J+L,1, mxDOUBLE_CLASS, mxREAL); 
    
    output=(mxGetPr(plhs[0]));
    
    switch(nlhs)
    {

	case 0: 
	case 1: 
	output[0] = get_loglik(x);
	break;

	case 2: 
	output[0]=get_loglikg(x,grad);
    	memcpy(mxGetPr(plhs[1]), grad,sizeof(double)*n);
	break;

	case 3: 
	output[0]=get_loglikg(x,grad);
        std_err = new double [J+L];
	cout << "Starting..." << endl;
	get_stderr(x,std_err);
	cout << "VICTORY " << endl;
    	memcpy(mxGetPr(plhs[1]), grad,sizeof(double)*n);
    	memcpy(mxGetPr(plhs[2]), std_err,sizeof(double)*(J+L));
	break;
    }
		
    delete x,grad,std_err;
    return;

} 
