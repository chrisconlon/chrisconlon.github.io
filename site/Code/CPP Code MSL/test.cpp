#include <iostream>
#include <mat.h>
#include "armadillo"
#include "armahelp.h"
#include "globals.h"
#include "msllcpp.h"
#include "knitro.h"
#include "omp.h"

const double MYEPS =0.000001;

using namespace arma;
using namespace std;

#define DSNAME "ds3"
#define THETANAME "t0"
#define RESULTS_NAME "knitroresults.mat"

mxArray *openDataset(const char *file, const char *varname){
	MATFile *mymatfile;
	mxArray *pa;
    
	printf("Reading file %s...\n\n", file);
	mymatfile = matOpen(file, "r");
	if (mymatfile == NULL) {
		printf("Error reopening file %s\n", file);
		exit(1);
	}
	pa = matGetVariable(mymatfile, varname);
	if (pa == NULL) {
		printf("Error: cannot find element %s in file\n", varname);
		exit(1);
	}
    
	if (matClose(mymatfile) != 0) {
		printf("Error closing file %s\n",file);
		exit(1);
	}
	return(pa);
}

double NumericalDerivative(double *x, int k){
	double *y1,*y2;
        double g;
       
	 if(k < 0 | k > NPARAM-1){
		cout << "Cannot differentiate w.r.t. invalid parameter number" << k << endl;
		exit(1);
	}

        y1=new double [NPARAM];
	y2=new double [NPARAM];  
         
        memcpy(y1, x, sizeof(double)*NPARAM);
        memcpy(y2, x, sizeof(double)*NPARAM);

        y1[k] += MYEPS;
        y2[k] -= MYEPS;

 	g= (get_loglik(y1)-get_loglik(y2)) / (2.0 * MYEPS);

    	delete y1;
	delete y2;
	
	return g;

}


int main(int argc, char **argv)
{
    mxArray *dataset, *theta;
    double ll0, ll1;
    double *x, *grad,*std_vec;
    int i,j,n;
    
    double c0, c1, c2;
    double cputime, cputime2;
    double g2;
    double myeps = .000001;
    double *y1, *y2;

    if (argc <2){
        printf("Usage: msleval <matfile> <outfile>");
        return(0);
    }
    
    // Load dataset structure
    if (!mxIsStruct(dataset=openDataset(argv[1],"ds3"))){
        printf("Error: dataset is not valid structure \n");
        exit(1);
    }
    cout << "verifying dataset...";
    VerifyDataset(dataset);
    cout << "done!" << endl;
    
    if(mxGetM(theta=openDataset(argv[1],THETANAME))!=NPARAM){
        cout << "Parameter Vector must have NFE+K+J entries: " << NPARAM << " not " << mxGetM(theta) << endl;
        exit(1);
    }
    
    cout << "n parameters " << NPARAM << endl;
    x = new double [NPARAM];
    grad = new double [NPARAM];
    std_vec= new double [NPARAM];

    memcpy(x, mxGetPr(theta), sizeof(double)*NPARAM);

    double TQ = as_scalar(accu(Q));
    printf("Total Sales: %8.4f\n",TQ);

    c0=omp_get_wtime();
    ll0 = get_loglik(x);
    c1=omp_get_wtime();
    printf("Log Likelihood: %8.4f in %6.4f seconds\n",ll0, (c1-c0));
    ll0 = get_loglikg(x,std_vec);
    get_loglikg(x,grad);
    c2=omp_get_wtime();
    cout << "function AND gradient in " << (c2 - c1) << " seconds" << endl;


    c1=omp_get_wtime();
    get_stderr(x,std_vec);
    c2=omp_get_wtime();
    cout << "Hessian in " << (c2 - c1) << " seconds" << endl;

/*
    for(i=0;i< J+L+5;i++){
	g2=NumericalDerivative(x,i);
        cout << "grad " << i+1 <<":\t"  << grad[i] << "\t" << g2<< endl;
    }
	for(i = NFE+J+L-2; i < NFE+J+L; i++){
	g2=NumericalDerivative(x,i);
        cout << "grad " << i+1 << ":\t"  << grad[i] << "\t" << g2<< endl;
    }
*/
    if(argc==3){
    	writemat(x,NPARAM,argv[2]);
    }
//    delete grad;
//    delete x;
}



