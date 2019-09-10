#include <iostream>
#include <mat.h>
#include "armadillo"
#include "armahelp.h"
#include "globals.h"
#include "msllcpp.h"
#include "knitro.h"

using namespace arma;
using namespace std;


#define DSNAME "ds3"
#define THETANAME "t0"

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

void getProblemSizes (int *const n, int *const m, int *const nnzJ, int *const nnzH)
{
    *n = J + L + NFE ;
    *m = 0;
    *nnzJ = 0;
    *nnzH = 0;
    
    return;
}

void  getProblemData (int    * const  objType,       /*-- SCALAR */
                      int    * const  objGoal,       /*-- SCALAR */
                      double * const  xLoBnds,       /*-- ARRAY LENGTH n */
                      double * const  xUpBnds,       /*-- ARRAY LENGTH n */
                      double * const  xInitial,      /*-- ARRAY LENGTH n */
                      int    * const  cType,         /*-- ARRAY LENGTH m */
                      double * const  cLoBnds,       /*-- ARRAY LENGTH m */
                      double * const  cUpBnds,       /*-- ARRAY LENGTH m */
                      int    * const  jacIndexVars,  /*-- ARRAY LENGTH nnzJ */
                      int    * const  jacIndexCons,  /*-- ARRAY LENGTH nnzJ */
                      int    * const  hessRows,      /*-- ARRAY LENGTH nnzH */
                      int    * const  hessCols,      /*-- ARRAY LENGTH nnzH */
                      int    * const  objFnType,     /*-- MIP SCALAR */
                      int    * const  xType,         /*-- MIP ARRAY LENGTH n */
                      int    * const  cFnType)      /*-- MIP ARRAY LENGTH m */
{
    
    int i;
    *objType = 0;
    *objGoal = KTR_OBJGOAL_MAXIMIZE;

     for(i=0; i<J; i++){
        xUpBnds[i] =  0.0;
        xLoBnds[i] =  -20.0;
     }
     for(i=J; i<J+L; i++){
        xUpBnds[i] =  10.0;
        xLoBnds[i] =  0.0;
     }
     for(i=J+L; i<J+L+NFE; i++){
        xUpBnds[i] =  30.0;
        xLoBnds[i] =  -30.0;
     }

//xLoBnds[J+L+5007-1] = 1.0;
//xUpBnds[J+L+5007-1] = 1.0;

xLoBnds[J+L] = 0.0;
xUpBnds[J+L] = 0.0;


return;
}


int callbackEvalGA (const int evalRequestCode,
                    const int n,
                    const int m,
                    const int nnzJ,
                    const int nnzH,
                    const double *const x,
                    const double *const lambda,
                    double *const obj,
                    double *const c,
                    double *const objGrad,
                    double *const jac,
                    double *const hessian,
                    double *const hessVector, void *userParams)
{
    int i;
    
    if (evalRequestCode != KTR_RC_EVALGA)
    {
        printf ("*** callbackEvalGA incorrectly called with eval code %d\n",
                evalRequestCode);
        return (-1);
    }
    
    *obj = get_loglikg(x, objGrad);
    return (0);
}  

int callbackEvalFC (const int evalRequestCode,
                    const int n,
                    const int m,
                    const int nnzJ,
                    const int nnzH,
                    const double *const x,
                    const double *const lambda,
                    double *const obj,
                    double *const c,
                    double *const objGrad,
                    double *const jac,
                    double *const hessian,
                    double *const hessVector, void *userParams)
{
    if (evalRequestCode != KTR_RC_EVALFC)
    {
        printf ("*** callbackEvalFC incorrectly called with eval code %d\n",
                evalRequestCode);
        return (-1);
    } 
    *obj = get_loglik (x);
    return (0);
}     

int main(int argc, char **argv)
{
    int result,i;
    
    int  nStatus;
    int  nHessOpt;
    
    KTR_context  *kc;
    int          n, m, nnzJ, nnzH, objGoal, objType;
    int          *cType;
    int          *jacIndexVars, *jacIndexCons, *hessRows, *hessCols;
    double       obj, *x, *lambda;
    double       *xLoBnds, *xUpBnds, *xInitial, *cLoBnds, *cUpBnds;
    
    mxArray *dataset, *theta;
	
    if (argc <2){
		printf("Usage: msleval <matfile> <outfile>");
		return(0);
	}
    
    // Load dataset structure
    if (!mxIsStruct(dataset=openDataset(argv[1],"ds3"))){
		printf("Error: dataset is not valid structure \n");
		exit(1);
	}
     VerifyDataset(dataset);
    
    // Load starting values for parameters
    if(mxGetM(theta=openDataset(argv[1],THETANAME))!=NFE+L+J){
        cout << "Parameter Vector must have NFE+K+J entries"<< endl;
        exit(1);
    }
    

    getProblemSizes (&n, &m, &nnzJ, &nnzH);

    x 		= new double[n];
    xLoBnds 	= new double[n];
    xUpBnds 	= new double[n];
    xInitial 	= new double[n];
    cType 	= new int[m];
    cLoBnds 	= new double[m];
    cUpBnds 	= new double[m];

    jacIndexVars = new int[nnzH];
    jacIndexCons = new int[nnzJ];
    hessRows     = new int[nnzH];
    hessCols     = new int[nnzH]; 
    
    cout << n << " parameters " << endl;

    getProblemData (&objType, &objGoal, xLoBnds, xUpBnds, xInitial,
                    cType, cLoBnds, cUpBnds,
                    jacIndexVars, jacIndexCons,
                    hessRows, hessCols, NULL, NULL,NULL);
    
    
    xInitial=mxGetPr(theta);
    memcpy(x, xInitial, sizeof(double)*n);    
    lambda = (double *) malloc ((m+n) * sizeof(double));
    cout << "Log likelihood: "<< get_loglik(x) << endl;

    kc = KTR_new();
    if (kc == NULL)
    {
        cout << "Failed to find a Ziena license." << endl;
        return( -1 );
    }
    
    if (KTR_load_param_file (kc, "knitro.opt") != 0)
        exit( -1 );
    if (KTR_get_int_param_by_name (kc, "hessopt", &nHessOpt) != 0)
        exit( -1 );
    
    if (KTR_set_func_callback (kc, &callbackEvalFC) != 0)
        exit( -1 );
	if (KTR_set_grad_callback (kc, &callbackEvalGA) != 0)
        exit( -1 );
    nStatus = KTR_init_problem (kc, n, objGoal, objType,
                                xLoBnds, xUpBnds,
                                m, cType, cLoBnds, cUpBnds,
                                nnzJ, jacIndexVars, jacIndexCons,
                                nnzH, hessRows, hessCols, xInitial, NULL);
/*    
    delete[] xLoBnds;
    delete[] xUpBnds;
    delete[] xInitial;
    delete[] cType;
    delete[] cLoBnds;
    delete[] cUpBnds;
    delete[] jacIndexVars;
    delete[] jacIndexCons;
    delete[] hessRows;
    delete[] hessCols;
*/    
    nStatus = KTR_solve (kc, x, lambda, 0, &obj, NULL, NULL, NULL, NULL, NULL, NULL);
    
    printf ("\n\n");
    if (nStatus != 0){
        printf ("KNITRO failed to solve the problem, final status = %d\n",nStatus);
        writemat(x,n,argv[2]);
    }else{
        printf ("KNITRO successful, feasibility violation    = %e\n",
                KTR_get_abs_feas_error (kc));
        printf ("                   KKT optimality violation = %e\n",
                KTR_get_abs_opt_error (kc));
        writemat(x,n,argv[2]);
    }
    KTR_free (&kc);
    
    delete[] x;
    delete[] lambda;
     
    return 0;  
}
