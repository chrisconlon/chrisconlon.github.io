#include "globals.h"
#include <iostream>
#include <armadillo>
#include <mat.h>

using namespace arma;
using namespace std;

mat Q,V,X,AV;
mat X0;
cube Z;
vec W;
uvec FE;
mwSize N,J,L,NS,M,NFE,NPARAM;
field<uvec> AVS, AV0 ;

#define DSNAME "ds3"
#define THETANAME "t0"
#define RESULTS_NAME "knitroresults.mat"

int VerifyDataset(const mxArray *mystruct){
    const mwSize *dim_array;
    mwSize elements,number_of_dims, nsubs;
    
    mxArray *array_q, *array_w, *array_v, *array_x, *array_z, *array_av, *array_mktid;
    int t, k,j;
    mwIndex *jc,*ir;
    double *qmat;
        
    /* Check data type of input argument */
    if (!(mxIsStruct(mystruct) )) {
        cout << "Input 1 must be dataset structure.";
        exit(1);
    }
    if((array_q=mxGetField(mystruct, 0, "qmat"))==NULL){
        cout << "Include the Q Matrix" ;
        exit(1);
    }
    if((array_w=mxGetField(mystruct, 0, "w"))==NULL){
        cout << "Include the quadrature weights ";
        exit(1);
    }
    if((array_v=mxGetField(mystruct, 0, "v"))==NULL){
        cout << "Include the node matrix " ;
        exit(1);
    }
    if((array_z=mxGetField(mystruct, 0, "Z"))==NULL){
        printf("Include the Z matrix (outer product of X_jl nu_il ");
        exit(1);
    }
    if((array_x=mxGetField(mystruct, 0, "Xmat"))==NULL){
        printf("Include the characteristic matrix ");
        exit(1);
    }
    if((array_av=mxGetField(mystruct, 0, "avail"))==NULL){
        printf("Include the availability matrix ");
        exit(1);
    }
    
    N = mxGetN(array_av);
    J = mxGetM(array_q)-1;
    L = mxGetN(array_x);
    dim_array=mxGetDimensions(array_z);
    NS = mxGetM(array_w);
    
    Q=zeros<mat>(J+1,N);
    AV=zeros<mat>(J,N);
    
    W=vec(mxGetPr(array_w), mxGetM(array_w));
    V=mat(mxGetPr(array_v),mxGetM(array_v),mxGetN(array_v));
    X=mat(mxGetPr(array_x),mxGetM(array_x),mxGetN(array_x));
    Z=cube(mxGetPr(array_z), dim_array[0], dim_array[1], dim_array[2]);

    if((array_mktid=mxGetField(mystruct, 0, "mktid"))==NULL){
        printf("No fixed effects included, assuming 1:1 \n");
    	NFE=N;
	FE = counting_vector(N,0);
    }else{
	vec FETemp=vec(mxGetPr(array_mktid), mxGetM(array_mktid))-1;
	FE=conv_to<uvec>::from(FETemp);
	NFE = max(FE)+1;
	cout << "Number of FE: " << NFE << " length " << FE.n_rows << endl;
    }
    
    ir=mxGetIr(array_q);
    jc=mxGetJc(array_q);
    qmat =mxGetPr(array_q);
    if(!mxIsSparse(array_q)){
	cout << "Q Matrix Must be Sparse" << endl;
	exit(1);
    }
    if(!mxIsSparse(array_av)){
	cout << "Q Matrix Must be Sparse" << endl;
	exit(1);
    }


  
   NPARAM = NFE + J + L; 
   cout << "N= " << N << ", NFE " << NFE <<  " , J= " << J << endl;

    // Nonsparse Q
    for(t=0;t<N;t++){
        for(k=jc[t]; k < jc[t+1]; k++){
            j=ir[k];
            Q(j,t)= qmat[k];
        }
    }
    cout << "Total Sales " << accu(Q) << endl;
        
    // Nonsparse AV
    ir=mxGetIr(array_av);
    jc=mxGetJc(array_av);
    for(t=0;t<N;t++){
        for(k=jc[t]; k < jc[t+1]; k++){
            j=ir[k];
            AV(j,t)= 1;
        }
    }
    
    cout << "Z has " << Z.n_cols << " columns" << endl;
    cout << "Z has " << Z.n_rows << " rows" << endl;
    cout << "Z has " << Z.n_slices << " slices" << endl;
    
    AVS=field<uvec>(N);
    AV0=field<uvec>(N);
    uvec::fixed<1> ending;
    ending(0) = J;
    
    X0 = X;
    X0.insert_rows(J,1);
    
    for(t=0; t<N;t++){
        AVS(t) = find(AV.col(t));
        AV0(t) = find(AV.col(t));
        AVS(t).insert_rows(AVS(t).n_rows,ending);
    }
}



void writemat(double *data, int n, const char *file){
    MATFile *pmat;
    mxArray *pa1;
    int status;

    printf("Creating file %s...\n\n", file);
    pmat = matOpen(file, "w");
    pa1 = mxCreateDoubleMatrix(n,1,mxREAL);
    memcpy((void *)(mxGetPr(pa1)), (void *)data, n * sizeof(double));
    status = matPutVariable(pmat, "t1", pa1);
    mxDestroyArray(pa1);
    if (matClose(pmat) != 0) {
        printf("Error closing file %s\n",file);
        return;
    }
    printf("Done\n");
    return;
}
