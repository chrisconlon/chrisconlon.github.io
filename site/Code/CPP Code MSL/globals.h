#ifndef GLOBAL_H
#define GLOBAL_H

#define LL_GRAD 1
#define LL_ONLY 0

#include "armadillo"
#include <mat.h>
#include "armahelp.h"

using namespace arma;

extern mat Q,V,X,AV;
extern mat X0;
extern cube Z;
extern vec W;
extern uvec FE;
extern mwSize N,J,L,NS,M,NFE,NPARAM;
extern field<uvec> AVS, AV0 ;
extern vec totalSales;

int VerifyDataset(const mxArray *mystruct);
void writemat(double *data, int n, const char *file);

#endif
