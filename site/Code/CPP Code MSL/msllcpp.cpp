#include <iostream>
#include "armadillo"

#include "globals.h"
#include "msllcpp.h"
#include "armahelp.h"

#define MYEPS .001

using namespace arma;
using namespace std;

mat create_utils(vec theta){
    vec dj      =theta.subvec( 0, J-1 );
    vec sigma   = theta.subvec(J,J+L-1);
    mat utils=mat(J,NS);        

#pragma omp parallel
{
#pragma omp for
    for(int i=0; i<NS;i++){
        utils.col(i) = exp(trans(Z.slice(i)) * sigma +dj);
    }
}
    return utils;
}


vec parallel_likelihood(const mat utils, const vec theta, double &ll, int mode){
    // extract parameters
    vec dj      =theta.subvec( 0, J-1 );
    vec sigma   =theta.subvec(J,J+L-1);
    vec xi      =theta.subvec(J+L,J+L+NFE-1);
    mat scores = zeros<mat>(N,J+1+L);
    double    ll0=0.0;
    vec newxi = xi.elem(FE);

    if (newxi.n_elem != N){
	cout << "Expanded FE vector of wrong length" << newxi.n_elem << " instead of " << N << endl;
    }
        
#pragma omp parallel
{
#pragma omp for reduction(+:ll0)
    //loop over markets
        for(int t=0; t<N; t++){
            mat pij=zeros<mat>(J+1,NS);
            vec qt = Q.col(t);
            
            // construct pij matrix for each t
            for(int i=0;i<NS;i++){
                double denom=1.0/(exp(newxi(t))+dot(AV.col(t),utils.col(i))); 
                vec pvec = (AV.col(t) % utils.col(i))*denom;                      
                pij.col(i).subvec(0,J-1) = pvec;
                pij(J,i) = 1.0-sum(pvec);
            }
            // compute p_j by integration and compute likelihood
            ll0 += as_scalar((qt.t()* pos_log(pij * W)));
	    if(mode==LL_GRAD){
		scores.row(t) = score_contribution(pij, qt);
	    }
        }
}
    ll = ll0;
    if(mode==LL_GRAD){
	return compute_gradient(scores);
    }else{
	return zeros<vec>(NFE+J+L);
    }
}

rowvec score_contribution(const mat pij, const vec qvec)
{
    rowvec grad_dj=zeros<rowvec>(J+1);
    rowvec grad_sigma=zeros<rowvec>(L);
    vec qw = element_divide(qvec,pij * W);
    
    // Loop over draws --do outer product of pijt's and sum up
    for(int i=0;i<NS;i++){
        vec p1 = pij.col(i);
        grad_dj -= trans(qw) * (p1 * trans(p1))  *W(i);

        // Loop over sigma_i characteristics
        for(int l=0; l < L; l++){
            double mysum = as_scalar(trans(p1) * X0.col(l));
            for(int j=0; j<J+1; j++){
                grad_sigma(l) += W(i) * qw(j)* p1(j) * V(i,l) * (X0(j,l) - mysum);
            }
        }
    }
    grad_dj += qvec.t();
    return join_rows(grad_dj,grad_sigma);
}

mat probs_only(const mat utils, const vec theta){
    // extract parameters
    vec dj      =theta.subvec( 0, J-1 );
    vec sigma   =theta.subvec(J,J+L-1);
    vec xi      =theta.subvec(J+L,J+L+N-1);
    vec newxi = xi.elem(FE);
    mat pjt = zeros<mat>(J+1,N);
    
#pragma omp parallel
    {
#pragma omp for
        //loop over markets
        for(int t=0; t<N; t++){
            mat pij=zeros<mat>(J+1,NS);
            vec qt = Q.col(t);
            
            // construct pij matrix for each t
            for(int i=0;i<NS;i++){
                double denom=1.0/(exp(newxi(t))+dot(AV.col(t),utils.col(i))); 
                vec pvec = (AV.col(t) % utils.col(i))*denom;                      
                pij.col(i).subvec(0,J-1) = pvec;
                pij(J,i) = 1.0-sum(pvec);
            }
            
            // compute p_j by integration and compute likelihood
            vec pj = pij * W;
            pjt.col(t) = pj;
        }
    }
    return pjt;
}
double get_loglik(const double *theta)
{
    double ll;
    vec x=vec(theta,NFE+L+J);
    mat utils=create_utils(x);
    
    parallel_likelihood(utils,x,ll,LL_ONLY);
    return ll;
}

double get_loglikg(const double *theta, double *grad)
{
    double ll;
    vec x=vec(theta,NFE+L+J);
    mat utils=create_utils(x);

    vec gradient=parallel_likelihood(utils,x,ll,LL_GRAD);
    memcpy(grad,gradient.memptr(),sizeof(double)*(NFE+J+L));

    return ll;
}

double get_stderr(const double *theta, double *sevec)
{
    double ll;
    vec x=vec(theta,NFE+L+J);
    mat utils=create_utils(x);
    vec std_vector=compute_se(utils,x);
	cout << std_vector << endl;
    memcpy(sevec,std_vector.memptr(),sizeof(double)*(NFE+J+L));
    return ll;
}

mat get_probs(const double *theta)
{
    
    vec x=vec(theta,NFE+L+J);
    mat utils=create_utils(x);    
    return probs_only(utils,x);
}

vec compute_gradient(const mat scores){
    vec temp_vec = trans(sum(scores));
    vec grad_dj = temp_vec.subvec(0,J-1);
    vec grad_sigma = temp_vec.subvec(J+1,J+L);
    vec grad_xi2 = accu_vec(FE,scores.col(J));
    return join_cols(join_cols(grad_dj,grad_sigma),grad_xi2);
}

mat score_contribution2(const mat pij, const vec qvec, mat &A, vec &B)
{
    mat scores;
    mat score_dj=zeros<mat>(J+1,J+1);
    mat score_sigma=zeros<mat>(J+1,L);
    mat BHHH=zeros<mat>(J+L+1,J+L+1);
    uvec sold = find(qvec > MYEPS);
    vec pj = pij * W;
    
    // Loop over draws --do outer product of pijt's and sum up
    for(int i=0;i<NS;i++){
        vec p1 = pij.col(i);
        score_dj -=  (p1 * trans(p1))  *W(i);
 	for(int l=0; l < L; l++){
           double mysum = as_scalar(trans(p1) * X0.col(l));
     	   for(int ii=0; ii< sold.n_elem; ii++){
     		int j=sold(ii);
                score_sigma(j,l) += W(i) * p1(j) * V(i,l) * (X0(j,l) - mysum);
            }
        }
    }

    //fix own products on diagonal
    score_dj.diag() += pj;

     // normalize scores by 1/pj for entries with positive sales
     for(int i=0; i< sold.n_elem; i++){
     	int j=sold(i);
	score_dj.row(j) /= pj(j);
	score_sigma.row(j) /= pj(j);
     }

     scores=join_rows(score_dj,score_sigma);


     // construct outer product of scores
     // only use contributions with positive sales 
     for(int i=0; i< sold.n_elem; i++){
     	int j=sold(i);
	vec myvec = trans(scores.row(j));
	BHHH += trimatu(myvec * myvec.t() * qvec(j));
      }

     // extract xi_t components
     vec xi2=BHHH.col(J).subvec(0,J-1);
     vec xi1=trans(BHHH.row(J).subvec(J+1,J+L));
     vec xipart=join_cols(xi2, xi1);
     xipart.resize(J+L+1);
     xipart(J+L) = BHHH.at(J,J);

     // remove xi_t from the Hessian
     BHHH.shed_row(J);
     BHHH.shed_col(J);

     A=BHHH;
     B=xipart;
}

vec compute_se(const mat utils, const vec theta){
    // extract parameters
    vec dj      =theta.subvec( 0, J-1 );
    vec sigma   =theta.subvec(J,J+L-1);
    vec xi	=theta.subvec(J+L,J+L+NFE-1);
    vec newxi = xi.elem(FE);
   
    mat HA = zeros<mat>(J+L,J+L);
    mat HB = zeros<mat>(J+L,N);
    vec D = zeros<vec>(N);

    if (newxi.n_elem != N){
        cout << "Expanded FE vector of wrong length" << newxi.n_elem << " instead of " << N << endl;
    }

    //loop over markets
 	for(int t=0; t<N; t++){
            mat pij=zeros<mat>(J+1,NS);
            vec qt = Q.col(t); 
            // construct pij matrix for each t
            for(int i=0;i<NS;i++){
                double denom=1.0/(exp(newxi(t))+dot(AV.col(t),utils.col(i)));
         	vec pvec = (AV.col(t) % utils.col(i))*denom;
                pij.col(i).subvec(0,J-1) = pvec;
                pij(J,i) = 1.0-sum(pvec);
            }
            // Compute the Hessian Contribution of each market
   	   mat A = zeros<mat>(J+L,J+L);
           vec B = zeros<vec>(J+L+1);
	   score_contribution2(pij, qt,A,B);
           HA+=A;
           D(t)=1.0/B(J+L);
           HB.col(t) = B.subvec(0,J+L-1);
        }
        mat HA2= symmatu(HA);
	mat part2=HB*diagmat(D)*HB.t();
	vec stderr=sqrt(diagvec(inv(symmatu(HA2+part2))));	
	return stderr;
}



double sum_FE(const double *theta)
{
        vec x = vec(theta,NFE+L+J);
        double mysum=sum(x.subvec(J+L,J+L+NFE-1))/NFE;
        return mysum;
}

void get_jacobian(double *jac)
{
        vec y = zeros<vec>(J+L+N,1);
        for(int i=J+L; i<J+L+NFE; i++){
                y(i) = 1.0/NFE;
        }
        memcpy(jac,y.memptr(),sizeof(double)*(NFE+J+L));
}

















