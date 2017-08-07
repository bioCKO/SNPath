/* logistic-ridge.c
 * 3/2008
 * C module called by an R wrapper function for logistic
 * ridge algorithm.
 */

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

#define SIGN(a) ((a)>0?1:((a)<0?-1:0))

double compute_loglik(double *X, double *Y, double *W, double *BETA, int *N, int *P);

double compute_pllik(double *X, double*Y, double *W, double *BETA, int *N, int *P, int *INTERCEPT, int *GRP, int *NGRP, int *UG, double *LAMBDA, int *EQSNP){
    double temp1, temp2, temp3, llik;
    int g,j;
    temp1 = 0.0;
    llik = - compute_loglik(X, Y, W, BETA, N, P);
	for (g=*INTERCEPT; g<*NGRP; g++){
		temp2 = 0.0;
		temp3 = 0.0;
		for (j=0; j<*P; j++){
			if (GRP[j] == UG[g]){
				temp3 += 1.0;
				temp2 += fabs(BETA[j]);
			}
		}
		temp2 *= temp2;
		if (*EQSNP==1) {
			temp1 +=temp2;
    	} else{
			temp1 += sqrt(temp3)*temp2;
		}
 	}
	llik += *LAMBDA*temp1;

    return llik;
}


double compute_beta0(double *X, double *Y, double *W, double *BETA, int *N, int *P, int *INTERCEPT, int *GRP, int *NGRP, int *UG, double *LAMBDA, int *EQSNP, double * TOL){
	double beta0, dStep, oldLlik, newLlik;
	double dDenominator, dNumerator, delta, diff;
	double temp1, temp3;
	double tol;
	int i,j,n,p;
	n = *N;
	double temp2[n];
	p = *P;
	tol =*TOL;
    dStep=0.0;
    delta = 1.0;
	newLlik = compute_pllik(X, Y, W, BETA, N, P, INTERCEPT, GRP, NGRP, UG, LAMBDA, EQSNP);

	for (i=0; i<n; i++) {
		temp2[i]=0.0;
        for (j=1; j<p; j++)
            temp2[i] += X[i*p+j]*BETA[j];
        temp2[i] *= Y[i];
    }
	do {
	    oldLlik = newLlik;

        dNumerator = 0.0;
        dDenominator = 0.0;
	    for (i=0; i<n; i++){
	        dNumerator += W[i]*Y[i]/(1+exp(BETA[0]*Y[i]+temp2[i]));

        diff = fabs(BETA[0]*Y[i]+temp2[i])-fabs(delta);
		if (diff<=0.0)
			temp1 = 0.25;
		else
		    temp1 = 1.0/(2.0+exp(diff)+exp(-diff));
			dDenominator += W[i]*temp1;
		}
        dStep = dNumerator/dDenominator;

        if (dStep < -delta)
			dStep = -delta;
		else if (dStep > delta)
			dStep = delta;

	    BETA[0] += dStep;
	    temp1 = 2.0*fabs(dStep);
		temp3 = delta/2.0;
		delta = (temp1>temp3?temp1:temp3);

    	newLlik = compute_pllik(X, Y, W, BETA, N, P, INTERCEPT, GRP, NGRP, UG, LAMBDA, EQSNP);


    } while (fabs(dStep) >= 1.490116e-8*fabs(BETA[0])+tol/3 );
    beta0 = BETA[0];
	return beta0;
}


void compute_loglik_c(double *X, double *Y, double *W, double *BETA, int*N, int *P, double *LLIK){
	*LLIK = compute_loglik(X, Y, W, BETA, N, P);
}

/* logistic.ridge.c
	X - input data 0/1 vector of length N*P
	N - number of observations
	P - number of variables
	LAMBDA - ridge penalty
	TOL - tolerance for termination
	MAXITER - maximum number of iterations (unlimited if MAXITER<=0)
	BETA - output beta
	ITER - number of iterations performed before convergence */
void logistic_ridge_c(double *X, double *Y, double *W, int * GRP, int * NGRP, int *UG, int *N, int *P, double
*LAMBDA, double *TOL, int *MAXITER, int *INTERCEPT, int *EQSNP, double *BETA, double *LLIK) {
  int n, p, nMaxIter;
  int ngrp;
  double lambdag, tol;
  /*double *delta, *b2;/* step interval bounds */
  double dOldLlik,dNewLlik, gOldLlik, gNewLlik;
  double dStep;  /* \Delta */
  double dNumerator;   /* numerator of dStep */
  double dDenominator; /* denominator of dStep */
  int i,j,k,g,iter,giter;
  double diff,temp1,temp2, temp3;
  double Bj;

  /* initialization */
  p = *P;
  n = *N;
  ngrp = *NGRP;

  tol = *TOL;
  nMaxIter = (*MAXITER>0?*MAXITER:INT_MAX);
  iter = 0;
  giter =0;

  double b2[p]; /*(double *) malloc(p*1*sizeof(double));*/
  double delta[p]; /*(double *) malloc(p*1*sizeof(double));*/
  double kg[ngrp];

  double dotprod[n];


  for (i=0; i<p; i++)
	BETA[i]=0.0;

  kg[0]=0.0;
  for (g=*INTERCEPT; g<ngrp; g++){
	  if (*EQSNP==1) {
		  kg[g]=1.0;
 	  } else {
	    kg[g]=0.0;
	    for (i=*INTERCEPT;i<p;i++)
		  if(GRP[i]==UG[g]) kg[g]=kg[g]+1.0;
  }}

  /* main loop */
  dNewLlik = compute_pllik(X, Y, W, BETA, N, P, INTERCEPT, GRP, NGRP, UG, LAMBDA, EQSNP);

  do {

	dOldLlik = dNewLlik;

    for (i=0; i<p; i++){
		delta[i]=1.0;
    }

    if (*INTERCEPT==1) BETA[0]= compute_beta0(X, Y, W, BETA, N, P, INTERCEPT, GRP, NGRP, UG, LAMBDA, EQSNP, TOL);



    for (i=0;i<n;i++){
		dotprod[i]=0.0;
		for (k=0; k<p; k++) {
			dotprod[i] += X[i*p+k]*BETA[k];
		}
	    dotprod[i] *= Y[i];
	}



    for (g=*INTERCEPT; g<ngrp; g++){

	  lambdag = *LAMBDA * sqrt(kg[g]);


      for (i=p; i>=0; i--)
  	       b2[i]=BETA[i];

      gNewLlik = compute_pllik(X, Y, W, b2, N, P, INTERCEPT, GRP, NGRP, UG, LAMBDA, EQSNP);

      do {

    	gOldLlik = gNewLlik;

	    for (j=0; j<p; j++){
		  dStep =0.0;
		  dNumerator=0.0;
		  dDenominator= 2.0* lambdag;



		  if (GRP[j]== UG[g]){
		  /* compute tentative step */

          Bj=0.0;
          for (i=0; i<p; i++) {
  	           if (GRP[i]==UG[g]) Bj += fabs(b2[i]);
	      }

		  for (i=0; i<n; i++) {

			  dNumerator += W[i]*Y[i]*X[i*p+j]/(1+exp(dotprod[i]));

			  diff = fabs(dotprod[i])-fabs(delta[j] * X[i*p+j]);
			  if (diff<=0.0)
			    temp1 = 0.25;
			  else
			    temp1 = 1.0/(2.0+exp(diff)+exp(-diff));
			  dDenominator += W[i]*X[i*p+j]*X[i*p+j] * temp1;
		    }

		    if (b2[j]==0.0) {
			  dStep = dNumerator/dDenominator - 2.0*Bj*lambdag/dDenominator; /* L1 penalty */
			  if (dStep<=0.0){
				  dStep = dNumerator/dDenominator + 2.0*Bj*lambdag/dDenominator; /* L1 penalty */
				  if (dStep>=0.0)
				    dStep = 0.0;
			  }
		    } else {
			  temp1 = 2*Bj*lambdag/dDenominator; /* L1 penalty */
			  dStep = dNumerator/dDenominator-SIGN(b2[j])*temp1;
			  if ((dStep+b2[j])*SIGN(b2[j]) < 0.0)
			    dStep = -b2[j];    /* negative direction */
		    }

		    /* update beta and delta */
	        if (dStep < -delta[j])
	              dStep = -delta[j];
	        else if (dStep > delta[j])
	              dStep = delta[j];

		    b2[j] += dStep;
		    temp1 = 2.0*fabs(dStep);
		    temp2 = delta[j]/2.0;
		    delta[j] = (temp1>temp2?temp1:temp2);


	        for (i=0;i<n;i++){
					dotprod[i] += Y[i]*X[i*p+j]*dStep;
			}
          }
	    }

        gNewLlik = compute_pllik(X, Y, W, b2, N, P, INTERCEPT, GRP, NGRP, UG, LAMBDA, EQSNP);

      } while ((++ giter<nMaxIter) && (fabs(gNewLlik-gOldLlik) >= tol*fabs(gOldLlik)));

      for (i=0; i<p; i++)	 BETA[i]=b2[i];

    }

    dNewLlik = compute_pllik(X, Y, W, BETA, N, P, INTERCEPT, GRP, NGRP, UG, LAMBDA, EQSNP);

  } while ((++ iter<nMaxIter) && (fabs(dNewLlik-dOldLlik) >= tol*fabs(dOldLlik)));


   *LLIK = dNewLlik;

}


double compute_loglik(double *X, double *Y, double *W, double *BETA, int *N, int *P) {
  int n,p;
  n = *N;
  p = *P;
  double llik=0.0;
  double sum;
  int i,j;
  for (i=0; i<n; i++) {
	sum = 0.0;
	for (j=0; j<p; j++){
	  sum += X[i*p+j]*BETA[j];
	}
	llik -= W[i]*log(exp(-Y[i]*sum)+1);
  }
  return llik;
}

