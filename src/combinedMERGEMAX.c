
#include <math.h>
#include "survS.h"
#include "survproto.h"

/* This code is used in the merging process in the Boomerang. A large part of this code was taken from C code available in the R CRAN package 'survival' on November 19, 2014. I have modified the code for use in the Boomerang algorithm. For the freely available Boomerang library, please visit www.github.com/sedevlin.
 */


double cpeftnM(double *time, int *status, double *xvec, double *currentlevels, int nused, int nvar,
                  double *offset, double *weights, int *strata,double eps,double toler) {
   
    int i,j,k, person;
    double  **cmat, **imat;  /*ragged arrays */
    double  wtave;
    double *a, *newbeta;
    double *a2, **cmat2;
    double *scale;
    double  denom=0, zbeta, risk;
    double  temp, temp2;
    int     ndead;  /* number of death obs at a time point */
    double  tdeath=0;  /* ndead= total at a given time point, tdeath= all */

    double  newlk=0;
    double  dtime, d2;
    double  deadwt;  /*sum of case weights for the deaths*/
    double  efronwt; /* sum of weighted risk scores for the deaths*/
    int     halving;    /*are we doing step halving at the moment? */
    int     nrisk;   /* number of subjects in the current risk set */
    double  *maxbeta;
    
    /** hard selection of a few features*/
    int method=1;
    int doscale=1;
    int maxiter = 20;

    /* returned objects */
    SEXP imat2, means2, u2, loglik2;
    double *u, *loglik, *means;
    SEXP sctest2, flag2, iter2;
    double *sctest;
    int *flag, *iter;
    int nprotect;  /* number of protect calls I have issued */

    nprotect =0;
 
    SEXP beta2;
    PROTECT(beta2= allocVector(REALSXP, nvar));
    nprotect++;
    double *beta;
    beta = REAL(beta2);
    for (i=0; i< nvar; i++) {
        beta[i]=0;
    }
    
    SEXP covar2;
    double **covar;
    covar2 =  PROTECT(allocVector(REALSXP, nused*nvar));
    covar= dmatrix(REAL(covar2), nused, nvar);
   	nprotect++;
 
    
    SEXP covart2;
    double **covart;
    covart2 =  PROTECT(allocVector(REALSXP, nused*nvar));
    covart= dmatrix(REAL(covart2), nused, nvar);
   	nprotect++;
    
    
    int markers;
    for (person=0; person<nused; person++) {
        for (markers=0; markers<nvar; markers++) {
            if(xvec[person] == currentlevels[markers]) covar[markers][person]=1;
            if(xvec[person] != currentlevels[markers]) covar[markers][person]=0;
            if(xvec[person] == currentlevels[markers]) covart[markers][person]=1;
            if(xvec[person] != currentlevels[markers]) covart[markers][person]=0;
        }
    }

    /* end of addin    */
    SEXP lpout2;
    double *lpout;
    lpout2 =    PROTECT(allocVector(REALSXP, nused));
    lpout = REAL(lpout2);
    nprotect++;
    
    
    PROTECT(imat2 = allocVector(REALSXP, nvar*nvar));
    nprotect++;
    imat = dmatrix(REAL(imat2),  nvar, nvar);
    a = (double *) R_alloc(2*nvar*nvar + 5*nvar, sizeof(double));
    newbeta = a + nvar;
    a2 = newbeta + nvar;
    maxbeta = a2 + nvar;
    scale = maxbeta + nvar;
    cmat = dmatrix(scale + nvar,   nvar, nvar);
    cmat2= dmatrix(scale + nvar +nvar*nvar, nvar, nvar);

    /* 
    ** create output variables
    */ 
  
    PROTECT(means2 = allocVector(REALSXP, nvar));
    means = REAL(means2);
    PROTECT(u2 = allocVector(REALSXP, nvar));
    u = REAL(u2);
    PROTECT(loglik2 = allocVector(REALSXP, 2)); 
    loglik = REAL(loglik2);
    PROTECT(sctest2 = allocVector(REALSXP, 1));
    sctest = REAL(sctest2);
    PROTECT(flag2 = allocVector(INTSXP, 1));
    flag = INTEGER(flag2);
    PROTECT(iter2 = allocVector(INTSXP, 1));
    iter = INTEGER(iter2);
    nprotect += 6;
    
    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable.
    */
    tdeath=0; temp2=0;
    for (i=0; i<nused; i++) {
	temp2 += weights[i];
	tdeath += weights[i] * status[i];
    }	
    for (i=0; i<nvar; i++) {
	temp=0;
	for (person=0; person<nused; person++) 
	    temp += weights[person] * covar[i][person];
	temp /= temp2;
	means[i] = temp;
	for (person=0; person<nused; person++) covar[i][person] -=temp;
	if (doscale==1) {  /* and also scale it */
	    temp =0;
	    for (person=0; person<nused; person++) {
		temp += weights[person] * fabs(covar[i][person]);
	    }
	    if (temp > 0) temp = temp2/temp;   /* scaling */
	    else temp=1.0; /* rare case of a constant covariate */
	    scale[i] = temp;
	    for (person=0; person<nused; person++)  covar[i][person] *= temp;
	    }
	}
    if (doscale==1) {
	for (i=0; i<nvar; i++) beta[i] /= scale[i]; /*rescale initial betas */
	}
    else {
	for (i=0; i<nvar; i++) scale[i] = 1.0;
	}

    /*
    ** do the initial iteration step
    */
    strata[nused-1] =1;
    loglik[1] =0;
    for (i=0; i<nvar; i++) {
	u[i] =0;
	a2[i] =0;
	for (j=0; j<nvar; j++) {
	    imat[i][j] =0 ;
	    cmat2[i][j] =0;
	    }
	}

    for (person=nused-1; person>=0; ) {
	if (strata[person] == 1) {
	    nrisk =0 ;  
	    denom = 0;
	    for (i=0; i<nvar; i++) {
		a[i] = 0;
		for (j=0; j<nvar; j++) cmat[i][j] = 0;
		}
	    }

	dtime = time[person];
	ndead =0; /*number of deaths at this time point */
	deadwt =0;  /* sum of weights for the deaths */
	efronwt=0;  /* sum of weighted risks for the deaths */
	while(person >=0 &&time[person]==dtime) {
	    /* walk through the this set of tied times */
	    nrisk++;
	    zbeta = offset[person];    /* form the term beta*z (vector mult) */
	    for (i=0; i<nvar; i++)
		zbeta += beta[i]*covar[i][person];
	    risk = exp(zbeta) * weights[person];
	    denom += risk;

	    /* a is the vector of weighted sums of x, cmat sums of squares */
	    for (i=0; i<nvar; i++) {
		a[i] += risk*covar[i][person];
		for (j=0; j<=i; j++)
		    cmat[i][j] += risk*covar[i][person]*covar[j][person];
	        }

	    if (status[person]==1) {
		ndead++;
		deadwt += weights[person];
		efronwt += risk;
		loglik[1] += weights[person]*zbeta;

		for (i=0; i<nvar; i++) 
		    u[i] += weights[person]*covar[i][person];
		if (method==1) { /* Efron */
		    for (i=0; i<nvar; i++) {
			a2[i] +=  risk*covar[i][person];
			for (j=0; j<=i; j++)
			    cmat2[i][j] += risk*covar[i][person]*covar[j][person];
		        }
		    }
	        }
	    
	    person--;
	    if (strata[person]==1) break;  /*ties don't cross strata */
	    }


	if (ndead >0) {  /* we need to add to the main terms */
	    if (method==0) { /* Breslow */
		loglik[1] -= deadwt* log(denom);
	   
		for (i=0; i<nvar; i++) {
		    temp2= a[i]/ denom;  /* mean */
		    u[i] -=  deadwt* temp2;
		    for (j=0; j<=i; j++)
			imat[j][i] += deadwt*(cmat[i][j] - temp2*a[j])/denom;
		    }
		}
	    else { /* Efron */
		/*
		** If there are 3 deaths we have 3 terms: in the first the
		**  three deaths are all in, in the second they are 2/3
		**  in the sums, and in the last 1/3 in the sum.  Let k go
		**  from 0 to (ndead -1), then we will sequentially use
		**     denom - (k/ndead)*efronwt as the denominator
		**     a - (k/ndead)*a2 as the "a" term
		**     cmat - (k/ndead)*cmat2 as the "cmat" term
		**  and reprise the equations just above.
		*/
		for (k=0; k<ndead; k++) {
		    temp = (double)k/ ndead;
		    wtave = deadwt/ndead;
		    d2 = denom - temp*efronwt;
		    loglik[1] -= wtave* log(d2);
		    for (i=0; i<nvar; i++) {
			temp2 = (a[i] - temp*a2[i])/ d2;
			u[i] -= wtave *temp2;
			for (j=0; j<=i; j++)
			    imat[j][i] +=  (wtave/d2) *
				((cmat[i][j] - temp*cmat2[i][j]) -
					  temp2*(a[j]-temp*a2[j]));
		        }
		    }
		
		for (i=0; i<nvar; i++) {
		    a2[i]=0;
		    for (j=0; j<nvar; j++) cmat2[i][j]=0;
		    }
		}
	    }
	}   /* end  of accumulation loop */
    loglik[0] = loglik[1]; /* save the loglik for iter 0 */
    
    /*
    ** Use the initial variance matrix to set a maximum coefficient
    **  (The matrix contains the variance of X * weighted number of deaths)
    */
    for (i=0; i<nvar; i++) 
	maxbeta[i] = 20* sqrt(imat[i][i]/tdeath);

    /* am I done?
    **   update the betas and test for convergence
    */
    for (i=0; i<nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
	a[i] = u[i];

    *flag= cholesky2(imat, nvar, toler);
    chsolve2(imat,nvar,a);        /* a replaced by  a *inverse(i) */

    temp=0;
    for (i=0; i<nvar; i++)
	temp +=  u[i]*a[i];
    *sctest = temp;  /* score test */

    /*
    **  Never, never complain about convergence on the first step.  That way,
    **  if someone HAS to they can force one iter at a time.
    */
    for (i=0; i<nvar; i++) {
	newbeta[i] = beta[i] + a[i];
	}
    if (maxiter==0) {
	chinv2(imat,nvar);
	for (i=0; i<nvar; i++) {
	    beta[i] *= scale[i];  /*return to original scale */
	    u[i] /= scale[i];
	    imat[i][i] *= scale[i]*scale[i];
	    for (j=0; j<i; j++) {
		imat[j][i] *= scale[i]*scale[j];
		imat[i][j] = imat[j][i];
		}
	    }
	goto finish;
    }

    /*
    ** here is the main loop
    */
    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (*iter=1; *iter<= maxiter; (*iter)++) {
	newlk =0;
	for (i=0; i<nvar; i++) {
	    u[i] =0;
	    for (j=0; j<nvar; j++)
		imat[i][j] =0;
	    }

	/*
	** The data is sorted from smallest time to largest
	** Start at the largest time, accumulating the risk set 1 by 1
	*/
	for (person=nused-1; person>=0; ) {
	    if (strata[person] == 1) { /* rezero temps for each strata */
		denom = 0;
		nrisk =0;
		for (i=0; i<nvar; i++) {
		    a[i] = 0;
		    for (j=0; j<nvar; j++) cmat[i][j] = 0;
		    }
		}

	    dtime = time[person];
	    deadwt =0;
	    ndead =0;
	    efronwt =0;
	    while(person>=0 && time[person]==dtime) {
		nrisk++;
		zbeta = offset[person];
		for (i=0; i<nvar; i++)
		    zbeta += newbeta[i]*covar[i][person];
		risk = exp(zbeta) * weights[person];
		denom += risk;

		for (i=0; i<nvar; i++) {
		    a[i] += risk*covar[i][person];
		    for (j=0; j<=i; j++)
		    cmat[i][j] += risk*covar[i][person]*covar[j][person];
		    }

		if (status[person]==1) {
		    ndead++;
		    deadwt += weights[person];
		    newlk += weights[person] *zbeta;
		    for (i=0; i<nvar; i++) 
			u[i] += weights[person] *covar[i][person];
		    if (method==1) { /* Efron */
			efronwt += risk;
			for (i=0; i<nvar; i++) {
			    a2[i] +=  risk*covar[i][person];
			    for (j=0; j<=i; j++)
				cmat2[i][j] += risk*covar[i][person]*covar[j][person];
			    }   
		        }
	  	    }
		
		person--;
		if (strata[person]==1) break; /*tied times don't cross strata*/
	        }

	    if (ndead >0) {  /* add up terms*/
		if (method==0) { /* Breslow */
		    newlk -= deadwt* log(denom);
		    for (i=0; i<nvar; i++) {
			temp2= a[i]/ denom;  /* mean */
			u[i] -= deadwt* temp2;
			for (j=0; j<=i; j++)
			    imat[j][i] +=  (deadwt/denom)*
				(cmat[i][j] - temp2*a[j]);
		        }
    		    }
		else  { /* Efron */
		    for (k=0; k<ndead; k++) {
			temp = (double)k / ndead;
			wtave= deadwt/ ndead;
			d2= denom - temp* efronwt;
			newlk -= wtave* log(d2);
			for (i=0; i<nvar; i++) {
			    temp2 = (a[i] - temp*a2[i])/ d2;
			    u[i] -= wtave*temp2;
			    for (j=0; j<=i; j++)
				imat[j][i] +=  (wtave/d2)*
				    ((cmat[i][j] - temp*cmat2[i][j]) -
				    temp2*(a[j]-temp*a2[j]));
    		            }
    		        }

		    for (i=0; i<nvar; i++) { /*in anticipation */
			a2[i] =0;
			for (j=0; j<nvar; j++) cmat2[i][j] =0;
		        }
	            }
		}
	    }   /* end  of accumulation loop  */

	/* am I done?
	**   update the betas and test for convergence
	*/
	*flag = cholesky2(imat, nvar, toler);

	if (fabs(1-(loglik[1]/newlk))<= eps && halving==0) { /* all done */
	    loglik[1] = newlk;
	    chinv2(imat, nvar);     /* invert the information matrix */
	    for (i=0; i<nvar; i++) {
		beta[i] = newbeta[i]*scale[i];
		u[i] /= scale[i];
		imat[i][i] *= scale[i]*scale[i];
		for (j=0; j<i; j++) {
		    imat[j][i] *= scale[i]*scale[j];
		    imat[i][j] = imat[j][i];
		    }
	    }
	    goto finish;
	}

	if (*iter== maxiter) break;  /*skip the step halving calc*/

	if (newlk < loglik[1])   {    /*it is not converging ! */
		halving =1;
		for (i=0; i<nvar; i++)
		    newbeta[i] = (newbeta[i] + beta[i]) /2; /*half of old increment */
		}
	else {
	    halving=0;
	    loglik[1] = newlk;
	    chsolve2(imat,nvar,u);
	    j=0;
	    for (i=0; i<nvar; i++) {
		beta[i] = newbeta[i];
		newbeta[i] = newbeta[i] +  u[i];
		if (newbeta[i] > maxbeta[i]) newbeta[i] = maxbeta[i];
		else if (newbeta[i] < -maxbeta[i]) newbeta[i] = -maxbeta[i];
	        }
	    }
	}   /* return for another iteration */

    /*
    ** We end up here only if we ran out of iterations 
    */
    loglik[1] = newlk;
    chinv2(imat, nvar);
    for (i=0; i<nvar; i++) {
	beta[i] = newbeta[i]*scale[i];
	u[i] /= scale[i];
	imat[i][i] *= scale[i]*scale[i];
	for (j=0; j<i; j++) {
	    imat[j][i] *= scale[i]*scale[j];
	    imat[i][j] = imat[j][i];
	    }
	}
    *flag = 1000;

    
finish:

    for (person=0; person<nused; person++){
        temp=0;
        temp2=0;
        for (i=0; i<nvar; i++) {
            temp2 += beta[i]*means[i];
            temp += covart[i][person]*beta[i];
        }
        lpout[person]= temp-temp2;
        }
    
     
    double finalCPE;
    double CPE, tempCPE, bxjxi,bxixj,denomji,denomij,Scale1;
    CPE = 0;
    Scale1 = 1.0/nused;
    for(i=0; i<(nused-1); i++) {
        tempCPE = 0;
        for(j=(i+1); j<nused; j++) {
            bxjxi = lpout[j] - lpout[i];
            bxixj = 0 - bxjxi;
            denomji = 2 + expm1(bxjxi);
            denomij = 2 + expm1(bxixj);
            tempCPE += 1.0*(bxjxi <= 0)/denomji + 1.0*(bxixj < 0)/denomij;
        }
        CPE += Scale1 * tempCPE;
    }
    
   finalCPE  = 2.0*CPE/(nused - 1);
    
    unprotect(nprotect);
    return(finalCPE);
    }


/**************************************  starting file ************************************/
/**************************************  starting file ************************************/



SEXP combinedMERGEMAX(SEXP existClass, SEXP uniqueExist,
                       SEXP time2,   SEXP status2,
                       SEXP offset2, SEXP weights2,
                       SEXP strata2,   SEXP eps2, SEXP toler2){
    
    int i,j,k;
    int protected=0;
    int currentnum,currentnumpre;
    
    R_len_t nused;
    double  eps, toler;
    
    /* vector inputs */
    double *time, *weights, *offset;
    int *status, *strata;
    
    eps  = asReal(eps2);     /* convergence criteria */
    toler = asReal(toler2);  /* tolerance for cholesky */
    nused = LENGTH(offset2);
    
    time = REAL(time2);
    status = INTEGER(status2);
    weights = REAL(weights2);
    offset= REAL(offset2);
    strata = INTEGER(strata2);
    
    int *existClassPtr, *uniqueExistPtr;
    existClassPtr = INTEGER(existClass);
    uniqueExistPtr= INTEGER(uniqueExist);
    
    currentnumpre=length(uniqueExist);
    
    int numptions=0;
    numptions= 0.5*currentnumpre*(currentnumpre-1);   /* n choose 2 */
    
    int *lowgroupPtr;
    SEXP lowgroup;
    PROTECT(lowgroup = allocVector(INTSXP, numptions));
    protected++;
    lowgroupPtr=INTEGER(lowgroup);
    
    int *highgroupPtr;
    SEXP highgroup;
    PROTECT(highgroup = allocVector(INTSXP, numptions));
    protected++;
    highgroupPtr=INTEGER(highgroup);

    SEXP finaldesign;
    double *finaldesignPtr;
    PROTECT(finaldesign = allocVector(REALSXP, nused*numptions));
    protected++;
    finaldesignPtr=REAL(finaldesign);

    int ll=0;
    for(i=0; i<(currentnumpre-1); i++){
        for(j=(i+1); j < currentnumpre; j++){
            for(k=0; k <nused; k++){
                if(existClassPtr[k] < uniqueExistPtr[i]){
                    finaldesignPtr[k+(ll*nused)]=existClassPtr[k];
                }
                if((existClassPtr[k] == uniqueExistPtr[i]) ||  (existClassPtr[k] == uniqueExistPtr[j])){
                    finaldesignPtr[k+(ll*nused)]=uniqueExistPtr[i];
                }
                if((existClassPtr[k] > uniqueExistPtr[i]) &&  (existClassPtr[k] < uniqueExistPtr[j])){
                    finaldesignPtr[k+(ll*nused)]=existClassPtr[k];
                }
                if(existClassPtr[k] > uniqueExistPtr[j]){
                    finaldesignPtr[k+(ll*nused)]=existClassPtr[k]-1;
                }
                }
            lowgroupPtr[ll] =uniqueExistPtr[i];
            highgroupPtr[ll]=uniqueExistPtr[j];
            ll++;
        }
    }
    
    currentnum=currentnumpre-2;
    int totalvalid;
    totalvalid=numptions;
    
    double *currentlevels;
    SEXP currentlevels2;
    PROTECT(currentlevels2 = allocVector(REALSXP, currentnum+1));
    protected++;
    currentlevels=REAL(currentlevels2);
    
    for(i=0; i<(currentnum+1); i++){
        currentlevels[i]=i+1;
    }
    
    double **xmatrix;
    SEXP xmatrix2;
    PROTECT(xmatrix2 = finaldesign);
    protected++;
    xmatrix = dmatrix(REAL(xmatrix2), nused, totalvalid);
    
    SEXP xvecCOL2;
    PROTECT(xvecCOL2 = allocVector(REALSXP, nused));
    protected++;
    double *xvecCOL;
    xvecCOL = REAL(xvecCOL2);
    
    SEXP currentCPE2;
    PROTECT(currentCPE2 = allocVector(REALSXP, totalvalid));
    protected++;
    double *currentCPE;
    currentCPE = REAL(currentCPE2);
    
    double templist;
    for (i=0; i< totalvalid; i++) {
        templist=0;
        for(j=0; j <nused; j++){
           xvecCOL[j]=xmatrix[i][j];
        }
        templist = cpeftnM(time, status,xvecCOL, currentlevels,  nused,  currentnum,offset, weights,  strata,  eps, toler);
        currentCPE[i]=templist;
         }
    
    int maxcolumn=0;
    double maxCPE=0;
    int maxhigh=0;
    int maxlow=0;
    for (i=0; i< totalvalid; i++) {
          if( currentCPE[i] >maxCPE){
              maxCPE=currentCPE[i];
              maxcolumn= i;
              maxhigh=highgroupPtr[i];
              maxlow=lowgroupPtr[i];
          }
      }
    
    for(j=0; j <nused; j++){
        xvecCOL[j]=xmatrix[maxcolumn][j];
    }
    
    SEXP outdel2;
    PROTECT(outdel2 = allocVector(REALSXP, 3));
    protected++;
    double *outdel;
    outdel =REAL(outdel2);
    outdel[0]=maxCPE;
    outdel[1]=maxlow;
    outdel[2]=maxhigh;
    
    SEXP rlist,rlistnames;
    PROTECT(rlist= allocVector(VECSXP, 2));
    SET_VECTOR_ELT(rlist, 0, outdel2);
    SET_VECTOR_ELT(rlist, 1, xvecCOL2);

    PROTECT(rlistnames = allocVector(STRSXP, 2));
    SET_STRING_ELT(rlistnames, 0, mkChar("metrics"));
    SET_STRING_ELT(rlistnames, 1, mkChar("newriskgrp"));

    setAttrib(rlist, R_NamesSymbol, rlistnames);
    
    unprotect(protected+2);
    return(rlist);

}

