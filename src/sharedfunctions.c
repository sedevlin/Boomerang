#include <math.h>
#include "survS.h"
#include "survproto.h"




void chinv2(double **matrix , int n)
{
    register double temp;
    register int i,j,k;
    
    /*
     ** invert the cholesky in the lower triangle
     **   take full advantage of the cholesky's diagonal of 1's
     */
    for (i=0; i<n; i++){
        if (matrix[i][i] >0) {
            matrix[i][i] = 1/matrix[i][i];   /*this line inverts D */
            for (j= (i+1); j<n; j++) {
                matrix[j][i] = -matrix[j][i];
                for (k=0; k<i; k++)     /*sweep operator */
                    matrix[j][k] += matrix[j][i]*matrix[i][k];
            }
        }
    }
    
    /*
     ** lower triangle now contains inverse of cholesky
     ** calculate F'DF (inverse of cholesky decomp process) to get inverse
     **   of original matrix
     */
    for (i=0; i<n; i++) {
        if (matrix[i][i]==0) {  /* singular row */
            for (j=0; j<i; j++) matrix[j][i]=0;
            for (j=i; j<n; j++) matrix[i][j]=0;
        }
        else {
            for (j=(i+1); j<n; j++) {
                temp = matrix[j][i]*matrix[j][j];
                if (j!=i) matrix[i][j] = temp;
                for (k=i; k<j; k++)
                    matrix[i][k] += temp*matrix[j][k];
            }
        }
    }
}




int cholesky2(double **matrix, int n, double toler)
{
    double temp;
    int  i,j,k;
    double eps, pivot;
    int rank;
    int nonneg;
    
    nonneg=1;
    eps =0;
    for (i=0; i<n; i++) {
        if (matrix[i][i] > eps)  eps = matrix[i][i];
        for (j=(i+1); j<n; j++)  matrix[j][i] = matrix[i][j];
    }
    eps *= toler;
    
    rank =0;
    for (i=0; i<n; i++) {
        pivot = matrix[i][i];
        if (pivot < eps) {
            matrix[i][i] =0;
            if (pivot < -8*eps) nonneg= -1;
        }
        else  {
            rank++;
            for (j=(i+1); j<n; j++) {
                temp = matrix[j][i]/pivot;
                matrix[j][i] = temp;
                matrix[j][j] -= temp*temp*pivot;
                for (k=(j+1); k<n; k++) matrix[k][j] -= temp*matrix[k][i];
            }
        }
    }
    return(rank * nonneg);
}



void chsolve2(double **matrix, int n, double *y)
{
    register int i,j;
    register double temp;
    
    /*
     ** solve Fb =y
     */
    for (i=0; i<n; i++) {
        temp = y[i] ;
        for (j=0; j<i; j++)
            temp -= y[j] * matrix[i][j] ;
        y[i] = temp ;
    }
    /*
     ** solve DF'z =b
     */
    for (i=(n-1); i>=0; i--) {
        if (matrix[i][i]==0)  y[i] =0;
        else {
            temp = y[i]/matrix[i][i];
            for (j= i+1; j<n; j++)
                temp -= y[j]*matrix[j][i];
            y[i] = temp;
        }
    }
}

    
    
    
    double **dmatrix(double *array, int ncol, int nrow)
    {
        
        int i;
        double **pointer;
        
        pointer = (double **) ALLOC(nrow, sizeof(double *));
        for (i=0; i<nrow; i++) {
            pointer[i] = array;
            array += ncol;
        }
        return(pointer);
    }
    
    
    
    void coxcpeOnly(const int *ROW, const double *xbeta, double * result){
        int i, j;
        double CPE, tempCPE, bxjxi,bxixj,denomji,denomij,Scale1;
        CPE = 0;
        Scale1 = 1.0/(*ROW);
        for(i=0; i<((*ROW)-1); i++) {
            tempCPE = 0;
            for(j=(i+1); j<(*ROW); j++) {
                bxjxi = xbeta[j] - xbeta[i];
                bxixj = 0 - bxjxi;
                denomji = 2 + expm1(bxjxi);
                denomij = 2 + expm1(bxixj);
                tempCPE += 1.0*(bxjxi <= 0)/denomji + 1.0*(bxixj < 0)/denomij;
            }
            CPE += Scale1 * tempCPE;
        }
        *result = 2.0*CPE/((*ROW) - 1);
    }
    
    
    void cpeOnlyNoTies(const int *ROW, const double *xbeta, double * result){
        int i,j,N;
        double CPE,tempCPE,bxjxi,bxixj,denomji,denomij,Scale1;
        
        CPE = 0;
        Scale1 = 1.0/(*ROW);
        N = 0;
        for(i=0; i<((*ROW)-1); i++) {
            tempCPE = 0;
            for(j=(i+1); j<(*ROW); j++) {
                if(xbeta[j] != xbeta[i]){
                    N = N + 1;
                    bxjxi = xbeta[j] - xbeta[i];
                    bxixj = 0 - bxjxi;
                    denomji = 2 + expm1(bxjxi);
                    denomij = 2 + expm1(bxixj);
                    tempCPE += 1.0*(bxjxi < 0)/denomji + 1.0*(bxixj < 0)/denomij;  /* [1] Kn(beta.hat) */
                }
            }
            CPE += Scale1 * tempCPE;
        }
        *result = CPE/N*(*ROW);
    }
    
    

    
