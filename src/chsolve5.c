/*
** Solve the equation Ab = y, where the cholesky decomposition of A and y
**   are the inputs.
**
** Input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**        y[n] contains the right hand side
**
**  y is overwriten with b
**
** This differs from chsolve2 only in the flag.
**   We recieved the cholesky A= LDL' where L is lower triangular, this
** is solveed in 3 stages:  L' a = y,  Db =a,  Lc = b.
** If flag=0 we do all three, if 1 we do 1 and sqrt(D)b =a,
**    if 2 we do ssqrt(D)b=a and stage 3.
** These latter support the backsolve routine.
**
**  Terry Therneau
*/
#include <math.h>

void chsolve5(double **matrix, int n, double *y, int flag) {	
    int i,j;
    double temp;

    /*
    ** solve L'z =y, 
    */
    if (flag <2) {
	for (i=0; i<n; i++) {
	    temp = y[i] ;
	    for (j=0; j<i; j++)
		temp -= y[j] * matrix[i][j] ;
	    y[i] = temp ;
	}
    }
    if (flag>0) {
	/*
	** solve D^{1/2}b =z
	*/
	for (i=0; i<n; i++) {
	    if (matrix[i][i]<=0) y[i]=0;
	    else y[i] /= sqrt(matrix[i][i]);
	}
    }
    else { /* divide by D */
	for(i=0; i<n; i++) {
	    if (matrix[i][i]==0) y[i]=0;
	    else y[i] /= matrix[i][i];
	}
    }	
    if (flag != 1) {
	/*
	** solve Lz =y
	*/
	for (i=(n-1); i>=0; i--) {
	    temp = y[i];
	    for (j= i+1; j<n; j++)
		temp -= y[j]*matrix[j][i];
	    y[i] = temp;
	}
    }
}	
