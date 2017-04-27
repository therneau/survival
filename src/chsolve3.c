/* $Id: chsolve3.c 11166 2008-11-24 22:10:34Z therneau $ */
/*
** Solve the equation Ab = y, where the cholesky decomposition of A and y
**   are the inputs.
** This is a specialized form for the frailty problem.  The matrix C in this
**   case has C[1:m, 1:m] diagonal and  C[(m+1):n, 1:n)] is dense. 
**
** arguments are:
**     n         the size of the matrix to be factored
**     m         the size of the diagonal upper portion
**     diag      the diagonal upper portion
**     **matrix, which contains the chol decomp of the dense portion
**     y[n] contains the right hand side
**
**  y is overwriten with b
**
**  Terry Therneau
*/
#include "survS.h"
#include "survproto.h"

void chsolve3(double **matrix, int n, int m, double *diag, double *y)
     {
     int i,j, n2;
     double temp;

     n2 = n-m;
     /*
     ** solve Fb =y   (the diagonal portion is unchanged)
     */
     for (i=0; i<n2; i++) {
	  temp = y[i+m];
	  for (j=0; j<m; j++) temp -= y[j]   * matrix[i][j];
	  for (j=0; j<i; j++) temp -= y[j+m] * matrix[i][j+m] ;
	  y[i+m] = temp ;
	  }
     /*
     ** solve DF'z =b
     */
     /* dense portion */
     for (i=(n2-1); i>=0; i--) {
	  if (matrix[i][i+m]==0)  y[i+m] =0;
	  else {
	      temp = y[i+m]/matrix[i][i+m];
	      for (j= i+1; j<n2; j++)
		   temp -= y[j+m]*matrix[j][i+m];
	      y[i+m] = temp;
	      }
	  }
     /* diag portion */
     for (i=(m-1); i>=0; i--) {
	 if (diag[i] == 0)  y[i] =0;
	 else {
	     temp = y[i] / diag[i];
	     for (j=0; j<n2; j++)
		  temp -= y[j+m]*matrix[j][i];
	     y[i] = temp;
	     }
         }
     }


