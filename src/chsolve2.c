/*  $Id: chsolve2.c 11376 2009-12-14 22:53:57Z therneau $
**
** Solve the equation Ab = y, where the cholesky decomposition of A and y
**   are the inputs.
**
** Input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**        y[n] contains the right hand side
**
**  y is overwriten with b
**
**  Terry Therneau
*/
#include "survS.h"
#include "survproto.h"

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
