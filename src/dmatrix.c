/* $Id: dmatrix.c 11525 2012-12-07 17:20:39Z therneau $
**
** set up ragged arrays, with #of columns and #of rows,
**  where nrow (second arg) is what R thinks are columns
**  but C thinks are rows.
*/
#include "survS.h"
#include "survproto.h"

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
