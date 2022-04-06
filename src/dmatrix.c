/*
** set up the indices so that C code can use x[i][j] notation for R
**  matrices.  Remember that R sees matrices in column order and C in
**  row order, so every reference in the C code will be x[col][row].
**
** array = pointer to the data
** nrow, ncol = number of rows and colums, from R's point of view.
**
**  Sometime in 2015-2018 R allowed for matrices whose total number of
**  elements is > 2^31, although the numer of columns and rows must be <2^31.
**  On a 64 bit machine the array variable is of type *double
**  which is a 64 bit integer; using "array += nrow" rather than
**  "pointer[i] = array + i*nrow" is critical to avoiding an integer
**   overflow.
*/
#include "survS.h"
#include "survproto.h"

double **dmatrix(double *array, int nrow, int ncol)
    {
    int i;
    double **pointer;

    pointer = (double **) ALLOC(ncol, sizeof(double *));
    for (i=0; i<ncol; i++) {
	pointer[i] = array;
	array += nrow;
	}
    return(pointer);
    }

int **imatrix(int *array, int nrow, int ncol)
    {
    int i;
    int **pointer;

    pointer = (int **) ALLOC(ncol, sizeof(int *));
    for (i=0; i<ncol; i++) {
	pointer[i] = array;
	array += nrow;
	}
    return(pointer);
    }
