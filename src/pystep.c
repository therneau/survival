/*  $Id: pystep.c 11166 2008-11-24 22:10:34Z therneau $ */
/*
** Returns the amount of time that will be spent in the current "cell",
**  along with the index of the cell (treating a multi-way array as linear).
** This is a basic calculation in all of the person-years work.
**
** Input
**      nc:         number of categories
**      data[nc]    start points, for the data values
**      fac[nc]     1: category is a factor, 0: it is continuous
**                  >=2: special handling for "years" dim of US rate tables
**      dims[nc]    the extent of each category
**      cuts[nc,dims+1] ragged array, containing the start for each interval
**      step        the amount of time remaining for the subject.
**      edge        if =0, then the cuts contain +1 obs, and we are strict
**                      about out-of-range cells.  If it is a 1, then the
**                      table is assummed to extend infinitly at the edges.
**
** Output
**      *index   linear index into the array
**      if *index == -1, then the returned amount of time is "off table";
**  if one of the dimensions has fac >1 --
**      *index2   second index for linear interpolation
**      *wt       a number between 0 and 1, amount of wt for the first index
**                this will be 1 if none of the dims have fac >1
**
** Return value     amount of time in indexed cell.
*/
#include "survS.h"
#include "survproto.h"

double pystep(int nc,        int  *index,  int  *index2,   double *wt, 
	      double *data,  Sint *fac,    Sint *dims,     double **cuts, 
	      double step,   int  edge)
    {
    int i,j;
    double maxtime;
    double shortfall;
    double temp;
    int kk, dtemp;


    kk=1;
    *index =0;  *index2=0;
    *wt =1;
    shortfall =0;
    maxtime = step;
    for (i=0; i<nc; i++) {
	if (fac[i]==1) *index += (data[i]-1) * kk;
	else {
	    if (fac[i]>1) dtemp = 1 + (fac[i]-1)*dims[i];
	    else          dtemp = dims[i];
	    for (j=0; j<dtemp; j++) if (data[i] < cuts[i][j]) break;

	    if (j==0) {  /* less than first cut */
		temp = cuts[i][j] - data[i];  /* time to next cutpoint */
		if (edge==0 && temp > shortfall) {
		    if (temp > step) shortfall = step;
		    else             shortfall = temp;
		    }
		if (temp < maxtime)  maxtime = temp;
		}
	    else if (j==dtemp){  /*bigger than last cutpoint */
		if (edge==0) {
		    temp = cuts[i][j] - data[i];  /* time to upper limit */
		    if (temp <=0) shortfall = step;
		    else if (temp < maxtime) maxtime = temp;
		    }
		if (fac[i] >1) j = dims[i] -1;   /*back to normal indices */
		else  j--;
		}
	    else {
		temp = cuts[i][j] - data[i];  /* time to next cutpoint */
		if (temp < maxtime)  maxtime = temp;
		j--;
		if (fac[i] >1) { /*interpolate the year index */
		    *wt = 1.0 - (j%fac[i])/ (double)fac[i];
		    j /= fac[i];
		    *index2 = kk;
		    }
		}
	    *index += j*kk;
	    }
	kk *= dims[i];
	}

    *index2 += *index;
    if (shortfall ==0) return(maxtime);
    else  {
	*index = -1;
	return(shortfall);
	}
    }
