/*
** Data sets can be created where there are useless rows: ones that are not
**  ever at risk at one of the event times.  
** It improves the accuracy of the accumulation routines if these are left
**  out.  This routine marks them.
*/
#include "survS.h"
#include "survproto.h"

int *norisk(int n, double *time1, double *time2, double *status, 
	    int *sort1, int *sort2, int *strata) {
    int i, j, istrat, p1, p2;
    int ndeath;
    double dtime;
    int *notused;

    notused = ALLOC(n, sizeof(int));

    /*
    ** Algorithm: keep a running total of deaths.
    **   When a subject enters their counter is set to  current total,
    **   when they leave, set it to 0/1 depending on whether the current
    **   total is larger than the saved count.
    ** Reset the grand total at the start of each stratum.
    ** This logic would be simpler if we walked forward in time, but the
    **  caller will have always sorted the data by (strata, reverse time, status)
    */
    ndeath=0;
    istrat=0;
    j =0;   /* tracks time1 */
    for (i=0; i<n; i++) {
	p2 = sort2[i];
	dtime = time2[p2];   
	
	if (i== strata[istrat]) {
	    /* first obs of a new stratum */
	    for (; j<i; j++) { /* finish off old stratum */
		p1 = sort1[j];
		if (ndeath > notused[p1]) notused[p1] =1; 
		else notused[p1] =0;
	    }	
	    ndeath=0;
	    istrat++;
	    }
	else {
	    for (; j<i && time1[sort1[j]] >= dtime; j++) {
		p1 = sort1[j];
		if (ndeath > notused[p1]) notused[p1] =1; 
		else notused[p1] =0;
	    }	
	} 
	ndeath += status[p2];
	notused[p1] = ndeath;
    }	

    for (; j<n; j++) {
	p1 = sort2[j];
	
	if (ndeath > notused[p1]) notused[p1] =1; 
	else notused[p1] =0;
	}

    return(notused);
}

	    
	
