/* 
** Compute the Kalbfleisch-Prentice estimate of a Cox model survival curve.
** This is not the default computation, so I don't worry about making the
**  simple bisection solution below at all fast.
** See survival:Kalbfleisch-Prentice in the methods document
*/
#include "survS.h"
#include "survproto.h"

void agsurv4(int   *ndeath,   double *risk,    double *wt,
             int   *sn,        double *denom,   double *km) 
{
    int i,j,k, l;
    int n;  /* number of unique death times */
    double sumt, guess, inc;    
    
    n = *sn;
    j =0;
    for (i=0; i<n; i++) {
        if (ndeath[i] ==0) km[i] =1;
        else if (ndeath[i] ==1) { /* not a tied death */
            km[i] = pow(1- wt[j]*risk[j]/denom[i], 1/risk[j]);
            }
        else { /* biscection solution */
            guess = .5;
            inc = .25;
            for (l=0; l<35; l++) { /* bisect it to death */
                sumt =0;
                for (k=j; k<(j+ndeath[i]); k++) {
                    sumt +=  wt[k]*risk[k]/(1-pow(guess, risk[k]));
                }
            if (sumt < denom[i])  guess += inc;
            else          guess -= inc;
            inc = inc/2;
            }
            km[i] = guess;
        }
        j += ndeath[i];
    }
}
