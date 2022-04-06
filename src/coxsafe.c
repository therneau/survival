/*
** A very few pathologic cases can cause the Newton Raphson iteration
**  path in coxph to generate a horrific argument to exp().  Since all these
**  calls to exp result in (essentially) relative risks we choose a
**  fixed value of LARGE on biological grounds: any number less than 
**  1/(population of the earth) is essentially a zero, that is, an exponent
**  outside the range of +-23.  
** A sensible numeric limit would be log(.Machine$double.xmax) which is 
**  about 700, perhaps divided by 2 or log(n) to keep a few more bits.
**  However, passing this down the R calling chain to the c-routine is a lot
**  more hassle than I want to implement for this very rare case.
**
** Actually, the argument does not have to get large enough to have any
**  single exponential overflow.  In (start, stop] data we keep a running 
**  sum of scores exp(x[i]*beta), which involves both adding subjects in and
**  subtracting them out.  An outlier x value that enters and then leaves can
**  erase all the digits of accuracy.  Most machines have about 16 digits of
**  accuracy and exp(21) uses up about 9 of them, leaving enough that the
**  routine doesn't fall on it's face.  (A user data set with outlier that
**  got exp(54) and a overlarge first beta on the first iteration led to this 
**  paragraph.)  When beta-hat is infinite and x well behaved, the loglik 
**  usually converges before xbeta gets to 15, so this protection should not
**  harm the iteration path of even edge cases; only fix those that truely
**  go astray. 
**
** The truncation turns out not to be necessary for small values, since a risk
**  score of exp(-50) or exp(-1000) or 0 all give essentially the same effect.
** We only cut these off enough to avoid underflow.
*/

#define LARGE 22
#define SMALL -200

double coxsafe(double x) {
    if (x< SMALL) return(SMALL);
    if (x> LARGE) return(LARGE);
    return (x);
    }
