# survival
![Survival package for R](man/figures/logo.png)
 
This is the source code for the "survival" package in R.  It gets posted to the
comprehensive R archive (CRAN) at intervals, each such posting preceded a
throrough test. (I run the test suite for all 600+ packages that depend on
survival.)  In general, each new push to CRAN will update the second term of
the version number, e.g. 2.40-5 to 2.41-0.  Updates only to the github source
increment after the dash.  (If an error is found in the process of CRAN
submission then the published CRAN version may be x.yy-1 or even x.yy-2 or 3.) 

The vignette2 directory contains material that is not posted to CRAN.
The file "tutorial.Rnw", for instance, requires data from
the mstate package, survival is a recommended package, and such packages can 
only depend on other recommended packages.  (This allows for a consistent 
distribution bundle.)  The sas.Rnw vignette has a discussion of compute time and
takes too long to run, etc.

A large portion of the source is found in the noweb directory, and is based on
the literate programming ideas of Knuth. The reason is that it allows more
complete documentation of the methods. I can have things like blocks of
equations, and find having the "real" equations side by side with the code makes
it much easier to get it right.  Anyone who wants to study the methods is 
advised to perform "make code.pdf" in the noweb directory and then look at the 
relevant portion of that pdf file.  Any file in the R or src directories that
starts with an "automatically generated ..." comment should NOT be modified
directly, instead work with the noweb source.  (You will need to have the noweb
package loaded in order to run the Makefile.)

You should be able to install this using the following R code:
library(devtools); install_github("therneau/survival")

Note that good practice would be to make derived files such as R/tmerge.R
"on the fly" using a configure script; that way there would not be a 
danger of someone trying to modify the derived file rather than the actual
source (noweb/tmerge.Rnw).  However, I was not able to create a configure
file that worked reliably on all platforms, and voted for usability rather than
purity.
