# survival
Survival package for R
This is the source code for the "survival" package in R.  It gets posted to the comprehensive R archive (CRAN) at intervals, each such posting preceded a throrough test.  In general, each new push to CRAN will update the second term of
the version number, e.g. 2.40-5 to 2.41-0.  Updates only to the github source increment after the dash.  (If an error is found in the process of CRAN submission
then the published CRAN version may be x.yy-1 or even x.yy-2.) 

The vignette "tutorial.Rnw" is not posted to CRAN, since it requires data from
the mstate package, survival is a recommended package, and such packages can 
only depend on other recommended packages.  (This allows for a consistent 
distribution bundle.)  It is included here in the vignette2 directory.  

A large portion of the source is found in the noweb directory, and is based
on the literate programming ideas of Knuth.  After making a github copy and 
before building the R package you need to do "make fun" in the noweb directory.
This will populate the R directory with a number of files whose first line is
   # Automatically generated from the noweb directory
See the noweb/Readme file for more information.

