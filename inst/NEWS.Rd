\name{NEWS}
\title{NEWS file for the survival package}
\section{Changes in version 2.37-7}{
  \itemize{
    \item Remove a dependency on the survey package from the adjusted
    survival curves vignette, at the request of CRAN.  (The base +
    required bundle needs to be capable of a stand-alone build.)
    
    \item Fix error in calcuation of the y-axis range for survival curve
    plots whenever the "fun" argument could produce infinite values,
    e.g., complimentary log-log plots transform 1 to -Inf.  Pointed out
    by Eva Boj del Val.  (Add finite=TRUE to range() call).
}}

\section{Changes in version 2.37-6}{
  \itemize{
    \item The plot for competing risk curves could have a spurious
    segment.  (Found within 3 hours of submitting 2.37-5 to CRAN.)

    \item The lines method for survexp objects was defaulting to a
    step function, restore the documented default of a connected
    line.
    \item Add a levels method for tcut objects. 14Jan2014
  }
}  
\section{Changes in version 2.37-5}{
  \itemize{
    \item Add vignette on adjusted survival curves.
    \item Add vignette concerning "type 3" tests.
    
    \item Make the tt() function invisible outside of a coxph formula.
    There was a complaint about conflicts with another package, and
    there is not really a good reason to have it be a global name.
    An R-devel discussion just over 1 year ago showed how to accomplish this.
   
    \item The modeling routines are set in two parts, e.g., coxph sets
    up the model and coxph.fit does the work.  Export more of the ".fit"
    routines to make it easier for other packages to build on top of
    this one.
    
    \item Updates to the model.matrix and model.frame logic for
    coxph.  A note from F Harrell showed that I was not correctly
    dealing with the "assign" attribute when there are strata * factor
    interactions.  This led to cleanup in other cases that I had missed
    but which never had proven fatal.  Also added support for tt() terms
    to the stand alone model.matrix and model.frame functions.
    (Residuals for tt models are still not available, but this was a
    necessary first step to that end.) 26Dec13
    
    \item The Surv function now remembers attributes of the input
    variables that were passed to it; they are saved as
    "inputAttributes".  This allows the rms package, for instance, to
    retain labels and units through the call.
    
    \item Update summary.coxph.penal to produce an object, which in turn
    has a print method, i.e., make it a "standard" summary function.
    
    \item Add a logLik method for coxph and survfit objects.
    
    \item Allow for Inf as the end of the time interval, for interval
    censored data in the Surv function.
    
    \item The predict.coxph function would fail if it had both a newdata
    and a collapse argument.  Pointed out by Julian Bothe.  25Sep13
    
    \item Survexp can now produce expecteds based on a stratified Cox
    model.  Add the 'individual.s' and 'individual.h' options to return
    indivudual survival and cumulative hazard estimates, respectively.
    The result of survfit now (sometimes) includes the cumulative
    hazard.  This will be expanded. 29Jul13
    
    \item Change code in the coxpenal.fit routine: the use of a vector of
    symbols as arguments to my .C calls was confusing to a new CRAN
    consistency check.  Both the old and new are legal R; but the old
    was admittedly an unusual construction and it was simpler to change it.
    
    \item Fix a bug in survfit.coxph pointed out by Chris Andrews, whose
    root cause was incorrect curve labels when the id option is used. 27Jun13
    
    \item Add rsurvreg routine.
    
    \item Change survfit.coxph routine so that it detects whether
    newdata contains or does not contain strata variables, and acts
    accordingly. If newdata does containe strata then  the output will
    contain only those data-value and strata combinations specified by the
    user. Retain strata levels in the coxph routine for use
    in the survfit routine, to correctly reconstruct strata levels.
    Warn about curves with interactions.  18Ju13

    \item Add a dim method for survival curves.
    
    \item For competing risks curves that use the istate option, the
    plotted curves now start with the correct (initial) prevalence of
    each state.  22May13
    
    \item The survreg function failed with the "robust=T"
    option. Pointed out by Jon Peck.  Test case added.  6May13
    
    \item Kazuki Yoshida pointed out that rep() had no method for Surv
    objects.  This caused the survSplit routine to fail if the data
    frame contained a Surv object.  3May13
    
    \item Per a request from Milan Bouchet-Valet fix an issue in survfit
    that arose when the OutDec option is set to ',': it did not
    correctly convert times back from character to numeric.

    \item The plot.survfit function now obeys "cex" for the size of the
      marks used for censored observations.
}}

\section{Changes in version 2.37-4}{
  \itemize{
    \item Subscripting error in predict.coxph for type=expected, se=T,
    strata in the model, newdata, and multiple strata in the new data
    set.  Pointed out by Chris Andrews.  The test program has been
    tweaked to include multiple strata in newdata.
}}
    

\section{Changes in version 2.37-3}{
  \itemize{
    \item Minor flaw in [.survfit.  If "fit" had multiple curves, and
    fit$surv was a matrix, and one of those curves had only a single observation
    time, fit[i,] would collapse columns when "i" selected that curve,
    though it shouldn't.  

    \item Changed all of the .C and .Call statements to make use of
    "registered native routines", per R-core request. Add file src/init.c 

    \item Error in plot.survfit pointed out by K Hoggart -- the "+"
    signs for censored observations were printing one survival time to
    the left of the proper spot.  Eik Vettorazi found another error if
    mark.time is a vector of numerics.  These are the results of merging
    the code for plot, lines and points due to some discrepancies between
    them, plus not having any graphical checks in the test suite.
    
    \item Repair an error in using double subscripts for the survfitms
    objects.
    
    \item Add the US population data set, with yearly totals by age and sex 
    for 2000 onward.  It is named uspop2, since there is already a "uspop" 
    data set containing decennial totals from 1790 to 1970.
    
    \item Not all combinations of strata Y/N and CI Y/N worked in the
    quantile.survfit function, pointed out by Daniel Wallschlaeger
    (missing a function argument in one if-else combination).
    Added a new test routine that verifies all paths.

    \item The first example in predict.survreg help file needed to have
    \code{I(age^2)} instead of \code{age^2} in the model: R ignores the
    second form. (I'm almost sure this worked at one time, perhaps in Splus).
    It also needed different plot symbols to actually match the
    referenced figure.  Pointed out by Evan Newell.

    \item Fix a long-standing problem with cch pointed out by Ornulf Borgan
    leading to incorrect standard errors.  A check in the underlying
    coxph routines to deal with out of bounds exponents, added in
    version 2.36-6, interacted badly with the -100 offset used in cch.
    It only affected models using (start, stop) survival times. 
}}

\section{Changes in version 2.37-2}{
  \itemize{
    \item Two bugs were turned up by running tests for all the
    packages that depend on survival (158 of them).
}}
\section{Changes in version 2.37-1}{
  \itemize{
    \item Add a new multi-state type to the Surv object.  Update the
    survfit routine to work with it.  The major change is addition
    of a proper variance for this case.
    More functionality is planned.
    
    \item Remove the fr_colon.R test program.  It tests an ability that
    has been superseded by coxme, on a numerically touchy data set, and
    it was slow besides.  For several other tests that produce warning
    messages and are supposed to produce said messages, add extra
    comments to that effect so testers will know it is expected.
    
    \item The code has had several "if.R" clauses to accomodate Splus vs
    R differences, which are mostly class vs oldClass.  These are now being
    removed as I encounter them; since our institution no longer uses
    Splus I can no longer test the clauses' validity.
    
    \item The fast subsets routine coxexact.fit incorrectly returned the
    linear predictor vector in the (internal) sorted order rather than
    data set order.  Pointed out by Tatsuki Koyama, affecting the result
    of a clogit call.  6Nov2012
    
    \item Jason Law pointed out that the sample data set "rats" is from
    the paper by Mantel et.al, but the documentation was for a data
    set from Gail, Santner and Brown.  Added the Gail data as rats2 and
    fixed the documentation for rats.
    
    \item For predict.coxph with type="terms", use "sample" as the
    default value for the reference option.  For all others the default
    remains "strata", the current value.  Type terms are nearly always passed
    forward for further manipulation and per strata centering can mess
    things up: termplot() for instance will no longer show a smooth
    function if the results are recentered within strata.
        
    \item Fix bug in summary.aareg, which was unhappy (without cause) if
    the maxtime option was used for a fit that did not include the
    dfbeta option.  Pointed out by Asa Johannesen.
    
    \item The coxph fitting functions would report an error for a null model
      (no X variables) if init was specified as numeric(0) rather than NULL.
    
    \item Update the description and citation files to use the new
    "person" function described in the R Journal.  Also add the
    ByteCompile directive per suggestion of R core.

    \item Allow an ordinary vector as the left hand side of survConcordance.

    \item Update anova.coxphlist to reject models with a robust variance.
 
    \item The survfit function had an undocumented backwards-compatability 
    that allows the newdata argument to be a vector with no names.  An
    example from Damon Krstajic showed that this does not work when the
    original model has a matrix in the formula.  Removed the feature.
    (This is for survfit.coxph.)  Also clarified the code and its 
    documentation about what is found where -- environments, formulas,
    and the arguments of eval, which fixes a problem pointed out by xxx
    where the result of a Surv call is used in the coxph formula.

    \item Fix an issue in summary.survfit pointed out by Frank Harrell.  The
    strata variable for the output always had its labels in sorted order,
    even when a factor creating the survival curves was otherwise.  (This was
    due to a call to factor() in the code.)  The print routine would then
    list curves in sorted order, which might well be contrary to the user's
    wishes.  The curves were numerically correct.

    \item Add the anova.coxmelist function to the namespace so that it is
    visible.  If someone has a list of models the first of which was
    a coxph fit and the list includes coxme fits, then anova.coxph will
    be the function called by R, and it will call anova.coxmelist.

    \item Fix a bug pointed out by Yi Zhang and Mickael Hartweg.
    If a coxph model used an offset, then a predicted survival curve
    that used newdata (and the offset variable of course) would be
    wrong, e.g. survival values > 1.
    A simple misplaced parenthesis was the cause.
    A recent paper by Langholz shows how to get absolute survival from
    case-control data using an offset, which seems to have suddenly made this
    feature popular.

    \item Per further interaction with Yi Zhang, a few items were
    missing from the S3methods in the NAMESPACE file: as.matrix.Surv,
    model.matrix.coxph, model.matrix.survreg, model.frame.survreg.
}}

\section{Changes in version 2.36-14}{
  \itemize{
    \item A supposedly cosmetic change to coxph in the last release
	caused formulas with a "." on the right hand side to fail.  Fix this
	and add a case with "." to the test suite.
      }
    }

\section{Changes in version 2.36-13}{
  \itemize{
    \item Add the anova.coxmelist function.  This is in the survival
    package rather than in coxme since "anova(fit1, fit2)" is valid when
    fit1 is a coxph and fit2 a coxme object, a case which will cause this
    function to be called by way of anova.coxph.
    
    \item More work on "predvars" handling for the pspline function,
    when used in predict calls.  Add a new test of this to the suite,
    and the makepredictcall method to the namespace.
    Fixes a bug pointed out by C Crowson.  

    \item Deprecate the "robust" option of coxph.  When there are
    multiple observations per subject it is almost surely the wrong
    thing to do, while adding a "cluster(id)" term does the correct
    thing. When there is only one obs per subject both methods work
    correctly.

    \item Add documentation of the output structure to the aareg help
    file.

    \item Change ratetableDate so that it still allows use of chron
    objects, but doesn't need the chron library.  This eliminates a
    warning messge from the package checks, but is also a reasonable
    support strategy for a moribund package.  (Some of the local users keep
    datasets for a long long time.)

    \item Fix a bug in summary.survfit for a multiple-strata survival
    object.  If one of the curves had no data after application of the
    times argument, an output label was the wrong length.
    
    \item Fix a bug pointed out by Charles Berry: predict for a Cox
    model which has strata, and the strata is a factor with not all its
    levels represented in the data.  I had a mistake in the subscripting
    logic: number of groups is not equal to max(as.integer(strata)).
    
    \item Changes to avoid overflow in the exponent made in 2.36-6
    caused failure for one special usage: in case-cohort designs a dummy
    offset of -100 could be added to some observations.  This was being
    rounded away.  The solution is to 1: have coxsafe not truncate small
    exponents and 2: do not recenter user provided offset values.

    \item Fix bug in survfit.coxph.  Due to an indexing error I would
    sometimes create a huge scratch vector midway through the
    calculations (size = max value of "id"); the final result was always
    correct however.  Data set provided by Cindy Crowson which had a
    user id in the billions.

    \item Fix bug pointed out by Nicholas Horton: predictions of
    type expected, with newdata, from a Cox model without a strata
    statement would fail with "x not found".  A misplaced parenthesis
    from an earlier update caused it to not recreate the X matrix even
    though it was needed later.  Also add some further information to
    the predict manual page to clarify an issue with frailty terms.
    }}

\section{Changes in version 2.36-12}{
  \itemize{
    \item Fix a bug in the new fast subsets code.  The test suite had no
    examples of strata + lots of tied times, so of course that's the
    case where I had an indexing error.  Add a test case using the
    clogit function, which exercises this.

    \item Further memory tuning for survexp.
    }}

    \section{Changes in version 2.36-11}{
  \itemize{
    \item Make survexp more efficient.  The X matrix was being modified
    in several places, leading to multiple copies of the data.  When the
    data set was large this would lead to a memory shortage.

    \item Cause anova.coxph to call anova.coxme when a list of models
    has both coxph and coxme objects.

    \item Add the quantile.survfit function.  This allows a user to
    extract arbitrary quantiles from a fitted curve (and std err).

    \item Fix an error in predict.coxph.  When the model had a strata
    and the newdata and reference="sample" arguments were used, it
    would (incorrectly) ask for a strata variable in the new data set.

    \item Incorporate the fast subsets algorithm of Gail et al, when 
    using coxph with the "exact" option.  The speed increase is profound 
    though at the cost of some memory.    Reflect this in the
    documentation for the clogit routine.
    Note that the fast computation is not yet implemented for
    (start,stop) coxph models.

    \item Change the C routine used by coxph.fit from .C to .Call
    semantics to improve memory efficiency, in particular fewer copies
    of the X matrix.

    \item Add scaling to the above routine.  This was prompted by a user
    who had some variables with a 0-1 range and others that were
    0 - 10^7, resulting in 0 digits of accuracy in the variance matrix.
    (Economics data).

    \item Comment out some code sections that are specific to Splus.
    This reduced the number of "function not found" warnings from R CMD check.
}}

\section{Changes in version 2.36-10}{
  \itemize{
    \item 30 Sept 2011: The na.action argument was being ignored in
    predict.coxph; pointed out by Cindy Crowson.

    \item The log-likelihood for survreg was incorrect when there are
    case weights in the model.  The error is a fixed constant for any
    given data set, so had no impact on tests or inferences.  The error
    and correction were pointed out by Robert Kusher.
    
    \item A variable name was incorrect in survpenal.fit.  This was in a
    program path that had never been traversed until Carina Salt used
    survreg with a psline(..., method='aic') call, leading to a
    "variable not found" message.

    \item Punctuation error in psline made it impossible for a user to
    specify the boundary.knots argument.  Pointed out by Brandon
    Stewart.

    \item Add an "id" variable to the output of survobrien.

    \item The survfitCI routine would fail for a curve with only one
    jump point (a matrix collapsed into a vector).

    \item Fix an error in survfit.coxph when the coxph model has both a
    strata by covariate interaction and a cluster statement.  The
    cluster term was not dropped from the Terms object as it should have
    been, led to a spurious "variable not found" error.   Pointed out
    by Eva Bouguen.

    \item If a coxph model with penalized terms (frailty, pspline) also
    had a redundant covariate, the linear predictor would be returned as
    NA.  Pointed out by Pavel Krivitsky.
}}

\section{Changes in version 2.36-9}{
  \itemize{
    \item Due to a mistake in my script that submits to CRAN, the fix in
    2.36-8 below was actually not propogated to the CRAN submission.

    \item Fix an error in the Cauchy example found in the survreg.distributions
    help page, pointed out by James Price.

    \item Update the coxph.getdata routine to use the model.frame.coxph
    and model.matrix.coxph methods.

    \item Add the concordance statistic to the printout for penalized models.
}}

\section{Changes in version 2.36-8}{
  \itemize{
    \item Unitialized variable in calcuation of the variance of the
    concordance.  Found on platform cross-checking by Brian Ripley.
    \item Changed testci to use a fixed file of results from cmprsk
    rather than invoking that package on-the-fly.  Suggested by the CRAN
    maintainers.
  }
  }
\section{Changes in version 2.36-7}{
  \itemize{
    \item Due to changes in R 2.13 default printout, the results of many of
    the test programs change in trivial way (one more or fewer digits).
    Update the necessary test/___.Rout.save files.  Per the core team's
    suggestion the dependency for the package is marked as >=2.13.
  }}
  
\section{Changes in version 2.36-6}{
  \itemize{
    \item An example from A Drummond caused iteration failure in
    coxph: x=c(1,1,1,0,1, rep(0,35)), time=1:40, status=1.  The first
    iteration overshoots the solution and lands on an almost perfectly
    linear part of the loglik surface, which made the second iteration
    go to a huge number and exp() overflows.  A sanity check routine
    coxsafe is now invoked on all values of the linear predictor.
    
    \item 1 April: Fix minor bug in survfit.  For left censored data
    where all the left censored are on the very left, it would give a
    spurious warning message when trying to create a 0 row matrix that
    it didn't need or use.  Pointed out by Steve Su.
    
    \item 31 March 2011: One of the plots in the r_sas test was wrong
    (it's been a long time since I visually checked these).  The error
    was in predict.survreg; it had not taken into account a change in
    R2.7.1: the intercept attribute is reset to 1 whenever one
    subscripts a terms object, leading to incorrect results for a model
    with "-1" in the formula and a strata(): the intercept returned
    when removing the strata.  I used this opportunity to move most of
    the logic into model.frame.survreg and model.matrix.survreg
    functions.  Small change to the model.frame.coxph and
    model.matrix.coxph functions due to a better understanding of
    xlevels processing.
    
    \item Round off error issue in survfit: it used both unique(time)
    and table(time), and the resulting number of unique values is not
    guarranteed to be the same for times that differ by a tiny amount.
    Now times are coverted to a factor first.  Peter Savicky from the R
    core team provided a nice discussion of the issue and helped me
    clarify how best to deal with it.  The prior fix of first rounding
    to 15 digits was good enough for almost every data set -- except the
    one found by a local user just last week.
    
    \item Round off error in print.survfit pointed out by Micheal Faye.
    If a survival value was .5 in truth, but .5- eps due to round off
    the printed median was wrong. But it was ok for .5+eps. Simple
    if-then logic error.
    
    \item Re-fix a bug in survfit.  It uses both unique and table in
    various places, which do not round the same; I had added a
    pre-rounding step to the code.  A data set from Fan Chun showed that
    I didn't round quite enough. But the prior rounding did work for a
    time of 2 vs (sqrt(2))^2: this bug is very hard to produce.  I now
    use as.numeric(as.character(factor(x))), which induces exactly the
    same rounding as table, since it is the same compuation path.

    \item Further changes to pspline.  The new Boundary.knots argument
    allows a user to set the boundary knots inside the range of data.
    Code for extrapolation outside that range was needed,
    essentially a copy of the code found in ns() for the same issue.
    Also added a psplineinverse function, which may be useful with
    certain tt() calls in coxph.
    
    \item 10 Mar 2011: Add the capablilty for time-dependent
    transformations to coxph, along with a small vignette describing use
    of the feature.  This code is still incompletely incorporated in
    that the models work but other methods (residuals, predict, etc) are
    not yet defined.
    
    \item 8 Mar 2011: Expand the survConcordance function.  The function
    now correctly handles strata and  time dependent covariates, and
    computes a standard error for the estimate.  All computation is based
    on a balanced binary tree strucure, which leads to computation in
    \eqn{O(n \log_2(n))}{O(n log(n))} time.
    The \code{coxph} function now adds concordance to its output, and
    \code{summary.coxph} displays the result.

    \item 8 Mar 2011: Add the "reference" option to predict.coxph, a
    feature and need pointed out by Stephen Bond.

    \item 4 Mar 2011: Add a makepredictcall method for pspline(), which
    in turn required addition of a Boundary.knots argument to the
    function.
    
    \item 25 Feb 2011: Bug in pyears pointed out by Norm Phillips.  If a
    subject started out with "off table" time, their age was not
    incremented by that amount as they moved forward to the next "in
    table" cell of the result.  This could lead to using the wrong
    expected rate from the rate table.
  }
}

\section{Changes in version 2.36-5}{
  \itemize{
    \item 20 Feb 2011: Update survConcordance to correctly handle case
    weights, time dependent covariates, and strata.

    \item 18 Feb 2011: Bug in predict.coxph found by a user (1 day after
    36-4!).  If the coxph call had a subset and predict used newdata,
    the subset clause was "remembered" in the newdata construction,
    which is not appropriate.
}}

\section{Changes in version 2.36-4}{
\itemize{
  \item 17 Feb 2011: Fix to predict.coxph.  A small typo that only was
  exercised if the coxph model had x=T.  Discovered via induced error in
  the rankhazard package.  Added lines to the test suite to test for
  this in the future.

  \item Removed some files from test and src that are no longer needed.

  \item Update the configure script per suggestion from Kurt H.
  }}
\section{Changes in version 2.36-3}{
  \itemize{

    \item 13 Feb 2011: Add the rmap argument to pyears, as was done for
    survexp, and update the manual pages and examples. Fix one last bug
    in predict.coxph (na.action use).  Passes all the
    tests for inclusion on the next R release.

    \item 8 Feb 2011: Change the name of the new survfit.coxph.fit routine to
     survfitcoxph.fit; R was mistaking it for a survfit method.  Fix
     errors in predict.coxph when there is a newdata argument, including
     adding yet another test program.

   \item 1 Feb 2011: Fix bugs in coxph and survreg pointed out by Heinz
    Tuechler and dtdenes@cogpsyphy.hu, independently, that were the same
    wrong line in both programs.  With interactions, a non-penalized
    term could be marked as penalized due to a mismatched vector length,
    leading to a spurious error message later in the code.
  
  \item 1 Feb 2011: Update survfit.coxph to handle the case of a strata by
    covariate interaction.  All prior releases of the code did this
    wrong, but it is a very rare case (found by Frank Harrell).  Added a
    new test routine coxsurv4.  Also found a bug in [.survfit; for a
    curve with both strata and multiple columns, as produced by
    survfit.coxph, it could drop the n.censored item when subscripting.
    A minor issue was fixed in coxph: when iter=0 the output coefficient
    vector should be equal to the input even when the variance is
    singular. 

  \item 30 Jan 2011: Move the noweb files to a top level directory, out of
    inst/.  They don't need to be copied to binary installs.
    
  \item 22 Jan 2011: Convert the Changelog files to the new inst/NEWS.Rd
    format.
  
  \item 1 Jan 202011: The match.ratetable would fail when passed a data
    frame with a character variable.  This was pointed out by Heinz
    Tuechler, who also did most of the legwork to find it.  It was
    triggered by the first few lines of tests/jasa.R (expect <- ....)
    when options(stringsAsFactors=FALSE) is set. 
 }
}

\section{Changes in version 2.36-2}{
\itemize{
  \item 20 Dec 2010: Add more test cases for survfit.coxph,
    which led to significant updates in the code. 

  \item 18 Nov 2010: Add nevent to the coxph output and printout in
    response to a long standing user request.

  \item 14 Dec 2010: Add an as.matrix method for Surv objects. 

  \item 11 Nov 2010: The prior changes broke 5 packages: the dependencies form a bigger
test suite than mine!  1. Survival curve for a coxph model with sparse
frailty fit; fixed and added a new test case.  2. survexp could fail if 
called from within a function due to a scoping error.  3. "Tsiatis" was once
a valid type (alias for 'aalen')  for survfit.coxph; now removed from the 
documentation but the code needed to be backwards compatable.  The other two
conflicts were fixed in the packages that call survival.  There are still
issues with the rms package which I am working out with Frank H.
 }
}

\section{Changes in version 2.36-1}{
 \itemize{
  \item{27 Oct 2010: Finish corrections and test to the new code. It now passes
the checks.  The predict.coxph routine now does strata and standard errors
correctly, factors propogate through to predictions, and numerous small 
errors are addressed.  Predicted survival curves for a Cox model has been
rewritten in noweb and expanded. Change the version number to 2.36-1.}

  \item{17 Oct 2010: Per a request from Frank Harrell (interaction with his library),
survfit.coxph no longer reconstructs the model frame unless it really needs
it: in some cases the 'x' and 'y' matrices may be sufficient, and may be
saved in the result.  Add an argument "mf" to model.matrix.coxph for more 
efficient interaction when a parent routine has already recovered the model 
frame.
  In general, we are trying to make use of model.matrix.coxph in many of the
routines, so that the logic contained there (remove cluster() calls, pull
out strata, how to handle intercepts) need not be replicated in multiple 
places.}

  \item{12 Oct 2010: Fix a bug in the modified lower limits for survfit (Dory & Korn).
A logical vector was being inadvertently converted to numeric.  Pointed out
by Andy Mugglin.  A new case was added to the test suite. }
 }
}

\section{Changes in version 2.35}{ 
 \itemize{
  \item{15 July 2010: Add a coxph method for the logLik function.  This is used by
the AIC function and was requested by a user.}

  \item{29 July 2010: Fix 2 bugs in pyears.  The check for a US rate table was off
(minor effect on calculations), and there was a call to julian which assumed
that the origin argument could be a vector.  }

  \item{21 July 2010: Fix a problem pointed out by a user: calling survfit with almost
tied times, e.g., c(2, sqrt(2)^2), could lead to an inconsistent result.  Some
parts of the code saw these as 2 unique values per the unique() function, some
as a single value using the results of table().  We now pre-round the input
times to one less decimal digit than the max from .Machine$double.digits.
  Also added the noweb.R processing function from the coxme package, so that
the noweb code can be extracted "on the fly" during installation using
commands in the configure and cleanup scripts.  }

  \item{11 July 2010: A rewrite of the majority of the survfit.coxph code.  The primary
benefits are 1: finally tracked down and eliminated the bug for standard errors
of case weights + Cox survival + Efron method; 2: the individual=TRUE and FALSE
options now use the same underlying code for curves, before there were some
options valid only for one or the other; 3: code was rewritten using noweb 
with a considerable increase in documentation; 4: during the verification
process some errors were found in the test suite and corrected, e.g., a 
typo in my book led to failure of an all.equal test in book4.R.  Similar
to the rewrite for survfit several years ago, the new code has far less use
of .C to help transparency.}

  \item{21 May 2010: Fix bug in summary.survfit.  For a survival curve from a Cox model
with start,stop data, the 'times' argument would generate an error.}

  \item{24 May 2010: Fix an annoyance in summary.survfit.  When the survival data 
had an event or censor at time 0 and summary is called with a times argument,
then my constructed call to approx() would have duplicate x values.  The answer
was always right, but approx has begun to print a bothersome warning message.  
A small change to the constructed argument vector avoids it.}

  \item{7 April 2010: Minor bug pointed out by Fredrik Lundgren.  In survfit if the
method was KM (default) and error = Tsiatis an error message results.  Simple
fix: code went down the wrong branch.}

  \item{24 Feb 2010: Serious bug pointed out by Kevin Buhr.  In Surv(time1, time2,stat)
if there were i) missing values in time1 and/or time2, ii) illegal value
sets with time1 >=time2, and iii) all the instances of ii do not preceed all
the instances of i, then the wrong observation (not the illegal) will be 
thrown out.  Repaired, and a new test added.
  Minor updates to 3 test files: survreg2, testci, ratetable.}

  \item{8 Feb 2010: Bug pointed out by Heinz Tuechler -- if a subscript was dropped from
a rate table the 'type' attribute got dropped, e.g. survexp.usr[,1,,].}

  \item{26 Jan 2010: At the request of Alex Bokov, added the xmax, xscale, and fun
arguments to points.survfit.}

  \item{26 Jan 2010: Fix bug pointed out by Thomas Lumley -- with case weights <1 a Cox 
model with (start, stop) input would inappropriately decide it needed to do
step halving to find a solution, eventually failing to converge.  It was
treating a loglik >0 as an indication of failure, but such values arise for
small case weights. Let L(w) be the loglik for a data set where everyone is
given a weight of w, then L(w)= wL(1) - d log(w) where d=number of deaths in
the data.  For small enough w positivity of L(w) is certain.}

  \item{25 Jan 2010: Fix bug in summary.ratetable pointed out by Heinze Tuechler.  Added
a call to the function to the test suite as well.}

  \item{15 Dec 2009: Two users pointed out a bug that crept into survreg() with a
cluster statement, when a t(x)%*%x was replaced with crossprod.  A trivial
fix, but in response I added another test that more formally checks the
dfbeta residuals and found a major oversight for the case of multiple
strata.  }

  \item{14 Dec 2009: 1.Fix bug in frailty.xxx, if there is a missing value in the levels
it gets counted by "length(unique(x))" (frailty is called before NA
removal.)  2.SurvfitCI had an incorrect CI with case weights, and 3. in
survreg a call to resid instead of residuals.survreg, before the class
was attached.}

  \item{11 Nov 2009: The 'type' argument does not make sense for plot.survfit.  (If 
type='p', should one plot the tops of the step function, the bottoms, or
both?).  Make it explicitly disallowed in response to an R-help query,
rather than the confusing error message that currently arose.}

  \item{28 Oct 2009: The basehaz function would reorder the labels of the strata
factor.  Not a bug really, but a "why do this?"  Unintended consequence of
a character -> factor conversion.}

  \item{1 Oct 2009: Fix a bug pointed out by Ben Domingue.  There was one if-then-else
path into step-halving in the frailty.controldf routine that would refer to
a non-existent variable.  A very rarely followed path, obviously, and with
the obvious fix. The mathematics of the update was fine.}

  \item{30 Sep 2009: For coxph and model.matrix.coxph, re-attach the
attributues lost from the X matrix when the intercept is removed,
i.e., X <- X[,1].  In particular, some downstream libraries depend on
the assign attribute. 
  For predict.coxph remove an earlier edit so that a single variable model
+ type='terms' returns a matrix, not a vector.  This is expected by the 
termplot() function.  It led to a whole lot of changes in the test suite
results, though, due to more "matrix" printouts.}

  \item{4 Sep 2009: Added a model.matrix.coxph and model.frame.coxph methods.  The
model.matrix.default function ceased to work for coxph models sometime
between R 2.9 and 2.9.2 (best guess).  This wasn't picked up in the test
suite but rather by failure of 3 packages that depend on survival.  Also added
a test.  Update CRAN since this broke other's packages.}

  \item{20 Aug 2009: One more fix to predict.coxph.  It needed to use 
delete.response(Terms) rather than Terms, so as to not look for (unnecessarily)
the response variable when the newdata argment is used.  Pointed out by Michael
Conklin.}

  \item{17 Aug 2009: Small bug in survfit.coxph.null pointed out by Frank Harrell.  The
'n' component would be missing if the input data included strata, i.e., the
initial model had used x=TRUE.  He also pointed out the fix.}

  \item{10 June 2009: Fix an error pointed out by Nick Reich, who was the first to use
interval censored data + user defined distribution in survreg, jointly.  There
was no test case and creating one uncovered several errors (but only for this
combination).  All the error cases led to catastrophic failure, highlighting
the extreme rarity of a user requesting this combination.}

  \item{2 June 2009: Surv(time1, time2, status, type='interval') would fail for an NA
status code.  Pointed out by Achim Zeilus.}

  \item{22 May 2009: Allow single subscripts to rate tables, e.g. survexp[1:10: .  
Returns a simple vector of values.  The str() function does this to print out
a short summary.  Problem pointed out by Heinz Tuechler.}

  \item{21 May 2009: Create a test case for factor variables/newdata/predict for coxph 
and survreg.  This led to a set of minor fixes; the code is now in line with
the R standard for model functions.  One consequence is that model.frame.coxph
and model.frame.survreg are no longer needed, so have been removed.}

  \item{20 May 2009: The manual page for survfit was confusing, since it tries to 
document both the standard KM (formula method) and the coxph method.  I've
split them out so that now survfit documents only the basic method and points 
a user the appropriate specialized page.}

  \item{1 May 2009: The anova.coxph function was incorrect for models with a strata
term.  Fixed this, and made chisquare tests the default.}

  \item{22 April 2009: The coxph code had an override to iter and eps, making both of
them more strict for a penalized model.  However, the overall default values
have changed over time, so that these lines actually decreased accuracy - the
opposite of their intent.  Removed the lines.  Also removed the iter.miss and
eps.miss components (on which this check depended) from coxph.control, which
makes that function match its documentation.}

 }
}

\section{Changes in version 2.34 and earlier}{
\subsection{Merge of the TMT source code tree with the Lumley code tree}{
 \itemize{

  \item Issues/decisions in remerging the Mayo and R code: 
  For most of routines, it was easier to start with the Lumley code and add
the Therneau fixes.  This is because Tom had expanded a lot of partial 
matches, e.g., fit$coef in the TT code vs fit$coefficients.  Routines with
substantial changes were, of course, a special case.
  The most common change is an is.R() construct to choose class vs oldClass.
  \item xtras.R: Move anova.coxph and anova.coxphlist to their own
  source files.
  The remainder of the code is R only.

  \item survsum: removed from package

  \item survreg.old: has been removed from the package

  \item survfit.s: 
   Depreciate the "formula with no ~1" option
   Mayo code for [ allows for reordering curves
   Separate out the R "basehaz" function as a separate source file

  \item survfit.km.s: The major change of did not get copied into R, so lots of
changes.  R had "new.time" and Splus 'start.time' for the same argument.  Allow
them both as synonyms.
   The output structure also changed: adapt the new one.  This is mostly some
name changes in the components, removing unneeded redundancies created by
a different programmer. 

  \item survfit.coxph.s:  TMT code finally fixed the "Can't (yet) to case weights" 
problem.  There must have been 10 years been the intent and execution. 

  \item survexp.s:  Add "bareterms" function from R, which replaces a prior use of
  terms.inner (in Splus but not R). 

  \item survdiff.s: R code had the old (incorrect) expected <- sum(1-offset), since
corrected to sum(-log(offset)) . 

  \item{summary.coxph.s: This was a mess, since Tom and I had independently made the
 addition of a print.summary.coxph function.  Below, TMT means that it was the
 choice in the Splus code, TL means that it was the choice in R
	1. Put the coef=T argument in the print function, not summary (TMT)
	2. Change the output's name from coef to coefficients (suggestion of
  Peter Dalgaard).  Also change one column name to Pr(>|z|) for R.
	3. Remove last vestiges of a reference to the 'icc' component (TMT)
	4. Do not include score, rscore, naive.var in the result (TL)
	5. Do include loglik in the result (TMT)
	6. Compute the test statistics (loglik, Wald, etc) in the summary
   function rather than in the print.summary function (TL)
	7. Remove the digits option from summary, it belongs in print.summary.
   (neither)}

  \item{strata.s: R code added a sep argument, this is ok
     R changed the character string NA to as.character(NA).  Not okay
          1. won't work with Splus, 
          2. This is a label, designed for printing, 
	          and so it should be a character string.  }

  \item{residuals.coxph.s: R had added type='partial'.  (Which I'm not very partial to,
  from their statistical properties.  But they are legal, and I assume that
  someone requested them).}

  \item{print.survfit.s: Rewritten as a part of the general survival rewrite.  Created
the function 'survmean' which does most of the work, and is shared by print and
summary, so that the values from 'print' are now available.
   Fix the minmin function: min(NULL) gives NA in Splus, which is the right
answer for a non-estimable median, but Inf in R.  Explicitly deal with this
case, and add a bunch of comments.
   R had the print.rmean option, this has been expanded to a more general
rmean option that allows setting the cutoff point.
   R added a print.n option with 3 choices, my code includes all 3 in the 
output.  }

  \item{lines.survfit.s: 
    The S version has a new block of code for guessing "firstx" more
intellegently when it is missing.  (Or, one hopes is is more intellegent!)}

  \item{coxph.control.s: 
   The R code had tighter tolerances (eps= 1e-9) than Splus (1e-4) and
a higher iterationn count (20 vs 10). 
   Set eps to 1e-8 and iter to 15, mostly bending to the world.  The tighter
iteration is defensible, but I still maintain that a Cox model that takes >10 
iterations is not going to finish if you give it 100.  The likelihood surface
is almost perfectly quadratic near the minimum.  (Not true for survreg by the
way).}

  \item{: In Surv, the Mayo code creates NA's out of invalid status values or
start,stop pairs, rather than a stop and error message.  This is to
allow for example
       coxph(Surv(time1,time2, status).... , subset=(goodlines))
succeed, when "goodlines" is the subset with correct values.}
}
}

\subsection{Older changes}{
 \itemize{
  \item{25SepO7: How embarrassing -- someone pointed out that I had Dave Harrington's
name spelled wrong in the options to survfit.coxph!}

  \item{9Jul07: In a model with offsets, survreg mistakenly omitted the offset
from the returned linear.predictor component.}

  \item{10May07: Change summary.coxph so that it returns an object of class
summary.coxph, and add a print method for that object.}

  \item{22Jun06: Update match.ratetable, so that more liberal matches are now
allowed.  For instance,  'F', 'f', 'female', 'fem', 'FEMA', etc are 
now all considered matches to the dimname "female" in survexp.us.}

  \item{26Apr06: Fix bug in summary.survfit, pointed out by Bob Treder.  With
the times option, the value of n.risk would be wrong for "in between"
times; e.g., the data had events and/or censoring at times 10,
20,... and we asked for printout at time 15.  It should give n.risk at time
20, it was returning the value at time 10. Interestingly, the code had
a very careful treatment of this case, along with an example in the
comments, and the "the right answer is" part of the comment was wrong!
So the code correctly computed an incorrect answer.  Added another
test case to the test suite, survtest2.}

  \item{21Apr06: Fix problem in [.survfit, pointed out by Thomas Lumley.  If
fit <- survfit(Surv(time, status) ~ ph.ecog, lung), then fit[2:1] did
not reorder the output correctly.  I had never tested putting the 
subscripts in non-increasing order.}

  \item{7Feb06: Fix a problem in the coxph iteration (coxfit2.c, coxfit5, agfit3,
agfit5, agexact).  It will likely never catch anyone again, even if I 
didn't fix it.  In a particular data set, beta overshot and step halving
was invoked.  During step halving, a loglik happened to occur that was
within eps of the prior step's loglik --- and the routine decided, erroneously,
that it had converged!  (A nice quadratic curve, a first guess b1 to the 
left of  the desired max of the curve.  The next guess b2 overshot and ends 
up with a lower loglik, on the right side of the max.  Back up to the 
midpoint of b1 and b2, and this guess, still to the right of the max (still
too large) has EXACTLY the same value of y as b1 did, but on the other side
of the max from b1.  "Last two guesses give the same answer, I'm done" said
the routine).}

  \item{27Sep05: Found and fixed a nasty bug in survfit.  When method='fh2' and
there were multiple groups I had a subscripting bug, leading to vectors
that were supposed to be the same length, but weren't, passed into C.
The resulting curves were obviously wrong -- survival precipitously drops 
to zero.}

  \item{5May05: Add the drop=F arg to one subscripting selection in survfit.coxph.
	        temp <- (matrix(surv$y, ncol=3))[ntime,,drop=F] 
If you selected only 1 time point (1 row) in the final output, the code
would fail.  Pointed out by Cindy Crowson.}

  \item{18Apr05: Bug in survfit.turnbull.  The strata variable was not being
filled in (number of points per curve).  So if multiple curves were
generated at once, i.e., with something on the right hand side of ~ in
the formula, all the downstream print/plot functions would not work
with the result.}

  \item{8Feb05: Fix small typo in is.ratetable, introduced on 24Nov04: (Today
was the first time I added to the standard library, and thus ended
up using the non-verbose mode.)}

  \item{8Feb05: Add the data.frame argument to pyears.  This causes the output
to contain a dataframe rather than a set of arrays.  It is useful for
further processing of the data using Poisson regression.}

  \item{7Feb05: Modified print.ratetable to be more useful.  It now tells 
about the ratetable, rather than printing all of its values.}

  \item{8Dec04: Fix a small bug in survfit.turnbull. If there are people left 
censored before the first
time point of any other kind (interval, exact, or right censored),
the the plotted height of the curve from "rightmost left censoring time"
to "leftmost event time", that is the flat tail on the left, was at
the wrong height.
  Added another test to testreg/reliability.s for this.}

  \item{24Nov04: Change is.ratetable to give longer messages}
}
}
}
