# There are two major differences with model.frame.coxph
# Less complex: we don't have to worry about tt() terms, since coxphms models 
#  don't allow them.  Most of the work is recreating the model.frame call

model.frame.coxphms <- function(formula, data, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "subset", "weights",
                          "id", "cluster", "istate"), 
                        names(dots), 0)] 
    # If nothing has changed and the coxph object had a model component,
    #   simply return it.
    if (length(nargs) ==0  && missing(data) && !is.null(formula$model)) 
        return(formula$model)
    
    # Otherwise, rebuild the original call to model.frame
    Terms <- terms(formula)
    fcall <- formula$call
    indx <- match(c("formula", "data", "weights", "subset", "na.action",
                    "cluster", "id", "istate"),
              names(fcall), nomatch=0) 
    if (indx[1] ==0) stop("The coxph call is missing a formula!")

    temp <- fcall[c(1,indx)]  # only keep the arguments we wanted
    temp[[1]] <- quote(stats::model.frame)  # change the function called
    temp$xlev <- formula$xlevels  # this will turn strings to factors
    temp$formula <- Terms   #keep the predvars attribute
    if (is.list(fcall$formula))
        temp$na.action <- na.pass # defer missings until later

    # Now, any arguments that were on this call overtake the ones that
    #  were in the original call.  
    if (length(nargs) >0)
        temp[names(nargs)] <- nargs

    # The documentation for model.frame implies that the environment arg
    #  to eval will be ignored, but if we omit it there is a problem.
    if (is.null(environment(formula$terms))) 
        mf <- eval(temp, parent.frame())
    else mf <- eval(temp, environment(formula$terms), parent.frame())

    if (!is.null(attr(formula$terms, "dataClasses")))
        .checkMFClasses(attr(formula$terms, "dataClasses"), mf)

    if (inherits(model.response(mf), "Surv2")) 
        mf <- surv2counting(mf) # this is Surv2 style data, modifiy it

    # Here is a subtle bit.  If the coxph call did not save the
    #  model frame (and so we had to rebuild it above), it will
    #  nevertheless have saved the na.action attribute based on
    #  rows to be deleted, using the multimiss logic/routine below.
    # Thus we don't have to repeat that logic, just apply it. However,
    #  that attribute applies to the saved data frame, which is in
    #  counting process form, i.e., do this after the surv2data transform
    # Exception!: if the data argument is present, then this is a user call
    #  to build a new model.frame from a new data set, using the same formula
    #  and etc.  The most common reason for this is to get new predicted
    #  values, and this routine will like have been called by model.matrix.
    # In this case we don't want to do anything at all about missing values:
    #  if the user put in a data set of length n he/she should get back something
    #  of that length, and if an input is NA they should get back an NA
    #
    if (!is.null(formula$na.action)) mf[-formula$na.action,]
    else mf
}

# The complex NA manipuation needed for multistate coxph models that have
#  a list of formulas.  See the section on multistate missing in the
#  methods document
multimiss <- function(mf, y, id, istate, tmap) {
    # first off: any row missing y, id, or weight will be lost
    miss0 <- is.na(y) | is.na(id) # y and id are not optional
    if (!is.null(istate)) miss0 <- miss0 | is.na(istate)
    if (!is.null(istate)) miss0 <- miss0 | is.na(istate)
    weights <- model.weights(mf)
    if (!is.null(weights)) miss0 <- miss0 | is.na(weights)
    cluster <- model.extract(mf, "cluster")
    if (!is.null(cluster)) miss0 <- miss0 | is.na(cluster)

    # the starting state for each transition
    from  <- as.numeric(sub(":[0-9]+$", "", colnames(tmap)))

    # Now look at covariates.  The first row of tmap descibes baseline 
    #  hazards and can be ignored. The remainder is one row per term.
    # If someone has a term like sex:trt in the model but no main effect for trt
    #  then the rows in tmap won't properly map to columns in mf wrt missing,
    #  so we have to create tmap2 with one row for mf column
    # Row 1 of tmap is "(Baseine)" which we can ignore
    tmap2 <- matrix(0, ncol(mf), ncol(tmap))  # will have 1 for "mf col was used"
    temp <- sapply(strsplit(rownames(tmap)[-1], ":"),    
                   function(x) match(x, colnames(mf)))
    for (i in seq(along.with=temp)[-1]) {  # skip the (Baseline) row
        tmap2[temp[[i]], tmap[i,] >0] <- 1
    }

    # create a missing indicator for each column of mf
    termiss <- sapply(mf, is.na)  # will be a matrix, mf[1] is the response
    # matrix that shows if row i of the data contibutes a missing to transition j
    makemiss <- termiss %*% tmap2 # positive only if missing and used
    omit <- miss0 | (rowSums(makemiss) >0)
    omit
}
