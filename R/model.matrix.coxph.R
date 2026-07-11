# In internal use "data" will often be an already derived model frame.
#  We detect this via it having a terms attribute.
model.matrix.coxph <- function(object, data=NULL, 
                               contrast.arg=object$contrasts, ...) {
    # 
    # If the object has an "x" component, return it, unless a new
    #   data set is given
    if (is.null(data) && !is.null(object[['x']])) 
        return(object[['x']]) #don't match "xlevels"

    Terms <- delete.response(object$terms)
    if (is.null(data)) mf <- stats::model.frame(object)
    else {
        if (is.null(attr(data, "terms")))
            mf <- stats::model.frame(Terms, data, xlev=object$xlevels)
        else mf <- data  #assume "data" is already a model frame
    }

    cluster <- attr(Terms, "specials")$cluster
    if (length(cluster)) {
        temp <- untangle.specials(Terms, "cluster")
        dropterms <- temp$terms
    }
    else dropterms <- NULL
    
    strats <- attr(Terms, "specials")$strata
    hasinteractions <- FALSE
    if (length(strats)) {
        stemp <- untangle.specials(Terms, 'strata', 1)
        if (length(stemp$vars)==1) strata.keep <- mf[[stemp$vars]]
        else strata.keep <- strata(mf[,stemp$vars], shortlabel=TRUE)
        istrat <- as.integer(strata.keep)

        for (i in stemp$vars) {  #multiple strata terms are allowed
            # The factors attr has one row for each variable in the frame, one
            #   col for each term in the model.  Pick rows for each strata
            #   var, and find if it participates in any interactions.
            if (any(attr(Terms, 'order')[attr(Terms, "factors")[i,] >0] >1))
                hasinteractions <- TRUE  
        }
        if (!hasinteractions) dropterms <- c(dropterms, stemp$terms) 
    } else istrat <- NULL


    if (length(dropterms)) {
        Terms2 <- Terms[-dropterms]
        X <- model.matrix(Terms2, mf, constrasts.arg=contrast.arg)
        # we want to number the terms wrt the original model matrix
        temp <- attr(X, "assign")
        shift <- sort(dropterms)
        for (i in seq(along.with=shift))
            temp <- temp + 1*(shift[i] <= temp)
        attr(X, "assign") <- temp 
    }
    else X <- model.matrix(Terms, mf, contrasts.arg=contrast.arg)

    # drop the intercept after the fact, and also drop strata if necessary
    Xatt <- attributes(X) 
    if (hasinteractions) adrop <- c(0, untangle.specials(Terms, "strata")$terms)
    else adrop <- 0
    xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
    X <- X[, !xdrop, drop=FALSE]
    attr(X, "assign") <- Xatt$assign[!xdrop]
    attr(X, "contrasts") <- Xatt$contrasts
    if (length(istrat)>0) attr(X, "strata") <- strata.keep
    X
}
model.frame.coxph <- function(formula, ...) {
    Call <- match.call()
    newargs <- Call[match(c("data", "na.action", "subset", "weights",
                          "id", "cluster", "istate"), 
                        names(Call), nomatch=0)] 
    # If nothing has changed and the coxph object had a model component,
    #   simply return it.
    if (length(newargs)==0  && !is.null(formula$model)) 
        return(formula$model)
    else {
        # Rebuild the original call to model.frame. Ignore na.action
        #  and cluster: we won't use them
        newdata <- !is.null(Call$data)  # was there a "data" arg
        Terms <- terms(formula)
        fcall <- formula$call
        indx <- match(c("formula", "data", "weights", "subset", "id", "istate",
                  "cluster"), names(fcall), nomatch=0) 
        if (indx[1] ==0) stop("The coxph call is missing a formula!")
   
        temp <- fcall[c(1,indx)]  # only keep the arguments we wanted
        temp[[1]] <- quote(stats::model.frame)  # change the function called
        
        temp$na.action <- quote(stats::na.pass) # defer NA processing
        temp$xlev <- formula$xlevels  # this will turn strings to factors
        temp$formula <- Terms   #keep the predvars attribute
        # Now, any arguments that were on this call overtake the ones that
        #  were in the original call.  
        if (length(newargs) >0) 
            temp[names(newargs)] <- newargs

        # Make "tt" visible for coxph formulas, 
        if (!is.null(attr(temp$formula, "specials")$tt)) {
            coxenv <- new.env(parent= environment(temp$formula))
            assign("tt", function(x) x, envir=coxenv)
            environment(temp$formula) <- coxenv
        }

        # The documentation for model.frame implies that the environment arg
        #  to eval will be ignored, but if we omit it there is a problem.
        if (is.null(environment(formula$terms))) 
            mf <- eval(temp, parent.frame())
        else mf <- eval(temp, environment(formula$terms), parent.frame())

        if (!is.null(attr(formula$terms, "dataClasses")))
            .checkMFClasses(attr(formula$terms, "dataClasses"), mf)
       
        # None of the functions that would call model.frame (survfit, predict
        #  or residuals) allow for tt(), so we don't need to expand it.
        # But we do need to deal with timeline data
        Y <- model.response(mf)
        id <- model.extract(mf, "id")
        if (inherits(Y, "Surv2") || (!is.null(id) && any(duplicated(id)) && 
                                 attr(Y, 'type') %in% c("right", "mright"))) {
            # timeline data, convert to regular
            mf <- surv2counting(mf)
            Y <- model.response(mf)
            id <- model.extract(mf, "id")
        }   
        
        if (!newdata) {
            # remove na.action rows when recreating original data
            # for new data from a user, leave them in
            if (!is.null(formula$na.action) && length(formula$na.action)>0) {
                mf <- mf[-formula$na.action,, drop=FALSE]
                attr(mf, "na.action") <- formula$na.action
            }
        }        
        # should I also do aeqSurv here? I think yes
        timefix <- formula$timefix  #if timefix was a variable we want to eval
        if (is.null(timefix) || timefix) mf[[1]] <- aeqSurv(mf[[1]])
        mf
    }
}
