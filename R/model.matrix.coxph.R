# In internal use "data" will often be an already derived model frame.
#  We detect this via it having a terms attribute.
model.matrix.coxph <- function(object, data=NULL, 
                               contrast.arg=object$contrasts, ...) {
    
    # If the object has an "x" component, return it, unless a new
    #   data set is given
    if (is.null(data) && !is.null(object[['x']])) 
        return(object[['x']]) #don't match "xlevels"

    Terms <- delete.response(object$terms)
    if (is.null(data)) mf <- model.frame(object)
    else {
        if (is.null(attr(data, "terms")))
            mf <- model.frame(Terms, data, xlev=object$xlevels)
        else mf <- data  #assume "data" is already a model frame
    }

    attr(Terms,"intercept")<- 1  #Cox model always has \Lambda_0
    strats <- attr(Terms, "specials")$strata
    cluster<- attr(Terms, "specials")$cluster
    dropx <- NULL
    if (length(cluster)) {
        tempc <- untangle.specials(Terms, 'cluster', 1:10)
        ord <- attr(Terms, 'order')[tempc$terms]
        if (any(ord>1)) stop ("Cluster can not be used in an interaction")
        dropx <- tempc$terms
        }
    if (length(strats)) {
        temp <- untangle.specials(Terms, 'strata', 1)
        dropx <- c(dropx, temp$terms)
        }

    # I need to keep the intercept in the model when creating the
    #   model matrix (so factors generate correct columns), then
    #   remove it.
    if (length(dropx)) {
        newTerms <- Terms[-dropx]
        # The predvars attribute, if present, is lost when we
        #  subscript.  The attribute is a Call, so has one more element
        #  than term wrt subscripting, i.e., the called function "list"
        if (!is.null(attr(terms, "predvars"))) 
            attr(newTerms, "predvars") <- attr(terms, "predvars")[-(dropx+1)]
         X <- model.matrix(newTerms, mf, contrasts=contrast.arg)
        }
    else {
        newTerms <- Terms
        X <- model.matrix(Terms, mf, contrasts=contrast.arg)
        }

    # Save attributes that are removed by subscripting, then
    #  put them back on.  Dim and dimnames are correctly changed
    #  by the subscripting.
    Xatt <- attributes(X)
    X <- X[,-1,drop=F]
    Xatt$dim <- attr(X, 'dim')
    Xatt$dimnames <- attr(X, 'dimnames')
    Xatt$assign <- Xatt$assign[-1]
    attributes(X) <- Xatt
    X
}

#  This function is somewhat confusing to read: the first argument of the
# generic model.frame is "formula", so we have to use the same label.  
# However, our first arg is actually a coxph object, from which we want
# to extract the formula!
#  There are two separate paths:
#     1. no new arguments
#          -- if there is a "model" component that that's what we want
#          -- if not, recall all the original arguments, maybe replace some
#     2. there is a new data argument
#          -- only remember the prior formula
model.frame.coxph <- function(formula, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset", "weights"), 
                        names(dots), 0)] 
    # If nothing has changed and the coxph object had a model component,
    #   simply return it.
    if (length(nargs) ==0  && !is.null(formula$model)) formula$model
    else {
        # Rebuild the original call to model.frame
        fcall <- formula$call
        indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                  names(fcall), nomatch=0) 
        if (indx[1] ==0) stop("The coxph call is missing a formula!")

        temp <- fcall[c(1,indx)]  # only keep the arguments we wanted
        temp[[1]] <- as.name('model.frame')  # change the function called
        temp$xlev <- formula$xlevels
        temp$formula <- terms(formula)  #keep the predvars attribute
        # Now, any arguments that were on this call overtake the ones that
        #  were in the original call.  
        if (length(nargs) >0)
            temp[names(nargs)] <- nargs

        # The documentation for model.frame implies that the environment arge
        #  to mf is ignored, but if we omit it there is a problem.
        if (is.null(environment(formula$terms))) 
            mf <- eval(temp, parent.frame())
        else mf <- eval(temp, environment(formula$terms), parent.frame())

       if (!is.null(attr(formula$terms, "dataClasses")))
	   .checkMFClasses(attr(formula$terms, "dataClasses"), mf)
        mf
    }         
}   

        
