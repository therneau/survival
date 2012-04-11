# The mf argument is mostly for internal calls, when the model frame
#  has already been constructed.  
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
#
model.frame.coxph <- function(formula, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)] 

    # If nothing has changed and the coxph object had a model component,
    #   simply return it.
    if (length(nargs) ==0  && !is.null(formula$model)) formula$model
    else {
        # First, find out what arguments existed in the original call
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

        if (is.null(environment(formula$terms)))
                eval(temp, parent.frame())
            else eval(temp, environment(formula$terms), parent.frame())
        # In the line just above the third argument is ignored since the 
        #  second arg is an environment.  But we mimic model.frame.lm by
        #  including it.
    }
}   

        
