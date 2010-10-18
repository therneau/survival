# The mf argument is mostly for internal calls, when the model frame
#  has already been constructed.  
model.matrix.coxph <- function(object, data=NULL,contrast.arg=object$contrasts,
                               mf, ...){
    if (is.null(data) && missing(mf) && !is.null(object[['x']])) 
        object[['x']] #don't match "xlevels"
    else {
        Terms <- object$terms
        if (missing(mf)) {
            if (is.null(data)) mf <- model.frame(object, ...)
            else mf <- model.frame(object, data=data, ...)
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

        if (length(dropx)) {
            # I need to keep the intercept in the model when creating the
            #   model matrix (so factors generate correct columns), then
            #   remove it.
            newTerms <- Terms[-dropx]
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
    }

model.frame.coxph <- function(formula, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset", "drop.unused.levels",
                          "xlev"), names(dots), 0)]

    if (length(nargs) ==0  && !is.null(formula$model)) formula$model
    else {
        fcall <- formula$call
        indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                  names(fcall), nomatch=0) 
        if (indx[1] ==0) stop("A formula argument is required")

        temp <- fcall[c(1,indx)]  # only keep the arguments we wanted
        temp[[1]] <- as.name('model.frame')  # change the function called
        temp$xlev <- formula$xlevels

        if (length(nargs) >0)
            temp[names(nargs)] <- nargs
        if (is.R()) {
            if (is.null(environment(formula$terms)))
                eval(temp, parent.frame())
            else eval(temp, environment(formula$terms))
            }
        else  eval(temp, sys.parent())
        }
    }

        
