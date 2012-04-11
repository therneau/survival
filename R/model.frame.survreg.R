model.frame.survreg <- function (formula, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), 
            names(dots), 0)]

    if (length(nargs) || is.null(formula$model)) {
        fcall <- formula$call
        indx <- match(c("formula", "data", "weights", "subset", 
            "na.action"), names(fcall), nomatch = 0)
        if (indx[1] == 0) 
            stop("The coxph call is missing a formula!")
        temp <- fcall[c(1, indx)]
        temp[[1]] <- as.name("model.frame")
        temp$xlev <- formula$xlevels
        if (length(nargs) > 0) 
            temp[names(nargs)] <- nargs
        if (is.null(environment(formula$terms))) 
            eval(temp, parent.frame())
        else eval(temp, environment(formula$terms), parent.frame())
     }
    else formula$model
}

# 
model.matrix.survreg <- function(object, data,  ...) {
    if (missing(data) && !is.null(object[["x"]]))
        object[["x"]]
    else {
        Terms <- delete.response(object$terms)
        strats <- attr(Terms, "specials")$strata
        cluster<- attr(Terms, "specials")$cluster
        dropx <- NULL
        if (length(cluster)) {
            tempc <- untangle.specials(Terms, 'cluster', 1:10)
            dropx <- tempc$terms
        }
        
        if (length(strats)) {
            temp <- untangle.specials(Terms, 'strata', 1)
            dropx <- c(dropx, temp$terms)
        }       

        if (length(dropx)) {
            newTerms <- Terms[-dropx]
            # R (version 2.7.1) adds intercept=T anytime you drop something
            attr(newTerms, 'intercept') <- attr(Terms, 'intercept')
            # The predvars attribute, if present, is lost when we
            #  subscript.  The attribute is a Call, so has one more element
            #  than term wrt subscripting, i.e., the called function "list"
            if (!is.null(attr(terms, "predvars"))) 
                attr(newTerms, "predvars") <- attr(terms, "predvars")[-(dropx+1)]
        }
        else newTerms <- Terms


        # Grab the model frame.  By using "newterms" for a new data set,
        #  we allow the new data to be missing things we don't need: y,
        #  strata, and cluster.  For the original data we can assume they
        #  are present.
        if (missing(data)) 
            mf <- model.frame(object, ...)
        else {
            if (is.null(attr(data, "terms")))
                mf <- model.frame(newTerms, data, xlev=object$xlevels)
            else mf <- data  #assume we were given a model frame     
        }
        model.matrix(newTerms, mf, contrasts.arg= object$contrasts)
    }
}
