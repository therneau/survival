model.frame.survfit <- function(formula, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "subset", "id", "cluster"), 
            names(dots), 0)]
    # if this was called without updating the data, and the model was saved
    #  then just return the model component
    if (length(nargs) ==0 && !is.null(formula$model)) formula$model
    else {
        fcall <- formula$call
        na.action <- getOption("na.action")
        if (is.character(na.action))
            na.action <- get(na.action)  # this is a temporary hack
        # create a copy of the call that has only the arguments we want,
        #  and use it to call model.frame()
        indx <- match(c('formula', 'data', 'weights', 'subset','na.action',
                       'istate', 'id', 'cluster', "etype"), names(fcall), 
                     nomatch=0)
        # The next error message is usually due to a typo
        #  eg survfit(wt=Surv(time, status) ~1) 
        if (indx[1]==0) stop("a formula argument is required")
        temp <- fcall[c(1, indx)]
        temp$xlev <- formula$xlevels
        if (length(nargs) > 0) 
            temp[names(nargs)] <- nargs
        
        temp[[1L]] <- quote(stats::model.frame)

        if (is.null(environment(formula$terms))) 
            mf <- eval(temp, parent.frame())
        else mf <- eval(temp, environment(formula$terms), parent.frame())

        n <- nrow(mf)
        if (n==0) stop("data set has no non-missing observations")
        mf
    }
}
