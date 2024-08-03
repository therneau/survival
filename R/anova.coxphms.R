# This is has arguments like anova.coxph: a call with a single model is
#  not allowed.  So: we avoid having separate anova.coxphms and anova.coxphmslist
#  functions.
# 
anova.coxphms <- function (object, ...,  test = c("score", "Wald", "PL")) {
    if (!inherits(object, "coxphms"))
        stop ("argument must be the fit of a multistate hazard  model")

    # All the ... args need to be coxphms fits. If any of them
    #  have a name attached, e.g., 'charlie=T' we assume a priori
    #  that they are illegal
    #
    dotargs <- list(...)
    named <- if (is.null(names(dotargs))) 
	           rep(FALSE, length(dotargs))
             else (names(dotargs) != "")
    if (any(named)) 
        warning(paste("The following arguments to anova.coxphms(..)", 
            "are invalid and dropped:", paste(deparse(dotargs[named]), 
                collapse = ", ")))

    allmod <- c(object, dotargs[!named])  # all the models
    is.coxms <- sapply(dotargs, function(x) inherits(x, "coxphms"))
    if (!all(is.coxms))
        stop("All arguments must be multistate hazard models")    

    ties <- sapply(allmod, function(x) x$method)
    if (any(ties != ties[1]))
        stop("all models must have the same ties option")
stop("anova not yet available for multistate")
    responses <- as.character(unlist(lapply(allmod, 
				     function(x) deparse(formula(x)[[2]]))))
    sameresp <- (responses == responses[1])
    if (!all(sameresp)) {
        allmod <- allmod[sameresp]
        warning(paste("Models with response", deparse(responses[!sameresp]), 
            "removed because response differs from", "model 1"))
    }
    nmodel <- length(allmod)
    if (nmodel < 2) stop("must have more than one model")

    # Check that they were all fit to the same data set.  This isn't perfect:
    #  we could have distince data sets with the same number of rows in the
    #  model frame (n), same number of unique id values, and same number of
    #  events.  But it is unlikely.  The main reason for this is that one fit
    #  had an extra variable that was missing on some subjects.
    ns <- sapply(object, function(x) c(x$n, x$n.id, x$nevent))
    if (!(is.matrix(ns))) # only occurs if a user messed with the object
        stop("at least one model is missing the n, n.id, or nevent element")
    if (any(apply(ns, 2, function(x) any(x != x[1])))) 
        stop("models were not all fit to the same dataset")

    # I can only handle models with the same stratification structure
    stest <- sapply(allmod, function(x) identical(x$map, allmod[[1]]$smap))
    if (!all(stest))
        stop("not all models have the same structure of baseline hazards")
    
    # Models must be in increasing order of complexity
    nvar <- sapply(allmod, function(x) length(coef(x)))
    if (any(diff(nvar) < 1)) 
        stop("models must be in increasing order of complexity")
    # do the more complex nesting via variable names
    for (i in 2:nmodel) {
        indx <- match(rownames(allmod[i-1][["cmap"]]),
                      rownames(allmod[i][["cmap"]]), nomatch=0)
        if (any(indx==0))
            stop(paste("model", i-1, "contains variables not in model", i))
    }
 
    # They all must have the same response
    for (i in 2:nmodel)
        if (!identical(allmod[[i]]$y, allmod[[i-1]]$y))
            stop("all models must have the same response")

    # Now for the real work
    test <- df <- pval <- double(nmodel -1)
    for (i in 2:nmodel)
    
    tfun <- function(x) paste(as.character(delete.response(terms(formula(x)))),
                              collapse=' ')
    variables <- lapply(object, tfun)
    dimnames(table) <- list(1:nmodel, 
			    c("loglik", "Chisq", "Df"))
    title <- paste("Analysis of Deviance Table\n Cox model: response is ",
		   responses[1]) 
    topnote <- paste(" Model ", format(1:nmodel), ": ", variables, 
		     sep = "", collapse = "\n")
    if (!is.null(test)) {
        table[['Pr(>|Chi|)']] <- pchisq(table$Chisq, table$Df, lower.tail=FALSE)
        }
    structure(table, heading = c(title, topnote), 
			  class = c("anova", "data.frame"))
} 

