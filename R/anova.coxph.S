# $Id$
# The anova function for a coxph object
anova.coxph <- function (object, ...,  test = 'Chisq') {
    if (!inherits(object, "coxph"))
        stop ("argument must be a cox model")

    # All the ... args need to be coxph or coxme fits.  If any of them
    #  have a name attached, e.g., 'charlie=T' we assume a priori
    #  that they are illegal
    #
    dotargs <- list(...)
    named <- if (is.null(names(dotargs))) 
	           rep(FALSE, length(dotargs))
             else (names(dotargs) != "")
    if (any(named)) 
        warning(paste("The following arguments to anova.coxph(..)", 
            "are invalid and dropped:", paste(deparse(dotargs[named]), 
                collapse = ", ")))
    dotargs <- dotargs[!named]

    if (length(dotargs) >0) {
        # Check that they are all cox or coxme models
        is.coxmodel <-unlist(lapply(dotargs, function(x) inherits(x, "coxph")))
        is.coxme <- unlist(lapply(dotargs, function(x) inherits(x, "coxme")))
        if (!all(is.coxmodel | is.coxme))
            stop("All arguments must be Cox models")
        
        if (any(is.coxme))
            return(anova.coxmelist(c(list(object), dotargs), test = test))
        else return(anova.coxphlist(c(list(object), dotargs), test = test))
    }

    #
    # I have one Cox model 
    #
    if (length(object$rscore)>0)
        stop("Can't do anova tables with robust variances")

    varlist <- attr(object$terms, "variables")
    termlist<-attr(object$terms,"term.labels")
    
    specials <- attr(object$terms, "specials")

    if (!is.null(specials$strata)) {
        termlist <- termlist[-(specials$strata -1)]
        }
    
    nmodel <- length(termlist)
    df <- integer(nmodel+1)
    loglik <- double(nmodel+1)
    df[nmodel+1] <- sum(!is.na(object$coefficients))
    loglik[nmodel+1] <- object$loglik[2]
    df[1] <- 0
    loglik[1] <- object$loglik[1]

    # Fit a series of Cox models, dropping terms one by one
    # To deal properly with missings I may need a subset statement
    #   
    temp <- object$na.action
    if (!is.null(temp) && (class(temp) %in% c('exclude', 'omit')) &&
        length(temp) >0)  tsub <- -as.vector(temp)
    else tsub <- 1:object$n    

    if (is.R()) {
        fenv <- environment(formula(object))
        assign('tsub', tsub, envir=fenv)
        }

    fit <- object
    for (i in (nmodel:1)[-nmodel]) {
        form <- paste(".~ .", termlist[i], sep=' - ')
        if (is.R()) fit <-update(fit, as.formula(form,env=fenv), 
                                 subset=tsub)
        else        fit <- update(fit, as.formula(form), subset=tsub)
        df[i] <- sum(!is.na(coef(fit)))
        loglik[i] <- fit$loglik[2]
        }

    table <- data.frame(loglik=loglik, Chisq=c(NA, 2*diff(loglik)), 
                        Df=c(NA, diff(df))) 

    if (nmodel == 0) #failsafe for a NULL model
             table <- table[1, , drop = FALSE]

    if (length(test) >0 && test[1]=='Chisq') {
        table[['Pr(>|Chi|)']] <- 1- pchisq(table$Chisq, table$Df)
        }
    row.names(table) <- c('NULL', termlist)

    title <- paste("Analysis of Deviance Table\n Cox model: response is ",
		   deparse(object$terms[[2]]),
		   "\nTerms added sequentially (first to last)\n", 
		   sep = "")
    if (is.R())
	    structure(table, heading = title, class = c("anova", "data.frame"))
    else    structure(table, heading= title, class='anova')
}
