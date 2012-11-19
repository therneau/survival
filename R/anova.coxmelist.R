# This is usually called from anova.coxph or anova.coxme, no by a user
#  It is included as part of the survival library because the anova list
# may contain both coxph and coxme objects; if the first arg of a list is a
# coxph object then this will be called via anova.coxph
anova.coxmelist <- function (object, test =  'Chisq' ,...) {
    if (!is.list(object)) stop("First argument must be a list")
    is.coxph <- sapply(object, function(x) inherits(x, "coxph"))
    is.coxme <- sapply(object, function(x) inherits(x, "coxme"))
    if (!all(is.coxph | is.coxme))
	     stop("Argument must be a list of coxme or coxph models")
    
    # Be fussy: all the models should have the same response variable
    # and the same number of events and subjects.
    # So that the coxme library doesn't need to be loaded, we don't make
    #  use of the formula.coxme method.  (Otherwise we'd get a circular
    #  dependency that survival needs coxme which needs survival which...)
    rfun <- function(x)
        if (inherits(x, "coxph")) formula(x) else x$call$formula
            
    responses <- as.character(sapply(object, 
				     function(x) deparse(rfun(x)[[2]])))
    sameresp <- (responses == responses[1])
    if (!all(sameresp)) {
        object <- object[sameresp]
        warning(paste("Models with response", deparse(responses[!sameresp]), 
            "removed because response differs from", "model 1"))
    }
    ns <- sapply(object, function(x) if (inherits(x, "coxph")) c(x$nevent, x$n)
                                         else x$n)
    # Will be a two row matrix with deaths and total n
    if (any(ns[2,] != ns[2,1])) 
        stop("models do not have the same size of dataset")
    if (any(ns[1,] != ns[1,1]))
        stop("models do not have the same number of events")

    nmodels <- length(object)
    if (nmodels == 1) # only one model remains
        return(anova(object[[1]], test = test))

    loglik <- unlist(lapply(object, function(x) x$loglik[2]))
    dffun <- function(x) {
        if (inherits(x, "coxph")) sum(!is.na(coef(x)))
        else x$df[1]
    }
    df <- sapply(object,dffun)

    table <- data.frame(loglik, Chisq= c(NA, abs(2*diff(loglik))), 
                        Df= abs(c(NA, diff(df))))

    tfun <- function(x) paste(deparse(x$call$formula[-2]), collapse=' ')
    variables <- lapply(object, tfun)
    dimnames(table) <- list(1:nmodels, 
			    c("loglik", "Chisq", "Df"))
    title <- paste("Analysis of Deviance Table\n Cox model: response is ",
		   responses[1]) 
    topnote <- paste(" Model ", format(1:nmodels), ": ", variables, 
		     sep = "", collapse = "\n")
    if (!is.null(test)) {
        table[['P(>|Chi|)']] <- 1-pchisq(table$Chisq, table$Df)
        }
    if (is.R()) structure(table, heading = c(title, topnote), 
			  class = c("anova", "data.frame"))
    else structure(table, heading = c(title, topnote), 
			  class = "anova")	     
} 

