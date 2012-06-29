# $Id$
# This is usually called from anova.coxph, not a user
#  It's first argument must be a list of coxph models
anova.coxphlist <- function (object, test =  'Chisq' ,...) {
    if (!is.list(object)) stop("First argument must be a list")
    if (!all(unlist(lapply(object, function(x) inherits(x, 'coxph')))))
	     stop("Argument must be a list of coxph models")
    if (any(sapply(object, function(x) !is.null(x$naive.var))))
        stop("Can't do anova tables with robust variances")
    
    responses <- as.character(unlist(lapply(object, 
				     function(x) deparse(formula(x)[[2]]))))
    sameresp <- (responses == responses[1])
    if (!all(sameresp)) {
        object <- object[sameresp]
        warning(paste("Models with response", deparse(responses[!sameresp]), 
            "removed because response differs from", "model 1"))
    }
    ns <- sapply(object, function(x) length(x$residuals))
    if (any(ns != ns[1])) 
        stop("models were not all fitted to the same size of dataset")
    nmodels <- length(object)
    if (nmodels == 1) # only one model remains
        return(anova.coxph(object[[1]], test = test))
    loglik <- unlist(lapply(object, function(x) x$loglik[2]))
    df <- unlist(lapply(object, function(x) sum(!is.na(coef(x)))))

    table <- data.frame(loglik, Chisq= c(NA, abs(2*diff(loglik))), 
                        Df= abs(c(NA, diff(df))))

    tfun <- function(x) paste(as.character(delete.response(terms(formula(x)))),
                              collapse=' ')
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

