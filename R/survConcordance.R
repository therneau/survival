# Automatically generated from the noweb directory
# Replaced by concordance.R; this code now deprecated and will eventually be
#  removed
survConcordance <- function(formula, data,
                            weights, subset, na.action) {
    Call <- match.call()  # save a copy of of the call, as documentation
    .Deprecated("concordance")
    m <- match.call(expand.dots=FALSE)
    m[[1L]] <-  quote(stats::model.frame)
    m$formula <- if(missing(data)) terms(formula, "strata")
                 else              terms(formula, "strata", data=data)
    m <- eval(m, sys.parent())
    Terms <- attr(m, 'terms')

    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) {
        if (is.numeric(Y) && is.vector(Y))  Y <- Surv(Y)
        else stop("left hand side of the formula  must be a numeric vector or a surival")
    }
    n <- nrow(Y)

    wt <- model.extract(m, 'weights')
    offset<- attr(Terms, "offset")
    if (length(offset)>0) stop("Offset terms not allowed")

    stemp <- untangle.specials(Terms, 'strata')
    if (length(stemp$vars)) {
        if (length(stemp$vars)==1) strat <- m[[stemp$vars]]
        else strat <- strata(m[,stemp$vars], shortlabel=TRUE)
        Terms <- Terms[-stemp$terms]
    }
    else strat <- NULL
    
    x <- model.matrix(Terms, m)[,-1, drop=FALSE]  #remove the intercept
    if (ncol(x) > 1) stop("Only one predictor variable allowed")

    count <- survConcordance.fit(Y, x, strat, wt)
    if (is.null(strat)) {
        concordance <- (count[1] + count[3]/2)/sum(count[1:3])
        std.err <- count[5]/(2* sum(count[1:3]))
        }
    else {
        temp <- colSums(count)
        concordance <- (temp[1] + temp[3]/2)/ sum(temp[1:3])
        std.err <- temp[5]/(2*sum(temp[1:3]))
        }

    fit <- list(concordance= concordance, stats=count, n=n, 
                std.err=std.err, call=Call)
    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action

    oldClass(fit) <- 'survConcordance'
    fit
}

print.survConcordance <- function(x, ...) {
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
        }
    omit <- x$na.action
    if(length(omit))
        cat("  n=", x$n, " (", naprint(omit), ")\n", sep = "")
    else cat("  n=", x$n, "\n")
    cat("Concordance= ", format(x$concordance), " se= ", format(x$std.err),
        '\n', sep='')
    print(x$stats)

    invisible(x)
    }
