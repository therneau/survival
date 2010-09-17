# $Id: survConcordance.S 11207 2009-02-09 04:25:50Z therneau $
# Compute the concordance between a survival time x
#  and a continuous variable y
survConcordance <- function(formula, data,
                            weights, subset, na.action) {
    call <- match.call()         # save a copy of of the call, as documentation

    m <- match.call(expand=FALSE)
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    Terms <- terms(formula)
    attr(Terms,"intercept")<- 0  # There is never an intercept

    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    if (attr(Y, 'type') != 'right') stop("Only works for right-censored data")

    time <- Y[,1]
    status <- Y[,2]
    n <- length(time)

    weights <- model.extract(m, 'weights')
    if (length(weights) > 0) warning("Case weights ignored")
    offset<- attr(Terms, "offset")
    if (length(offset)>0) stop("Offset terms not allowed")

    x <- model.matrix(Terms, m)
    if (ncol(x) > 1) stop("Only one predictor variable allowed")

    ord <- order(-time, status)  # longest times first
    time <- time[ord]
    status <- status[ord]
    x <- x[ord]
    #
    #  The C code uses a balanced binary tree of the unique X values (this
    #     idea is courtesy Brad Broom).
    #   
    x2 <- sort(unique(x))
    n2 <- length(x2)
    fit <- .C("survConcordance",
              as.integer(n),
              as.double(time),
              as.integer(status),
              as.double(x),
              as.integer(n2),
              as.double(x2),
              integer(2*n2),
              result=integer(5),
              copy=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE))
    agree <- fit$result[1]
    disagree <- fit$result[2]
    tie2  <- fit$result[3]
    ties <- fit$result[4]
    incomp <- fit$result[5]

    concord <- (agree + ties/2)/(agree + disagree + ties)

    fit <- list(concordance=concord, 
                stats=c(agree=agree, disagree=disagree, 
                        tied.x=ties, tied.time=tie2, incomparable= incomp),
                call=call)
    fit$n <- n

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
    cat("Concordance= ", format(x$concordance), ", Gamma= ", 
        format((x$stats[1] - x$stats[2])/(x$stats[1] + x$stats[2])),"\n",
        sep='')
    print(x$stats)

    invisible(x)
    }
