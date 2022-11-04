#
# This routine works in the opposite way to survSplit.
# Two data rows like this can be replaced with a single row that goes from
# 0 to 25
#   id    time1  time2   status  x1  x2
#   Jones   0      10       0     1   3
#   Jones  10      25       0     1   3
# Delete the first row and replace time1 on the second row.

survcondense <- function(formula, data, subset, id, start="tstart",
                              end="tstop", event="event") {
    Call <- match.call()
    if (missing(id)) stop("id is required")

    ss <- c("cluster")
    Terms <- if (missing(data)) terms(formula, specials=ss) else
                 terms(formula, specials=ss, data=data)
    if (length(attr(Terms, "specials")$cluster) >0)
        stop("function does not handle cluster terms")

    indx <- match(c("formula", "data", "id"), names(Call), nomatch=0)
    if (indx[1] ==0) stop("A formula argument is required")
    temp <- Call[c(1,indx)]  # only keep the arguments we wanted
    temp[[1L]] <- quote(stats::model.frame)   # change the function called
    mf <- eval(temp, parent.frame())

    Y <- model.response(mf)
    if (!is.Surv(Y)) stop("the response must be a Surv object")
    if (attr(Y, "type") != "counting")
        stop("invalid survival type")
    id <- model.extract(mf, "id")
    index <- order(id, Y[,2])  # sort by ending time
    
    nvar <- ncol(mf)
    n <- nrow(mf)
    Xdup <- sapply(2:nvar, function(i) {
        c(mf[index[-1],i] == mf[index[-n],i], FALSE)})
    Xdup <- apply(Xdup, 1, all)
    Ydup <- c(Y[index[-1], 1]== Y[index[-n], 2], FALSE)

    droprow <- unname(Xdup & Ydup)

    # There will often be clusters rows with droprow=TRUE, the start time
    #  for the first row of said cluster is moved to the first row after
    #  the cluster.
    temp1 <- rle(droprow)
    temp2 <- cumsum(temp1$lengths) # last row of each cluster of T or F

    # The last rle value is always a FALSE
    if (temp1$values[1]) { # sequence starts with a TRUE
        # j = last row of each TRUE block
        # k = first row of each TRUE block
        j <- temp2[seq(from= 1L, by=2L, to=length(temp2) -1L)] 
        k <- 1L + c(0, temp2[seq(from=2L, by=2L, to=length(temp2) -2L)])
        Y[j+1, 1] <- Y[k,1]
    } else { #sequence starts with a FALSE
        j <- temp2[seq(from=2, by=2, to= length(temp2) -1L)]
        k <- 1L + temp2[seq(from=1, by=2, to=length(temp2) -2L)]
        Y[j+1, 1] <- Y[k,1]
    }

    # thin the data rows, and remove the first colum and last= "(id)"
    newdata <- mf[!droprow, -c(1, ncol(mf)), drop=FALSE]
    row.names(newdata) <- NULL    # erase R's manufactured row names
    attr(newdata, "terms") <- NULL
    Y <- Y[!droprow,]

    # Create the response, similar to survSplit.R
    states <- attr(Y, "states")
    status <- Y[, 3]
    if (!is.null(states)) status <- factor(status, labels=c("censor", states))

    # Did the user hand me a Surv call with multiple variables, or a
    #  premade Surv object?
    if (inherits(formula[[2]], "call") && formula[[2]][[1]]== as.name("Surv")){
        # it was a call, figure out the names
        # The user might have used something like Surv(status=abc, time=fred),
        #  so use match.call to find "abc" and "fred".  But give up if there
        #  is anything complex.
        temp <- match.call(Surv, formula[[2]])
        if (missing(end) && !is.null(temp$time2) && is.name(temp$time2)) 
            end <- as.character(temp$time2)
        if (missing(event) && !is.null(temp$event) && is.name(temp$event))
            event <- as.character(temp$event)
        if (missing(start) && !is.null(temp$time) && is.name(temp$time))
            start <- as.character(temp$time)

        newdata[[start]] <- Y[,1]
        newdata[[end]]   <- Y[,2]
        newdata[[event]] <- status
    }
    else {
        if (!inherits(formula[[2]], "name"))
            stop("left hand side not recognized")
        temp <- as.character(formula[[2]])
        newdata[temp] <- Y
    }
    newdata
}

    
