#
# Routine to turn a Surv2 type dataset into a Surv type of data set
#  The task is pretty simple
#     1. An id with 10 rows will have 9 in the new data set, 1-9 contribute
#        covariates, and 2-10 the endpoints for those rows.
#     2. Missing covariates are filled in using last-value-carried forward.
#     3. A response, id, current state, and new data set are returned.
# If check=FALSE, it is being called by survcheck.  In that case don't fail
#   if there is a duplicate time, but rather let it worry about it.
#

surv2data <- function(mf, check=FALSE) {
    Terms <- terms(mf)
    y <- model.response(mf)
    if (!inherits(y, "Surv2")) stop("response must be a Surv2 object")
    n <- nrow(y)
    states <- attr(y, "states")

    id <- model.extract(mf, "id")
    if (length(id) != n) stop("id statement is required")

    # relax this later
    if (any(is.na(id)) || any(is.na(y[,1])))
        stop("id and time cannot be missing")
    
    isort <- order(id, y[,1])
    id2 <- id[isort]
    y2  <- y[isort,]
    first <- !duplicated(id2)
    last  <- !duplicated(id2, fromLast=TRUE)

    status <- y2[!first,2]
    y3 <- cbind(y2[!last,1], y2[!first,1], ifelse(is.na(status), 0, status))
    if (!is.null(states)) {
        if (all(y3[,1] ==0)) {
            y3 <- y3[,2:3]           
            attr(y3, "type") <- "mright"
        }
        else attr(y3, "type") <- "mcounting"
        attr(y3, "states") <- states
    }
    else {
        if (all(y3[,1] ==0)) {
            y3 <- y3[,2:3]
            attr(y3, "type") <- "right"
        }
        else attr(y3, "type") <- "counting"
    }
    class(y3) <- "Surv"

    if (!check && ncol(y3)==3 && any(y3[,1]==y3[,2])) 
        stop("duplicated time values for a single id")

    # We no longer need the last row of the data
    # tmerge3 expects the data in time within id order
    mf2 <- mf[isort,][!last,]  
    id3 <- id2[!last]
    # Use LVCF on all the variables in the data set, other than
    #  the response and id.  The response is always first and
    #  the id and cluster will be at the end
    fixup <- !(names(mf) %in% c("(id)", "(cluster)"))
    if (is.factor(id3)) idi <- as.integer(id3)
    else    idi <- match(id3, unique(id3))  
    for (i in (which(fixup))[-1]) {
        miss <- is.na(mf2[[i]])
        if (any(miss)) {
            k <- .Call(Ctmerge3, idi, miss)
            update <- (miss & k>0)  # some values will be updated
            if (any(update))  { # some successful replacements
                if (is.matrix(mf2[[i]])) mf2[[i]][update,] <-mf2[[i]][k[update],]
                else mf2[[i]][update] <- (mf2[[i]])[k[update]]
            }
        }}      
    
    # create the current state vector
    #  this is simply the non-lagged y -- except, we need to lag the state
    #  over missing values.  Do that with survcheck.
    istate <- y2[!last, 2]
    first <- !duplicated(id3)
    n <- nrow(y3)
    temp <- istate[first]
    if (all(is.na(temp) | (temp==0))) {
        # special case -- no initial state for anyone
        temp <- survcheck(y3 ~1, id=id3)
    } 
    else if (any(is.na(temp) | (temp==0)))
        stop("everyone or no one should have an initial state")
    else {
        # survcheck does not like missing istate values.  Only the initial
        #  one for each subject needs to be correct though
        itemp <- istate[match(id3, id3)]
        if (is.null(states)) 
            temp <- survcheck(y3~1, id=id3, istate= itemp)
        else temp <- survcheck(y3~1, id=id3, 
                               istate=factor(itemp, 1:length(states), states))
    }
     
    if (check) list(y=y3, id=id3, istate= temp$istate, mf= mf2, isort=isort,
                    last=last)
    else list(y=y3, id=id3, istate= temp$istate, mf= mf2)
}
   
# User callable version
Surv2data <- function(formula, data, subset, id){
    Call <- match.call()
    indx <- match(c("formula", "data", "weights", "subset", "na.action",
                    "cluster", "id", "istate"),
                  names(Call), nomatch=0)
    if (indx[1] ==0) stop("A formula argument is required")
    tform <- Call[c(1, indx)]  # only keep arguments we wanted
    tform$na.action <- stats::na.pass
    tform[[1L]] <- quote(stats::model.frame)
    mf <- eval(tform, parent.frame())
    temp <- surv2data(mf, check=FALSE)
    
    mf2 <- temp$mf
    mf2[['S2.y']] <- temp$y
    index <- match(as.character(Call$id), names(mf2), nomatch=0)
    if (index >0) mf2[[index]] <- temp$id
    else    mf2[['S2.id']] <- temp$id
    mf2[['S2.state']] <- temp$istate
    attr(mf2, "terms") <- NULL
    mf2
}

    
