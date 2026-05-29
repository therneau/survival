# The complex NA manipuation needed for multistate coxph models that have
#  a list of formulas.  See the section on multistate missing in the
#  methods document.  
multimiss <- function(mf, y, istate, tmap, states) {
    # mf= the model frame
    # y = Surv object (also the first column of mf)
    # istate: initial state, as created or vetted by survcheck
    # tmap = terms map, from parsecovar2
    # states = the canonical list of states
    # rows that are missing y, id, istate, weights, or cluster have already
    #  been removed by the parent

    # the starting and ending state for each transition  
    # they are numbered in the order of the states argument
    from  <- as.numeric(sub(":[0-9]+$", "", colnames(tmap)))
    to    <- as.numeric(sub("^[0-9]+:", "", colnames(tmap)))

    # the states attribute of y might have only a subset of states,
    #  make a version where they line up, 0= censored
    ystat <- c(0, match(attr(y, "states"), states))[1L + y[,ncol(y)]]
    istate <- as.numeric(istate)

    # Now look at covariates.  The first row of tmap descibes baseline 
    #  hazards and can be ignored. The remainder is one row per term.
    # If someone has a term like sex:trt in the model but no main effect for trt
    #  then the rows in tmap won't properly map to columns in mf wrt missing,
    #  so we have to create tmap2 with one row for each mf column
    # Row 1 of tmap is "(Baseine)" which we can ignore
    tmap <- tmap[-1,,drop=FALSE]  # ignore this row
    tmap2 <- matrix(0, ncol(mf), ncol(tmap))  # will have 1 for "mf col was used"
    temp <- sapply(strsplit(rownames(tmap), ":"),    
                   function(x) match(x, colnames(mf)))
    for (i in seq(along.with=temp)) { 
        tmap2[temp[[i]], tmap[i,] >0] <- 1
    }
    # create a missing indicator for each column of mf
    termiss <- sapply(mf, is.na)  # will be a matrix, mf[1] is the response

    # Every term is used somewhere, otherwise mf wouldn't have included it,
    #  which means that every row of tmap has at least one non-zero value.
    # A given row will be used in all linear predictors that match its istate
    #  a. any row that has at least one successful contribution is retained; that
    #   row is part of at least one denominator in the study
    #  b. a row that contributes NA to the 1:3 transition, say, and has ystat==3
    #   causes a reduction in the transitions count
    nstate <- length(states)
    tcount <- matrix(0, nstate, nstate,
                     dimnames= list(states, states))
    # the following two matrices have one column per transition
    ismiss <- (termiss %*% tmap2)        # creates a missing value
    ispart <- outer(istate, from, "==")  # is part of this linear predictor
    isused <- rowSums(ispart & !ismiss)  # contributes a value at least one LP
    for (i in 1:length(from)) {
        tcount[from[i], to[i]] <- sum(istate== from[i] & ystat== to[i] &
                                      ismiss[,i]>0)
    }    

    # If a subject's last row(s) are censored, and all those terminal censored
    #  rows are removed per above (isused==0), then the "(censored)" column
    #  of the transitions matrix also is reduced by 1.
    # a. doing this calculation is a bit of a PITA
    # b. the rows are not necessarity sorted, which makes it worse
    # c. no one really looks at those values anyway, we are usually interested
    #   in event counts
    # This final "dotting the 'i's and crossing the 't's" part of updating the
    #  transitions table is deferred to some rainy day when I don't have 
    #  anything better to do.
    list(omit= (isused ==0), count=tcount)
}
