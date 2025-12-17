# There are two major differences with model.frame.coxph
# Less complex: we don't have to worry about tt() terms, since coxphms models 
#  don't allow them
# More complex: missing value omission has to happen as a second step. If 
#  subject i is missing variable zed, but they participate in a transtion that
#  doesn't use zed, that observation needs to stay.  Only obs that will not
#  be used in any transition are tossed, e.g., the variable 'sex' is used in
#  all transitions then any row with missing sex will be omitted.  
# The stacker routine will eventually deal with this and create a stacked
#  data set that has no missing rows.
#
model.frame.coxphms <- function(formula, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset", "weights",
                          "id", "cluster", "istate"), 
                        names(dots), 0)] 
    # If nothing has changed and the coxph object had a model component,
    #   simply return it.
    if (length(nargs) ==0  && !is.null(formula$model)) return(formula$model)
    
    # Otherwise, rebuild the original call to model.frame
    Terms <- terms(formula)
    fcall <- formula$call
    indx <- match(c("formula", "data", "weights", "subset", "na.action",
                    "cluster", "id", "istate"),
              names(fcall), nomatch=0) 
    if (indx[1] ==0) stop("The coxph call is missing a formula!")

    temp <- fcall[c(1,indx)]  # only keep the arguments we wanted
    temp[[1]] <- quote(stats::model.frame)  # change the function called
    temp$xlev <- formula$xlevels  # this will turn strings to factors
    temp$formula <- Terms   #keep the predvars attribute
    temp$na.action <- na.pass # defer missings until later

    # Now, any arguments that were on this call overtake the ones that
    #  were in the original call.  
    if (length(nargs) >0)
        temp[names(nargs)] <- nargs

    # The documentation for model.frame implies that the environment arg
    #  to eval will be ignored, but if we omit it there is a problem.
    if (is.null(environment(formula$terms))) 
        mf <- eval(temp, parent.frame())
    else mf <- eval(temp, environment(formula$terms), parent.frame())

    if (!is.null(attr(formula$terms, "dataClasses")))
        .checkMFClasses(attr(formula$terms, "dataClasses"), mf)

    cmap  <- formula$cmap  # coefficient map, one row per coef, col per state
    # the starting state for each transition
    from  <- as.numeric(sub(":[0-9]+$", "", colnames(cmap)))
    #  In a multi-state not every variable might be used in every transition,
    # only data rows that have no use at all are removed. So for instance if sex
    # is missing in an obs, but that particular obs is at risk for at least one 
    # transition that doesn't depend on sex, said obs is retained in the data
    # frame. (The stacker routine will deal with partial missings, later.)
    #  The model frame has terms but our cmap matrix has covariates, e.g., mf
    # has a factor and cmap the dummy variables created from it.  In theory we
    # could use the terms structure to map each covariate to its term, but it
    # is easier to go the other way and use model.matrix to map from terms
    # to covariates.
    #  First though, any obs with missing Y, state, weight, id, or cluster is
    # out regardless
    Y <- model.response(mf)
    id <- model.extract(mf, "id")
    istate <- model.extract(mf, "istate")
    miss0 <- is.na(Y) | is.na(id) # Y and id are not optional
    if (!is.null(istate)) miss0 <- miss0 | is.na(istate)
    mcheck <-survcheck2(Y, id, istate) # use survcheck to curate istate
    istate <- mcheck$istate

    weights <- model.weights(mf)
    if (!is.null(weights)) miss0 <- miss0 | is.na(weights)
    cluster <- model.extract(mf, "cluster")
    if (!is.null(cluster)) miss0 <- miss0 | is.na(cluster)

    # now look at covariates
    # Not every obs contributes to every transition, delete those rows that
    #  will will be tossed for all transitions.  Variables with names like
    #  "ph(2:1)" in cmap will later be constructed, and are not in X matrix
    # 
    phvar <- grepl("^ph\\([0-9]+:[0-9]+\\)$", row.names(cmap))
    cmap2 <- cmap[!phvar,, drop=FALSE]
    Xna <- is.na(model.matrix(delete.response(Terms), data=mf))
    Xna <- Xna[,-1, drop=FALSE] # delete the intercept, has to be done last
    bad <- good <- rep(FALSE, nrow(Y))
    for (i in 1:ncol(cmap)) {
        j <- which(cmap2[,i] >0)
        atrisk <- (as.integer(istate) == from[i])
        bad[atrisk] <- bad[atrisk] | apply(Xna[atrisk, j, drop=FALSE], 1, any)
        good[atrisk]<- good[atrisk] | apply(!Xna[atrisk,j,drop=FALSE], 1, all)
    }

    # good[i] is TRUE if observation i makes a contribution to some transition
    # bad [i] is TRUE observation i contributes a missing for some transition

    # last, deal with any missing in the strata, which follow the same rules
    # as covariates. Strata in multistate models are unusual, so we don't need
    # to be efficient
    smap <- formula$smap
    if (nrow(smap) >1) {
        for (i in 2:nrow(smap)) {
            stemp <- is.na(mf[[rownames(smap)[i]]])
            for (j in 1:ncol(smap)) { # same cols as cmap, one per trans
                atrisk <- (as.integer(istate) == from[i])
                if (smap[2,j]>0) {
                    bad[atrisk] <- bad[atrisk] | stemp[atrisk]
                    good[atrisk] <-good[atrisk]| !stemp[atrisk]
                }
            }
        }
    }
    n.partially.used <- sum(good & bad & !miss0) # unused, but interesting
    omit <- (!good & bad) | miss0

    if (any(omit)) {
        temp <- setNames(seq(omit)[omit], attr(mf, "row.names")[omit])
        attr(temp, "class") <- "omit"
        mf <- mf[!omit,, drop=FALSE]
        attr(mf, "na.action") <- temp
    }
    mf
}

