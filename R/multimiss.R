# The complex NA manipuation needed for multistate coxph models that have
#  a list of formulas.  See the section on multistate missing in the
#  methods document
multimiss <- function(mf, y, id, istate, tmap) {
    # first off: any row missing y, id, or weight will be lost
    miss0 <- is.na(y) | is.na(id) # y and id are not optional
    if (!is.null(istate)) miss0 <- miss0 | is.na(istate)
    if (!is.null(istate)) miss0 <- miss0 | is.na(istate)
    weights <- model.weights(mf)
    if (!is.null(weights)) miss0 <- miss0 | is.na(weights)
    cluster <- model.extract(mf, "cluster")
    if (!is.null(cluster)) miss0 <- miss0 | is.na(cluster)

    # the starting state for each transition
    from  <- as.numeric(sub(":[0-9]+$", "", colnames(tmap)))

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
    # matrix that shows if row i of the data contibutes a missing to transition j
    makemiss <- termiss %*% tmap2 # positive only if missing and used
    # element i,j of makemiss =1 if objs i will make linear predictor j missing
    #  only rows that have no utility (make all LP missing) are tossed
    omit <- miss0 | (rowSums(makemiss) == ncol(makemiss))
    omit
}
