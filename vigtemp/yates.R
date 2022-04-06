#
# A function to compute Yate's contrasts
# For explanation of how it works, see the vignette 
#  "Cox models and type III tests", particularly appendix 1
#
yates <- function(fit, method=c("ATT", "STT", "NSTT")) {
    method <- match.arg(method)

    X <- model.matrix(fit)
    ax <- attr(X, "assign")
    xcon <- attr(X, "contrasts")
    Terms <- terms(fit)
    dfac   <- attr(Terms, "factors")     #which variables are in which term

    # Ignore any columns corresponding to NA coefficients
    #  by assigning them to "intercept"
    nacoef <- is.na(coef(fit))
    if (any(nacoef)) {
        ax[nacoef] <- 0
        X[,nacoef] <- 0
    }
    ux <- unique(ax)
    ux <- ux[ux!=0]  #no contrast is returned for the intercept
    if (length(ux)==0) return(NULL)  #all coefs are NA
    clist <- vector("list", length(ux))
    names(clist) <- dimnames(dfac)[[2]]

    # Create the default contrasts for each term
    dtemp <- diag(length(ax))
    dimnames(dtemp) <- list(NULL, dimnames(X)[[2]])
    for (i in 1:length(ux)) 
            clist[[i]] <- dtemp[ax== ux[i],]
        
    # Divide the X matrix into factors and non-factors
    # By definition a factor is something listed in attr(X, "contrasts")
    # Any factor's contrast is ajusted for other factor terms that
    #  contain it.
    # The matrix "included" has a row for each term, 1= term + extenders
    forder <- attr(Terms, "order")
    fvar <- dimnames(dfac)[[1]] %in% names(xcon)
    fterm <- colSums(dfac[!fvar,, drop=FALSE]) ==0  #terms involving only factors
    included <- diag(ncol(dfac)) #one row for each term, it + others
    for (i in which(fterm)) {
        temp <- dfac[,i] * dfac
        keep <- (colSums(temp) >0 & forder >=forder[i] & fterm)
        included[i,keep] <- 1
    }
    adjust <- which(rowSums(included) >1) #factors with an including interaction
    if (length(adjust) >0 && method != "NSTT") {
        # Use one of the more sophisticated methods
 
        if (method=="ATT") {
            # Averaging method
            for (i in adjust) {
                icol <- ax %in% which(included[i,] >0)
                xfac <- unique(X[, icol, drop=FALSE])
                afac <- ax[icol]
                           
                xtemp <- xfac[,afac==i, drop=FALSE]  #effect of interest
                index <- xmatch(xtemp)
                dtemp <- matrix(0., max(index), ncol(X),
                            dimnames=list(NULL, dimnames(X)[[2]]))
                for (j in 1:max(index)) 
                    dtemp[j, icol] <- colMeans(xfac[index==j,, drop=FALSE])
                clist[[i]] <- dtemp[-1,] - rep(dtemp[1,], each=nrow(dtemp)-1)
            }            
        } else {
        #SAS type III method, first build the D matrix
        }       
     }

    #perform the contrasts
    temp <- lapply(clist, contrast, fit)
    # pull out the SS and df, since that's what folks most often want
    list(ss = sapply(temp, function(x) x$ss),
         df = sapply(temp, function(x) x$df),
         detail = temp)
}

xmatch <- function(x) {
    # match each row of a matrix x to unique(x)
    # there may be a fast, clever, non loop way to do this,
    #  but these matrices are very small so it's not worth subtlety
    ux <- unique(x)
    result <- integer(nrow(x))
    for (i in 1:nrow(x)){
        for (j in 1:nrow(ux)) if (all(x[i,] == ux[j,])) result[i] <- j
    }
    result
}
qform <- function(beta, var) # quadratic form b' (V-inverse) b
    sum(beta * solve(var, beta))

contrast <- function(cmat, fit) {
    varmat <- vcov(fit)
    if (class(fit) == "lm") sigma2 <- summary(fit)$sigma^2              
    else sigma2 <- 1   # for the Cox model case

    beta <- coef(fit)
    if (!is.matrix(cmat)) cmat <- matrix(cmat, nrow=1)
    if (ncol(cmat) != length(beta)) stop("wrong dimension for contrast")

    estimate <- drop(cmat %*% beta)  #vector of contrasts
    ss <- qform(estimate, cmat %*% varmat %*% t(cmat)) *sigma2
    list(contrast=cmat, estimate=estimate, ss=ss, df=qr(cmat)$rank,
         var=drop(cmat %*% varmat %*% t(cmat)))
    }
