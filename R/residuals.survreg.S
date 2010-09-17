# Automatically generated from all.nw using noweb
# $Id$
#
#  Residuals for survreg objects
residuals.survreg <- function(object, type=c('response', 'deviance',
                      'dfbeta', 'dfbetas', 'working', 'ldcase',
                      'ldresp', 'ldshape', 'matrix'), 
                      rsigma =TRUE, collapse=FALSE, weighted=FALSE, ...) {
    type <-match.arg(type)
    n <- length(object$linear.predictors)
    Terms <- object$terms
    if(!inherits(Terms, "terms"))
            stop("invalid terms component of  object")

    # If there was a cluster directive in the model statment then remove
    #  it.  It does not correspond to a coefficient, and would just confuse
    #  things later in the code.
    cluster <- untangle.specials(Terms,"cluster")$terms
    if (length(cluster) >0 )
        Terms <- Terms[-cluster]

    strata <- attr(Terms, 'specials')$strata
    coef <- object$coefficients
    intercept <- attr(Terms, "intercept") 
    response  <- attr(Terms, "response")
    weights <- object$weights
    if (is.null(weights)) weighted <- FALSE

    if (is.character(object$dist)) 
                dd <- survreg.distributions[[object$dist]]
    else dd <- object$dist
    if (is.null(dd$itrans)) {
            itrans <- dtrans <-function(x)x
            }
    else {
            itrans <- dd$itrans
            dtrans <- dd$dtrans
            }
    if (!is.null(dd$dist))  dd <- survreg.distributions[[dd$dist]]
    deviance <- dd$deviance
    dens <- dd$density
    if (is.null(object$naive.var)) vv <- object$var
    else                           vv <- object$naive.var

    need.x <- is.na(match(type, c('response', 'deviance', 'working')))
    if (is.null(object$y) || !is.null(strata) || (need.x & is.null(object[['x']])))
        mf <- model.frame(object)

    y <- object$y
    if (is.null(y)) {
        y <- model.extract(mf, 'response')
        if (!is.null(dd$trans)) {
            tranfun <- dd$trans
            exactsurv <- y[,ncol(y)] ==1
            if (any(exactsurv)) logcorrect <-sum(log(dd$dtrans(y[exactsurv,1])))

            if (type=='interval') {
                if (any(y[,3]==3))
                        y <- cbind(tranfun(y[,1:2]), y[,3])
                else y <- cbind(tranfun(y[,1]), y[,3])
                }
            else if (type=='left')
                 y <- cbind(tranfun(y[,1]), 2-y[,2])
            else     y <- cbind(tranfun(y[,1]), y[,2])
            }
        else {
            if (type=='left') y[,2] <- 2- y[,2]
            else if (type=='interval' && all(y[,3]<3)) y <- y[,c(1,3)]
            }
        }

    if (!is.null(strata)) {
        temp <- untangle.specials(Terms, 'strata', 1)
        Terms2 <- Terms[-temp$terms]
        if (length(temp$vars)==1) strata.keep <- mf[[temp$vars]]
        else strata.keep <- strata(mf[,temp$vars], shortlabel=TRUE)
        strata <- as.numeric(strata.keep)
        nstrata <- max(strata)
        sigma <- object$scale[strata]
        }
    else {
        Terms2 <- Terms
        nstrata <- 1
        sigma <- object$scale
        }
            
    if (need.x) { 
       x <- object[['x']]  #don't grab xlevels component
       if (is.null(x)) 
            x <- model.matrix(Terms2, mf, contrasts.arg=object$contrasts)
        }
    if (type=='response') {
        yhat0 <- deviance(y, sigma, object$parms)
        rr <-  itrans(yhat0$center) - itrans(object$linear.predictor)
        }
    else {
        status <- y[,ncol(y)]
        eta <- object$linear.predictors
        z <- (y[,1] - eta)/sigma
        dmat <- dens(z, object$parms)
        dtemp<- dmat[,3] * dmat[,4]    #f'
        if (any(status==3)) {
            z2 <- (y[,2] - eta)/sigma
            dmat2 <- dens(z2, object$parms)
            }
        else {
            dmat2 <- dmat   #dummy values
            z2 <- 0
            }

        tdenom <- ((status==0) * dmat[,2]) +  #right censored
                  ((status==1) * 1 )       +  #exact
                  ((status==2) * dmat[,1]) +  #left
                  ((status==3) * ifelse(z>0, dmat[,2]-dmat2[,2], 
                                             dmat2[,1] - dmat[,1])) #interval
        g <- log(ifelse(status==1, dmat[,3]/sigma, tdenom))  #loglik
        tdenom <- 1/tdenom
        dg <- -(tdenom/sigma) *(((status==0) * (0-dmat[,3])) +    #dg/ eta
                                ((status==1) * dmat[,4]) +     
                                ((status==2) * dmat[,3]) +      
                                ((status==3) * (dmat2[,3]- dmat[,3])))

        ddg <- (tdenom/sigma^2) *(((status==0) * (0- dtemp)) +  #ddg/eta^2
                                  ((status==1) * dmat[,5]) +
                                  ((status==2) * dtemp) +
                                  ((status==3) * (dmat2[,3]*dmat2[,4] - dtemp))) 

        ds  <- ifelse(status<3, dg * sigma * z,
                                tdenom*(z2*dmat2[,3] - z*dmat[,3]))
        dds <- ifelse(status<3, ddg* (sigma*z)^2,
                                tdenom*(z2*z2*dmat2[,3]*dmat2[,4] -
                                        z * z*dmat[,3] * dmat[,4]))
        dsg <- ifelse(status<3, ddg* sigma*z,
                      tdenom *(z2*dmat2[,3]*dmat2[,4] - z*dtemp))
        deriv <- cbind(g, dg, ddg=ddg- dg^2, 
                       ds = ifelse(status==1, ds-1, ds), 
                       dds=dds - ds*(1+ds), 
                       dsg=dsg - dg*(1+ds))
        if (type=='deviance') {
            yhat0 <- deviance(y, sigma, object$parms)
            rr <- (-1)*deriv[,2]/deriv[,3]  #working residuals
            rr <- sign(rr)* sqrt(2*(yhat0$loglik - deriv[,1]))
            }

        else if (type=='working') rr <- (-1)*deriv[,2]/deriv[,3]

        else if (type=='dfbeta' || type== 'dfbetas' || type=='ldcase') {
            score <- deriv[,2] * x  # score residuals
            if (rsigma) {
                if (nstrata > 1) {
                    d4 <- matrix(0., nrow=n, ncol=nstrata)
                    d4[cbind(1:n, strata)] <- deriv[,4]
                    score <- cbind(score, d4)
                    }
                else score <- cbind(score, deriv[,4])
                }
            rr <- score %*% vv
            if (type=='dfbetas') rr <- rr %*% diag(1/sqrt(diag(vv)))
            if (type=='ldcase')  rr<- rowSums(rr*score)
            }

        else if (type=='ldresp') {
            rscore <-  deriv[,3] *  (x * sigma)
            if (rsigma) {
                if (nstrata >1) {
                    d6 <- matrix(0., nrow=n, ncol=nstrata)
                    d6[cbind(1:n, strata)] <- deriv[,6]*sigma
                    rscore <- cbind(rscore, d6)
                    }
                else rscore <- cbind(rscore, deriv[,6] * sigma)
                }
            temp <-  rscore %*% vv
            rr <- rowSums(rscore * temp)
            }

        else if (type=='ldshape') {
            sscore <- deriv[,6] *x
            if (rsigma) {
                if (nstrata >1) {
                    d5 <- matrix(0., nrow=n, ncol=nstrata)
                    d5[cbind(1:n, strata)] <- deriv[,5]
                    sscore <- cbind(sscore, d5)
                    }
                else sscore <- cbind(sscore, deriv[,5])
                }
            temp <- sscore %*% vv
            rr <- rowSums(sscore * temp)
            }

        else {  #type = matrix
            rr <- deriv
            }
        }
    #case weights
    if (weighted) rr <- rr * weights

    #Expand out the missing values in the result
    if (!is.null(object$na.action)) {
        rr <- naresid(object$na.action, rr)
        if (is.matrix(rr)) n <- nrow(rr)
        else               n <- length(rr)
        }

    # Collapse if desired
    if (!missing(collapse)) {
        if (length(collapse) !=n) stop("Wrong length for 'collapse'")
        rr <- drop(rowsum(rr, collapse))
        }

    rr
    }
