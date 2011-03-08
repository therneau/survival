# $Id 
summary.coxph <- function(object,  conf.int = 0.95, scale = 1, ...) {
    cox<-object
    beta <- cox$coefficients
    if (is.null(cox$coefficients)) {   # Null model
	return(object)  #The summary method is the same as print in this case
	}
    nabeta <- !(is.na(beta))          #non-missing coefs
    beta2 <- beta[nabeta]
    if(is.null(beta) | is.null(cox$var))
         stop("Input is not valid")
    se <- sqrt(diag(cox$var))
    if (!is.null(cox$naive.var)) nse <- sqrt(diag(cox$naive.var))

    rval<-list(call=cox$call,fail=cox$fail, na.action=cox$na.action,
                n=cox$n, loglik=cox$loglik)
    if (!is.null(cox$nevent)) rval$nevent <- cox$nevent

    if (is.null(cox$naive.var)) {
        tmp <- cbind(beta, exp(beta), se, beta/se,
                     1 - pchisq((beta/ se)^2, 1))
        dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
                                             "se(coef)", "z", "Pr(>|z|)"))
        }
    else {
        tmp <- cbind(beta, exp(beta), nse, se, beta/se,
                     1 - pchisq((beta/ se)^2, 1))
        dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
                               "se(coef)", "robust se", "z", "Pr(>|z|)"))
        }
    rval$coefficients <- tmp
     
    if (conf.int) {
        z <- qnorm((1 + conf.int)/2, 0, 1)
        beta <- beta * scale
        se <- se * scale
        tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se),
            exp(beta + z * se))
        dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
            paste("lower .", round(100 * conf.int, 2), sep = ""),
            paste("upper .", round(100 * conf.int, 2), sep = "")))
        rval$conf.int <- tmp
        }

    df <- length(beta2)
    logtest <- -2 * (cox$loglik[1] - cox$loglik[2])
    rval$logtest <- c(test=logtest,
                      df=df,
                      pvalue=1 - pchisq(logtest, df))
    rval$sctest <- c(test=cox$score,
                     df=df,
                     pvalue=1 - pchisq(cox$score, df))
    rval$rsq<-c(rsq=1-exp(-logtest/cox$n),
                 maxrsq=1-exp(2*cox$loglik[1]/cox$n))
    rval$waldtest<-c(test=as.vector(round(cox$wald.test, 2)),
                      df=df,
                      pvalue=1 - pchisq(as.vector(cox$wald.test), df))
    if (!is.null(cox$rscore))
         rval$robscore<-c(test=cox$rscore,
                          df=df,
                          pvalue=1 - pchisq(cox$rscore, df))
    rval$used.robust<-!is.null(cox$naive.var)

    if (!is.null(cox$concordance)) {
        if (is.matrix(cox$concordance)) temp <- colSums(cox$concordance)
        else temp <- cox$concordance
        rval$concordance <- c("concordance"= (temp[1] + temp[3]/2)/
                              sum(temp[1:3]), "se"= temp[5]/(2*sum(temp[1:3])))
    }
        

    if (is.R()) class(rval)    <-"summary.coxph"
    else        oldClass(rval) <- "summary.coxph"
    rval
    }
