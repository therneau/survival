vcov.coxph<-function (object, ...) {
        rval<-object$var
        dimnames(rval)<-list(names(coef(object)),names(coef(object)))
        rval
    }

vcov.survreg<-function (object, ...) {
        object$var
    }


extractAIC.coxph.penal<- function(fit,scale,k=2,...){
        edf<-sum(fit$df)
        loglik <- fit$loglik[length(fit$loglik)]
        c(edf, -2 * loglik + k * edf)
    }


labels.survreg <- function(object, ...) attr(object,"term.labels")

rep.Surv <- function(x, ...) {
    indx <- rep(1:nrow(x), ...)
    x[indx,]
}

