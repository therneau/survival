# The vcov and extractAIC methods are not defined in Splus, so they
#  do not need survival methods.
if (is.R()) {
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
}
