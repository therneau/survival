## conditional logistic regression
##
## case~exposure+strata(matching)
##

clogit<-function(formula,data,method=c("exact","approximate"),
                 na.action=getOption("na.action"),subset=NULL,
                 control=coxph.control()){
    
    mf<-match.call()
    mf[[1]]<-as.name("model.frame")
    mf$method<-mf$control<-NULL
    mfn<-mf

    mfn$na.action<-"I"
    mfn$subset<-NULL
    nrows<-NROW(eval(mfn,parent.frame()))

    mf<-eval(mf,parent.frame())
    
    coxcall<-match.call()
    coxcall[[1]]<-as.name("coxph")
    newformula<-formula
    newformula[[2]]<-substitute(Surv(rep(1,nn),case),list(case=formula[[2]],nn=nrows))
    environment(newformula)<-environment(formula)
    coxcall$formula<-newformula
    coxcall$method<-switch(match.arg(method),exact="exact","breslow")

    coxcall<-eval(coxcall,sys.frame(sys.parent()))
    coxcall$userCall<-sys.call()
    
    class(coxcall)<-c("clogit","coxph")
    coxcall
}


print.clogit<-function(x,...){

    x$call<-x$userCall
    NextMethod()

}
