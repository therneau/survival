### Suite of programs for case-cohort analysis
### Main program

cch <- function(formula, data=sys.parent(), subcoh, id, stratum=NULL, cohort.size, 
                method=c("Prentice", "SelfPrentice", "LinYing","I.Borgan","II.Borgan"),
                robust=FALSE){
    call <- match.call()
    
    if (is.data.frame(data)){
        if (inherits(id,"formula"))
            id<-model.frame(id,data,na.action=na.fail)[,1]
        if (inherits(subcoh,"formula"))
            subcoh<-model.frame(subcoh,data,na.action=na.fail)[,1]
        if (inherits(stratum,"formula"))
            stratum<-model.frame(stratum,data,na.action=na.fail)[,1]
    }

    ## Check id, subcoh and cohort.size variables
    if(length(id)!=length(unique(id)))
        stop("Multiple records per id not allowed")
    if (is.logical(subcoh))
        subcoh <- as.numeric(subcoh)
    tt <- table(subcoh)
    if(min(charmatch(names(tt), c("0","1"), 0))==0)
        stop("Permissible values for subcohort indicator are 0/1 or TRUE/FALSE")
    if(length(id)>sum(cohort.size))
        stop("Number of records greater than cohort size")
    nn <- cohort.size	
    method<-match.arg(method)
    stratified<-method %in% c("I.Borgan","II.Borgan")
    if (!is.null(stratum))
        stratum<-factor(stratum)
    if (stratified){
        if (robust)
            warning("`robust' not implemented for stratified analysis.")
        if (is.null(stratum))
            stop("method (",method,") requires 'stratum'")
        if (length(cohort.size)!=length(levels(stratum)))
            stop("cohort.size and stratum do not match")
        if (!(all(levels(stratum) %in% names(cohort.size))))
            warning("stratum levels and names(cohort.size) do not agree")
        subcohort.sizes<-table(stratum)
    } else if(!stratified) {
        if (!(method =="LinYing"))
            warning("`robust' ignored for  method (",method,")")
        if (!is.null(stratum))
            warning("'stratum' ignored for method (",method,")")
        if (length(cohort.size)!=1)
            stop("cohort size must be a scalar for unstratified analysis")
        subcohort.sizes<-length(id)
    }
    if (any(subcohort.sizes>cohort.size))
        stop("Population smaller than sample in some strata")
    ## Evaluate model formula
    m <- match.call(expand.dots=FALSE)
    m$method <- m$cohort.size <- m$id <- m$subcoh <- m$stratum <-m$robust<- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m,sys.parent())
    Terms <- attr(m,"terms")
    Y <- model.extract(m, "response")
    if(!inherits(Y, "Surv"))
        stop("Response must be a survival object")
    type <- attr(Y, "type")
    itype<-charmatch(type,c("right","counting"),nomatch=0)
    cens<-switch(itype+1,
                 stop(paste("Cox model doesn't support \"", type, "\" survival data", sep = "")),
                 Y[,2],
                 Y[,3])
    if (any(!subcoh & !cens))
        stop(sum(!subcoh & !cens),"censored observations not in subcohort")
    cc<-cens+1-subcoh
    texit<-switch(itype+1, stop(), Y[,1], Y[,2])
    tenter<-switch(itype+1, stop(), rep(0,length(texit)), Y[,1])
    X <- model.matrix(Terms, m)
    X <- X[,2:ncol(X)]
    fitter <- get(method)
    if (stratified)        
        out<-fitter(tenter=tenter, texit=texit, cc=cc, id=id, X=X,
                    stratum=as.numeric(stratum), stratum.sizes=cohort.size)
    else
        out<-fitter(tenter=tenter, texit=texit, cc=cc, id=id, X=X, ntot=nn, robust=robust)
    out$method <- method
    names(out$coefficients) <- dimnames(X)[[2]]
    if(!is.null(out$var))
        dimnames(out$var) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
    if(!is.null(out$naive.var))
        dimnames(out$naive.var) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
    out$call <- call
    out$cohort.size <- cohort.size
    out$stratified<-stratified
    if (stratified){
        out$stratum<-stratum
        out$subcohort.size <-subcohort.sizes
    } else {
        out$subcohort.size <- tt[2]
    }
    class(out) <- "cch"
    out
}

### Subprograms

Prentice <- function(tenter, texit, cc,  id, X, ntot,robust){
    eps <- 0.00000001
    cens <- as.numeric(cc>0) # Censorship indicators
    subcoh <- as.numeric(cc<2) # Subcohort indicators

    ## Calculate Prentice estimate
    ent2 <- tenter
    ent2[cc==2] <- texit[cc==2]-eps
    fit1 <- coxph(Surv(ent2,texit,cens)~X,eps=eps,x=TRUE)

    ## Calculate Prentice estimate and variance
    nd <- sum(cens) # Number of failures
    nc <- sum(subcoh) # Number in subcohort
    ncd <- sum(cc==1) #Number of failures in subcohort
    X <- as.matrix(X)
    aent <- c(tenter[cc>0],tenter[cc<2])
    aexit <- c(texit[cc>0],texit[cc<2])
    aX <- rbind(as.matrix(X[cc>0,]),as.matrix(X[cc<2,]))
    aid <- c(id[cc>0],id[cc<2])
    dum <- rep(-100,nd)
    dum <- c(dum,rep(0,nc))
    gp <- rep(1,nd)
    gp <- c(gp,rep(0,nc))
    fit <- coxph(Surv(aent,aexit,gp)~aX+offset(dum)+cluster(aid),eps=eps,x=TRUE,
                 iter.max=35,init=fit1$coefficients)
    db <- resid(fit,type="dfbeta")
    db <- as.matrix(db)
    db <- db[gp==0,]
    fit$phase2var<-(1-(nc/ntot))*t(db)%*%(db)
    fit$naive.var <- fit$naive.var+fit$phase2var
    fit$var<-fit$naive.var
    
    fit$coefficients <- fit$coef <- fit1$coefficients
    fit
}

SelfPrentice <- function(tenter, texit, cc,  id, X, ntot,robust){
    eps <- 0.00000001
    cens <- as.numeric(cc>0) # Censorship indicators
    subcoh <- as.numeric(cc<2) # Subcohort indicators

    ## Calculate Self-Prentice estimate and variance
    nd <- sum(cens) # Number of failures
    nc <- sum(subcoh) # Number in subcohort
    ncd <- sum(cc==1) #Number of failures in subcohort
    X <- as.matrix(X)
    aent <- c(tenter[cc>0],tenter[cc<2])
    aexit <- c(texit[cc>0],texit[cc<2])
    aX <- rbind(as.matrix(X[cc>0,]),as.matrix(X[cc<2,]))
    aid <- c(id[cc>0],id[cc<2])
    dum <- rep(-100,nd)
    dum <- c(dum,rep(0,nc))
    gp <- rep(1,nd)
    gp <- c(gp,rep(0,nc))
    fit <- coxph(Surv(aent,aexit,gp)~aX+offset(dum)+cluster(aid),eps=eps,x=TRUE)
    db <- resid(fit,type="dfbeta")
    db <- as.matrix(db)
    db <- db[gp==0,,drop=FALSE]
    fit$phase2var<-(1-(nc/ntot))*t(db)%*%(db)
    fit$naive.var <- fit$naive.var+fit$phase2var
    fit$var<-fit$naive.var
    fit
}

LinYing <- function(tenter, texit, cc,  id, X, ntot,robust){
    eps <- 0.000000001
    cens <- as.numeric(cc>0) # Censorship indicators
    subcoh <- as.numeric(cc<2) # Subcohort indicators
    nd <- sum(cens) # Number of failures
    nc <- sum(subcoh) # Number in subcohort
    ncd <- sum(cc==1) #Number of failures in subcohort

    ## Calculate Lin-Ying estimate and variance
    offs <- rep((ntot-nd)/(nc-ncd),length(texit))
    offs[cc>0] <- 1
    loffs <- log(offs)
    fit <- coxph(Surv(tenter, texit, cens)~X+offset(loffs)+cluster(id),
                 eps=eps,x=TRUE)
    db <- resid(fit,type="dfbeta")
    db <- as.matrix(db)
    db0 <- db[cens==0,,drop=FALSE]
    dbm <- apply(db0,2,mean)
    db0 <- sweep(db0,2,dbm)
    fit$phase2var<-(1-(nc-ncd)/(ntot-nd))*crossprod(db0)
    fit$naive.var <- fit$naive.var+fit$phase2var
    if (robust)
        fit$var<- crossprod(db,db/offs)+fit$phase2var
    else
        fit$var<-fit$naive.var
    fit
}

I.Borgan <- function(tenter, texit, cc,  id, X, stratum, stratum.sizes){
  eps <- 0.00000001
  nobs <- length(texit)
  idx <- 1:length(nobs)
  jj <- max(stratum)
  nn <- stratum.sizes  ## Cohort stratum sizes
  n <- table(stratum)  ## Sample stratum sizes
  d <- table(stratum[cc>0]) ## Failures in each stratum
  tt <- table(cc,stratum)
  cens <- as.numeric(cc>0) ## Failure indicators
  subcoh <- as.numeric(cc<2) ## Subcohort indicators
  nd <- sum(cens) ## Number of failures
  nc <- sum(subcoh) ## Number in subcohort
  ncd <- sum(as.numeric(cc==1)) #Number of failures in subcohort
  m0 <- tt[1,] ## Subcohort stratum sizes (noncases only)
  if (ncd>0) m <- m0+tt[2,] else m <- m0 #Subcohort stratum sizes
  X <- as.matrix(X)
  kk <- ncol(X) ## Number of variables
  wt <- as.vector(nn/m) ## Weights for Estimator I
  stratum <- c(stratum[cc>0],stratum[cc<2])
  w <- wt[stratum]
  ent <- c(tenter[cc > 0], tenter[cc < 2])
  exit <- c(texit[cc > 0], texit[cc < 2])
  X <- rbind(as.matrix(X[cc > 0,  ]), 
             as.matrix(X[cc < 2,  ]))
  id <- c(id[cc > 0], id[cc < 2])
  dum <- rep(-100, nd)
  dum <- c(dum, rep(0, nc))
  gp <- rep(1, nd)
  gp <- c(gp, rep(0, nc))
  w[gp==1] <- 1
  fit <- coxph(Surv(ent,exit,gp)~X+offset(dum)+cluster(id),
               weights=w, eps=eps,x=T, iter.max=25)  
  score <- resid(fit, type = "score", weighted=F)
  sc <- resid(fit, type="score", collapse=id, weighted=T)
  score <- as.matrix(score)
  score <- score[gp == 0,,drop=F]
  st <- stratum[gp==0]
  sto <- st %o% rep(1,kk)
  Index <- col(score)
  tscore <- tapply(score,list(sto,Index),mean)
  pscore <- tapply(score,list(sto,Index))
  score <- score-tscore[pscore]
  delta <- matrix(0,kk,kk)
  opt <- NULL
  for (j in 1:jj) {
    temp <- t(score[st==j,])%*%score[st==j,]/(m[j]-1)
    delta <- matrix(delta+(wt[j]-1)*nn[j]*temp,kk,kk)
    if(is.null(opt)) 
      opt <- nn[j]*sqrt(diag(fit$naive.var %*% temp %*% fit$naive.var)) 
    else
      opt <- rbind(opt,nn[j]*sqrt(diag(fit$naive.var %*% temp %*% fit$naive.var))) 
  }
  z <- apply(opt,2,sum)
  fit$opt <- sweep(opt,2,z,FUN="/")
  fit$phase2var<-fit$naive.var%*%delta%*%fit$naive.var
  fit$naive.var <- fit$naive.var+fit$phase2var
  fit$var<-fit$naive.var
  fit$delta <- delta
  fit$sc <- sc
  fit
}

II.Borgan <- function(tenter, texit, cc,  id, X, stratum, stratum.sizes){
  eps <- 0.00000001
  jj <- max(stratum)
  nn <- stratum.sizes  ## Cohort stratum sizes
  n <- table(stratum)  ## Sample stratum sizes
  d <- table(stratum[cc>0]) ## Failures in each stratum
  tt <- table(cc,stratum)
  cens <- as.numeric(cc>0) ## Failure indicators
  subcoh <- as.numeric(cc<2) ## Subcohort indicators
  nd <- sum(cens) ## Number of failures
  nc <- sum(subcoh) ## Number in subcohort
  ncd <- sum(as.numeric(cc==1)) #Number of failures in subcohort
  m0 <- tt[1,] ## Subcohort stratum sizes (controls only)
  if (ncd>0) m <- m0+tt[2,] else m <- m0 #Subcohort stratum sizes
  X <- as.matrix(X)
  kk <- ncol(X) ## Number of variables
  nn0 <- nn-as.vector(d) #Noncases in cohort
  wt <- as.vector(nn0/m0)
  w <- wt[stratum]
  w[cens==1] <- 1
  fit <- coxph(Surv(tenter,texit,cens)~X+cluster(id),
               weights=w,eps=eps,x=T, iter.max=25)  ## Borgan Estimate II
  score <- resid(fit, type = "score", weighted=F)
  sc <- resid(fit,type="score", collapse=id, weighted=T)
  score <- as.matrix(score)
  score <- score[cens == 0,,drop=F] ## Scores for controls
  st <- stratum[cens==0] ## Stratum indicators for controls
  sto <- st %o% rep(1,kk)
  Index <- col(score)
  tscore <- tapply(score,list(sto,Index),mean) ## Within stratum control score means
  pscore <- tapply(score,list(sto,Index))
  score <- score-tscore[pscore] ## Subtract off within stratum score means
  delta <- matrix(0,kk,kk)
  opt <- NULL
  for (j in 1:jj) {
    temp <- t(score[st==j,])%*%score[st==j,]/(m0[j]-1) ## Borgan equation (19)
    delta <- delta+(wt[j]-1)*nn0[j]*temp ## Borgan equation (17)
    if(is.null(opt)) 
      opt <- nn0[j]*sqrt(diag(fit$naive.var %*% temp %*% fit$naive.var)) 
    else
      opt <- rbind(opt,nn0[j]*sqrt(diag(fit$naive.var %*% temp %*% fit$naive.var))) 
  }
  z <- apply(opt,2,sum)
  fit$opt <- sweep(opt,2,z,FUN="/")
  fit$phase2var<-fit$naive.var %*% delta %*% fit$naive.var
  fit$naive.var <- fit$naive.var+fit$phase2var
  fit$var<-fit$naive.var
  fit$delta <- delta
  fit$sc <- sc
  fit
}

## Methods

vcov.cch<-function(object,...) object$var

"print.cch"<- function(x,...)
{
    ## produces summary from an x of the class "cch"
    call<-x$call
    coef <- coef(x)
    method <- x$method
    se <- sqrt(diag(vcov(x)))
    Z<- abs(coef/se)
    p<- pnorm(Z)
    cohort.size<-x$cohort.size
    subcohort.size<-x$subcohort.size
    coefficients <- matrix(0, nrow = length(coef), ncol = 4)
    dimnames(coefficients) <- list(names(coef), c("Value", 
                                                  "SE", "Z", "p"))
    coefficients[, 1] <- coef
    coefficients[, 2] <- se
    coefficients[, 3] <- Z
    coefficients[, 4] <- 2*(1-p)

 
    if (x$stratified){
        cat("Exposure-stratified case-cohort analysis,", x$method, "method.\n")
        m<-rbind(subcohort=x$subcohort.size, cohort=x$cohort.size)
        prmatrix(m,quote=FALSE)
    } else{
        cat("Case-cohort analysis,")
        cat("x$method,", x$method,"\n with subcohort of",
            x$subcohort.size,"from cohort of", x$cohort.size,"\n\n")
    }
    cat("Call: "); print(x$call)
    cat("\nCoefficients:\n")
    print(coefficients)
    invisible(x)
}


"summary.cch"<-function(object,...)
{
    ## produces summary from an object of the class "cch"
    call<-object$call
	coef <- coef(object)
	method <- object$method
	se <- sqrt(diag(vcov(object)))
      Z<- abs(coef/se)
      p<- pnorm(Z)
      cohort.size<-object$cohort.size
      subcohort.size<-object$subcohort.size
    coefficients <- matrix(0, nrow = length(coef), ncol = 4)
    dimnames(coefficients) <- list(names(coef), c("Value", 
                                                  "SE", "Z", "p"))
    coefficients[, 1] <- coef
    coefficients[, 2] <- se
    coefficients[, 3] <- Z
    coefficients[, 4] <- 2*(1-p)
    structure(list(call=call, method=method, cohort.size=cohort.size,
                   subcohort.size=subcohort.size, coefficients = coefficients,
                   stratified=object$stratified), 
              class = "summary.cch")
}

print.summary.cch <- function(x,digits=3,...){

    if (x$stratified){
        cat("Exposure-stratified case-cohort analysis,", x$method, "method.\n")
        m<-rbind(subcohort=x$subcohort.size, cohort=x$cohort.size)
        prmatrix(m,quote=FALSE)
    } else{
        cat("Case-cohort analysis,")
        cat("x$method,", x$method,"\n with subcohort of",
            x$subcohort.size,"from cohort of", x$cohort.size,"\n\n")
    }
    cat("Call: "); print(x$call)
    cat("\nCoefficients:\n")
    output<-cbind(Coef=x$coefficients[,1],HR=exp(x$coefficients[,1]),
                  "(95%"=exp(x$coefficients[,1]-1.96*x$coefficients[,2]),
                  "CI)"=exp(x$coefficients[,1]+1.96*x$coefficients[,2]),
                  "p"=x$coefficients[,4]
                  )
    print(round(output,3))
    invisible(x)
}

