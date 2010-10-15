pyears2html <- function(object, 
                   filename=paste(deparse(substitute(object)),'.html',sep=''),
                   title=NULL, exp=T, r=F, rr=T, ci.rr = F, 
                   conf.level= 0.95, digits=2) {

  if(class(object) != 'pyears') 
      stop('Your object does not have class "pyears"')

  nvar <- length(dim(object$n)) #number of categorization variables
  if(nvar>4) 
    stop('Currently this function can handle a maximum of 4 variables in the pyears model')
  
  ## If no expected rates are found in the object, then automatically set
  ## exp, rr, and ci.rr options to FALSE
  foundExpected <- !(is.na(match('expected', names(object))))
  if(!foundExpected) {
    exp <- F
    rr <- F
    ci.rr <- F
  }
  
  ## create an initial table that includes a legend.
  legend <- '<b>Cell Content:</b><p>N Subjects<br># Person-years<br>#Observed Events'
  if(exp) legend <- paste(legend,'<br># Expected Events')
  if(r) legend <- paste(legend,'<br>Obs/Exp')
  if(rr) legend <- paste(legend,'<br>Risk-Ratio')
  if(ci.rr) legend <- paste(legend,'<br>RR ',conf.level*100,'% CI',sep='')
  
  html.table(legend, file=filename, main=title)

  ## add in some spacing before introducing the next table
  sink(filename, append=T)
  cat('<p>&nbsp;</p>')
  sink()

  ## IF THERE IS ONLY 1 VARIABLE IN THE PYEARS CALL
  
  if(nvar==1) {

    ## 1 variable (vectors) ##
    n.mat <- c(object$n,' ') ## don't know actual N overall
    p.mat <- c(object$pyears, sum(object$pyears, na.rm=T))
    obs.mat <- c(object$event, sum(object$event, na.rm=T))
    r.mat <- obs.mat/p.mat

    if(exp){
      exp.mat <- c(object$expected, sum(object$expected, na.rm=T))
      rr.mat <- obs.mat/exp.mat
    }

    if(ci.rr){
      ## stick together the confidence interval
      tmp.ci <- cipoisson(obs.mat, exp.mat, p=conf.level)
      lower.mat <- round(tmp.ci[,1], digits)
      upper.mat <- round(tmp.ci[,2], digits)
      ci.mat <- paste(lower.mat,'-', upper.mat, sep='')
    }

    ## create a matrix that is a character string
    tmp <- matrix(paste(n.mat,'<br>',
                        round(p.mat,digits),'<br>',
                        obs.mat,'<br>',
                        if(r) paste(round(r.mat,digits),'<br>'),
                        if(exp) paste(round(exp.mat,digits),'<br>'),
                        if(rr) paste(round(rr.mat,digits),'<br>'),
                        if(ci.rr) ci.mat,
                        sep=''), nrow=length(n.mat), ncol=1)

    dimnames(tmp) <- list(c(names(object$n), 'Total'), NULL)
    html.table(tmp, file=filename, append=T)
    
  }

  if(nvar==2) {

    ## 2 variables (everything is a matrix) ##
    n.mata <- cbind(object$n,' ')
    n.mat <- rbind(n.mata,' ')
    p.mata <- cbind(object$pyears,rowSums(object$pyears))
    p.mat <- rbind(p.mata,colSums(p.mata))
    obs.mata <- cbind(object$event,rowSums(object$event))
    obs.mat <- rbind(obs.mata,colSums(obs.mata))
    r.mat <- obs.mat/p.mat
    
    if(exp){
      exp.mata <- cbind(object$expected,rowSums(object$expected))
      exp.mat <- rbind(exp.mata,colSums(exp.mata))
      rr.mat <- obs.mat/exp.mat
    }

    if(ci.rr){
      ## stick together the confidence interval
      tmp.ci <- cipoisson(obs.mat, exp.mat, p=conf.level)
      lower.mat <- matrix(round(tmp.ci[,1], digits), nrow=nrow(n.mat))
      upper.mat <- matrix(round(tmp.ci[,2], digits), nrow=nrow(n.mat))
      ci.mat <- paste(lower.mat,'-', upper.mat, sep='')
    }
    
    tmp <- matrix(paste(n.mat,'<br>',
                        round(p.mat,digits),'<br>',
                        obs.mat,'<br>',
                        if(r) paste(round(r.mat,digits),'<br>'),
                        if(exp) paste(round(exp.mat,digits),'<br>'),
                        if(rr) paste(round(rr.mat,digits),'<br>'),
                        if(ci.rr) ci.mat,
                        sep=''), nrow=nrow(n.mat), ncol=ncol(n.mat))

    dimnames(tmp) <- list(c(dimnames(object$n)[[1]], 'Total'),
                          c(dimnames(object$n)[[2]], 'Total'))

    html.table(tmp, file=filename, append=T)
  }
  

  if(nvar == 3) {
    tmpall <- NA
    
    ## loop over the last variable
    for(i in 1:(dim(object$n)[3])) {
      dimtitle <- dimnames(object$n)[[3]] [i]

      n <- object$n[,,i]
      pyears <- object$pyears[,,i] 
      event <- object$event[,,i]

      ## 2 variables (everything is a matrix) ##
      n.mata <- cbind(n,' ')
      n.mat <- rbind(n.mata,' ')
      p.mata <- cbind(pyears,rowSums(pyears))
      p.mat <- rbind(p.mata,colSums(p.mata))
      obs.mata <- cbind(event,rowSums(event))
      obs.mat <- rbind(obs.mata,colSums(obs.mata))
      r.mat <- obs.mat/p.mat
      
      if(exp){
        expected <- object$expected[,,i]
        exp.mata <- cbind(expected,rowSums(expected))
        exp.mat <- rbind(exp.mata,colSums(exp.mata))
        rr.mat <- obs.mat/exp.mat
      }

      if(ci.rr){
        ## stick together the confidence interval
        tmp.ci <- cipoisson(obs.mat, exp.mat, p=conf.level)
        lower.mat <- matrix(round(tmp.ci[,1], digits), nrow=nrow(n.mat))
        upper.mat <- matrix(round(tmp.ci[,2], digits), nrow=nrow(n.mat))
        ci.mat <- paste(lower.mat,'-', upper.mat, sep='')
      }

      
      tmp <- matrix(paste(n.mat,'<br>',
                          round(p.mat,digits),'<br>',
                          obs.mat,'<br>',
                          if(r) paste(round(r.mat,digits),'<br>'),
                          if(exp) paste(round(exp.mat,digits),'<br>'),
                          if(rr) paste(round(rr.mat,digits),'<br>'),
                          if(ci.rr) ci.mat,
                          sep=''), nrow=nrow(n.mat), ncol=ncol(n.mat))

      dimnames(tmp) <- list(c(dimnames(n)[[1]], 'Total'),
                            c(dimnames(n)[[2]], 'Total'))
      
      html.table(tmp, file=filename, main=paste('Group=',dimtitle), append=T)
      tmpall <- list(tmpall, tmp)
    }
    tmp <- tmpall
  }


  if(nvar == 4) {

    tmpall <- NA
    ## loop over the last variable
    for(i in 1:(dim(object$n)[4])) {
      dimtitle1 <- dimnames(object$n)[[4]] [i]
      ## loop over second to last variable
      for(j in 1:(dim(object$n)[3])) {
        dimtitle2 <- dimnames(object$n)[[3]] [j]

        n <- object$n[,,j,i]; pyears <- object$pyears[,,j,i] 
        event <- object$event[,,j,i]; 

        ## 2 variables (everything is a matrix) ##
        n.mata <- cbind(n,' ')
        n.mat <- rbind(n.mata,' ')
        p.mata <- cbind(pyears,rowSums(pyears))
        p.mat <- rbind(p.mata,colSums(p.mata))
        obs.mata <- cbind(event,rowSums(event))
        obs.mat <- rbind(obs.mata,colSums(obs.mata))
        r.mat <- obs.mat/p.mat
        
        if(exp){
          expected <- object$expected[,,j,i]
          exp.mata <- cbind(expected,rowSums(expected))
          exp.mat <- rbind(exp.mata,colSums(exp.mata))
          rr.mat <- obs.mat/exp.mat
        }

        if(ci.rr){
          ## stick together the confidence interval
          tmp.ci <- cipoisson(obs.mat, exp.mat, p=conf.level)
          lower.mat <- matrix(round(tmp.ci[,1], digits), nrow=nrow(n.mat))
          upper.mat <- matrix(round(tmp.ci[,2], digits), nrow=nrow(n.mat))
          ci.mat <- paste(lower.mat,'-', upper.mat, sep='')
        }

        
        tmp <- matrix(paste(n.mat,'<br>',
                            round(p.mat,digits),'<br>',
                            obs.mat,'<br>',
                            if(r) paste(round(r.mat,digits),'<br>'),
                            if(exp) paste(round(exp.mat,digits),'<br>'),
                            if(rr) paste(round(rr.mat,digits),'<br>'),
                            if(ci.rr) ci.mat,
                            sep=''), nrow=nrow(n.mat), ncol=ncol(n.mat))

        dimnames(tmp) <- list(c(dimnames(n)[[1]], 'Total'),
                              c(dimnames(n)[[2]], 'Total'))
        
        html.table(tmp, file=filename, main=paste('Group 1=',dimtitle1, '   Group 2=',
                                         dimtitle2), append=T)
        tmpall <- list(tmpall, tmp)
      }
    }
    tmp <- tmpall
  }

  invisible(tmp)
}
