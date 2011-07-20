#
# The test for survival proposed by Peter O'Brien
#
survobrien <- function(formula, data= sys.frame(sys.parent()),
                       transform=function(r) (r-.5)/(.5+length(r)-r)) {

    unprotect<-function(Terms){
        namei<-as.name("I")
        namefactor<-as.name("factor")
        unprotect1<-function(aterm){
            if (is.name(aterm))
                return(aterm)
            if (is.name(aterm[[1]])){
                if (aterm[[1]]==namei | aterm[[1]]==namefactor)
                    return(aterm[[2]])
                if (length(aterm)==1)
                    return(aterm)
                for (i in 2:length(aterm)){
                    aterm[[i]]<-unprotect1(aterm[[i]])
                }
            }
            aterm
        }
        if (length(Terms)==3)
            Terms[[3]]<-unprotect1(Terms[[3]])
        else
            Terms[[2]]<-unprotect1(Terms[[2]])
        terms.formula(Terms)
    }
    
    
  
    m <- model.frame(formula, data, na.action= function(x) x )
    n <- nrow(m)
    Terms <- attr(m, 'terms')

    y <- model.extract(m, 'response')
    if (!inherits(y, "Surv")) stop ("Response must be a survival object")
    if (attr(y, 'type') != 'right') stop("Can only handle right censored data")

    # Figure out which are the continuous predictor variables
    factors <- unlist(lapply(m, is.factor))
    protected <- unlist(lapply(m, function(x) inherits(x, "AsIs")))
    keepers <- factors | protected
    cont <- ((seq(keepers))[!keepers]) [-1]
    if (length(cont)==0) stop ("No continuous variables to modify")
    else {
	temp <- (names(m))[-1]      #ignore the response variable
	protected <- protected[-1]
	if (any(protected)| any(factors)) {
###       temp[protected] <- (attr(terms.inner(Terms), 'term.labels'))[protected]
          temp <- attr(unprotect(Terms),'term.labels')
        }
	kname <- temp[keepers[-1]]
      }

    ord <- order(y[,1])
    x <- as.matrix(m[ord, cont, drop=FALSE])
    time <- y[ord,1]
    status <- y[ord,2]
    nvar <- length(cont)

    nline <- 0
    for (i in unique(time[status==1])) nline <- nline + sum(time >=i)
    start <- stop <- event <- double(nline)
    xx <- matrix(double(nline*nvar), ncol=nvar, 
		 dimnames=list(NULL, dimnames(x)[[2]]))
    ltime <- 0
    j<- 1
    keep.index <- NULL

    for (i in unique(time[status==1])) {
	who <- (time >=i)
	nrisk <- sum(who)

	temp <- apply(x[who,,drop=FALSE], 2, rank)
	deaths <- (status[who]==1 & time[who]==i)

	k <- seq(from=j, length.out=nrisk)
	start[k] <- ltime
	stop[k] <-  i
	event[k] <- deaths
	xx[k,] <- transform(temp)
	j <- j + nrisk
	ltime <- i
	keep.index <- c(keep.index, ord[who])
	}

    if (any(keepers)){
	temp <- m[keep.index, keepers, drop=FALSE]
	names(temp) <- kname
	data.frame(start, stop, event, temp, xx, id=keep.index)
        }
    else  data.frame(start, stop, event, xx, id=keep.index)
    }
