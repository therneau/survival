#  $Id: coxpenal.df.S 11166 2008-11-24 22:10:34Z therneau $
#
# degrees of freedom computation, based on Bob Gray's paper
#
#  hmat = right hand slice of cholesky of H
#  hinv = right hand slice of cholesky of H-inverse
#  fdiag= diagonal of D-inverse
#  assign.list: terms information
#  ptype= 1 or 3 if a sparse term exists, 2 or 3 if a non-sparse exists
#  nvar = # of non-sparse terms
#  pen1 = the penalty matrix (diagonal thereof) for the sparse terms
#  pen2 = the penalty matrix for the non-sparse terms
#  sparse = indicates which term is the sparse one
coxpenal.df <- function(hmat, hinv, fdiag, assign.list, ptype, nvar,
			pen1, pen2, sparse) {

    if (ptype ==1 & nvar==0) {  #only sparse terms
	hdiag <- 1/fdiag
	list(fvar2=(hdiag-pen1)*fdiag^2, df=sum((hdiag-pen1)*fdiag),
	     fvar = fdiag, trH=sum(fdiag))
        }
    
    else if (ptype==2) {  # only dense ones
	hmat.full <- t(hmat) %*% (ifelse(fdiag==0, 0,1/fdiag)* hmat)
	hinv.full <- hinv %*% (fdiag* t(hinv))
	if (length(pen2)==length(hmat.full)) imat <- hmat.full - pen2
	else                                 imat <- hmat.full - diag(pen2)

	var <- hinv.full %*% imat %*% hinv.full

	if (length(assign.list)==1)
		list(var2=var, df=sum(imat * hinv.full), 
		              trH=sum(diag(hinv.full)), var=hinv.full)
	else {
	    df <- trH <- NULL
	    d2 <- diag(hinv.full)
	    for (i in assign.list) {
		temp <- coxph.wtest(hinv.full[i,i], var[i,i])$solve
		if (is.matrix(temp)) df <- c(df, sum(diag(temp)))
		else                 df <- c(df, sum(temp))
		trH<- c(trH, sum(d2[i]))
	        }
	    list(var2=var, df=df, trH=trH, var = hinv.full)
	    }
        }

    else {
	# sparse terms + other vars
	nf <- length(fdiag) - nvar
       	nr1 <- 1:nf
	nr2 <- (nf+1):(nf+nvar)

	d1 <- fdiag[nr1]
	d2 <- fdiag[nr2]
	temp <- t(hinv[nr1,])
	temp2<- t(hinv[nr2,,drop=FALSE])
	A.diag <- d1 + c(rep(1,nvar) %*% (temp^2*d2))
	B  <- hinv[nr1,] %*% (d2 * temp2)
	C  <- hinv[nr2,] %*% (d2 * temp2)  #see notation in paper
	var2 <- C - t(B) %*% (pen1 * B)

	if (ptype==3) {
	    #additional work when we have penalties on both the sparse term
	    #  and on  non-sparse terms
	    hmat.22 <- t(hmat) %*%(ifelse(fdiag==0, 0,1/fdiag)* hmat)
	    temp <- C - coxph.wtest(hmat.22, diag(nvar))$solve
	    if (nvar==1) {
		var2 <- var2 - C*pen2*C    # C will be 1 by 1
		temp2 <- c(temp*pen2)
		}
	    else if (length(pen2) == nvar) {
		var2 <- var2 - C %*% (pen2 * C)  #diagonal penalty
		temp2 <- sum(diag(temp) * pen2)
		}
	    else {
		var2 <- var2 - C %*% matrix(pen2,nvar) %*% C
		temp2 <- sum(diag(temp * pen2))
		}
	    }
	else temp2 <- 0  #temp2 contains trace[B'A^{-1}B P2], this line: P2=0

	df <- trH <- NULL
	cdiag <- diag(C)

	for (i in 1:length(assign.list)) {
	    if (sparse==i){
		df <- c(df, nf - (sum(A.diag * pen1) + temp2))
		trH <- c(trH, sum(A.diag))
	        }
	    else {
		j <- assign.list[[i]] 
		temp <- coxph.wtest(C[j,j], var2[j,j])$solve
		if (is.matrix(temp)) df <- c(df, sum(diag(temp)))
		else                 df <- c(df, sum(temp))
		trH <- c(trH, sum(cdiag[j]))
	        }
	    }
	list(var=C, df=df, trH=trH, fvar=A.diag, var2=var2)
	}
    }

