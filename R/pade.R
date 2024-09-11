# This is a trimmed copy of expm:::expm2.R, optimized for transition matrices
# The row sums of A are 0, with the diagonal the only negative element, and
#    transitions away from the diagonal are normally small.  
#  I don't need balancing, don't need Matrix objects
#  The matrices are usually very nice.
#  If I need derivatives, then deriv is an (n,n, k) array, the first
#  slice contains dA/(d beta1), elementwise,  and etc for the second, third etc.
#
pade <- function(A, deriv) {
    n <- nrow(A)
    nA <- max(colSums(abs(A)))
    I <- diag(n)
    if (!missing(deriv)) {
        if (!is.array(deriv) || length(dim(deriv))!=3 ||
            any(dim(deriv)[1:2] != n)) stop("invalid derivative array")
        nderiv <- dim(deriv)[3]
        dmat <- array(0., dim=c(n, n, nderiv))
    }
    else nderiv <- 0

    ## If the norm is small enough, use the Pade Approximation (PA) directly
    if (nA <= 2.1) {
    	t <- c(0.015, 0.25, 0.95, 2.1)
    	## the minimal m for the PA :
    	l <- which.max(nA <= t)
        nterm <- l+1  # for reporting: number of terms used

    	## Calculate PA
    	C <- rbind(c(120,60,12,1,0,0,0,0,0,0),
    		   c(30240,15120,3360,420,30,1,0,0,0,0),
    		   c(17297280,8648640,1995840,277200,25200,1512,56,1,0,0),
    		   c(17643225600,8821612800,2075673600,302702400,30270240,
    		     2162160,110880,3960,90,1))
        A2 <- A %*% A
    	P <- I
    	U <- C[l,2]*I
    	V <- C[l,1]*I

        for (k in 1:l) {
    	    P <- P %*% A2
    	    U <- U + C[l,(2*k)+2]*P
    	    V <- V + C[l,(2*k)+1]*P
    	}
    	U <- A %*% U
    	X <- solve(V-U,V+U)

        if (nderiv > 0) {
            for (dk in 1:nderiv) {
                P <- I
                U <- C[l,2]*I
                V <- C[l,1]*I
                dA <- deriv[,,dk]
                dA2 <- A %*% dA + dA %*% A
                dU <- dV <- dP <-  0*I
                for (k in 1:l) {
                    dP <- dP %*% A2 + P %*% dA2
                    P <- P %*% A2
           	    U <- U + C[l,(2*k)+2]*P
                    dU <- dU + C[l,(2*k)+2]* dP
                    V <- V + C[l,(2*k)+1]*P
                    dV <- dV + C[l,(2*k)+1]*dP
                }
                dU <- dA %*%U + A %*% dU
                U <- A %*% U
                dmat[,,dk] <- -solve(V-U, dV-dU)%*% X + solve(V-U, dV + dU)
            }
        }
    }
    else {
        # This section should be rarer than rare
        s <- log2(nA/5.4)
        B <- A
        ## Scaling
        if (s > 0) {
            s <- ceiling(s)
            B <- B/(2^s)
            if (nderiv>0) deriv <- deriv/(2^s)
        }
        nterm <- max(0, ceiling(s)) + 6

	## Calculate PA
	c. <- c(64764752532480000,32382376266240000,7771770303897600,
                1187353796428800, 129060195264000,10559470521600, 670442572800,
                33522128640, 1323241920, 40840800,960960,16380, 182,1)
	B2 <- B %*% B
	B4 <- B2 %*% B2
	B6 <- B2 %*% B4

	U <- B %*% (B6 %*% (c.[14]*B6 + c.[12]*B4 + c.[10]*B2) +
		    c.[8]*B6 + c.[6]*B4 + c.[4]*B2 + c.[2]*I)
	V <- B6 %*% (c.[13]*B6 + c.[11]*B4 + c.[9]*B2) +
	    c.[7]*B6 + c.[5]*B4 + c.[3]*B2 + c.[1]*I

	X <- solve(V-U,V+U)

        if (nderiv > 0) {
            for (dk in 1:nderiv) {
                dB <- deriv[,,dk]
                dB2 <- B %*% dB + dB %*% B
                dB4 <- B2%*% dB2 + dB2 %*% B2
                dB6 <- B2 %*%dB4 + dB2 %*% B4
                dU <- dB %*% (B6 %*% (c.[14]*B6 + c.[12]*B4 + c.[10]*B2) +
                            c.[8]*B6 + c.[6]*B4 + c.[4]*B2 + c.[2]*I) +
                       B %*% (dB6 %*% (c.[14]*B6 + c.[12]*B4 + c.[10]*B2) +
                               B6 %*% (c.[14]*dB6 + c.[12]*dB4 + c.[10]*dB2) +
                              c.[8]*dB6 + c.[6]*dB4 + c.[4]*dB2)
                dV <-  dB6 %*% (c.[13]*B6 + c.[11]*B4 + c.[9]*B2) +
                        B6 %*% (c.[13]*dB6 + c.[11]*dB4 + c.[9]*dB2) +
                            c.[7]*dB6 + c.[5]*dB4 + c.[3]*dB2

                dX <- -solve(V-U, dV-dU)%*% X + solve(V-U, dV + dU)
            }
            if (s > 0) for (t in 1:s) {
                dX <- X %*% dX + dX %*% X
                X <- X %*% X
            }
            dmat[,,dk] <- dX
        }
        else if (s>0) for (t in 1:s) {
                 X <- X %*% X
        }
        
    }
    if (nderiv > 0) list(P=X, dmat=dmat)
    else X
}

