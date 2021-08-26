# Natural spline, parameterized so that the coefficients are the values of the
#  spline at the knots.
#
nsk <- function(x, df=NULL, knots=NULL, intercept=FALSE, b=.05,
                Boundary.knots = quantile(x, c(b, 1-b), na.rm=TRUE)) {

    if ((is.logical(Boundary.knots) && !Boundary.knots[1]) ||
        length(Boundary.knots) == 0)  kx <- knots
    else if (length(knots) == 0) kx <- Boundary.knots
    else {
        if (Boundary.knots[2] <= max(knots)) Boundary.knots <- Boundary.knots[1]
        if (Boundary.knots[1] >= min(knots)) Boundary.knots <- Boundary.knots[-1]
        kx <- sort(c(knots, Boundary.knots))
        }

    if (length(kx) ==2) { # nothing left but a line
        if (intercept) cbind(1.0, x) else as.matrix(x)
    }
    else { 
        j <- c(1, length(kx))
        bknot <- kx[j]
        iknot <- kx[-j]
        if (length(iknot) ==0) 
            basis <- ns(x, df=df, intercept=intercept,
                        Boundary.knots = bknot)
        else basis <- ns(x, df=df, knots= iknot, intercept=intercept,
                       Boundary.knots = bknot)
        iknot <- attr(basis, "knots")
        kx <- c(bknot[1], iknot, bknot[2])
        kbasis <- ns(kx, df=df, knots=iknot, intercept=intercept,
                     Boundary.knots = bknot)         
        # 
        # We know that gamma = Kbasis *beta = yhat at knots
        # then (basis* K-inverse) gamma =  basis * beta
        #  inverting a 3x3 or 4x4 matrix is not compute issue
        #
        if (intercept) ibasis <- basis %*% solve(kbasis)
        else ibasis <- (cbind(1, basis) %*% solve(cbind(1, kbasis)))[, -1]
 
        attributes(ibasis) <- attributes(basis)  # retain the attributes of ns
        class(ibasis) <- c("nsk", class(basis))
        ibasis
    }
}

makepredictcall.nsk <- function(var, call)
{
    if(as.character(call)[1L] == "nsk" || 
       (is.call(call) && identical(eval(call[[1L]]), nsk))) {
	at <- attributes(var)[c("knots", "Boundary.knots", "intercept")]
	call <- call[1L:2L]
	call[names(at)] <- at
    }
    call
}

# the predict method is inherited from ns



