# Natural spline, parameterized so that the coefficients are the values of the
#  spline at the knots.
#
nsk <- function(x, df=NULL, knots=NULL, intercept=FALSE, b=.05,
                Boundary.knots = quantile(x, c(b, 1-b), na.rm=TRUE)) {

    basis <- ns(x, df=df, knots=knots, intercept=intercept,
                Boundary.knots = Boundary.knots)
    knots <- attr(basis, "knots")

    kx <- sort(c(knots, Boundary.knots))
    kbasis <- ns(kx, df=df, knots=knots, intercept=intercept,
                Boundary.knots = Boundary.knots)   
        
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



