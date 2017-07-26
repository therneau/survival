#
# Code for the generalized cholesky  A = LDL', where L is lower triangular
#   with 1's on the diagonal, and D is diagonal.
# The decompostions exists for any square symmetric matrix.  
# If A is positive definite, then all elements of D will be positve.
# If A is not full rank, then 0's on the diagonal of D signal the redundant 
#   columns.  
#
#  Gchol already exists as an exported class in bdsmatrix, and I can't change 
#   that routine since other libraries depend on it.  I can't redefine
#   gchol here since coxme needs both survival and bdsmatrix, which will lead
#   to a conflict.  However, I really only need the routine for the 
#   Yates function, so set it up as a gchol2 class.
#
as.matrix.gchol2 <- function(x, ...) {
    diag(x) <- 1
    class(x) <- "matrix"
    }

gchol2 <- function(x,  tolerance = sqrt(.Machine$double.eps)) {
    d <- dim(x)
    if (d[1] != d[2]) 
        stop("Cholesky decomposition requires a square matrix")
    storage.mode(x) <- "double"  # failsafe
    temp <- .C("gchol", as.integer(d[1]),
               x =   as.double(x),
               rank= as.double(tolerance))
    
    attr(temp, "dimnames") <- dimnames(x)
    attr(temp, dim) <- d
    class(temp) <- "gchol2"
    temp
}   

solve.gchol <- function(a, b, full=TRUE, ...) {
    if (full) flag<-0 else flag<-1

    d <- dim(a)
    if (missing(b)) {
	# Return the inverse of the original matrix, for which a is the chol
	temp <- .C("gchol_inv", as.integer(d), 
                   x=as.double(a),
                   as.integer(flag))$x
	matrix(temp, d[1])
	}

    else {  # solve for right-hand side
	if (length(b) == d[1]) {
	    temp <- .C("gchol_solve", as.integer(d),
		                      as.double(a@.Data),
		                      y=as.double(b),
                                      as.integer(flag))$y
	    temp
	    }
	else {
	    if (!is.matrix(b) || nrow(b) != d[1]) 
		stop("number or rows of b must equal number of columns of a")
	    new <- b
	    for (i in 1:ncol(b)) {
		new[,i] <- .C("gchol_solve", as.integer(d),
			                     as.double(a@.Data),
		                             y=as.double(b[,i]),
                                             as.integer(flag))$y
		}
	    new
	    }
	}
    }
