\name{gchol2}
\alias{gchol2}
\title{
  Generalized Cholesky decompostion
}
\description{
Perform the generalized Cholesky decompostion of a real symmetric matrix.
}
\usage{
gchol2(x, tolerance=sqrt(.Machine$double.eps))
}
\arguments{
\item{x}{
the symmetric matrix to be factored
}
\item{tolerance}{
the numeric tolerance for detection of singular columns in x.
}
}
\value{
an object of class \code{gchol} containing the
generalized Cholesky decompostion.
It has the appearance of a lower triangular matrix.
}
\details{
A symmetric matrix A can be decomposed as LDL', where L is a lower
triangular matrix with 1's on the diagonal, L' is the transpose of
L, and D is diagonal.
The inverse of L is also lower-triangular, with 1's on the diagonal.  
If all elements of D are positive, then A must be symmetric
positive definite (SPD), and the solution can be reduced the usual
Cholesky decomposition U'U where U is upper triangular and
U = sqrt(D) L'.

The main advantage of the generalized form is that it admits
of matrices that are not of full rank: D will contain zeros 
marking the redundant columns, and the rank of A is the
number of non-zero columns.  If all elements of D are zero or
positive, then A is a non-negative definite (NND) matrix.
The generalized form also has the (quite minor) numerical advantage
of not requiring square roots during its calculation.
To extract the components of the decompostion, use the 
\code{diag} and \code{as.matrix}
functions.

This method is heavily used in the survival code.  One purpose of this
routine is to expose the underlying C routine.
}
\note{
The \code{bdsmatrix} library implements an equivalent \code{gchol}
function, but using the S4 class structure.
}
\seealso{\code{\link{solve.gchol}}
}
\examples{
# Create a matrix that is symmetric, but not positive definite
#   The matrix temp has column 6 redundant with cols 1-5
smat <- matrix(1:64, ncol=8)
smat <- smat + t(smat) + diag(rep(20,8))  #smat is 8 by 8 symmetric
temp <-  smat[c(1:5, 5:8), c(1:5, 5:8)]
ch1  <- gchol2(temp)

print(as.matrix(ch1), digits=4)   # print out L
print(diag(ch1))        # Note the zero at position 6

ginv <- solve(ch1)    # generalized inverse
diag(ginv)            # also has column 6 marked as singular
}
\keyword{array}
