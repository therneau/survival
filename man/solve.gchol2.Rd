\name{solve.gchol2}
\alias{solve.gchol2}
\title{
  Solve a matrix equation using the generalized Cholesky decompostion
}
\description{
This function solves the equation Ax=b for x, given b and the
generalized Cholesky decompostion of A.
If only the first argument is given, then a G-inverse of A is
returned.
}
\usage{
\method{solve}{gchol2}(a, b, full=TRUE, ...)
}
\arguments{
\item{a}{
a generalized cholesky decompostion of a matrix, as
returned by the \code{gchol} function.
}
\item{b}{
a numeric vector or matrix, that forms the right-hand side of the equation.
}
\item{full}{
solve the problem for the full (orignal) matrix, or for the cholesky
matrix.
}
\item{...}{other arguments are ignored}
}
\value{
if argument \code{b} is not present, the inverse of
\code{a} is returned, otherwise the solution to 
matrix equation.
}
\details{
A symmetric matrix A can be decomposed as LDL', where L is a lower
triangular matrix with 1's on the diagonal, L' is the transpose of
L, and D is diagonal.  
This routine solves either the original problem Ay=b 
(\code{full} argument) or the subproblem sqrt(D)L'y=b.
If \code{b} is missing it returns the inverse of
A or L, respectively.
}
\seealso{
gchol2
}
\examples{
# Create a matrix that is symmetric, but not positive definite
#   The matrix temp has column 6 redundant with cols 1-5
smat <- matrix(1:64, ncol=8)
smat <- smat + t(smat) + diag(rep(20,8))  #smat is 8 by 8 symmetric
temp <-  smat[c(1:5, 5:8), c(1:5, 5:8)]
ch1  <- gchol(temp)

ginv <- solve(ch1, full=FALSE)  # generalized inverse of ch1
tinv <- solve(ch1, full=TRUE)   # generalized inverse of temp
all.equal(temp \%*\% tinv \%*\% temp, temp)
}
\keyword{array}
