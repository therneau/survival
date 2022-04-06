\name{nsk}
\alias{nsk}
\title{
  Natural splines with knot heights as the basis.
}
\description{
Create the design matrix for a natural spline, such that the coefficient
of the resulting fit are the values of the function at the knots.
}
\usage{
nsk(x, df = NULL, knots = NULL, intercept = FALSE, b = 0.05, 
    Boundary.knots = quantile(x, c(b, 1 - b), na.rm = TRUE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the predictor variable.  Missing values are allowed.
}
  \item{df}{
    degrees of freedom. One can supply df rather than knots; ns() then
    chooses df - 1 - intercept knots at suitably chosen quantiles of x
    (which will ignore missing values).
    The default, df = NULL, sets the number of inner knots as length(knots).
}
  \item{knots}{
    breakpoints that define the spline. The default is no knots;
    together with the natural boundary conditions this results in a
    basis for linear regression on x.
    Typical values are the mean or median for one knot, quantiles for
    more 
    knots. See also Boundary.knots.
}
  \item{intercept}{
    if TRUE, an intercept is included in the basis; default is FALSE
}
  \item{b}{default placement of the boundary knots.  A value of
    \code{bs=0} will replicate the default behavior of \code{ns}.
}
\item{Boundary.knots}{
  boundary points at which to impose the natural boundary conditions and
  anchor the B-spline basis.  Beyond these points the function is
  assumed to be linear.
  If both knots and Boundary.knots are supplied, the basis parameters do
  not depend on x. Data can extend beyond Boundary.knots
}
}
\details{
  The \code{nsk} function behaves identically to the \code{ns} function,
  with two exceptions.  The primary one is that the returned basis is
  such that coefficients correspond to the value of the fitted function
  at the knot points.  If \code{intercept = FALSE}, there will be k-1
  coefficients corresponding to the k knots, and they will be the difference
  in predicted value between knots 2-k and knot 1.
  The primary advantage to the basis is that the coefficients are
  directly interpretable.  A second is that tests for the linear and
  non-linear components are simple contrasts.

  The second differnce with \code{ns} is one of opinion with respect to
  the default position for the boundary knots.  The default here is
  closer to that found in the \code{rms::rcs} function.

  This function is a trial if a new idea, it's future inclusion in the
  package is not yet guarranteed.
}
\value{
  A matrix of dimension length(x) * df where either df was supplied or,
  if knots were supplied, df = length(knots) + 1 + intercept.
  Attributes are returned that correspond to the arguments to kns,
  and explicitly give the knots, Boundary.knots etc for use by predict.kns().
}
\note{
  A thin flexible metal or wooden strip is called a spline, and is the
  traditional method for laying out a smooth curve, e.g., for a ship's
  hull or an airplane wing.  Pins are put into a board and the strip is
  passed through them, each pin is a 'knot'.

  A mathematical spline is a piecewise function between each knot.  A
  linear spline will be a set of connected line segments, a quadratic
  spline is a set of connected local quadratic functions, constrained to
  have a continuous first derivative, a cubic spline is cubic between
  each knot, constrained to have continuous first and second
  derivatives, and etc.  Mathematical splines are not an exact
  representation of natural splines: being a physical object the wood or
  metal strip will have continuous derivatives of all orders.  Cubic
  splines are commonly used because they are sufficiently smooth to
  look natural to the human eye.

  If the mathematical spline is further constrained to be linear beyond
  the end knots, this is often called a 'natural spline', due to the
  fact that a wooden or metal spline will also be linear beyond the last
  knots.  Another name for the same object is a 'restricted cubic
  spline', since it is achieved in code by adding further constraints.
  Given a vector of data points and a set of knots, it is possible to
  create a basis matrix X with one column per knot, such that ordinary
  regression of X on y will fit the cubic spline function, hence these
  are also called 'regression splines'. (One of these three labels is no
  better or worse than another, in our opinion).

  Given a basis matrix X with k columns, the
  matrix Z= XT for any k by k nonsingular matrix T is is also a basis
  matrix, and will result in identical predicted values, but a new set
  of coefficients gamma = (T-inverse) beta in place of beta.  One can
  choose the basis functions so that X is easy to construct, to
  make the regression numerically stable, to make tests easier, or based
  on other considerations.
  It seems as though every spline library returns a different basis
  set, which unfortunately makes fits difficult to compare between packages.
  This is yet one more basis set, chosen to make the coefficients
  more interpretable. 
}

\seealso{
\code{\link[splines]{ns}}
}
\examples{
# make some dummy data
tdata <- data.frame(x= lung$age, y = 10*log(lung$age-35) + rnorm(228, 0, 2))
fit1 <- lm(y ~ -1 + nsk(x, df=4, intercept=TRUE) , data=tdata)
fit2 <- lm(y ~ nsk(x, df=3), data=tdata)

# the knots (same for both fits)
knots <- unlist(attributes(fit1$model[[2]])[c('Boundary.knots', 'knots')])
sort(unname(knots))
unname(coef(fit1))  # predictions at the knot points

unname(coef(fit1)[-1] - coef(fit1)[1])  # differences: yhat[2:4] - yhat[1]
unname(coef(fit2))[-1]                  # ditto

\dontrun{
plot(y ~ x, data=tdata)
points(sort(knots), coef(fit1), col=2, pch=19)
coef(fit)[1] + c(0, coef(fit)[-1])
}
}
\keyword{ smooth }
