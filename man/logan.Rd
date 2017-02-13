\name{logan}
\docType{data}
\alias{logan}
\title{Data from the 1972-78 GSS data used by Logan}
\usage{data(logan)}
\description{
Intergenerational occupational mobility data with covariates.
}
\format{
  A data frame with 838 observations on the following 4 variables.
  \describe{
    \item{occupation}{subject's occupation, a factor with levels
      \code{farm}, \code{operatives}, \code{craftsmen}, \code{sales},
       and \code{professional}}
    \item{focc}{father's occupation}
    \item{education}{total years of schooling, 0 to 20}
    \item{race}{levels of \code{non-black} and \code{black}}
  }
}
\source{
  General Social Survey data, see the web site for detailed information
  on the variables.
  \url{http://www3.norc.org/GSS+Website}.
}
\references{
  Logan, John A. (1983). A Multivariate Model for Mobility Tables.
  \cite{American Journal of Sociology} 89: 324-349.}
\keyword{datasets}
