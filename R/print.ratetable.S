## $Id: print.ratetable.S 11166 2008-11-24 22:10:34Z therneau $

print.ratetable <- function(x, ...) {
  cat ("Rate table with dimension(s):", attr(x, 'dimid'), "\n")
  attributes(x) <- attributes(x)[c("dim", "dimnames")]
  NextMethod()
}
