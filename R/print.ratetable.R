
print.ratetable <- function(x, ...) {
    if (is.null(attr(x, 'dimid')))
        cat ("Rate table with dimension(s):", names(dimnames(x)), "\n")
    else  cat ("Rate table with dimension(s):", attr(x, 'dimid'), "\n")
  attributes(x) <- attributes(x)[c("dim", "dimnames")]
  NextMethod()
}
