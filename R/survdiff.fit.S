survdiff.fit <- function(y, x, strat, rho=0) {
    #
    # This routine is almost always called from survdiff
    #  If called directly, remember that it does no error checking
    #
    n <- length(x)
    if (ncol(y) !=2) stop ("Invalid y matrix")
    if (nrow(y) !=n | length(x) !=n) stop("Data length mismatch")

    ngroup <- length(unique(x))
    if (ngroup <2) stop ("There is only 1 group")
    if (is.category(x)) x <- as.numeric(x)
    else x <- match(x, unique(x))

    if (missing(strat)) strat <- rep(1,n)
    else strat <- as.numeric(as.factor(strat))
    nstrat <- length(unique(strat))
    if (length(strat) !=n) stop("Data length mismatch")

    ord <- order(strat, y[,1], -y[,2])
    strat2 <- c(1*(diff(strat[ord])!=0), 1)

    xx <- .C(Csurvdiff2, as.integer(n),
		   as.integer(ngroup),
		   as.integer(nstrat),
		   as.double(rho),
		   as.double(y[ord,1]),
		   as.integer(y[ord,2]),
		   as.integer(x[ord]),
		   as.integer(strat2),
		   observed = double(ngroup*nstrat),
		   expected = double(ngroup*nstrat),
		   var.e    = double(ngroup * ngroup),
		   double(ngroup), double(n))

    if (nstrat==1)  list(expected = xx$expected,
			 observed = xx$observed,
			 var      = matrix(xx$var.e, ngroup, ngroup))
    else            list(expected = matrix(xx$expected, ngroup),
			 observed = matrix(xx$observed, ngroup),
			 var      = matrix(xx$var.e, ngroup, ngroup))
    }
