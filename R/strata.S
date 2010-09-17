# $Id: strata.S 11166 2008-11-24 22:10:34Z therneau $
# Create a strata variable, possibly from many objects
#
strata <- function(..., na.group=FALSE, shortlabel=FALSE, sep=', ') {
    # First, grab a copy of the call, which will be used to manufacture
    #  labels for unlabeled arguments
    # Then get the arguments as a list
    words <- as.character((match.call())[-1])
    allf <- list(...)
    # If there is only one argument, and it itself is a list, use
    #  it instead
    if(length(allf) == 1 && is.list(ttt <- unclass(allf[[1]]))) allf <- ttt
    nterms <- length(allf)

    # Keep the names of named args as their label, what was typed otherwise
    if (is.null(names(allf))) argname <- words[1:nterms]
    else argname <- ifelse(names(allf) == '', words[1:nterms], names(allf))

    # Process the first argument
    what <- allf[[1]]
    if(is.null(levels(what)))
	    what <- factor(what)
    levs <- unclass(what) - 1
    wlab <- levels(what)
    if (na.group && any(is.na(what))){
	# add "NA" as a level
	levs[is.na(levs)] <- length(wlab)
	wlab <- c(wlab, "NA")
	}

    if (shortlabel) labs <- wlab
    else            labs <- paste(argname[1], wlab, sep='=')

    # Now march through the other variables, if any
    for (i in (1:nterms)[-1]) {
	what <- allf[[i]]
	if(is.null(levels(what)))
		what <- factor(what)
	wlab <- levels(what)
	wlev <- unclass(what) - 1
	if (na.group && any(is.na(wlev))){
	    wlev[is.na(wlev)] <- length(wlab)
	    wlab <- c(wlab, "NA")
	    }
	if (!shortlabel) wlab <- format(paste(argname[i], wlab, sep='='))
	levs <- wlev + levs*(length(wlab))
	labs <- paste(rep(labs, rep(length(wlab), length(labs))),
		      rep(wlab, length(labs)), sep=sep)
	}
    levs <- levs + 1
    ulevs <- sort(unique(levs[!is.na(levs)]))
    levs <- match(levs, ulevs)
    labs <- labs[ulevs]

    factor(levs, labels=labs)
    }
