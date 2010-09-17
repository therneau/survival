# SCCS $Id: aareg.taper.S 11166 2008-11-24 22:10:34Z therneau $
#
# Do running averages of an information matrix
# 
aareg.taper <- function(taper, imat, nevent) {
    dd <- dim(imat)
    if (length(taper)==0 || !is.numeric(taper) || any(taper <=0)) 
        stop("Invalid taper vector")

    ntaper <- length(taper)
    ntime <- dd[3]
    if (ntaper > ntime) {
        taper <- taper[1:ntime]
        ntaper <- ntime
        }

    #
    # Turn imat into an array: 1 row per coef, one col per time
    #  and then scale it by the number of events to get a variance
    # (coxph.detail returns imat = var(X) * nevents)
    #
    imat <- matrix(as.vector(imat), ncol=dd[3])
    imat <- imat / rep(nevent, rep(dd[1]*dd[2], dd[3]))

    if (ntaper >1) {
        smoother <- matrix(0., ntime, ntime)
        tsum <- cumsum(rev(taper))
        for (i in 1:ntaper) 
            smoother[1:i, i] <- taper[seq(to=ntaper, length=i)]/tsum[i]
        if (ntaper < ntime) {
            for (i in (ntaper+1):ntime)
                smoother[seq(to=i, length=ntaper),i] <- taper/tsum[ntaper]
            }
        imat <- imat %*% smoother
        }
    array(imat, dim=dd)
    }

            
