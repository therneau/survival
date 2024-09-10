library(survival)

visits <- function(sfit, debug=FALSE) {
    # Compute the expected number of visits by time t, from a mutlistate
    #  survfit object
    from <- as.numeric(sub("\\:[0-9]*$", "", colnames(sfit$n.transition)))
    to   <- as.numeric(sub("^[0-9]*\\:", "", colnames(sfit$n.transition)))
    nstate <- length(sfit$states)
    ntime  <- length(sfit$time)
    if (is.null(sfit$p0)) p0 <- c(1, rep(0, nstate-1))
    else p0 <- sfit$p0

    cmat <- model.matrix(~ factor(to, levels=1:nstate) -1) # collapse matrix
    # Do the following per strata (rows in the stratum)
    tfun <- function(irow, p0) {
        n <- length(irow)
        temp <- fit$n.transition[irow,]
        lambda <- temp/ifelse(fit$n.risk[irow,from]==0, 1, fit$n.risk[irow,from])
        increment <- rbind(p0[from], fit$pstate[irow[-n], from]) * lambda
        if (debug) browser()
        apply(increment %*% cmat, 2, cumsum)
    }

    if (is.null(sfit$strata)) vmat <- tfun(1:ntime, p0)
    else {
        nstrat <- length(sfit$strata)
        rowlist <- split(1:length(sfit$time), rep(1:nstrat, sfit$strata))
        if (!is.matrix(p0)) p0 <- matrix(rep(p0, each=nstrat, nrow=nstrat))
        vmat <- matrix(0, ntime, nstate)
        for (i in 1:nstrat) {
            vmat[rowlist[[i]],] <- tfun(rowlist[[i]], p0[i,])
        }
    }  
    colnames(vmat) <- sfit$states
    vmat
}

tdata <- myeloid  # temporary working copy
tied <- with(tdata, (!is.na(crtime) & !is.na(txtime) & crtime==txtime))
tdata$crtime[tied] <- tdata$crtime[tied] -1
mdata <- tmerge(tdata[,1:2], tdata,  id=id,  death= event(futime, death),
                sct = event(txtime), cr = event(crtime), 
                relapse = event(rltime),
                priorcr = tdc(crtime), priortx = tdc(txtime))
temp <- with(mdata, cr + 2*sct  + 4*relapse + 8*death)
table(temp)
mdata$event <- factor(temp, c(0,1,2,4,8),
                       c("none", "CR", "SCT", "relapse", "death"))
fit <- survfit(Surv(tstart, tstop, event) ~ trt, mdata, id=id)
test <- visits(fit)
round(test[cumsum(fit$strata),],2)



