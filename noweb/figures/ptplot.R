#
# Plot the parse tree for an object
#
ptplot <- function(x, depth=0, ...) {
    if (class(x) == 'call' || class(x)=='formula') {
        temp <- lapply(x[-1], ptplot, depth=depth+1)
        ypos <- lapply(temp, function(x) x$pos)
        offset <- c(0, cumsum(unlist(lapply(ypos, max))))
        if (length(ypos) > 1) {
            for (i in 2:length(ypos)) 
                ypos[[i]] <- ypos[[i]] + offset[i]
            }
            
        mypos = mean(unlist(lapply(ypos, function(x) x[1])))
        rlist <- list(pos= c(mypos, unlist(ypos)),
             depth=c(depth, unlist(lapply(temp, function(x) x$depth))),
             string=c(paste(class(x), deparse(x[[1]]), sep=': '), 
                      unlist(lapply(temp, function(x) x$string))),
             connect.n = c(length(temp),
                           unlist(lapply(temp, function(x) x$connect.n))),
             connect.y = c(unlist(lapply(ypos, function(x) x[1])) - mypos,
                           unlist(lapply(temp, function(x) x$connect.y))))
        }

    else if (class(x) == '(') {
        temp <- ptplot(x[[2]], depth+1)
        rlist <- list(pos=c(temp$pos[1], temp$pos),
                      depth= c(depth, temp$depth),
                      string=c( '(: (', temp$string),
                      connect.n = c(1, temp$connect.n),
                      connect.y = c(0, temp$connect.y))
        }

    else if (is.recursive(x)) {
        rlist <- list(pos=1,
                      depth= depth,
                      string= paste("List:", names(x), collapse=' '),
                      connect.n=0)
        }
    else rlist <- list(pos=1,
                       depth= depth,
                       string= ifelse(is.name(x)|| is.numeric(x), as.character(x),
                                      paste(class(x), x, sep=': ')),
                       connect.n=0)

    if (depth>0) return(rlist)  # recur the function
    else {
        #
        # plot the results
        #
        frame()
        pdepth <- -1 * rlist$depth  # plot larger depths lower on the graph
        par(usr=c(range(rlist$pos), range(pdepth))+ c(-.8,.8,-.5, .5))
        text(rlist$pos, pdepth, rlist$string, ...)
        
        j <- 0
        for (i in 1:length(rlist$pos)) {
            k <- rlist$connect.n[i]
            if (k>0) {
                segments(rep(rlist$pos[i],k), 
                         rep(pdepth[i], k) -.2,
                         rlist$connect.y[j+ 1:k] + rlist$pos[i],
                         rep(pdepth[i], k) -.8, ...)
               j <- j+k
                }
            }
        invisible(rlist)
        }
    }


