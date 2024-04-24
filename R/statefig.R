# Automatically generated from the noweb directory
statefig <- function(layout, connect, margin=.03, box=TRUE,
                     cex=1, col=1, lwd=1, lty=1, bcol= col,
                     acol=col, alwd = lwd, alty= lty, offset=0) {
    # set up an empty canvas
    frame();  # new environment
    par(usr=c(0,1,0,1))
    if (!is.numeric(layout))
        stop("layout must be a numeric vector or matrix")
    if (!is.matrix(connect) || nrow(connect) != ncol(connect))
        stop("connect must be a square matrix")
    nstate <- nrow(connect)
    dd <- dimnames(connect)
    if (!is.null(dd[[1]])) statenames <- dd[[1]]
    else if (is.null(dd[[2]])) 
        stop("connect must have the state names as dimnames")
    else statenames <- dd[[2]]

    # expand out all of the graphical parameters.  This lets users
    #  use a vector of colors, line types, etc
    narrow <- sum(connect!=0) 
    acol <- rep(acol, length.out=narrow)
    alwd <- rep(alwd, length.out=narrow)
    alty <- rep(alty, length.out=narrow)

    bcol <- rep(bcol, lengt.outh=nstate)
    lty  <- rep(lty, length.out=nstate)
    lwd  <- rep(lwd, length.out=nstate)
    
    col <- rep(col, length.out=nstate)  # text colors
 
    if (is.matrix(layout) && ncol(layout)==2 && nrow(layout) > 1) {
        # the user provided their own
        if (any(layout <0) || any(layout >1))
            stop("layout coordinates must be between 0 and 1")
        if (nrow(layout) != nstate)
            stop("layout matrix should have one row per state")
        cbox <- layout
    }
    else {
        if (any(layout <=0 | layout != floor(layout)))
            stop("non-integer number of states in layout argument")
        space <- function(n) (1:n -.5)/n   # centers of the boxes
        if (sum(layout) != nstate) stop("number of boxes != number of states")
        cbox <- matrix(0, ncol=2, nrow=nstate)  #coordinates will be here
        n <- length(layout)
     
        ix <- rep(seq_along(layout), layout) 
        if (is.vector(layout) || ncol(layout)> 1) { #left to right     
            cbox[,1] <- space(n)[ix]
            for (i in 1:n) cbox[ix==i,2] <- 1 -space(layout[i])
        } else { # top to bottom
            cbox[,2] <- 1- space(n)[ix]
            for (i in 1:n) cbox[ix==i,1] <- space(layout[i])
        }
    }
    text(cbox[,1], cbox[,2], statenames, cex=cex, col=col)  # write the labels
    textwd <- strwidth(statenames, cex=cex)
    textht <- strheight(statenames, cex=cex)
    temp <- par("pin")   #plot region in inches
    dx <- margin * temp[2]/mean(temp)  # extra to add in the x dimension
    dy <- margin * temp[1]/mean(temp)  # extra to add in y

    if (box) {
        drawbox <- function(x, y, dx, dy, lwd, lty, col) {
            lines(x+ c(-dx, dx, dx, -dx, -dx),
                  y+ c(-dy, -dy, dy, dy, -dy), lwd=lwd, lty=lty, col=col)
        }
        for (i in 1:nstate) 
            drawbox(cbox[i,1], cbox[i,2], textwd[i]/2 + dx, textht[i]/2 + dy,
                    col=bcol[i], lwd=lwd[i], lty=lty[i])
        dx <- 2*dx; dy <- 2*dy   # move arrows out from the box
        }
    arrow2 <- function(...) arrows(..., angle=20, length=.1)
    doline <- function(x1, x2, d, delta1, delta2, lwd, lty, col) {
        if (d==0 && x1[1] ==x2[1]) { # vertical line
            if (x1[2] > x2[2]) # downhill
                arrow2(x1[1], x1[2]- delta1[2], x2[1], x2[2] + delta2[2],
                       lwd=lwd, lty=lty, col=col)
            else arrow2(x1[1], x1[2]+ delta1[2], x2[1], x2[2] - delta2[2],
                        lwd=lwd, lty=lty, col=col)
        }
        else if (d==0 && x1[2] == x2[2]) {  # horizontal line
            if (x1[1] > x2[1])  # right to left
                arrow2(x1[1]-delta1[1], x1[2], x2[1] + delta2[1], x2[2],
                       lwd=lwd, lty=lty, col=col)
            else arrow2(x1[1]+delta1[1], x1[2], x2[1] - delta2[1], x2[2],
                        lwd=lwd, lty=lty, col=col)
        }
        else {
            temp <- phi(x1[1], x1[2], x2[1], x2[2], d, delta1, delta2)
            if (d==0) {        
                arrow2(temp$center[1] + temp$r*cos(temp$angle[1]),
                       temp$center[2] + temp$r*sin(temp$angle[1]),
                       temp$center[1] + temp$r*cos(temp$angle[2]),
                       temp$center[2] + temp$r*sin(temp$angle[2]),
                       lwd=lwd, lty=lty, col=col)
            }
            else {
                # approx the curve with 21 segments
                #  arrowhead on the last one
                phi <- seq(temp$angle[1], temp$angle[2], length.out=21)
                lines(temp$center[1] + temp$r*cos(phi),
                      temp$center[2] + temp$r*sin(phi), lwd=lwd, lty=lty, col=col)
                arrow2(temp$center[1] + temp$r*cos(phi[20]),
                       temp$center[2] + temp$r*sin(phi[20]),
                       temp$center[1] + temp$r*cos(phi[21]),
                       temp$center[2] + temp$r*sin(phi[21]),
                       lwd=lwd, lty=lty, col=col)
            }
        }
    }
    k <- 1
    for (j in 1:nstate) {
        for (i in 1:nstate) {
            if (i != j && connect[i,j] !=0) {
                if (connect[i,j] == 2-connect[j,i] && offset!=0) {
                    #add an offset
                    toff <- c(cbox[j,2] - cbox[i,2], cbox[i,1] - cbox[j,1])
                    toff <- -offset *toff/sqrt(sum(toff^2))
                    doline(cbox[i,]+toff, cbox[j,]+toff, connect[i,j]-1,
                           delta1 = c(textwd[i]/2 + dx, textht[i]/2 + dy),
                           delta2 = c(textwd[j]/2 + dx, textht[j]/2 + dy),
                           lty=alty[k], lwd=alwd[k], col=acol[k])
                    }
                else doline(cbox[i,], cbox[j,], connect[i,j]-1,
                            delta1 = c(textwd[i]/2 + dx, textht[i]/2 + dy),
                            delta2 = c(textwd[j]/2 + dx, textht[j]/2 + dy),
                            lty=alty[k], lwd=alwd[k], col=acol[k])
                k <- k +1
            }
        }
    }
    
    dimnames(cbox) <- list(statenames, c("x", "y"))
    invisible(cbox)
}
statefigx <- function(x, C, r, a1, a2) {
    temp <-(x - C[1])/r
    if (abs(temp) >1) return(NULL)  # no intersection of the arc and x
    phi <- acos(temp)  # this will be from 0 to pi
    pi <- 3.1415926545898   # in case someone has a variable "pi" 
    if (x > C[1]) phi <-  c(phi, pi - phi)
    else          phi <- -c(phi, pi - phi)
    # Add reflection about the X axis, in both forms
    phi <- c(phi, -phi, 2*pi - phi) 
    amax <- max(a1, a2)
    amin <- min(a1, a2)
    phi[phi<amax & phi > amin]
}
statefigy <-  function(y, C, r, a1, a2) {
    pi <- 3.1415926545898   # in case someone has a variable named "pi" 
    amax <- max(a1, a2)
    amin <- min(a1, a2)
    temp <-(y - C[2])/r
    if (abs(temp) >1) return(NULL)  # no intersection of the arc and y
    phi <- asin(temp)  # will be from -pi/2 to pi/2
    phi <- c(phi, sign(phi)*pi -phi)  # reflect about the vertical
    phi <- c(phi, phi + 2*pi)
    phi[phi<amax & phi > amin]
}
phi <- function(x1, y1, x2, y2, d, delta1, delta2) {
    # d = height above the line
    theta <- atan2(y2-y1, x2-x1)    # angle from center to center
    if (abs(d) < .001) d=.001       # a really small arc looks like a line

    z <- sqrt((x2-x1)^2 + (y2 - y1)^2) /2 # half length of chord
    ab <- c((x1 + x2)/2, (y1 + y2)/2)      # center of chord
    r  <- abs(z*(1 + d^2)/ (2*d))
    if (d >0) C  <- ab + (r - d*z)* c(-sin(theta), cos(theta)) # center of arc
    else      C  <- ab + (r + d*z)* c( sin(theta), -cos(theta))

    a1 <- atan2(y1-C[2], x1-C[1])   # starting angle
    a2 <- atan2(y2-C[2], x2-C[1])   # ending angle
    if (abs(a2-a1) > pi) {
        # a1= 3 and a2=-3, we don't want to include 0
        # nor for a1=-3 and a2=3
        if (a1>0) a2 <- a2 + 2 *pi 
        else a1 <- a1 + 2*pi
    }
    if (d > 0) { #counterclockwise
        phi1 <- min(statefigx(x1 + delta1[1], C, r, a1, a2),
                    statefigx(x1 - delta1[1], C, r, a1, a2),
                    statefigy(y1 + delta1[2], C, r, a1, a2),
                    statefigy(y1 - delta1[2], C, r, a1, a2), na.rm=TRUE)
        phi2 <- max(statefigx(x2 + delta2[1], C, r, a1, a2),
                    statefigx(x2 - delta2[1], C, r, a1, a2),
                    statefigy(y2 + delta2[2], C, r, a1, a2),
                    statefigy(y2 - delta2[2], C, r, a1, a2), na.rm=TRUE)
    }
    else { # clockwise
        phi1 <- max(statefigx(x1 + delta1[1], C, r, a1, a2),
                    statefigx(x1 - delta1[1], C, r, a1, a2),
                    statefigy(y1 + delta1[2], C, r, a1, a2),
                    statefigy(y1 - delta1[2], C, r, a1, a2), na.rm=TRUE)
        phi2 <- min(statefigx(x2 + delta2[1], C, r, a1, a2),
                    statefigx(x2 - delta2[1], C, r, a1, a2),
                    statefigy(y2 + delta2[2], C, r, a1, a2),
                    statefigy(y2 - delta2[2], C, r, a1, a2), na.rm=TRUE)
    }

    list(center=C, angle=c(phi1, phi2), r=r)
}
