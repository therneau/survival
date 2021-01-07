# bounded links for pseudovalues

blogit <- function(edge=.05) {
    new <- make.link("logit")
    new$linkfun <- function(mu) { 
        x <- (pmax(edge, pmin(mu, 1-edge)))
        log(x/(1-x))
    }
    new$name <- "blogit"
    new
}

bcloglog <- function(edge=.05) {
    new <- make.link("cloglog")
    new$linkfun <- function(mu) { 
        x <- (pmax(edge, pmin(mu, 1-edge)))
        log(-log(1-x))
    }
    new$name <- "bcloglog"
    new
}

bprobit <- function(edge=.05) {
    new <- make.link("probit")
    new$linkfun <- function(mu) { 
        x <- (pmax(edge, pmin(mu, 1-edge)))
        qnorm(x)
    }
    new$name <- "probit"
    new
}
