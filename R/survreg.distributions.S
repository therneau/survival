# $Id: survreg.distributions.S 11198 2009-02-02 05:01:22Z therneau $
#
# Create the survreg.distributions object
#
survreg.distributions <- list(
'extreme' = list(
    name = "Extreme value",
    variance = function(parm) pi^2/6,
    init  = function(x, weights, ...) {
	mean <- sum(x*weights)/ sum(weights)
	var  <- sum(weights*(x-mean)^2)/ sum(weights)
	c(mean + .572, var/1.64)
	},
    deviance= function(y, scale, parms) {
	status <- y[,ncol(y)]
	width <- ifelse(status==3,(y[,2] - y[,1])/scale, 1)
	temp <- width/(exp(width)-1)
	center <- ifelse(status==3, y[,1] - log(temp), y[,1])
	temp3 <- (-temp) + log(1- exp(-exp(width)))
	best <- ifelse(status==1, -(1+log(scale)),
				    ifelse(status==3, temp3, 0))
	list(center=center, loglik=best) 
	},
    density = function(x,parms) {
	w <- exp(x)
	ww <- exp(-w)
	cbind(1-ww, ww, w*ww, (1-w), w*(w-3) +1)
	},
    quantile = function(p,parms) log(-log(1-p))
    ),

logistic = list(
    name  = "Logistic",
    variance = function(parm) pi^2/3,
    init  = function(x, weights, ...) {
	mean <- sum(x*weights)/ sum(weights)
	var  <- sum(weights*(x-mean)^2)/ sum(weights)
	c(mean, var/3.2)
	},
    deviance= function(y, scale, parms) {
	status <- y[,ncol(y)]
	width <- ifelse(status==3,(y[,2] - y[,1])/scale, 0)
	center <- y[,1] - width/2
	temp2 <- ifelse(status==3, exp(width/2), 2) #avoid a log(0) message
	temp3 <- log((temp2-1)/(temp2+1))
	best <- ifelse(status==1, -log(4*scale),
				    ifelse(status==3, temp3, 0))
	list(center=center, loglik=best) 
	},
    density = function(x, parms) {
	w <- exp(x)
	cbind(w/(1+w), 1/(1+w), w/(1+w)^2, (1-w)/(1+w), (w*(w-4) +1)/(1+w)^2)
	},
    quantile = function(p, parms) log(p/(1-p))
    ),

gaussian = list(
    name  = "Gaussian",
    variance = function(parm) 1,
    init  = function(x, weights, ...) {
	mean <- sum(x*weights)/ sum(weights)
	var  <- sum(weights*(x-mean)^2)/ sum(weights)
	c(mean, var)
	},
    deviance= function(y, scale, parms) {
	status <- y[,ncol(y)]
	width <- ifelse(status==3,(y[,2] - y[,1])/scale, 0)
	center <- y[,1] - width/2
	temp2 <- log(1 - 2*pnorm(width/2))
	best <- ifelse(status==1, -log(sqrt(2*pi)*scale),
				ifelse(status==3, temp2, 0))
	list(center=center, loglik=best) 
	},
    density = function(x, parms) {
	cbind(pnorm(x), pnorm(-x), dnorm(x), -x, x^2-1)	
	},
    quantile = function(p, parms) qnorm(p)
    ),

weibull = list(
    name  = "Weibull",
    dist  = 'extreme',
    trans = function(y) log(y),
    dtrans= function(y) 1/y ,
    itrans= function(x) exp(x)
    ),

exponential = list(
    name  = "Exponential",
    dist  = 'extreme',
    trans = function(y) log(y),
    dtrans= function(y) 1/y,
    scale =1,
    itrans= function(x) exp(x)
    ),

rayleigh = list(
    name  = "Rayleigh",
    dist  = 'extreme',
    trans = function(y) log(y),
    dtrans= function(y) 1/y,
    itrans= function(x) exp(x),
    scale =0.5
    ),

loggaussian = list(
    name  = "Log Normal",
    dist  = 'gaussian',
    trans = function(y) log(y),
    itrans= function(x) exp(x),
    dtrans= function(y) 1/y
    ),

lognormal = list(
    name  = "Log Normal",
    dist  = 'gaussian',
    trans = function(y) log(y),
    itrans= function(x) exp(x),
    dtrans= function(y) 1/y
    ),

loglogistic = list(
    name = "Log logistic",
    dist = 'logistic',
    trans = function(y) log(y),
    dtrans= function(y) 1/y ,
    itrans= function(x) exp(x)
    ),

t = list(
    name  = "Student-t",
    variance = function(df) df/(df-2),
    parms = c(df=4),
    init  = function(x, weights, df) {
	if (df <=2) stop ("Degrees of freedom must be >=3")
	mean <- sum(x*weights)/ sum(weights)
	var  <- sum(weights*(x-mean)^2)/ sum(weights)
	c(mean, var*(df-2)/df)
	},
    deviance= function(y, scale, parms) {
	status <- y[,ncol(y)]
	width <- ifelse(status==3,(y[,2] - y[,1])/scale, 0)
	center <- y[,1] - width/2
	temp2 <- log(1 - 2*pt(width/2, df=parms))
	best <- ifelse(status==1, -log(dt(0, df=parms)*scale),
				ifelse(status==3, temp2, 0))
	list(center=center, loglik=best) 
	},
    density = function(x, df) {
	cbind(pt(x, df), pt(-x, df), dt(x,df),
	      -(df+1)*x/(df+x^2), 
	      (df+1)*(x^2 *(df+3)/(df+x^2) - 1)/(df +x^2))
	},
    quantile = function(p, df) qt(p, df)
  )
)







