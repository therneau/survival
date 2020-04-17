# Automatically generated from the noweb directory
survfit.coxph <-
  function(formula, newdata, se.fit=TRUE, conf.int=.95, individual=FALSE,
            stype=2, ctype, 
            conf.type=c("log", "log-log", "plain", "none", "logit", "arcsin"),
            censor=TRUE, start.time, id, influence=FALSE,
            na.action=na.pass, type, ...) {

      Call <- match.call()
      Call[[1]] <- as.name("survfit")  #nicer output for the user
      object <- formula     #'formula' because it has to match survfit

      Terms  <- terms(object)
      robust <- !is.null(object$naive.var)   # did the coxph model use robust var?

      if (!is.null(attr(object$terms, "specials")$tt))
          stop("The survfit function can not process coxph models with a tt term")

      if (!missing(type)) {  # old style argument
          if (!missing(stype) || !missing(ctype))
              warning("type argument ignored")
          else {
              temp1 <- c("kalbfleisch-prentice", "aalen", "efron",
                         "kaplan-meier", "breslow", "fleming-harrington",
                         "greenwood", "tsiatis", "exact")
              
              survtype <- match(match.arg(type, temp1), temp1)
              stype <- c(1,2,2,1,2,2,2,2,2)[survtype]
              if (stype!=1) ctype <-c(1,1,2,1,1,2,1,1,1)[survtype]
          }
      }
      if (missing(ctype)) {
          # Use the appropriate one from the model
          temp1 <- match(object$method, c("exact", "breslow", "efron"))
          ctype <- c(1,1,2)[temp1]
      }
      else if (!(ctype %in% 1:2)) stop ("ctype must be 1 or 2")

      if (!se.fit) conf.type <- "none"
      else conf.type <- match.arg(conf.type)

      tfac <- attr(Terms, 'factors')
      temp <- attr(Terms, 'specials')$strata 
      has.strata <- !is.null(temp)
      if (has.strata) {
          stangle = untangle.specials(Terms, "strata")  #used multiple times, later
          # Toss out strata terms in tfac before doing the test 1 line below, as
          #  strata end up in the model with age:strat(grp) terms or *strata() terms
          #  (There might be more than one strata term)
          for (i in temp) tfac <- tfac[,tfac[i,] ==0]  # toss out strata terms
      }
      if (any(tfac >1))
          stop("not able to create a curve for models that contain an interaction without the lower order effect")

      Terms <- object$terms
      n <- object$n[1]
      if (!has.strata) strata <- NULL
      else strata <- object$strata

      missid <- missing(id) # I need this later, and setting id below makes
                            # "missing(id)" always false
      if (!missid & !missing(individual))
          warning("the `id' option supersedes `individual'")

      if (!missid) individual <- TRUE
      else if (missid && individual) id <- rep(0,n)  #dummy value
      else id <- NULL

      if (individual & missing(newdata)) {
          stop("the id and/or individual options only make sense with new data")
      }
      if (has.strata) {
          temp <- attr(Terms, "specials")$strata
          factors <- attr(Terms, "factors")[temp,]
          strata.interaction <- any(t(factors)*attr(Terms, "order") >1)
      }
      coxms <- inherits(object, "coxphms")
      if (coxms || is.null(object$y) || is.null(object[['x']]) ||
          !is.null(object$call$weights) || !is.null(object$call$id) ||
          (has.strata && is.null(object$strata)) ||
          !is.null(attr(object$terms, 'offset'))) {
          
          mf <- stats::model.frame(object)
          }
      else mf <- NULL  #useful for if statements later
      position <- NULL
      Y <- object[['y']]
      if (is.null(mf)) {
          weights <- object$weights  # let offsets/weights be NULL until needed
          offset <- NULL
          X <- object[['x']]
      }
      else {
          weights <- model.weights(mf)
          offset <- model.offset(mf)
          X <- model.matrix.coxph(object, data=mf)
          if (is.null(Y) || coxms) {
              Y <- model.response(mf)
              if (is.null(object$timefix) || object$timefix) Y <- aeqSurv(Y)
          }
          oldid <- model.extract(mf, "id")
          if (length(oldid) && ncol(Y)==3) position <- survflag(Y, oldid)
          else position <- NULL
          if (!coxms && (nrow(Y) != object$n[1])) 
              stop("Failed to reconstruct the original data set")
          if (has.strata) {
              if (length(strata)==0) {
                  if (length(stangle$vars) ==1) strata <- mf[[stangle$vars]]
                  else strata <- strata(mf[, stangle$vars], shortlabel=TRUE)
              }
          }

      }
      if (FALSE) {
      if (!is.null(mf)){
          y2 <- object[['y']]
          if (!is.null(y2)) {
              if (ncol(y2) != ncol(Y) || length(y2) != length(Y))
                  stop("Could not reconstruct the y vector")
          }
      }
      }
      type <- attr(Y, 'type')
      if (!type %in% c("right", "counting", "mright", "mcounting"))
          stop("Cannot handle \"", type, "\" type survival data")
      if (type=="right" || type== "mright") individual <- FALSE

      if (!missing(start.time)) {
          if (!is.numeric(start.time) || length(start.time) > 1)
              stop("start.time must be a single numeric value")
          # Start the curves after start.time
          # To do so, remove any rows of the data with an endpoint before that
          #  time.
          if (ncol(Y)==3) {
              keep <- Y[,2] > start.time
              Y[keep,1] <- pmax(Y[keep,1], start.time)
          }
          else keep <- Y[,1] > start.time
          if (!any(Y[keep, ncol(Y)]==1)) 
              stop("start.time argument has removed all endpoints")
          Y <- Y[keep,,drop=FALSE]
          X <- X[keep,,drop=FALSE]
          if (!is.null(offset)) offset <- offset[keep]
          if (!is.null(weights)) weights <- weights[keep]
          if (!is.null(strata))  strata <- strata[keep]
          if (length(id) >0 ) id <- id[keep]
          if (length(position) >0) position <- position[keep]
          n <- nrow(Y)
      }
      if (length(object$means) ==0) { # a model with only an offset term
          # Give it a dummy X so the rest of the code goes through
          #  (This case is really rare)
          # se.fit <- FALSE
          X <- matrix(0., nrow=n, ncol=1)
          if (is.null(offset)) offset <- rep(0, n)
          xcenter <- mean(offset)
          coef <- 0.0
          varmat <- matrix(0.0, 1, 1)
          risk <- rep(exp(offset- mean(offset)), length=n)
      }
      else {
          varmat <- object$var
          beta <- ifelse(is.na(object$coefficients), 0, object$coefficients)
          if (is.null(offset)) xcenter <- sum(object$means * beta)
          else xcenter <- sum(object$means * beta)+ mean(offset)
          if (!is.null(object$frail)) {
             keep <- !grepl("frailty(", dimnames(X)[[2]], fixed=TRUE)
             X <- X[,keep, drop=F]
          }
              
          if (is.null(offset)) risk <- c(exp(X%*% beta - xcenter))
          else     risk <- c(exp(X%*% beta + offset - xcenter))
      }
      if (missing(newdata)) {
          # If the model has interactions, print out a long warning message.
          #  People may hate it, but I don't see another way to stamp out these
          #  bad curves without backwards-incompatability.  
          # I probably should complain about factors too (but never in a strata
          #   or cluster term).
          if (any(attr(Terms, "order") > 1) )
              warning("the model contains interactions; the default curve based on columm means of the X matrix is almost certainly not useful. Consider adding a newdata argument.")
          
          if (length(object$means)) {
              mf2 <- as.list(object$means)   #create a dummy newdata
              names(mf2) <- names(object$coefficients)
              mf2 <- as.data.frame(mf2)
              x2 <- matrix(object$means, 1)
          }
          else { # nothing but an offset
              mf2 <- data.frame(X=0)
              x2 <- 0
          }
          offset2 <- 0
          found.strata <- FALSE  
      }
      else {
          if (!is.null(object$frail))
              stop("Newdata cannot be used when a model has frailty terms")

          Terms2 <- Terms 
          if (!individual)  Terms2 <- delete.response(Terms)
          if (is.vector(newdata, "numeric")) {
              if (individual) stop("newdata must be a data frame")
              if (is.null(names(newdata))) {
                  stop("Newdata argument must be a data frame")
              }
              newdata <- data.frame(as.list(newdata), stringsAsFactors=FALSE)
          }
          if (has.strata) {
              found.strata <- TRUE
              tempenv <- new.env(, parent=emptyenv())
              assign("strata", function(..., na.group, shortlabel, sep)
                  list(...), envir=tempenv)
              assign("list", list, envir=tempenv)
              for (svar in stangle$vars) {
                  temp <- try(eval(parse(text=svar), newdata, tempenv),
                              silent=TRUE)
                  if (!is.list(temp) || 
                      any(unlist(lapply(temp, class))== "function"))
                      found.strata <- FALSE
              }
              
              if (!found.strata) {
                  ss <- untangle.specials(Terms2, "strata")
                  Terms2 <- Terms2[-ss$terms]
              }
          }

          tcall <- Call[c(1, match(c('id', "na.action"), 
                                       names(Call), nomatch=0))]
          tcall$data <- newdata
          tcall$formula <- Terms2
          tcall$xlev <- object$xlevels[match(attr(Terms2,'term.labels'),
                                             names(object$xlevels), nomatch=0)]
          tcall[[1L]] <- quote(stats::model.frame)
          mf2 <- eval(tcall)
      }
      if (has.strata && found.strata) { #pull them off
          temp <- untangle.specials(Terms2, 'strata')
          strata2 <- strata(mf2[temp$vars], shortlabel=TRUE)
          strata2 <- factor(strata2, levels=levels(strata))
          if (any(is.na(strata2)))
              stop("New data set has strata levels not found in the original")
          # An expression like age:strata(sex) will have temp$vars= "strata(sex)"
          #  and temp$terms = integer(0).  This does not work as a subscript
          if (length(temp$terms) >0) Terms2 <- Terms2[-temp$terms]
      }
      else strata2 <- factor(rep(0, nrow(mf2)))

      if (!robust) cluster <- NULL
      if (individual) {
          if (missing(newdata)) 
              stop("The newdata argument must be present when individual=TRUE")
          if (!missid) {  #grab the id variable
              id2 <- model.extract(mf2, "id")
              if (is.null(id2)) stop("id=NULL is an invalid argument")
              }
          else id2 <- rep(1, nrow(mf2))
          
          x2 <- model.matrix(Terms2, mf2)[,-1, drop=FALSE]  #no intercept
          if (length(x2)==0) stop("Individual survival but no variables")

          offset2 <- model.offset(mf2)
          if (length(offset2) ==0) offset2 <- 0
                       
          y2 <- model.extract(mf2, 'response')
          if (attr(y2,'type') != type)
              stop("Survival type of newdata does not match the fitted model")
          if (attr(y2, "type") != "counting")
              stop("Individual=TRUE is only valid for counting process data")
          y2 <- y2[,1:2, drop=F]  #throw away status, it's never used
      }
      else if (missing(newdata)) {
          if (has.strata && strata.interaction)
              stop ("Models with strata by covariate interaction terms require newdata")
          offset2 <- 0
          if (length(object$means)) {
              x2 <- matrix(object$means, nrow=1, ncol=ncol(X))
          } else {
              # model with only an offset and no new data: very rare case 
              x2 <- matrix(0.0, nrow=1, ncol=1)   # make a dummy x2
          }
      } else {
          offset2 <- model.offset(mf2)
          if (length(offset2) >0) offset2 <- offset2 
          else offset2 <- 0
          x2 <- model.matrix(Terms2, mf2)[,-1, drop=FALSE]  #no intercept
      }
      if (missing(newdata)) risk2 <- 1
      else {
          if (length(object$means)) 
              risk2 <- exp(c(x2 %*% beta) + offset2 - xcenter)
          else risk2 <- exp(offset2 - xcenter)
      }
      if (individual) {
          result <- coxsurv.fit(ctype, stype, se.fit, varmat, cluster, 
                                 Y, X, weights, risk, position, strata, oldid,
                                 y2, x2, risk2, strata2, id2)
      }
      else {
          result <- coxsurv.fit(ctype, stype, se.fit, varmat, cluster, 
                                 Y, X, weights, risk, position, strata, oldid,
                                 y2, x2, risk2)
          if (has.strata && found.strata) {
              if (is.matrix(result$surv)) {
                  nr <- nrow(result$surv)  #a vector if newdata had only 1 row
                  indx1 <- split(1:nr, rep(1:length(result$strata), result$strata))
                  rows <- indx1[as.numeric(strata2)]  #the rows for each curve

                  indx2 <- unlist(rows)  #index for time, n.risk, n.event, n.censor
                  indx3 <- as.integer(strata2) #index for n and strata

                  for(i in 2:length(rows)) rows[[i]] <- rows[[i]]+ (i-1)*nr #linear subscript
                  indx4 <- unlist(rows)   #index for surv and std.err
                  temp <- result$strata[indx3]
                  names(temp) <- row.names(mf2)
                  new <- list(n = result$n[indx3],
                              time= result$time[indx2],
                              n.risk= result$n.risk[indx2],
                              n.event=result$n.event[indx2],
                              n.censor=result$n.censor[indx2],
                              strata = temp,
                              surv= result$surv[indx4],
                              cumhaz = result$cumhaz[indx4])
                  if (se.fit) new$std.err <- result$std.err[indx4]
                  result <- new
              }
          }
      }
      if (!censor) {
          kfun <- function(x, keep){ if (is.matrix(x)) x[keep,,drop=F] 
                                    else if (length(x)==length(keep)) x[keep]
                                    else x}
          keep <- (result$n.event > 0)
          if (!is.null(result$strata)) {
              temp <- factor(rep(names(result$strata), result$strata),
                             levels=names(result$strata))
              result$strata <- c(table(temp[keep]))
              }
          result <- lapply(result, kfun, keep)
          }
      result$logse = TRUE   # this will migrate further in

      if (se.fit && conf.type != "none") {
          ci <- survfit_confint(result$surv, result$std.err, logse=result$logse,
                                conf.type, conf.int)
          result <- c(result, list(lower=ci$lower, upper=ci$upper, 
                                   conf.type=conf.type, conf.int=conf.int))
      }

      if (!missing(start.time)) result$start.time <- start.time

      result$call <- Call
      class(result) <- c('survfitcox', 'survfit')
      result
      }
