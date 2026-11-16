library(survival)

# A full-rank covariance must not lose a term when its units change (#353).
v <- matrix(c(4, 1, 1, 9), 2)
b <- cbind(c(1, 2), c(-1, 3))
answer <- matrix(c(7, 7, -12, 13), 2) / 35
for (s in list(c(1, 1), c(1e-7, 1), c(1, 1e7))) {
    z <- coxph.wtest(v * outer(s, s), b * s)
    stopifnot(z$df == 2,
              isTRUE(all.equal(z$test, c(3/5, 51/35))),
              isTRUE(all.equal(z$solve * s, answer)))
}

# Singular covariance matrices retain the ordered generalized solution.
for (s in list(c(1, 1), c(1e-7, 1), c(1, 1e7))) {
    z <- coxph.wtest(matrix(1, 2, 2) * outer(s, s), c(2, 2) * s)
    stopifnot(z$df == 1, isTRUE(all.equal(z$test, 4)),
              isTRUE(all.equal(z$solve * s, c(2, 0))))
    z <- coxph.wtest(diag(c(1, 0)) * outer(s, s), c(2, 0) * s)
    stopifnot(z$df == 1, isTRUE(all.equal(z$test, 4)),
              isTRUE(all.equal(z$solve * s, c(2, 0))))
}
z <- coxph.wtest(matrix(0, 2, 2), c(0, 0))
stopifnot(z$df == 0, z$test == 0, identical(z$solve, c(0, 0)))
z <- coxph.wtest(matrix(c(1, 1-1e-12, 1-1e-12, 1), 2), c(1, 1))
stopifnot(z$df == 1, isTRUE(all.equal(z$test, 1)))

# Vector, one-column matrix, scalar and empty cases keep their return shapes.
z <- coxph.wtest(v, b[, 1])
stopifnot(is.null(dim(z$solve)), isTRUE(all.equal(z$solve, answer[, 1])))
z <- coxph.wtest(v, b[, 1, drop=FALSE])
stopifnot(is.null(dim(z$solve)), isTRUE(all.equal(z$solve, answer[, 1])))
dimnames(v) <- list(c("a", "b"), c("a", "b"))
z <- coxph.wtest(v, b[, 1])
stopifnot(is.null(names(z$solve)))
z <- coxph.wtest(4, 2)
stopifnot(z$df == 1, z$test == 1, z$solve == 0.5)
z <- coxph.wtest(4, matrix(c(2, 4), 1))
stopifnot(identical(dim(z$solve), c(1L, 2L)),
          identical(as.vector(z$solve), c(0.5, 1)))
z <- coxph.wtest(matrix(numeric(0), 0, 0), numeric(0))
stopifnot(z$df == 0, identical(z$test, numeric(0)), z$solve == 0)
stopifnot(inherits(try(coxph.wtest(diag(2), c(1, Inf)), silent=TRUE),
                   "try-error"),
          inherits(try(coxph.wtest(diag(2), 1:3), silent=TRUE), "try-error"))

# Cox Wald and robust score tests must be invariant, including counting data.
d <- na.omit(lung[, c("time", "status", "age", "sex")])
d$start <- 0
d$id <- seq_len(nrow(d))
for (counting in c(FALSE, TRUE)) for (robust in c(FALSE, TRUE)) {
    form <- if (counting) Surv(start, time, status) ~ age + sex else
        Surv(time, status) ~ age + sex
    fit <- coxph(form, d, robust=robust, id=id)
    beta <- coef(fit)
    variance <- vcov(fit)
    reference <- sum(beta * solve(variance, beta))
    for (k in c(1e-8, 1e5)) {
        scaled <- d
        scaled$age <- scaled$age * k
        fit2 <- coxph(form, scaled, robust=robust, id=id)
        stopifnot(isTRUE(all.equal(fit2$wald.test, reference)),
                  isTRUE(all.equal(fit2$linear.predictors, fit$linear.predictors)))
        if (robust) stopifnot(isTRUE(all.equal(fit2$rscore, fit$rscore)))
    }
}

# An aliased coefficient remains NA and does not add a test degree of freedom.
d$copy <- d$age
fit <- coxph(Surv(time, status) ~ age + sex, d)
fit2 <- coxph(Surv(time, status) ~ age + sex + copy, d)
stopifnot(is.na(coef(fit2)["copy"]), summary(fit2)$waldtest["df"] == 2,
          isTRUE(all.equal(fit2$wald.test, fit$wald.test)))
