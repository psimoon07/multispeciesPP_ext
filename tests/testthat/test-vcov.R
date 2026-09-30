## inverse.hessian=TRUE returns the covariance matrix of the coefficients: the matrix
## whose diagonal std.errs is already taken from, on the standardized scale
## (normalized.vcov) and mapped to the original covariate scale (vcov).

vcov.data <- function(shift=0) {
    land <- sim.landscape(n.cells=3000, seed=31)
    coefs <- rbind(sp1=c(-2.0, 0.8, -0.5), sp2=c(-2.5, -0.4, 0.9))
    lam <- sim.lambda(land, coefs)
    set.seed(13)
    PO <- sim.PO(land, lam, gamma=c(sp1=1.2, sp2=1.5), delta=0.6)
    PA <- sim.PA(land, lam, sites=sample(nrow(land), 800), quadrat.size=3)
    move <- function(d) { d$x1 <- d$x1 + shift; d }
    list(land=move(land), PO=lapply(PO, move), PA=move(PA))
}

vcov.fit <- function(d, ...) {
    multispeciesPP(~x1+x2, ~z1, PA=d$PA, PO=d$PO, BG=d$land, quadrat.size=3,
                   region.size=nrow(d$land)*d$land$cell.area[1], ...)
}

test_that("the covariance is only returned on request", {
    fit <- vcov.fit(vcov.data())
    expect_null(fit$vcov)
    expect_null(fit$normalized.vcov)
})

test_that("asking for the covariance does not change the fit", {
    d <- vcov.data()
    expect_equal(vcov.fit(d, inverse.hessian=TRUE)$all.coef, vcov.fit(d)$all.coef)
})

test_that("the diagonal of the covariance is the squared standard errors", {
    fit <- vcov.fit(vcov.data(), inverse.hessian=TRUE)
    expect_equal(dimnames(fit$vcov), list(names(fit$all.coef), names(fit$all.coef)))
    expect_equal(fit$vcov, t(fit$vcov))
    expect_equal(sqrt(diag(fit$normalized.vcov)), fit$normalized.std.errs)

    ## on the original scale std.errs ignores the shift of the intercepts, vcov does not
    shifted <- grepl(":(\\(Intercept\\)|isPO)$", names(fit$all.coef))
    expect_equal(sqrt(diag(fit$vcov))[!shifted], fit$std.errs[!shifted])
})

test_that("the covariance follows the coefficients when a covariate is shifted", {
    ## x1 + 3 changes nothing but the intercepts: alpha' = alpha - 3*beta_x1.  The
    ## variance of alpha' in the shifted fit must then equal the variance of that
    ## contrast in the unshifted fit, which needs the off-diagonal terms to be right.
    fit <- vcov.fit(vcov.data(), inverse.hessian=TRUE)
    moved <- vcov.fit(vcov.data(shift=3), inverse.hessian=TRUE)

    for(sp in c("sp1","sp2")) {
        a <- paste0(sp,":(Intercept)")
        b <- paste0(sp,":x1")
        expect_equal(moved$all.coef[[a]], fit$all.coef[[a]] - 3*fit$all.coef[[b]])
        expect_equal(moved$vcov[a,a],
                     fit$vcov[a,a] + 9*fit$vcov[b,b] - 6*fit$vcov[a,b])
        expect_equal(moved$vcov[b,b], fit$vcov[b,b])
    }
})
