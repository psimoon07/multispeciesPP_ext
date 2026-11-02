## Survey-source detection offsets: source.PA / source.RL / source.ref.  Each
## non-reference source s gets a species-specific is<s> coefficient that enters the
## linear predictor of that source's PA or RL rows only, so exp(is<s>) is the share
## of its nominal exposure the source effectively searches, relative to source.ref.

src.setup <- function(p.ZK=0.2, p.RL=0.1, seed=2024) {
    land <- sim.landscape(n.cells=4000, seed=7)
    coefs <- rbind(sp1=c(-2.0, 0.8, -0.5), sp2=c(-2.5, -0.4, 0.9))
    lam <- sim.lambda(land, coefs)
    set.seed(seed)
    PO <- sim.PO(land, lam, gamma=c(sp1=1.2, sp2=1.5), delta=0.6)
    ## ST is the reference: it searches its nominal area.  ZK and RL search only a
    ## share p of theirs, but are given to the model with their nominal exposure.
    st <- sim.PA(land, lam, sites=sample(nrow(land), 900), quadrat.size=2)
    zk <- sim.PA(land, lam, sites=sample(nrow(land), 900), quadrat.size=8*p.ZK)
    zk$quadrat.size <- 8
    PA <- rbind(st, zk); rownames(PA) <- paste0("pa", seq_len(nrow(PA)))
    rl <- sim.RL(land, lam, plots=sample(nrow(land), 1500), weight=20*p.RL)
    rl$BG.RL$quadrat.size <- 20
    list(land=land, coefs=coefs, PO=PO, PA=PA,
         source.PA=rep(c("ST","ZK"), each=900),
         RL=rl$RL, BG.RL=rl$BG.RL, region.size=nrow(land))
}

src.fit <- function(d, ...) {
    suppressMessages(suppressWarnings(
        multispeciesPP(~x1+x2, ~z1, PA=d$PA, PO=d$PO, BG=d$land,
                       RL=d$RL, BG.RL=d$BG.RL,
                       quadrat.size=d$PA$quadrat.size,
                       quadrat.size.RL=d$BG.RL$quadrat.size,
                       region.size=d$region.size, ...)))
}

test_that("source offsets recover the detection shares and the true lambda", {
    d <- src.setup()
    fit <- src.fit(d, source.PA=d$source.PA, source.RL="RL", source.ref="ST",
                   inverse.hessian=TRUE)

    expect_equal(fit$source.ref, "ST")
    expect_equal(fit$source.levels, c("RL","ZK"))
    expect_equal(rownames(fit$species.coef),
                 c("(Intercept)","x1","x2","isRL","isZK","isPO"))

    expect_lt(max(abs(fit$species.coef["isZK",] - log(0.2))), 0.25)
    expect_lt(max(abs(fit$species.coef["isRL",] - log(0.1))), 0.25)
    expect_lt(max(abs(fit$species.coef["(Intercept)",] - d$coefs[,1])), 0.15)
    expect_lt(max(abs(fit$species.coef["x1",] - d$coefs[,2])), 0.1)
    expect_lt(max(abs(fit$species.coef["x2",] - d$coefs[,3])), 0.1)
    expect_lt(max(abs(fit$species.coef["isPO",] - c(1.2, 1.5))), 0.2)
    expect_lt(abs(fit$bias.coef - 0.6), 0.1)

    ## offsets are unstandardized indicators: no rescaling, no mean shift, so the
    ## normalized and original-scale values and standard errors coincide
    is.src <- grepl(":is(RL|ZK)$", names(fit$all.coef))
    expect_equal(fit$all.coef[is.src], fit$normalized.all.coef[is.src])
    expect_equal(sqrt(diag(fit$vcov))[is.src], fit$std.errs[is.src])
})

test_that("ignoring the detection shares biases alpha, the offsets remove that bias", {
    d <- src.setup()
    naive <- src.fit(d)
    fit <- src.fit(d, source.PA=d$source.PA, source.RL="RL", source.ref="ST")
    err <- function(f) max(abs(f$species.coef["(Intercept)",] - d$coefs[,1]))
    expect_gt(err(naive), 0.5)
    expect_lt(err(fit), 0.15)
})

test_that("predict() leaves the source offsets out of lambda", {
    d <- src.setup()
    fit <- src.fit(d, source.PA=d$source.PA, source.RL="RL", source.ref="ST")
    nd <- d$land[1:20,]
    manual <- cbind(1, nd$x1, nd$x2) %*% fit$species.coef[c("(Intercept)","x1","x2"),]
    expect_equal(unname(predict(fit, nd)), unname(manual))
    ## with bias, isPO is added but still no source offset
    expect_equal(unname(predict(fit, nd, bias=TRUE)),
                 unname(manual + c(nd$z1 * fit$bias.coef) +
                        rep(fit$species.coef["isPO",], each=20)))
})

test_that("a source labelling only the reference gives the fit without sources", {
    d <- src.setup()
    d$PA <- d$PA[d$source.PA == "ST",]
    plain <- src.fit(d)
    same <- src.fit(d, source.PA=rep("ST", nrow(d$PA)), source.ref="ST")
    expect_equal(same$all.coef, plain$all.coef)
    expect_equal(same$source.levels, character(0))
})

test_that("unlabelled blocks count as the reference source", {
    ## only RL labelled: PA (both halves) is the reference
    d <- src.setup(p.ZK=1)
    fit <- src.fit(d, source.RL="RL", source.ref="ST")
    expect_equal(fit$source.levels, "RL")
    expect_lt(max(abs(fit$species.coef["isRL",] - log(0.1))), 0.25)
})

test_that("source arguments are validated", {
    d <- src.setup()
    expect_error(src.fit(d, source.PA=d$source.PA), "source.ref")
    expect_error(src.fit(d, source.PA=d$source.PA[-1], source.ref="ST"), "every row of PA")
    expect_error(src.fit(d, source.RL=c("RL","RL"), source.ref="ST"), "single source label")
    expect_error(src.fit(d, source.PA=d$source.PA, source.RL="RL", source.ref="XX"),
                 "labels no rows")
    expect_error(src.fit(d, source.PA=sub("ZK","PO",d$source.PA), source.ref="ST"),
                 "reserved")
})
