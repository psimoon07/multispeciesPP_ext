multispeciesPP <-
function(sdm.formula, bias.formula, PA=NULL, PO=NULL, BG,
             species=union(names(PO),names(RL)),
             species.PA=species,species.PO=species,
         RL=NULL, BG.RL=NULL, species.RL=names(RL),
         quadrat.size=1,region.size=1,quadrat.size.RL=1,
         start=NULL,inverse.hessian=FALSE,
         penalty.l2.sdm=.1,penalty.l2.bias=.1,
         penalty.l2.intercept=1E-4,  # should be small but nonzero, protects against singularities
         weights=rep(1,n.species*nrow(x)),
         control=list()) {
    control <- do.call("glm.control", control)

    ## RL detections and RL quadrature must be supplied together: the detections
    ## give the sum term of the point-process likelihood and the quadrature rows
    ## give the integral term.  One without the other is not a likelihood.
    if(is.null(RL) != is.null(BG.RL))
        stop("RL and BG.RL must both be supplied, or both be NULL")

    if(is.null(PA)) species.PA <- character(0)
    if(is.null(PO)) species.PO <- character(0)
    if(is.null(RL)) species.RL <- character(0)
    species <- union(union(species.PO,species.PA),species.RL)
    if(length(species) == 0)
        stop("no species have data; supply at least one of PA, PO or RL (and 'species' if it cannot be inferred)")

    use.PA <- !is.null(PA)
    use.RL <- !is.null(BG.RL)

    if(use.PA && !(length(quadrat.size) %in% c(1L,nrow(PA))))
        stop("quadrat.size must be a scalar or a numeric vector of length nrow(PA)")
    if(use.RL && !(length(quadrat.size.RL) %in% c(1L,nrow(BG.RL))))
        stop("quadrat.size.RL must be a scalar or a numeric vector of length nrow(BG.RL)")

    ## alpha is shared across all sources, so every exposure must be measured the
    ## same way.  A mismatch still converges; it just returns the wrong alpha.
    if(use.RL || length(quadrat.size) > 1)
        message("multispeciesPP: region.size, quadrat.size and quadrat.size.RL must all be ",
                "expressed in the same space-time units (e.g. m^2 * day); lambda is an ",
                "intensity per unit of that measure.")

    sdm.formula <- update(sdm.formula,~.+1)   # ensure intercept for sdm model
    bias.formula <- update(bias.formula,~.-1) # get rid of intercept for bias model


    sdm.mf <- model.frame(sdm.formula,data=BG)
    bias.mf <- model.frame(bias.formula,data=BG)

    ## Extract data and construct standardized matrices.  Note that model.matrix automatically omits NAs!
    ## Note that we need the BG data to get our standardization, whether or not we have any PO data.
    sdm.BG.model.matrix <- model.matrix(terms(sdm.mf), BG)
    sdm.means <- c(0,apply(sdm.BG.model.matrix[,-1,drop=FALSE],2,mean))
    sdm.BG.model.matrix <- sweep(sdm.BG.model.matrix,2,sdm.means,"-")
    sdm.sds <- c(1,apply(sdm.BG.model.matrix[,-1,drop=FALSE],2,sd))
    sdm.BG.model.matrix <- sweep(sdm.BG.model.matrix,2,sdm.sds,"/")
    sdm.standardize <- function(mat) sweep(sweep(mat,2,sdm.means,"-"),2,sdm.sds,"/")

    bias.BG.model.matrix <- model.matrix(terms(bias.mf), BG)
    bias.means <- apply(bias.BG.model.matrix,2,mean)
    bias.BG.model.matrix <- sweep(bias.BG.model.matrix,2,bias.means,"-")
    bias.sds <- apply(bias.BG.model.matrix,2,sd)
    bias.BG.model.matrix <- sweep(bias.BG.model.matrix,2,bias.sds,"/")
    bias.standardize <- function(mat) sweep(sweep(mat,2,bias.means,"-"),2,bias.sds,"/")

    BG.good.rows <- intersect(rownames(sdm.BG.model.matrix),
                              rownames(bias.BG.model.matrix))

    if(use.PA) {
        sdm.PA.model.matrix <- sdm.standardize(model.matrix(terms(sdm.mf), PA))
        PA.good.rows <- rownames(sdm.PA.model.matrix)
    } else {
        sdm.PA.model.matrix <- sdm.BG.model.matrix[0,,drop=FALSE]
        PA.good.rows <- character(0)
    }

    ## PA.good.rows is a character vector of *rownames*, and can be a strict subset of
    ## rownames(PA) because model.matrix drops rows with NA covariates.  match() aligns
    ## a per-row quadrat.size to it whether or not the vector carries names.
    if(length(quadrat.size) > 1)
        quadrat.size <- quadrat.size[match(PA.good.rows,rownames(PA))]

    if(length(species.PO) > 0) {
        sdm.PO.model.matrices <- lapply(as.list(species.PO),
                                        function(sp) sdm.standardize(model.matrix(terms(sdm.mf), PO[[sp]])))
        names(sdm.PO.model.matrices) <- species.PO
        bias.PO.model.matrices <- lapply(as.list(species.PO),
                                         function(sp) bias.standardize(model.matrix(terms(bias.mf), PO[[sp]])))
        names(bias.PO.model.matrices) <- species.PO

        PO.good.rows <- lapply(as.list(species.PO),
                               function(sp) intersect(rownames(sdm.PO.model.matrices[[sp]]),
                                                      rownames(bias.PO.model.matrices[[sp]])))
            names(PO.good.rows) <- species.PO
    }

    ## RL detections need only the sdm covariates: there is no bias term for an
    ## exhaustively censused source, so no bias model matrix is built for them.
    if(length(species.RL) > 0) {
        sdm.RL.model.matrices <- lapply(as.list(species.RL),
                                        function(sp) sdm.standardize(model.matrix(terms(sdm.mf), RL[[sp]])))
        names(sdm.RL.model.matrices) <- species.RL
        RL.good.rows <- lapply(sdm.RL.model.matrices, rownames)
        names(RL.good.rows) <- species.RL
    }

    if(use.RL) {
        sdm.BG.RL.model.matrix <- sdm.standardize(model.matrix(terms(sdm.mf), BG.RL))
        BG.RL.good.rows <- rownames(sdm.BG.RL.model.matrix)
    } else {
        sdm.BG.RL.model.matrix <- sdm.BG.model.matrix[0,,drop=FALSE]
        BG.RL.good.rows <- character(0)
    }
    if(length(quadrat.size.RL) > 1)
        quadrat.size.RL <- quadrat.size.RL[match(BG.RL.good.rows,rownames(BG.RL))]

    n.species <- length(species)
    p.sdm <- ncol(sdm.BG.model.matrix)-1  # not including intercept
    p.bias <- ncol(bias.BG.model.matrix)

    ## Compress all PO examples into margins
    sdm.margins.ab <- matrix(0,n.species,p.sdm+1,dimnames=list(species,colnames(sdm.BG.model.matrix)))
    sdm.margins.gamma <- matrix(0,n.species,1,dimnames=list(species,"isPO"))
    bias.margins <- matrix(0,1,p.bias,dimnames=list(NULL,colnames(bias.BG.model.matrix)))
    for(sp in species.PO) {
        k <- match(sp,species)
        sdm.margins.ab[k,] <- colSums(sdm.PO.model.matrices[[sp]][PO.good.rows[[sp]],,drop=FALSE])
        sdm.margins.gamma[k,] <- length(PO.good.rows[[sp]])
        bias.margins <- bias.margins + colSums(bias.PO.model.matrices[[sp]][PO.good.rows[[sp]],,drop=FALSE])
    }

    ## Compress all RL detections into margins too.  The margin rows use the linear()
    ## pseudo-family, whose score contribution for a y=1 row is exactly the row itself,
    ## so summing RL into sdm.margins.ab is *identical* to carrying separate RL margin
    ## rows -- and it leaves every downstream row index untouched.  RL contributes
    ## nothing to sdm.margins.gamma (no gamma_RL) or to bias.margins (no delta).
    sdm.margins.ab.RL <- matrix(0,n.species,p.sdm+1,dimnames=list(species,colnames(sdm.BG.model.matrix)))
    for(sp in species.RL) {
        k <- match(sp,species)
        sdm.margins.ab.RL[k,] <- colSums(sdm.RL.model.matrices[[sp]][RL.good.rows[[sp]],,drop=FALSE])
    }
    sdm.margins.ab <- sdm.margins.ab + sdm.margins.ab.RL

    if(length(species.RL) > 0) {
        no.detections <- species.RL[vapply(species.RL,
                                           function(sp) all(sdm.margins.ab.RL[match(sp,species),] == 0) &&
                                                        !(sp %in% species.PA) && !(sp %in% species.PO),
                                           logical(1))]
        if(length(no.detections) > 0)
            warning("species with RL as their only data source but no RL detections: ",
                    paste(no.detections,collapse=", "),
                    ".  Their alpha is determined by penalty.l2.intercept alone.")
    }

    abcd.from.all.coef <- function(all.coef) {
        sdm.coef <- matrix(all.coef[1:(n.species*(p.sdm+2))],p.sdm+2,n.species)
        alpha <- sdm.coef[1,]
        beta <- t(sdm.coef[2:(p.sdm+1),,drop=FALSE])
        gamma <- sdm.coef[p.sdm+2,]
        delta <- all.coef[-(1:(n.species*(p.sdm+2)))]

        names(alpha) <- names(gamma) <- species
        colnames(beta) <- colnames(sdm.margins.ab)[-1]
        rownames(beta) <- species
        names(delta) <- colnames(bias.BG.model.matrix)
        return(list(alpha=alpha,beta=beta,gamma=gamma,delta=delta))
    }

    all.coef.from.abcd <- function(alpha,beta,gamma,delta) {
        c(rbind(alpha,beta,gamma),delta)
    }

    n.PA <- length(PA.good.rows)
    n.BG <- length(BG.good.rows)
    n.RL <- length(BG.RL.good.rows)
    subsamp.PA.offset <- 0
    subsamp.BG.offset <- 0
    n.sites <- n.BG+n.PA+n.RL

    ## Combine matrices into X and Z.  Row blocks, in order:
    ##   n.species sdm margins | 1 bias margin | n.PA | n.BG | n.RL | penalty
    ## The last column of x is the isPO indicator: 1 on BG rows only.  RL rows get 0,
    ## because an exhaustive census has no bias intercept -- RL's linear predictor is
    ## exactly alpha_k + beta_k' X(u), the same coefficients every other source uses.
    x <- cbind(rbind(sdm.margins.ab,    # extra rows for the positive obs
                     0,
                     sdm.PA.model.matrix[PA.good.rows,,drop=FALSE],
                     sdm.BG.model.matrix[BG.good.rows,,drop=FALSE],
                     sdm.BG.RL.model.matrix[BG.RL.good.rows,,drop=FALSE]),
               c(sdm.margins.gamma,rep(c(0,1,0),c(1+n.PA,n.BG,n.RL))))

    ## add pseudodata for penalty on (non-intercept) beta and delta
    x <- rbind(x,
               diag(sqrt(c(penalty.l2.intercept, rep(penalty.l2.sdm,p.sdm), penalty.l2.intercept))),
               matrix(0,p.bias,p.sdm+2))

    z <- rbind(matrix(0,n.species,p.bias),
               bias.margins,
               matrix(0,n.PA,p.bias),
               bias.BG.model.matrix[BG.good.rows,,drop=FALSE],
               matrix(0,n.RL,p.bias),   # no bias covariates enter the RL block
               matrix(0,p.sdm+2,p.bias),
               sqrt(penalty.l2.bias/n.species)*diag(p.bias))

    y <- rep(0,nrow(x)*n.species)
    offset <- rep(0,nrow(x)*n.species)
    for(k in 1:n.species) {
        yk <- rep(0,nrow(x))
        yk[1:n.species] <- 1*(1:n.species == k)     # sdm margin only active for species k
        yk[1 + n.species] <- 1*(1 == k)             # bias margin only active for species no. 1

        if(n.PA > 0) {
            if(species[k] %in% species.PA) {
                yk[1 + n.species + (1:n.PA)] <- PA[PA.good.rows,species[k]]  # PA sites
            } else {
                yk[1 + n.species + (1:n.PA)] <- NA
            }
        }
        if(species[k] %in% species.PO) {
            yk[1 + n.species + n.PA + (1:n.BG)] <- 0             # BG sites
        } else {
            yk[1 + n.species + n.PA + (1:n.BG)] <- NA
        }
        ## RL quadrature sites.  y is 0 whether or not that plot-visit detected species k:
        ## a non-detection enters the likelihood only through the -integral(lambda) term,
        ## which is exactly what these rows carry.  Species without RL data are masked out
        ## with NA, or their alpha would absorb an integral with no matching detections.
        if(n.RL > 0) {
            if(species[k] %in% species.RL) {
                yk[1 + n.species + n.PA + n.BG + (1:n.RL)] <- 0
            } else {
                yk[1 + n.species + n.PA + n.BG + (1:n.RL)] <- NA
            }
        }
        yk[1 + n.species + n.sites + (1:(p.sdm + 2 + p.bias))] <- 0  # penalty

        y[(k-1)*nrow(x) + 1:nrow(x)] <- yk

        offk <- rep(0,nrow(x))
        if(n.PA > 0) offk[1 + n.species + (1:n.PA)] <- log(quadrat.size)
        offk[1 + n.species + n.PA + (1:n.BG)] <- log(region.size) - log(n.BG)  # BG sites
        if(n.RL > 0) offk[1 + n.species + n.PA + n.BG + (1:n.RL)] <- log(quadrat.size.RL)
        offset[(k-1)*nrow(x) + 1:nrow(x)] <- offk
    }

    block.index <- function(start.row,n.rows) {
        if(n.rows == 0) return(integer(0))
        (start.row + 1):(start.row + n.rows) + rep((0:(n.species-1))*nrow(x),each=n.rows)
    }
    which.PA <- block.index(1 + n.species, n.PA)
    which.BG <- block.index(1 + n.species + n.PA, n.BG)
    which.RL <- block.index(1 + n.species + n.PA + n.BG, n.RL)

    ## Initialize coefficients and fit model
    if(is.null(start)) {
        start.alpha <- start.gamma <- rep(0,n.species)
        total.quadrat.size <- if(n.PA > 0) sum(rep(quadrat.size,length.out=n.PA)) else 0
        total.quadrat.size.RL <- if(n.RL > 0) sum(rep(quadrat.size.RL,length.out=n.RL)) else 0
        for(k in 1:n.species) {
            if((species[k] %in% species.PA) &&
               sum(!is.na(PA[PA.good.rows,species[k]])>0))
                start.alpha[k] <- log((1+sum(PA[PA.good.rows,species[k]],na.rm=TRUE))/total.quadrat.size)
            else if((species[k] %in% species.RL) && total.quadrat.size.RL > 0)
                ## without this, a species with no PA data starts at lambda = 1 per unit
                ## of exposure, which for plot-sized weights means a wildly divergent
                ## first Poisson IRLS step
                start.alpha[k] <- log((1+length(RL.good.rows[[species[k]]]))/total.quadrat.size.RL)
            if(species[k] %in% species.PO)  start.gamma[k] <- log1p(sdm.margins.gamma[k,]) - start.alpha[k]- log(region.size)

        }
        start <- all.coef.from.abcd(start.alpha,
                                    matrix(0,p.sdm,n.species),
                                    start.gamma,
                                    rep(0,p.bias))
    }
    fit <- block.glm.fit(x,z,y,weights=weights,start=start,offset=offset,
                         families=list(linear(),binomial(link="cloglog"),poisson(),gaussian()),
                         row.families=rep(rep(1:4,c(1+n.species,n.PA,n.BG+n.RL,p.sdm+p.bias+2)),n.species),
                         control=control,inverse.hessian=inverse.hessian)

    ## Postprocess fit
    all.coef <- fit$coefficients
    eta <- fit$linear.predictors
    mu <- fit$fitted.values
    names(all.coef)[1:(n.species*(p.sdm+2))] <-
        paste(rep(species,each=p.sdm+2),c(colnames(sdm.BG.model.matrix)[1:(p.sdm+1)],"isPO"),sep=":")
    names(all.coef)[-(1:(n.species*(p.sdm+2)))] <- paste("isPO:",colnames(bias.BG.model.matrix),sep="")

    std.errs <- fit$fit$std.errs
    names(std.errs) <- names(all.coef)

    species.coef <- matrix(all.coef[1:(n.species*(p.sdm+2))],p.sdm+2,n.species,
                           dimnames=list(c(colnames(sdm.margins.ab),"isPO"),species))
    bias.coef <- all.coef[-(1:(n.species*(p.sdm+2)))]
    names(bias.coef) <- colnames(bias.BG.model.matrix)

    n.PA.rows <- if(use.PA) nrow(PA) else 0
    PA.rownames <- if(use.PA) dimnames(PA)[[1]] else NULL
    fit.PA <- linear.fit.PA <- matrix(NA,n.PA.rows,length(species),dimnames=list(PA.rownames,species))
    if(n.PA > 0) {
        linear.fit.PA[PA.good.rows,] <- eta[which.PA]
        fit.PA[PA.good.rows,] <- mu[which.PA]
    }

    fit.BG <- linear.fit.BG <-
        bias.fit.BG <- linear.bias.fit.BG <- matrix(NA,nrow(BG),length(species),
                                                            dimnames=list(dimnames(BG)[[1]],species))
    linear.fit.BG[BG.good.rows,] <- matrix(eta[which.BG],ncol=n.species) + log(n.BG) - log(region.size)
    fit.BG[BG.good.rows,] <- matrix(mu[which.BG],ncol=n.species) * n.BG / region.size
    linear.bias.fit.BG[BG.good.rows,] <- c(bias.BG.model.matrix[BG.good.rows,,drop=FALSE] %*% bias.coef)
    bias.fit.BG[BG.good.rows,] <- exp(linear.bias.fit.BG[BG.good.rows,])

    ## lambda_k at the RL quadrature points, on a per-unit-exposure scale
    n.BG.RL.rows <- if(use.RL) nrow(BG.RL) else 0
    BG.RL.rownames <- if(use.RL) dimnames(BG.RL)[[1]] else NULL
    fit.RL <- linear.fit.RL <- matrix(NA,n.BG.RL.rows,length(species),
                                      dimnames=list(BG.RL.rownames,species))
    if(n.RL > 0) {
        qs.RL <- rep(quadrat.size.RL,length.out=n.RL)
        linear.fit.RL[BG.RL.good.rows,] <- matrix(eta[which.RL],ncol=n.species) - log(qs.RL)
        fit.RL[BG.RL.good.rows,] <- matrix(mu[which.RL],ncol=n.species) / qs.RL
    }

    fitted.sdm.margins.gamma <- colSums(fit.BG[BG.good.rows,,drop=FALSE])*region.size/n.BG
    fitted.bias.margins <- colSums(t(fit.BG[BG.good.rows,species.PO,drop=FALSE]) %*%
                                   bias.BG.model.matrix[BG.good.rows,,drop=FALSE]*region.size/n.BG)

    ## RL touches neither gamma nor delta (isPO column is 0 and z is 0 on its rows), so
    ## both score checks below are unchanged by the RL block.
    score.check.gamma <- fitted.sdm.margins.gamma - sdm.margins.gamma + penalty.l2.intercept*species.coef[p.sdm+2,]
    score.check.gamma <- score.check.gamma[species %in% species.PO]

    score.check.bias <- fitted.bias.margins - bias.margins + penalty.l2.bias*bias.coef

    ##stopifnot(mean(score.check.ab^2) < f.tol)
    if(length(score.check.gamma)>0) stopifnot(mean((score.check.gamma/fit$deviance)^2) < control$epsilon)
    stopifnot(mean((score.check.bias/fit$deviance)^2) < control$epsilon)

    sd.normalizer <- c(rep(c(sdm.sds,1),n.species),bias.sds)
    unstandardized.coef <- all.coef/sd.normalizer
    gamma.adjust <- sum(unstandardized.coef[-(1:(n.species*(p.sdm+2)))] * bias.means)
    for(k in 1:n.species) {
        jk <- (p.sdm+2)*(k-1) + 1:(p.sdm+1)
        coef.block <- unstandardized.coef[jk]
        unstandardized.coef[jk[1]] <- coef.block[1] - sum(coef.block[-1] * sdm.means[-1])       # adjust alpha(k)
        unstandardized.coef[jk[1]+p.sdm+1] <- unstandardized.coef[jk[1]+p.sdm+1] - gamma.adjust # adjust gamma(k)
    }
    unstandardized.species.coef <- matrix(unstandardized.coef[1:(n.species*(p.sdm+2))],p.sdm+2,n.species,
                                          dimnames=list(c(colnames(sdm.margins.ab),"isPO"),species))
    unstandardized.bias.coef <- unstandardized.coef[-(1:(n.species*(p.sdm+2)))]
    names(unstandardized.bias.coef) <- colnames(bias.BG.model.matrix)



    ## terms(*.mf) carry 'predvars': the knots/centres of data-dependent terms
    ## (ns(), pbs(), poly(), scale()) as fixed on BG.  predict() must rebuild the
    ## design from these, not from the bare formula, or it re-derives them from newdata.
    tr <- list(sdm.formula=sdm.formula,bias.formula=bias.formula,
               sdm.terms=terms(sdm.mf),bias.terms=terms(bias.mf),
               normalized.species.coef=species.coef,normalized.bias.coef=bias.coef,
               normalized.all.coef=all.coef,normalized.std.errs=std.errs,
               all.coef=unstandardized.coef,
               std.errs=std.errs/sd.normalizer, #### NOTE STD ERRS FOR INTERCEPTS ARE WRONG!!!
               species.coef=unstandardized.species.coef,
               bias.coef=unstandardized.bias.coef,
               species.PA=species.PA,species.PO=species.PO,species.RL=species.RL,
               sdm.margins.ab.RL=sdm.margins.ab.RL,
               linear.fit.PA=linear.fit.PA,fit.PA=fit.PA,
               linear.bias.fit.BG=linear.bias.fit.BG,bias.fit.BG=bias.fit.BG,
               linear.fit.BG=linear.fit.BG,fit.BG=fit.BG,
               linear.fit.RL=linear.fit.RL,fit.RL=fit.RL)

    ## Covariance of the coefficients: the inverse Hessian of the penalized likelihood
    ## at the last IRLS step, i.e. the matrix whose diagonal std.errs already reports.
    if(inverse.hessian) {
        normalized.vcov <- fit$fit$inverse.hessian
        dimnames(normalized.vcov) <- list(names(all.coef),names(all.coef))

        ## all.coef -> unstandardized.coef is a linear map (divide by the sds, then shift
        ## alpha(k) and gamma(k) by the means, exactly as above); A is that map, so the
        ## covariance on the original scale is A V A'.  Unlike std.errs, this is right
        ## for the intercepts too.
        jd <- n.species*(p.sdm+2) + seq_len(p.bias)
        A <- diag(length(all.coef))
        for(k in 1:n.species) {
            jk <- (p.sdm+2)*(k-1) + 1:(p.sdm+1)
            A[jk[1],jk[-1]] <- -sdm.means[-1]        # alpha(k)
            A[jk[1]+p.sdm+1,jd] <- -bias.means       # gamma(k)
        }
        A <- sweep(A,2,sd.normalizer,"/")
        vcov <- A %*% normalized.vcov %*% t(A)
        dimnames(vcov) <- dimnames(normalized.vcov)

        tr$normalized.vcov <- normalized.vcov
        tr$vcov <- vcov
    }

    class(tr) <- c("multispeciesPP","list")
    tr
}
