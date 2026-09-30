predict.multispeciesPP <-
function (object, newdata,
              sdm=TRUE, bias=FALSE,
              species=colnames(object$fit.PA),
              dispersion = NULL, terms = NULL,
              na.action = na.pass, ...)
{
    na.act <- object$na.action
    object$na.action <- NULL
    pred <- matrix(0,nrow(newdata),length(species),dimnames=list(rownames(newdata),species))
    if(sdm) {
        ## fits from before sdm.terms was stored fall back to the bare formula
        sdm.terms <- if(is.null(object$sdm.terms)) object$sdm.formula else object$sdm.terms
        sdm.mat <-  model.matrix(sdm.terms,newdata)
        good.rows <- row.names(sdm.mat)
        pred[!(row.names(pred) %in% good.rows),] <- NA
        pred[good.rows,] <- sdm.mat %*% object$species.coef[-nrow(object$species.coef),species]
    }
    if(bias) {
        bias.terms <- if(is.null(object$bias.terms)) object$bias.formula else object$bias.terms
        bias.mat <- model.matrix(bias.terms,newdata)
        good.rows <- rownames(bias.mat)
        pred[!(row.names(pred) %in% good.rows),] <- NA
        ## the bias term is shared by all species: one value per row, recycled down the
        ## columns.  (Adding the n x 1 product directly fails for more than one species,
        ## and bias.mat[good.rows,] drops to a vector for a single bias covariate.)
        pred[good.rows,] <- pred[good.rows,,drop=FALSE] +
            c(bias.mat[good.rows,,drop=FALSE] %*% object$bias.coef)
        pred <- pred + rep(object$species.coef["isPO",species],each=nrow(pred))
    }
    pred
}
