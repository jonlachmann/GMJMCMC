#' @export
predict.gmjmcmc <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), ...) {
  merged <- merge_results(list(object))
  return(predict.gmjmcmc_merged(merged, x, link, quantiles))
}

#' New idea for a more streamlined function...
#' Produces slightly different results from the fun above since this is using all lo.models too.
#' @inheritParams predict.gmjmcmc_merged
#' @param pop The population to use.
predict.gmjmcmc.2 <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), pop = 1, ...) {

  mmodel <- lapply(object[1:8], function (x) x[[pop]])

  # Precalculate the features for the new data (c(0,1...) is because precalc features thinks there is an intercept and y col).
  x.precalc <- precalc.features(cbind(0, 1, x), mmodel$populations)[, -1]
  return(predict.mjmcmc(mmodel, x.precalc, link, quantiles))
}

#' Predict using a BGNLM model.
#'
#' @param object The model to use.
#' @param x The new data to use for the prediction, a matrix where each row is an observation.
#' @param link The link function to use
#' @param quantiles The quantiles to calculate credible intervals for the posterior moddes (in model space).
#' @param ... Not used.
#'
#' @export
predict.gmjmcmc_merged <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), ...) {
  x <- as.matrix(x)
  preds <- list()
  for (i in seq_along(object$results)) {
    preds[[i]] <- list()
    for (j in seq_along(object$results[[i]]$populations)) {
      # Select the models and features to predict from at this iteration
      models <- object$results[[i]]$models[[j]]
      features <- object$results[[i]]$populations[[j]]
      model.probs <- object$results[[i]]$model.probs[[j]]

      # Precalculate the features for the new data (c(0,1...) is because precalc features thinks there is an intercept and y col).
      x.precalc <- precalc.features(cbind(0, 1, x), features)[, -1]

      yhat <- matrix(0, nrow=nrow(x), ncol=length(models))
      for (k in seq_along(models)) {
        # Models which have 0 weight are skipped since they may also be invalid, and would not influence the predictions.
        if (models[[k]]$crit == -.Machine$double.xmax) next
        yhat[, k] <- link(x.precalc[, c(TRUE, models[[k]]$model), drop=FALSE] %*% models[[k]]$coefs)
      }

      mean.pred <- rowSums(yhat %*% diag(as.numeric(model.probs)))
      pred.quant <- apply(yhat, 1, weighted.quantiles, weights=model.probs, prob=quantiles)

      preds[[i]][[j]] <- list(mean=mean.pred, quantiles=pred.quant, weight=object$results[[i]]$pop.weights[j])
    }
  }

  aggr <- list()
  aggr$mean <- 0 * preds[[1]][[1]]$mean
  aggr$quantiles <- 0 * preds[[1]][[1]]$quantiles
  for (i in seq_along(preds)) {
    for (j in seq_along(preds[[i]])) {
      aggr$mean <- aggr$mean + preds[[i]][[j]]$mean * object$results[[i]]$pop.weights[j]
      aggr$quantiles <- aggr$quantiles + preds[[i]][[j]]$quantiles * object$results[[i]]$pop.weights[j]
    }
  }

  return(list(aggr=aggr, preds=preds))
}

#' @export
predict.mjmcmc <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), ...) {
  # Select the models and features to predict from at this iteration
  models <- c(object$models, object$lo.models)[object$model.probs.idx]

  yhat <- matrix(0, nrow=nrow(x), ncol=length(models))
  for (k in seq_along(models)) {
    # Models which have 0 weight are skipped since they may also be invalid, and would not influence the predictions.
    if (models[[k]]$crit == -.Machine$double.xmax) next
    yhat[, k] <- link(x[, c(TRUE, models[[k]]$model), drop=FALSE] %*% models[[k]]$coefs)
  }

  mean.pred <- rowSums(yhat %*% diag(as.numeric(object$model.probs)))
  pred.quant <- apply(yhat, 1, weighted.quantiles, weights = object$model.probs, prob = quantiles)

  return(list(mean = mean.pred, quantiles = pred.quant))
}

#' @export
predict.mjmcmc_parallel <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), ...) {
  max.crits <- numeric()
  for (i in seq_along(object)) {
    max.crits <- c(max.crits, object[[i]]$best.crit)
  }
  max.crit <- max(max.crits)
  result.weights <- exp(max.crits - max.crit) / sum(exp(max.crits - max.crit))

  preds <- list()
  for (i in seq_along(object)) {
    preds[[i]] <- predict.mjmcmc(object[[i]], x, link, quantiles)
  }

  aggr <- list()
  aggr$mean <- 0 * preds[[1]]$mean
  aggr$quantiles <- 0 * preds[[1]]$quantiles
  for (i in seq_along(preds)) {
    aggr$mean <- aggr$mean + preds[[i]]$mean * result.weights[i]
    aggr$quantiles <- aggr$quantiles + preds[[i]]$quantiles * result.weights[i]
  }

  return(list(aggr = aggr, preds = preds))
}

#' @export
predict.gmjmcmc_parallel <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), ...) {
  merged <- merge_results(object, ...)
  predict.gmjmcmc_merged(merged, x, link, quantiles)
}

#' Calculate weighted quantiles
#'
#' @param values The values to use
#' @param weights The weights of the values
#' @param prob The probabilities of the quantiles to use
#'
#' @return Weighted quantiles
weighted.quantiles <- function (values, weights, prob = c(0.025, 0.975)) {
  ordered <- order(values)
  P <- cumsum(weights[ordered])

  iv <- integer(length(prob))
  for (i in seq_along(iv)) {
    iv[i] <- which.max(P >= prob[i])
  }
  {values[ordered]}[iv]
}
