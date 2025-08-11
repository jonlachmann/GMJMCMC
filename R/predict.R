#' Predict responses from a BGNLM model
#'
#' This function generates predictions from a fitted \code{bgnlm_model} object given a new dataset.
#'
#' @param object A fitted \code{bgnlm_model} object obtained from the BGNLM fitting procedure. 
#'              It should contain the estimated coefficients in \code{model$coefs}.
#' @param x A \code{data.frame} containing the new data for which predictions are to be made. 
#'          The variables in \code{x} must match the features used in the model.
#' @param link A link function to apply to the linear predictor. 
#'             By default, it is the identity function \code{function(x)\{x\}}, 
#'             but it can be any function such as \code{plogis} for logistic regression models.
#' @param x_train Training design matrix to be provided when imputations are to be made from them
#' @param ...  Additional arguments to pass to prediction function.
#'
#' @return A numeric vector of predicted values for the given data \code{x}. 
#'         These predictions are calculated as \eqn{\hat{y} = \text{link}(X \beta)}, 
#'         where \eqn{X} is the design matrix and \eqn{\beta} are the model coefficients.
#'
#' @examples
#' \dontrun{
#' # Example with simulated data
#' set.seed(123)
#' x_train <- data.frame(PlanetaryMassJpt = rnorm(100), RadiusJpt = rnorm(100))
#' model <- list(
#'   coefs = c(Intercept = -0.5, PlanetaryMassJpt = 0.2, RadiusJpt = -0.1),
#'   class = "bgnlm_model"
#' )
#' class(model) <- "bgnlm_model"
#'
#' # New data for prediction
#' x_new <- data.frame(PlanetaryMassJpt = c(0.1, -0.3), RadiusJpt = c(0.2, -0.1))
#'
#' # Predict using the identity link (default)
#' preds <- predict.bgnlm_model(model, x_new)
#' }
#'
#' @export
predict.bgnlm_model <- function(object, x, link = function(x) {x}, x_train = NULL, ... ) {
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x_test = x, x_train = x_train)
  x <- data.frame(x)
  if (object$needs.precalc) {
    if (object$intercept) {
      x <- cbind(1, x)
    }
    
    if(length(object$features)==0)
    {
      warning("MPM has no featres included! All posteriors below 0.5! Baseline only used.")
      x.precalc <-  model.matrix(~1, data = x)
    } else precalc <- precalc.features(list(x = x, y = NULL, fixed = object$fixed), object$features)
    
    if (dim(precalc$x)[2]>length(object$coefs[object$coefs!=0])) {
      precalc$x <- as.matrix(precalc$x[,-1])
    }
    
    yhat <- link(precalc$x %*% object$coefs[object$coefs != 0])
  } else {
    
    if(length(object$coefs)==1)
    {
      warning("MPM has no featres included! All posteriors below 0.5! Baseline only used.")
      x.precalc <-  model.matrix(~1, data = x)
    }
    else{
      x.precalc <- model.matrix(
        as.formula(paste0("~I(", paste0(names(object$coefs)[-1][object$coefs[-1]!=0], collapse = ")+I("), ")")),
        data = x
      )
    }
    
    if (dim(x.precalc)[2]<length(object$coefs[object$coefs!=0])) {
      x.precalc <- cbind(1,x.precalc)
    } else if (dim(x.precalc)[2]>length(object$coefs[object$coefs!=0])) {
      x.precalc <- as.matrix(x.precalc[,-1])
    }
    yhat <- link(x.precalc %*% object$coefs[object$coefs!=0])
  }
  return(yhat)
}


#' Predict using a gmjmcmc result object.
#'
#' @inheritParams predict.gmjmcmc_merged
#' @return A list containing aggregated predictions and per model predictions.
#' \item{aggr}{Aggregated predictions with mean and quantiles.}
#' \item{preds}{A list of lists containing individual predictions per model per population in object.}
#' 
#' @examples
#' result <- gmjmcmc(
#'  x = matrix(rnorm(600), 100),
#'  y = matrix(rnorm(100), 100),
#'  P = 2,
#'  transforms = c("p0", "exp_dbl")
#' )
#' preds <- predict(result, matrix(rnorm(600), 100))
#' 
#' 
#' @export
predict.gmjmcmc <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975),  pop = NULL,tol =  0.0000001, x_train = NULL, ...) {
  transforms.bak <- set.transforms(object$transforms)
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x, x_train)
  
  merged <- merge_results(list(object), data = list(x = x, object$fixed), populations = pop, tol = tol)
  set.transforms(transforms.bak)
  return(predict.gmjmcmc_merged(merged, x, link, quantiles))
}

#' New idea for a more streamlined function...
#' Produces slightly different results from the fun above since this is using all lo.models too.
#' @inheritParams predict.gmjmcmc_merged
#' @param pop The population to use.
#' @noRd
predict.gmjmcmc.2 <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), pop = 1, x_train = NULL, ...) {
  transforms.bak <- set.transforms(object$transforms)
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x, x_train)
  
  mmodel <- lapply(object[1:8], function (x) x[[pop]])
  
  # Precalculate the features for the new data (c(0,1...) is because precalc features thinks there is an intercept and y col).
  x.precalc <- precalc.features(cbind(0, 1, x), mmodel$populations)[, -1]
  set.transforms(transforms.bak)
  return(predict.mjmcmc(mmodel, x.precalc, link, quantiles))
}

#' Predict using a merged gmjmcmc result object.
#'
#' @param object The model to use.
#' @param x The new data to use for the prediction, a matrix where each row is an observation.
#' @param link The link function to use
#' @param quantiles The quantiles to calculate credible intervals for the posterior modes (in model space).
#' @param pop The population to plot, defaults to last
#' @param tol The tolerance to use for the correlation when finding equivalent features, default is 0.0000001
#' @param x_train Training design matrix to be provided when imputations are to be made from them
#' 
#' @param ... Not used.
#' @return A list containing aggregated predictions and per model predictions.
#' \item{aggr}{Aggregated predictions with mean and quantiles.}
#' \item{preds}{A list of lists containing individual predictions per model per population in object.}
#'
#' @examples
#' result <- gmjmcmc.parallel(
#'  runs = 1,
#'  cores = 1,
#'  x = matrix(rnorm(600), 100),
#'  y = matrix(rnorm(100), 100),
#'  P = 2,
#'  transforms = c("p0", "exp_dbl")
#' )
#' preds <- predict(result, matrix(rnorm(600), 100))
#'
#' @export
predict.gmjmcmc_merged <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), pop = NULL, tol = 0.0000001, x_train = NULL, ...) {
  
  
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x, x_train)
  if (object$intercept) {
    x <- cbind(1, x)
  }
  
  
  transforms.bak <- set.transforms(object$transforms)
  if (!is.null(pop))
    object <- merge_results(object$results.raw, pop, 2, tol, data = list(x = x, fixed = object$fixed))
  
  preds <- list()
  for (i in seq_along(object$results)) {
    preds[[i]] <- list()
    for (j in seq_along(object$results[[i]]$populations)) {
      # Select the models and features to predict from at this iteration
      models <- object$results[[i]]$models[[j]]
      features <- object$results[[i]]$populations[[j]]
      model.probs <- object$results[[i]]$model.probs[[j]]
      
      # Precalculate the features for the new data
      x.precalc <- precalc.features(list(x = x, fixed = object$fixed), features)$x
      
      yhat <- matrix(0, nrow = nrow(x), ncol = length(models))
      for (k in seq_along(models)) {
        # Models which have 0 weight are skipped since they may also be invalid, and would not influence the predictions.
        if (models[[k]]$crit == -.Machine$double.xmax) next
        yhat[, k] <- link(x.precalc[, c(rep(TRUE, object$fixed), models[[k]]$model), drop=FALSE] %*% models[[k]]$coefs)
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
  set.transforms(transforms.bak)
  return(list(aggr = aggr, preds = preds))
}

#' Predict using a mjmcmc result object.
#'
#' @inheritParams predict.gmjmcmc_merged
#' @return A list containing aggregated predictions.
#' \item{mean}{Mean of aggregated predictions.}
#' \item{quantiles}{Quantiles of aggregated predictions.}
#' 
#' @examples
#' result <- mjmcmc(
#' x = matrix(rnorm(600), 100),
#' y = matrix(rnorm(100), 100),
#' loglik.pi = gaussian.loglik)
#' preds <- predict(result, matrix(rnorm(600), 100))
#' 
#' @export
predict.mjmcmc <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), x_train = NULL, ...) {
  # Select the models and features to predict from at this iteration
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x, x_train)
  
  
  if (object$intercept) {
    x <- cbind(1, x)
  }
  
  
  models <- object$models[object$model.probs.idx]
  
  yhat <- matrix(0, nrow = nrow(x), ncol = length(models))
  for (k in seq_along(models)) {
    # Models which have 0 weight are skipped since they may also be invalid, and would not influence the predictions.
    if (models[[k]]$crit == -.Machine$double.xmax) next
    yhat[, k] <- link(x[, c(rep(TRUE, object$fixed), models[[k]]$model), drop=FALSE] %*% models[[k]]$coefs)
  }
  
  mean.pred <- rowSums(yhat %*% diag(as.numeric(object$model.probs)))
  pred.quant <- apply(yhat, 1, weighted.quantiles, weights = object$model.probs, prob = quantiles)
  
  return(list(mean = mean.pred, quantiles = pred.quant))
}

#' Predict using a mjmcmc result object from a parallel run.
#'
#' @inheritParams predict.gmjmcmc_merged
#' @return A list containing aggregated predictions.
#' \item{mean}{Mean of aggregated predictions.}
#' \item{quantiles}{Quantiles of aggregated predictions.}
#' 
#' @examples
#' result <- mjmcmc.parallel(runs = 1, 
#' cores = 1, 
#' x = matrix(rnorm(600), 100),
#' y = matrix(rnorm(100), 100), 
#' loglik.pi = gaussian.loglik)
#' preds <- predict(result, matrix(rnorm(600), 100))
#' 
#' @export
predict.mjmcmc_parallel <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), x_train = NULL, ...) {
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x, x_train)
  
  max.crits <- sapply(object$chains, function (x) x$best.crit)
  max.crit <- max(max.crits)
  result.weights <- exp(max.crits - max.crit) / sum(exp(max.crits - max.crit))
  
  preds <- lapply(object$chains, predict.mjmcmc,x, link, quantiles)
  
  aggr <- list()
  aggr$mean <- 0 * preds[[1]]$mean
  aggr$quantiles <- 0 * preds[[1]]$quantiles
  for (i in seq_along(preds)) {
    aggr$mean <- aggr$mean + preds[[i]]$mean * result.weights[i]
    aggr$quantiles <- aggr$quantiles + preds[[i]]$quantiles * result.weights[i]
  }
  
  return(list(aggr = aggr, preds = preds))
}

#' Predict using a gmjmcmc result object from a parallel run.
#'
#' @inheritParams predict.gmjmcmc_merged
#' @param ... Additional arguments to pass to merge_results.
#' @return A list containing aggregated predictions and per model predictions.
#' \item{aggr}{Aggregated predictions with mean and quantiles.}
#' \item{preds}{A list of lists containing individual predictions per model per population in object.}
#' 
#' @examples
#' result <- gmjmcmc.parallel(
#'  runs = 1,
#'  cores = 1,
#'  x = matrix(rnorm(600), 100),
#'  y = matrix(rnorm(100), 100),
#'  P = 2,
#'  transforms = c("p0", "exp_dbl")
#' )
#' preds <- predict(result, matrix(rnorm(600), 100))
#' 
#' @export
predict.gmjmcmc_parallel <- function (object, x, link = function(x) x, quantiles = c(0.025, 0.5, 0.975), x_train = NULL, ...) {
  transforms.bak <- set.transforms(object$transforms)
  if(is.null(x_train))
    x <- impute_x(object, x)
  else
    x <- impute_x_pred(object, x, x_train)
  merged <- merge_results(object,data = cbind(1, x), ...)
  results <- predict.gmjmcmc_merged(merged, x, link, quantiles)
  set.transforms(transforms.bak)
  return(results)
}

#' Calculate weighted quantiles
#'
#' @param values The values to use
#' @param weights The weights of the values
#' @param prob The probabilities of the quantiles to use
#'
#' @return Weighted quantiles
#' @noRd
weighted.quantiles <- function (values, weights, prob = c(0.025, 0.975)) {
  ordered <- order(values)
  P <- cumsum(weights[ordered])
  
  iv <- integer(length(prob))
  for (i in seq_along(iv)) {
    iv[i] <- which.max(P >= prob[i])
  }
  {values[ordered]}[iv]
}

impute_x <- function (object, x) {
  if (!is.null(attr(object, which = "imputed"))) {
    df <- data.frame(x)
    na.matr <- data.frame(1 * (is.na(df)))
    cm <- colMeans(na.matr)
    na.matr <- na.matr[, attr(object, which = "imputed")]
    names(na.matr) <- paste0("mis_", names(na.matr))
    for (i in which(cm != 0)){
      med <- median(df[[i]], na.rm = TRUE)
      if(is.na(med))
        stop("No data for imputation in test set, provide x_train in predict!")
      df[[i]][is.na(df[[i]])] <- med
    }
    return(as.matrix(data.frame(df,na.matr)))
  }
  return(as.matrix(x))
}


impute_x_pred <- function (object, x_test, x_train) {
  if (!is.null(attr(object, which = "imputed"))) {
    df <- data.frame(x_test)
    x_train <- data.frame(x_train)
    na.matr <- data.frame(1 * (is.na(df)))
    cm <- colMeans(na.matr)
    na.matr <- na.matr[, attr(object, which = "imputed")]
    names(na.matr) <- paste0("mis_", names(na.matr))
    for (i in which(cm != 0)){
      med <- median(x_train[[i]], na.rm = TRUE)
      if(is.na(med))
      {
        warning("One or more missing in test columns do not have any data in x_train, test set will be used for imputations!")
        med <-  median(df[[i]], na.rm = TRUE)
      }
      df[[i]][is.na(df[[i]])] <- med
    }
    return(as.matrix(data.frame(df,na.matr)))
  }
  return(as.matrix(x_test))
}
