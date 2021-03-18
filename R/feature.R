# Title     : Features for use in GMJMCMC
# Objective : Define the features and get a useful way to turn them into a string
# Created by: jonlachmann
# Created on: 2021-02-10

# transform = 0 is a multiplication of the (two) features listed in f
# transform > 0 is a transform of the type denoted by the transform variable
# f is a list of features involved in the feature
# alpha is the coefficient before each feature listed in f,
# and also possibly one more for an intercept

# A feature matrix has the structure
# | TRANSFORM   WIDTH ALPHA0 |
# | DEPTH       FEAT1 ALPHA1 |
# | NA          FEAT2 ALPHA2 |
# | NA          FEAT3 ALPHA3 |
# ...

#' Create method for "feature" class
#'
#' @param transform A numeric denoting the transform type
#' @param features A list of features to include
#' @param transform A numeric vector denoting the alphas to use
#'
#' @export create.feature
create.feature <- function (transform, features, alphas=NULL) {
  # Given no alphas, assume no intercept and unit coefficients
  if (is.null(alphas)) alphas <- c(0, rep(1, length(features)))
  if (length(alphas) != (length(features) + 1)) stop("Invalid alpha/feature count")
  # Calculate the depth and width of the new feature
  if (transform == 0) {
    depth <- 1 + depth.feature(features[[1]]) + depth.feature(features[[2]])
    width <- 2 + width.feature(features[[1]]) + width.feature(features[[2]])
  }
  else {
    depth <- 0
    width <- length(features)
    for (i in 1:(length(features))) {
      locdepth <- depth.feature(features[[i]])
      locwidth <- width.feature(features[[i]])
      width <- width + locwidth
      if (locdepth > depth) depth <- locdepth
    }
    depth <- depth + 1
  }

  # Generate the new feature matrix
  newFeature <- list(matrix(c(transform, depth, rep(NA,length(alphas)-2),
                      width, 1:(length(features)), alphas), length(alphas)))
  feature <- append(newFeature, features, 0)
  class(feature) <- "feature"
  return(feature)
}

#' Update alphas on a feature
#'
#' @param feature The feature to be updated
#' @param alphas The alphas that will be used
#' @param recurse If we are recursing, to note the number of alphas used
update.alphas <- function (feature, alphas, recurse=FALSE) {
  feat <- feature[[length(feature)]]
  alpha <- 0
  # This is a more complex feature
  if (is.matrix(feat)) {
    # Adjust intercept if it is not multiplication
    if (feat[1,1] > 0 && nrow(feat) > 2) {
      alpha <- alpha + 1
      feat[1,3] <- alphas[alpha]
    }
    for (i in 2:nrow(feat)) {
      # Multiplication does not have alphas, and a zero intercept is no intercept
      if (feat[1,1] > 0 && feat[1,3] != 0) {
        alpha <- alpha + 1
        feat[i,3] <- alphas[alpha]
      }
      # If we have a nested feature, recurse into it
      if (is.list(feature[[feat[i,2]]])) {
        recur <- update.alphas(feature[[feat[i,2]]], alphas[alpha+1:length(alphas)], TRUE)
        feature[[feat[i,2]]] <- recur$feature
        alpha <- alpha + recur$alpha
      }
    }
  }
  feature[[length(feature)]] <- feat
  if (recurse) return(list(alpha=alpha, feature=feature))
  else return(feature)
}

#' Print method for "feature" class
#'
#' @param feature An object of class "feature"
#' @param transforms The list of transforms to use
#' @param dataset Set the regular covariates as columns in a dataset
#' @param alphas Print a "?" instead of actual alphas to prepare the output for alpha estimation
#'
#' @method print feature
#' @export
print.feature <- function (feature, transforms, dataset=F, alphas=F) {
  fString <- ""
  feat <- feature[[length(feature)]]
  # This is a more complex feature
  if (is.matrix(feat)) {
    # Assume that we are not doing multiplication
    op <- "+"
    # Add the outer transform is there is one
    if (feat[1,1] > 0) fString <- paste0(fString, transforms[[feat[1,1]]], "(")
    # If g = 0, we are doing multiplication
    else {
      op <- "*"
      fString <- paste0(fString, "(")
    }
    for (j in 1:nrow(feat)) {
      # No plus or multiplication sign on the last one
      if (j == nrow(feat)) op <- ""
      # If this is an intercept just add it in
      if (j == 1 && feat[j,3] != 0) {
        if (!alphas) fString <- paste0(fString, feat[j,3], op)
        else fString <- paste0(fString, "?", op)
      }
      # Otherwise this is a feature or covariate, do a recursive conversion
      if (j != 1) {
        # Process alphas, which are only present if there is more than one term in the feature
        # this implies that the feature is not a multiplication (i.e. only one _term_).
        if (nrow(feat) > 2 && feat[1,1] > 0) {
          if (alphas) fString <- paste0(fString, "?*")
          else fString <- paste0(fString, feat[j,3], "*")
        }
        fString <- paste0(fString, print.feature(feature[[feat[j,2]]], transforms, dataset, alphas), op)
      }
    }
    fString <- paste0(fString, ")")
  }
  # This is a plain covariate
  else if (is.numeric(feat)) {
    if (dataset) fString <- paste0("data[,", feat+1, "]")
    else fString <- paste0("x", feat)
  } else stop("Invalid feature structure")
  return(fString)
}

# A function to get the depth of a feature
depth.feature <- function (feature) {
  feat <- feature[[length(feature)]]
  if (is.matrix(feat)) return(feat[2,1])
  else if (is.numeric(feat)) return(0)
  else stop("Invalid feature structure")
}

# A function to get the width of a feature
width.feature <- function (feature) {
  feat <- feature[[length(feature)]]
  if (is.matrix(feat)) return(feat[1,2])
  else if (is.numeric(feat)) return(1)
  else stop("Invalid feature structure")
}

# A function to get the complexity measures of a list of features
complex.features <- function (features) {
  featcount <- length(features)
  width <- rep(NA, featcount)
  depth <- rep(NA, featcount)
  for (i in 1:featcount) {
    width[i] <- width.feature(features[[i]])
    depth[i] <- depth.feature(features[[i]])
  }
  return(list(width=width, depth=depth))
}

#' Function to generate a function string for a model consisting of features
#'
#' @param transform A numeric denoting the transform type
#' @param features A list of features to include
#' @param transform A numeric vector denoting the alphas to use
#'
#' @export create.feature
model.function <- function (model, features, transforms, link) {
  modelstring <- paste0(sapply(features[model], print.feature, transforms, alphas=T), collapse="+")
  modelfun <- set_alphas(modelstring)
  modelfun$formula <- paste0(link, "(", modelfun$formula, ")")
  return(modelfun)
}