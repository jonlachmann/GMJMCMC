# Title     : Features for use in GMJMCMC
# Objective : Define the features and get a useful way to turn them into a string
# Created by: jonlachmann
# Created on: 2021-02-10

# transform = 0 is a multiplication of the (two) features listed in f
# transform > 0 is a transform of the type denoted of the features listed in f
# f is a list of features involved in the feature
# alpha is the coefficient before each feature listed in f,
# and also possibly one more for an intercept

createFeature <- function (transform, features, alphas=NULL) {
  # Given no alphas, assume no intercept and unit coefficients
  if (is.null(alphas)) alphas <- c(0, rep(1, length(features)))
  if (length(alphas) != (length(features) + 1)) stop("Invalid alpha/feature count")

  # Generate the new feature matrix
  newFeature <- list(matrix(c(transform, rep(NA,length(alphas)-1),
                      NA, 1:(length(features)), alphas), length(alphas)))
  feature <- append(newFeature, features, 0)
  class(feature) <- "feature"
  return(feature)
}

print.feature <- function (feature) {
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
      if (is.na(feat[j,2]) && feat[j,3] != 0) fString <- paste0(fString, feat[j,3], op)
      # Otherwise this is a feature or covariate, do a recursive conversion
      if (!is.na(feat[j,2])) {
        if (feat[j,3] == 1) fString <- paste0(fString, featureToString(feature[[feat[j,2]]]), op)
        else fString <- paste0(fString, feat[j,3], "*", featureToString(feature[[feat[j,2]]]), op)
      }
    }
    fString <- paste0(fString, ")")
  }
  # This is a plain covariate
  else if (is.numeric(feat)) {
    fString <- paste0("x", feat)
  } else stop("Invalid feature structure")
  return(fString)
}

test <- createFeature(1, 1)

transforms <- list("sin", "cos")

test2 <- createFeature(2, list(test, test3), c(2,1,2))
test3 <- createFeature(2, list(1,2,3,6))

print(test3)
