# Validation script for fbms.mlik.master
# Tests all supported (family, prior) combinations with required parameters

set.seed(42)  # Ensure reproducibility

# Generate synthetic data
gen_data <- function(family) {
  n <- 50
  p <- 3
  x <- cbind(1, matrix(rnorm(n * p), n, p))  # Include intercept
  beta <- c(1, -1, 0.5, -0.5)
  
  if (family == "gaussian") {
    y <- x %*% beta + rnorm(n, mean = 0, sd = 1000)
  } else if (family == "binomial") {
    prob <- 1 / (1 + exp(-x %*% beta + rnorm(n, mean = 0, sd = 1)))
    y <- rbinom(n, 1, prob)
  } else if (family == "poisson") {
    lambda <- exp(x %*% beta + rnorm(n, mean = 0, sd = 1))
    y <- rpois(n, lambda)
  } else if (family == "gamma") {
    shape <- 2
    rate <- exp(-x %*% beta + rnorm(n, mean = 0, sd = 1))
    y <- rgamma(n, shape = shape, rate = rate)
  } else {
    stop("Unsupported family")
  }
  
  list(y = as.vector(y), x = x)
}

# Define prior lists
glm_and_gaussian_priors <- c("ZS-adapted", "beta.prime", "EB-local", "g-prior", "hyper-g", "hyper-g-n", 
                             "intrinsic", "Jeffreys", "uniform", "benchmark", "robust", "Jeffreys-BIC", 
                             "CH", "tCCH", "TG")
gaussian_only_priors <- c("ZS-null", "ZS-full", "hyper-g-laplace", "AIC", "JZS","EB-global")

glm_priors <- glm_and_gaussian_priors
gaussian_priors <- c(glm_and_gaussian_priors, gaussian_only_priors)

families <- c("gaussian", "binomial", "poisson", "gamma")

# Required parameters for priors
prior_params <- list(
  "g-prior" = list(g = 50,a = 50),
  "hyper-g" = list(a = 3),
  "hyper-g-n" = list(a = 3),
  "ZS-null" = list(a = 3),
  "ZS-full" = list(a = 500),
  "hyper-g-laplace" = list(a = 3),
  "AIC" = list(a = 3),
  "JZS" = list(a = 3),
  "EB-global" = list(a = 3),
  "EB-local" = list(a = 3),
  "CH" = list(a = 1, b = 2, s = 1),
  "tCCH" = list(a = 1, b = 2, s = 0, rho = 1, v = 1, k = 1),
  "TG" = list(a = 2, s = 1)
)

# Testing loop
for (family in families) {
  priors <- if (family == "gaussian") gaussian_priors else glm_priors
  data <- gen_data(family)
  
  cat("\n===== Testing family:", family, "=====")
  
  for (prior in priors) {
    
    print(prior)
    
    params <- list(family = family,beta_prior = list(type = prior))
    
    params_old <- list(family = family,prior_beta = prior)
    
    
    # Add required parameters if applicable
   
    if (prior %in% names(prior_params)) {
      params$beta_prior <- c(params$beta_prior, prior_params[[prior]])
    }
    
    if (prior %in% names(prior_params)) {
      params_old <- c(params_old, prior_params[[prior]])
    }
    
    # Run the model
    tryCatch({
      
      set.seed(1)
      result <- fbms.mlik.master(data$y, data$x, model = c(T, rep(TRUE, ncol(data$x) - 1)), 
                                 complex = list(oc = 1), params = params)
      set.seed(1)
      result.null <- fbms.mlik.master(data$y, data$x, model = c(T, T, rep(FALSE, ncol(data$x) - 2)), 
                        complex = list(oc = 1), params = params)
      set.seed(1)
      result.old <- fbms.mlik.master_old(data$y, data$x, model = c(T, rep(TRUE, ncol(data$x) - 1)), 
                                       complex = list(oc = 1), params = params_old)#
      set.seed(1)
      result.null.old <- fbms.mlik.master_old(data$y, data$x, model = c(T, T, rep(FALSE, ncol(data$x) - 2)), 
                       complex = list(oc = 1), params = params_old)
      
      
      crit_rounded <- round(result$crit - result.null$crit - result.old$crit + result.null.old$crit, 8) 
      coefs_mean <- round(mean(result$coefs) - mean(result.null$coefs) - mean(result.old$coefs) + mean(result.null.old$coefs), 8)
      
      cat(sprintf("\nPrior: %-15s -> crit: %8.4f, mean(coefs): %8.4f", prior, crit_rounded, coefs_mean))
      
      print("Finished")
      
    }, error = function(e) {
      cat(sprintf("\nPrior: %-15s -> ERROR: %s", prior, conditionMessage(e)))
    })
  }
}
