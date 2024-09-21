
#' rmclapply: Cross-platform rmclapply/parLapply function
#'
#' This function applies a function in parallel to a list or vector (`X`) using multiple cores.
#' It automatically chooses the appropriate parallel method based on the operating system.
#' On Windows, it uses `parLapply`, while on Linux/macOS it uses `rmclapply`.
#'
#' @param X A vector or list to which the function will be applied.
#' @param FUN The function to be applied to each element of `X`.
#' @param ... Additional arguments to pass to `FUN`.
#' @param mc.cores Number of cores to use for parallel processing. If not provided, defaults to `detectCores() - 1`.
#'
#' @return A list with the results of applying `FUN` to each element of `X`.
#'
#' @examples
#' # Define a function
#' my_function <- function(x) {
#'   return(x^2)
#' }
#' 
#' # Apply rrmclapply with 2 cores
#' result <- rmclapply(1:10, my_function, mc.cores = 2)
#' print(result)
#'
#' @export
rmclapply <- function(X, FUN, mc.cores = NULL,varlist, ...) {
  # Check the operating system
  os_type <- "windows"#.Platform$OS.type
  
  # Use provided cores or default to detectCores() - 1
  if (is.null(mc.cores)) {
    mc.cores <- detectCores() - 1
  }
  
  if (os_type == "windows") {
    # For Windows, use parLapply
    # Set the future plan
    cl <- makeCluster(mc.cores)
    
    clusterExport(cl, varlist = "FUN", envir = environment())
    
    # Capture additional arguments
    args <- list(...)
    clusterExport(cl, varlist = names(args), envir = environment())
    
    # Call parLapply with the captured arguments
    result <- parLapply(cl, X, function(x) {
      do.call(FUN, c(list(x), args))
    })
    
    stopCluster(cl)
    
    

  } else {
    # For other OS (Linux/macOS), use rmclapply
    result <- mclapply(X, FUN, mc.cores = mc.cores, ... )  # Use the provided or default core count
  }
  
  return(result)
}




#' Run multiple mjmcmc runs in parallel, merging the results before returning.
#' @param runs The number of runs to run
#' @param cores The number of cores to run on
#' @param ... Further params passed to mjmcmc.
#' @return Merged results from multiple mjmcmc runs
#' 
#' @examples
#' result <- mjmcmc.parallel(runs = 1, cores = 1, matrix(rnorm(600), 100), gaussian.loglik)
#' summary(result)
#' plot(result)
#' 
#' @export
mjmcmc.parallel <- function (runs = 2, cores = getOption("mc.cores", 2L), ...) {
  results <- rmclapply(seq_len(runs), function (x) { mjmcmc(...) }, mc.cores = cores)
  class(results) <- "mjmcmc_parallel"
  return(results)
}


#' Run multiple gmjmcmc (Genetically Modified MJMCMC) runs in parallel returning a list of all results.
#' @param runs The number of runs to run
#' @param cores The number of cores to run on
#' @param merge.options A list of options to pass to the [merge_results()] function run after the
#' @inheritParams gmjmcmc
#' @param ... Further params passed to mjmcmc.
#' @return Results from multiple gmjmcmc runs
#' 
#' @examples
#' result <- gmjmcmc.parallel(
#'  runs = 1,
#'  cores = 1,
#'  list(populations = "best", complex.measure = 2, tol = 0.0000001),
#'  matrix(rnorm(600), 100),
#'  P = 2,
#'  gaussian.loglik,
#'  loglik.alpha = gaussian.loglik.alpha,
#'  c("p0", "exp_dbl")
#' )
#' 
#' summary(result)
#' 
#' plot(result)
#' 
#' 
#' @export
gmjmcmc.parallel <- function (runs = 2, cores = getOption("mc.cores", 2L), merge.options = list(populations = "best", complex.measure = 2, tol = 0.0000001), data, loglik.pi = gaussian.loglik, loglik.alpha = gaussian.loglik.alpha(), transforms, ...) {
  options("gmjmcmc-transformations" = transforms)
  
  #to fix
  os_type <- "linux"
  
  if (os_type == "windows") {
    # For Windows, use parLapply
    # Set the future plan
    # Combine additional arguments into a list
    extra_args <- list(...)
    
    # Prepare arguments for gmjmcmc
    gmjmcmc_args <- c(list(data = data, 
                           loglik.pi = loglik.pi, 
                           loglik.alpha = loglik.alpha, 
                           transforms = transforms), 
                           extra_args)
    
    #list2env(extra_args)
    
    # Set up the cluster
    cl <- makeCluster(cores)
    clusterExport(cl, varlist = ls(envir = environment()), envir = environment())
    #clusterExport(cl, varlist = names(gmjmcmc_args), envir = environment())
    clusterExport(cl, varlist = "gmjmcmc_args", envir = environment())
    # Run gmjmcmc in parallel using parLapply
    results <- parLapply(cl, seq_len(runs), function(x) {
      # Create a new environment for each worker
     
      # Call gmjmcmc using do.call in the local environment
      result <- tryCatch({
        do.call(gmjmcmc, gmjmcmc_args,envir = new.env())
      }, error = function(e) {
        cat("Error in iteration", x, ":", e$message, "\n")
        e  # Return NULL in case of an error
      })
      return(result)
    })
    
    print(results)
    # Stop the cluster
    stopCluster(cl)
  }else
  {
    results <- mclapply(seq_len(runs), function (x) {
      gmjmcmc(data = data, loglik.pi = loglik.pi, loglik.alpha = loglik.alpha, transforms = transforms, ...)
    }, mc.cores = cores)
  }
  
  class(results) <- "gmjmcmc_parallel"
  merged <- merge_results(results, merge.options$populations, merge.options$complex.measure, merge.options$tol, data = data)
  return(merged)
}
