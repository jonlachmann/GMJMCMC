#' rmclapply: Cross-platform mclapply/forking hack for Windows
#'
#' This function applies a function in parallel to a list or vector (`X`) using multiple cores.
#' On Linux/macOS, it uses `mclapply`, while on Windows it uses a hackish version of parallelism.
#' The Windows version is based on `parLapply` to mimic forking following Nathan VanHoudnos.
#' @param ... Additional arguments to pass to `FUN`.
#' @param mc.cores Number of cores to use for parallel processing. Defaults to `detectCores()`.
#' @param verbose Should a message be shown for Windows users? Defaults to FALSE.
#'
#'
#' @return A list of results, with one element for each element of `X`.
#' @export
#'
#' @examples
#' # Simple example: Square numbers from 1 to 10 using parallel processing
#' result <- rmclapply(1:10, function(x) x^2, mc.cores = 2)
#' print(result)
#'
#' # Example with a custom function and more cores
#' my_func <- function(x) { Sys.sleep(1); return(x * 2) }
#' result <- rmclapply(1:5, my_func, mc.cores = 2, verbose = TRUE)
#' print(result)
rmclapply <- function(..., verbose=FALSE, mc.cores=NULL) {
  
  if(is.null(mc.cores) ) {
    size.of.list <- length(list(...)[[1]])
    mc.cores <- min(size.of.list, detectCores())
  }
  
  if (Sys.info()[['sysname']] == 'Windows' & mc.cores > 1) {
      if (verbose) {
        message("Using parallelization hack for Windows with parLapply.")
      }
      
 
    ## N.B. setting outfile to blank redirects output to
    ##      the master console, as is the default with
    ##      mclapply() on Linux / Mac
    cl <- makeCluster( mc.cores, outfile="" )
    
    ## Find out the names of the loaded packages
    loaded.package.names <- c(
      ## Base packages
      sessionInfo()$basePkgs,
      ## Additional packages
      names( sessionInfo()$otherPkgs ))
    
    tryCatch( {
      
      ## Copy over all of the objects within scope to
      ## all clusters.
      this.env <- environment()
      while( identical( this.env, globalenv() ) == FALSE ) {
        clusterExport(cl,
                      ls(all.names=TRUE, env=this.env),
                      envir=this.env)
        this.env <- parent.env(environment())
      }
      clusterExport(cl,
                    ls(all.names=TRUE, env=globalenv()),
                    envir=globalenv())
      
      ## Load the libraries on all the clusters
      ## N.B. length(cl) returns the number of clusters
      parLapply( cl, 1:length(cl), function(xx){
        lapply(loaded.package.names, function(yy) {
          require(yy , character.only=TRUE)})
      })
      
      ## Run the lapply in parallel
      return( parLapply( cl, ...) )
    }, finally = {
      ## Stop the cluster
      stopCluster(cl)
    })
    
    ## Warn the user if they are using Windows
    if( Sys.info()[['sysname']] == 'Windows' & verbose == TRUE){
      message(paste(
        "\n",
        "   *** Microsoft Windows detected ***\n",
        "   \n",
        "   For technical reasons, the MS Windows version of mclapply()\n",
        "   is implemented as a serial function instead of a parallel\n",
        "   function.",
        "   \n\n",
        "   As a quick hack, we replace this serial version of mclapply()\n",
        "   with a wrapper to parLapply() for this R session. Please see\n\n",
        "     http://www.stat.cmu.edu/~nmv/2014/07/14/implementing-mclapply-on-windows \n\n",
        "   for details.\n\n"))
    }
  }else
  {
    return(mclapply(..., mc.cores = mc.cores))
  }
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
  
  results <- rmclapply(seq_len(runs), function (x) {
    gmjmcmc(data = data, loglik.pi = loglik.pi, loglik.alpha = loglik.alpha, transforms = transforms, ...)
  }, mc.cores = cores)
  
  class(results) <- "gmjmcmc_parallel"
  merged <- merge_results(results, merge.options$populations, merge.options$complex.measure, merge.options$tol, data = data)
  return(merged)
}
