#' Break up existing chunks of registry ids
#' 
#' @param reg batchtools registry object
#' @param n.chunks Number of chunks to be submitted
#' @param chunk.size Size of chunks to be submitted
#' @return Chunked ids to be passed to batchtools::submitJobs
#' @export
chunk_registry <- function(reg, n.chunks = NULL, chunk.size = NULL) {
  ids <- batchtools::getJobTable(reg = reg)
  ids$chunk = batchtools::chunk(ids$job.id, 
                                n.chunks = n.chunks, 
                                chunk.size = chunk.size)
  return(ids)
}



#' Conditionally clear a registry if requested
#' 
#' @param clear_existing Whether to clear registry
#' @param registry Batchtools registry to clear (or not)
#' @return None (if registry empty), or error message (if not empty and not
#' requested to clear). Used for side effects.
#' @export
process_clear_registry <- function(clear_existing, registry) {
  if (clear_existing) {
    batchtools::clearRegistry(reg = registry)
  } else {
    if (NROW(batchtools::getJobTable(reg = registry)) > 0) {
      stop("Registry has existing jobs but clear_existing is FALSE")  
    }
  }
}



#' Function to submit many semicompeting risk scenario simulation study 
#' replicates using a batchtools registry
#' 
#' @param registry batchtools registry object
#' @param scenario Scenario vector (\code{scenario = 1:3} does all 3)
#' @param seed Seed vector (\code{seed = 1:R} does R replicates)
#' @param clear_existing Whether to clear existing registry first
#' @param n Values for sample sizes
#' @param iter Number of MCMC iterations
#' @param chains Number of MCMC chains
#' @param chunk.size How many jobs should be chunked together
#' @param time_each Number of minutes for each job
#' @param memory Memory to allocate at the smallest n
#' @param max.concurrent.jobs Maximum number of jobs at the same time
#' @param submit Whether to actually submit the jobs
#' @param ... Additional arguments to pass to Stan
#' @return None; jobs will be submitted and update in registry
#' @export
submit_scenario_jobs <- function(registry, scenario, seed, 
                                 clear_existing = FALSE,
                                 n = 6000,
                                 iter = 2000, chains = 4,
                                 init_r = 0.5,
                                 shared_beta = TRUE,
                                 init = get_init_truth(scenario = scenario,
                                                       chains = chains,
                                                       with_randomness = TRUE),
                                 eval_t = c(30, 90),
                                 parallelize_chains = FALSE,
                                 chunk.size = 1,
                                 time_each = 120,
                                 memory = 1500,
                                 max.concurrent.jobs = 4000,
                                 submit = FALSE,
                                 ...) {
  
  process_clear_registry(clear_existing, registry)
  
  prog_opt <- getOption("batchtools.progress")
  options(batchtools.progress = FALSE)
  
  # Calculate truths only once
  truths <- summarize_scenario_truths(scenarios = scenario, eval_t = eval_t)
  
  # Make job
  args <- data.table::CJ(seed = seed,
                         scenario = scenario, 
                         n = n)
  batchtools::batchMap(fun = run_replicate_with_oc, 
                       args = args, 
                       more.args = list(shared_beta = shared_beta,
                                        eval_t = eval_t, 
                                        truths = truths,
                                        init = init, init_r = init_r,
                                        iter = iter, chains = chains, 
                                        mc.cores = ifelse(parallelize_chains, 
                                                          chains, 
                                                          1),
                                        ...), 
                       reg = registry)
  
  walltime <- 60 * time_each * chunk.size
  if (submit) {
    batchtools::submitJobs(ids = chunk_registry(reg = registry,
                                                chunk.size = chunk.size),
                           reg = registry,
                           resources = list(walltime = walltime,
                                            memory = memory,
                                            max.concurrent.jobs = 
                                              max.concurrent.jobs,
                                            ncpus = ifelse(parallelize_chains, 
                                                           chains, 
                                                           1),
                                            ntasks = 1))
    
    # Reset option
    options(batchtools.progress = prog_opt)    
  }
  else {
    message("Jobs ready to be submitted.")
  }
}
