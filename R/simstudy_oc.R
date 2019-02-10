#' Function to make a data set of the operating characterstics of a simulation
#' scenario replicate
#' 
#' @param res Result loaded from a batchtools registry
#' @return Long data set columns containing bias, MSE, and 0/1 interval
#' coverage for TV-SACE and RM-SACE at each evaluated time point
#' @export
extract_oc <- function(res) {
  bias <- res$oc$bias
  mse <- bias
  mse[, -1] <- mse[, -1]^2
  cov <- res$oc$ci_coverage
  rs <- reshape2::melt(bias, id = "eval_t", variable.name = "qty", 
                       value.name = "bias")
  rs$mse <- reshape2::melt(mse, id = "eval_t", variable.name = "qty", 
                           value.name = "mse")$mse
  rs$cov <- reshape2::melt(cov, id = "eval_t", variable.name = "qty", 
                           value.name = "cov")$cov
  return(rs)
}



#' Extract job parameters for Done jobs (that would be processed using
#' reduceResultList)
#' 
#' @param reg Batchtools registry with simulation scenarios
#' @return Length-K list of data frames with job parameters, where K is number
#' of successfully completed jobs
#' @export
extract_job_pars <- function(reg) {
  lapply(batchtools::getJobPars(ids = batchtools::findDone(reg = reg),
                                reg = reg)$job.pars,
         FUN = as.data.frame)
}



#' Extract operating characteristics from a registry along with job submission
#' parameters (e.g., seed, n, ...)
#' 
#' @param reg Batchtools registry
#' @return Data frame with columns for bias, mse, and coverage, as well as 
#' job submission parameters
#' @export 
extract_oc_with_pars <- function(reg) {
  # Make data frame of output and combine with job parameters (seed, scenario,
  # etc)
  res <- as.data.frame(data.table::rbindlist(
         Map(cbind, 
             extract_job_pars(reg = reg),
             batchtools::reduceResultsList(fun = extract_oc, reg = reg))
         ))
  
  # Reorder to have seed column be first if it exists
  # Add a false 1st column in case it does not
  if (!("seed" %in% names(res))) {
    res$seed <- NA
  } 
  if (names(res)[1] != "seed") {
    which_seed <- which(names(res) == "seed")
    res <- res[, c(which_seed, (1:NCOL(res))[-which_seed])]
  } 
  return(res)
}



#' Function to summarize operating characteristics across
#' all jobs in a registry
#' 
#' @param reg Batchtools registry containing SCR replicates
#' @return Data table with mean bias, mse, and coverage for
#' each quantity at each time point
#' @export 
get_registry_oc <- function(reg) {
  prior_opt <- getOption("batchtools.progress")
  options(batchtools.progress = FALSE)
  extracted <- extract_oc_with_pars(reg = reg)
  grouping_vars <- names(extracted)[!(names(extracted) %in% 
                                        c("bias", "mse", "cov", "seed"))]
  res <- extracted %>% 
    dplyr::group_by(.dots = grouping_vars) %>%
    dplyr::summarize(bias = mean(bias, na.rm = TRUE),
                     mse = mean(mse, na.rm = TRUE),
                     cov = mean(cov, na.rm = TRUE))
  res$qty <- factor(res$qty, 
                    levels = c("frac_aa", "tv_sace", "rm_sace"), 
                    labels = c("Fraction always-alive at $t$",
                               "$TV\\text{-}SACE(t,t)$",
                               "$RM\\text{-}SACE(t,t)$"))
  res <- as.data.frame(res)
  if (all(c("scenario", "n", "qty", "eval_t") %in% names(res))) {
    res <- res[order(res$scenario, res$n, res$qty, res$eval_t), ]
    res[, c("eval_t", "qty")] <- res[, c("qty", "eval_t")]
    colnames(res)[colnames(res) %in% c("qty", "eval_t")] <- 
      c("qty", "eval_t")
  } 
  options(batchtools.progress = prior_opt)
  return(res)
}



#' Change column name(s) if they exist
#' 
#' @param df Data frame
#' @param old_name Character ot vector of old names 
#' @param old_name Character ot vector of corresponding new names 
#' @return Data frame with renamed columns
#' @export 
change_name_if <- function(df, old_name, new_name) {
  if ((length(old_name) > 1) && (length(old_name) == length(new_name))) {
    for (i in 1:length(old_name)) {
      df <- change_name_if(df, old_name[i], new_name[i])
    }
  } else {
    if (old_name %in% colnames(df)) {
      colnames(df)[colnames(df) == old_name] <- new_name
    }  
  }
  return(df)
}



#' Function to apply some pretty to the operating characteristics table
#' 
#' Used in interior of \code{\link{make_oc_table}}
#' 
#' @param tab Table to prettify
#' @return Data frame with less duplicate text, etc
#' @export
prettify_oc_table <- function(tab) {

  # Drop sample size if constant across replicates
  if (("n" %in% colnames(tab)) && (length(unique(tab$n)) == 1)) {
    # n <- formatC(tab$n[1], big.mark = ",")
    tab <- tab[ , (colnames(tab) != "n")]
  } 
  
  # Cut down on duplicate text
  tab$qty <- as.character(tab$qty)
  tab$scenario[duplicated(tab$scenario)] <- ""
  min_t <- min(tab$eval_t)
  tab$qty[tab$eval_t != min_t] <- ""
  
  # Better column headings
  tab <- change_name_if(tab, 
                        old_name = c("n", "scenario", "qty", "eval_t",
                                     "bias", "mse", "cov"),
                        new_name = c("$n$", "Scenario", "Quantity", "$t$",
                                     "Bias", "MSE", "Coverage"))
  return(tab)
}
  


#' Make operating characteristics table
#' 
#' @param regs List of batchtools registries
#' @inheritDotParams knitr::kable
#' @return
#' @export
make_oc_table <- function(regs, prettify = TRUE, label = "tab:simoc", ...) {
  tab <- as.data.frame(data.table::rbindlist(lapply(regs, 
                                                    FUN = get_registry_oc)))
  if (prettify) {
    tab <- prettify_oc_table(tab)
  }
  tab <- knitr::kable(tab, linesep = "", ...)
  if (is.null(label)) {
    tab <- gsub(x = tab, pattern = "\\label{tab:}", "", fixed = TRUE)
  } else {
    tab <- gsub(x = tab, pattern = "\\label{tab:}", 
                paste0("\\label{", label, "}"), 
                fixed = TRUE)
  }  
  return(tab)
}

