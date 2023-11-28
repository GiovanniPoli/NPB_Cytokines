library(foreach)
library(doParallel)
library(iterators)
library(parallel)
library(doSNOW)
library(rstudioapi)

directory_source = setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory_source = setwd(dirname(rstudioapi::getSourceEditorContext()$path))

setwd(directory_source)

source("Scripts/functions_v3.0.R")


scenarios  = c("Scenario_I","Scenario_II","Scenario_III","Scenario_IV")
models     = c("Mod_v3","Mod_v3_ExFDR","Mod_hierarchical","Mod_hierarchical_ExFDR")
var_levels = c("low","mid","high")

### Cycle folder 
scenario = scenarios[1]
model    = models[1]
var      = var_levels[1]


folder_to_check = paste(scenario, model, var , sep = "/")
results_list    = list.files(folder_to_check)

size_to_check =   length(results_list)


start_time = Sys.time()

iterations <- size_to_check
n.cores = 50
n.cores = min(n.cores,detectCores(logical = TRUE)-2,iterations)
cl <- makeCluster(n.cores)
registerDoSNOW(cl)

cat("Detected cores:",detectCores(logical = TRUE),
    "\nUsing:",min(n.cores,detectCores(logical = TRUE)-2),
    "\nCheck:", getDoParWorkers(),
    "\nSimulations:",iterations)

pb = txtProgressBar(max = iterations, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)


results = foreach(file = results_list, .options.snow = opts) %dopar%{
  
  stored = readRDS(paste(folder_to_check, file , sep = "/"))
  PRIORS = stored$PRIORS
  
  if(model == "Mod_v3" | model == "Mod_v3_ExFDR"){
    first_sample = cytokines.model_gibbs.sampler_v3.5( data = stored$data,
                                        PRIORS = stored$PRIORS,
                                        starting_values = stored$start,
                                        missing  = FALSE,
                                        verbouse = FALSE,
                                        sample   = 1,
                                        burn.in  = 100,
                                        thining  = 25,
                                        seed = stored$seed_chain)
    
    
    
  }else if(model == "Mod_hierarchical_ExFDR" | model == "Mod_hierarchical"){
    first_sample = cytokines.model_gibbs.sampler_hierarchical( data = stored$data,
                                                       PRIORS = stored$PRIORS,
                                                       starting_values = stored$start,
                                                       missing  = FALSE,
                                                       verbouse = FALSE,
                                                       sample   = 1,
                                                       burn.in  = 100,
                                                       thining  = 25,
                                                       seed = stored$seed_chain)

    
  }
  
  return(abs( first_sample[[1]]$sigma0 - stored$chain[[1]]$sigma0) < 1e-6)
}

close(pb)
stopCluster(cl)
time_end = Sys.time()
