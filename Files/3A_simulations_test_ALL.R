library(foreach)
library(doParallel)
library(iterators)
library(parallel)
library(doSNOW)
library(lme4)
library(rstudioapi)

directory_source = setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory_source = setwd(dirname(rstudioapi::getSourceEditorContext()$path))

setwd(directory_source)

source("Scripts/functions_v3.0.R")


scenarios  = c("Scenario_I","Scenario_II","Scenario_III","Scenario_IV")
models     = c("Mod_v3","Mod_v3_FDR","Mod_hierarchical_ExFDR","Mod_hierarchical")
var_levels = c("low","mid","high")

### Cycle folder 


start_time = Sys.time()

iterations <- 4
n.cores = 4
n.cores = min(n.cores,detectCores(logical = TRUE)-2,iterations)
cl <- makeCluster(n.cores)
registerDoSNOW(cl)

cat("Detected cores:",detectCores(logical = TRUE),
    "\nUsing:",min(n.cores,detectCores(logical = TRUE)-2),
    "\nCheck:", getDoParWorkers(),
    "\nSimulations:",iterations)

results = foreach(scenario = scenarios) %dopar% {
  
  check = TRUE 
  
  for(var in var_levels){
    for(model in models){

      folder_to_check = paste(scenario, model, var , sep = "/")
      results_list    = list.files(folder_to_check)
      
      size_to_check =   length(results_list)

      for(file in results_list){

        stored = readRDS(paste(folder_to_check, file , sep = "/"))
        PRIORS = stored$PRIORS
        
        if(model == "Mod_v3" | model == "Mod_v3_FDR"){
          
          first_sample = cytokines.model_gibbs.sampler_v3.5( data = stored$data,
                                                             PRIORS = stored$PRIORS,
                                                             starting_values = stored$start,
                                                             missing  = FALSE,
                                                             verbouse = FALSE,
                                                             sample   = 1,
                                                             burn.in  = 100,
                                                             thining  = 25,
                                                             seed = stored$seed_chain)
          
          check =  abs(first_sample[[1]]$sigma0 - stored$chain[[1]]$sigma0) < 1e-6 & check
          
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
          
          check = abs(first_sample[[1]]$sigma0 - stored$chain[[1]]$sigma0) < 1e-6 & check
          
        }      
      }
    }
  }
  return(check)
}

end_time = Sys.time() - start_time
close(pb)
