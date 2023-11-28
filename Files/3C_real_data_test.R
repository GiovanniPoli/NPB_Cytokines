# Simulation Results Scripts 21/11/2023
# Giovanni Poli

script.start = Sys.time()

#### Functions & Packages -----

library(mclust)
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

#### DATA IMPORT            ----
cytokines_data <- read_csv("Data/Cytokine_2020.csv", 
                           col_types = cols(state = readr::col_factor(levels = c("healthy", "diseased")),
                                            location = readr::col_factor(levels = c("mucosa", "submucosa","serum")),
                                            id = readr::col_factor(levels = c("IBD1","IBD2","BD5","BD6","BD7","BD8",
                                                                              "BD9","BD10","BD11","BD12","BD20","BD21"))))

lower_limit = c(10,1246,160	,305 ,0.55,	12	,0.73,2.01,	2.26,	6.86,	3.27,2.33,	12	,6.57,	9.69,	20  ,17	,20	,8.84,	7.42,9.23,	2.4	,7.57,1.89,	3.61,	5.1 ,	7.2)
upper_limit = c(41000,5104000,653700,312450,2250,12675,3000,	8250,9250,	28100,13400,9550,	49500,26900,	39700,82500,	68500,	82000,	36200,30400,	37800,	9850,31000,	7750,	14800,5225,	29500)

labels_variable = colnames(cytokines_data)[2:28]
names(upper_limit) = names(lower_limit) = labels_variable

labels_location    = levels(cytokines_data$location)
labels_id          = levels(cytokines_data$id)
labels_state       = levels(cytokines_data$state)

cytokines_data = reshape2::melt(cytokines_data,  id.vars = c("id", "state", "location"))
cytokines_data = na.omit(cytokines_data)

cytokines_data$censored =   cytokines_data$value == 0

cytokines_data$value = log(cytokines_data$value)
cytokines_data$value[cytokines_data$censored] = log(lower_limit)[cytokines_data$variable][cytokines_data$censored]

keep = ls()


#### ACTUAL MODELS          -----

SETTINGS = list(##  Model 4 & EB criteria
  
  list(seed =  1, version = "v3",             start = "singleton"),
  list(seed =  2, version = "v3",             start = "mclust"   ), 
  list(seed = 11, version = "v4_1",           start = "singleton"),
  list(seed = 12, version = "v4_1",           start = "mclust"   ),
  list(seed = 11, version = "v4_2",           start = "singleton"),
  list(seed = 12, version = "v4_2",           start = "mclust"   ),
  list(seed = 21, version = "v3_equal",       start = "singleton"),
  list(seed = 22, version = "v3_equal",       start = "mclust"   ),
  list(seed = 31, version = "v3_quarter",     start = "singleton"),
  list(seed = 32, version = "v3_quarter",     start = "mclust"   ),
  list(seed = 41, version = "v3_doble",       start = "singleton"),
  list(seed = 42, version = "v3_doble",       start = "mclust"   ),
  list(seed = 51, version = "v3_diffuse",     start = "singleton"),
  list(seed = 52, version = "v3_diffuse",     start = "mclust"   ),
  list(seed = 61, version = "hierarchical",  start = "means"    ),
  list(seed = 62, version = "hierarchical",  start = "zeros"    ),
  
  list(seed =  3, version = "v3",             start = "singleton"),
  list(seed =  4, version = "v3",             start = "mclust"   ), 
  list(seed = 13, version = "v4_1",           start = "singleton"),
  list(seed = 14, version = "v4_1",           start = "mclust"   ),
  list(seed = 13, version = "v4_2",           start = "singleton"),
  list(seed = 14, version = "v4_2",           start = "mclust"   ),
  list(seed = 23, version = "v3_equal",       start = "singleton"),
  list(seed = 24, version = "v3_equal",       start = "mclust"   ),
  list(seed = 33, version = "v3_quarter",     start = "singleton"),
  list(seed = 34, version = "v3_quarter",     start = "mclust"   ),
  list(seed = 43, version = "v3_doble",       start = "singleton"),
  list(seed = 44, version = "v3_doble",       start = "mclust"   ),
  list(seed = 53, version = "v3_diffuse",     start = "singleton"),
  list(seed = 54, version = "v3_diffuse",     start = "mclust"   ),
  list(seed = 63, version = "hierarchical",  start = "means"    ),
  list(seed = 64, version = "hierarchical",  start = "zeros"    ),
  #
  #### M values 
  #
  list(seed = 71, version = "v3_M_1_1_4",             start = "singleton"),
  list(seed = 72, version = "v3_M_1_1_4",             start = "mclust"   ),
  list(seed = 73, version = "v3_M_1_1_4",             start = "singleton"),
  list(seed = 74, version = "v3_M_1_1_4",             start = "mclust"   ),
  #
  list(seed = 81, version = "v3_M_10_2",             start = "singleton"),
  list(seed = 82, version = "v3_M_10_2",             start = "mclust"   ),
  list(seed = 83, version = "v3_M_10_2",             start = "singleton"),
  list(seed = 84, version = "v3_M_10_2",             start = "mclust"   ),
  #
  list(seed = 81, version = "v3_M_50_10",             start = "singleton"),
  list(seed = 82, version = "v3_M_50_10",             start = "mclust"   ),
  list(seed = 83, version = "v3_M_50_10",             start = "singleton"),
  list(seed = 84, version = "v3_M_50_10",             start = "mclust"   ),
  
  ### n_sigma
  
  list(seed = 91, version = "v3_ns_1",             start = "singleton"),
  list(seed = 92, version = "v3_ns_1",             start = "mclust"   ),
  list(seed = 93, version = "v3_ns_1",             start = "singleton"),
  list(seed = 94, version = "v3_ns_1",             start = "mclust"   ),
  
  list(seed = 101, version = "v3_ns_15",             start = "singleton"),
  list(seed = 102, version = "v3_ns_15",             start = "mclust"   ),
  list(seed = 103, version = "v3_ns_15",             start = "singleton"),
  list(seed = 104, version = "v3_ns_15",             start = "mclust"   ),
  
  
  list(seed = 111, version = "v3_ns_30",           start = "singleton"),
  list(seed = 112, version = "v3_ns_30",           start = "mclust"   ),
  list(seed = 113, version = "v3_ns_30",           start = "singleton"),
  list(seed = 114, version = "v3_ns_30",           start = "mclust"   ),
  
  ## Pi prior
  
  list(seed = 115, version = "v3_pi_5_15",            start = "singleton"),
  list(seed = 116, version = "v3_pi_5_15",             start = "mclust"   ),
  list(seed = 117, version = "v3_pi_5_15",             start = "singleton"),
  list(seed = 118, version = "v3_pi_5_15",             start = "mclust"   ),
  
  list(seed = 119, version = "v3_pi_1_1",              start = "singleton"),
  list(seed = 120, version = "v3_pi_1_1",              start = "mclust"   ),
  list(seed = 121, version = "v3_pi_1_1",              start = "singleton"),
  list(seed = 122, version = "v3_pi_1_1",              start = "mclust"   ),
  
  list(seed = 123, version = "v3_pi_15_5",             start = "singleton"),
  list(seed = 124, version = "v3_pi_15_5",             start = "mclust"   ),
  list(seed = 125, version = "v3_pi_15_5",             start = "singleton"),
  list(seed = 126, version = "v3_pi_15_5",             start = "mclust"   )
)

v3_ids = c("v3","v4_1","v4_2","v3_equal",   
           "v3_quarter","v3_doble","v3_diffuse",
           "v3_M_1_1_4","v3_M_10_2","v3_M_50_10","v3_ns_1",     
           "v3_ns_15","v3_ns_30","v3_pi_5_15","v3_pi_1_1",
           "v3_pi_15_5")
v4_ids = c("v4_1","v4_2")
vhier  = "hierarchical"

iterations <- length(SETTINGS)
n.cores = 4
n.cores = min(n.cores,detectCores(logical = TRUE)-2,iterations)
cl <- makeCluster(n.cores)
registerDoSNOW(cl)

cat("Detected cores:",detectCores(logical = TRUE),
   "\nUsing:",min(n.cores,detectCores(logical = TRUE)-2),
   "\nCheck:", getDoParWorkers(),
   "\nSimulations:",iterations)
 
results = foreach(it = SETTINGS) %dopar% {
  
  OBJECT = readRDS(paste0("Chains_real_data/",it$version,"/Cytokines_model_chain_",
             it$version,"_",it$start,"_",it$seed,".rds"))
  
  PRIORS          = OBJECT$PRIORS
  starting_values = OBJECT$start
  first_sigma20   = OBJECT$chain.MCMC[[1]]$sigma0 # Use Sigma^2_0 for v3 and hier 
  first_M         = OBJECT$chain.MCMC[[1]]$M      # Use M for v4
  remove(OBJECT)
  
  if(it$version %in% v3_ids){
    
    chain_to_check = cytokines.model_gibbs.sampler_v3.5( data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                         sample = 1, burn.in = 500, thining = 50,  seed = it$seed,
                                                         missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
    return( abs(chain_to_check[[1]]$sigma0 - first_sigma0) < 1e-6 )
    
  }else if(it$version %in% v4_ids){
    
    chain_to_check = cytokines.model_gibbs.sampler_v4.0( data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                         sample = 1, burn.in = 500, thining = 50,  seed = it$seed,
                                                         missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
    return( abs(chain_to_check[[1]]$M - first_M) < 1e-6 )
  
  }else if(it$version %in% vhier){
    
    chain_to_check = cytokines.model_gibbs.sampler_hierarchical( data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                                 sample = 1, burn.in = 500, thining = 50,  seed = it$seed,
                                                                 missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
    return( abs(chain_to_check[[1]]$sigma0 - first_sigma0) < 1e-6 )
  }else{
    cat("error")
  }
  
}

stopCluster(cl)
end.time = Sys.time() - script.start
cat("Script execution time:", end.time ,units(end.time))
