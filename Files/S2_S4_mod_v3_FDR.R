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

script.start = Sys.time()

xi_mean =  c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,10,10,10)
xi_var  =  rep(4,30)
xi_par = cbind("mean" = xi_mean, "var" = xi_var)
rownames(xi_par) = paste0("comp.",1:30)
tau_max = lambda_max = sigma_max = 5

M_par = c("shape"= 5,"rate"= 2)
pi_par = matrix(c(1.8,0.2), ncol = 2, nrow = 3, byrow = TRUE)

colnames(pi_par) = c("shape1","shape2")
rownames(pi_par) = paste0("loc.", 1:3)

PRIORS_all = list(
  "xi"     = xi_par,
  "tau"    = tau_max,
  "lambda" = lambda_max,
  "nu0"     = 5,
  "sigma0" = sigma_max,
  "M"      = M_par,
  "pi"     = pi_par)


remove(xi_mean,xi_var,xi_par, tau_max,lambda_max, sigma_max, M_par,pi_par)



settings = c(
  lapply(501:600, function(x) list("seed_data" = as.numeric(x), "var_lv" = "low",  "seed_chain"=   55000+x-500,
                                   "min" = 0.25, "max"= 0.75,    "eff.min" = .1, "eff.max" = 1)),
  lapply(601:700, function(x) list("seed_data" = as.numeric(x), "var_lv" = "mid", "seed_chain"=  56000+x-600,
                                   "min" = 0.75, "max"= 1.50,    "eff.min" = .1, "eff.max" = 1)),
  lapply(701:800, function(x) list("seed_data" = as.numeric(x), "var_lv" = "high",   "seed_chain"=  57000+x-700,
                                   "min" = 1.50, "max"= 3.00,   "eff.min" = .1, "eff.max" = 1)),
  lapply(801:900, function(x) list("seed_data" = as.numeric(x), "var_lv" = "low",   "seed_chain"=  58000+x-800,
                                   "min" = 0.25, "max"= 0.75,   "eff.min" = 0, "eff.max" = 0)),
  lapply(901:1000, function(x) list("seed_data" = as.numeric(x), "var_lv" = "mid", "seed_chain"= 59000+x-900,
                                   "min" = 0.75, "max"= 1.50,   "eff.min" = 0, "eff.max" = 0)),
  lapply(1001:1100, function(x) list("seed_data" = as.numeric(x), "var_lv" = "high", "seed_chain"= 510000+x-1000,
                                   "min" = 1.50, "max"= 3.00,   "eff.min" = 0, "eff.max" = 0))
)

iterations <- length(settings)
n.cores = 100
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

foreach(it = settings, .options.snow = opts, .packages = c("lme4","reshape2")) %dopar%{
  
  SIMULATION = simulate_data_heteroscedastic_no_groups(sd_min = it$min,   sd_max = it$max,
                                                       min. = it$eff.min, max.   = it$eff.max,
                                                       seed = it$seed_data)
  
  data = SIMULATION$simulated_data
  
  obj_name = paste0("_singleton_",it$var_lv, "_",it$seed_data,".rds")
  
  obbserved_eff = sapply(paste0("comp.",1:30), function(com) {
    HS = data[data$state=="healthy"  & data$variable == com,"value"]
    DS = data[data$state=="diseased" & data$variable == com,"value"]
    tapply(DS-HS, rep(paste0("loc.",1:3),10), mean)})
  
  PRIORS = PRIORS_all
  PRIORS$effects =  t(apply(obbserved_eff, 1, function(row) c("mean" = mean(row),"var" = var(row)/4)))
  
  mu_start = t(sapply(paste0("comp.",1:30), function(com) {
    HS = data[data$state=="healthy"  & data$variable == com,"value"]
    tapply(HS, rep(paste0("loc.",1:3),10), mean)}) )
  
  
  mod = lmer("value~-1+variable*location+state*location-state-location+(1|id)", data = data)
  delta_start = ranef(mod)$id[,"(Intercept)"]
  names(delta_start) = paste0("sub.",1:10)
  
  rho_start = paste0("latent.group.",1:30)
  names(rho_start) = paste0("comp.",1:30)
  
  theta_start = sapply(unique(rho_start), function(x){
    tapply(data$value[data$variable    %in% names(rho_start)[rho_start==x] & data$state == "diseased"] -
             data$value[data$variable    %in% names(rho_start)[rho_start==x] & data$state == "healthy"],
           data$location[data$variable %in% names(rho_start)[rho_start==x] & data$state == "healthy"],
           mean)},  USE.NAMES = TRUE, simplify = FALSE)
  
  
  
  pi_start = c("loc.1" = 0.5, "loc.2" = 0.5, "loc.3" = 0.5)
  xi_start = rowMeans(mu_start)
  tau_start = rep(.5,30)
  lambda_start = .5
  sigma_start = rep(it$min^2,30)
  names(sigma_start) = names(tau_start) = paste0("comp.",1:30)
  sigma0_start = it$min^2
  M_start = 1
  
  starting_values = list( theta.star =  theta_start,
                          rho        =  rho_start,
                          
                          pi         =  pi_start,
                          
                          mu         =  mu_start,
                          delta      =  delta_start, 
                          
                          xi         =  xi_start,
                          tau        =  tau_start,
                          
                          sigma      =  sigma_start, 
                          lambda     =  lambda_start,
                          
                          sigma0     = sigma0_start,
                          M          = M_start)
  
  
  
  t0 = Sys.time()
  chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = data, PRIORS = PRIORS, starting_values = starting_values,
                                                     sample = 2500, thining = 25, burn.in = 100, seed = it$seed_chain,
                                                     verbouse = FALSE, missing = FALSE)
  
  if(it$eff.max!=0){
    
    saveRDS(list("PRIORS"          = PRIORS,
                 "start"           = starting_values,
                 "data"            = data,
                 "chain"           = chain_to_save,
                 "seed_simulation" = it$seed_data,
                 "seed_chain"      = it$seed_chain,
                 "time"            = Sys.time() - t0),
            file = paste0("Scenario_II/Mod_v3_ExFDR/",it$var_lv,"/sim_s2",obj_name))
  }else{
    saveRDS(list("PRIORS"          = PRIORS,
                 "start"           = starting_values,
                 "data"            = data,
                 "chain"           = chain_to_save,
                 "seed_simulation" = it$seed_data,
                 "seed_chain"      = it$seed_chain,
                 "time"            = Sys.time() - t0),
            file = paste0("Scenario_IV/Mod_v3_ExFDR/",it$var_lv,"/sim_s4",obj_name))      
  }
  
  return(NULL)
}

close(pb)
stopCluster(cl) 
end.time = Sys.time() - script.start
cat("Script execution time:", end.time ,units(end.time))