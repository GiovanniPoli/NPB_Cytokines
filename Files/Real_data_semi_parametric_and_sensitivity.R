# Simulation Results Scripts 06/07/2022
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

#### PRIORS                 ----
xi_par = cbind( (log(upper_limit) + 2*log(lower_limit)) /3 , (( log(upper_limit) - log(lower_limit)) /2/qnorm(.975))^2 )
colnames(xi_par) = c("mean","var")

tau_max    = lambda_max = sigma_max = 5
sigma_fix  = .04
M_par = c("shape" = 5,"rate" = 2)

pi_par = matrix(c(1.8,.2), ncol=2, nrow = 3, byrow = TRUE)
colnames(pi_par) = c("shape1" , "shape2")
rownames(pi_par) = labels_location

nu0 = 5

obbserved_eff = sapply(labels_variable, function(com) {
  index.eff = cytokines_data[cytokines_data$state=="healthy"  & cytokines_data$variable == com,"location"]
  HS = cytokines_data[cytokines_data$state=="healthy"  & cytokines_data$variable == com,"value"]
  DS = cytokines_data[cytokines_data$state=="diseased" & cytokines_data$variable == com,"value"]
  tapply(DS-HS, index.eff, mean)})


eff_par = cbind("mean" = apply(obbserved_eff, 1, mean), "var" = apply(obbserved_eff, 1, var)/4)

PRIORS.v3 = list(xi      =  xi_par,
                 tau     =  tau_max,
                 lambda  =  lambda_max,
                 M       =  M_par,
                 pi      =  pi_par,
                 sigma0  =  sigma_max,
                 nu0     =  nu0,
                 effects =  eff_par)

PRIORS.v4 = list(xi      =  xi_par,
                 tau     =  tau_max,
                 lambda  =  lambda_max,
                 M       =  M_par,
                 pi      =  pi_par,
                 sigma   =  sigma_fix,
                 effects =  eff_par)

PRIORS.hierarchical = list(xi      =  xi_par,
                            tau     =  tau_max,
                            lambda  =  lambda_max,
                            pi      =  pi_par,
                            sigma0  =  sigma_max,
                            nu0     =  nu0,
                            effects =  eff_par)


#### STARTING VALUES        -----

J = length(labels_variable) 

mu_start = t(sapply(labels_variable, function(com) {
  index.eff = cytokines_data[cytokines_data$state=="healthy"  & cytokines_data$variable == com,"location"]
  HS = cytokines_data[cytokines_data$state=="healthy"  & cytokines_data$variable == com,"value"]
  tapply(HS, index.eff, mean)}) )


mod = lmer("value~-1+variable*location+state*location-state-location+(1|id)", data = cytokines_data )
delta_start = ranef(mod)$id[,"(Intercept)"]
names(delta_start) = labels_id

rho_start_singleton = paste0("latent.group.",1:J)
rho_start_mclust    = paste0("latent.group.",Mclust(t(obbserved_eff),G = 3:9,"VII")$classification)
names(rho_start_singleton) = names(rho_start_singleton) = labels_variable

theta_start_singleton = sapply(unique(rho_start_singleton), function(x){
  tapply(cytokines_data$value[cytokines_data$variable    %in% names(rho_start_singleton)[rho_start_singleton==x] & cytokines_data$state == "diseased"] -
           cytokines_data$value[cytokines_data$variable    %in% names(rho_start_singleton)[rho_start_singleton==x] & cytokines_data$state == "healthy"],
         cytokines_data$location[cytokines_data$variable %in% names(rho_start_singleton)[rho_start_singleton==x] & cytokines_data$state == "healthy"],
         mean)},  USE.NAMES = TRUE, simplify = FALSE)
theta_start_mclust = sapply(unique(rho_start_mclust), function(x){
  tapply(cytokines_data$value[cytokines_data$variable    %in% names(rho_start_singleton)[rho_start_singleton==x] & cytokines_data$state == "diseased"] -
           cytokines_data$value[cytokines_data$variable    %in% names(rho_start_singleton)[rho_start_singleton==x] & cytokines_data$state == "healthy"],
         cytokines_data$location[cytokines_data$variable %in% names(rho_start_singleton)[rho_start_singleton==x] & cytokines_data$state == "healthy"],
         mean)},  USE.NAMES = TRUE, simplify = FALSE)


beta_start_means = t(simplify2array(theta_start_singleton))
rownames(beta_start_means) = labels_variable

pi_start = c(0.5,0.5,0.5)
names(pi_start) = labels_location
xi_start = rowMeans(mu_start)
tau_start = rep(.5,J)
lambda_start = .5

sigma_start = rep(mod@devcomp$cmp["sigmaREML"],J)
sigma_start_v4 = (rep(sigma_fix,J))^2 - 1e-5 
names(sigma_start) = names(sigma_start_v4) = names(tau_start) = labels_variable
sigma0_start = c(mod@devcomp$cmp["sigmaREML"], use.names = FALSE)
M_start = 1

starting_values_v3 = list( theta.star =  theta_start_singleton,
                           rho        =  rho_start_singleton,
                           
                           pi         =  pi_start,
                           
                           mu         =  mu_start,
                           delta      =  delta_start, 
                           
                           xi         =  xi_start,
                           tau        =  tau_start,
                           lambda     =  lambda_start,
                           
                           sigma      = sigma_start, 
                           sigma0     = sigma0_start,
                           
                           M          = M_start)


starting_values_v4 = list( theta.star =  theta_start_singleton,
                           rho        =  rho_start_singleton,
                           
                           pi         =  pi_start,
                           
                           mu         =  mu_start,
                           delta      =  delta_start, 
                           
                           xi         =  xi_start,
                           tau        =  tau_start,
                           lambda     =  lambda_start,
                           
                           sigma      =  sigma_start_v4, 
                           
                           M          = M_start)

starting_values_hierarchical = list( 
  beta =  beta_start_means,
  
  pi         =  pi_start,
  
  mu         =  mu_start,
  delta      =  delta_start, 
  
  xi         =  xi_start,
  tau        =  tau_start,
  lambda     =  lambda_start,
  
  sigma      = sigma_start, 
  sigma0     = sigma0_start)

remove(list = ls()[! ls() %in% c("PRIORS.hierarchical", "starting_values_hierarchical",
                                 "PRIORS.v3"           , "starting_values_v3",
                                 "PRIORS.v4"           , "starting_values_v4", 
                                 "theta_start_mclust"  , "rho_start_mclust", keep)])

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

iterations <- length(SETTINGS)
n.cores = 12
n.cores = min(n.cores,detectCores(logical = TRUE)-2,iterations)
cl <- makeCluster(n.cores)
registerDoSNOW(cl)

cat("Detected cores:",detectCores(logical = TRUE),
    "\nUsing:",min(n.cores,detectCores(logical = TRUE)-2),
    "\nCheck:", getDoParWorkers(),
    "\nSimulations:",iterations)

foreach(it = SETTINGS) %dopar%{
  
  if(it$version == "v3"){
    
    
    PRIORS = PRIORS.v3
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_equal"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$effects[,"var"] = PRIORS.v3$effects[,"var"] * 4
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_quarter"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$effects[,"var"] = PRIORS.v3$effects[,"var"] / 4
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_doble"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$effects[,"var"] = PRIORS.v3$effects[,"var"] * 16
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_diffuse"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$effects[,"var"] = 4
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v4_1"){
    
    PRIORS = PRIORS.v4
    PRIORS$sigma = sqrt(.5)
    starting_values = starting_values_v4
    
    if(it$start=="mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    chain_to_save = cytokines.model_gibbs.sampler_v4.0(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v4_2"){
    
    PRIORS = PRIORS.v4
    PRIORS$sigma = sqrt(.1)
    starting_values = starting_values_v4
    
    if(it$start=="mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    chain_to_save = cytokines.model_gibbs.sampler_v4.0(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "hierarchical"){
    
    PRIORS = PRIORS.hierarchical
    starting_values = starting_values_hierarchical
    
    if(it$start == "zeros"){
      starting_values$beta[,] = 0
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_hierarchical(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                               sample = 2500, burn.in = 500, thining = 50, seed = it$seed,
                                                               missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_M_1_1_4"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$M = c("shape" = 1, "rate" = .25)
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_M_10_2"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$M = c("shape" = 10, "rate" = 2)
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_M_50_10"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$M = c("shape" = 50, "rate" = 10)
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_ns_1"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$nu0 = 1
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_ns_15"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$nu0 = 15
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_ns_30"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$nu0 = 30
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_pi_5_15"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$pi[,] = cbind(rep(.5, 3), rep(1.5, 3))
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_pi_1_1"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$pi[,] = cbind(rep(1, 3), rep(1, 3))
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }else if(it$version == "v3_pi_15_5"){
    
    
    PRIORS = PRIORS.v3
    PRIORS$pi[,]   = cbind(rep(1.5, 3), rep(.5, 3))
    starting_values = starting_values_v3
    
    if(it$start == "mclust"){
      starting_values$theta.star = theta_start_mclust
      starting_values$rho[labels_variable] = rho_start_mclust 
    }
    
    
    chain_to_save = cytokines.model_gibbs.sampler_v3.5(data = cytokines_data, PRIORS = PRIORS, starting_values = starting_values,
                                                       sample = 2500, burn.in = 500, thining = 50,  seed = it$seed,
                                                       missing = TRUE, lower_censure = log(lower_limit), verbouse = FALSE)
    
  }
  
  saveRDS(list("PRIORS"        = PRIORS,
               "start"         = starting_values,
               "chain.MCMC"    = chain_to_save,
               "seed"          = it$seed),
          file = paste0("Chains_real_data/",it$version,"/Cytokines_model_chain_",it$version,"_",it$start,"_",it$seed,".rds"))
  
  return(NULL)
  
}

stopCluster(cl) 
end.time = Sys.time() - script.start
cat("Script execution time:", end.time ,units(end.time))