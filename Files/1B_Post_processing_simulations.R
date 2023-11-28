################################################
######## Simulations results 28/06/2022 ########

# DatD:   30/06/2022
# Author: Giovanni Poli


##### Functions & Library      -----                                                                                             ----

library(ROCR)
library(devtools)
library(mcclust.ext)
library(BiocManager)
library(limma)
library(statmod)
library(fossil)
library(foreach)
library(doParallel)
library(iterators)
library(parallel)
library(doSNOW)
library(lme4)
library(reshape2)
library(sn)
library(rstudioapi)

directory_source = setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory_source = setwd(dirname(rstudioapi::getSourceEditorContext()$path))

setwd(directory_source)

source("Scripts/functions_v3.0.R")
source("Scripts/My_GG_Scripts.R")

PATHS =  list()
pos = 0 

true.eff = matrix( c(rep(c(0,0,1),20), rep(c(1,1,1),10)), ncol=3, byrow= TRUE)
true.clust  = c(rep(1,10),rep(2,10),rep(3,10))
colnames(true.eff) = paste0("loc.",1:3)
rownames(true.eff) = names(true.clust) = paste0("comp.",1:30)


##### For loop to read results -----
t0 = Sys.time()

for(scenario in paste0("Scenario_",c("I","II","III","IV"))){
  models = list.files(paste(directory_source,scenario, sep ="/", collape = "/"))
  for(mod in models){
    var_levels = list.files(paste(directory_source,scenario,mod, sep ="/", collape = "/"))
    for(var in var_levels){
      pos = pos + 1 
      PATHS[[pos]] = list("path"     = paste(directory_source,scenario,mod,var, sep ="/", collape = "/"),
                          "scenario" = scenario,
                          "model"    = mod,
                          "variance" = var)     
    }
  }
}

iterations <- length(PATHS)
n.cores = 6
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


results = foreach(it = PATHS, .options.snow = opts, .packages = c(
  "ROCR",        #
  "devtools",    #
  "mcclust.ext", #
  "BiocManager", #
  "limma",       #
  "statmod")) %dopar%{
  
  DF.ret = data.frame()
  
  for(mod in list.files(it$path)){
    
    sim = readRDS(paste(it$path,mod,  collapse = "/", sep = ""))
    
    if( it$model == "Mod_hierarchical"){ 
    
      
      post_prob   = Reduce("+", lapply(sim$chain, function(x) x$beta == 0) ) / length(sim$chain)
      pred_sp     = prediction(c(post_prob),    c(true.eff))
      AUC.iter    = performance(pred_sp,   "auc")@y.values[[1]]
      

      DF.ret     = rbind(DF.ret,    c(it$model,                     # Model  
                                      it$scenario,                                # Scenario
                                      it$variance,                                    # Ranges of Variance Scale
                                      AUC.iter,                                   # AUC  
                                      NA,                                         # Binder Loss Groups
                                      NA,                                         # VI Groups
                                      NA,                                         # ARI (Binder Loss)                         
                                      NA,                                         # ARI (VI)                 
                                      NA,                                         # Latent groups
                                      sum(post_prob < .1) /90,                    # FPRs3 .1 
                                      sum(post_prob < .2) /90,                    # FPRs3 .2 
                                      sum(post_prob < .3) /90,                    # FPRs3 .3 
                                      sum(post_prob < .4) /90,                    # FPRs3 .4 
                                      sum(post_prob < .5) /90,                    # FPRs3 .5
                                      sim$seed_simulation,
                                      sim$seed_chain
      ))
      
      
      ## NOTE: When analyzing the chain of this model, the frequentist analyses are also performed. 
      
      data = sim$data

      #### LIMMA
      
      Array  = matrix(data[,"value"], ncol=60, nrow=30, byrow= TRUE)
      fix    = paste0(data$state[1:60],data$location[1:60])
      block  = data$id[1:60]
      design = model.matrix(~0 + fix)
      
      ### Duplicate for same subject replications
      rand = duplicateCorrelation(object = Array,
                                  design = design,
                                  block  = block)
      
      fit = lmFit(Array, design, correlation = rand$consensus.correlation)
      
      ### Contrast matrix for inflammed vs healthy 
      
      cont.dif = makeContrasts( loc.1 = fixdiseasedloc.1 -  fixhealthyloc.1,
                                loc.2 = fixdiseasedloc.2 -  fixhealthyloc.2,
                                loc.3 = fixdiseasedloc.3 -  fixhealthyloc.3,
                                levels = design)
      
      fit2 = contrasts.fit(fit, cont.dif)
      fit3 = eBayes(fit2, robust = TRUE)
      
      post_probs_limma = fit3$p.value
      
      rownames(post_probs_limma) = paste0("comp.",1:30)
      
      pred_limma = prediction(c(post_probs_limma), c(true.eff))
      AUC.limma = performance(pred_limma,"auc")@y.values[[1]]
      

      DF.ret     = rbind(DF.ret,    c("Limma",                                    # Model 
                                      it$scenario,                                # Scenario
                                      it$variance,                                    # Ranges of Variance Scale
                                      AUC.limma,                                  # AUC  
                                      NA,                                         # Binder Loss Groups
                                      NA,                                         # VI Groups
                                      NA,                                         # ARI (Binder Loss)                         
                                      NA,                                         # ARI (VI)                 
                                      NA,                                         # Latent groups
                                      sum(post_probs_limma < .1) /90,                    # FPRs3 .1 
                                      sum(post_probs_limma < .2) /90,                    # FPRs3 .2 
                                      sum(post_probs_limma < .3) /90,                    # FPRs3 .3 
                                      sum(post_probs_limma < .4) /90,                    # FPRs3 .4 
                                      sum(post_probs_limma < .5) /90,                    # FPRs3 .5
                                      sim$seed_simulation,
                                      NA
      ))
      
      ## Wilcoxon Rank Sum Test for all pair
      
      WK.pvalues = sapply(paste0("loc.",1:3), function(l) sapply(paste0("comp.",1:30), function(c){
                           test = wilcox.test(x = data$value[ data$location == l & data$variable == c & data$state == "diseased"],
                           y = data$value[ data$location == l & data$variable == c & data$state == "healthy"],
                           paired = TRUE)
        return(test$p.value)
        }  ))
      
      pred_WK = prediction(c(WK.pvalues), c(true.eff))
      AUC.WK  = performance(pred_WK,"auc")@y.values[[1]]
      
      
      DF.ret     = rbind(DF.ret,    c("Wilcoxon Rank",                      # Model 
                                      it$scenario,                          # Scenario      
                                      it$variance,                          # Ranges of Variance Scale
                                      AUC.WK,                                     # AUC  
                                      NA,                                         # Binder Loss Groups
                                      NA,                                         # VI Groups
                                      NA,                                         # ARI (Binder Loss)                         
                                      NA,                                         # ARI (VI)                 
                                      NA,                                         # Latent groups
                                      sum(WK.pvalues < .1) /90,                   # FPRs3 .1 
                                      sum(WK.pvalues < .2) /90,                   # FPRs3 .2 
                                      sum(WK.pvalues < .3) /90,                   # FPRs3 .3 
                                      sum(WK.pvalues < .4) /90,                   # FPRs3 .4 
                                      sum(WK.pvalues < .5) /90,                   # FPRs3 .5
                                      sim$seed_simulation,
                                      NA
      ))

      
      colnames(DF.ret) = c("Model",
                           "Scenario",
                           "Errors_Scale",
                           "AUC",
                           "groups_BL",
                           "groups_VI",
                           "ARI_BL",
                           "ARI_VI",
                           "H_Ex",
                           "FDR_1",
                           "FDR_2",
                           "FDR_3",
                           "FDR_4",
                           "FDR_5",
                           "seed_data",
                           "seed_simulation")
    }else if(it$model == "Mod_hierarchical_ExFDR"){

      post_prob   = Reduce("+", lapply(sim$chain, function(x) x$beta == 0) ) / length(sim$chain)
      pred_sp     = prediction(c(post_prob),    c(true.eff))
      AUC.iter    = performance(pred_sp,   "auc")@y.values[[1]]
      
      
      DF.ret     = rbind(DF.ret,    c(it$model,                       # Model  
                                      it$scenario,                                # Scenario     
                                      it$variance,                                # Ranges of Variance Scale
                                      AUC.iter,                                   # AUC  
                                      NA,                                         # Binder Loss Groups
                                      NA,                                         # VI Groups
                                      NA,                                         # ARI (Binder Loss)                         
                                      NA,                                         # ARI (VI)                 
                                      NA,                                         # Latent groups
                                      sum(post_prob < .1) /90,                    # FPRs3 .1 
                                      sum(post_prob < .2) /90,                    # FPRs3 .2 
                                      sum(post_prob < .3) /90,                    # FPRs3 .3 
                                      sum(post_prob < .4) /90,                    # FPRs3 .4 
                                      sum(post_prob < .5) /90,                    # FPRs3 .5
                                      sim$seed_simulation,
                                      sim$seed_chain
      ))
      
        colnames(DF.ret) = c("Model",
                           "Scenario",
                           "Errors_Scale",
                           "AUC",
                           "groups_BL",
                           "groups_VI",
                           "ARI_BL",
                           "ARI_VI",
                           "H_Ex",
                           "FDR_1",
                           "FDR_2",
                           "FDR_3",
                           "FDR_4",
                           "FDR_5",
                           "seed_data",
                           "seed_simulation")
      
    }else{
      ### SEMI PARAMETRIC POST PROCESSIONG
      # Mass_point group c(0,0,0) from the multiple groups of DP as one
      
      chain = Mass_point_group(sim$chain)
      Coocorence. = Coocorence(chain) 
      
      cl_binder = minbinder.ext(Coocorence.)
      cl_VI     = minVI(Coocorence.)
      

      post_prob   = post.inference(sim$chain,  plot = FALSE)$zero.probs
      pred_sp     = prediction(c(post_prob),    c(true.eff))
      AUC.iter    = performance(pred_sp,   "auc")@y.values[[1]]
      
      H.it        = mean(sapply(chain, function(x) length(unique(x$rho))))
      
      DF.ret     = rbind(DF.ret,    c(it$model,                                   # Model 
                                      it$scenario,                                # Scenario
                                      it$variance,                                # Ranges of Variance Scale
                                      AUC.iter,                                   # AUC  
                                      length(unique(cl_binder$cl)),               # Binder Loss Groups
                                      length(unique(cl_VI$cl)),                   # VI Groups
                                      NA,                                         # ARI (Binder Loss)                         
                                      NA,                                         # ARI (VI)                 
                                      H.it,                                       # Latent groups
                                      sum(post_prob < .1) /90,                    # FPRs3 .1 
                                      sum(post_prob < .2) /90,                    # FPRs3 .2 
                                      sum(post_prob < .3) /90,                    # FPRs3 .3 
                                      sum(post_prob < .4) /90,                    # FPRs3 .4 
                                      sum(post_prob < .5) /90,                    # FPRs3 .5
                                      sim$seed_simulation,
                                      sim$seed_chain
      ))
      
         colnames(DF.ret) = c("Model",
                           "Scenario",
                           "Errors_Scale",
                           "AUC",
                           "groups_BL",
                           "groups_VI",
                           "ARI_BL",
                           "ARI_VI",
                           "H_Ex",
                           "FDR_1",
                           "FDR_2",
                           "FDR_3",
                           "FDR_4",
                           "FDR_5",
                           "seed_data",
                           "seed_simulation")
    }
  }
  return(DF.ret)
  }
results = Reduce(rbind, results)
close(pb)
stopCluster(cl)

# saveRDS(results, "results.RDS")
# results = readRDS("results.RDS")

##### FIGURE 1                 -----

results = transform(results,
          Model         = factor(Model),
          Scenario      = factor(Scenario),
          Errors_Scale  = factor(Errors_Scale),
          AUC           = as.numeric(AUC),
          groups_BL     = as.numeric(groups_BL),
          groups_VI     = as.numeric(groups_VI),
          ARI_BL        = as.numeric(ARI_BL),
          ARI_VI        = as.numeric(ARI_VI),
          H_Ex          = as.numeric(H_Ex),
          FDR_1         = as.numeric(FDR_1),
          FDR_2         = as.numeric(FDR_2), 
          FDR_3         = as.numeric(FDR_3),    
          FDR_4         = as.numeric(FDR_4), 
          FDR_5         = as.numeric(FDR_5), 
          seed_data     = as.numeric(seed_data), 
          seed_simulation  = as.numeric(seed_simulation))

summary(results)

# Scenario IV: NULL case has no AUC
df.plot.1_1 = results[results$Scenario != "Scenario_IV",]

df.plot.1_1$Scenario =   factor( x      = as.character(df.plot.1_1$Scenario),
                                 levels = c("Scenario_I","Scenario_II","Scenario_III"),
                                 labels = c("Scenario I","Scenario II","Scenario III"))


df.plot.1_1$Model = factor(     x = df.plot.1_1$Model,
                           levels = c("Mod_v3",
                                      "Mod_v3_ExFDR",
                                      "Mod_hierarchical",
                                      "Mod_hierarchical_ExFDR",
                                      "Limma",
                                      "Wilcoxon Rank"),
                           labels = c("Semi-Parametric (Unif.)", 
                                     "Semi-Parametric (FDR)",
                                     "Hierarchical (Unif.)",
                                     "Hierarchical (FDR)",     
                                     "Limma",                
                                     "Wilcoxon Rank"))


df.plot.1_1$Errors_Scale = factor(x = df.plot.1_1$Errors_Scale,
                                levels = c("low","mid","high"),
                                labels  = c("low","mid","high"))

p1 <- ggplot(df.plot.1_1, aes(y=AUC, fill = Model, x = Model, col = Model)) + 
      facet_grid(Scenario ~ Errors_Scale) +  
      geom_violin(alpha=0.1)  +
      geom_boxplot(width=0.2, fill = "white") + theme_bw() +
      scale_y_percent("",values = c(10:5)/10, breaks= c(10:5)/10) +
      coord_cartesian(ylim =  c(5,10)/10) +
      scale_color_discrete("") +
      scale_fill_discrete("")  + 
      xlab(NULL) + ylab("AUC") + ggtitle("AUC") + 
      theme(legend.position = "top", text = element_text(family = "serif"),
            axis.text.x = element_blank())

p1


df.plot.1_2 = results[results$Scenario == "Scenario_IV", c("Model","Errors_Scale","Scenario",
                                                           "FDR_1","FDR_2","FDR_3","FDR_4","FDR_5")]



df.plot.1_2     = melt(df.plot.1_2, id=c("Errors_Scale","Scenario", "Model"))

df.plot.1_2$Errors_Scale = factor(     x = df.plot.1_2$Errors_Scale,
                                  levels = c("low","mid","high"),
                                  labels = c("low","mid","high"))

df.plot.1_2$variable =  factor(x  = df.plot.1_2$variable,
                             levels = c("FDR_1", "FDR_2", "FDR_3", "FDR_4", "FDR_5"),
                             labels = c("10%", "20%", "30%", "40%", "50%"))


df.plot.1_2$Model = factor(   x = df.plot.1_2$Model,
                         levels = c("Mod_v3",
                                    "Mod_v3_ExFDR",
                                    "Mod_hierarchical",
                                    "Mod_hierarchical_ExFDR",
                                    "Limma",
                                    "Wilcoxon Rank"),
                         labels = c("Semi-Parametric (Unif.)", 
                                    "Semi-Parametric (FDR)",
                                    "Hierarchical (Unif.)",
                                    "Hierarchical (FDR)",     
                                    "Limma",                
                                    "Wilcoxon Rank"))


ordered_lv = c(
           "FDR_1Semi-Parametric (Unif.)",
           "FDR_2Semi-Parametric (Unif.)",
           "FDR_3Semi-Parametric (Unif.)",
           "FDR_4Semi-Parametric (Unif.)",
           "FDR_5Semi-Parametric (Unif.)",
           "Blank1",
           "FDR_1Semi-Parametric (FDR)",
           "FDR_2Semi-Parametric (FDR)",
           "FDR_3Semi-Parametric (FDR)",
           "FDR_4Semi-Parametric (FDR)",
           "FDR_5Semi-Parametric (FDR)",
           "Blank2",
           "FDR_1Hierarchical (Unif.)",
           "FDR_2Hierarchical (Unif.)",
           "FDR_3Hierarchical (Unif.)",
           "FDR_4Hierarchical (Unif.)",
           "FDR_5Hierarchical (Unif.)",
           "Blank3",
           "FDR_1Hierarchical (FDR)",
           "FDR_2Hierarchical (FDR)",
           "FDR_3Hierarchical (FDR)",
           "FDR_4Hierarchical (FDR)",
           "FDR_5Hierarchical (FDR)",
           "Blank4",
           "FDR_1Limma",
           "FDR_2Limma",
           "FDR_3Limma",
           "FDR_4Limma",
           "FDR_5Limma",
           "Blank5",
           "FDR_1Wilcoxon Rank",
           "FDR_2Wilcoxon Rank",
           "FDR_3Wilcoxon Rank",
           "FDR_4Wilcoxon Rank",
           "FDR_5Wilcoxon Rank") 


df.plot.1_2$x.f = paste0(df.plot.1_2$variable,df.plot.1_2$x.f$Model) 
df.plot.1_2$x.f = factor(x = df.plot.1_2$x.f,
                       levels = ordered_lv,
                       labels = ordered_lv)
df.plot.1_2$Scenario = factor(x = df.plot.1_2$Scenario,
                            levels = "Scenario_IV",
                            labels = "Scenario IV")                   

p2 <- ggplot(df.plot.1_2, aes(y = value, x = variable, col = Model)) + 
  facet_grid( Scenario ~ Errors_Scale) +  
  geom_boxplot() + theme_bw() +
  scale_y_percent("",values = c(0,2.5,5,7.5,10)/10, breaks= c(0,2.5,5,7.5,10)/10) +
  scale_color_discrete("") +
  scale_fill_discrete("")  +
  theme(legend.position = "top", text = element_text(family = "serif") ) +
  ylab("FDR") + xlab(NULL) + ggtitle("False discovery rate by thresholds") 
p2

legend = get_legend(p1)

plot.export = grid.arrange(legend, p1 + theme(legend.position = "none"), p2 + theme(legend.position = "none"),
                            ncol = 1, nrow = 3, heights  = c(.5,3,1.2))

ggsave("FIG1.PDF", plot = plot.export,
       width  = 21,
       height = 29, 
       units  = "cm",
       dpi    = 800)

##### CREATE TABLE 2           -----
TAB_2 = list()
  
for(SC in c("Scenario_I","Scenario_II","Scenario_III","Scenario_IV")){
    T1 = sapply(unique(results$Errors_Scale), function(V){
      data = results[results$Scenario == SC & results$Errors_Scale == V,]
      tapply(data[,"H_Ex"], data$Model, mean)})[c(4,5),c(2,3,1)]
    T2 = sapply(unique(results$Errors_Scale), function(V){
      data = results[results$Scenario == SC & results$Errors_Scale == V,]
      tapply(data[,"groups_VI"], data$Model, mean)})[c(4,5),c(2,3,1)]
    T3 = sapply(unique(results$Errors_Scale), function(V){
      data = results[results$Scenario == SC & results$Errors_Scale == V,]
      tapply(data[,"groups_BL"], data$Model, mean)})[c(4,5),c(2,3,1)]
    
    ret = round(cbind(T1,T2,T3),2)
    TAB_2[[SC]] = ret
}

TAB_2  = lapply(TAB_2, function(x){ colnames(x) = paste0(rep(c("Low",   "Mid","High"),        3),
                                                         rep(c("E[H|y]","VI" ,"BL"),   each = 3))
                                    rownames(x) = c("Semi-Parametric (Unif.)","Semi-Parametric (FDR)")
                                    return(x)})

TAB_2
